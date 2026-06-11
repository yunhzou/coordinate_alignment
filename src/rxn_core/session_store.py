"""MongoDB/GridFS-backed session persistence for the local web app."""
from __future__ import annotations

import json
import mimetypes
import shutil
from datetime import datetime, timezone
from pathlib import Path


def _utcnow():
    return datetime.now(timezone.utc)


def _json_ready(value):
    if isinstance(value, datetime):
        return value.isoformat()
    if isinstance(value, dict):
        return {str(k): _json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(v) for v in value]
    return value


def _safe_session_id(value):
    text = "".join(
        ch if ch.isalnum() or ch in "._-" else "_"
        for ch in str(value or "")
    ).strip("._-")
    return text or "session"


class MongoSessionStore:
    """Store webapp sessions in MongoDB and run artifacts in GridFS."""

    def __init__(self, uri="mongodb://localhost:27017",
                 db_name="rxn_core", *, server_timeout_ms=1500):
        try:
            import gridfs
            from pymongo import MongoClient
        except Exception as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                "Mongo session persistence requires pymongo. Install the "
                "optional sessions dependency first."
            ) from exc
        self.client = MongoClient(
            uri, serverSelectionTimeoutMS=int(server_timeout_ms))
        self.client.admin.command("ping")
        self.db = self.client[db_name]
        self.sessions = self.db["webapp_sessions"]
        self.fs = gridfs.GridFS(self.db, collection="webapp_session_files")
        self.sessions.create_index([("updated_at", -1)])

    def available(self):
        try:
            self.client.admin.command("ping")
            return True
        except Exception:
            return False

    def list_sessions(self, limit=50):
        docs = self.sessions.find(
            {},
            {
                "state": False,
                "artifact_ids": False,
            },
        ).sort("updated_at", -1).limit(int(limit))
        return [_session_summary(doc) for doc in docs]

    def save_session(self, name, state, *, job_dir=None, session_id=None):
        now = _utcnow()
        session_id = _safe_session_id(session_id or f"{name}-{now.timestamp()}")
        existing = self.sessions.find_one({"_id": session_id}) or {}
        for item in existing.get("artifact_ids", []):
            file_id = item.get("file_id")
            if file_id is not None:
                try:
                    self.fs.delete(file_id)
                except Exception:
                    pass

        artifact_ids = []
        if job_dir is not None and Path(job_dir).exists():
            artifact_ids = self._store_artifacts(session_id, Path(job_dir))

        doc = {
            "_id": session_id,
            "name": str(name or "session"),
            "state": _json_ready(state or {}),
            "artifact_ids": artifact_ids,
            "updated_at": now,
            "created_at": existing.get("created_at", now),
            "artifact_count": len(artifact_ids),
        }
        self.sessions.replace_one({"_id": session_id}, doc, upsert=True)
        return _session_summary(doc)

    def load_session(self, session_id, *, restore_root=None):
        doc = self.sessions.find_one({"_id": str(session_id)})
        if not doc:
            raise KeyError(f"session not found: {session_id}")
        restored = None
        if restore_root is not None:
            restored = self.restore_artifacts(session_id, restore_root)
        return {
            **_session_summary(doc),
            "state": _json_ready(doc.get("state") or {}),
            "restored": restored,
        }

    def delete_session(self, session_id):
        doc = self.sessions.find_one({"_id": str(session_id)})
        if not doc:
            return False
        for item in doc.get("artifact_ids", []):
            file_id = item.get("file_id")
            if file_id is not None:
                try:
                    self.fs.delete(file_id)
                except Exception:
                    pass
        self.sessions.delete_one({"_id": str(session_id)})
        return True

    def restore_artifacts(self, session_id, restore_root):
        doc = self.sessions.find_one({"_id": str(session_id)})
        if not doc:
            raise KeyError(f"session not found: {session_id}")
        job_id = f"session_{_safe_session_id(session_id)}"
        target_root = Path(restore_root) / job_id
        if target_root.exists():
            shutil.rmtree(target_root)
        target_root.mkdir(parents=True, exist_ok=True)
        paths = []
        for item in doc.get("artifact_ids", []):
            rel_path = item.get("path")
            file_id = item.get("file_id")
            if not rel_path or file_id is None:
                continue
            out = target_root / rel_path
            out.parent.mkdir(parents=True, exist_ok=True)
            handle = self.fs.get(file_id)
            out.write_bytes(handle.read())
            paths.append(str(rel_path))
        viewer_rel = _find_first(paths, suffix="/view.html") or _find_first(
            paths, exact="view.html")
        stage_rel = _find_first(paths, exact="rp_stage.json")
        return {
            "job_id": job_id,
            "artifact_paths": paths,
            "viewer_url": (
                f"/artifacts/{job_id}/{viewer_rel}" if viewer_rel else None),
            "stage_url": (
                f"/artifacts/{job_id}/{stage_rel}" if stage_rel else None),
        }

    def _store_artifacts(self, session_id, job_dir):
        artifact_ids = []
        for path in sorted(Path(job_dir).rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(job_dir).as_posix()
            content_type = (
                mimetypes.guess_type(path.name)[0] or
                ("application/json" if path.suffix == ".json"
                 else "application/octet-stream")
            )
            file_id = self.fs.put(
                path.read_bytes(),
                filename=f"{session_id}/{rel}",
                metadata={
                    "session_id": session_id,
                    "path": rel,
                    "content_type": content_type,
                },
            )
            artifact_ids.append({
                "path": rel,
                "file_id": file_id,
                "content_type": content_type,
                "bytes": int(path.stat().st_size),
            })
        return artifact_ids


def _find_first(paths, *, suffix=None, exact=None):
    for path in paths:
        if exact is not None and path == exact:
            return path
        if suffix is not None and path.endswith(suffix):
            return path
    return None


def _session_summary(doc):
    return {
        "id": str(doc.get("_id")),
        "name": doc.get("name"),
        "created_at": _json_ready(doc.get("created_at")),
        "updated_at": _json_ready(doc.get("updated_at")),
        "artifact_count": int(doc.get("artifact_count", 0)),
    }


def dumps_session_state(state):
    """Stable JSON helper used by tests and debugging scripts."""
    return json.dumps(_json_ready(state), sort_keys=True)
