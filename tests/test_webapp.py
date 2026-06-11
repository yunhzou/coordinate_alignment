import numpy as np
import pytest

from rxn_core import webapp


class FakeSessionStore:
    def __init__(self):
        self.docs = {}
        self.saved_job_dir = None

    def list_sessions(self, limit=50):
        return [
            {
                "id": key,
                "name": doc["name"],
                "artifact_count": doc.get("artifact_count", 0),
            }
            for key, doc in sorted(self.docs.items())
        ][:limit]

    def save_session(self, name, state, *, job_dir=None, session_id=None):
        session_id = session_id or "fake-session"
        self.saved_job_dir = job_dir
        doc = {
            "id": session_id,
            "name": name,
            "state": state,
            "artifact_count": 1 if job_dir else 0,
        }
        self.docs[session_id] = doc
        return {k: v for k, v in doc.items() if k != "state"}

    def load_session(self, session_id, *, restore_root=None):
        doc = self.docs[session_id]
        return {
            **doc,
            "restored": {
                "job_id": "session_fake-session",
                "viewer_url": "/artifacts/session_fake-session/views/demo/view.html",
                "stage_url": "/artifacts/session_fake-session/rp_stage.json",
            },
        }

    def delete_session(self, session_id):
        return self.docs.pop(session_id, None) is not None


def test_parse_xyz_allows_blank_comment_line():
    elements, coords = webapp.parse_xyz_text(
        "2\n\nH 0 0 0\nH 0 0 0.74\n")

    assert elements == ["H", "H"]
    assert coords.shape == (2, 3)


def test_parse_wbo_supports_three_atom_xtb_triples_and_matrix():
    triples = webapp.parse_wbo_text("1 2 1.0\n2 3 2.0\n", 3)
    assert triples[0, 1] == pytest.approx(1.0)
    assert triples[1, 2] == pytest.approx(2.0)

    matrix = webapp.parse_wbo_text("0 1 0\n1 0 1\n0 1 0\n", 3)
    assert matrix.tolist() == [
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 0.0],
    ]


def test_infer_single_bond_wbo_from_xyz_distance():
    wbo = webapp.infer_single_bond_wbo(
        ["H", "H"],
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]),
    )

    assert wbo[0, 1] == pytest.approx(1.0)


def test_preview_payload_from_smiles_preserves_formal_wbo():
    pytest.importorskip("rdkit")

    data = webapp.preview_payload({
        "mode": "smiles",
        "reactant": {"smiles": "[O:1][H:2]"},
        "product": {"smiles": "[O-:1].[H+:2]"},
    })

    assert data["reactant"]["elements"] == ["O", "H"]
    assert data["reactant"]["wbo"][0][1] == pytest.approx(1.0)
    assert data["reactant"]["atom_maps"] == {"0": 1, "1": 2}


def test_webapp_2d_result_has_branch_degeneracy_toggle():
    html = webapp.APP_HTML

    assert 'id="show2dDegeneracy"' in html
    assert "branchColorGroups" in html
    assert "branchRDegMap" in html
    assert "branchProductDegMap" in html


def test_webapp_rp_api_writes_existing_viewer(tmp_path):
    pytest.importorskip("rdkit")

    app = webapp.AAMWebApp(tmp_path)
    result = app.api_rp({
        "name": "smiles_deprot",
        "mode": "smiles",
        "reactant": {"smiles": "[O:1][H:2]"},
        "product": {"smiles": "[O-:1].[H+:2]"},
        "config": {"n_seeds": 1},
    })

    assert result["rp"]["stage"] == "rp"
    assert result["rp"]["mechanisms"]
    assert result["viewer_url"].endswith("/view.html")
    assert (tmp_path / result["job_id"] / "rp_stage.json").exists()


def test_webapp_session_api_uses_store_and_job_artifacts(tmp_path):
    store = FakeSessionStore()
    app = webapp.AAMWebApp(tmp_path, session_store=store)
    job_dir = tmp_path / "job1"
    job_dir.mkdir()
    (job_dir / "rp_stage.json").write_text("{}")

    saved = app.api_save_session({
        "id": "session1",
        "name": "demo",
        "job_id": "job1",
        "state": {"name": "demo", "lastRun": {"job_id": "job1"}},
    })
    listed = app.api_list_sessions()
    loaded = app.api_load_session("session1")
    deleted = app.api_delete_session("session1")

    assert saved["session"]["id"] == "session1"
    assert store.saved_job_dir == job_dir
    assert listed["sessions"][0]["name"] == "demo"
    assert loaded["state"]["lastRun"]["job_id"] == "job1"
    assert loaded["restored"]["viewer_url"].endswith("/view.html")
    assert deleted["deleted"] is True


def test_subgraph_payload_matches_drawn_formal_graphs():
    result = webapp.run_subgraph_payload({
        "query": {
            "type": "drawn",
            "graph": {
                "nodes": [
                    {"element": "O", "x": 0.0, "y": 0.0},
                    {"element": "C", "x": 1.0, "y": 0.0},
                ],
                "edges": [{"i": 0, "j": 1, "order": 1.0}],
            },
        },
        "target": {
            "type": "drawn",
            "graph": {
                "nodes": [
                    {"element": "C", "x": 0.0, "y": 0.0},
                    {"element": "O", "x": 1.0, "y": 0.0},
                    {"element": "C", "x": 2.0, "y": 0.0},
                ],
                "edges": [
                    {"i": 0, "j": 1, "order": 1.0},
                    {"i": 1, "j": 2, "order": 1.0},
                ],
            },
        },
        "iso_tol": 0.1,
    })

    assert result["n_matches"] >= 1
    assert any(match["mapping"].get("0") == 1 for match in result["matches"])
