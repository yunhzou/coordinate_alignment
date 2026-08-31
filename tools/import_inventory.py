#!/usr/bin/env python3
"""Normalize inventory workbooks and resolve structures through PubChem."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from zipfile import ZipFile

from rdkit import Chem


XML_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
PUBCHEM_PROPERTIES = "SMILES,ConnectivitySMILES,InChIKey,IUPACName"
NORMALIZED_FIELDS = (
    "source_file", "barcode", "product_name", "cas_number", "vendor",
    "product_number", "stock_number", "location", "location_path",
    "amount_remaining", "unit", "amount_kg", "storage_requirement",
    "pubchem_cid", "smiles", "connectivity_smiles", "inchikey",
    "iupac_name", "resolution_method", "resolution_query",
    "resolution_status",
)


def _column_index(reference):
    letters = re.match(r"[A-Z]+", reference).group()
    value = 0
    for letter in letters:
        value = value * 26 + ord(letter) - 64
    return value - 1


def _xlsx_rows(path):
    with ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(XML_NS + "si"):
                shared.append("".join(
                    node.text or "" for node in item.iter(XML_NS + "t")))
        root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        for row in root.findall(".//" + XML_NS + "row"):
            values = {}
            for cell in row.findall(XML_NS + "c"):
                node = cell.find(XML_NS + "v")
                value = "" if node is None else (node.text or "")
                if cell.get("t") == "s" and value:
                    value = shared[int(value)]
                values[_column_index(cell.get("r"))] = value.strip()
            yield values


def _read_workbook(path):
    rows = iter(_xlsx_rows(path))
    header = next(rows)
    columns = {name: index for index, name in header.items()}
    for row in rows:
        yield {
            "source_file": path.name,
            "barcode": row[columns["Barcode#"]],
            "product_name": row[columns["ProductName"]],
            "cas_number": row.get(columns["CAS#"], ""),
            "vendor": row.get(columns["Vendor"], ""),
            "product_number": row.get(columns["ProductNo."], ""),
            "stock_number": row.get(columns["StockNumber"], ""),
            "location": row.get(columns["Location"], ""),
            "location_path": row.get(columns["Location Path"], ""),
            "amount_remaining": row.get(columns["AmountRemaining"], ""),
            "unit": row.get(columns["U-O-M"], ""),
            "amount_kg": row.get(columns["Amt Rem (kg)"], ""),
            "storage_requirement": row.get(columns["StorageReq."], ""),
        }


def _clean_name(name):
    name = re.sub(r"\s*,?\s*\d+(?:\.\d+)?\s*%.*$", "", name).strip()
    return name.rstrip(" -")


def _lookup_key(row):
    if row["cas_number"]:
        return "cas", row["cas_number"]
    return "name", _clean_name(row["product_name"])


class PubChemResolver:
    def __init__(self, cache_path, requests_per_second=4.0):
        self.cache_path = cache_path
        self.interval = 1.0 / requests_per_second
        self.lock = threading.Lock()
        self.next_request = 0.0
        self.cache = {}
        if cache_path.exists():
            with cache_path.open(encoding="utf-8") as stream:
                for line in stream:
                    record = json.loads(line)
                    if not (record["status"].startswith("http_")
                            or record["status"] == "request_error"):
                        self.cache[(record["method"], record["query"])] = record

    def _wait(self):
        with self.lock:
            now = time.monotonic()
            delay = max(0.0, self.next_request - now)
            self.next_request = max(now, self.next_request) + self.interval
        if delay:
            time.sleep(delay)

    def resolve(self, key):
        if key in self.cache:
            return self.cache[key]
        method, query = key
        self._wait()
        encoded = urllib.parse.quote(query, safe="")
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{encoded}/property/{PUBCHEM_PROPERTIES}/JSON"
        )
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                properties = json.load(response)["PropertyTable"]["Properties"]
            if len(properties) != 1:
                status = "ambiguous"
                compound = {}
            else:
                status = "resolved"
                compound = properties[0]
        except urllib.error.HTTPError as error:
            status = "not_found" if error.code == 404 else f"http_{error.code}"
            compound = {}
        except (TimeoutError, urllib.error.URLError):
            status = "request_error"
            compound = {}
        record = {
            "method": method,
            "query": query,
            "status": status,
            "cid": compound.get("CID", ""),
            "smiles": compound.get("SMILES", ""),
            "connectivity_smiles": compound.get("ConnectivitySMILES", ""),
            "inchikey": compound.get("InChIKey", ""),
            "iupac_name": compound.get("IUPACName", ""),
        }
        return record

    def store(self, record):
        key = record["method"], record["query"]
        if key in self.cache:
            return
        self.cache[key] = record
        with self.cache_path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record, separators=(",", ":")) + "\n")


def _canonical_smiles(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"PubChem returned invalid SMILES: {smiles!r}")
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def _write_inventory(rows, resolutions, output):
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED_FIELDS)
        writer.writeheader()
        for row in rows:
            method, query = _lookup_key(row)
            resolution = resolutions[(method, query)]
            output_row = dict(row)
            output_row.update({
                "pubchem_cid": resolution["cid"],
                "smiles": (_canonical_smiles(resolution["smiles"])
                           if resolution["smiles"] else ""),
                "connectivity_smiles": resolution["connectivity_smiles"],
                "inchikey": resolution["inchikey"],
                "iupac_name": resolution["iupac_name"],
                "resolution_method": method,
                "resolution_query": query,
                "resolution_status": resolution["status"],
            })
            writer.writerow(output_row)


def _write_structure_bank(rows, output):
    groups = defaultdict(list)
    for row in rows:
        if row["resolution_status"] == "resolved":
            groups[row["smiles"]].append(row)
    fields = (
        "SMILES", "Inventory ID", "Container Count", "Product Names",
        "CAS Numbers", "PubChem CIDs", "Vendors",
    )
    with gzip.open(output, "wt", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for index, (smiles, items) in enumerate(sorted(groups.items()), 1):
            writer.writerow({
                "SMILES": smiles,
                "Inventory ID": f"INVENTORY-{index:06d}",
                "Container Count": len(items),
                "Product Names": " | ".join(sorted({
                    item["product_name"] for item in items})),
                "CAS Numbers": " | ".join(sorted({
                    item["cas_number"] for item in items
                    if item["cas_number"]})),
                "PubChem CIDs": " | ".join(sorted({
                    item["pubchem_cid"] for item in items})),
                "Vendors": " | ".join(sorted({
                    item["vendor"] for item in items if item["vendor"]})),
            })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    input_paths = sorted(Path(args.input_dir).glob("*.xlsx"))
    rows = [row for path in input_paths for row in _read_workbook(path)]
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    resolver = PubChemResolver(output_dir / "pubchem_cache.jsonl")
    keys = sorted({_lookup_key(row) for row in rows})
    missing = [key for key in keys if key not in resolver.cache]
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for index, record in enumerate(executor.map(resolver.resolve, missing), 1):
            resolver.store(record)
            if index % 100 == 0 or index == len(missing):
                print(f"resolved {index}/{len(missing)} new identifiers", flush=True)

    resolutions = {key: resolver.cache[key] for key in keys}
    normalized_path = output_dir / "inventory_normalized.csv"
    _write_inventory(rows, resolutions, normalized_path)
    with normalized_path.open(encoding="utf-8") as stream:
        normalized_rows = list(csv.DictReader(stream))
    bank_path = output_dir / "inventory_structure_bank.csv.gz"
    _write_structure_bank(normalized_rows, bank_path)
    status_counts = defaultdict(int)
    for row in normalized_rows:
        status_counts[row["resolution_status"]] += 1
    print(json.dumps({
        "containers": len(rows),
        "unique_queries": len(keys),
        "status_counts": dict(status_counts),
        "normalized_inventory": str(normalized_path),
        "structure_bank": str(bank_path),
    }, indent=2))


if __name__ == "__main__":
    main()
