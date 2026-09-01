#!/usr/bin/env python3
"""Build a graph-centered PDF explaining fragment detection and assembly."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from reportlab.lib.colors import HexColor, white
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas

from rxn_core.fragment_matching import FragmentDetectionConfig, detect_fragments
from rxn_core.retrosynthesis import assemble_fragment_cover
from rxn_core.smiles import smiles_to_weighted_graph


PAGE = landscape(A4)
WIDTH, HEIGHT = PAGE
INK = HexColor("#172033")
MUTED = HexColor("#526273")
LIGHT = HexColor("#F3F6F9")
LINE = HexColor("#D7E0E8")
BLUE = HexColor("#2878D0")
ORANGE = HexColor("#E98524")
PURPLE = HexColor("#7857C6")
GREEN = HexColor("#159A6B")
RED = HexColor("#D64045")
CYAN = HexColor("#2AA7A1")
ALERT = HexColor("#B73522")
ALERT_BG = HexColor("#FFF0E8")
GOLD = HexColor("#F4B400")

RD_BLUE = (0.16, 0.47, 0.82)
RD_ORANGE = (0.93, 0.49, 0.13)
RD_PURPLE = (0.49, 0.34, 0.78)
RD_GREEN = (0.08, 0.62, 0.39)
RD_RED = (0.86, 0.18, 0.20)
RD_CYAN = (0.10, 0.66, 0.64)


def _graph(smiles):
    return smiles_to_weighted_graph(smiles, expand_hydrogens=True)


def _detect(source_id, smiles, target_smiles, config):
    return detect_fragments(
        _graph(smiles),
        _graph(target_smiles),
        source_id=source_id,
        config=config,
    )


def _assembly(target_smiles, sources, config):
    results = [
        _detect(source_id, smiles, target_smiles, config)
        for source_id, smiles in sources
    ]
    assembled = assemble_fragment_cover(
        _graph(target_smiles),
        tuple(candidate for result in results for candidate in result.candidates),
        maximum_precursors=len(sources),
        assembly_limit=100,
        require_attachment_bonds=False,
    )
    if not assembled.assemblies:
        raise RuntimeError("example did not produce a complete assembly")
    return results, assembled.assemblies[0]


def _rdkit_png(
        smiles, atom_colors, *, cut_bonds=(), formed_bonds=(),
        width=700, height=370, atom_indices=False):
    molecule = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(molecule)
    heavy_atom_count = molecule.GetNumAtoms()
    colors = {
        int(atom): color for atom, color in atom_colors.items()
        if int(atom) < heavy_atom_count
    }
    bond_colors = {}
    highlighted_bonds = []
    for bond in molecule.GetBonds():
        left = bond.GetBeginAtomIdx()
        right = bond.GetEndAtomIdx()
        if left in colors and right in colors and colors[left] == colors[right]:
            highlighted_bonds.append(bond.GetIdx())
            bond_colors[bond.GetIdx()] = colors[left]
    for left, right in cut_bonds:
        if left >= heavy_atom_count or right >= heavy_atom_count:
            continue
        bond = molecule.GetBondBetweenAtoms(int(left), int(right))
        if bond is not None:
            highlighted_bonds.append(bond.GetIdx())
            bond_colors[bond.GetIdx()] = RD_RED
    for left, right in formed_bonds:
        if left >= heavy_atom_count or right >= heavy_atom_count:
            continue
        bond = molecule.GetBondBetweenAtoms(int(left), int(right))
        if bond is not None:
            highlighted_bonds.append(bond.GetIdx())
            bond_colors[bond.GetIdx()] = RD_GREEN

    drawer = Draw.MolDraw2DCairo(width, height)
    options = drawer.drawOptions()
    options.padding = 0.08
    options.bondLineWidth = 2.2
    options.highlightBondWidthMultiplier = 7
    options.atomHighlightsAreCircles = True
    options.addAtomIndices = atom_indices
    drawer.DrawMolecule(
        molecule,
        highlightAtoms=sorted(colors),
        highlightAtomColors=colors,
        highlightBonds=sorted(set(highlighted_bonds)),
        highlightBondColors=bond_colors,
        highlightAtomRadii={atom: 0.27 for atom in colors},
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _source_colors(candidate, color):
    return {atom: color for atom in candidate.retained_atoms}


def _target_colors_by_source(assembly, color_by_source):
    output = {}
    for candidate in assembly.candidates:
        color = color_by_source[candidate.source_id]
        output.update({
            atom: color for atom in candidate.covered_target_atoms
        })
    return output


def _fragment_colors(candidate, colors):
    source = {}
    target = {}
    mapping = dict(candidate.mapping)
    for fragment, color in zip(candidate.retained_fragments, colors):
        for atom in fragment:
            source[atom] = color
            target[mapping[atom]] = color
    return source, target


class Report:
    def __init__(self, output):
        self.output = output
        self.pdf = canvas.Canvas(str(output), pagesize=PAGE)
        self.page_number = 0

    def new_page(self, title, subtitle=None):
        if self.page_number:
            self.pdf.showPage()
        self.page_number += 1
        self.pdf.setFillColor(INK)
        self.pdf.setFont("Helvetica-Bold", 23)
        self.pdf.drawString(38, HEIGHT - 46, title)
        if subtitle:
            self.pdf.setFillColor(MUTED)
            self.pdf.setFont("Helvetica", 10.5)
            self.pdf.drawString(40, HEIGHT - 63, subtitle)
        self.pdf.setStrokeColor(LINE)
        self.pdf.line(38, HEIGHT - 73, WIDTH - 38, HEIGHT - 73)

    def finish_page(self):
        self.pdf.setFillColor(MUTED)
        self.pdf.setFont("Helvetica", 8)
        self.pdf.drawString(38, 18, "Building-block recommendation workflow")
        self.pdf.drawRightString(WIDTH - 38, 18, f"Page {self.page_number}")

    def save(self):
        self.pdf.save()

    def text(self, x, y, text, *, width=340, size=10.5, leading=14,
             color=INK, bold=False):
        self.pdf.setFillColor(color)
        self.pdf.setFont("Helvetica-Bold" if bold else "Helvetica", size)
        words = text.split()
        lines = []
        line = ""
        for word in words:
            trial = word if not line else f"{line} {word}"
            if stringWidth(trial, "Helvetica", size) <= width:
                line = trial
            else:
                lines.append(line)
                line = word
        if line:
            lines.append(line)
        for item in lines:
            self.pdf.drawString(x, y, item)
            y -= leading
        return y

    def bullets(self, x, y, items, *, width=350, size=10.2, leading=14):
        for item in items:
            self.pdf.setFillColor(BLUE)
            self.pdf.circle(x + 3, y + 3, 2.2, fill=1, stroke=0)
            y = self.text(
                x + 13, y, item, width=width - 13, size=size,
                leading=leading)
            y -= 5
        return y

    def box(self, x, y, width, height, title, body="", *, color=BLUE):
        self.pdf.setFillColor(LIGHT)
        self.pdf.setStrokeColor(LINE)
        self.pdf.roundRect(x, y, width, height, 9, fill=1, stroke=1)
        self.pdf.setFillColor(color)
        self.pdf.roundRect(x, y + height - 8, width, 8, 4, fill=1, stroke=0)
        self.pdf.setFillColor(INK)
        self.pdf.setFont("Helvetica-Bold", 12)
        self.pdf.drawString(x + 13, y + height - 29, title)
        if body:
            self.text(
                x + 13, y + height - 48, body,
                width=width - 26, size=9.3, leading=12, color=MUTED)

    def arrow(self, x1, y1, x2, y2, *, color=MUTED, width=1.8):
        self.pdf.setStrokeColor(color)
        self.pdf.setFillColor(color)
        self.pdf.setLineWidth(width)
        self.pdf.line(x1, y1, x2, y2)
        angle = 6
        if abs(x2 - x1) >= abs(y2 - y1):
            direction = 1 if x2 > x1 else -1
            self.pdf.line(x2, y2, x2 - direction * angle, y2 + 4)
            self.pdf.line(x2, y2, x2 - direction * angle, y2 - 4)
        else:
            direction = 1 if y2 > y1 else -1
            self.pdf.line(x2, y2, x2 + 4, y2 - direction * angle)
            self.pdf.line(x2, y2, x2 - 4, y2 - direction * angle)

    def image(self, png, x, y, width, height):
        self.pdf.drawImage(
            ImageReader(io.BytesIO(png)), x, y, width, height,
            preserveAspectRatio=True, anchor="c", mask="auto")

    def label(self, x, y, text, color):
        self.pdf.setFillColor(color)
        self.pdf.circle(x + 5, y + 4, 5, fill=1, stroke=0)
        self.pdf.setFillColor(INK)
        self.pdf.setFont("Helvetica", 9.5)
        self.pdf.drawString(x + 16, y, text)


def _cover_page(report):
    report.new_page(
        "Purpose: Find Building Blocks for P_target",
        "Use a reaction-bank pool to identify structural foundations for retrosynthesis")
    report.box(
        42, 410, 758, 72, "Purpose",
        "Given P_target, find feasible building blocks whose conserved molecular pieces can serve as the foundation for retrosynthesis.",
        color=BLUE)

    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 14)
    report.pdf.drawString(42, 382, "General flow")
    report.box(
        48, 290, 175, 70, "Reaction bank",
        "Pool of candidate building blocks.", color=CYAN)
    report.box(
        48, 205, 175, 65, "P_target",
        "Desired product structure.", color=ORANGE)
    report.box(
        305, 245, 185, 95, "Geometric matching",
        "Find coherent R fragments in P_target and assemble complete target coverage.",
        color=PURPLE)
    report.box(
        580, 245, 215, 95, "Candidate foundations",
        "Rank feasible building-block sets and retain distinct construction patterns.",
        color=GREEN)
    report.arrow(225, 325, 300, 300, color=CYAN, width=2.3)
    report.arrow(225, 237, 300, 280, color=ORANGE, width=2.3)
    report.arrow(495, 292, 575, 292, color=GREEN, width=2.5)

    report.pdf.setFillColor(ALERT_BG)
    report.pdf.setStrokeColor(ALERT)
    report.pdf.roundRect(42, 82, 758, 90, 10, fill=1, stroke=1)
    report.pdf.setFillColor(ALERT)
    report.pdf.roundRect(42, 162, 758, 10, 5, fill=1, stroke=0)
    report.pdf.setFont("Helvetica-Bold", 13)
    report.pdf.drawString(58, 140, "CURRENT SCOPE: PURE MOLECULAR GEOMETRY")
    report.text(
        58, 118,
        "The workflow identifies structurally compatible building blocks. It does not yet prove chemical reactivity or construct a complete multi-step synthesis route.",
        width=720, size=10.3, leading=14, color=INK)
    report.finish_page()


def _architecture_page(report):
    report.new_page(
        "Objective and general approach",
        "Recommend available starting materials whose conserved pieces can account for a desired product")
    y = 380
    report.box(42, y, 150, 95, "Inventory",
               "Available compounds represented as molecular graphs.", color=CYAN)
    report.box(226, y, 175, 95, "Detect useful pieces",
               "Find coherent regions of each compound that occur in the target.", color=BLUE)
    report.box(435, y, 155, 95, "Combine pieces",
               "Select complementary candidates that account for the target.", color=PURPLE)
    report.box(624, y, 175, 95, "Recommend",
               "Present a small, diverse set of plausible precursor sets.", color=ORANGE)
    report.arrow(194, y + 48, 222, y + 48)
    report.arrow(403, y + 48, 431, y + 48)
    report.arrow(592, y + 48, 620, y + 48)

    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 14)
    report.pdf.drawString(42, 330, "What the workflow is trying to preserve")
    report.bullets(44, 305, [
        "Large, recognizable molecular pieces are preferred over scattered atom-by-atom coincidences.",
        "Each product atom is assigned to one proposed source, giving the recommendation clear provenance.",
        "Different ways of dividing the target remain visible as distinct construction patterns.",
        "Every 2D color comes from the computed atom correspondence rather than manual annotation.",
    ], width=745)
    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(42, 82, 757, 105, 10, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 13)
    report.pdf.drawString(60, 160, "Important scope")
    report.text(60, 138,
        "This is a structural recommendation system. It proposes coherent atom sources and implied bond edits; it does not by itself prove kinetic feasibility, select catalysts, or guarantee a one-step laboratory reaction.",
        width=715, size=10.5, leading=15)
    report.finish_page()


def _detection_page(report):
    report.new_page(
        "Stage 1: detect useful molecular pieces",
        "Compare each available compound with the target and retain coherent shared regions")
    xs = (80, 330, 580)
    titles = (
        "A. Find the main shared piece", "B. Account for the remainder",
        "C. Record a colored candidate")
    bodies = (
        "Locate the largest connected region of the source that appears in the desired product.",
        "Separate material that is not retained, while allowing another coherent source piece to match elsewhere.",
        "Store which source atoms correspond to which target atoms. These mappings drive the colors.",
    )
    colors = (CYAN, PURPLE, ORANGE)
    for x, title, body, color in zip(xs, titles, bodies, colors):
        report.box(x, 355, 185, 135, title, body, color=color)
    for left, right in zip(xs, xs[1:]):
        report.arrow(left + 188, 422, right - 3, 422)

    report.pdf.setFont("Helvetica-Bold", 14)
    report.pdf.setFillColor(INK)
    report.pdf.drawString(42, 325, "Graph interpretation")
    report.pdf.setStrokeColor(LINE)
    report.pdf.setLineWidth(2)
    node_positions = [(95, 245), (155, 245), (215, 245), (275, 245)]
    for index, (x, y) in enumerate(node_positions):
        if index:
            report.pdf.line(node_positions[index - 1][0] + 13, y, x - 13, y)
        report.pdf.setFillColor(CYAN if index < 3 else LIGHT)
        report.pdf.setStrokeColor(INK)
        report.pdf.circle(x, y, 13, fill=1, stroke=1)
    report.pdf.setStrokeColor(RED)
    report.pdf.setLineWidth(4)
    report.pdf.line(228, 245, 262, 245)
    report.label(50, 205, "main shared piece", CYAN)
    report.label(220, 205, "boundary of that piece", RED)
    report.label(365, 205, "unused source material", PURPLE)
    report.label(595, 205, "another recognized piece", ORANGE)
    report.text(42, 165,
        "The important idea is fragment-level conservation. A source can contribute more than one coherent piece, but atoms that do not belong in the product are not forced into arbitrary product positions.",
        width=745, size=10.5, leading=15)
    report.finish_page()


def _co2_page(report, config):
    target = "O=C(O)Cc1ccccc1"
    result = _detect("carbon dioxide", "O=C=O", target, config)
    candidate = max(result.candidates, key=lambda item: item.retained_size)
    source_colors, target_colors = _fragment_colors(
        candidate, (RD_CYAN, RD_ORANGE, RD_PURPLE))
    source_png = _rdkit_png(
        "O=C=O", source_colors, cut_bonds=candidate.boundary_bonds,
        width=500, height=280)
    target_png = _rdkit_png(
        target, target_colors, width=700, height=330)

    report.new_page(
        "Example 3: one building block contributes multiple pieces",
        "CO2 is recognized within the carboxyl group of phenylacetic acid")
    report.box(42, 405, 315, 80, "R: carbon dioxide",
               "The carbonyl piece is recognized first. The second oxygen is then retained as another useful piece of the same source.", color=CYAN)
    report.box(485, 405, 315, 80, "P_target: phenylacetic acid",
               "Only the colored carboxyl atoms are owned by this candidate. The aromatic side of P remains available to other sources.", color=ORANGE)
    report.image(source_png, 55, 225, 285, 170)
    report.image(target_png, 440, 215, 375, 190)
    report.arrow(360, 322, 425, 322, color=CYAN, width=2.5)
    report.arrow(360, 298, 425, 298, color=ORANGE, width=2.5)
    report.pdf.setFillColor(MUTED)
    report.pdf.setFont("Helvetica", 8)
    report.pdf.drawCentredString(392, 337, "R fragments mapped into P_target")
    report.label(60, 205, "initial retained fragment", CYAN)
    report.label(240, 205, "additional target-owned fragment", ORANGE)
    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(42, 82, 758, 95, 10, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 11)
    report.pdf.drawString(58, 152, "What the colors show")
    report.bullets(58, 130, [
        "The cyan carbonyl is the main connected match.",
        "The orange oxygen is recognized as a second useful piece from the same R.",
        "Together, the colored atoms account for the CO2 contribution; the uncolored target remains available to other building blocks.",
    ], width=720, size=9.2, leading=11)
    report.finish_page()


def _assembly_page(report):
    report.new_page(
        "Stage 2: assemble a complete product explanation",
        "Combine complementary pieces so every target atom has one proposed source")
    report.box(45, 405, 185, 90, "Candidate index",
               "Collect useful target regions found across the inventory.", color=PURPLE)
    report.box(328, 405, 185, 90, "Construction patterns",
               "Explore distinct ways of dividing the target among sources.", color=BLUE)
    report.box(611, 405, 185, 90, "Recommendation",
               "Show the best precursor sets for each useful partition.", color=ORANGE)
    report.arrow(235, 450, 322, 450)
    report.arrow(518, 450, 605, 450)

    report.pdf.setFont("Helvetica-Bold", 13)
    report.pdf.setFillColor(INK)
    report.pdf.drawString(48, 345, "A target atom has exactly one owner")
    centers = [(115, 270), (185, 270), (255, 270), (325, 270),
               (395, 270), (465, 270), (535, 270), (605, 270)]
    module_colors = (BLUE, BLUE, BLUE, BLUE, ORANGE, ORANGE, ORANGE, ORANGE)
    for index, ((x, y), color) in enumerate(zip(centers, module_colors)):
        if index:
            report.pdf.setStrokeColor(GREEN if index == 4 else INK)
            report.pdf.setLineWidth(4 if index == 4 else 1.5)
            report.pdf.line(centers[index - 1][0] + 14, y, x - 14, y)
        report.pdf.setFillColor(color)
        report.pdf.setStrokeColor(INK)
        report.pdf.circle(x, y, 14, fill=1, stroke=1)
        report.pdf.setFillColor(white)
        report.pdf.setFont("Helvetica-Bold", 8)
        report.pdf.drawCentredString(x, y - 3, str(index))
    report.label(75, 220, "atoms owned by source 1", BLUE)
    report.label(270, 220, "new bond between owners", GREEN)
    report.label(500, 220, "atoms owned by source 2", ORANGE)
    report.bullets(48, 175, [
        "Selected pieces do not overlap: each product atom has one proposed source.",
        "Together, the selected pieces explain the complete product graph.",
        "A bond joining two colors is a proposed connection between source-derived modules.",
        "The result is a recommendation for review, not a claim that reaction conditions are known.",
    ], width=720, size=10.2, leading=13)
    report.finish_page()


def _suzuki_page(report, config):
    target = "Clc1ccc(-c2ccccc2)cc1"
    sources = (
        ("bromobenzene", "Brc1ccccc1"),
        ("4-chlorophenylboronic acid", "OB(O)c1ccc(Cl)cc1"),
    )
    _results, assembly = _assembly(target, sources, config)
    selected = {candidate.source_id: candidate for candidate in assembly.candidates}
    source_images = []
    for (source_id, smiles), color in zip(sources, (RD_BLUE, RD_ORANGE)):
        candidate = selected[source_id]
        source_images.append(_rdkit_png(
            smiles,
            _source_colors(candidate, color),
            cut_bonds=candidate.boundary_bonds,
            width=570,
            height=300,
        ))
    color_by_source = {
        "bromobenzene": RD_BLUE,
        "4-chlorophenylboronic acid": RD_ORANGE,
    }
    target_png = _rdkit_png(
        target,
        _target_colors_by_source(assembly, color_by_source),
        formed_bonds=assembly.formed_bonds,
        width=900,
        height=360,
    )

    report.new_page(
        "Example 1: two building blocks cover P_target",
        "Complete 2D molecules are shown; identical colors are assigned by the computed atom mapping")
    report.image(source_images[0], 65, 330, 280, 145)
    report.image(source_images[1], 495, 325, 290, 155)
    report.pdf.setFont("Helvetica-Bold", 10.5)
    report.pdf.setFillColor(BLUE)
    report.pdf.drawCentredString(205, 485, "R1: bromobenzene")
    report.pdf.setFillColor(ORANGE)
    report.pdf.drawCentredString(640, 485, "R2: 4-chlorophenylboronic acid")
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 24)
    report.pdf.drawCentredString(421, 387, "+")
    report.arrow(285, 329, 220, 278, color=BLUE, width=2.5)
    report.arrow(557, 329, 340, 278, color=ORANGE, width=2.5)
    report.pdf.setFont("Helvetica-Bold", 11)
    report.pdf.drawCentredString(255, 265, "P_target: 4-chlorobiphenyl")
    report.image(target_png, 55, 75, 400, 185)
    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(475, 82, 325, 165, 9, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 11)
    report.pdf.drawString(492, 222, "What the mapping shows")
    report.label(492, 194, "blue atoms come from R1", BLUE)
    report.label(492, 166, "orange atoms come from R2", ORANGE)
    report.label(492, 138, "red source groups are not retained", RED)
    report.label(492, 110, "green is the new bond joining the pieces", GREEN)
    report.finish_page()


def _alternative_patterns_page(report, config):
    target = "O=C(O)Cc1ccccc1"
    source_smiles = {
        "phenylacetyl chloride": "O=C(Cl)Cc1ccccc1",
        "water": "O",
        "benzyl bromide": "BrCc1ccccc1",
        "formic acid": "O=CO",
        "benzene": "c1ccccc1",
        "acetic acid": "CC(=O)O",
    }
    detections = {
        source_id: _detect(source_id, smiles, target, config)
        for source_id, smiles in source_smiles.items()
    }
    candidates = tuple(
        candidate
        for detection in detections.values()
        for candidate in detection.candidates
    )
    search = assemble_fragment_cover(
        _graph(target),
        candidates,
        maximum_precursors=2,
        assembly_limit=1_000,
        require_attachment_bonds=False,
    )
    requested = (
        (
            "Pattern A",
            "carbonyl framework retained",
            ("phenylacetyl chloride", "water"),
        ),
        (
            "Pattern B",
            "benzyl and carboxyl modules",
            ("benzyl bromide", "formic acid"),
        ),
        (
            "Pattern C",
            "aryl and acetic-acid modules",
            ("benzene", "acetic acid"),
        ),
    )
    rows = []
    for label, description, source_ids in requested:
        matches = [
            assembly for assembly in search.assemblies
            if frozenset(assembly.precursor_ids) == frozenset(source_ids)
        ]
        if not matches:
            raise RuntimeError(f"missing alternative assembly: {source_ids}")
        assembly = min(matches, key=lambda item: (
            len(item.broken_bonds),
            len(item.formed_bonds),
            tuple(candidate.covered_target_atoms
                  for candidate in item.candidates),
        ))
        selected = {
            candidate.source_id: candidate
            for candidate in assembly.candidates
        }
        color_by_source = {
            source_ids[0]: RD_BLUE,
            source_ids[1]: RD_ORANGE,
        }
        source_pngs = []
        for source_id in source_ids:
            candidate = selected[source_id]
            source_pngs.append(_rdkit_png(
                source_smiles[source_id],
                _source_colors(candidate, color_by_source[source_id]),
                cut_bonds=candidate.boundary_bonds,
                width=500,
                height=260,
            ))
        target_png = _rdkit_png(
            target,
            _target_colors_by_source(assembly, color_by_source),
            formed_bonds=assembly.formed_bonds,
            width=700,
            height=280,
        )
        rows.append((
            label, description, source_ids, source_pngs, target_png,
        ))

    report.new_page(
        "Example 2: alternative building-block patterns",
        "One best representative is shown for each distinct structural division of phenylacetic acid")
    row_y = (365, 225, 85)
    for (label, description, source_ids, source_pngs, target_png), y in zip(
            rows, row_y):
        report.pdf.setFillColor(LIGHT)
        report.pdf.setStrokeColor(LINE)
        report.pdf.roundRect(42, y, 758, 120, 9, fill=1, stroke=1)
        report.pdf.setFillColor(INK)
        report.pdf.setFont("Helvetica-Bold", 11)
        report.pdf.drawString(56, y + 91, label)
        report.text(
            56, y + 73, description, width=105, size=8.5,
            leading=11, color=MUTED)
        report.image(source_pngs[0], 165, y + 22, 145, 83)
        report.pdf.setFillColor(BLUE)
        report.pdf.setFont("Helvetica-Bold", 8.5)
        report.pdf.drawCentredString(237, y + 104, f"R1: {source_ids[0]}")
        report.pdf.setFillColor(INK)
        report.pdf.setFont("Helvetica-Bold", 15)
        report.pdf.drawCentredString(321, y + 57, "+")
        report.image(source_pngs[1], 330, y + 22, 135, 83)
        report.pdf.setFillColor(ORANGE)
        report.pdf.setFont("Helvetica-Bold", 8.5)
        report.pdf.drawCentredString(397, y + 104, f"R2: {source_ids[1]}")
        report.arrow(473, y + 68, 525, y + 68, color=BLUE, width=2)
        report.arrow(473, y + 51, 525, y + 51, color=ORANGE, width=2)
        report.image(target_png, 535, y + 14, 245, 96)
        report.pdf.setFillColor(INK)
        report.pdf.setFont("Helvetica-Bold", 8.5)
        report.pdf.drawCentredString(657, y + 104, "P_target: phenylacetic acid")
    report.bullets(48, 62, [
        "Every row completely covers the same P_target, but divides it between R1 and R2 differently.",
        "The workflow returns one representative per pattern; chemical feasibility is reviewed separately.",
    ], width=745, size=8.3, leading=9)
    report.finish_page()


def _complex_page(report, config):
    target = "Brc1cccc(-c2ccc(N(c3ccccc3)c3ccccc3)cc2)c1"
    sources = (
        ("1,3-dibromobenzene", "Brc1cccc(Br)c1"),
        ("triarylamine BPin", "CC1(C)OB(c2ccc(N(c3ccccc3)c3ccccc3)cc2)OC1(C)C"),
    )
    _results, assembly = _assembly(target, sources, config)
    selected = {candidate.source_id: candidate for candidate in assembly.candidates}
    left = selected["1,3-dibromobenzene"]
    right = selected["triarylamine BPin"]
    left_png = _rdkit_png(
        sources[0][1], _source_colors(left, RD_BLUE),
        cut_bonds=left.boundary_bonds, width=600, height=320)
    right_png = _rdkit_png(
        sources[1][1], _source_colors(right, RD_ORANGE),
        cut_bonds=right.boundary_bonds, width=900, height=360)
    color_by_source = {
        "1,3-dibromobenzene": RD_BLUE,
        "triarylamine BPin": RD_ORANGE,
    }
    target_png = _rdkit_png(
        target, _target_colors_by_source(assembly, color_by_source),
        formed_bonds=assembly.formed_bonds, width=1000, height=390)

    report.new_page(
        "Example 4: the same workflow on a larger target",
        "R1 and R2 are mapped into a complete P_target using the same color convention")
    report.image(left_png, 35, 325, 220, 145)
    report.image(right_png, 225, 300, 330, 185)
    report.image(target_png, 570, 295, 255, 190)
    report.arrow(550, 393, 590, 393, color=BLUE, width=2.4)
    report.arrow(550, 371, 590, 371, color=ORANGE, width=2.4)
    report.pdf.setFillColor(BLUE)
    report.pdf.setFont("Helvetica-Bold", 9.5)
    report.pdf.drawCentredString(145, 475, "R1: 1,3-dibromobenzene")
    report.pdf.setFillColor(ORANGE)
    report.pdf.drawCentredString(390, 475, "R2: triarylamine BPin")
    report.pdf.setFillColor(INK)
    report.pdf.drawCentredString(697, 475, "P_target: assembled product")
    report.label(48, 270, "retained smaller aryl module", BLUE)
    report.label(275, 270, "retained triarylamine module", ORANGE)
    report.label(585, 270, "new inter-module bond", GREEN)
    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(42, 75, 758, 155, 10, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 12)
    report.pdf.drawString(58, 204, "What this example demonstrates")
    report.bullets(58, 180, [
        "The large triarylamine system remains one coherent orange contribution from R2.",
        "The smaller blue aryl module comes from R1.",
        "Unused bromine and the BPin unit are uncolored because they are not retained in P_target.",
        "Together, the blue and orange regions cover the complete target; green marks their new bond.",
    ], width=720, size=9.6, leading=12)
    report.finish_page()


def _multistep_page(report, config):
    aniline = "Nc1ccccc1"
    bromobenzene = "Brc1ccccc1"
    diphenylamine = "c1ccc(Nc2ccccc2)cc1"
    triphenylamine = "c1ccc(N(c2ccccc2)c2ccccc2)cc1"

    _step_one_results, step_one = _assembly(
        diphenylamine,
        (("aniline", aniline), ("bromobenzene", bromobenzene)),
        config,
    )
    step_one_selected = {
        candidate.source_id: candidate for candidate in step_one.candidates
    }
    step_one_colors = {
        "aniline": RD_BLUE,
        "bromobenzene": RD_ORANGE,
    }
    aniline_png = _rdkit_png(
        aniline,
        _source_colors(step_one_selected["aniline"], RD_BLUE),
        cut_bonds=step_one_selected["aniline"].boundary_bonds,
        width=500,
        height=280,
    )
    first_bromobenzene_png = _rdkit_png(
        bromobenzene,
        _source_colors(step_one_selected["bromobenzene"], RD_ORANGE),
        cut_bonds=step_one_selected["bromobenzene"].boundary_bonds,
        width=500,
        height=280,
    )
    intermediate_colors = _target_colors_by_source(
        step_one, step_one_colors)
    intermediate_png = _rdkit_png(
        diphenylamine,
        intermediate_colors,
        formed_bonds=step_one.formed_bonds,
        width=760,
        height=300,
    )

    _step_two_results, step_two = _assembly(
        triphenylamine,
        (("diphenylamine", diphenylamine),
         ("bromobenzene", bromobenzene)),
        config,
    )
    step_two_selected = {
        candidate.source_id: candidate for candidate in step_two.candidates
    }
    intermediate_candidate = step_two_selected["diphenylamine"]
    final_colors = {
        target_atom: intermediate_colors[source_atom]
        for source_atom, target_atom in intermediate_candidate.mapping
        if source_atom in intermediate_colors
    }
    final_colors.update({
        target_atom: RD_PURPLE
        for target_atom in step_two_selected[
            "bromobenzene"].covered_target_atoms
    })
    intermediate_as_source_png = _rdkit_png(
        diphenylamine,
        intermediate_colors,
        cut_bonds=intermediate_candidate.boundary_bonds,
        formed_bonds=step_one.formed_bonds,
        width=760,
        height=300,
    )
    second_bromobenzene_png = _rdkit_png(
        bromobenzene,
        _source_colors(
            step_two_selected["bromobenzene"], RD_PURPLE),
        cut_bonds=step_two_selected["bromobenzene"].boundary_bonds,
        width=500,
        height=280,
    )
    final_png = _rdkit_png(
        triphenylamine,
        final_colors,
        formed_bonds=step_two.formed_bonds,
        width=850,
        height=320,
    )

    report.new_page(
        "Harder Example: Two Steps with Multiple R",
        "The same geometric workflow is applied recursively; bromobenzene is used twice")
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 12)
    report.pdf.drawString(45, 480, "Step 1: construct the intermediate")
    report.image(aniline_png, 45, 330, 165, 125)
    report.image(first_bromobenzene_png, 235, 330, 165, 125)
    report.image(intermediate_png, 515, 320, 280, 140)
    report.pdf.setFont("Helvetica-Bold", 9)
    report.pdf.setFillColor(BLUE)
    report.pdf.drawCentredString(127, 462, "R1: aniline")
    report.pdf.setFillColor(ORANGE)
    report.pdf.drawCentredString(317, 462, "R2: bromobenzene, copy 1")
    report.pdf.setFillColor(INK)
    report.pdf.drawCentredString(655, 462, "P1: diphenylamine")
    report.pdf.setFont("Helvetica-Bold", 17)
    report.pdf.drawCentredString(220, 385, "+")
    report.arrow(415, 401, 495, 401, color=BLUE, width=2)
    report.arrow(415, 381, 495, 381, color=ORANGE, width=2)

    report.pdf.setStrokeColor(LINE)
    report.pdf.line(42, 305, 800, 305)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 12)
    report.pdf.drawString(45, 280, "Step 2: expand the intermediate to P_target")
    report.image(intermediate_as_source_png, 45, 145, 280, 125)
    report.image(second_bromobenzene_png, 350, 145, 150, 125)
    report.image(final_png, 585, 135, 210, 140)
    report.pdf.setFont("Helvetica-Bold", 9)
    report.pdf.setFillColor(INK)
    report.pdf.drawCentredString(185, 270, "R1: P1 intermediate")
    report.pdf.setFillColor(PURPLE)
    report.pdf.drawCentredString(425, 270, "R3: bromobenzene, copy 2")
    report.pdf.setFillColor(INK)
    report.pdf.drawCentredString(690, 270, "P_target: triphenylamine")
    report.pdf.setFont("Helvetica-Bold", 17)
    report.pdf.drawCentredString(337, 205, "+")
    report.arrow(510, 221, 565, 221, color=BLUE, width=2)
    report.arrow(510, 201, 565, 201, color=ORANGE, width=2)
    report.arrow(510, 181, 565, 181, color=PURPLE, width=2)

    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(42, 55, 758, 66, 9, fill=1, stroke=0)
    report.bullets(58, 98, [
        "Both stages independently achieve complete explicit-H target coverage without a branch-cap hit.",
        "The direct foundation search also recovers aniline plus two copies of bromobenzene; a route planner is needed to introduce P1 as the intermediate.",
    ], width=720, size=8.8, leading=10)
    report.finish_page()


def _large_star_page(report, config):
    target = (
        "c1ccc(-c2ccc(-c3cc(-c4ccc(-c5ccccc5)cc4)cc"
        "(-c4ccc(-c5ccccc5)cc4)c3)cc2)cc1"
    )
    core_smiles = "Brc1cc(Br)cc(Br)c1"
    arm_smiles = "CN1CC(=O)OB(c2ccc(-c3ccccc3)cc2)OC(=O)C1"
    target_graph = _graph(target)
    core_result = _detect(
        "1,3,5-tribromobenzene", core_smiles, target, config)
    arm_result = _detect(
        "4-biphenylboronic acid MIDA ester",
        arm_smiles,
        target,
        config,
    )
    search = assemble_fragment_cover(
        target_graph,
        core_result.candidates + arm_result.candidates,
        maximum_precursors=4,
        assembly_limit=1_000,
        allow_repeated_precursors=True,
        require_attachment_bonds=False,
    )
    if not search.assemblies:
        raise RuntimeError("large repeated-arm example did not assemble")
    assembly = search.assemblies[0]
    core = next(
        candidate for candidate in assembly.candidates
        if candidate.source_id == "1,3,5-tribromobenzene")
    arms = sorted(
        (candidate for candidate in assembly.candidates
         if candidate.source_id ==
         "4-biphenylboronic acid MIDA ester"),
        key=lambda candidate: candidate.covered_target_atoms,
    )
    arm_colors = (RD_ORANGE, RD_PURPLE, RD_CYAN)
    core_png = _rdkit_png(
        core_smiles,
        _source_colors(core, RD_BLUE),
        cut_bonds=core.boundary_bonds,
        width=500,
        height=300,
    )
    arm_pngs = [
        _rdkit_png(
            arm_smiles,
            _source_colors(candidate, color),
            cut_bonds=candidate.boundary_bonds,
            width=650,
            height=300,
        )
        for candidate, color in zip(arms, arm_colors)
    ]
    target_colors = {
        atom: RD_BLUE for atom in core.covered_target_atoms
    }
    for candidate, color in zip(arms, arm_colors):
        target_colors.update({
            atom: color for atom in candidate.covered_target_atoms
        })
    target_png = _rdkit_png(
        target,
        target_colors,
        formed_bonds=assembly.formed_bonds,
        width=1_050,
        height=450,
    )

    report.new_page(
        "Really Difficult Example: Seven-Ring Star Target",
        "One central building block plus three repeated biphenyl-arm building blocks")
    centers = (100, 285, 470, 655)
    report.image(core_png, 25, 335, 150, 120)
    for index, png in enumerate(arm_pngs):
        report.image(png, 185 + index * 185, 330, 180, 130)
    labels = (
        ("R1: tribromobenzene core", BLUE),
        ("R2: biphenyl arm, copy 1", ORANGE),
        ("R3: biphenyl arm, copy 2", PURPLE),
        ("R4: biphenyl arm, copy 3", CYAN),
    )
    for center, (label, color) in zip(centers, labels):
        report.pdf.setFillColor(color)
        report.pdf.setFont("Helvetica-Bold", 8.5)
        report.pdf.drawCentredString(center, 470, label)
    endpoints = (155, 245, 335, 425)
    for center, endpoint, color in zip(
            centers, endpoints, (BLUE, ORANGE, PURPLE, CYAN)):
        report.arrow(
            center, 325, endpoint, 300,
            color=color, width=2)

    report.image(target_png, 45, 90, 460, 200)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 10.5)
    report.pdf.drawCentredString(
        275, 292, "P_target: complete seven-ring product")
    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(525, 92, 275, 195, 9, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 11)
    report.pdf.drawString(542, 262, "Computed result")
    report.bullets(542, 238, [
        "42 heavy atoms and 72 explicit atoms are fully covered.",
        "Four R copies are selected: one core and three identical arm structures.",
        "Three inter-module bonds are formed; all searches finish below branch cap 100.",
        "This gives the foundations, but does not choose the order or chemistry of the three coupling stages.",
    ], width=240, size=8.8, leading=11)
    report.finish_page()


def _extreme_star_page(report, result_path):
    payload = json.loads(Path(result_path).read_text())
    target_smiles = payload["target"]["smiles"]
    assembly = payload["assembly"]["assemblies"][0]
    source_smiles = {
        detection["source_id"]: detection["representation"]
        for detection in payload["detections"]
    }
    detection_candidates = {
        (
            detection["source_id"],
            tuple(tuple(pair) for pair in candidate["mapping"]),
        ): candidate
        for detection in payload["detections"]
        for candidate in detection["candidates"]
    }
    core = next(
        candidate for candidate in assembly["candidate_mappings"]
        if candidate["source_id"] == "1,3,5-tribromobenzene")
    arms = sorted(
        (candidate for candidate in assembly["candidate_mappings"]
         if candidate["source_id"] == "large OLED boronic ester"),
        key=lambda candidate: candidate["covered_target_atoms"],
    )
    colors = (RD_ORANGE, RD_PURPLE, RD_CYAN)

    def source_png(candidate, color, width):
        mapping = tuple(tuple(pair) for pair in candidate["mapping"])
        detection_candidate = detection_candidates[
            (candidate["source_id"], mapping)]
        return _rdkit_png(
            source_smiles[candidate["source_id"]],
            {source: color for source, _target in mapping},
            cut_bonds=tuple(
                tuple(bond)
                for bond in detection_candidate["boundary_bonds"]),
            width=width,
            height=320,
        )

    core_png = source_png(core, RD_BLUE, 500)
    arm_pngs = [
        source_png(candidate, color, 800)
        for candidate, color in zip(arms, colors)
    ]
    target_colors = {
        atom: RD_BLUE for atom in core["covered_target_atoms"]
    }
    for candidate, color in zip(arms, colors):
        target_colors.update({
            atom: color for atom in candidate["covered_target_atoms"]
        })
    target_png = _rdkit_png(
        target_smiles,
        target_colors,
        formed_bonds=tuple(
            tuple(bond) for bond in assembly["formed_bonds"]),
        width=1_300,
        height=600,
    )

    report.new_page(
        "Extreme Size Test: 183 Explicit Atoms",
        "A 108-heavy-atom target assembled from one core and three repeated OLED building blocks")
    centers = (92, 280, 468, 656)
    report.image(core_png, 22, 335, 140, 125)
    for index, png in enumerate(arm_pngs):
        report.image(png, 170 + index * 188, 325, 190, 140)
    labels = (
        ("R1: central core", BLUE),
        ("R2: large arm, copy 1", ORANGE),
        ("R3: large arm, copy 2", PURPLE),
        ("R4: large arm, copy 3", CYAN),
    )
    for center, (label, color) in zip(centers, labels):
        report.pdf.setFillColor(color)
        report.pdf.setFont("Helvetica-Bold", 8.5)
        report.pdf.drawCentredString(center, 470, label)
    endpoints = (155, 245, 335, 425)
    for center, endpoint, color in zip(
            centers, endpoints, (BLUE, ORANGE, PURPLE, CYAN)):
        report.arrow(center, 320, endpoint, 294, color=color, width=2)

    report.image(target_png, 35, 72, 500, 220)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 10.5)
    report.pdf.drawCentredString(
        285, 296, "P_target: complete 108-heavy-atom product")
    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(555, 78, 245, 215, 9, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 11)
    report.pdf.drawString(572, 268, "Saved computed result")
    report.bullets(572, 244, [
        "Complete coverage: 9 + 58 + 58 + 58 explicit atoms.",
        "Three product bonds join the four source-derived modules.",
        "Detection takes 147.6 seconds; final assembly takes 0.004 seconds.",
        "Maximum branch count is four, with no cap hit: source scanning, not branching, is the bottleneck.",
    ], width=210, size=8.7, leading=11)
    report.finish_page()


def _next_questions_page(report):
    report.new_page(
        "Open Questions: From Geometry to Retrosynthesis",
        "The geometric candidates are a foundation; two capabilities are still required")

    report.pdf.setFillColor(LIGHT)
    report.pdf.setStrokeColor(GREEN)
    report.pdf.roundRect(42, 410, 758, 72, 10, fill=1, stroke=1)
    report.pdf.setFillColor(GREEN)
    report.pdf.setFont("Helvetica-Bold", 13)
    report.pdf.drawString(58, 456, "CURRENT OUTPUT")
    report.text(
        58, 434,
        "Structurally compatible building-block sets that cover P_target from a pure geometrical perspective.",
        width=720, size=10.5, leading=14)

    report.pdf.setFillColor(ALERT_BG)
    report.pdf.setStrokeColor(ALERT)
    report.pdf.roundRect(42, 205, 360, 175, 10, fill=1, stroke=1)
    report.pdf.setFillColor(ALERT)
    report.pdf.roundRect(42, 368, 360, 12, 5, fill=1, stroke=0)
    report.pdf.setFont("Helvetica-Bold", 14)
    report.pdf.drawString(58, 344, "1. Chemical reactivity rules")
    report.pdf.setFont("Helvetica-Bold", 11)
    report.pdf.drawString(58, 320, "Can the proposed transformation actually occur?")
    report.bullets(60, 294, [
        "Which bonds can be formed or broken by known chemistry?",
        "Are functional groups, reagents, catalysts, and conditions compatible?",
        "Will chemo-, regio-, and stereoselectivity be acceptable?",
    ], width=320, size=9.3, leading=12)

    report.pdf.setFillColor(HexColor("#FFF8DC"))
    report.pdf.setStrokeColor(GOLD)
    report.pdf.roundRect(440, 205, 360, 175, 10, fill=1, stroke=1)
    report.pdf.setFillColor(GOLD)
    report.pdf.roundRect(440, 368, 360, 12, 5, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 14)
    report.pdf.drawString(456, 344, "2. Multi-step route construction")
    report.pdf.setFont("Helvetica-Bold", 11)
    report.pdf.drawString(456, 320, "How do the foundations become a complete route?")
    report.bullets(458, 294, [
        "Use accepted building blocks or intermediates as new targets.",
        "Connect validated single steps into a retrosynthesis path.",
        "Control branching and rank complete routes by chemical feasibility.",
    ], width=320, size=9.3, leading=12)

    y = 85
    report.box(48, y, 160, 70, "Geometric candidates",
               "Current workflow", color=PURPLE)
    report.box(250, y, 160, 70, "Chemical validation",
               "Apply reactivity rules", color=ALERT)
    report.box(452, y, 160, 70, "Multi-step expansion",
               "Construct route paths", color=GOLD)
    report.box(654, y, 140, 70, "Proposed routes",
               "Rank for review", color=GREEN)
    report.arrow(211, y + 35, 246, y + 35, color=MUTED, width=2)
    report.arrow(413, y + 35, 448, y + 35, color=MUTED, width=2)
    report.arrow(615, y + 35, 650, y + 35, color=MUTED, width=2)
    report.finish_page()


def _inventory_page(report):
    report.new_page(
        "Inventory integration and interpretation",
        "The structure bank provides candidates; graph matching supplies atom-level evidence")
    report.box(42, 400, 170, 90, "3,059 containers",
               "Physical inventory rows retained with barcode and location metadata.", color=CYAN)
    report.box(237, 400, 170, 90, "2,769 resolved",
               "Containers with a PubChem structure that parses in RDKit.", color=BLUE)
    report.box(432, 400, 170, 90, "1,919 structures",
               "Unique stereochemical SMILES in the deduplicated search bank.", color=PURPLE)
    report.box(627, 400, 170, 90, "290 unresolved",
               "Not-found or ambiguous identifiers remain excluded and auditable.", color=ORANGE)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 14)
    report.pdf.drawString(42, 350, "Recommended use")
    report.bullets(44, 325, [
        "Compare every resolved inventory structure with the desired target.",
        "Collect compounds that contribute substantial coherent pieces.",
        "Combine complementary pieces into several distinct construction patterns.",
        "Review the leading recommendations with stereochemical and chemical knowledge.",
    ], width=750)
    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(42, 92, 758, 120, 10, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 12)
    report.pdf.drawString(58, 185, "Reading the colors")
    report.text(58, 163,
        "Blue, orange, and purple identify source ownership. Red marks a source boundary bond cut. Green marks a product bond connecting different owners. Uncolored atoms remain in a source residual or are not supplied by that candidate. Hydrogens participate explicitly in the computation, although routine C-H labels are hidden in the skeletal drawings for legibility.",
        width=718, size=10, leading=14)
    report.text(58, 112,
        "Recommendations favor a small number of source structures that retain a large fraction of their atoms with few structural changes. This ordering is explainable but is not a substitute for reaction conditions or mechanistic judgment.",
        width=718, size=9.5, leading=13, color=MUTED)
    report.finish_page()


def _references_page(report):
    report.new_page("How to interpret the result")
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 14)
    report.pdf.drawString(42, 475, "What the report supports")
    report.bullets(44, 450, [
        "A visual explanation of which molecular pieces can be conserved in the target.",
        "A transparent proposal for which inventory compounds supply each product region.",
        "A compact set of structurally distinct precursor recommendations.",
        "A starting point for a chemist to assess reactivity, conditions, selectivity, and availability.",
    ], width=745)
    report.pdf.setFont("Helvetica-Bold", 14)
    report.pdf.setFillColor(INK)
    report.pdf.drawString(42, 320, "References")
    refs = [
        "1. PubChem PUG REST: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest",
        "2. RDKit documentation: https://www.rdkit.org/docs/",
        "3. Project design specification: RETROSYNTHESIS_DESIGN.md in the coordinate_alignment repository.",
    ]
    y = 290
    for ref in refs:
        y = report.text(55, y, ref, width=730, size=10, leading=14)
        y -= 8
    report.pdf.setFillColor(LIGHT)
    report.pdf.roundRect(42, 82, 758, 92, 10, fill=1, stroke=0)
    report.pdf.setFillColor(INK)
    report.pdf.setFont("Helvetica-Bold", 11)
    report.pdf.drawString(58, 148, "Figure provenance")
    report.text(58, 128,
        "Regenerate this PDF with tools/build_fragment_matching_report.py. Molecular highlights are recomputed from the current detector and assembly implementation each time.",
        width=718, size=9.7, leading=13)
    report.finish_page()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--large-star-result",
        default="reports/large_star_183_result.json",
    )
    args = parser.parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    config = FragmentDetectionConfig(
        minimum_fragment_size=1,
        iso_tolerance=0.5,
        branch_limit=500,
        candidate_limit=500,
    )
    report = Report(output)
    _cover_page(report)
    _suzuki_page(report, config)
    _alternative_patterns_page(report, config)
    _co2_page(report, config)
    _complex_page(report, config)
    _multistep_page(report, config)
    _large_star_page(report, config)
    _extreme_star_page(report, args.large_star_result)
    _next_questions_page(report)
    report.save()
    print(output.resolve())


if __name__ == "__main__":
    main()
