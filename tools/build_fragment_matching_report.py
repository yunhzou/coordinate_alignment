#!/usr/bin/env python3
"""Build a graph-centered PDF explaining fragment detection and assembly."""
from __future__ import annotations

import argparse
import io
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
        "Finding Building Blocks for a Target Molecule",
        "A visual explanation of the single-step precursor recommendation workflow")
    report.pdf.setFillColor(HexColor("#10243B"))
    report.pdf.roundRect(38, 90, 766, 380, 18, fill=1, stroke=0)
    report.pdf.setFillColor(white)
    report.pdf.setFont("Helvetica-Bold", 16)
    report.pdf.drawString(68, 430, "Objective")
    report.pdf.setFont("Helvetica", 12)
    report.pdf.drawString(
        68, 405,
        "Given P_target and a bank of compounds, find small sets of building blocks whose")
    report.pdf.drawString(
        68, 386,
        "conserved molecular pieces can be combined to account for the complete target.")
    report.pdf.setFont("Helvetica-Bold", 16)
    report.pdf.drawString(68, 346, "Workflow")
    y = 205
    report.box(60, y, 155, 105, "1. Search the bank",
               "Compare every building block with P_target.", color=CYAN)
    report.box(250, y, 155, 105, "2. Find useful pieces",
               "Retain coherent fragments that occur in the target.", color=BLUE)
    report.box(440, y, 155, 105, "3. Assemble coverage",
               "Combine complementary pieces until P_target is covered.", color=PURPLE)
    report.box(630, y, 155, 105, "4. Show alternatives",
               "Return the best candidate from each distinct pattern.", color=ORANGE)
    report.arrow(218, y + 52, 246, y + 52, color=white, width=2)
    report.arrow(408, y + 52, 436, y + 52, color=white, width=2)
    report.arrow(598, y + 52, 626, y + 52, color=white, width=2)
    report.pdf.setFillColor(white)
    report.pdf.setFont("Helvetica", 10)
    report.pdf.drawString(68, 150, "The next slides show the result directly as R1 + R2 to P_target mappings.")
    report.pdf.drawString(68, 130, "Matching colors come from computed atom correspondences, not manual annotation.")
    report.pdf.drawString(68, 110, "The output is a structural recommendation; chemical feasibility still requires review.")
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
    report.save()
    print(output.resolve())


if __name__ == "__main__":
    main()
