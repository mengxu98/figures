from html import escape
from pathlib import Path
import xml.etree.ElementTree as ET

WIDTH, HEIGHT = 1440, 970
FONT = "Arial, Helvetica, sans-serif"
PALETTE = {
    "ink": "#221815", "panel": "#E9EAF0", "edge": "#727171",
    "border": "#A0A0A0", "white": "#FFFFFF", "gene": "#F14F0D",
    "red": "#C9242B", "blue": "#2F4E9A", "green": "#59B34E",
    "pale_blue": "#D7DFEE", "pale_green": "#DAEBCF", "pale_red": "#FBDCD7",
    "heat0": "#2F4E9A", "heat1": "#8E97CA", "heat2": "#BFC2E2",
    "heat3": "#E3E5F3", "heat4": "#F6F6FB", "heat5": "#FEF3F1",
    "heat6": "#FADAD5", "heat7": "#F0A498", "heat8": "#DD6F61",
    "heat9": "#C9242B",
}


def complement(value):
    return "#" + "".join(f"{255 - int(value[i:i + 2], 16):02X}" for i in (1, 3, 5))


def draw(theme):
    colors = {k: complement(v) if theme == "dark" else v for k, v in PALETTE.items()}
    out = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}" '
        'role="img" aria-labelledby="title description">',
        '<title id="title">inferCSN: input, sparse regression, evidence and networks</title>',
        '<desc id="description">A four-panel schematic matching the multiCSN homepage. '
        'Expression matrices, optional supplied pseudotime and gene subsets enter '
        'Greedy-L0 sparse regression. Edge-deletion BIC evidence is transformed into '
        'signed ordinal weights. Static and pseudotime-informed networks illustrate '
        'alternative inference modes. All matrices, cells and networks are schematic.</desc>',
        '<defs>',
    ]
    for key in ("ink", "edge", "red", "blue", "green"):
        out.append(f'<marker id="arrow-{key}" viewBox="0 0 12 12" refX="11" refY="6" '
                   'markerWidth="11" markerHeight="11" markerUnits="userSpaceOnUse" '
                   f'orient="auto"><path d="M 1 1 L 11 6 L 1 11 L 4 6 Z" '
                   f'fill="{colors[key]}"/></marker>')
    out.append('</defs>')

    def color(key):
        return colors.get(key, key)

    def text(x, y, label, size=24, anchor="start", fill="ink"):
        out.append(f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" '
                   f'fill="{color(fill)}" text-anchor="{anchor}">{escape(label)}</text>')

    def rect(x, y, w, h, fill="panel", radius=0, stroke=None):
        border = f' stroke="{color(stroke)}" stroke-width="2"' if stroke else ""
        out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" '
                   f'fill="{color(fill)}"{border}/>')

    def circle(x, y, radius, fill, stroke=None, width=2):
        border = f' stroke="{color(stroke)}" stroke-width="{width}"' if stroke else ""
        out.append(f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{color(fill)}"{border}/>')

    def line(x1, y1, x2, y2, fill="ink", width=2, arrow=False, dash=False):
        marker = f' marker-end="url(#arrow-{fill})"' if arrow else ""
        pattern = ' stroke-dasharray="7 6"' if dash else ""
        out.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
                   f'stroke="{color(fill)}" stroke-width="{width}"{marker}{pattern}/>')

    def path(d, fill="ink", width=2, arrow=False):
        marker = f' marker-end="url(#arrow-{fill})"' if arrow else ""
        out.append(f'<path d="{d}" fill="none" stroke="{color(fill)}" '
                   f'stroke-width="{width}" stroke-linejoin="round"{marker}/>')

    def start(name):
        out.append(f'<g id="{name}">')

    def end():
        out.append('</g>')

    def matrix(x, y, rows, cols, cw, ch, variant=0):
        for row in range(rows):
            for col in range(cols):
                distance = abs(col - (row + variant) % cols)
                value = (row * 2 + col * 3 + variant) % 7
                value = min(9, value + 3) if distance <= 2 else value % 6
                rect(x + cw * col, y + ch * row, cw, ch, f"heat{value}")

    def regulator(x, y, side=25):
        rect(x - side / 2, y - side / 2, side, side, "ink")

    def target(x, y, radius=18):
        circle(x, y, radius, "gene")

    def network(top, lagged=False):
        xs, targets = [1174, 1270, 1366], [1220, 1320]
        ry, ty = top + 88, top + 186
        edges = [(0, 0, False), (1, 1, False), (2, 0, True), (2, 1, False)]
        if lagged:
            edges = [(0, 1, False), (1, 0, True), (1, 1, False), (2, 1, True)]
        for a, b, negative in edges:
            line(xs[a], ry + 15, targets[b], ty - 21, "edge", 2.1, True, negative)
        for x in xs:
            regulator(x, ry, 27)
        for x in targets:
            target(x, ty)

    start("a-expression-input")
    text(10, 32, "a", 39)
    rect(32, 38, 312, 380, radius=24)
    for x, body, nucleus in [(99, "pale_blue", "blue"),
                             (186, "pale_green", "green"),
                             (272, "pale_red", "red")]:
        circle(x, 89, 26, body, nucleus, 1.4)
        circle(x - 4, 85, 9, nucleus)
        path(f'M {x + 3} 102 q 5 -8 10 -1 q 4 5 8 -3', nucleus, 1.7)
    text(189, 143, "Single-cell expression", 25, "middle")
    line(189, 155, 189, 184, arrow=True)
    text(207, 212, "Genes", 23, "middle")
    matrix(95, 229, 7, 10, 21, 17, 1)
    text(61, 293, "n", 24, "middle")
    line(61, 232, 61, 267, width=1.5)
    line(61, 304, 61, 347, width=1.5, arrow=True)
    line(44, 365, 332, 365, dash=True)
    text(188, 403, "Cells × genes", 25, "middle")
    end()

    start("a-supplied-pseudotime")
    rect(32, 449, 312, 245, radius=24)
    text(188, 488, "Supplied pseudotime", 26, "middle")
    path("M 66 565 C 111 565 141 548 169 554 C 199 525 232 517 304 525", "blue", 2.5)
    path("M 169 554 C 202 569 230 594 305 593", "green", 2.5)
    for x, y in [(69,565),(100,562),(133,554),(168,554),(199,533),(236,522),(274,521),(301,524)]:
        circle(x, y, 6.5, "blue")
    for x, y in [(198,570),(231,587),(268,592),(302,593)]:
        circle(x, y, 6.5, "green")
    line(82, 627, 294, 627, width=1.5, arrow=True)
    text(82, 649, "t", 22, "middle")
    text(287, 649, "t + lag", 22, "middle")
    text(188, 680, "Optional vector / branch matrix", 19, "middle")
    end()

    start("a-gene-subsets")
    rect(32, 739, 312, 177, radius=24)
    text(188, 778, "Gene subsets (optional)", 25, "middle")
    text(50, 823, "Regulators", 23)
    text(50, 882, "Targets", 23)
    for x in [222, 263, 304]:
        regulator(x, 815, 20)
        target(x, 874, 11)
    end()

    start("input-connections")
    for y in [265, 568, 827]:
        line(344, y, 368, y)
    line(368, 265, 368, 827)
    line(368, 380, 415, 380, arrow=True)
    end()

    start("b-sparse-regression")
    text(376, 32, "b", 39)
    rect(404, 38, 660, 140)
    text(425, 84, "Greedy-L0 sparse-regression engine", 30)
    text(734, 139, "BIC(S) = n log[max{RSS(S) / n, 10⁻¹²}] + |S| log n", 24, "middle")
    rect(424, 227, 620, 320, radius=25)
    text(734, 267, "One sparse model per target gene", 29, "middle")
    text(480, 305, "Target", 22, "middle")
    text(711, 305, "Candidate regulators", 23, "middle")
    text(936, 305, "β", 28, "middle")
    matrix(466, 325, 8, 1, 27, 20, 3)
    text(522, 414, "≈", 37, "middle")
    matrix(556, 325, 8, 10, 29, 20, 0)
    text(881, 414, "·", 38, "middle")
    for i in range(10):
        fill = "red" if i in (1, 7) else "blue" if i == 4 else "white"
        rect(924, 325 + i * 16, 27, 16, fill)
    for d in ["M 459 325 H 452 V 485 H 459", "M 500 325 H 507 V 485 H 500",
              "M 548 325 H 541 V 485 H 548", "M 853 325 H 860 V 485 H 853",
              "M 916 325 H 909 V 485 H 916", "M 958 325 H 965 V 485 H 958"]:
        path(d, width=1.5)
    text(734, 528, "BIC-guided add / drop / swap and refit", 24, "middle")
    line(734, 554, 734, 601, arrow=True, width=2.5)
    end()

    start("c-edge-evidence")
    text(376, 600, "c", 39)
    rect(404, 614, 660, 302)
    text(734, 655, "Deletion evidence and signed ranking", 29, "middle")
    text(510, 694, "Fitted model", 23, "middle")
    text(893, 694, "Edge removed", 23, "middle")
    for center, remove in [(510, False), (893, True)]:
        for i, x in enumerate([center - 53, center, center + 53]):
            removed = remove and i == 2
            line(x, 735, center, 790, "edge", 2, not removed, removed)
            regulator(x, 722, 23)
        target(center, 810, 17)
    line(600, 770, 789, 770, arrow=True)
    text(695, 751, "Remove + refit", 23, "middle")
    line(907, 757, 925, 775, "red", 3)
    line(907, 775, 925, 757, "red", 3)
    text(734, 861, "ΔBIC = BIC(edge removed) − BIC(fitted)", 25, "middle")
    text(734, 899, "Signed weight = sign(β) × ordinal evidence rank", 23, "middle")
    line(1064, 781, 1113, 781, arrow=True)
    end()

    start("d-network-output")
    text(1096, 32, "d", 39)
    rect(1120, 38, 300, 878, "none", stroke="border")
    text(1270, 80, "Network outputs", 29, "middle")
    rect(1132, 100, 276, 263)
    text(1270, 137, "Static network", 26, "middle")
    line(1132, 152, 1408, 152, dash=True)
    text(1270, 181, "Regulators", 21, "middle")
    network(122)
    text(1270, 352, "Target genes", 22, "middle")
    rect(1132, 382, 276, 263)
    text(1270, 419, "Pseudotime-informed", 23, "middle")
    line(1132, 434, 1408, 434, dash=True)
    text(1270, 464, "Regulators (t)", 21, "middle")
    network(405, lagged=True)
    text(1270, 635, "Targets (t + lag)", 22, "middle")
    rect(1132, 664, 276, 240)
    text(1270, 702, "Signed edge weights", 24, "middle")
    line(1132, 718, 1408, 718, dash=True)
    for x, label in [(1148, "regulator"), (1256, "target"), (1342, "weight")]:
        text(x, 750, label, 19)
    for y, r, g, w in [(793,"r₁","g₁","+w₁"),(836,"r₂","g₂","−w₂"),(879,"r₃","g₁","+w₃")]:
        text(1187, y, r, 25, "middle")
        text(1279, y, g, 25, "middle")
        text(1370, y, w, 25, "middle")
    end()

    start("legend")
    line(955, 945, 998, 945, "edge", 2, True)
    text(1008, 953, "Positive", 21)
    line(1200, 945, 1243, 945, "edge", 2, True, True)
    text(1253, 953, "Negative", 21)
    end()
    out.append('</svg>')
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    destination = Path(__file__).resolve().parent
    for theme, filename in [("light", "inferCSN.svg"), ("dark", "inferCSN-dark.svg")]:
        svg = draw(theme)
        root = ET.fromstring(svg)
        assert not root.findall('.//{http://www.w3.org/2000/svg}image')
        assert root.findall('.//{http://www.w3.org/2000/svg}text')
        (destination / filename).write_text(svg, encoding="utf-8")
        print(filename)
