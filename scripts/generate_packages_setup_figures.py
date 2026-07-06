from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
IMG_DIR = ROOT / "img"

RESOURCE_PNG = IMG_DIR / "packages_setup_resources_cs4.png"
WORKFLOW_PNG = IMG_DIR / "packages_setup_workflow_map_cs4.png"

RESOURCE_W, RESOURCE_H = 2000, 960
WORKFLOW_W, WORKFLOW_H = 2000, 2860

BG = "#F7FAFD"
WHITE = "#FFFFFF"
NAVY = "#08254D"
INK = "#17202A"
MUTED = "#4B6178"
LINE = "#CBD7E2"
BLUE = "#1D4ED8"
TEAL = "#0F766E"
GREEN = "#059669"
CORAL = "#D55E00"
PURPLE = "#7E22CE"
GOLD = "#B7791F"
ROSE = "#BE123C"
SLATE = "#475569"

PALE_BLUE = "#EAF2FB"
PALE_TEAL = "#E8F7F2"
PALE_GREEN = "#ECFDF5"
PALE_CORAL = "#FFF0E8"
PALE_PURPLE = "#F4ECFA"
PALE_GOLD = "#FFF7DF"
PALE_SLATE = "#F1F5F9"

CORE_COLOR = TEAL
REBUILD_COLOR = BLUE
SWEEP_COLOR = PURPLE

FONT_REG_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/Library/Fonts/Arial.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
]
FONT_BOLD_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/Library/Fonts/Arial Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
]


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = FONT_BOLD_CANDIDATES if bold else FONT_REG_CANDIDATES
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size)
    return ImageFont.load_default(size=size)


H1 = font(62, True)
H2 = font(42, True)
H3 = font(30, True)
BODY = font(25)
BODY_BOLD = font(25, True)
SMALL = font(21)
SMALL_BOLD = font(21, True)
TINY = font(17)
TINY_BOLD = font(17, True)


def rgb(hex_value: str) -> tuple[int, int, int]:
    hex_value = hex_value.lstrip("#")
    return tuple(int(hex_value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(hex_value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*rgb(hex_value), alpha)


def tint(hex_value: str, amount: float = 0.9) -> str:
    r, g, b = rgb(hex_value)
    rr = round(r * (1 - amount) + 255 * amount)
    gg = round(g * (1 - amount) + 255 * amount)
    bb = round(b * (1 - amount) + 255 * amount)
    return f"#{rr:02x}{gg:02x}{bb:02x}"


def text_size(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def wrap(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont, max_width: int) -> list[str]:
    lines: list[str] = []
    for para in text.split("\n"):
        words = para.split()
        if not words:
            lines.append("")
            continue
        line = words[0]
        for word in words[1:]:
            trial = f"{line} {word}"
            if text_size(draw, trial, fnt)[0] <= max_width:
                line = trial
            else:
                lines.append(line)
                line = word
        lines.append(line)
    return lines


class Audit:
    def __init__(self) -> None:
        self.items: list[dict[str, object]] = []
        self.failures: list[dict[str, object]] = []

    def add(self, name: str, box: tuple[int, int, int, int], parent: tuple[int, int, int, int] | None = None) -> None:
        self.items.append({"name": name, "box": list(box)})
        if parent:
            x0, y0, x1, y1 = box
            px0, py0, px1, py1 = parent
            if x0 < px0 or y0 < py0 or x1 > px1 or y1 > py1:
                self.failures.append({"type": "bounds", "name": name, "box": list(box), "parent": list(parent)})

    @staticmethod
    def overlap(a: tuple[int, int, int, int], b: tuple[int, int, int, int], pad: int = 0) -> bool:
        ax0, ay0, ax1, ay1 = a
        bx0, by0, bx1, by1 = b
        return not (ax1 + pad <= bx0 or bx1 + pad <= ax0 or ay1 + pad <= by0 or by1 + pad <= ay0)

    def no_overlap(self, group: str, boxes: list[tuple[str, tuple[int, int, int, int]]], pad: int = 0) -> None:
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                if self.overlap(boxes[i][1], boxes[j][1], pad=pad):
                    self.failures.append(
                        {
                            "type": "overlap",
                            "group": group,
                            "a": boxes[i][0],
                            "a_box": list(boxes[i][1]),
                            "b": boxes[j][0],
                            "b_box": list(boxes[j][1]),
                        }
                    )


def draw_text_block(
    draw: ImageDraw.ImageDraw,
    audit: Audit,
    name: str,
    x: int,
    y: int,
    text: str,
    fnt: ImageFont.FreeTypeFont,
    max_width: int,
    fill: str,
    parent: tuple[int, int, int, int] | None = None,
    line_gap: int = 7,
) -> tuple[int, tuple[int, int, int, int]]:
    start = y
    max_x = x
    for idx, line in enumerate(wrap(draw, text, fnt, max_width)):
        w, h = text_size(draw, line or "Ag", fnt)
        draw.text((x, y), line, font=fnt, fill=fill)
        audit.add(f"{name}-{idx}", (x, y, x + w, y + h), parent)
        max_x = max(max_x, x + w)
        y += h + line_gap
    full = (x, start, max_x, max(start, y - line_gap))
    audit.add(name, full, parent)
    return y, full


def shadow_rect(
    image: Image.Image,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str,
    outline: str,
    width: int = 2,
    shadow_alpha: int = 16,
    shadow_blur: int = 18,
    shadow_offset: tuple[int, int] = (0, 9),
) -> None:
    layer = Image.new("RGBA", image.size, (0, 0, 0, 0))
    ld = ImageDraw.Draw(layer)
    x0, y0, x1, y1 = box
    ox, oy = shadow_offset
    ld.rounded_rectangle((x0 + ox, y0 + oy, x1 + ox, y1 + oy), radius=radius, fill=(15, 23, 42, shadow_alpha))
    layer = layer.filter(ImageFilter.GaussianBlur(shadow_blur))
    image.alpha_composite(layer)
    draw = ImageDraw.Draw(image, "RGBA")
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, text: str, color: str, fnt: ImageFont.FreeTypeFont = TINY_BOLD) -> tuple[int, int, int, int]:
    w, h = text_size(draw, text, fnt)
    box = (x, y, x + w + 28, y + h + 16)
    draw.rounded_rectangle(box, radius=(h + 16) // 2, fill=tint(color, 0.86), outline=color, width=2)
    draw.text((x + 14, y + 8), text, font=fnt, fill=color)
    return box


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 6, dashed: bool = False) -> None:
    x0, y0 = start
    x1, y1 = end
    if dashed:
        segments = 10
        for i in range(segments):
            if i % 2 == 0:
                sx = x0 + (x1 - x0) * i / segments
                sy = y0 + (y1 - y0) * i / segments
                ex = x0 + (x1 - x0) * (i + 1) / segments
                ey = y0 + (y1 - y0) * (i + 1) / segments
                draw.line((sx, sy, ex, ey), fill=color, width=width)
    else:
        draw.line((x0, y0, x1, y1), fill=color, width=width)
    if x1 >= x0:
        head = [(x1, y1), (x1 - 24, y1 - 13), (x1 - 24, y1 + 13)]
    else:
        head = [(x1, y1), (x1 + 24, y1 - 13), (x1 + 24, y1 + 13)]
    draw.polygon(head, fill=color)


def draw_cube(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    pts = [(x + 48, y + 8), (x + 94, y + 32), (x + 94, y + 84), (x + 48, y + 110), (x + 2, y + 84), (x + 2, y + 32)]
    draw.line(pts + [pts[0]], fill=color, width=5, joint="curve")
    draw.line((x + 48, y + 8, x + 48, y + 58), fill=color, width=4)
    draw.line((x + 2, y + 32, x + 48, y + 58, x + 94, y + 32), fill=color, width=4, joint="curve")
    draw.text((x + 26, y + 62), "R", font=font(24, True), fill=color)
    draw.text((x + 55, y + 62), "Py", font=font(20, True), fill=color)


def draw_repo(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    draw.rounded_rectangle((x + 6, y + 10, x + 106, y + 104), radius=12, outline=color, width=5)
    for yy in [36, 58, 80]:
        draw.line((x + 26, y + yy, x + 88, y + yy), fill=color, width=5)
    draw.ellipse((x + 18, y + 28, x + 36, y + 46), fill=color)
    draw.line((x + 36, y + 37, x + 52, y + 37, x + 52, y + 75, x + 66, y + 75), fill=color, width=4)
    draw.ellipse((x + 64, y + 66, x + 82, y + 84), fill=color)


def draw_archive(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    top = (x + 10, y + 8, x + 108, y + 38)
    bottom = (x + 10, y + 78, x + 108, y + 108)
    draw.ellipse(top, fill="#EAF2FB", outline=color, width=5)
    draw.line((x + 10, y + 23, x + 10, y + 93), fill=color, width=5)
    draw.line((x + 108, y + 23, x + 108, y + 93), fill=color, width=5)
    draw.arc((x + 10, y + 42, x + 108, y + 84), start=0, end=180, fill="#9DB7D6", width=4)
    draw.ellipse(bottom, fill="#F3F7FC", outline=color, width=5)
    draw.rounded_rectangle((x + 32, y + 42, x + 90, y + 92), radius=10, fill="white", outline="#9DB7D6", width=3)
    draw.line((x + 43, y + 54, x + 79, y + 54, x + 43, y + 82, x + 81, y + 82), fill=color, width=8, joint="curve")


def draw_cluster(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    for row in range(3):
        yy = y + row * 32
        draw.rounded_rectangle((x, yy, x + 118, yy + 22), radius=7, outline=color, width=4)
        for i in range(3):
            draw.ellipse((x + 12 + i * 18, yy + 7, x + 20 + i * 18, yy + 15), fill=color)
    draw.line((x + 60, y + 97, x + 60, y + 122), fill=color, width=4)
    draw.line((x + 36, y + 122, x + 84, y + 122), fill=color, width=4)


def draw_cells(draw: ImageDraw.ImageDraw, x: int, y: int, scale: float = 1.0) -> None:
    colors = [TEAL, BLUE, "#F2C14E", "#7B61B3", "#E76F51", "#2A9D8F"]
    coords = [(0, 20), (50, 0), (103, 24), (152, 6), (44, 70), (112, 82)]
    for idx, (dx, dy) in enumerate(coords):
        cx = x + int(dx * scale)
        cy = y + int(dy * scale)
        r = int(26 * scale)
        draw.ellipse((cx, cy, cx + 2 * r, cy + 2 * r), fill=tint(colors[idx], 0.25), outline=colors[idx], width=max(2, int(4 * scale)))
        draw.ellipse((cx + int(11 * scale), cy + int(11 * scale), cx + int(33 * scale), cy + int(33 * scale)), fill=tint(colors[idx], 0.55), outline=colors[idx], width=max(1, int(2 * scale)))


def draw_matrix_icon(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    cell = 15
    for r in range(5):
        for c in range(5):
            fill = tint(color, 0.25 + 0.10 * ((r + c) % 3))
            draw.rectangle((x + c * (cell + 4), y + r * (cell + 4), x + c * (cell + 4) + cell, y + r * (cell + 4) + cell), fill=fill)
    draw.line((x + 112, y, x + 112, y + 95), fill=color, width=4)
    for i, col in enumerate([TEAL, CORAL, PURPLE, GOLD]):
        draw.rounded_rectangle((x + 132 + i * 24, y + 18, x + 146 + i * 24, y + 76), radius=5, fill=col)


def draw_file_stack(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    for offset in [18, 9, 0]:
        draw.rounded_rectangle((x + offset, y + offset, x + 104 + offset, y + 78 + offset), radius=9, fill=WHITE, outline=color, width=4)
        draw.line((x + 22 + offset, y + 28 + offset, x + 80 + offset, y + 28 + offset), fill=color, width=4)
        draw.line((x + 22 + offset, y + 48 + offset, x + 70 + offset, y + 48 + offset), fill=color, width=4)


def draw_doc_icon(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    draw.rounded_rectangle((x + 8, y + 0, x + 100, y + 110), radius=10, fill=WHITE, outline=color, width=5)
    draw.line((x + 31, y + 28, x + 80, y + 28), fill=color, width=5)
    draw.line((x + 31, y + 51, x + 80, y + 51), fill=color, width=5)
    draw.line((x + 31, y + 74, x + 67, y + 74), fill=color, width=5)
    draw.rectangle((x + 22, y + 26, x + 25, y + 29), fill=color)


def draw_resource_card(
    image: Image.Image,
    audit: Audit,
    box: tuple[int, int, int, int],
    title: str,
    subtitle: str,
    body: str,
    color: str,
    fill: str,
    icon: str,
    badges: list[str],
    dashed: bool = False,
) -> None:
    draw = ImageDraw.Draw(image, "RGBA")
    outline = color
    if dashed:
        shadow_rect(image, box, 26, fill, "#D7CBE9", width=2, shadow_alpha=10)
        draw = ImageDraw.Draw(image, "RGBA")
        x0, y0, x1, y1 = box
        for i in range(x0 + 18, x1 - 18, 34):
            draw.line((i, y0, min(i + 18, x1 - 18), y0), fill=outline, width=4)
            draw.line((i, y1, min(i + 18, x1 - 18), y1), fill=outline, width=4)
        for i in range(y0 + 18, y1 - 18, 34):
            draw.line((x0, i, x0, min(i + 18, y1 - 18)), fill=outline, width=4)
            draw.line((x1, i, x1, min(i + 18, y1 - 18)), fill=outline, width=4)
    else:
        shadow_rect(image, box, 26, fill, outline, width=3, shadow_alpha=14)
        draw = ImageDraw.Draw(image, "RGBA")

    x0, y0, x1, y1 = box
    inner = (x0 + 32, y0 + 32, x1 - 32, y1 - 32)
    if icon == "docker":
        draw_cube(draw, x0 + 34, y0 + 48, color)
    elif icon == "repo":
        draw_repo(draw, x0 + 38, y0 + 48, color)
    elif icon == "archive":
        draw_archive(draw, x0 + 35, y0 + 50, color)
    elif icon == "cluster":
        draw_cluster(draw, x0 + 35, y0 + 50, color)

    tx = x0 + 170
    cursor, title_box = draw_text_block(draw, audit, f"{title}-title", tx, y0 + 40, title, H3, x1 - tx - 28, color, inner, 6)
    cursor += 8
    cursor, sub_box = draw_text_block(draw, audit, f"{title}-sub", tx, cursor, subtitle, SMALL_BOLD, x1 - tx - 28, color, inner, 5)
    cursor = max(cursor + 18, y0 + 166)
    cursor, body_box = draw_text_block(draw, audit, f"{title}-body", x0 + 34, cursor, body, BODY, x1 - x0 - 68, INK, inner, 8)
    bx = x0 + 34
    by = y1 - 42
    badge_boxes = []
    for label in badges:
        b = badge(draw, bx, by, label, color)
        badge_boxes.append((label, b))
        bx = b[2] + 12
    audit.no_overlap(title, [("title", title_box), ("subtitle", sub_box), ("body", body_box), *badge_boxes], pad=4)


@dataclass(frozen=True)
class WorkflowStep:
    n: int
    title: str
    body: str
    color: str
    fill: str
    icon: str
    badges: tuple[str, ...] = ()


CORE_STEPS = [
    WorkflowStep(1, "Open the prepared PBMC data", "Use the included teaching data with cells, genes, condition, donor, cell type, and UMAP coordinates.", TEAL, PALE_TEAL, "cells", ("R tab", "Python tab")),
    WorkflowStep(2, "Inspect model-resolution evidence", "Read included sweep summaries and diagnostic tables before interpreting latent patterns.", GOLD, PALE_GOLD, "files", ("included summaries",)),
    WorkflowStep(3, "Load selected CoGAPS outputs", "Use lightweight R-primary exports; Python tabs use parallel saved outputs where data work is repeated.", BLUE, PALE_BLUE, "matrix", ("R primary", "Python parallel")),
    WorkflowStep(4, "Check diagnostics", "Confirm that trace, runtime, and stability evidence support using the saved model for interpretation.", PURPLE, PALE_PURPLE, "trace", ("included diagnostics",)),
    WorkflowStep(5, "Summarize activity", "Compare pattern activity across condition, donor, and annotated PBMC cell type.", GREEN, PALE_GREEN, "bars", ("cell activity",)),
    WorkflowStep(6, "Interpret genes and direction", "Pair top gene weights with expression direction so weights are not mistaken for up- or down-regulation.", CORAL, PALE_CORAL, "direction", ("gene weights", "directionality")),
]

REBUILD_STEPS = [
    WorkflowStep(1, "Get source and large files", "Use GEO plus the Zenodo/local large-file archive when rebuilding inputs or full model objects.", REBUILD_COLOR, PALE_BLUE, "archive", ("optional",)),
    WorkflowStep(2, "Recreate preprocessing", "Rebuild the analysis-ready matrix using the recorded filtering and HVG-preparation choices.", REBUILD_COLOR, PALE_BLUE, "filter", ("full local",)),
    WorkflowStep(3, "Fit the preselected model", "Run the chosen R CoGAPS workflow locally, then recreate the lightweight exports used by learners.", REBUILD_COLOR, PALE_BLUE, "matrix", ("R primary",)),
]


def draw_step_icon(draw: ImageDraw.ImageDraw, x: int, y: int, icon: str, color: str) -> None:
    if icon == "cells":
        draw_cells(draw, x - 5, y + 5, 0.48)
    elif icon == "files":
        draw_file_stack(draw, x, y, color)
    elif icon == "matrix":
        draw_matrix_icon(draw, x, y + 6, color)
    elif icon == "trace":
        draw.line((x + 8, y + 88, x + 108, y + 88), fill=color, width=4)
        draw.line((x + 8, y + 88, x + 8, y + 16), fill=color, width=4)
        pts = [(x + 10, y + 78), (x + 30, y + 58), (x + 47, y + 66), (x + 72, y + 32), (x + 105, y + 38)]
        draw.line(pts, fill=color, width=5, joint="curve")
        for px, py in pts:
            draw.ellipse((px - 6, py - 6, px + 6, py + 6), fill=WHITE, outline=color, width=3)
    elif icon == "bars":
        for i, h in enumerate([38, 78, 54, 94]):
            draw.rounded_rectangle((x + 18 + i * 24, y + 105 - h, x + 32 + i * 24, y + 105), radius=5, fill=color)
        draw.line((x + 10, y + 106, x + 126, y + 106), fill=color, width=4)
    elif icon == "direction":
        draw.line((x + 34, y + 92, x + 34, y + 20), fill=CORAL, width=6)
        draw.polygon([(x + 34, y + 10), (x + 16, y + 32), (x + 52, y + 32)], fill=CORAL)
        draw.line((x + 88, y + 16, x + 88, y + 88), fill=BLUE, width=6)
        draw.polygon([(x + 88, y + 98), (x + 70, y + 76), (x + 106, y + 76)], fill=BLUE)
    elif icon == "archive":
        draw_archive(draw, x, y, color)
    elif icon == "filter":
        draw.line((x + 10, y + 16, x + 110, y + 16, x + 70, y + 58, x + 70, y + 98, x + 50, y + 108, x + 50, y + 58, x + 10, y + 16), fill=color, width=5, joint="curve")
    else:
        draw_doc_icon(draw, x, y, color)


def draw_workflow_step(
    image: Image.Image,
    audit: Audit,
    step: WorkflowStep,
    box: tuple[int, int, int, int],
    lane_color: str,
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(image, "RGBA")
    shadow_rect(image, box, 22, step.fill, step.color, width=3, shadow_alpha=12)
    draw = ImageDraw.Draw(image, "RGBA")
    x0, y0, x1, y1 = box
    inner = (x0 + 28, y0 + 28, x1 - 28, y1 - 28)
    draw.rounded_rectangle((x0, y0, x0 + 25, y1), radius=22, fill=lane_color)
    draw.rectangle((x0 + 12, y0, x0 + 34, y1), fill=lane_color)
    draw.ellipse((x0 + 45, y0 + 36, x0 + 91, y0 + 82), fill=step.color)
    n = str(step.n)
    nw, nh = text_size(draw, n, BODY_BOLD)
    draw.text((x0 + 45 + (46 - nw) / 2, y0 + 36 + (46 - nh) / 2 - 1), n, font=BODY_BOLD, fill=WHITE)
    draw_step_icon(draw, x1 - 152, y0 + 48, step.icon, step.color)
    tx = x0 + 116
    tw = x1 - tx - 180
    cursor, title_box = draw_text_block(draw, audit, f"workflow-{step.title}-title", tx, y0 + 34, step.title, H3, tw, step.color, inner, 6)
    cursor += 8
    cursor, body_box = draw_text_block(draw, audit, f"workflow-{step.title}-body", tx, cursor, step.body, SMALL, tw, INK, inner, 7)
    bx = tx
    by = y1 - 58
    badge_boxes = []
    for label in step.badges:
        b = badge(draw, bx, by, label, step.color)
        badge_boxes.append((label, b))
        bx = b[2] + 10
    audit.no_overlap(step.title, [("title", title_box), ("body", body_box), *badge_boxes], pad=2)
    return box


def generate_resource_figure() -> None:
    image = Image.new("RGBA", (RESOURCE_W, RESOURCE_H), rgb(BG) + (255,))
    audit = Audit()
    draw = ImageDraw.Draw(image, "RGBA")

    draw.text((70, 48), "What setup provides", font=H1, fill=NAVY)
    draw.text((72, 118), "Use the Docker image and GitHub-safe files for the main lesson; add large files only for full local reproduction.", font=BODY, fill=MUTED)

    main_band = (64, 190, 1292, 788)
    shadow_rect(image, main_band, 30, "#FBFEFD", "#B8DCD6", width=3, shadow_alpha=10)
    draw = ImageDraw.Draw(image, "RGBA")
    draw.text((110, 222), "Core learner path", font=H2, fill=CORE_COLOR)
    draw.text((112, 272), "enough to work through the case study", font=SMALL_BOLD, fill=MUTED)

    draw_resource_card(
        image,
        audit,
        (108, 338, 596, 690),
        "Docker image",
        "software environment",
        "R, Bioconductor CoGAPS, Python, Quarto, plotting tools, and OpenMP-enabled runtime support.",
        BLUE,
        WHITE,
        "docker",
        ["tools"],
    )
    draw_resource_card(
        image,
        audit,
        (706, 338, 1248, 690),
        "GitHub repository",
        "case-study evidence",
        "Quarto pages, scripts, manifests, teaching data, model summaries, lightweight CoGAPS outputs, and figures.",
        TEAL,
        WHITE,
        "repo",
        ["main files", "included evidence"],
    )
    arrow(draw, (606, 514), (690, 514), "#9DB0C5", width=8)
    draw.text((170, 725), "Most learners only need these two pieces.", font=SMALL_BOLD, fill=CORE_COLOR)

    draw_resource_card(
        image,
        audit,
        (1360, 214, 1934, 498),
        "Zenodo / local archive",
        "optional large files",
        "Dense inputs and full model objects support deeper inspection or local rebuilding.",
        REBUILD_COLOR,
        PALE_BLUE,
        "archive",
        ["full reproduction"],
    )
    draw_resource_card(
        image,
        audit,
        (1360, 548, 1934, 820),
        "HPC sweep provenance",
        "separate audit path",
        "Included summaries support rank review. Rerunning the sweep is separate HPC work.",
        SWEEP_COLOR,
        PALE_PURPLE,
        "cluster",
        ["not main/rebuild"],
        dashed=True,
    )
    arrow(draw, (1320, 434), (1280, 530), REBUILD_COLOR, width=5, dashed=True)
    draw.text((1284, 486), "optional", font=TINY_BOLD, fill=REBUILD_COLOR)
    draw.text((1354, 82), "Optional layers", font=H2, fill=NAVY)
    draw.text((1356, 132), "use when you want to rebuild or audit beyond the core path", font=SMALL, fill=MUTED)

    footer = (64, 884, 1934, 928)
    draw.rounded_rectangle(footer, radius=16, fill=WHITE, outline="#D8E4EF", width=2)
    draw.text((92, 895), "Setup rule:", font=SMALL_BOLD, fill=NAVY)
    draw.text((224, 895), "Docker gives the tools; GitHub gives the main lesson files; Zenodo/local files and HPC sweeps are optional extensions.", font=SMALL, fill=MUTED)

    save_with_qa(image, audit, RESOURCE_PNG, "Case Study 4 setup resources, Pillow vector draft")


def generate_workflow_figure() -> None:
    image = Image.new("RGBA", (WORKFLOW_W, WORKFLOW_H), rgb(BG) + (255,))
    audit = Audit()
    draw = ImageDraw.Draw(image, "RGBA")

    draw.text((80, 54), "Choose the setup path you need", font=H1, fill=NAVY)
    draw.text((82, 126), "The main lesson analyzes included evidence. Full local reproduction rebuilds selected files. The HPC sweep is separate.", font=BODY, fill=MUTED)

    core_panel = (70, 214, 1230, 2686)
    rebuild_panel = (1292, 214, 1932, 1638)
    sweep_panel = (1292, 1752, 1932, 2458)
    shadow_rect(image, core_panel, 28, "#FBFEFD", "#B8DCD6", width=3, shadow_alpha=12)
    shadow_rect(image, rebuild_panel, 28, "#FAFCFF", "#B8D0EA", width=3, shadow_alpha=12)
    shadow_rect(image, sweep_panel, 28, "#FCFAFF", "#DAC6F0", width=3, shadow_alpha=12)
    draw = ImageDraw.Draw(image, "RGBA")

    draw.text((112, 248), "Core learner path", font=H2, fill=CORE_COLOR)
    draw.text((114, 300), "Use included teaching data and saved CoGAPS evidence.", font=SMALL_BOLD, fill=MUTED)
    badge(draw, 942, 252, "no model sweep", CORE_COLOR, SMALL_BOLD)

    y = 370
    h = 300
    gap = 52
    boxes = []
    for step in CORE_STEPS:
        box = (130, y, 1168, y + h)
        boxes.append((f"core-{step.n}", draw_workflow_step(image, audit, step, box, CORE_COLOR)))
        y += h + gap
    draw = ImageDraw.Draw(image, "RGBA")
    for _, box1 in boxes[:-1]:
        idx = boxes.index((_, box1))
        box2 = boxes[idx + 1][1]
        arrow(draw, ((box1[0] + box1[2]) // 2, box1[3] + 8), ((box2[0] + box2[2]) // 2, box2[1] - 12), "#9DB0C5", width=7)


    draw.text((1332, 248), "Full local reproduction", font=H2, fill=REBUILD_COLOR)
    draw.text((1334, 300), "Optional rebuilding with source and large files.", font=SMALL_BOLD, fill=MUTED)
    y = 382
    rebuild_boxes = []
    for step in REBUILD_STEPS:
        box = (1332, y, 1890, y + 314)
        rebuild_boxes.append((f"rebuild-{step.n}", draw_workflow_step(image, audit, step, box, REBUILD_COLOR)))
        y += 348
    draw = ImageDraw.Draw(image, "RGBA")
    for idx, (_, box1) in enumerate(rebuild_boxes[:-1]):
        box2 = rebuild_boxes[idx + 1][1]
        arrow(draw, ((box1[0] + box1[2]) // 2, box1[3] + 8), ((box2[0] + box2[2]) // 2, box2[1] - 12), rgba(REBUILD_COLOR, 210), width=6)
    arrow(draw, (1330, 1282), (1212, 1188), REBUILD_COLOR, width=6, dashed=True)
    draw.text((1358, 1422), "rejoins core at selected outputs", font=TINY_BOLD, fill=REBUILD_COLOR)

    draw.text((1332, 1788), "Separate HPC sweep", font=H2, fill=SWEEP_COLOR)
    draw.text((1334, 1840), "Not part of the main lesson or full local reproduction.", font=SMALL_BOLD, fill=MUTED)
    draw_cluster(draw, 1370, 1908, SWEEP_COLOR)
    draw_file_stack(draw, 1560, 1906, SWEEP_COLOR)
    arrow(draw, (1496, 1970), (1544, 1970), SWEEP_COLOR, width=6)
    draw_text_block(draw, audit, "sweep-body", 1340, 2070, "The full rank sweep was run on high-performance computing. The case study uses included summary tables so learners can evaluate rank evidence without launching jobs.", BODY, 530, INK, sweep_panel, 8)
    badge(draw, 1340, 2320, "HPC only", SWEEP_COLOR, SMALL_BOLD)
    badge(draw, 1488, 2320, "summary tables included", SWEEP_COLOR, SMALL_BOLD)

    footer = (70, 2720, 1932, 2802)
    draw.rounded_rectangle(footer, radius=20, fill=WHITE, outline="#D8E4EF", width=2)
    draw.text((104, 2742), "Reading guide:", font=SMALL_BOLD, fill=NAVY)
    draw.text((268, 2742), "Follow the left path for the case study. Use the upper-right path only to rebuild selected files. Treat the lower-right HPC sweep as provenance.", font=SMALL, fill=MUTED)

    save_with_qa(image, audit, WORKFLOW_PNG, "Case Study 4 setup workflow and path-choice map, Pillow vector draft")


def save_with_qa(image: Image.Image, audit: Audit, png_path: Path, rendering: str) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    image.convert("RGB").save(png_path, quality=95)
    qa = {
        "output": str(png_path),
        "size": list(image.size),
        "rendering": rendering,
        "audit_item_count": len(audit.items),
        "audit_failures": audit.failures,
        "status": "pass" if not audit.failures else "review",
    }
    print(json.dumps(qa, indent=2))


def main() -> None:
    generate_resource_figure()
    generate_workflow_figure()


if __name__ == "__main__":
    main()
