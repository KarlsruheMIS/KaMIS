#!/usr/bin/env python3
"""Generate the GitHub social preview image (1280x640) for KaMIS."""

from PIL import Image, ImageDraw, ImageFont
import math
import os
import random

random.seed(42)

W, H = 1280, 640
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "social-preview.png")

canvas = Image.new("RGB", (W, H), (15, 23, 42))
draw = ImageDraw.Draw(canvas)

# Subtle radial gradient
cx, cy = W // 2, H // 2
for r in range(500, 0, -2):
    frac = r / 500
    c = int(15 + (30 - 15) * (1 - frac))
    c2 = int(23 + (45 - 23) * (1 - frac))
    c3 = int(42 + (72 - 42) * (1 - frac))
    draw.ellipse([cx - r * 1.6, cy - r, cx + r * 1.6, cy + r], fill=(c, c2, c3))

# Graph vertices
positions = [
    (180, 200), (280, 130), (380, 180), (320, 300), (200, 340),
    (140, 260), (420, 280), (260, 220), (350, 350), (160, 150),
    (440, 150), (500, 240), (480, 350), (100, 310), (230, 400),
    (370, 420), (500, 400), (540, 300), (120, 420), (300, 450),
    (60, 180), (550, 130), (80, 450), (560, 450), (450, 480),
    (170, 490), (340, 510), (500, 500),
]

# Edges between nearby vertices
edges = []
for i in range(len(positions)):
    for j in range(i + 1, len(positions)):
        x1, y1 = positions[i]
        x2, y2 = positions[j]
        dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        if dist < 160:
            edges.append((i, j))

# Greedy independent set
adj = {i: set() for i in range(len(positions))}
for i, j in edges:
    adj[i].add(j)
    adj[j].add(i)

independent_set = set()
excluded = set()
for v in sorted(range(len(positions)), key=lambda v: len(adj[v])):
    if v not in excluded:
        independent_set.add(v)
        excluded.add(v)
        excluded.update(adj[v])

# Colors
edge_color = (50, 65, 90)
vertex_color = (100, 120, 160)
is_color = (59, 180, 120)

# Draw edges
for i, j in edges:
    x1, y1 = positions[i]
    x2, y2 = positions[j]
    draw.line([(x1, y1), (x2, y2)], fill=edge_color, width=1)

# Glow for independent set vertices
for i in independent_set:
    x, y = positions[i]
    for r in range(20, 4, -1):
        alpha_frac = 1 - (r - 4) / 16
        gc = tuple(
            int(is_color[k] * 0.15 * alpha_frac + (1 - 0.15 * alpha_frac) * [15, 23, 42][k])
            for k in range(3)
        )
        draw.ellipse([x - r, y - r, x + r, y + r], fill=gc)

# Draw vertices
for i, (x, y) in enumerate(positions):
    if i in independent_set:
        r = 8
        draw.ellipse([x - r, y - r, x + r, y + r], fill=is_color, outline=(80, 220, 150), width=2)
    else:
        r = 5
        draw.ellipse([x - r, y - r, x + r, y + r], fill=vertex_color, outline=(130, 150, 190), width=1)

# Fonts
try:
    font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 80)
    font_sub = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 28)
    font_tag = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
    font_ver = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
except OSError:
    font_title = font_sub = font_tag = font_ver = ImageFont.load_default()

text_x = 660

# Title
draw.text((text_x, 180), "KaMIS", fill=(240, 245, 255), font=font_title)

# Version pill
ver_text = "v3.2"
vbbox = draw.textbbox((0, 0), ver_text, font=font_ver)
vw = vbbox[2] - vbbox[0]
title_bbox = draw.textbbox((text_x, 180), "KaMIS", font=font_title)
pill_x = title_bbox[2] + 18
pill_y = 210
draw.rounded_rectangle([pill_x, pill_y, pill_x + vw + 18, pill_y + 32], radius=14, fill=is_color)
draw.text((pill_x + 9, pill_y + 3), ver_text, fill=(255, 255, 255), font=font_ver)

# Separator
draw.line([(text_x, 290), (text_x + 500, 290)], fill=(60, 80, 110), width=2)

# Subtitle
draw.text((text_x, 310), "Karlsruhe Maximum", fill=(180, 195, 220), font=font_sub)
draw.text((text_x, 348), "Independent Sets", fill=(180, 195, 220), font=font_sub)

# Tagline
draw.text((text_x, 420), "Fast solvers for the maximum independent", fill=(110, 130, 160), font=font_tag)
draw.text((text_x, 450), "set problem on large sparse graphs", fill=(110, 130, 160), font=font_tag)

# Decorative dots
for i, dx in enumerate(range(0, 300, 40)):
    dot_x = text_x + dx
    dot_y = 540
    if i % 2 == 0:
        draw.ellipse([dot_x - 5, dot_y - 5, dot_x + 5, dot_y + 5], fill=is_color)
    else:
        draw.ellipse([dot_x - 3, dot_y - 3, dot_x + 3, dot_y + 3], fill=(60, 80, 110))

canvas.save(OUT_PATH, "PNG", quality=95)
print(f"Saved {OUT_PATH}")
