"""
Lectin H-Score Pipeline - Algorithm Figure Generator
Generates a publication-ready workflow figure (300 DPI, PDF + PNG)
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

# ── Color palette (colorblind-safe) ──────────────────────────────────────────
C_INPUT    = "#4E79A7"   # blue
C_PROC     = "#F28E2B"   # orange
C_ORG      = "#76B7B2"   # teal
C_STAT     = "#59A14F"   # green
C_SUM      = "#EDC948"   # yellow
C_DOC      = "#B07AA1"   # purple
C_OUT      = "#E15759"   # red
C_ARROW    = "#4D4D4D"
C_BG       = "#FAFAFA"
WHITE      = "#FFFFFF"
DARK       = "#2D2D2D"

FIG_W, FIG_H = 14, 20
FONT_FAMILY  = "DejaVu Sans"

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), facecolor=C_BG)
ax.set_xlim(0, 14)
ax.set_ylim(0, 20)
ax.axis("off")
fig.patch.set_facecolor(C_BG)


# ── Helper functions ──────────────────────────────────────────────────────────

def draw_box(ax, x, y, w, h, label, sublabel=None,
             color=C_PROC, alpha=0.92, fontsize=10, radius=0.25,
             icon=None, text_color=WHITE):
    """Draw a rounded rectangle with centered text."""
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=1.4, edgecolor=DARK, facecolor=color, alpha=alpha,
        zorder=3
    )
    ax.add_patch(box)

    # Icon + main label
    full_label = f"{icon}  {label}" if icon else label
    ax.text(x, y + (0.15 if sublabel else 0), full_label,
            ha="center", va="center", fontsize=fontsize,
            fontweight="bold", color=text_color, fontfamily=FONT_FAMILY,
            zorder=4, wrap=True)

    if sublabel:
        ax.text(x, y - 0.35, sublabel,
                ha="center", va="center", fontsize=8,
                color=text_color, alpha=0.88, fontfamily=FONT_FAMILY,
                zorder=4, style="italic")


def draw_diamond(ax, x, y, w, h, label, color=C_ORG, fontsize=9, text_color=WHITE):
    """Draw a diamond shape (decision/branch node)."""
    dx, dy = w / 2, h / 2
    diamond = plt.Polygon(
        [[x, y + dy], [x + dx, y], [x, y - dy], [x - dx, y]],
        closed=True, linewidth=1.4, edgecolor=DARK,
        facecolor=color, alpha=0.92, zorder=3
    )
    ax.add_patch(diamond)
    ax.text(x, y, label, ha="center", va="center",
            fontsize=fontsize, fontweight="bold",
            color=text_color, fontfamily=FONT_FAMILY, zorder=4)


def arrow(ax, x1, y1, x2, y2, label=None, color=C_ARROW, lw=1.8):
    """Draw a straight arrow with optional label."""
    ax.annotate("",
                xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color,
                                lw=lw, mutation_scale=14),
                zorder=2)
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mx + 0.12, my, label, fontsize=7.5, color=color,
                fontfamily=FONT_FAMILY, va="center", zorder=5)


def section_label(ax, x, y, text, color=DARK):
    ax.text(x, y, text, fontsize=8, color=color, fontfamily=FONT_FAMILY,
            alpha=0.7, ha="left", va="center", style="italic", zorder=5)


# ═══════════════════════════════════════════════════════════════════════════════
# TITLE
# ═══════════════════════════════════════════════════════════════════════════════
ax.text(7, 19.4, "Lectin Histochemistry H-Score Quantification Pipeline",
        ha="center", va="center", fontsize=13, fontweight="bold",
        color=DARK, fontfamily=FONT_FAMILY)
ax.text(7, 19.0, "DAB staining analysis · Color deconvolution · Statistical comparison",
        ha="center", va="center", fontsize=9, color=DARK,
        fontfamily=FONT_FAMILY, alpha=0.65, style="italic")

# Horizontal rule under title
ax.plot([0.5, 13.5], [18.72, 18.72], color=DARK, lw=0.8, alpha=0.3, zorder=2)


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 1 — INPUT
# ═══════════════════════════════════════════════════════════════════════════════
section_label(ax, 0.35, 18.3, "① INPUT")
draw_box(ax, 7, 18.0, 6.5, 0.75,
         "Lectin-Stained Histology Images",
         sublabel="kesitler/  ·  .tif / .png / .jpg  ·  WGA · SNA · MAL",
         color=C_INPUT, fontsize=10)

arrow(ax, 7, 17.62, 7, 17.18)


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 2 — COLOR DECONVOLUTION (analyze_dab.py)
# ═══════════════════════════════════════════════════════════════════════════════
section_label(ax, 0.35, 16.95, "② analyze_dab.py")

# Main process box
draw_box(ax, 7, 16.55, 10.5, 0.85,
         "Color Deconvolution  (Ruifrok & Johnston, 2001)",
         sublabel="RGB → Optical Density  ·  Stain matrix unmixing  ·  DAB channel extraction",
         color=C_PROC, fontsize=10)

arrow(ax, 7, 16.12, 7, 15.62)

# OD sub-step
draw_box(ax, 7, 15.28, 10.5, 0.68,
         "Optical Density (OD) Calculation",
         sublabel="OD = −log₁₀(RGB / 255)   ·   Lambert–Beer law",
         color=C_PROC, fontsize=9.5)

arrow(ax, 7, 14.94, 7, 14.42)

# Classification box
draw_box(ax, 7, 14.08, 10.5, 0.68,
         "Intensity Classification",
         sublabel="Negative (<0.22 OD)  ·  Weak 1+ (0.22–0.35)  ·  Moderate 2+ (0.35–0.50)  ·  Strong 3+ (≥0.50)",
         color=C_PROC, fontsize=9.5)

arrow(ax, 7, 13.74, 7, 13.22)

# H-Score box
draw_box(ax, 7, 12.88, 10.5, 0.68,
         "H-Score Calculation  (0 – 300 scale)",
         sublabel="H = (1 × %Weak) + (2 × %Moderate) + (3 × %Strong)",
         color=C_PROC, fontsize=9.5)

# Output annotations right side
ax.annotate("", xy=(12.6, 12.35), xytext=(12.6, 16.12),
            arrowprops=dict(arrowstyle="-", color=C_PROC,
                            lw=1.2, linestyle="dashed"), zorder=2)
ax.text(12.75, 14.2, "analiz_cikti/\n*.png overlays\n+\nlektin_hscore\n_sonuclari.csv",
        fontsize=7.2, color=C_PROC, fontfamily=FONT_FAMILY,
        va="center", ha="left", linespacing=1.4, zorder=5)

arrow(ax, 7, 12.54, 7, 12.0)


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 3 — DATA ORGANIZATION
# ═══════════════════════════════════════════════════════════════════════════════
section_label(ax, 0.35, 11.8, "③ DATA ORGANIZATION")
draw_box(ax, 7, 11.55, 10.5, 0.75,
         "Region-Specific CSV Files",
         sublabel="lektin_hscore_{PROJECT}_{REGION}_{TISSUE}_{LECTIN}.csv   ·   Alsancak · Pasaport · Urla",
         color=C_ORG, fontsize=10)

# Branch arrows to three parallel modules
arrow(ax, 4.5, 11.17, 2.5,  10.55)   # left
arrow(ax, 7,   11.17, 7,    10.55)   # center
arrow(ax, 9.5, 11.17, 11.5, 10.55)  # right


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 4 — THREE PARALLEL ANALYSIS MODULES
# ═══════════════════════════════════════════════════════════════════════════════
section_label(ax, 0.35, 10.4, "④ ANALYSIS MODULES  (parallel)")

# ── Left: Statistical analysis ─────────────────────────────
draw_box(ax, 2.5, 9.95, 4.5, 1.08,
         "final_stats_and_plots.py",
         sublabel="Kruskal–Wallis H-test\nMann–Whitney U (pairwise)\nBar + jitter plots (SE)",
         color=C_STAT, fontsize=9)
arrow(ax, 2.5, 9.41, 2.5, 8.85)

# ── Center: Summary table ──────────────────────────────────
draw_box(ax, 7, 9.95, 4.5, 1.08,
         "generate_summary.py",
         sublabel="Mean ± SD per group\nRegion / Tissue / Lectin\nhierarchy",
         color=C_SUM, fontsize=9, text_color=DARK)
arrow(ax, 7, 9.41, 7, 8.85)

# ── Right: Methods document ───────────────────────────────
draw_box(ax, 11.5, 9.95, 4.5, 1.08,
         "save_as_word.py",
         sublabel="Materials & Methods\nColor deconvolution\nStatistics section",
         color=C_DOC, fontsize=9)
arrow(ax, 11.5, 9.41, 11.5, 8.85)


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE 5 — OUTPUTS
# ═══════════════════════════════════════════════════════════════════════════════
section_label(ax, 0.35, 8.72, "⑤ OUTPUTS")

draw_box(ax, 2.5, 8.35, 4.5, 0.92,
         "istatistik_analiz/",
         sublabel="plot_*.png\nPublication figures\np-values annotated (*//**//**/**)",
         color=C_OUT, fontsize=8.5)

draw_box(ax, 7, 8.35, 4.5, 0.92,
         "mean_sd_ozet.txt",
         sublabel="Hierarchical summary\nMean ± SD · n values\nManuscript-ready",
         color=C_OUT, fontsize=8.5)

draw_box(ax, 11.5, 8.35, 4.5, 0.92,
         "materyal_metot.docx",
         sublabel="Methods section\nPeer-review ready\nWord format",
         color=C_OUT, fontsize=8.5)

# Connect all three outputs to final result box
arrow(ax, 2.5,  7.88, 4.2, 7.30)
arrow(ax, 7,    7.88, 7,   7.30)
arrow(ax, 11.5, 7.88, 9.8, 7.30)


# ═══════════════════════════════════════════════════════════════════════════════
# FINAL RESULT BOX
# ═══════════════════════════════════════════════════════════════════════════════
draw_box(ax, 7, 6.88, 10.5, 0.80,
         "Publication-Ready H-Score Dataset",
         sublabel="Quantified DAB staining · Statistical comparisons · Figures · Documentation",
         color=C_INPUT, fontsize=10.5)


# ═══════════════════════════════════════════════════════════════════════════════
# LEGEND — Stain matrix & OD thresholds reference table
# ═══════════════════════════════════════════════════════════════════════════════
legend_x, legend_y = 0.6, 5.9
legend_w, legend_h = 12.8, 5.5

# Background
legend_bg = FancyBboxPatch(
    (legend_x, legend_y - legend_h), legend_w, legend_h,
    boxstyle="round,pad=0,rounding_size=0.2",
    linewidth=1, edgecolor="#BBBBBB", facecolor=WHITE, alpha=0.8, zorder=2
)
ax.add_patch(legend_bg)

ax.text(7, 5.7, "Technical Reference",
        ha="center", va="center", fontsize=9, fontweight="bold",
        color=DARK, fontfamily=FONT_FAMILY, zorder=5)

# ── OD Threshold table ──────────────────────────────────────
col_headers = ["Category", "OD Range", "Score", "Color"]
col_x       = [1.8, 3.6, 5.2, 6.5]
rows = [
    ("Negative",  "< 0.22",      "—",  "#D3D3D3"),
    ("Weak (1+)", "0.22 – 0.35", "1",  "#FFC0CB"),
    ("Moderate (2+)", "0.35 – 0.50", "2", "#FF8C00"),
    ("Strong (3+)",   "≥ 0.50",      "3", "#B22222"),
]

row_y = 5.3
ax.text(col_x[0], row_y, "OD Intensity Thresholds",
        fontsize=8, fontweight="bold", color=DARK,
        fontfamily=FONT_FAMILY, va="center", zorder=5)
row_y -= 0.35
for h, cx in zip(col_headers, col_x):
    ax.text(cx, row_y, h, fontsize=7.5, fontweight="bold", color=DARK,
            fontfamily=FONT_FAMILY, va="center", zorder=5)
row_y -= 0.28
ax.plot([1.4, 7.4], [row_y + 0.1, row_y + 0.1], color="#BBBBBB", lw=0.8, zorder=4)
for cat, od, score, clr in rows:
    row_y -= 0.32
    ax.text(col_x[0], row_y, cat,   fontsize=7.5, color=DARK, fontfamily=FONT_FAMILY, va="center", zorder=5)
    ax.text(col_x[1], row_y, od,    fontsize=7.5, color=DARK, fontfamily=FONT_FAMILY, va="center", zorder=5)
    ax.text(col_x[2], row_y, score, fontsize=7.5, color=DARK, fontfamily=FONT_FAMILY, va="center", zorder=5)
    swatch = FancyBboxPatch((col_x[3] - 0.1, row_y - 0.1), 0.55, 0.22,
                             boxstyle="round,pad=0", linewidth=0.5,
                             edgecolor=DARK, facecolor=clr, zorder=5)
    ax.add_patch(swatch)

# ── Stain matrix ────────────────────────────────────────────
mx, my = 8.5, 5.3
ax.text(mx, my, "Color Deconvolution Stain Matrix (ImageJ standard)",
        fontsize=8, fontweight="bold", color=DARK,
        fontfamily=FONT_FAMILY, va="center", zorder=5)
matrix_rows = [
    ("DAB",         "[0.650", "0.704", "0.286]"),
    ("Haematoxylin","[0.268", "0.570", "0.776]"),
    ("Residual",    "[0.711", "0.423", "0.561]"),
]
col2 = [8.5, 9.85, 10.85, 11.75]
row2_y = my - 0.32
for ch_name, r, g, b in matrix_rows:
    row2_y -= 0.32
    ax.text(col2[0], row2_y, ch_name, fontsize=7.5, color=DARK,
            fontfamily=FONT_FAMILY, va="center", zorder=5)
    for val, cx in zip([r, g, b], col2[1:]):
        ax.text(cx, row2_y, val, fontsize=7.5, color=DARK,
                fontfamily=FONT_FAMILY, va="center", zorder=5)

# H-Score formula
form_y = row2_y - 0.5
ax.text(8.5, form_y,
        "H-Score  =  (1 × %Weak)  +  (2 × %Moderate)  +  (3 × %Strong)   ∈  [0, 300]",
        fontsize=8, color=C_PROC, fontweight="bold",
        fontfamily=FONT_FAMILY, va="center", zorder=5)

# Statistics section
stat_y = form_y - 0.4
ax.text(8.5, stat_y,
        "Statistics:  Kruskal–Wallis H-test  →  Post-hoc Mann–Whitney U  "
        "·  * p<0.05  ** p<0.01  *** p<0.001",
        fontsize=7.5, color=DARK, fontfamily=FONT_FAMILY,
        va="center", zorder=5, alpha=0.85)


# ═══════════════════════════════════════════════════════════════════════════════
# STAGE LEGEND (color key)
# ═══════════════════════════════════════════════════════════════════════════════
leg_items = [
    (C_INPUT, "Input / Final output"),
    (C_PROC,  "Image processing"),
    (C_ORG,   "Data organization"),
    (C_STAT,  "Statistical analysis"),
    (C_SUM,   "Summary / Tables"),
    (C_DOC,   "Documentation"),
    (C_OUT,   "File outputs"),
]
lx, ly = 0.65, 0.85
ax.text(lx, ly, "Stage legend:", fontsize=7.5, color=DARK,
        fontfamily=FONT_FAMILY, va="center", fontweight="bold", zorder=5)
for i, (clr, lbl) in enumerate(leg_items):
    px = lx + 1.2 + i * 1.83
    swatch = FancyBboxPatch((px, ly - 0.12), 0.28, 0.24,
                             boxstyle="round,pad=0", linewidth=0.5,
                             edgecolor=DARK, facecolor=clr, zorder=5, alpha=0.92)
    ax.add_patch(swatch)
    ax.text(px + 0.35, ly, lbl, fontsize=7, color=DARK,
            fontfamily=FONT_FAMILY, va="center", zorder=5)

# Footer
ax.text(13.5, 0.22,
        "github.com/erdal-b-dot/lectin-hscore-pipeline",
        ha="right", fontsize=7, color=DARK, alpha=0.5,
        fontfamily=FONT_FAMILY, zorder=5, style="italic")


# ═══════════════════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════════════════
plt.tight_layout(pad=0.3)

out_pdf = "/Users/erdalbalcan/Desktop/lectin_hscore_pipeline_figure.pdf"
out_png = "/Users/erdalbalcan/Desktop/lectin_hscore_pipeline_figure.png"

fig.savefig(out_pdf, dpi=300, bbox_inches="tight", facecolor=C_BG)
fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor=C_BG)

print(f"Saved:\n  {out_pdf}\n  {out_png}")
plt.show()
