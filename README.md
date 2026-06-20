# Lectin H-Score Pipeline

Pixel-area based H-Score quantification pipeline for lectin histochemistry images (DAB staining).  
Developed for comparative analysis of lectin binding patterns (WGA, SNA, MAL) across multiple regions and tissue types.

## Pipeline

| Script | Description |
|--------|-------------|
| `analyze_dab.py` | Color deconvolution (Ruifrok & Johnston 2001) + H-Score calculation + overlay visualization |
| `final_stats_and_plots.py` | Kruskal-Wallis + Mann-Whitney U post-hoc + significance-annotated bar plots |
| `generate_summary.py` | Mean ± SD summary table per region/tissue/lectin |
| `save_as_word.py` | Generates Materials & Methods section as Word document |

## Method

H-Score is calculated from pixel-area proportions (0–300 scale):

```
H-Score = (1 × %Weak) + (2 × %Moderate) + (3 × %Strong)
```

DAB intensity thresholds (OD units):
- Negative: < 0.22
- Weak (1+): 0.22 – 0.35
- Moderate (2+): 0.35 – 0.50
- Strong (3+): ≥ 0.50

## Requirements

```bash
pip install scikit-image numpy pandas matplotlib seaborn scipy python-docx
```

## Folder Structure

```
project/
├── kesitler/              # Input images (.jpg/.tif/.png)
├── BOLGE1_SONUCLAR/       # Per-lectin/tissue CSVs
├── BOLGE2_SONUCLAR/
├── analiz_cikti/          # Auto-created: H-score overlay PNGs
└── istatistik_analiz/     # Auto-created: stats + plots
```

## Reference

Ruifrok AC, Johnston DA (2001). Quantification of histochemical staining by color deconvolution. *Anal Quant Cytol Histol* 23(4):291–299.
