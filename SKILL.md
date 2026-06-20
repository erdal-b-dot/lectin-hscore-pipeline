---
name: "lectin-hscore-pipeline"
description: "Lektin histokimyası görüntülerinde piksel-alan tabanlı DAB renk dekonvolüsyonu ve H-Skor (0–300) hesaplama pipeline'ı. Çok bölge / çok lektin (WGA, SNA, MAL) karşılaştırmalı analizleri, Kruskal-Wallis + Mann-Whitney U istatistikleri ve anlamlılık yıldızlı bar grafikleri üretir. Midye (Mytilus), balık, deniz canlıları veya diğer dokulardaki lektin bağlanma örüntülerini sayısal olarak karşılaştırmak için kullanılır."
license: "CC-BY-4.0"
---

# Lectin H-Score Pipeline

## Overview

Bu skill, lektin histokimyası preparatlarında DAB (3,3′-diaminobenzidin) boyanma yoğunluğunu **piksel-alan tabanlı H-Skor** yöntemiyle sayısal olarak ölçer. Hücre segmentasyonu gerektirmez; bunun yerine her pikselin Optik Dansite (OD) değerini Ruifrok & Johnston (2001) renk dekonvolüsyonu ile hesaplar ve üç yoğunluk kategorisine (1+/2+/3+) ayırır.

**Pipeline 4 aşamadan oluşur:**
1. `analyze_dab.py` — DAB dekonvolüsyonu + H-Skor + overlay görselleştirme (görüntü başına CSV satırı)
2. Bölge/doku/lektin bazlı CSV'lerin `{BOLGE}_SONUCLAR/` klasörlerine taşınması
3. `final_stats_and_plots.py` — Tüm CSV'lerin birleştirilmesi, Kruskal-Wallis + post-hoc bar grafikleri
4. `generate_summary.py` — Mean ± SD özet tablosu

**Kullanılan kütüphaneler:** `scikit-image`, `numpy`, `pandas`, `matplotlib`, `seaborn`, `scipy`

---

## When to Use

- Lektin histokimyası preparatlarında (WGA, SNA, MAL veya diğer lektinler) DAB pozitifliğini sayısal olarak ölçmek
- Birden fazla bölge (lokasyon, grup, tedavi) ve birden fazla lektin tipini aynı anda karşılaştırmak
- Hücre sınırlarının belirsiz olduğu (mukus, bağ doku, gonad asini) preparatlarda piksel-tabanlı analiz yapmak
- Yayına hazır istatistiksel bar grafikleri (Mean ± SE + jitter + anlamlılık yıldızı) üretmek
- Materyal & Metot bölümü için standart metin oluşturmak

**Yerine `ihc-quantification` kullanın** eğer: hücre bazlı sayım yapmanız gerekiyorsa (Ki-67, CD marker), veya QuPath/TMA çıktısı istiyorsanız.

---

## Prerequisites

```bash
pip install scikit-image numpy pandas matplotlib seaborn scipy python-docx
```

**Klasör yapısı:**
```
proje_klasoru/
├── analyze_dab.py
├── final_stats_and_plots.py
├── generate_summary.py
├── save_as_word.py          # isteğe bağlı
├── kesitler/                # Ham görüntüler buraya (.jpg/.tif/.png)
├── BOLGE1_SONUCLAR/         # CSV'ler buraya taşınır (per lektin/doku)
├── BOLGE2_SONUCLAR/
├── analiz_cikti/            # Otomatik oluşturulur (H-skor overlay PNG'ler)
└── istatistik_analiz/       # Otomatik oluşturulur (istatistik + grafik)
```

**CSV adlandırma kuralı** (`final_stats_and_plots.py`'nin parse edebilmesi için):
```
lektin_hscore_{PROJE}_{BOLGE}_{DOKU}_{LEKTIN}.csv
# Örnek:
lektin_hscore_midye_ALSANCAK_GENEL_WGA.csv
lektin_hscore_midye_URLA_GONAD_SNA.csv
lektin_hscore_midye_URLA_GONAD_DISI_MAL.csv   # GONAD_DISI/GONAD_ERKEK gibi alt gruplar da desteklenir
```

---

## Quick Start

Tek görüntü üzerinde H-Skor hesapla:

```python
import numpy as np
from skimage import io, color

def imagej_color_deconvolution(rgb_img):
    """ImageJ-uyumlu renk dekonvolüsyonu (Ruifrok & Johnston 2001)."""
    rgb_float = rgb_img.astype(float) + 1.0
    od = -np.log10(rgb_float / 255.0)
    stain_matrix = np.array([
        [0.650, 0.704, 0.286],   # Hematoksilen
        [0.268, 0.570, 0.776],   # DAB
        [0.711, 0.423, 0.561]    # Artık
    ])
    stain_matrix /= np.linalg.norm(stain_matrix, axis=1)[:, np.newaxis]
    inverse_matrix = np.linalg.inv(stain_matrix)
    return np.dot(od, inverse_matrix)

img = io.imread("kesitler/ornek.jpg")
if img.shape[-1] == 4:
    img = color.rgba2rgb(img)

deconvolved = imagej_color_deconvolution(img)
dab_od = np.maximum(deconvolved[:, :, 1], 0)   # Kanal 1 = DAB

# Eşik değerleri (OD birimleri)
BASE_THRESH = 0.22
weak_mask   = (dab_od >= 0.22) & (dab_od < 0.35)
mod_mask    = (dab_od >= 0.35) & (dab_od < 0.50)
strong_mask = (dab_od >= 0.50)
total_px    = dab_od.size

p_weak   = weak_mask.sum()   / total_px * 100
p_mod    = mod_mask.sum()    / total_px * 100
p_strong = strong_mask.sum() / total_px * 100

h_score = (1 * p_weak) + (2 * p_mod) + (3 * p_strong)
print(f"H-Skor: {h_score:.2f}  |  Zayıf: {p_weak:.1f}%  Orta: {p_mod:.1f}%  Güçlü: {p_strong:.1f}%")
```

---

## Workflow

### Aşama 1 — DAB Analizi ve H-Skor Overlay (`analyze_dab.py`)

Tüm görüntüleri `kesitler/` klasöründen okur, H-Skor hesaplar, overlay PNG kaydeder, sonuçları CSV'ye yazar.

```python
import os
import numpy as np
import pandas as pd
from skimage import io, color
import matplotlib.pyplot as plt

def imagej_color_deconvolution(rgb_img):
    rgb_float = rgb_img.astype(float) + 1.0
    od = -np.log10(rgb_float / 255.0)
    stain_matrix = np.array([
        [0.650, 0.704, 0.286],
        [0.268, 0.570, 0.776],
        [0.711, 0.423, 0.561]
    ])
    stain_matrix /= np.linalg.norm(stain_matrix, axis=1)[:, np.newaxis]
    return np.dot(od, np.linalg.inv(stain_matrix))

def analyze_h_score(image_path, base_thresh=0.22):
    img = io.imread(image_path)
    if img.shape[-1] == 4:
        img = color.rgba2rgb(img)
    dab_od = np.maximum(imagej_color_deconvolution(img)[:, :, 1], 0)

    weak_mask   = (dab_od >= base_thresh) & (dab_od < 0.35)
    mod_mask    = (dab_od >= 0.35)        & (dab_od < 0.50)
    strong_mask = (dab_od >= 0.50)
    total_pos   = dab_od >= base_thresh

    total_px = dab_od.size
    p_weak   = weak_mask.sum()   / total_px * 100
    p_mod    = mod_mask.sum()    / total_px * 100
    p_strong = strong_mask.sum() / total_px * 100
    p_total  = total_pos.sum()   / total_px * 100
    h_score  = (1 * p_weak) + (2 * p_mod) + (3 * p_strong)
    mean_od  = dab_od[total_pos].mean() if total_pos.any() else 0

    return {'mean_od': mean_od, 'p_total': p_total,
            'p_weak': p_weak, 'p_mod': p_mod, 'p_strong': p_strong,
            'h_score': h_score, 'dab_od': dab_od}

def main():
    input_folder  = 'kesitler'
    output_plots  = 'analiz_cikti'
    BASE_THRESH   = 0.22
    os.makedirs(output_plots, exist_ok=True)
    results = []

    for filename in os.listdir(input_folder):
        if not filename.lower().endswith(('.tif', '.png', '.jpg', '.jpeg')):
            continue
        path = os.path.join(input_folder, filename)
        data = analyze_h_score(path, BASE_THRESH)
        if data is None:
            continue

        results.append({
            'File_Name': filename,
            'Mean_OD': round(data['mean_od'], 4),
            'Total_Positive_Area_Pct': round(data['p_total'], 2),
            'Weak_1_Plus_Area_Pct':    round(data['p_weak'], 2),
            'Moderate_2_Plus_Area_Pct':round(data['p_mod'], 2),
            'Strong_3_Plus_Area_Pct':  round(data['p_strong'], 2),
            'H_Score': round(data['h_score'], 2)
        })

        # H-skor yoğunluk haritası overlay
        h_map = np.zeros_like(data['dab_od'])
        h_map[data['dab_od'] >= 0.22] = 1
        h_map[data['dab_od'] >= 0.35] = 2
        h_map[data['dab_od'] >= 0.50] = 3

        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1); plt.imshow(io.imread(path)); plt.title('Orijinal'); plt.axis('off')
        plt.subplot(1, 2, 2); plt.imshow(h_map, cmap='YlOrRd')
        plt.title(f"H-Skor: {round(data['h_score'], 1)}"); plt.axis('off')
        plt.savefig(os.path.join(output_plots, f'hscore_{filename}.png')); plt.close()
        print(f"{filename}: H-Skor = {round(data['h_score'], 2)}")

    pd.DataFrame(results).to_csv('lektin_hscore_sonuclari.csv', index=False, sep=';', decimal=',')
    print("Tamamlandı → lektin_hscore_sonuclari.csv")

if __name__ == "__main__":
    main()
```

**Çalıştırma:**
```bash
python analyze_dab.py
```

---

### Aşama 2 — CSV'leri Sınıflandır

`analyze_dab.py` çıktısını grup/doku/lektin'e göre ayır ve ilgili `_SONUCLAR/` klasörlerine koy.  
**CSV adı kuralına dikkat et** (yukarıdaki Prerequisites'e bak).

---

### Aşama 3 — İstatistik ve Grafikler (`final_stats_and_plots.py`)

Tüm `_SONUCLAR/` klasörlerini okur, birleştirir, Kruskal-Wallis + Mann-Whitney U yapar, anlamlılık yıldızlı bar grafikleri üretir.

```python
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

base_path     = '/path/to/proje_klasoru'    # ← DEĞİŞTİR
sonuclar_dirs = ['BOLGE1_SONUCLAR', 'BOLGE2_SONUCLAR', 'BOLGE3_SONUCLAR']
output_dir    = os.path.join(base_path, 'istatistik_analiz')
os.makedirs(output_dir, exist_ok=True)

all_data = []
for sdir in sonuclar_dirs:
    for file in os.listdir(os.path.join(base_path, sdir)):
        if not file.endswith('.csv'):
            continue
        parts = file.replace('.csv', '').split('_')
        # Parse: lektin_hscore_{PROJE}_{BOLGE}_{DOKU}_{LEKTIN}
        region = parts[3].upper()
        name_upper = file.upper()
        if 'GONAD_DISI' in name_upper:
            tissue, lectin = 'GONAD_DISI', parts[-1].upper()
        elif 'GONAD_ERKEK' in name_upper:
            tissue, lectin = 'GONAD_ERKEK', parts[-1].upper()
        else:
            tissue, lectin = parts[4].upper(), parts[5].upper()

        df = pd.read_csv(os.path.join(base_path, sdir, file), sep=';', decimal=',')
        df['Region'] = region
        df['Tissue'] = tissue
        df['Lectin'] = lectin
        all_data.append(df)

full_df = pd.concat(all_data, ignore_index=True)
full_df['Tissue_Group'] = full_df['Tissue'].apply(lambda x: 'GONAD' if 'GONAD' in x else x)

def get_star(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return None

def annotate_significance(ax, sub_df, group_col, val_col):
    groups     = [t.get_text() for t in ax.get_xticklabels()]
    y_max      = sub_df[val_col].max()
    y_shift    = y_max * 0.1
    current_y  = y_max + y_shift
    for g1, g2 in combinations(groups, 2):
        d1 = sub_df[sub_df[group_col] == g1][val_col].dropna()
        d2 = sub_df[sub_df[group_col] == g2][val_col].dropna()
        if len(d1) > 0 and len(d2) > 0:
            _, p = mannwhitneyu(d1, d2)
            star = get_star(p)
            if star:
                i1, i2 = groups.index(g1), groups.index(g2)
                ax.plot([i1, i1, i2, i2],
                        [current_y - y_shift*0.2, current_y, current_y, current_y - y_shift*0.2],
                        lw=1.5, c='black')
                ax.text((i1+i2)/2, current_y + y_shift*0.1, star,
                        ha='center', va='bottom', fontsize=12, fontweight='bold')
                current_y += y_shift * 1.5

for tissue in full_df['Tissue_Group'].unique():
    for lectin in full_df['Lectin'].unique():
        sub = full_df[(full_df['Tissue_Group'] == tissue) & (full_df['Lectin'] == lectin)]
        if sub.empty or sub['Region'].nunique() < 2:
            continue
        order = sorted(sub['Region'].unique())
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.set_style("whitegrid")
        sns.barplot(data=sub, x='Region', y='H_Score', order=order,
                    capsize=.1, errorbar='se', palette='muted', ax=ax)
        sns.stripplot(data=sub, x='Region', y='H_Score', order=order,
                      color='black', alpha=0.4, jitter=True, ax=ax)
        ax.set_title(f'H-Skor: {tissue} — {lectin}', fontsize=14)
        ax.set_ylabel('H-Skor (Mean ± SE)', fontsize=12)
        annotate_significance(ax, sub, 'Region', 'H_Score')
        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 1.3)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'plot_{tissue}_{lectin}_stars.png'))
        plt.close()
        print(f"Kaydedildi: {tissue} — {lectin}")
```

---

### Aşama 4 — Mean ± SD Özet (`generate_summary.py`)

```python
import pandas as pd, os

base_path  = '/path/to/proje_klasoru'    # ← DEĞİŞTİR
csv_path   = os.path.join(base_path, 'istatistik_analiz', 'aggregated_hscore_data.csv')
output_txt = os.path.join(base_path, 'istatistik_analiz', 'mean_sd_ozet.txt')

df = pd.read_csv(csv_path, sep=';', decimal=',')
df['Tissue_Group'] = df['Tissue'].apply(lambda x: 'GONAD' if 'GONAD' in x else x)

summary = df.groupby(['Region', 'Tissue_Group', 'Lectin'])['H_Score'].agg(
    ['mean', 'std', 'count']).reset_index()

lines = ["H-SKOR ÖZET (Mean ± SD)\n" + "="*50 + "\n"]
for region in sorted(df['Region'].unique()):
    lines.append(f"BÖLGE: {region}")
    for tissue in sorted(df['Tissue_Group'].unique()):
        lines.append(f"  Doku: {tissue}")
        for lectin in sorted(df['Lectin'].unique()):
            row = summary[(summary['Region']==region) &
                          (summary['Tissue_Group']==tissue) &
                          (summary['Lectin']==lectin)]
            if not row.empty:
                m, s, n = row[['mean','std','count']].values[0]
                s = s if pd.notnull(s) else 0.0
                lines.append(f"    - {lectin:3}: {m:6.2f} ± {s:5.2f}  (n={int(n)})")
        lines.append("")

text = "\n".join(lines)
with open(output_txt, 'w', encoding='utf-8') as f:
    f.write(text)
print(text)
```

---

## Key Parameters

| Parametre | Varsayılan | Açıklama |
|-----------|-----------|----------|
| `base_thresh` | `0.22` | Pozitif kabul alt eşiği (OD birimi). Düşük boyanmalı örneklerde `0.15`'e indir. |
| Zayıf (1+) aralığı | `0.22 – 0.35` | H-Skor katkısı: ×1. Artır → daha seçici zayıf sınır. |
| Orta (2+) aralığı | `0.35 – 0.50` | H-Skor katkısı: ×2 |
| Güçlü (3+) eşiği | `≥ 0.50` | H-Skor katkısı: ×3. Çok az 3+ görünüyorsa `0.45`'e indir. |
| Stain matrix | ImageJ-uyumlu 3×3 | Kendi lektininiz için ImageJ → Colour Deconvolution ile kalibrasyon yapın. |
| `sep` / `decimal` | `';'` / `','` | Excel-TR formatı. İngilizce Excel için `','` / `'.'` olarak değiştirin. |
| İstatistik testi | Mann-Whitney U | Parametrik dağılım varsayımı yok; `kruskal` fonksiyonu ile genel F testi eklenebilir. |

---

## Common Recipes

### Tek görüntü için debug görselleştirme

```python
import matplotlib.pyplot as plt
import numpy as np
from skimage import io, color

def debug_image(image_path, base_thresh=0.22):
    img = io.imread(image_path)
    if img.shape[-1] == 4:
        img = color.rgba2rgb(img)

    # DAB kanalı
    rgb_float  = img.astype(float) + 1.0
    od         = -np.log10(rgb_float / 255.0)
    sm         = np.array([[0.650,0.704,0.286],[0.268,0.570,0.776],[0.711,0.423,0.561]])
    sm        /= np.linalg.norm(sm, axis=1)[:, np.newaxis]
    dab_od     = np.maximum(np.dot(od, np.linalg.inv(sm))[:,:,1], 0)

    h_map = np.zeros_like(dab_od)
    h_map[dab_od >= base_thresh] = 1
    h_map[dab_od >= 0.35]        = 2
    h_map[dab_od >= 0.50]        = 3

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].imshow(img);                    axes[0].set_title("Orijinal")
    axes[1].imshow(dab_od, cmap='hot');     axes[1].set_title("DAB OD Haritası")
    axes[2].imshow(h_map,  cmap='YlOrRd'); axes[2].set_title("H-Skor Kategorileri")
    for ax in axes: ax.axis('off')
    plt.tight_layout(); plt.show()

    print(f"OD min={dab_od.min():.3f}  max={dab_od.max():.3f}  mean={dab_od.mean():.3f}")
    print(f"Pozitif piksel oranı: {(dab_od>=base_thresh).mean()*100:.1f}%")

debug_image("kesitler/ornek.jpg")
```

### Eşik optimizasyonu (OD histogramı)

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_od_histogram(dab_od, thresholds=(0.22, 0.35, 0.50)):
    flat = dab_od.flatten()
    flat = flat[flat > 0.01]   # sıfır pikselleri dışla
    plt.figure(figsize=(8, 4))
    plt.hist(flat, bins=200, color='saddlebrown', alpha=0.7)
    colors = ['green', 'orange', 'red']
    labels = ['Negatif eşik (0+)', 'Zayıf eşik (1+→2+)', 'Güçlü eşik (2+→3+)']
    for t, c, l in zip(thresholds, colors, labels):
        plt.axvline(t, color=c, linestyle='--', label=l)
    plt.xlabel('DAB OD'); plt.ylabel('Piksel sayısı')
    plt.title('DAB OD Dağılımı — Eşik Kontrolü')
    plt.legend(); plt.tight_layout(); plt.show()
```

### Aggregated CSV oluşturma (Aşama 2→3 köprüsü)

```python
import pandas as pd, os

# Tüm SONUCLAR CSV'lerini aggregated dosyaya birleştir
base_path     = '/path/to/proje_klasoru'
sonuclar_dirs = ['BOLGE1_SONUCLAR', 'BOLGE2_SONUCLAR', 'BOLGE3_SONUCLAR']
all_dfs = []

for sdir in sonuclar_dirs:
    dpath = os.path.join(base_path, sdir)
    for f in os.listdir(dpath):
        if not f.endswith('.csv'): continue
        parts = f.replace('.csv','').split('_')
        region = parts[3].upper()
        name_u = f.upper()
        if   'GONAD_DISI'  in name_u: tissue, lectin = 'GONAD_DISI',  parts[-1].upper()
        elif 'GONAD_ERKEK' in name_u: tissue, lectin = 'GONAD_ERKEK', parts[-1].upper()
        else:                          tissue, lectin = parts[4].upper(), parts[5].upper()
        df = pd.read_csv(os.path.join(dpath, f), sep=';', decimal=',')
        df['Region'] = region; df['Tissue'] = tissue; df['Lectin'] = lectin
        all_dfs.append(df)

merged = pd.concat(all_dfs, ignore_index=True)
out    = os.path.join(base_path, 'istatistik_analiz', 'aggregated_hscore_data.csv')
merged.to_csv(out, index=False, sep=';', decimal=',')
print(f"Aggregated: {len(merged)} satır → {out}")
```

### Materyal & Metot Word belgesi oluşturma (`save_as_word.py`)

```python
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH

def create_methods_doc(output_path, lektin_listesi="WGA, SNA ve MAL", bolgeler="3 bölge"):
    doc = Document()
    doc.add_heading('MATERYAL VE METOT', 0).alignment = WD_ALIGN_PARAGRAPH.CENTER

    sections = [
        ("Lektin Histokimyası Görüntü Analizi ve H-Skor Hesaplamaları",
         f"Lektin ({lektin_listesi}) uygulaması sonrası elde edilen histolojik kesitlerdeki "
         "DAB boyanma yoğunluğu, bilgisayar destekli görüntü analiz yöntemleriyle değerlendirilmiştir."),
        ("Dijital Görüntü İşleme ve Renk Dekonvolüsyonu",
         "RGB görüntüler Python (v3.13) ve scikit-image kütüphanesi ile analiz edilmiştir. "
         "Ruifrok & Johnston (2001) renk dekonvolüsyonu algoritmasıyla DAB ve Hematoksilen "
         "spektral olarak ayrılmış; Lambert-Beer yasasına göre Optik Dansite (OD) değerleri hesaplanmıştır."),
        ("H-Skor Analizi",
         "Pozitif boyanma alt eşiği 0.22 OD olarak belirlenmiştir. Yoğunluk kategorileri: "
         "Zayıf (1+): 0.22–0.35, Orta (2+): 0.35–0.50, Güçlü (3+): ≥0.50. "
         "H-Skor = (1 × %Zayıf) + (2 × %Orta) + (3 × %Güçlü) formülüyle 0–300 ölçeğinde hesaplanmıştır."),
        ("İstatistiksel Analizler",
         f"Gruplar arası farklılıklar (bölgeler: {bolgeler}) Kruskal-Wallis H-Testi ile değerlendirilmiş; "
         "anlamlı farklılıklar Mann-Whitney U post-hoc testi ile ikili karşılaştırılmıştır (p < 0.05). "
         "Grafikler Mean ± SE olarak sunulmuş; * p<0.05, ** p<0.01, *** p<0.001.")
    ]
    for head, text in sections:
        doc.add_heading(head, level=1)
        p = doc.add_paragraph(text)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    doc.save(output_path)
    print(f"Word belgesi kaydedildi: {output_path}")

create_methods_doc("materyal_metot.docx",
                   lektin_listesi="WGA, SNA ve MAL",
                   bolgeler="Alsancak, Pasaport, Urla")
```

---

## Expected Outputs

| Dosya / Klasör | İçerik |
|----------------|--------|
| `lektin_hscore_sonuclari.csv` | Her görüntü için: H_Score, Mean_OD, %Weak, %Mod, %Strong |
| `analiz_cikti/hscore_*.png` | Orijinal + YlOrRd yoğunluk haritası (görüntü başına) |
| `istatistik_analiz/aggregated_hscore_data.csv` | Tüm bölge/doku/lektin verisi birleşik |
| `istatistik_analiz/plot_{DOKU}_{LEKTIN}_stars.png` | Bar + jitter + anlamlılık yıldızları |
| `istatistik_analiz/mean_sd_ozet.txt` | Region × Tissue × Lectin: Mean ± SD (n=) |
| `materyal_metot.docx` | Yayına hazır M&M bölümü (Word) |

**Her doku × lektin kombinasyonu için 1 grafik:** `GENEL×WGA`, `GENEL×SNA`, `GENEL×MAL`, `GONAD×WGA`, `GONAD×SNA`, `GONAD×MAL` → 6 PNG (veya alt grup varsa daha fazla).

---

## Troubleshooting

| Problem | Neden | Çözüm |
|---------|-------|-------|
| Tüm H-Skor değerleri 0 | `base_thresh` çok yüksek | OD histogramına bak; eşiği `0.15`'e düşür |
| Overlay haritası tamamen sarı | Güçlü eşik fazla düşük | `0.50` yerine `0.60`'a çıkar ve `debug_image()` ile kontrol et |
| CSV parse hatası | Dosya adı kuralı bozulmuş | `_{BOLGE}_{DOKU}_{LEKTIN}` formatını kontrol et; boşluk/Türkçe karakter olmamalı |
| `aggregated_hscore_data.csv` bulunamadı | `generate_summary.py` önce çalıştırılmamış | Önce `final_stats_and_plots.py` çalıştır → aggregated otomatik orada değil; aggregation recipe'sini çalıştır |
| Grafiklerde yıldız yok | Tüm p > 0.05 | Gerçekten anlamsız olabilir; veya `mannwhitneyu` için yeterli n yok (n≥3 gerekli) |
| Renk dekonvolüsyonu yanlış | Farklı protokol/boyama | ImageJ → Colour Deconvolution eklentisiyle kendi stain vektörlerini kalibre et, `stain_matrix`'i güncelle |
| `decimal=','` CSV Excel'de bozuk | Farklı bölge ayarı | `sep='\t', decimal='.'` dene veya Excel'de "Veriyi Al" ile nokta ondalık seç |

---

## References

- Ruifrok AC, Johnston DA (2001). "Quantification of histochemical staining by color deconvolution." *Anal Quant Cytol Histol* 23(4):291-299. — Stain matrix ve OD hesaplama yöntemi.
- Taylor CR, Bhatt DL (2011). "H-score methodology in IHC." *Histopathology* — H-Skor formülü standardı.
- [scikit-image color deconvolution](https://scikit-image.org/docs/stable/auto_examples/color_exposure/plot_ihc_color_separation.html) — HED dekonvolüsyon tutorial.
- Mann HB, Whitney DR (1947). "On a test of whether one of two random variables is stochastically larger than the other." *Ann Math Stat* 18:50-60.
- Kruskal WH, Wallis WA (1952). "Use of ranks in one-criterion variance analysis." *J Am Stat Assoc* 47:583-621.
