"""
Eşik Karşılaştırma: 0.22 vs 0.28 vs 0.30
Tüm lektin görüntüleri (WGA, SNA, MAL) üzerinde çalışır.
Bölge / Doku / Lektin bilgisini klasör adından otomatik çeker.
"""

import os
import re
import unicodedata
import numpy as np
import pandas as pd
from skimage import io, color
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.dirname(os.path.abspath(__file__))
THRESHOLDS = [0.22, 0.28, 0.30]
REGIONS = ['ALSANCAK', 'PASAPORT', 'URLA']
IMG_EXTS = ('.jpg', '.jpeg', '.tif', '.tiff', '.png')


def color_deconvolution(rgb_img):
    rgb_float = rgb_img.astype(float) + 1.0
    od = -np.log10(rgb_float / 255.0)
    stain_matrix = np.array([
        [0.650, 0.704, 0.286],
        [0.268, 0.570, 0.776],
        [0.711, 0.423, 0.561]
    ])
    stain_matrix /= np.linalg.norm(stain_matrix, axis=1)[:, np.newaxis]
    inv = np.linalg.inv(stain_matrix)
    return np.dot(od, inv)


def compute_hscore(image_path, base_thresh):
    try:
        img = io.imread(image_path)
        if img.ndim == 2:
            img = np.stack([img]*3, axis=-1)
        if img.shape[-1] == 4:
            img = color.rgba2rgb(img)
            img = (img * 255).astype(np.uint8)

        dab_od = np.maximum(color_deconvolution(img)[:, :, 1], 0)

        total_pixels = dab_od.size
        weak   = (dab_od >= base_thresh) & (dab_od < 0.35)
        mod    = (dab_od >= 0.35)        & (dab_od < 0.50)
        strong = (dab_od >= 0.50)
        pos    = dab_od >= base_thresh

        p_weak   = np.sum(weak)   / total_pixels * 100
        p_mod    = np.sum(mod)    / total_pixels * 100
        p_strong = np.sum(strong) / total_pixels * 100
        p_total  = np.sum(pos)    / total_pixels * 100
        h_score  = 1*p_weak + 2*p_mod + 3*p_strong
        mean_od  = float(np.mean(dab_od[pos])) if np.any(pos) else 0.0

        return dict(p_weak=p_weak, p_mod=p_mod, p_strong=p_strong,
                    p_total=p_total, h_score=h_score, mean_od=mean_od,
                    dab_od=dab_od)
    except Exception as e:
        print(f"  HATA: {os.path.basename(image_path)} -> {e}")
        return None


def parse_folder(folder_name):
    """'GENEL DOKU (LEKTòN-WGA)' -> ('GENEL DOKU', 'WGA')"""
    # macOS HFS+ NFD normalization: 'ò' stored as o + combining accent → NFC
    folder_nfc = unicodedata.normalize('NFC', folder_name)
    m = re.search(r'LEKT.N[-_](\w+)', folder_nfc, re.IGNORECASE)
    lectin = m.group(1).upper() if m else None

    if 'GONAD' in folder_name.upper():
        tissue = 'GONAD'
    elif 'GENEL' in folder_name.upper():
        tissue = 'GENEL DOKU'
    else:
        tissue = folder_name

    return tissue, lectin


def collect_images():
    """Lektin boyamalı tüm görüntüleri toplar."""
    records = []
    for region in REGIONS:
        region_dir = os.path.join(BASE, region)
        if not os.path.isdir(region_dir):
            continue
        for folder in os.listdir(region_dir):
            folder_path = os.path.join(region_dir, folder)
            if not os.path.isdir(folder_path):
                continue
            tissue, lectin = parse_folder(folder)
            if lectin is None:
                continue  # HE, AB vb. atla
            for f in os.listdir(folder_path):
                if f.lower().endswith(IMG_EXTS) and not f.startswith('.'):
                    records.append(dict(
                        region=region,
                        tissue=tissue,
                        lectin=lectin,
                        path=os.path.join(folder_path, f),
                        filename=f
                    ))
    return records


def run_comparison():
    images = collect_images()
    print(f"Toplam lektin görüntüsü: {len(images)}")
    if not images:
        print("Görüntü bulunamadı!")
        return

    rows = []
    for i, rec in enumerate(images, 1):
        print(f"[{i}/{len(images)}] {rec['region']} | {rec['tissue']} | {rec['lectin']} | {rec['filename']}")
        for thr in THRESHOLDS:
            res = compute_hscore(rec['path'], thr)
            if res:
                rows.append(dict(
                    Region=rec['region'],
                    Tissue=rec['tissue'],
                    Lectin=rec['lectin'],
                    File=rec['filename'],
                    Threshold=thr,
                    H_Score=round(res['h_score'], 2),
                    Total_Pos_Pct=round(res['p_total'], 2),
                    Weak_Pct=round(res['p_weak'], 2),
                    Mod_Pct=round(res['p_mod'], 2),
                    Strong_Pct=round(res['p_strong'], 2),
                    Mean_OD=round(res['mean_od'], 4)
                ))

    df = pd.DataFrame(rows)

    # ── CSV kaydet ──
    out_csv = os.path.join(BASE, 'threshold_comparison.csv')
    df.to_csv(out_csv, index=False, sep=';', decimal=',')
    print(f"\nCSV kaydedildi: {out_csv}")

    # ── Özet tablo ──
    summary = (df.groupby(['Threshold', 'Lectin', 'Tissue'])['H_Score']
                 .agg(['mean', 'std', 'count'])
                 .round(2)
                 .reset_index())
    summary.columns = ['Threshold', 'Lectin', 'Tissue', 'Mean_H', 'SD_H', 'N']
    out_sum = os.path.join(BASE, 'threshold_comparison_summary.csv')
    summary.to_csv(out_sum, index=False, sep=';', decimal=',')
    print(f"Özet CSV kaydedildi: {out_sum}")

    # ── Karşılaştırma Grafikleri ──
    plot_comparison(df, summary)
    print("Grafikler kaydedildi.")
    print("\n=== ÖZET (Tüm görüntüler ortalaması) ===")
    print(summary.to_string(index=False))


def plot_comparison(df, summary):
    lectins = sorted(df['Lectin'].unique())
    tissues = sorted(df['Tissue'].unique())
    colors  = {0.22: '#e74c3c', 0.28: '#e67e22', 0.30: '#2ecc71'}
    labels  = {0.22: 'Eşik 0.22 (mevcut)', 0.28: 'Eşik 0.28', 0.30: 'Eşik 0.30'}

    # ── FIG 1: Lektin × Doku ortalama H-skoru ──
    fig, axes = plt.subplots(len(tissues), len(lectins),
                             figsize=(5*len(lectins), 4*len(tissues)),
                             sharey=False)
    if len(tissues) == 1:
        axes = [axes]
    if len(lectins) == 1:
        axes = [[ax] for ax in axes]

    for ti, tissue in enumerate(tissues):
        for li, lectin in enumerate(lectins):
            ax = axes[ti][li]
            sub = summary[(summary['Lectin']==lectin) & (summary['Tissue']==tissue)]
            if sub.empty:
                ax.set_visible(False)
                continue
            x = np.arange(len(THRESHOLDS))
            bars = ax.bar(x, sub.set_index('Threshold').reindex(THRESHOLDS)['Mean_H'],
                          color=[colors[t] for t in THRESHOLDS],
                          yerr=sub.set_index('Threshold').reindex(THRESHOLDS)['SD_H'],
                          capsize=4, width=0.55, edgecolor='white')
            ax.set_xticks(x)
            ax.set_xticklabels([str(t) for t in THRESHOLDS])
            ax.set_title(f'{lectin} — {tissue}', fontsize=10, fontweight='bold')
            ax.set_xlabel('Alt Eşik (OD)')
            ax.set_ylabel('Ortalama H-Skoru')
            ax.set_ylim(0, max(sub['Mean_H'].max()*1.35, 10))
            for bar, thr in zip(bars, THRESHOLDS):
                val = sub[sub['Threshold']==thr]['Mean_H'].values
                if val.size:
                    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
                            f'{val[0]:.1f}', ha='center', va='bottom', fontsize=8)

    patches = [mpatches.Patch(color=colors[t], label=labels[t]) for t in THRESHOLDS]
    fig.legend(handles=patches, loc='upper center', ncol=3, fontsize=9,
               bbox_to_anchor=(0.5, 1.02))
    fig.suptitle('Eşik Karşılaştırması: Ortalama H-Skoru\n(hata çubukları = ±SD)',
                 y=1.05, fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(BASE, 'threshold_comparison_hscore.png'),
                dpi=150, bbox_inches='tight')
    plt.close()

    # ── FIG 2: Bölge bazlı H-skoru dağılımı (boxplot) ──
    fig2, axes2 = plt.subplots(1, len(lectins), figsize=(5*len(lectins), 5), sharey=False)
    if len(lectins) == 1:
        axes2 = [axes2]

    for li, lectin in enumerate(lectins):
        ax = axes2[li]
        sub = df[df['Lectin'] == lectin]
        data_by_thresh = [sub[sub['Threshold']==t]['H_Score'].values for t in THRESHOLDS]
        bp = ax.boxplot(data_by_thresh, patch_artist=True,
                        labels=[str(t) for t in THRESHOLDS])
        for patch, thr in zip(bp['boxes'], THRESHOLDS):
            patch.set_facecolor(colors[thr])
            patch.set_alpha(0.7)
        ax.set_title(f'{lectin}', fontsize=11, fontweight='bold')
        ax.set_xlabel('Alt Eşik (OD)')
        ax.set_ylabel('H-Skoru')

    fig2.suptitle('Eşik Karşılaştırması: H-Skoru Dağılımı (Tüm görüntüler)',
                  fontsize=12, fontweight='bold')
    plt.tight_layout()
    fig2.savefig(os.path.join(BASE, 'threshold_comparison_boxplot.png'),
                 dpi=150, bbox_inches='tight')
    plt.close()

    # ── FIG 3: Pozitif alan % değişimi ──
    fig3, axes3 = plt.subplots(1, len(lectins), figsize=(5*len(lectins), 4))
    if len(lectins) == 1:
        axes3 = [axes3]

    for li, lectin in enumerate(lectins):
        ax = axes3[li]
        sub = df[df['Lectin'] == lectin]
        means = [sub[sub['Threshold']==t]['Total_Pos_Pct'].mean() for t in THRESHOLDS]
        ax.plot([str(t) for t in THRESHOLDS], means, 'o-', color='#2c3e50', linewidth=2)
        for x, y in enumerate(means):
            ax.annotate(f'{y:.1f}%', (x, y), textcoords="offset points",
                        xytext=(0, 8), ha='center', fontsize=9)
        ax.set_title(f'{lectin} — Pozitif Alan %', fontsize=10, fontweight='bold')
        ax.set_xlabel('Alt Eşik (OD)')
        ax.set_ylabel('Ortalama Pozitif Alan (%)')
        ax.set_ylim(0, max(means)*1.4 if max(means) > 0 else 10)

    fig3.suptitle('Eşik Yükseldikçe Pozitif Alan Değişimi',
                  fontsize=12, fontweight='bold')
    plt.tight_layout()
    fig3.savefig(os.path.join(BASE, 'threshold_comparison_posarea.png'),
                 dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    run_comparison()
