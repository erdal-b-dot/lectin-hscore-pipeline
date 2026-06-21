"""
Lektin H-Skoru İstatistiksel Analiz ve Grafik Scripti
Eşik: 0.28 OD (threshold_comparison.csv kaynağı)

Bölgeler : ALSANCAK · PASAPORT · URLA
Dokular  : Hepatopankreas (GENEL DOKU) · Gonad
Lektinler: WGA · SNA · MAL

Testler:
  - Kruskal-Wallis (3 bölge arası genel fark)
  - Mann-Whitney U post-hoc (tüm ikili kombinasyonlar)
  - Bonferroni düzeltmesi (3 karşılaştırma için α=0.05/3≈0.0167)

Çıktılar (istatistik_analiz_0.28/ klasörü):
  - Her lektin × doku için bar+jitter grafik (PNG, 300 DPI)
  - Tüm lektinleri bir arada gösteren panel grafik
  - istatistik_ozet.csv  — tüm test sonuçları
  - ozet_tablo.csv       — ortalama ± SD özet
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations
import warnings

warnings.filterwarnings('ignore')

# ── Sabitler ────────────────────────────────────────────────────────────────
THRESHOLD    = 0.28
BASE_PATH    = os.path.dirname(os.path.abspath(__file__))
INPUT_CSV    = os.path.join(BASE_PATH, 'threshold_comparison.csv')
OUTPUT_DIR   = os.path.join(BASE_PATH, 'istatistik_analiz_0.28')

REGION_ORDER  = ['ALSANCAK', 'PASAPORT', 'URLA']
LECTIN_ORDER  = ['WGA', 'SNA', 'MAL']
TISSUE_MAP    = {'GENEL DOKU': 'Hepatopankreas', 'GONAD': 'Gonad'}
TISSUE_ORDER  = ['Hepatopankreas', 'Gonad']

# Renk paleti (bölge bazlı)
REGION_COLORS = {
    'ALSANCAK': '#2E86AB',
    'PASAPORT': '#E84855',
    'URLA'    : '#3BB273',
}

ALPHA_BONF = 0.05 / 3  # Bonferroni düzeltmeli eşik (3 ikili karşılaştırma)

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Veri yükle ───────────────────────────────────────────────────────────────
print(f"Veri okunuyor (eşik={THRESHOLD})...")
raw = pd.read_csv(INPUT_CSV, sep=';', decimal=',')
df  = raw[raw['Threshold'] == THRESHOLD].copy()
df['Tissue'] = df['Tissue'].map(TISSUE_MAP).fillna(df['Tissue'])

print(f"  Toplam kayıt : {len(df)}")
print(f"  Bölgeler     : {sorted(df['Region'].unique())}")
print(f"  Dokular      : {sorted(df['Tissue'].unique())}")
print(f"  Lektinler    : {sorted(df['Lectin'].unique())}\n")


# ── İstatistik fonksiyonları ─────────────────────────────────────────────────
def star(p, bonf=False):
    thr = ALPHA_BONF if bonf else 0.05
    if   p < 0.001:       return '***'
    elif p < 0.01:        return '**'
    elif p < thr:         return '*'
    return 'ns'

def run_stats(sub, group_col='Region', val_col='H_Score'):
    groups = [sub[sub[group_col]==g][val_col].dropna().values
              for g in REGION_ORDER if g in sub[group_col].values]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return None, []

    # Kruskal-Wallis
    try:
        kw_stat, kw_p = kruskal(*groups)
    except Exception:
        kw_stat, kw_p = np.nan, np.nan

    # Mann-Whitney U post-hoc (Bonferroni düzeltmeli)
    pairs = []
    regions_present = [r for r in REGION_ORDER if r in sub[group_col].values]
    for r1, r2 in combinations(regions_present, 2):
        d1 = sub[sub[group_col]==r1][val_col].dropna().values
        d2 = sub[sub[group_col]==r2][val_col].dropna().values
        if len(d1) > 0 and len(d2) > 0:
            _, p = mannwhitneyu(d1, d2, alternative='two-sided')
            pairs.append({'group1': r1, 'group2': r2, 'p_raw': p,
                          'p_bonf': min(p * 3, 1.0),
                          'star': star(min(p*3, 1.0), bonf=True)})

    kw_result = {'kw_stat': round(kw_stat, 3), 'kw_p': round(kw_p, 4)}
    return kw_result, pairs


# ── Grafik çizim fonksiyonu ──────────────────────────────────────────────────
def draw_significance(ax, pair_results, region_order, y_max, y_range):
    """Anlamlı çiftlere köprü + yıldız çizer."""
    sig_pairs = [p for p in pair_results if p['star'] != 'ns']
    if not sig_pairs:
        return y_max

    x_pos = {r: i for i, r in enumerate(region_order)}
    step   = y_range * 0.13
    y_cur  = y_max + y_range * 0.08

    for pair in sig_pairs:
        x1 = x_pos.get(pair['group1'])
        x2 = x_pos.get(pair['group2'])
        if x1 is None or x2 is None:
            continue
        # Köprü çizgisi
        ax.plot([x1, x1, x2, x2],
                [y_cur, y_cur + step*0.25, y_cur + step*0.25, y_cur],
                lw=1.2, c='#333333')
        ax.text((x1+x2)/2, y_cur + step*0.3, pair['star'],
                ha='center', va='bottom', fontsize=11,
                fontweight='bold', color='#333333')
        y_cur += step

    return y_cur + step * 0.5


def plot_lectin_tissue(sub, lectin, tissue, ax, show_legend=True):
    """Tek panel: bir lektin × doku kombinasyonu için bar+jitter grafik."""
    regions = [r for r in REGION_ORDER if r in sub['Region'].values]

    # --- Bar (ortalama ± SE) ---
    means = [sub[sub['Region']==r]['H_Score'].mean() for r in regions]
    sems  = [sub[sub['Region']==r]['H_Score'].sem()  for r in regions]
    x     = np.arange(len(regions))
    bars  = ax.bar(x, means, yerr=sems, capsize=5, width=0.55,
                   color=[REGION_COLORS[r] for r in regions],
                   edgecolor='white', linewidth=0.8,
                   error_kw=dict(elinewidth=1.2, ecolor='#444'))

    # --- Jitter (bireysel veri noktaları) ---
    rng = np.random.default_rng(42)
    for i, region in enumerate(regions):
        vals = sub[sub['Region']==region]['H_Score'].dropna().values
        jitter = rng.uniform(-0.18, 0.18, size=len(vals))
        ax.scatter(i + jitter, vals,
                   color=REGION_COLORS[region],
                   alpha=0.55, s=18, zorder=3,
                   edgecolors='white', linewidths=0.4)

    # --- İstatistik köprüleri ---
    kw_res, pairs = run_stats(sub)
    y_max   = sub['H_Score'].max()
    y_range = max(y_max * 0.5, 10)
    top_y   = draw_significance(ax, pairs, regions, y_max, y_range)

    # KW p değerini başlığa ekle
    kw_txt = ''
    if kw_res:
        p_val = kw_res['kw_p']
        kw_txt = f'  (KW p={p_val:.3f})' if p_val >= 0.001 else '  (KW p<0.001)'

    ax.set_title(f'{lectin} — {tissue}{kw_txt}',
                 fontsize=10, fontweight='bold', pad=6)
    ax.set_xticks(x)
    ax.set_xticklabels(regions, fontsize=9)
    ax.set_ylabel('H-Skoru (Ort ± SE)', fontsize=9)
    ax.set_ylim(0, max(top_y * 1.1, y_max * 1.45))
    ax.spines[['top','right']].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.4)
    ax.set_axisbelow(True)

    if show_legend:
        legend_els = [mpatches.Patch(color=REGION_COLORS[r], label=r)
                      for r in regions]
        ax.legend(handles=legend_els, fontsize=8, framealpha=0.7)

    return kw_res, pairs


# ── Ana analiz döngüsü ────────────────────────────────────────────────────────
all_stats = []
print("Grafikler üretiliyor...\n")

for tissue in TISSUE_ORDER:
    for lectin in LECTIN_ORDER:
        sub = df[(df['Tissue']==tissue) & (df['Lectin']==lectin)]
        if sub.empty:
            print(f"  [ATLA] {tissue} × {lectin}: veri yok")
            continue

        fig, ax = plt.subplots(figsize=(6, 5))
        kw_res, pairs = plot_lectin_tissue(sub, lectin, tissue, ax)
        plt.tight_layout()

        fname = f'{tissue.replace(" ","_")}_{lectin}_0.28.png'
        fig.savefig(os.path.join(OUTPUT_DIR, fname), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Kaydedildi: {fname}")

        # İstatistik kayıt
        for p in pairs:
            all_stats.append({
                'Tissue' : tissue,
                'Lectin' : lectin,
                'KW_stat': kw_res['kw_stat'] if kw_res else np.nan,
                'KW_p'   : kw_res['kw_p']   if kw_res else np.nan,
                'Group1' : p['group1'],
                'Group2' : p['group2'],
                'p_raw'  : round(p['p_raw'], 4),
                'p_bonf' : round(p['p_bonf'], 4),
                'Star'   : p['star'],
            })

# ── Panel grafik (tüm lektin × doku kombinasyonları) ────────────────────────
print("\nPanel grafik üretiliyor...")
fig, axes = plt.subplots(len(LECTIN_ORDER), len(TISSUE_ORDER),
                         figsize=(5*len(TISSUE_ORDER), 4.5*len(LECTIN_ORDER)),
                         squeeze=False)

for li, lectin in enumerate(LECTIN_ORDER):
    for ti, tissue in enumerate(TISSUE_ORDER):
        ax = axes[li][ti]
        sub = df[(df['Tissue']==tissue) & (df['Lectin']==lectin)]
        if sub.empty:
            ax.set_visible(False)
            continue
        plot_lectin_tissue(sub, lectin, tissue, ax, show_legend=(ti==0 and li==0))

fig.suptitle('Lektin H-Skoru — Bölge × Doku Karşılaştırması\n'
             f'(Eşik={THRESHOLD} OD, Mann-Whitney U + Bonferroni)',
             fontsize=13, fontweight='bold', y=1.01)

legend_els = [mpatches.Patch(color=REGION_COLORS[r], label=r) for r in REGION_ORDER]
fig.legend(handles=legend_els, loc='upper right', fontsize=9,
           bbox_to_anchor=(1.0, 1.0), framealpha=0.8)

plt.tight_layout()
panel_path = os.path.join(OUTPUT_DIR, 'PANEL_tum_lektinler_0.28.png')
fig.savefig(panel_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Panel kaydedildi: PANEL_tum_lektinler_0.28.png")

# ── Özet tablo (ortalama ± SD, n) ───────────────────────────────────────────
rows = []
for tissue in TISSUE_ORDER:
    for lectin in LECTIN_ORDER:
        for region in REGION_ORDER:
            sub = df[(df['Tissue']==tissue) & (df['Lectin']==lectin) &
                     (df['Region']==region)]['H_Score'].dropna()
            if len(sub) == 0:
                continue
            rows.append({
                'Tissue': tissue, 'Lectin': lectin, 'Region': region,
                'n'     : len(sub),
                'Mean'  : round(sub.mean(), 2),
                'SD'    : round(sub.std(),  2),
                'SE'    : round(sub.sem(),  2),
                'Min'   : round(sub.min(),  2),
                'Max'   : round(sub.max(),  2),
            })

summary_df = pd.DataFrame(rows)
summary_df.to_csv(os.path.join(OUTPUT_DIR, 'ozet_tablo.csv'),
                  index=False, sep=';', decimal=',')

# ── İstatistik özet CSV ──────────────────────────────────────────────────────
stats_df = pd.DataFrame(all_stats)
stats_df.to_csv(os.path.join(OUTPUT_DIR, 'istatistik_ozet.csv'),
                index=False, sep=';', decimal=',')

# ── Terminale özet yazdır ────────────────────────────────────────────────────
print("\n" + "="*60)
print("ÖZET TABLO (Ort ± SD)")
print("="*60)
print(summary_df.to_string(index=False))

print("\n" + "="*60)
print("İSTATİSTİK SONUÇLARI")
print("="*60)
if not stats_df.empty:
    print(stats_df.to_string(index=False))

print(f"\n✓ Tüm çıktılar: {OUTPUT_DIR}")
print(f"  • {len(LECTIN_ORDER)*len(TISSUE_ORDER)} bireysel grafik")
print(f"  • 1 panel grafik")
print(f"  • istatistik_ozet.csv + ozet_tablo.csv")
