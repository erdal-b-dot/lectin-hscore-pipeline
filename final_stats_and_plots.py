import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations
import numpy as np
import warnings

warnings.filterwarnings('ignore')

# ── Eşik (threshold) ──────────────────────────────────────────────────────
# 0.22 → 0.28 güncellendi (2026-06-21).
# Gerekçe: threshold_comparison.py ile yapılan n=421 karşılaştırmasında
# 0.22 OD'nin lektin histokimyasında arka plan sinyalini pozitif saydığı
# gösterildi. 0.28, WGA biyolojik sinyalini korurken SNA/MAL'daki
# non-spesifik boyanmayı ~%45-55 azaltıyor.
THRESHOLD = 0.28

# ── Veri kaynağı ──────────────────────────────────────────────────────────
# threshold_comparison.csv: tüm lektin görüntüleri (n=421) 0.22/0.28/0.30
# eşikleriyle analiz edilmiş. THRESHOLD=0.28 satırları kullanılır.
BASE_PATH = os.path.dirname(os.path.abspath(__file__))
COMPARISON_CSV = os.path.join(BASE_PATH, 'threshold_comparison.csv')
output_dir = os.path.join(BASE_PATH, 'istatistik_analiz')

print(f"Veri okunuyor (eşik={THRESHOLD})...")

if not os.path.exists(COMPARISON_CSV):
    print(f"HATA: {COMPARISON_CSV} bulunamadı.")
    print("Önce threshold_comparison.py çalıştırın.")
    exit(1)

raw = pd.read_csv(COMPARISON_CSV, sep=';', decimal=',')
full_df = raw[raw['Threshold'] == THRESHOLD].copy()
full_df = full_df.rename(columns={'H_Score': 'H_Score', 'Tissue': 'Tissue_Group'})
full_df['Tissue_Group'] = full_df['Tissue_Group'].apply(
    lambda x: 'GONAD' if 'GONAD' in str(x).upper() else x
)

print(f"Toplam kayıt: {len(full_df)} | Bölgeler: {sorted(full_df['Region'].unique())} | "
      f"Lektinler: {sorted(full_df['Lectin'].unique())}")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def get_star(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return None

def annotate_significance(ax, sub_df, group_col, val_col):
    groups = sub_df[group_col].unique()
    if len(groups) < 2: return
    
    # Get group positions on x-axis
    group_names = [t.get_text() for t in ax.get_xticklabels()]
    y_max = sub_df[val_col].max()
    y_shift = y_max * 0.1
    current_y = y_max + y_shift

    # Compare pairs
    for g1, g2 in combinations(group_names, 2):
        d1 = sub_df[sub_df[group_col] == g1][val_col].dropna()
        d2 = sub_df[sub_df[group_col] == g2][val_col].dropna()
        
        if len(d1) > 0 and len(d2) > 0:
            _, p = mannwhitneyu(d1, d2)
            star = get_star(p)
            
            if star:
                idx1 = group_names.index(g1)
                idx2 = group_names.index(g2)
                
                # Draw line and star
                lx, rx = idx1, idx2
                ly, ry = current_y, current_y
                ax.plot([lx, lx, rx, rx], [ly - y_shift*0.2, ly, ry, ry - y_shift*0.2], lw=1.5, c='black')
                ax.text((lx + rx) / 2, ly + y_shift*0.1, star, ha='center', va='bottom', color='black', fontsize=12, fontweight='bold')
                
                current_y += y_shift * 1.5 # Move up for next comparison

print("Generating plots with significance stars...")
for tissue in full_df['Tissue_Group'].unique():
    for lectin in full_df['Lectin'].unique():
        sub_df = full_df[(full_df['Tissue_Group'] == tissue) & (full_df['Lectin'] == lectin)]
        
        if sub_df.empty or sub_df['Region'].nunique() < 2:
            continue
            
        plt.figure(figsize=(10, 8))
        sns.set_style("whitegrid")
        
        # Determine order for consistent indexing
        order = sorted(sub_df['Region'].unique())
        
        ax = sns.barplot(data=sub_df, x='Region', y='H_Score', order=order, capsize=.1, errorbar='se', palette='muted')
        sns.stripplot(data=sub_df, x='Region', y='H_Score', order=order, color='black', alpha=0.4, jitter=True)
        
        plt.title(f'H-Score: {tissue} - {lectin}', fontsize=14)
        plt.ylabel('H-Score (Mean ± SE)', fontsize=12)
        
        # Add stars
        annotate_significance(ax, sub_df, 'Region', 'H_Score')
        
        # Set y limit higher to accommodate stars
        y_lim = ax.get_ylim()
        ax.set_ylim(y_lim[0], y_lim[1] * 1.3)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'plot_{tissue}_{lectin}_stars.png'))
        plt.close()
        print(f"Saved plot with stars: {tissue} - {lectin}")

print(f"\nCompleted! Check: {output_dir}")
