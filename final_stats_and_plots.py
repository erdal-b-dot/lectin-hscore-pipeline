import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations
import numpy as np
import warnings

warnings.filterwarnings('ignore')

# Directories
base_path = '/Users/erdalbalcan/midye_lektin_olcum'
sonuclar_dirs = ['ALSANCAK_SONUCLAR', 'PASAPORT_SONUCLAR', 'URLA_SONUCLAR']
output_dir = os.path.join(base_path, 'istatistik_analiz')

all_data = []

print("Reading files...")
for sdir in sonuclar_dirs:
    dir_path = os.path.join(base_path, sdir)
    if not os.path.exists(dir_path):
        continue
    
    files = [f for f in os.listdir(dir_path) if f.endswith('.csv')]
    for file in files:
        name_no_ext = file.replace('.csv', '')
        parts = name_no_ext.split('_')
        region = parts[3].upper()
        
        if 'GONAD_DISI' in name_no_ext.upper():
            tissue = 'GONAD_DISI'
            lectin = parts[-1].upper()
        elif 'GONAD_ERKEK' in name_no_ext.upper():
            tissue = 'GONAD_ERKEK'
            lectin = parts[-1].upper()
        else:
            tissue = parts[4].upper()
            lectin = parts[5].upper()
            
        file_path = os.path.join(dir_path, file)
        try:
            df = pd.read_csv(file_path, sep=';', decimal=',')
            df['Region'] = region
            df['Tissue'] = tissue
            df['Lectin'] = lectin
            all_data.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}")

if not all_data:
    print("No data found!")
    exit()

full_df = pd.concat(all_data, ignore_index=True)
full_df['Tissue_Group'] = full_df['Tissue'].apply(lambda x: 'GONAD' if 'GONAD' in x else x)

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
