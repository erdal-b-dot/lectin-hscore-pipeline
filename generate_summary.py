import pandas as pd
import os

# Paths
base_path = '/Users/erdalbalcan/midye_lektin_olcum'
csv_path = os.path.join(base_path, 'istatistik_analiz', 'aggregated_hscore_data.csv')
output_path = os.path.join(base_path, 'istatistik_analiz', 'mean_sd_ozet.txt')

if not os.path.exists(csv_path):
    print("Hata: Agregre veri dosyası bulunamadı.")
    exit()

# Read data
df = pd.read_csv(csv_path, sep=';', decimal=',')

# Ensure Tissue_Group exists
df['Tissue_Group'] = df['Tissue'].apply(lambda x: 'GONAD' if 'GONAD' in x else x)

# Group by Region, Tissue_Group, and Lectin (to be precise)
summary = df.groupby(['Region', 'Tissue_Group', 'Lectin'])['H_Score'].agg(['mean', 'std', 'count']).reset_index()

# Format the output string
output_str = "H-SKOR ÖZET İSTATİSTİKLERİ (Mean ± SD)\n"
output_str += "="*50 + "\n\n"

regions = sorted(df['Region'].unique())
tissues = sorted(df['Tissue_Group'].unique())
lectins = sorted(df['Lectin'].unique())

for region in regions:
    output_str += f"BÖLGE: {region}\n"
    output_str += "-"*20 + "\n"
    for tissue in tissues:
        output_str += f"  Doku: {tissue}\n"
        for lectin in lectins:
            row = summary[(summary['Region'] == region) & 
                          (summary['Tissue_Group'] == tissue) & 
                          (summary['Lectin'] == lectin)]
            if not row.empty:
                m = row['mean'].values[0]
                s = row['std'].values[0]
                n = row['count'].values[0]
                # If std is NaN (n=1), show 0.00
                s_val = s if pd.notnull(s) else 0.0
                output_str += f"    - {lectin:3}: {m:6.2f} ± {s_val:5.2f} (n={n})\n"
        output_str += "\n"
    output_str += "\n"

# Save to file
with open(output_path, 'w', encoding='utf-8') as f:
    f.write(output_str)

print(f"Özet dosya kaydedildi: {output_path}")
print("\n" + output_str)
