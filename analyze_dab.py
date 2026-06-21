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
    inverse_matrix = np.linalg.inv(stain_matrix)
    deconvolved = np.dot(od, inverse_matrix)
    return deconvolved

def analyze_h_score(image_path, base_thresh=0.28):
    try:
        img = io.imread(image_path)
        if img.shape[-1] == 4:
            img = color.rgba2rgb(img)

        deconvolved = imagej_color_deconvolution(img)
        dab_od = np.maximum(deconvolved[:, :, 1], 0)

        # Yoğunluk Sınıflandırması (Thresholds)
        # base_thresh altı negatif (arka plan) kabul ediliyor
        # 0.28 seçildi: eşik karşılaştırma analizinde (n=421) 0.22'nin
        # arka plan sinyalini pozitif saydığı gösterildi (threshold_comparison.py)
        weak_mask = (dab_od >= base_thresh) & (dab_od < 0.35)
        mod_mask = (dab_od >= 0.35) & (dab_od < 0.50)
        strong_mask = (dab_od >= 0.50)
        
        total_pos_mask = dab_od >= base_thresh
        
        # Alan Yüzdeleri
        total_pixels = dab_od.size
        p_weak = (np.sum(weak_mask) / total_pixels) * 100
        p_mod = (np.sum(mod_mask) / total_pixels) * 100
        p_strong = (np.sum(strong_mask) / total_pixels) * 100
        p_total = (np.sum(total_pos_mask) / total_pixels) * 100
        
        # H-Score Hesaplama (0-300)
        h_score = (1 * p_weak) + (2 * p_mod) + (3 * p_strong)
        
        # Ortalama OD (Sadece pozitif alanlarda)
        mean_od = np.mean(dab_od[total_pos_mask]) if np.any(total_pos_mask) else 0
        
        return {
            'mean_od': mean_od,
            'p_total': p_total,
            'p_weak': p_weak,
            'p_mod': p_mod,
            'p_strong': p_strong,
            'h_score': h_score,
            'dab_od': dab_od
        }

    except Exception as e:
        print(f"Hata: {image_path} -> {e}")
        return None

def main():
    input_folder = 'kesitler'
    output_plots = 'analiz_cikti'
    MY_THRESHOLD = 0.28  # güncellendi: 0.22 → 0.28 (threshold_comparison.py, 2026-06-21)
    
    if not os.path.exists(output_plots):
        os.makedirs(output_plots)

    results = []
    files = [f for f in os.listdir(input_folder) if f.lower().endswith(('.tif', '.png', '.jpg', '.jpeg'))]
    
    print(f"H-Score Analizi Başlıyor (Eşik: {MY_THRESHOLD})...")

    for filename in files:
        path = os.path.join(input_folder, filename)
        data = analyze_h_score(path, MY_THRESHOLD)
        
        if data:
            results.append({
                'File_Name': filename,
                'Mean_OD': round(data['mean_od'], 4),
                'Total_Positive_Area_Pct': round(data['p_total'], 2),
                'Weak_1_Plus_Area_Pct': round(data['p_weak'], 2),
                'Moderate_2_Plus_Area_Pct': round(data['p_mod'], 2),
                'Strong_3_Plus_Area_Pct': round(data['p_strong'], 2),
                'H_Score': round(data['h_score'], 2)
            })
            
            # Görselleştirme (H-Score Haritası)
            plt.figure(figsize=(10, 5))
            plt.subplot(1, 2, 1)
            plt.imshow(io.imread(path))
            plt.title('Orijinal')
            plt.axis('off')
            
            plt.subplot(1, 2, 2)
            # Yoğunluk haritası oluştur (0: Arka plan, 1: Zayıf, 2: Orta, 3: Güçlü)
            h_map = np.zeros_like(data['dab_od'])
            h_map[data['dab_od'] >= MY_THRESHOLD] = 1
            h_map[data['dab_od'] >= 0.35] = 2
            h_map[data['dab_od'] >= 0.50] = 3
            plt.imshow(h_map, cmap='YlOrRd')
            plt.title(f'H-Score: {round(data['h_score'], 1)}')
            plt.axis('off')
            
            plt.savefig(os.path.join(output_plots, f'hscore_{filename}.png'))
            plt.close()
            print(f"{filename} işlendi: H-Score = {round(data['h_score'], 2)}")

    # Sonuçları Kaydet (ASCII uyumlu ve Excel dostu)
    df = pd.DataFrame(results)
    df.to_csv('lektin_hscore_sonuclari.csv', index=False, sep=';', decimal=',')
    print("\nTamamlandı! 'lektin_hscore_sonuclari.csv' dosyasını Excel ile açabilirsiniz.")

if __name__ == "__main__":
    main()
