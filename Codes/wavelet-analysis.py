import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pywt
from scipy.stats import pearsonr
import os

def load_data(file_path):
    return pd.read_csv(file_path, sep=' ', parse_dates=['Time'])

def load_river_data(river_files):
    river_data = {}
    for i, file in enumerate(river_files, 1):
        river_data[f'R{i}'] = pd.read_csv(file, parse_dates=['Date'])
    return river_data

def load_tide_data(tide_file):
    return pd.read_csv(tide_file, parse_dates=['Time'])

def wavelet_analysis(data, location, river_data, tide_data):
    os.makedirs(f'wavelet_analysis_{location}', exist_ok=True)
    
    # Merge with river and tide data
    for river, discharge in river_data.items():
        data = data.merge(discharge, left_on='Time', right_on='Date', how='left')
        data = data.rename(columns={'Discharge': f'{river}_Discharge'})
    data = data.merge(tide_data, on='Time', how='left')
    
    variables = ['SAL', 'TEMP', 'V_x', 'V_y']
    depths = [f'{i*10}_{(i+1)*10}' for i in range(10)]
    
    # Wavelet parameters
    wavelet = 'cmor1.5-1.0'
    scales = np.arange(1, 129)
    
    for var in variables:
        # Wavelet transform for each depth
        for depth in depths:
            coef, freqs = pywt.cwt(data[f'{var}_{depth}'].dropna(), scales, wavelet)
            
            plt.figure(figsize=(12, 8))
            plt.imshow(np.abs(coef), extent=[0, len(data), 1, 128], aspect='auto',
                       cmap='viridis', interpolation='nearest')
            plt.colorbar(label='Magnitude')
            plt.title(f'Wavelet Transform of {var} at {depth.replace("_", "-")}% depth - {location}')
            plt.ylabel('Scale')
            plt.xlabel('Time')
            plt.savefig(f'wavelet_analysis_{location}/Wavelet_{var}_{depth}.png')
            plt.close()
        
        # Cross-wavelet transform between surface and bottom layers
        coef_surface, _ = pywt.cwt(data[f'{var}_0_10'].dropna(), scales, wavelet)
        coef_bottom, _ = pywt.cwt(data[f'{var}_90_100'].dropna(), scales, wavelet)
        cross_wavelet = coef_surface * np.conj(coef_bottom)
        
        plt.figure(figsize=(12, 8))
        plt.imshow(np.abs(cross_wavelet), extent=[0, len(data), 1, 128], aspect='auto',
                   cmap='viridis', interpolation='nearest')
        plt.colorbar(label='Magnitude')
        plt.title(f'Cross-Wavelet Transform of {var} between Surface and Bottom - {location}')
        plt.ylabel('Scale')
        plt.xlabel('Time')
        plt.savefig(f'wavelet_analysis_{location}/Cross_Wavelet_{var}.png')
        plt.close()
        
        # Wavelet coherence with river discharge and tidal height
        for factor in list(river_data.keys()) + ['Tide_Height']:
            coef_var, _ = pywt.cwt(data[f'{var}_0_10'].dropna(), scales, wavelet)
            coef_factor, _ = pywt.cwt(data[factor].dropna(), scales, wavelet)
            
            coherence = np.abs(coef_var * np.conj(coef_factor)) ** 2 / \
                        (np.abs(coef_var) ** 2 * np.abs(coef_factor) ** 2)
            
            plt.figure(figsize=(12, 8))
            plt.imshow(coherence, extent=[0, len(data), 1, 128], aspect='auto',
                       cmap='viridis', interpolation='nearest', vmin=0, vmax=1)
            plt.colorbar(label='Coherence')
            plt.title(f'Wavelet Coherence between {var} (Surface) and {factor} - {location}')
            plt.ylabel('Scale')
            plt.xlabel('Time')
            plt.savefig(f'wavelet_analysis_{location}/Coherence_{var}_{factor}.png')
            plt.close()
    
    # Global wavelet spectrum
    for var in variables:
        plt.figure(figsize=(10, 6))
        for depth in depths:
            coef, freqs = pywt.cwt(data[f'{var}_{depth}'].dropna(), scales, wavelet)
            power = (np.abs(coef) ** 2).mean(axis=1)
            plt.semilogx(1/freqs, power, label=f'{depth.replace("_", "-")}%')
        plt.title(f'Global Wavelet Spectrum of {var} - {location}')
        plt.xlabel('Period')
        plt.ylabel('Power')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'wavelet_analysis_{location}/Global_Wavelet_Spectrum_{var}.png')
        plt.close()
    
    # Scale-averaged wavelet power
    for var in variables:
        scale_power = np.abs(pywt.cwt(data[f'{var}_0_10'].dropna(), scales, wavelet)[0]) ** 2
        scale_avg_power = scale_power.mean(axis=0)
        
        plt.figure(figsize=(10, 6))
        plt.plot(data['Time'], scale_avg_power)
        plt.title(f'Scale-Averaged Wavelet Power of {var} (Surface) - {location}')
        plt.xlabel('Time')
        plt.ylabel('Power')
        plt.savefig(f'wavelet_analysis_{location}/Scale_Averaged_Power_{var}.png')
        plt.close()

# Load data
file_paths = [f'P{i}.txt' for i in range(1, n+1)]  
river_files = [f'R{i}.csv' for i in range(1, 6)]
tide_file = 'Tide.csv'

river_data = load_river_data(river_files)
tide_data = load_tide_data(tide_file)

# Run wavelet analysis 
for path in file_paths:
    data = load_data(path)
    location = data['Location'].iloc[0]
    wavelet_analysis(data, location, river_data, tide_data)

# Comparative analysis 
variables = ['SAL', 'TEMP', 'V_x', 'V_y']
wavelet = 'cmor1.5-1.0'
scales = np.arange(1, 129)

for var in variables:
    plt.figure(figsize=(12, 8))
    for path in file_paths:
        data = load_data(path)
        location = data['Location'].iloc[0]
        coef, freqs = pywt.cwt(data[f'{var}_0_10'].dropna(), scales, wavelet)
        power = (np.abs(coef) ** 2).mean(axis=1)
        plt.semilogx(1/freqs, power, label=location)
    plt.title(f'Global Wavelet Spectrum of Surface {var} across Locations')
    plt.xlabel('Period')
    plt.ylabel('Power')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'wavelet_analysis_comparison/Global_Wavelet_Spectrum_{var}_comparison.png')
    plt.close()
