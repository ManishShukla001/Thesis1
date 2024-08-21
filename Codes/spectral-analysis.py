import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import welch, coherence
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

def spectral_analysis(data, location, river_data, tide_data):
    os.makedirs(f'spectral_analysis_{location}', exist_ok=True)
    
    # Merge with river and tide data
    for river, discharge in river_data.items():
        data = data.merge(discharge, left_on='Time', right_on='Date', how='left')
        data = data.rename(columns={'Discharge': f'{river}_Discharge'})
    data = data.merge(tide_data, on='Time', how='left')
    
    variables = ['SAL', 'TEMP', 'V_x', 'V_y']
    depths = [f'{i*10}_{(i+1)*10}' for i in range(10)]
    
    # Sampling frequency (assuming hourly data)
    fs = 1 / 3600

    for var in variables:
        # Power Spectral Density for each depth
        plt.figure(figsize=(12, 8))
        for depth in depths:
            f, Pxx = welch(data[f'{var}_{depth}'].dropna(), fs=fs, nperseg=8192)
            plt.semilogx(f, Pxx, label=f'{int(depth.split("_")[0])}%-{int(depth.split("_")[1])}%')
        plt.title(f'Power Spectral Density of {var} at different depths - {location}')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Power Spectral Density')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'spectral_analysis_{location}/PSD_{var}.png')
        plt.close()
        
        # Coherence between surface and bottom layers
        f, Cxy = coherence(data[f'{var}_0_10'], data[f'{var}_90_100'], fs=fs, nperseg=8192)
        plt.figure(figsize=(10, 6))
        plt.semilogx(f, Cxy)
        plt.title(f'Coherence of {var} between surface and bottom layers - {location}')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Coherence')
        plt.grid(True)
        plt.savefig(f'spectral_analysis_{location}/Coherence_{var}.png')
        plt.close()
        
        # Correlation with river discharge and tidal height
        correlations = []
        for depth in depths:
            corr_river = {river: pearsonr(data[f'{var}_{depth}'], data[f'{river}_Discharge'])[0] 
                          for river in river_data.keys()}
            corr_tide = pearsonr(data[f'{var}_{depth}'], data['Tide_Height'])[0]
            correlations.append({**corr_river, 'Tide': corr_tide})
        
        corr_df = pd.DataFrame(correlations, index=depths)
        corr_df.to_csv(f'spectral_analysis_{location}/Correlations_{var}.csv')
        
        # Heatmap of correlations
        plt.figure(figsize=(10, 8))
        plt.imshow(corr_df.values.T, aspect='auto', cmap='coolwarm', vmin=-1, vmax=1)
        plt.colorbar(label='Correlation Coefficient')
        plt.title(f'Correlation of {var} with River Discharges and Tidal Height - {location}')
        plt.ylabel('Factor')
        plt.xlabel('Depth')
        plt.yticks(range(len(corr_df.columns)), corr_df.columns)
        plt.xticks(range(len(depths)), [f'{d.split("_")[0]}-{d.split("_")[1]}%' for d in depths], rotation=45)
        plt.savefig(f'spectral_analysis_{location}/Correlation_Heatmap_{var}.png')
        plt.close()
    
    # Cross-spectral analysis between variables
    for var1, var2 in [('SAL', 'TEMP'), ('V_x', 'V_y')]:
        f, Cxy = coherence(data[f'{var1}_0_10'], data[f'{var2}_0_10'], fs=fs, nperseg=8192)
        plt.figure(figsize=(10, 6))
        plt.semilogx(f, Cxy)
        plt.title(f'Coherence between {var1} and {var2} at surface - {location}')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Coherence')
        plt.grid(True)
        plt.savefig(f'spectral_analysis_{location}/Coherence_{var1}_{var2}.png')
        plt.close()

# Load data
file_paths = [f'P{i}.txt' for i in range(1, n+1)]  # Assuming n locations
river_files = [f'R{i}.csv' for i in range(1, 6)]
tide_file = 'Tide.csv'

river_data = load_river_data(river_files)
tide_data = load_tide_data(tide_file)

# Run spectral analysis for each location
for path in file_paths:
    data = load_data(path)
    location = data['Location'].iloc[0]
    spectral_analysis(data, location, river_data, tide_data)

# Comparative analysis across locations
variables = ['SAL', 'TEMP', 'V_x', 'V_y']
depths = [f'{i*10}_{(i+1)*10}' for i in range(10)]

for var in variables:
    plt.figure(figsize=(12, 8))
    for path in file_paths:
        data = load_data(path)
        location = data['Location'].iloc[0]
        f, Pxx = welch(data[f'{var}_0_10'].dropna(), fs=1/3600, nperseg=8192)
        plt.semilogx(f, Pxx, label=location)
    plt.title(f'Power Spectral Density of Surface {var} across Locations')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'spectral_analysis_comparison/PSD_{var}_comparison.png')
    plt.close()
