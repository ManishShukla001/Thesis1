import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from scipy.fft import fft, fftfreq
import utide


data = pd.read_csv('Tidalpropagation2.txt', delimiter='\t')


data['Time'] = pd.to_datetime(data['Time'], format='%d-%m-%Y %H:%M')


data.set_index('Time', inplace=True)


data = data.interpolate()

# Time series plots
plt.figure(figsize=(12, 6))
plt.plot(data.index, data['WS_Loc1'], label='Kochi')
plt.plot(data.index, data['WS_Loc2'], label='Thaneermukkom')
plt.xlabel('Time')
plt.ylabel('Water Surface Elevation')
plt.legend()
plt.title('a) Time Series of Water Surface Elevation')
plt.show()

# Time lag analysis using cross-correlation
corr = correlate(data['WS_Loc1'], data['WS_Loc2'], method='fft')
lag = np.argmax(corr) - (len(data['WS_Loc1']) - 1)
print(f"Time lag between Loc1 and Loc2: {lag} hours")

# Range and amplitude analysis
max_elev1, min_elev1 = data['WS_Loc1'].max(), data['WS_Loc1'].min()
max_elev2, min_elev2 = data['WS_Loc2'].max(), data['WS_Loc2'].min()
range_elev1, range_elev2 = max_elev1 - min_elev1, max_elev2 - min_elev2
print(f"Range of elevation at Loc1: {range_elev1}")
print(f"Range of elevation at Loc2: {range_elev2}")

# Harmonic analysis using FFT
fft_loc1 = fft(data['WS_Loc1'].values)
fft_loc2 = fft(data['WS_Loc2'].values)
freq_loc1 = fftfreq(len(data['WS_Loc1']), 1)  # Assuming data is hourly
freq_loc2 = fftfreq(len(data['WS_Loc2']), 1)  # Assuming data is hourly

plt.figure(figsize=(12, 6))
plt.plot(freq_loc1[1:], np.abs(fft_loc1[1:])*2/len(data['WS_Loc1']), label='Kochi')
plt.plot(freq_loc2[1:], np.abs(fft_loc2[1:])*2/len(data['WS_Loc2']), label='Thaneermukkom Barrage')
plt.xlabel('Frequency (1/hour)')
plt.ylabel('Amplitude')
plt.legend()
plt.title('b) Fourier Analysis of Water Surface Elevation')
plt.show()

# Phase difference
phase_diff = np.unwrap(np.angle(fft_loc1) - np.angle(fft_loc2))
plt.figure(figsize=(12, 6))
plt.plot(freq_loc1[1:], phase_diff[1:])
plt.xlabel('Frequency (1/hour)')
plt.ylabel('Phase Difference (radians)')
plt.title('c) Phase Difference Analysis')
plt.show()

# Tidal propagation speed
distance_km = 43872 / 1000  # Convert meters to kilometers
# Use the absolute value of the time lag
propagation_speed = distance_km / abs(lag)  # km/hour
print(f"Tidal Propagation Speed: {propagation_speed} km/hour")

# Energy dissipation
# Simplified approach: assume energy is proportional to amplitude squared
energy_dissipation = (range_elev1**2 - range_elev2**2) / range_elev1**2
print(f"Energy Dissipation: {energy_dissipation*100}%")

# Tidal constituent analysis using utide
def analyze_tidal_constituents(data, location, latitude):
    coef = utide.solve(data.index.to_pydatetime(), data[location].values, lat=latitude)
    return coef

latitude = 0.1  

constituents_loc1 = analyze_tidal_constituents(data, 'WS_Loc1', latitude)
constituents_loc2 = analyze_tidal_constituents(data, 'WS_Loc2', latitude)

print(f"Tidal Constituents for Loc1: {constituents_loc1}")
print(f"Tidal Constituents for Loc2: {constituents_loc2}")

# Compare main tidal constituents between locations
main_constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1']
def compare_constituents(constituents1, constituents2, main_constituents):
    comparison = []
    for constituent in main_constituents:
        indices1 = np.where(constituents1['name'] == constituent)[0]
        indices2 = np.where(constituents2['name'] == constituent)[0]
        if indices1.size > 0 and indices2.size > 0:
            index1 = indices1[0]
            index2 = indices2[0]
            amp1 = constituents1['A'][index1]
            amp2 = constituents2['A'][index2]
            phase1 = constituents1['g'][index1]
            phase2 = constituents2['g'][index2]
            comparison.append({
                'Constituent': constituent,
                'Loc1 Amplitude': amp1,
                'Loc2 Amplitude': amp2,
                'Loc1 Phase': phase1,
                'Loc2 Phase': phase2
            })
    return pd.DataFrame(comparison)

comparison_df = compare_constituents(constituents_loc1, constituents_loc2, main_constituents)
print(comparison_df)

# Residual analysis
def calculate_residuals(data, location, coef):
    tide = utide.reconstruct(data.index.to_pydatetime(), coef)
    residuals = data[location] - tide.h
    return residuals

residuals_loc1 = calculate_residuals(data, 'WS_Loc1', constituents_loc1)
residuals_loc2 = calculate_residuals(data, 'WS_Loc2', constituents_loc2)

plt.figure(figsize=(12, 6))
plt.plot(data.index, residuals_loc1, label='Residuals at Kochi')
plt.plot(data.index, residuals_loc2, label='Residuals at Thaneermukkom')
plt.xlabel('Time')
plt.ylabel('Residual Water Surface Elevation')
plt.legend()
plt.title('d) Residual Analysis of Water Surface Elevation')
plt.show()

# Export results to Excel
with pd.ExcelWriter('tidal_analysis_results.xlsx') as writer:
    data.to_excel(writer, sheet_name='Water_Surface_Elevation')
    pd.DataFrame({'Frequency': freq_loc1[1:], 'Loc1 Amplitude': np.abs(fft_loc1[1:])*2/len(data['WS_Loc1'])}).to_excel(writer, sheet_name='FFT_Loc1')
    pd.DataFrame({'Frequency': freq_loc2[1:], 'Loc2 Amplitude': np.abs(fft_loc2[1:])*2/len(data['WS_Loc2'])}).to_excel(writer, sheet_name='FFT_Loc2')
    pd.DataFrame({'Frequency': freq_loc1[1:], 'Phase Difference': phase_diff[1:]}).to_excel(writer, sheet_name='Phase_Difference')
    pd.DataFrame({'Time': data.index, 'Loc1 Residuals': residuals_loc1, 'Loc2 Residuals': residuals_loc2}).to_excel(writer, sheet_name='Residuals')
    comparison_df.to_excel(writer, sheet_name='Tidal_Constituent_Comparison')
