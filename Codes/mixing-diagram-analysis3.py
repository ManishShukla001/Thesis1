import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.spatial import ConvexHull
from matplotlib.colors import LinearSegmentedColormap
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
def mixing_diagram_analysis(data_list, river_data, tide_data):
    # Ensure output directory exists
    os.makedirs('mixing_diagram_analysis', exist_ok=True)
    # Combine all location data
    combined_data = pd.concat(data_list, keys=range(len(data_list)))
    combined_data = combined_data.reset_index(level=0).rename(columns={'level_0': 'Location_Index'})
    # Merge with river and tide data
    for river, discharge in river_data.items():
        combined_data = combined_data.merge(discharge, left_on='Time', right_on='Date', how='left')
        combined_data = combined_data.rename(columns={'Discharge': f'{river}_Discharge'})
    combined_data = combined_data.merge(tide_data, on='Time', how='left')
    depths = [f'{i*10}_{(i+1)*10}' for i in range(10)]
    variables = ['TEMP', 'TURB', 'V_x', 'V_y']
    # Basic Mixing Diagrams
    for var in variables:
        plt.figure(figsize=(12, 8))
        for depth in depths:
            sal_values = combined_data[f'SAL_{depth}']
            var_values = combined_data[f'{var}_{depth}']
            plt.scatter(sal_values, var_values, alpha=0.5, label=f'{depth}')
        plt.title(f'Mixing Diagram: {var} vs Salinity')
        plt.xlabel('Salinity')
        plt.ylabel(var)
        plt.legend(title='Depth Range', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'mixing_diagram_analysis/mixing_diagram_{var}.png')
        plt.close()
    # Vertical Profile Mixing Diagrams
    for var in variables:
        plt.figure(figsize=(12, 8))
        for time in combined_data['Time'].unique()[::24]:  # Plot every 24th time step
            data_slice = combined_data[combined_data['Time'] == time]
            sal_profile = [data_slice[f'SAL_{depth}'].mean() for depth in depths]
            var_profile = [data_slice[f'{var}_{depth}'].mean() for depth in depths]
            plt.plot(sal_profile, var_profile, '-o', alpha=0.5, label=time.strftime('%Y-%m-%d'))
        plt.title(f'Vertical Profile Mixing Diagram: {var} vs Salinity')
        plt.xlabel('Salinity')
        plt.ylabel(var)
        plt.legend(title='Date', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'mixing_diagram_analysis/vertical_profile_mixing_{var}.png')
        plt.close()
    # Mixing Line Analysis
    for var in variables:
        for depth in depths:
            sal_values = combined_data[f'SAL_{depth}']
            var_values = combined_data[f'{var}_{depth}']
            # Remove NaN values
            mask = ~np.isnan(sal_values) & ~np.isnan(var_values)
            sal_values = sal_values[mask]
            var_values = var_values[mask]
            if len(sal_values) > 1:
                slope, intercept, r_value, p_value, std_err = linregress(sal_values, var_values)
                plt.figure(figsize=(10, 6))
                plt.scatter(sal_values, var_values, alpha=0.5)
                plt.plot(sal_values, slope * sal_values + intercept, 'r', label='Mixing line')
                plt.title(f'Mixing Line: {var} vs Salinity at {depth} depth')
                plt.xlabel('Salinity')
                plt.ylabel(var)
                plt.legend()
                plt.text(0.05, 0.95, f'R² = {r_value**2:.4f}\np-value = {p_value:.4f}', 
                         transform=plt.gca().transAxes, verticalalignment='top')
                plt.grid(True)
                plt.savefig(f'mixing_diagram_analysis/mixing_line_{var}_{depth}.png')
                plt.close()
                # Save mixing line parameters
                with open(f'mixing_diagram_analysis/mixing_line_params_{var}_{depth}.txt', 'w') as f:
                    f.write(f'Slope: {slope}\n')
                    f.write(f'Intercept: {intercept}\n')
                    f.write(f'R-squared: {r_value**2}\n')
                    f.write(f'p-value: {p_value}\n')
                    f.write(f'Standard Error: {std_err}\n')
    # Conservative Mixing Analysis
    for var in variables:
        plt.figure(figsize=(12, 8))
        for depth in depths:
            sal_values = combined_data[f'SAL_{depth}']
            var_values = combined_data[f'{var}_{depth}']
            # Remove NaN values
            mask = ~np.isnan(sal_values) & ~np.isnan(var_values)
            sal_values = sal_values[mask]
            var_values = var_values[mask]
            if len(sal_values) > 2:
                # Calculate convex hull
                points = np.column_stack((sal_values, var_values))
                hull = ConvexHull(points)
                # Plot data points and convex hull
                plt.scatter(sal_values, var_values, alpha=0.5, label=depth)
                for simplex in hull.simplices:
                    plt.plot(points[simplex, 0], points[simplex, 1], 'r-')
        plt.title(f'Conservative Mixing Analysis: {var} vs Salinity')
        plt.xlabel('Salinity')
        plt.ylabel(var)
        plt.legend(title='Depth Range', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'mixing_diagram_analysis/conservative_mixing_{var}.png')
        plt.close()
    # Temporal variations in mixing
    for var in variables:
        plt.figure(figsize=(12, 8))
        for depth in depths:
            sal_values = combined_data[f'SAL_{depth}']
            var_values = combined_data[f'{var}_{depth}']
            times = combined_data['Time']
            scatter = plt.scatter(sal_values, var_values, c=times, cmap='viridis', alpha=0.5)
        plt.colorbar(scatter, label='Time')
        plt.title(f'Temporal Variations in Mixing: {var} vs Salinity')
        plt.xlabel('Salinity')
        plt.ylabel(var)
        plt.grid(True)
        plt.savefig(f'mixing_diagram_analysis/temporal_mixing_{var}.png')
        plt.close()
    # Influence of river discharge and tidal height on mixing
    for var in variables:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        for depth in depths:
            sal_values = combined_data[f'SAL_{depth}']
            var_values = combined_data[f'{var}_{depth}']
            river_discharge = combined_data['R1_Discharge']  # Assuming R1 is the main river
            tidal_height = combined_data['Tide_Height']
            scatter1 = ax1.scatter(sal_values, var_values, c=river_discharge, cmap='viridis', alpha=0.5)
            scatter2 = ax2.scatter(sal_values, var_values, c=tidal_height, cmap='viridis', alpha=0.5)
        fig.colorbar(scatter1, ax=ax1, label='River Discharge')
        fig.colorbar(scatter2, ax=ax2, label='Tidal Height')
        ax1.set_title(f'Influence of River Discharge on Mixing: {var} vs Salinity')
        ax2.set_title(f'Influence of Tidal Height on Mixing: {var} vs Salinity')
        ax1.set_xlabel('Salinity')
        ax2.set_xlabel('Salinity')
        ax1.set_ylabel(var)
        ax2.set_ylabel(var)
        ax1.grid(True)
        ax2.grid(True)
        plt.tight_layout()
        plt.savefig(f'mixing_diagram_analysis/discharge_tide_influence_{var}.png')
        plt.close()
    # Vertical structure of mixing
    for var in variables:
        plt.figure(figsize=(12, 8))
        depth_values = [int(d.split('_')[0]) for d in depths]
        for time in combined_data['Time'].unique()[::24]:  # Plot every 24th time step
            data_slice = combined_data[combined_data['Time'] == time]
            sal_profile = [data_slice[f'SAL_{depth}'].mean() for depth in depths]
            var_profile = [data_slice[f'{var}_{depth}'].mean() for depth in depths]
            plt.plot(sal_profile, depth_values, '-o', alpha=0.5, label=time.strftime('%Y-%m-%d'))
        plt.title(f'Vertical Structure of Mixing: Salinity vs Depth for {var}')
        plt.xlabel('Salinity')
        plt.ylabel('Depth (%)')
        plt.legend(title='Date', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.gca().invert_yaxis()  # Invert y-axis to show depth increasing downwards
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'mixing_diagram_analysis/vertical_structure_{var}.png')
        plt.close()
    # Mixing intensity analysis
    for var in variables:
        plt.figure(figsize=(12, 8))
        depth_values = [int(d.split('_')[0]) for d in depths]
        
        mixing_intensity = []
        for depth in depths:
            sal_values = combined_data[f'SAL_{depth}']
            var_values = combined_data[f'{var}_{depth}']
            # Calculate mixing intensity as the standard deviation of the variable
            intensity = np.std(var_values)
            mixing_intensity.append(intensity)
        plt.plot(mixing_intensity, depth_values, '-o')
        plt.title(f'Mixing Intensity Profile for {var}')
        plt.xlabel(f'Mixing Intensity (std dev of {var})')
        plt.ylabel('Depth (%)')
        plt.gca().invert_yaxis()  # Invert y-axis to show depth increasing downwards
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'mixing_diagram_analysis/mixing_intensity_{var}.png')
        plt.close()
    # Stratification analysis
    plt.figure(figsize=(12, 8))
    times = combined_data['Time'].unique()
    stratification = []
    for time in times:
        data_slice = combined_data[combined_data['Time'] == time]
        surface_sal = data_slice['SAL_0_10'].mean()
        bottom_sal = data_slice['SAL_90_100'].mean()
        stratification.append(bottom_sal - surface_sal)
    plt.plot(times, stratification)
    plt.title('Temporal Variation of Stratification')
    plt.xlabel('Time')
    plt.ylabel('Stratification (Bottom - Surface Salinity)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('mixing_diagram_analysis/stratification_timeseries.png')
    plt.close()
# Load data
file_paths = [f'P{i}.txt' for i in range(1, n+1)]  # Assuming n locations
river_files = [f'R{i}.csv' for i in range(1, 6)]
tide_file = 'Tide.csv'
data_list = [load_data(path) for path in file_paths]
river_data = load_river_data(river_files)
tide_data = load_tide_data(tide_file)
# Run mixing diagram analysis
mixing_diagram_analysis(data_list, river_data, tide_data)