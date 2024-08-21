import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.interpolate import griddata
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

def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def salinity_gradient_analysis(data_list, river_data, tide_data):
    # Ensure output directory exists
    os.makedirs('salinity_gradient_analysis', exist_ok=True)
    
    # Combine all location data
    combined_data = pd.concat(data_list, keys=range(len(data_list)))
    combined_data = combined_data.reset_index(level=0).rename(columns={'level_0': 'Location_Index'})
    
    # Merge with river and tide data
    for river, discharge in river_data.items():
        combined_data = combined_data.merge(discharge, left_on='Time', right_on='Date', how='left')
        combined_data = combined_data.rename(columns={'Discharge': f'{river}_Discharge'})
    combined_data = combined_data.merge(tide_data, on='Time', how='left')
    
    # Calculate distances between locations
    locations = combined_data.groupby('Location_Index').first()[['x', 'y']]
    distances = [0]
    for i in range(1, len(locations)):
        dist = calculate_distance(locations.iloc[i-1]['x'], locations.iloc[i-1]['y'],
                                  locations.iloc[i]['x'], locations.iloc[i]['y'])
        distances.append(distances[-1] + dist)
    
    depths = [int(d.split('_')[1]) for d in combined_data.columns if d.startswith('SAL_')]
    
    # Longitudinal salinity gradient
    for time in combined_data['Time'].unique():
        data_slice = combined_data[combined_data['Time'] == time]
        
        plt.figure(figsize=(12, 8))
        for depth in depths:
            sal_values = data_slice[f'SAL_0_{depth}']
            plt.plot(distances, sal_values, 'o-', label=f'{depth}% depth')
        
        plt.title(f'Longitudinal Salinity Profile at {time}')
        plt.xlabel('Distance along estuary (km)')
        plt.ylabel('Salinity')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'salinity_gradient_analysis/longitudinal_profile_{time.strftime("%Y%m%d%H%M")}.png')
        plt.close()
    
    # Vertical salinity gradient
    for loc_index in combined_data['Location_Index'].unique():
        data_slice = combined_data[combined_data['Location_Index'] == loc_index]
        
        plt.figure(figsize=(10, 12))
        for time in data_slice['Time'].iloc[::24]:  # Plot every 24th time step (assuming hourly data)
            sal_values = [data_slice.loc[data_slice['Time'] == time, f'SAL_0_{d}'].values[0] for d in depths]
            plt.plot(sal_values, depths, 'o-', label=time.strftime('%Y-%m-%d'))
        
        plt.title(f'Vertical Salinity Profile at Location {loc_index}')
        plt.xlabel('Salinity')
        plt.ylabel('Depth (%)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.gca().invert_yaxis()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'salinity_gradient_analysis/vertical_profile_location_{loc_index}.png')
        plt.close()
    
    # Salinity contour plot (distance vs depth)
    for time in combined_data['Time'].iloc[::24]:  # Plot every 24th time step
        data_slice = combined_data[combined_data['Time'] == time]
        
        X, Y = np.meshgrid(distances, depths)
        Z = np.array([[data_slice.loc[data_slice['Location_Index'] == i, f'SAL_0_{d}'].values[0] 
                       for i in range(len(locations))] for d in depths])
        
        plt.figure(figsize=(12, 8))
        contour = plt.contourf(X, Y, Z, levels=20, cmap='viridis')
        plt.colorbar(contour, label='Salinity')
        plt.title(f'Salinity Contour at {time}')
        plt.xlabel('Distance along estuary (km)')
        plt.ylabel('Depth (%)')
        plt.gca().invert_yaxis()
        plt.savefig(f'salinity_gradient_analysis/salinity_contour_{time.strftime("%Y%m%d%H%M")}.png')
        plt.close()
    
    # Correlation analysis
    correlations = {}
    for loc_index in combined_data['Location_Index'].unique():
        data_slice = combined_data[combined_data['Location_Index'] == loc_index]
        corr_data = {}
        
        for depth in depths:
            corr_rivers = {river: pearsonr(data_slice[f'SAL_0_{depth}'], data_slice[f'{river}_Discharge'])[0]
                           for river in river_data.keys()}
            corr_tide = pearsonr(data_slice[f'SAL_0_{depth}'], data_slice['Tide_Height'])[0]
            corr_data[depth] = {**corr_rivers, 'Tide_Height': corr_tide}
        
        correlations[loc_index] = pd.DataFrame(corr_data).T
    
    # Plot correlation heatmaps
    for loc_index, corr_df in correlations.items():
        plt.figure(figsize=(10, 8))
        plt.imshow(corr_df.values, aspect='auto', cmap='coolwarm', vmin=-1, vmax=1)
        plt.colorbar(label='Correlation Coefficient')
        plt.title(f'Correlation of Salinity with Rivers and Tide at Location {loc_index}')
        plt.ylabel('Depth (%)')
        plt.xlabel('Factor')
        plt.yticks(range(len(depths)), depths)
        plt.xticks(range(len(corr_df.columns)), corr_df.columns, rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f'salinity_gradient_analysis/correlation_heatmap_location_{loc_index}.png')
        plt.close()
        
        # Save correlation data
        corr_df.to_csv(f'salinity_gradient_analysis/correlation_data_location_{loc_index}.csv')
    
    # Time series of salinity gradient
    surface_gradient = []
    bottom_gradient = []
    times = []
    
    for time in combined_data['Time'].unique():
        data_slice = combined_data[combined_data['Time'] == time]
        surface_sal = data_slice['SAL_0_10']
        bottom_sal = data_slice['SAL_0_100']
        
        surface_grad = np.gradient(surface_sal, distances)
        bottom_grad = np.gradient(bottom_sal, distances)
        
        surface_gradient.append(np.mean(np.abs(surface_grad)))
        bottom_gradient.append(np.mean(np.abs(bottom_grad)))
        times.append(time)
    
    plt.figure(figsize=(12, 6))
    plt.plot(times, surface_gradient, label='Surface')
    plt.plot(times, bottom_gradient, label='Bottom')
    plt.title('Time Series of Mean Absolute Salinity Gradient')
    plt.xlabel('Time')
    plt.ylabel('Mean Absolute Salinity Gradient (per km)')
    plt.legend()
    plt.grid(True)
    plt.savefig('salinity_gradient_analysis/salinity_gradient_time_series.png')
    plt.close()

# Load data
file_paths = [f'P{i}.txt' for i in range(1, n+1)]  # Assuming n locations
river_files = [f'R{i}.csv' for i in range(1, 6)]
tide_file = 'Tide.csv'

data_list = [load_data(path) for path in file_paths]
river_data = load_river_data(river_files)
tide_data = load_tide_data(tide_file)

# Run salinity gradient analysis
salinity_gradient_analysis(data_list, river_data, tide_data)