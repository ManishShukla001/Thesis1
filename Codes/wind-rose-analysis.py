import pandas as pd
import matplotlib.pyplot as plt
from windrose import WindroseAxes
import numpy as np


df = pd.read_excel('P1_data.xlsx')


df['Date'] = pd.to_datetime(df['Date'])


df['Year'] = df['Date'].dt.year

# wind rose diagram
def create_windrose(data, year):
    wind_speed = data['v10']  # Assuming 'v10' is wind speed
    wind_dir = data['u10']    # Assuming 'u10' is wind direction

    # Create wind rose
    ax = WindroseAxes.from_ax()
    ax.bar(wind_dir, wind_speed, normed=True, opening=0.8, edgecolor='white')
    ax.set_legend()
    ax.set_title(f'Wind Rose Diagram - {year}')
    plt.savefig(f'wind_rose_{year}.png')
    plt.close()

# Filter data for 2009 and 2010
data_2009 = df[df['Year'] == 2009]
data_2010 = df[df['Year'] == 2010]

# Create wind rose diagrams
create_windrose(data_2009, 2009)
create_windrose(data_2010, 2010)

# Basic statistical analysis
def analyze_wind_data(data, year):
    wind_speed = data['v10']
    wind_dir = data['u10']
    
    print(f"\nWind Analysis for {year}:")
    print(f"Average Wind Speed: {wind_speed.mean():.2f}")
    print(f"Max Wind Speed: {wind_speed.max():.2f}")
    print(f"Min Wind Speed: {wind_speed.min():.2f}")
    print(f"Most Common Wind Direction: {wind_dir.mode().values[0]:.2f}")

analyze_wind_data(data_2009, 2009)
analyze_wind_data(data_2010, 2010)
