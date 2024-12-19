import pandas as pd
import os
from datetime import datetime
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="Workbook contains no default style, apply openpyxl's default")
# Define the directory where the .xlsx files are located
INPUT_DIR = Path('cdec')
OUTPUT_DIR = Path('hourly')
DATE_TIME_FORMAT = '%m/%d/%Y %H:%M'

COLUMN_MAPPING = {
    'atmospheric_pressure': ('Pressure_In', 'atmospheric_pressure_hourly.csv'),
    'wind_direction': ('Dir_degrees', 'wind_direction_hourly.csv'),
    'wind_speed': ('speed_mph', 'wind_speed_hourly.csv'),
    'relative_humidity': ('pcnt_hum', 'relative_hum_hourly.csv'),
    'air_temp': ('Temp_F', 'air_temp_avg_hourly.csv'),
}

# Get a list of all .xlsx files in that directory
xlsx_files = [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]
# WHAT ORDER TO FILES GO IN? THIS IS MISING THE DATA UP betweeen output sheets i.e. wind speed into temp....


# Loop over the list of .xlsx files
for i,xlsx_file in enumerate(xlsx_files):
    # Load the file into a DataFrame
    df = pd.read_excel(os.path.join(input_dir, xlsx_file))

    # Select only the 'DATE', 'TIME', and 'Value' columns
    df = df[['DATE TIME', 'VALUE']]
    
    # Rename the 'VALUE' column
    df = df.rename(columns={'VALUE': new_column_names[i]})

    # Change date format to 'month/day/year hour:min'
    df['DATE TIME'] = pd.to_datetime(df['DATE TIME'])

    df['DATE TIME'] = df['DATE TIME'].dt.strftime('%m/%d/%Y %H:%M')


    # Write the DataFram:w
    # e to a new CSV file in the 'hourly' directory
    output_file = os.path.join(output_dir, output_file_names[i])
   # output_file = os.path.join(output_dir, xlsx_file.replace('.xlsx', '_hourly.csv'))
    df.to_csv(output_file, index=False)


    # format Solar Rad csv

# Load the file into a DataFrame
Solar_df = pd.read_csv('cdec/CIMIS_Solar.csv', dtype={'Hour (PST)': str})

# Select only the 'DATE', 'TIME', and 'Value' columns
Solar_df = Solar_df[['Date', 'Hour (PST)', 'Sol Rad (W/sq.m)']]
Solar_df = Solar_df.rename(columns={'Sol Rad (W/sq.m)': 'Radiation_Wm2'})

Solar_df['Hour (PST)'] = Solar_df['Hour (PST)'].astype(str).str.split('.', expand=True)[0]

Solar_df['DATE TIME'] = Solar_df['Date'] + ' ' + Solar_df['Hour (PST)']

Solar_df=Solar_df[['DATE TIME', 'Radiation_Wm2']]

# Convert 'DateTime' to datetime format
Solar_df['DATE TIME'] = pd.to_datetime(Solar_df['DATE TIME'], format='%m/%d/%Y %H%M', errors='coerce')   


# Write the DataFrame to a new CSV file in the 'hourly' directory
output_file = os.path.join(output_dir, 'solar_rad_avg_hourly.csv')
Solar_df.to_csv(output_file, index=False)
