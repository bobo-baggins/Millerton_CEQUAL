import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class MetDataProcessor:
    def __init__(self, start_date='2023-10-01 00:00:00', end_date='2023-12-31 23:00:00'):
        self.start_date = start_date
        self.end_date = end_date
        self.setup_plotting()
        
    def setup_plotting(self):
        """Initialize plotting parameters"""
        sns.set_style("whitegrid", {"axes.facecolor": "1",'axes.edgecolor': '0.6','grid.color': '0.6'})
        sns.set_context({'grid.linewidth':'1'})
        plt.rcParams['figure.figsize'] = (10, 5)
        plt.rcParams['font.size'] = 12
        plt.rcParams['lines.linewidth'] = 1.5
        plt.rcParams['lines.linestyle'] = '-'
        plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
        plt.rcParams['axes.titlesize'] = plt.rcParams['font.size']
        plt.rcParams['legend.fontsize'] = 0.9*plt.rcParams['font.size']
        plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
        plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']

    def process_cdec_data(self):
        """Process CDEC data first"""
        from CDEC_Formatter import CDECFormatter
        formatter = CDECFormatter()
        formatter.process_all()
        logging.info("Completed CDEC data processing")

    def process_basic_met_data(self):
        """Process and QA/QC basic meteorological data"""
        # Create output directory if it doesn't exist
        Path('QAQC').mkdir(exist_ok=True)
        
        # Process air temperature
        df_temp = pd.read_csv('hourly/air_temp_avg_hourly.csv', index_col=0, parse_dates=True)
        df_temp['Temp_C'] = (df_temp.Temp_F-32)*(5/9)
        df_temp.to_csv('QAQC/air_temp_avg_hourly_QAQC.csv')
        
        # Process other meteorological variables
        variables = {
            'atmospheric_pressure_hourly.csv': 'Pressure_In',
            'solar_rad_avg_hourly.csv': 'Radiation_Wm2',
            'wind_direction_hourly.csv': 'Dir_degrees',
            'wind_speed_hourly.csv': 'speed_mph',
            'relative_hum_hourly.csv': 'pcnt_hum'
        }
        
        for filename, column in variables.items():
            df = pd.read_csv(f'hourly/{filename}', index_col=0, parse_dates=True)
            df = df.interpolate()
            df.to_csv(f'QAQC/{filename.replace(".csv", "_QAQC.csv")}')
            
        logging.info("Completed basic meteorological data processing")

    def calculate_dewpoint(self):
        """Calculate dewpoint temperature"""
        df = pd.read_csv('QAQC/air_temp_avg_hourly_QAQC.csv', index_col=0, parse_dates=True)
        dfrh = pd.read_csv('QAQC/relative_hum_hourly_QAQC.csv', index_col=0, parse_dates=True)
        
        df['rel_hum'] = dfrh.pcnt_hum
        RH = df.rel_hum.values * 0.01
        T = df.Temp_C.values
        Td = (112 + 0.9*T)*(RH**(1/8)) - 112 + 0.1*T
        df['DP_temp_C'] = Td
        df.to_csv('QAQC/dewpoint_temp_hourly_QAQC.csv')
        logging.info("Completed dewpoint temperature calculation")

    def calculate_cloud_cover(self):
        """Calculate cloud cover"""
        df = pd.read_csv('QAQC/solar_rad_avg_hourly_QAQC.csv', index_col=0, parse_dates=True)
        
        # Calculate solar parameters
        df['day_of_year'] = df.index.dayofyear + df.index.hour/24
        df['hour_of_day'] = df.index.hour
        
        # Solar calculations
        num_days = df.day_of_year.values
        declination_angle = -23.44 * np.cos(np.radians((360/365)*(num_days - 1 + 10)))
        hour = df.hour_of_day.values
        hour_angle = 15*(hour-12)
        
        # Calculate solar elevation angle
        lat = 32.32  # latitude
        sin_alpha = (np.sin(np.radians(lat)) * np.sin(np.radians(declination_angle)) +
                    np.cos(np.radians(lat)) * np.cos(np.radians(declination_angle)) *
                    np.cos(np.radians(hour_angle)))
        solar_elevation_angle = np.degrees(np.arcsin(sin_alpha))
        
        # Calculate cloud cover
        theta_p = solar_elevation_angle
        theta = np.zeros(len(theta_p))
        for t in range(1, len(theta_p)):
            theta[t] = (theta_p[t-1] + theta_p[t])/2
            
        clear_sky_insolation = 990*np.sin(np.radians(theta))-30
        
        R = df.Radiation_Wm2.values
        R0 = clear_sky_insolation
        cloud_cover = ((1/0.65)*(1-(R/R0)))**(1/2)
        
        # Process cloud cover data
        cloud_cover_new = np.nan_to_num(cloud_cover, nan=0.0)
        cloud_cover_new = np.clip(cloud_cover_new, 0, 1)
        
        df['cloud_cover'] = cloud_cover_new
        df = df[(df.index.hour >= 13) & (df.index.hour <= 15)]
        
        # Reindex and interpolate
        ix = pd.date_range(self.start_date, self.end_date, freq='h')
        df = df.reindex(ix)
        df.loc[self.start_date, 'cloud_cover'] = 1
        df.loc[self.end_date, 'cloud_cover'] = 1
        df['cloud_cover'] = df.cloud_cover.interpolate()
        df['cloud_cover_10'] = df.cloud_cover * 10
        
        df.to_csv('QAQC/cloud_cover_hourly_QAQC.csv')
        logging.info("Completed cloud cover calculation")

    def create_final_output(self):
        """Create final CE-QUAL formatted output"""
        from Final_Formatter import combine_met_files
        combine_met_files()
        logging.info("Created final CE-QUAL formatted output")

    def process_all(self):
        """Run all processing steps in order"""
        self.process_cdec_data()
        self.process_basic_met_data()
        self.calculate_dewpoint()
        self.calculate_cloud_cover()
        self.create_final_output()
        logging.info("Completed all meteorological data processing")

if __name__ == "__main__":
    processor = MetDataProcessor()
    processor.process_all()