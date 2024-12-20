import pandas as pd
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def combine_met_files(input_dir: Path = Path('QAQC'), 
                     output_file: Path = Path('2024_CEQUAL_met_inputs.csv'),
                     start_julian_day: float = 37529.0):  # October 1, 2024
    """
    Combine meteorological data files into a single CE-QUAL formatted file.
    
    Parameters:
        input_dir (Path): Directory containing the input CSV files
        output_file (Path): Path for the output combined CSV file
        start_julian_day (float): Starting Julian day for the dataset
    """
    try:
        # Dictionary to store our dataframes with their purpose
        dfs = {}
        
        # Read all CSV files from the input directory
        for csv_file in input_dir.glob('*.csv'):
            df = pd.read_csv(csv_file)
            df['DATE TIME'] = pd.to_datetime(df['DATE TIME'])
            df = df.set_index('DATE TIME')
            dfs[csv_file.stem] = df
            logging.info(f"Loaded {csv_file.name}")

        # Create a complete datetime index
        min_date = min(df.index.min() for df in dfs.values())
        max_date = max(df.index.max() for df in dfs.values())
        complete_index = pd.date_range(start=min_date, end=max_date, freq='H')

        # Initialize the final dataframe with the complete datetime index
        final_df = pd.DataFrame(index=complete_index)

        # Calculate Julian days starting from the specified start day
        hours_elapsed = (final_df.index - final_df.index[0]).total_seconds() / 3600
        final_df['JDAY'] = start_julian_day + (hours_elapsed / 24)

        # Map the input files to CE-QUAL variables (no unit conversion)
        mappings = {
            'air_temp_avg_hourly_QAQC': ('Temp_C', 'TAIR'),
            'dewpoint_temp_hourly_QAQC': ('DP_temp_C', 'TDEW'),
            'wind_speed_hourly_QAQC': ('speed_ms', 'WIND'),
            'atmospheric_pressure_hourly_QAQC': ('Pressure_In', 'PHI'),
            'cloud_cover_hourly_QAQC': ('cloud_cover', 'CLOUD'),
            'solar_rad_avg_hourly_QAQC': ('Radiation_Wm2', 'SRO')
        }

        # Process each variable
        for file_base, (input_col, output_col) in mappings.items():
            if file_base in dfs:
                final_df[output_col] = dfs[file_base][input_col]
                logging.info(f"Processed {file_base} into {output_col}")

        # Ensure all required columns are present and in correct order
        required_columns = ['JDAY', 'TAIR', 'TDEW', 'WIND', 'PHI', 'CLOUD', 'SRO']
        for col in required_columns:
            if col not in final_df.columns:
                logging.warning(f"Missing required column: {col}")
                final_df[col] = 0  # or another appropriate default value

        # Reorder columns to match CE-QUAL format
        final_df = final_df[required_columns]

        # Save the combined file
        final_df.to_csv(output_file)
        logging.info(f"Successfully created combined file: {output_file}")
        
        return final_df

    except Exception as e:
        logging.error(f"Error combining meteorological files: {str(e)}")
        raise

if __name__ == "__main__":
    # Example usage with explicit Julian day start
    combine_met_files(
        input_dir=Path('QAQC'),
        output_file=Path('2024_CEQUAL_met_inputs.csv'),
        start_julian_day=37529.0  # October 1, 2024
    )
