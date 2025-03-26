# This scripts post-processes model results and raw water temperature data from CDEC so that model performance can be evaluated
#FWQ is a sensor immeadately downstream of river release outlet

#Importing libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

#Reading in data
SJA_Results = pd.read_csv('../CEQUAL_model/EC_two_str4_seg31.csv', skiprows=2)
Weighted_Results = pd.read_csv('../CEQUAL_model/k2pW_two_str4_seg31.csv', skiprows=2)
K2P_Results = pd.read_csv('../CEQUAL_model/two_str4_seg31.csv', skiprows=2)
FWQ = pd.read_csv('FWQ_2024.csv')

#Converting Julian days to dates for Model Results
# Assuming the Julian days are for 2024 - adjust the year as needed
base_date = pd.to_datetime('2024-05-28 01:00')
SJA_Results['DATE TIME'] = pd.date_range(start=base_date, periods=len(SJA_Results), freq='1H')
Weighted_Results['DATE TIME'] = pd.date_range(start=base_date, periods=len(Weighted_Results), freq='1H')
K2P_Results['DATE TIME'] = pd.date_range(start=base_date, periods=len(K2P_Results), freq='1H')

#Converting date columns to datetime for sensor data
FWQ['DATE TIME'] = pd.to_datetime(FWQ['DATE TIME'])

#Remove -99 values before averaging
SJA_Results = SJA_Results[SJA_Results['T(C)'] != -99]
Weighted_Results = Weighted_Results[Weighted_Results['T(C)'] != -99]
K2P_Results = K2P_Results[K2P_Results['T(C)'] != -99]
FWQ = FWQ[FWQ['T(C)'] != -99]

#Averaging data to daily values
SJA_Results_Daily = SJA_Results.groupby(SJA_Results['DATE TIME'].dt.date).mean().round(2)
Weighted_Results_Daily = Weighted_Results.groupby(Weighted_Results['DATE TIME'].dt.date).mean().round(2)
K2P_Results_Daily = K2P_Results.groupby(K2P_Results['DATE TIME'].dt.date).mean().round(2)
FWQ_Daily = FWQ.groupby(FWQ['DATE TIME'].dt.date).mean().round(2)

#Saving data to csv, dropping the JDAY and second DATE TIME columns
SJA_Results_Daily = SJA_Results_Daily.drop(['JDAY', 'DATE TIME'], axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
Weighted_Results_Daily = Weighted_Results_Daily.drop(['JDAY', 'DATE TIME'], axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
K2P_Results_Daily = K2P_Results_Daily.drop(['JDAY', 'DATE TIME'], axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
FWQ_Daily = FWQ_Daily.drop('DATE TIME', axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})

#Save individual CSVs after renaming
SJA_Results_Daily.to_csv('SJA_Results_Daily.csv')
Weighted_Results_Daily.to_csv('Weighted_Results_Daily.csv')
K2P_Results_Daily.to_csv('K2P_Results_Daily.csv')
FWQ_Daily.to_csv('FWQ_Daily.csv')

#Create combined dataframe and save to CSV
combined_df = pd.DataFrame({
    'SJA_DAWT(C)': SJA_Results_Daily['DAWT(C)'],
    'Weighted_DAWT(C)': Weighted_Results_Daily['DAWT(C)'],
    'K2P_DAWT(C)': K2P_Results_Daily['DAWT(C)'],
    'FWQ_DAWT(C)': FWQ_Daily['DAWT(C)']
})
combined_df.to_csv('Combined_Daily_Temps.csv')

#Plotting data

# Set the figure size and DPI before plotting
plt.figure(figsize=(12, 8), dpi=100)  # Width: 12 inches, Height: 8 inches


plt.plot(SJA_Results_Daily['DAWT(C)'], label='2024 Actual Flows', color='red')
plt.plot(Weighted_Results_Daily['DAWT(C)'], label='K2P Weighted', color='blue')
plt.plot(K2P_Results_Daily['DAWT(C)'], label='K2P', color='purple')
plt.plot(FWQ_Daily['DAWT(C)'], label='FWQ Sensor Daily Avg.', linestyle=':', color='green')

# Calculate and plot temperature difference
#temp_diff = EC_Results_Daily['DAWT(C)'] - TCD_Results_Daily['DAWT(C)']
#plt.plot(temp_diff, label='EC-TCD Difference', color='purple', linestyle='--')

# Calculate max and average difference
#max_diff = temp_diff.max().round(2)
#avg_diff = temp_diff.mean().round(2)

# Calculate RMSE between EC and FWQ
# Align the indices first
common_indices = SJA_Results_Daily.index.intersection(FWQ_Daily.index)
sja_aligned = SJA_Results_Daily.loc[common_indices, 'DAWT(C)']
weighted_aligned = Weighted_Results_Daily.loc[common_indices, 'DAWT(C)']
k2p_aligned = K2P_Results_Daily.loc[common_indices, 'DAWT(C)']
fwq_aligned = FWQ_Daily.loc[common_indices, 'DAWT(C)']

# Calculate statistics for SJA
sja_rmse = np.sqrt(((sja_aligned - fwq_aligned) ** 2).mean()).round(2)
sja_r2 = round(r2_score(fwq_aligned, sja_aligned), 2)
mean_observed = fwq_aligned.mean()
sja_nse = (1 - np.sum((fwq_aligned - sja_aligned) ** 2) / 
          np.sum((fwq_aligned - mean_observed) ** 2)).round(2)

# Calculate statistics for Weighted K2P
weighted_rmse = np.sqrt(((weighted_aligned - fwq_aligned) ** 2).mean()).round(2)
weighted_r2 = round(r2_score(fwq_aligned, weighted_aligned), 2)
weighted_nse = (1 - np.sum((fwq_aligned - weighted_aligned) ** 2) / 
                np.sum((fwq_aligned - mean_observed) ** 2)).round(2)

# Calculate statistics for K2P
k2p_rmse = np.sqrt(((k2p_aligned - fwq_aligned) ** 2).mean()).round(2)
k2p_r2 = round(r2_score(fwq_aligned, k2p_aligned), 2)
k2p_nse = (1 - np.sum((fwq_aligned - k2p_aligned) ** 2) / 
           np.sum((fwq_aligned - mean_observed) ** 2)).round(2)

# Add text box with statistics for all scenarios
stats_text = (f'SJA Statistics:\n'
              f'R² : {sja_r2}\n'
              f'NSE: {sja_nse}\n'
              f'RMSE: {sja_rmse}°C\n\n'
              f'Weighted Statistics:\n'
              f'R² : {weighted_r2}\n'
              f'NSE: {weighted_nse}\n'
              f'RMSE: {weighted_rmse}°C\n\n'
              f'K2P Statistics:\n'
              f'R² : {k2p_r2}\n'
              f'NSE: {k2p_nse}\n'
              f'RMSE: {k2p_rmse}°C')

plt.text(0.02, 0.98, stats_text, 
         transform=plt.gca().transAxes, 
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'),
         verticalalignment='top')

# Add vertical dashed line at July 30th with label
july_30 = pd.to_datetime('2024-07-30').date()
plt.axvline(x=july_30, color='black', linestyle='--', alpha=0.7)

# Add boxed label with arrow
plt.annotate('K2P Data Available\nbeginning July 30th', 
            xy=(july_30, 19), # Point on the line
            xytext=(july_30 + pd.Timedelta(days=10), 19), # Text position (shifted right by 10 days)
            bbox=dict(facecolor='white', edgecolor='black', alpha=0.8),
            horizontalalignment='left',
            verticalalignment='top',
            arrowprops=dict(arrowstyle='-', color='black', linewidth=1))

# Set y-axis limits and ticks with 1.0 increments
plt.ylim(5, 20)
plt.yticks(np.arange(5, 20.5, 1))

# Add grid
plt.grid(True, linestyle='--', alpha=0.7, color='lightgray')

plt.xlabel('Date')
plt.ylabel('Avg. Water Temp. (°C)')
plt.legend()

# Save with high DPI and tight layout
plt.savefig('temperature_validation.png', 
            dpi=300,              # High DPI for clarity
            bbox_inches='tight',  # Removes extra whitespace
            format='png',         # Explicitly set format
            transparent=False,     # White background
            pad_inches=0.1)       # Small padding around the plot

plt.show()