# This scripts post-processes model results and raw water temperature data from CDEC so that model performance can be evaluated
#FWQ is a sensor immeadately downstream of river release outlet
#SJF is a sensor 1 mile downstream of river release outlet

#Importing libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Reading in data
EC_Results = pd.read_csv('../CEQUAL_model/two_str4_seg31.csv', skiprows=2)
#TCD_Results = pd.read_csv('../CEQUAL_model/two_str4_seg31.csv', skiprows=2)
FWQ = pd.read_csv('FWQ_2024.csv')
SJF = pd.read_csv('SJF_2024.csv')

#Converting Julian days to dates for Model Results
# Assuming the Julian days are for 2024 - adjust the year as needed
base_date = pd.to_datetime('2024-05-28 01:00')
EC_Results['DATE TIME'] = pd.date_range(start=base_date, periods=len(EC_Results), freq='1H')
#TCD_Results['DATE TIME'] = pd.date_range(start=base_date, periods=len(TCD_Results), freq='1H')

#Converting date columns to datetime for sensor data
FWQ['DATE TIME'] = pd.to_datetime(FWQ['DATE TIME'])
SJF['DATE TIME'] = pd.to_datetime(SJF['DATE TIME'])

#Remove -99 values before averaging
EC_Results = EC_Results[EC_Results['T(C)'] != -99]
#TCD_Results = TCD_Results[TCD_Results['T(C)'] != -99]
FWQ = FWQ[FWQ['T(C)'] != -99]
SJF = SJF[SJF['T(C)'] != -99]

#Averaging data to daily values
EC_Results_Daily = EC_Results.groupby(EC_Results['DATE TIME'].dt.date).mean().round(2)
#TCD_Results_Daily = TCD_Results.groupby(TCD_Results['DATE TIME'].dt.date).mean().round(2)
FWQ_Daily = FWQ.groupby(FWQ['DATE TIME'].dt.date).mean().round(2)
SJF_Daily = SJF.groupby(SJF['DATE TIME'].dt.date).mean().round(2)

#Saving data to csv, dropping the JDAY and second DATE TIME columns
EC_Results_Daily = EC_Results_Daily.drop(['JDAY', 'DATE TIME'], axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
#TCD_Results_Daily = TCD_Results_Daily.drop(['JDAY', 'DATE TIME'], axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
FWQ_Daily = FWQ_Daily.drop('DATE TIME', axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
SJF_Daily = SJF_Daily.drop('DATE TIME', axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})

#Save individual CSVs after renaming
EC_Results_Daily.to_csv('EC_Results_Daily.csv')
#TCD_Results_Daily.to_csv('TCD_Results_Daily.csv')
FWQ_Daily.to_csv('FWQ_Daily.csv')
#SJF_Daily.to_csv('SJF_Daily.csv')

#Create combined dataframe and save to CSV
combined_df = pd.DataFrame({
    'EC_DAWT(C)': EC_Results_Daily['DAWT(C)'],
    #'TCD_DAWT(C)': TCD_Results_Daily['DAWT(C)'],
    'FWQ_DAWT(C)': FWQ_Daily['DAWT(C)'],
    #'SJF_DAWT(C)': SJF_Daily['DAWT(C)']
})
combined_df.to_csv('Combined_Daily_Temps.csv')

#Plotting data

# Set the figure size and DPI before plotting
plt.figure(figsize=(12, 8), dpi=100)  # Width: 12 inches, Height: 8 inches


plt.plot(EC_Results_Daily['DAWT(C)'], label='2024 Actual Flows', color='red')
#plt.plot(TCD_Results_Daily['DAWT(C)'], label='Madera TCD Scenario', color='blue')
plt.plot(FWQ_Daily['DAWT(C)'], label='FWQ Sensor Daily Avg.', linestyle=':', color='green')
#plt.plot(SJF_Daily['DAWT(C)'], label='SJF')

# Calculate and plot temperature difference
#temp_diff = EC_Results_Daily['DAWT(C)'] - TCD_Results_Daily['DAWT(C)']
#plt.plot(temp_diff, label='EC-TCD Difference', color='purple', linestyle='--')

# Calculate max and average difference
#max_diff = temp_diff.max().round(2)
#avg_diff = temp_diff.mean().round(2)

# Calculate RMSE between EC and FWQ
# Align the indices first
common_indices = EC_Results_Daily.index.intersection(FWQ_Daily.index)
ec_aligned = EC_Results_Daily.loc[common_indices, 'DAWT(C)']
fwq_aligned = FWQ_Daily.loc[common_indices, 'DAWT(C)']
rmse = np.sqrt(((ec_aligned - fwq_aligned) ** 2).mean()).round(2)

# Calculate Nash-Sutcliffe Efficiency
mean_observed = fwq_aligned.mean()
nse = (1 - np.sum((fwq_aligned - ec_aligned) ** 2) / 
       np.sum((fwq_aligned - mean_observed) ** 2)).round(2)

# Add text box with statistics
plt.text(0.02, 0.98, f'RMSE : {rmse}°C\nNSE : {nse}', 
         transform=plt.gca().transAxes, 
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'),
         verticalalignment='top')

# Set y-axis limits and ticks with 0.5 increments
plt.ylim(0, 20)
plt.yticks(np.arange(0, 20.5, 1))

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