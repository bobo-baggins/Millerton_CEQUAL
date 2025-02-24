# This scripts post-processes model results and raw water temperature data from CDEC so that model performance can be evaluated
#FWQ is a sensor immeadately downstream of river release outlet
#SJF is a sensor 1 mile downstream of river release outlet

#Importing libraries
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg' if you prefer Qt
import matplotlib.pyplot as plt
import numpy as np

#Reading in data
Model_Results = pd.read_csv('../CEQUAL_outputs/2024/two_str4_seg31.csv', skiprows=2)
Old_Model = pd.read_fwf('../CEQUAL_outputs/2024/str_br1_old.opt', 
                        skiprows=2,
                        colspecs=[(1, 11),  # JDAY column
                                (45, 51)],  # 4th T(C) column
                        names=['JDAY', 'T(C)'])

Madera_TCD = pd.read_fwf('../CEQUAL_model/str_br1.opt', 
                        skiprows=2,
                        colspecs=[(1, 11),  # JDAY column
                                (45, 51)],  # 4th T(C) column
                        names=['JDAY', 'T(C)'])

FWQ = pd.read_csv('FWQ_2024.csv')
SJF = pd.read_csv('SJF_2024.csv')

print(Old_Model.head())

#Converting Julian days to dates for Model Results
# Assuming the Julian days are for 2024 - adjust the year as needed
base_date = pd.to_datetime('2024-05-28 01:00')
Model_Results['DATE TIME'] = pd.date_range(start=base_date, periods=len(Model_Results), freq='1h')
Madera_TCD['DATE TIME'] = pd.date_range(start=base_date, periods=len(Madera_TCD), freq='1h')
Old_Model['DATE TIME'] = pd.date_range(start=base_date, periods=len(Old_Model), freq='1h')
#Old_Model.to_csv('Old_Model_with_dates.csv')

#Converting date columns to datetime for sensor data
FWQ['DATE TIME'] = pd.to_datetime(FWQ['DATE TIME'])
SJF['DATE TIME'] = pd.to_datetime(SJF['DATE TIME'])

#Remove -99 values before averaging
Model_Results = Model_Results[Model_Results['T(C)'] != -99]
Madera_TCD = Madera_TCD[Madera_TCD['T(C)'] != -99]
Old_Model = Old_Model[Old_Model['T(C)'] != -99]
FWQ = FWQ[FWQ['T(C)'] != -99]
SJF = SJF[SJF['T(C)'] != -99]

#Averaging data to daily values
Model_Results_Daily = Model_Results.groupby(Model_Results['DATE TIME'].dt.date).mean().round(2)
Madera_TCD_Daily = Madera_TCD.groupby(Madera_TCD['DATE TIME'].dt.date).mean().round(2)
Old_Model_Daily = Old_Model.groupby(Old_Model['DATE TIME'].dt.date).mean().round(2)
FWQ_Daily = FWQ.groupby(FWQ['DATE TIME'].dt.date).mean().round(2)
SJF_Daily = SJF.groupby(SJF['DATE TIME'].dt.date).mean().round(2)

#Saving data to csv, dropping the JDAY and second DATE TIME columns
Model_Results_Daily = Model_Results_Daily.drop(['JDAY', 'DATE TIME'], axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
Madera_TCD_Daily = Madera_TCD_Daily.drop(['JDAY', 'DATE TIME'], axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
Old_Model_Daily = Old_Model_Daily.drop(['JDAY', 'DATE TIME'], axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
FWQ_Daily = FWQ_Daily.drop('DATE TIME', axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})
SJF_Daily = SJF_Daily.drop('DATE TIME', axis=1, errors='ignore').rename(columns={'T(C)': 'DAWT(C)'})

#Save to CSV after renaming
Model_Results_Daily.to_csv('Model_Results_Daily.csv')
Madera_TCD_Daily.to_csv('Madera_TCD_Daily.csv')
Old_Model_Daily.to_csv('Old_Model_Daily.csv')
FWQ_Daily.to_csv('FWQ_Daily.csv')
SJF_Daily.to_csv('SJF_Daily.csv')

#Calculating Stats

#Root Mean Square Error
Diff = Old_Model_Daily['DAWT(C)'] - FWQ_Daily['DAWT(C)']
RMSE = np.sqrt(np.mean(Diff**2))
print(f"RMSE: {RMSE}")  

#Plotting data
plt.rcParams['font.family'] = 'Arial'
#plt.plot(Model_Results_Daily['DAWT(C)'], label='New Model')
plt.plot(Old_Model_Daily['DAWT(C)'], label='Actual 2024 Flows')
plt.plot(Madera_TCD_Daily['DAWT(C)'], label='Madera TCD Scenario')
plt.plot(FWQ_Daily['DAWT(C)'], ls=':', label='FWQ Sensor')
#plt.plot(SJF_Daily['DAWT(C)'], label='SJF')
plt.ylim(5, 20)
plt.xlabel('Date', fontname='Arial')
plt.ylabel('Daily Avg. Water Temperature (°C)', fontname='Arial')
plt.title(f'RMSE: {RMSE:.2f}°C', fontname='Arial')
plt.legend()
plt.show()