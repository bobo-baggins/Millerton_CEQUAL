# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 15:48:22 2022
Inherited Novemeber 2023
Last update: July 2024

@author: jcohen
@editor: jHawkins
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import subprocess
import glob, os, shutil
import time
##################################  USER DEFINED VARIABLES #################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

#include file to update in excel: 
start_date = '6/2/2023 1:00' #start date of simulation (always start at 1:00)
end_date ='12/28/2023 23:00' #end date of simulation (always at at 23:00)
analog_years =  [2023]# available:  [1988,1989,1990,1994,2002,2007,2008,2013,2020,2021,2022,2023]
wait_time = 120 #seconds between simulations. Program will break if simulations take longer than wait time
#Code will excute after X amount of seconds, error if Model hasn't finished running by the end.
#Minimum wait time is 90 seconds.
make_output_folders = False #makes output folders, only if none exist
############################################################################################
############################################################################################
start_date = pd.to_datetime(start_date)
end_date = pd.to_datetime(end_date)
############################################################################################
############################################################################################

for analog_year in analog_years:
    print(analog_year)
    df_flow = pd.read_csv('flow_data/flow_data_base.csv', index_col = 0, parse_dates = True)
    df_temp = pd.read_csv('flow_data/flow_data_temp.csv', index_col = 0, parse_dates = True)
    df_flow = df_flow[(df_flow.index >= start_date) & (df_flow.index <= end_date)]
    df_temp = df_temp[(df_temp.index >= start_date) & (df_temp.index <= end_date)]

    JDAY_init = (pd.Timestamp(start_date)-pd.Timestamp('1-1-1921')).days + 1 #calulate initial Julinan Day for added data
    print(JDAY_init)
    JDAYS = np.arange(JDAY_init, JDAY_init + len(df_flow.index)) #create Julian Day array
    df_flow['JDAY'] = JDAYS #add Julian day array to flow dataframe

    
    ###################### create flow input files ###############################
    SPL_OUT = df_flow.SPL_OUT.values*0.028316847  # cfs to m3/s
    FKC_OUT = df_flow.FKC_OUT.values*0.028316847  # cfs to m3/s
    MC_OUT = df_flow.MC_OUT.values*0.028316847  # cfs to m3/s
    SJR_OUT = df_flow.SJR_OUT.values*0.028316847  # cfs to m3/s
    M_IN = df_flow.M_IN.values*0.028316847  # cfs to m3/s
    MIL_EVAP = df_flow.MIL_EVAP.values*-0.028316847  # cfs to m3/s
    JDAY = df_flow.JDAY.values *1.000
    Temp_predicted = df_temp['%s_Temp'%analog_year].values *1.000
    zero_filler = np.zeros(len(df_flow.index))*1.000 #fill-in for qin br2-4
    
    for i in range(0,len(JDAY)): #rounding to hundredths place
        JDAY[i] = '%0.2f'%JDAY[i]
        SPL_OUT[i] = '%0.2f'%SPL_OUT[i]
        FKC_OUT[i] = '%0.2f'%FKC_OUT[i]
        MC_OUT[i] = '%0.2f'%MC_OUT[i]
        SJR_OUT[i] = '%0.2f'%SJR_OUT[i]
        MIL_EVAP[i] = '%0.2f'%MIL_EVAP[i]
        M_IN[i] = '%0.2f'%M_IN[i]
        Temp_predicted[i] = '%0.2f'%Temp_predicted[i]
        zero_filler[i] = '%0.2f'%zero_filler[i]
    # changed 0.99 to 1.99 as this was causing an error during model run with input files being written missing the last day.
    JDAY_end = JDAY[-1]+1.99
    print(JDAY_end)
    ## create qot_br1.npt file
    with open('initial_files/mqot_br1_init.npt',"r") as f:
        lines = f.readlines()
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{SPL_OUT[i] : >8}{FKC_OUT[i] : >8}{MC_OUT[i] : >8}{SJR_OUT[i] : >8}" #make line for CEQUAL timestep
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}{0.0 : >8}{0.0 : >8}{0.0 : >8}\n")
    
    with open('mqot_br1.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    
    ##create mqdt_br1.npt
    with open('initial_files/mqdt_br1_init.npt',"r") as f:
        lines = f.readlines()
    
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{MIL_EVAP[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    
    with open('mqdt_br1.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    
    ## create mqin_br1.npt
    with open('initial_files/mqin_br1_init.npt',"r") as f:
        lines = f.readlines()
    
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{M_IN[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    
    with open('mqin_br1.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    
    ### create mqin_br2-4
    with open('initial_files/mqin_br2_init.npt',"r") as f:
        lines = f.readlines()
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{zero_filler[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    with open('mqin_br2.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    with open('initial_files/mqin_br3_init.npt',"r") as f:
        lines = f.readlines()
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{zero_filler[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    with open('mqin_br3.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    with open('initial_files/mqin_br4_init.npt',"r") as f:
        lines = f.readlines()
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{zero_filler[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    with open('mqin_br4.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    ### create mtin_br1-4
    with open('initial_files/mtin_br1_init.npt',"r") as f:
        lines = f.readlines()
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{Temp_predicted[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    with open('mtin_br1.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    with open('initial_files/mtin_br2_init.npt',"r") as f:
        lines = f.readlines()
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{Temp_predicted[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    with open('mtin_br2.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    with open('initial_files/mtin_br3_init.npt',"r") as f:
        lines = f.readlines()
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{Temp_predicted[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    with open('mtin_br3.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    with open('initial_files/mtin_br4_init.npt',"r") as f:
        lines = f.readlines()
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{Temp_predicted[i] : >8}"
        lines.append(l)
    lines.append(f"\n{JDAY_end : >8}{0.0 : >8}\n")
    with open('mtin_br4.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
    ###############################convert met data##########################
    df_met = pd.read_csv('met_data/%s_CEQUAL_met_inputs.csv'%analog_year,index_col = 0, parse_dates= True)
    df_met = df_met[(df_met.index >= start_date) & (df_met.index <= end_date)]
    #This creates the JDAY column for the met data. Commenting out for testing due to long run time
    hour_cycle = [0,0.04,0.08,0.13,0.17,0.21,0.25,0.29,0.34,0.38,0.42,0.46,0.50,0.55,0.59,0.63,0.67,0.71,0.76,0.80,0.84,0.88,0.92,0.97]
    JDAYS_hourly = [] #empty array for now Jdates added to original
    for d in range(JDAY_init, JDAY_init + len(df_flow.index)+1):
        for h in hour_cycle:
            JDAYS_hourly.append(d+h)
    JDAYS_hourly.pop(0)
    df_met['JDAY'] = JDAYS_hourly
    
    JDAY = df_met.JDAY.values*1.000  # cfs to m3/s
    TAIR = df_met.TAIR.values*1.000  # cfs to m3/s
    TDEW = df_met.TDEW.values*1.000  # cfs to m3/s
    WIND = df_met.WIND.values*1.000  # cfs to m3/s
    PHI = df_met.PHI.values*1.000  # cfs to m3/s
    CLOUD = df_met.CLOUD.values*1.000  # cfs to m3/s
    SRO = df_met.SRO.values *1.000
    
    for i in range(0,len(JDAY)): #rounding to hundredths place
        JDAY[i] = '%0.2f'%JDAY[i]
        TAIR[i] = '%0.2f'%TAIR[i]
        TDEW[i] = '%0.2f'%TDEW[i]
        WIND[i] = '%0.2f'%WIND[i]
        PHI[i] = '%0.2f'%PHI[i]
        CLOUD[i] = '%0.2f'%CLOUD[i]
        SRO[i] = '%0.2f'%SRO[i]
    
    with open('initial_files/mmet3_init.npt',"r") as f:
        lines = f.readlines()
    
    for i in range(0, len(JDAY)):
        l = f"\n{JDAY[i] : >8}{TAIR[i] : >8}{TDEW[i] : >8}{WIND[i] : >8}{PHI[i] : >8}{CLOUD[i] : >8}{SRO[i] : >8}" #make line for CEQUAL timestep
        lines.append(l)
    lines.append("\n")
    
    with open('mmet3.npt',"w") as update:
        update.writelines(lines)
    update.close()
    
       
    os.system("start w2_cvf_millerton_2_1_07.exe")
    time.sleep(wait_time)
    source_dir = r'../NewExNoAction2022'
    dest_dir =  r'../CEQUAL_outputs/%s'%analog_year  
    if make_output_folders == True:         
        os.mkdir(dest_dir)

    opt_files = glob.iglob(os.path.join("*.opt"))
    for file in opt_files:
        if os.path.isfile(file):
            shutil.copy2(file, dest_dir)      
            
    npt_files = glob.iglob(os.path.join("*.npt"))
    
    for file in npt_files:
        if os.path.isfile(file):
            shutil.copy2(file, dest_dir)      
        # os.mkdir(dest_dir)

    with open('../CEQUAL_outputs/%s/str_br1.opt'%analog_year,"r") as f:
        lines = f.readlines()
    lines.pop(0)
    lines.pop(0)
    for i,l in enumerate(lines):
        lines[i] = np.array(list(filter(None, np.array(l.split(' ')))), dtype = np.float32)
    ix = pd.date_range(start = start_date, end =end_date ,freq = 'h')

    #print("ix value:")
    #print(len(ix))
    #print(ix[0])
    # lengths between lines and ix are off by an hour.... or a value of 1.
    # Lines is missing an hour, so ix is clipped to match.
    ix=ix[1:]

    df = pd.DataFrame(index = ix, columns = ['JDAY', 'SPL_temp_C', 'FKC_temp_C', 'MC_temp_C', 'SJR_temp_C',
                                  'SPL_Q_m3s', 'FKC_Q_m3s', 'MC_Q_m3s', 'SJR_Q_m3s',
                                  'SPL_ELEVCL_m', 'FKC_ELEVCL_m', 'MC_ELEVCL_m', 'SJR_ELEVCL_m'])
    df[df.columns] = lines
    # dest_dir = r'../output_csvs/%s-analog-%s-exceedance'%(analog_year,exceedance)
    # os.mkdir(dest_dir)
    df.to_csv('../output_csvs/outflow_temps/%s_temp_analog-release_outputs.csv'%(analog_year))

    # can't locate parent directory for whatever reason... my virtual environment?
    #This is a quick fix..  .   .
   # df.to_csv('%s_analog-release_outputs.csv'%(analog_year))
