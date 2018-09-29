### As a sample test of model performance, this script calculates
###...Z Scores for rainfall depth at daily, monthly and multiyear resolutions 
###...and saves in textfiles.
### Change the filenames of observed and simulated data, and user inputs as required.


import numpy as np
import scipy.stats as st
import itertools as itr

####=== Data and User Inputs ===####
stn = '23034' ##Site or Station ID

data = np.loadtxt("data_site_"+str(stn)+".txt")  ## Observed daily rainfall timeseries
sim_data = np.loadtxt("DHMC_Simrain_Site_"+str(stn)+".txt") ## DHMC replicates of daily rainfall

ths = 0.3  ##Threshold rainfall (mm)
syear = 1979 ##Start year
nyear = 30 ##Data period in number of years
nrun = 100 ## Number of replicates in sim_data




##============= Calculate Observed Statistics ===========##
## Daily Variables
obs_d_rain = [[]*n for n in range(12)]
obs_d_mean = np.zeros(12)
obs_d_std = np.zeros(12)

## Monthly Variables
obs_m_rain = np.zeros((12,nyear))
obs_m_mean = np.zeros(12)
obs_m_std = np.zeros(12)

## Multi-year Variables
ob_arain = np.zeros(nyear)
obs_y_rains = np.zeros((10,nyear))
obs_y_mean = np.zeros(10)
obs_y_std = np.zeros(10)

for y in range(0,nyear):
    for m in range(0,12):
        for i in range (0,len(data)):
            if data[i,0]==y+syear and data[i,1]==m+1 and data[i,3]>=ths:
                obs_d_rain[m].append(data[i,3])
                obs_m_rain[m][y] += (data[i,3])
                ob_arain[y] += data[i,3]
## Observed Daily-Monthly 
for m in range(0,12):
    obs_d_mean[m] = np.mean(obs_d_rain[m])
    obs_d_std[m] = np.std(obs_d_rain[m])    
    obs_m_mean[m] = np.mean(obs_m_rain[m])
    obs_m_std[m] = np.std(obs_m_rain[m])

## Observed Multi-year
for k in range(0,nyear):
    for v in range(0,10):
        if v == 0:
            obs_y_rains[0][k] = ob_arain[k]
        if v > 0:
            obs_y_rains[v][k] = obs_y_rains[v-1][k] + ob_arain[k - v]        

for y in range(0,10):
    obs_y_mean[y] = np.mean(obs_y_rains[y])       
    obs_y_std[y] = np.std(obs_y_rains[y])


#######========= Calculate Simulated Statistics ==============#########
## Daily Variables
d_mean = [[]*n for n in range(12)]
d_std = [[]*n for n in range(12)]

## Monthly Variables
m_mean = [[]*n for n in range(12)]
m_std = [[]*n for n in range(12)]

## Multi-year Variables
y_mean = [[]*n for n in range(10)]
y_std = [[]*n for n in range(10)]

for n in range(0,nrun): #
    ###Variable for Amount Statistic
    rain = [[]*m for m in range(12)]
    mon_rain = np.zeros((12,nyear))
    ann_rain = np.zeros(nyear)
    yealy_stats = np.zeros((10,nyear))
    
    for y in range(0,nyear):
        for m in range(0,12):
            for i in range (0,len(sim_data)):
                if sim_data[i,0]==y+syear and sim_data[i,1]==m+1 and sim_data[i,(n+3)]>=ths:                   
                    rain[m].append(sim_data[i,(n+3)])
                    mon_rain[m][y] += (sim_data[i,(n+3)])
                    ann_rain[y] += (sim_data[i,(n+3)])

    ##Daily-Monthly Amount Stats
    for f in range(0,12):
        d_mean[f].append(np.mean(rain[f]))
        d_std[f].append(np.std(rain[f]))        
        m_mean[f].append(np.mean(mon_rain[f]))
        m_std[f].append(np.std(mon_rain[f]))       


    ##Multi-year Amount Stats
    for k in range(0,nyear):
        for v in range(0,10):
            if v==0:                
                yealy_stats[0][k] = ann_rain[k]
            if v > 0:
                yealy_stats[v][k] = yealy_stats[v-1][k] + ann_rain[k-v]
       
    for w in range(0,10):
        y_mean[w].append(np.mean(yealy_stats[w]))
        y_std[w].append(np.std(yealy_stats[w]))


###=========Simulated depth statistics=================###
#####Daily-Monthly
m_d_mn = [np.mean(s) for s in d_mean]
s_d_mn = [np.std(s) for s in d_mean]
m_d_sd = [np.mean(s) for s in d_std]
s_d_sd = [np.std(s) for s in d_std]

m_m_mn = [np.mean(s) for s in m_mean]
s_m_mn = [np.std(s) for s in m_mean]
m_m_sd = [np.mean(s) for s in m_std]
s_m_sd = [np.std(s) for s in m_std]

#####Multi-Year
m_y_mn = [np.mean(s) for s in y_mean]
s_y_mn = [np.std(s) for s in y_mean]
m_y_sd = [np.mean(s) for s in y_std]
s_y_sd = [np.std(s) for s in y_std]       

  
#######================= Z Score Calculation ============================########
####Daily and Monthly
z_d_mean = np.zeros(12)
z_d_std = np.zeros(12)
z_m_mean = np.zeros(12)
z_m_std = np.zeros(12)

for m in range(0,12):
    z_d_mean[m] = (obs_d_mean[m] - m_d_mn[m])/s_d_mn[m]
    z_d_std[m] = (obs_d_std[m] - m_d_sd[m])/s_d_sd[m]
    z_m_mean[m] = (obs_m_mean[m] - m_m_mn[m])/s_m_mn[m]
    z_m_std[m] = (obs_m_std[m] - m_m_sd[m])/s_m_sd[m]

####Multi-Year
z_y_mean = np.zeros(10)
z_y_std = np.zeros(10)

for y in range(0,10):
    z_y_mean[y] = (obs_y_mean[y] - m_y_mn[y])/s_y_mn[y]
    z_y_std[y] = (obs_y_std[y] - m_y_sd[y])/s_y_sd[y]

np.savetxt("ZScores_Depth_DailyMonthly_Site_"+str(stn)+".txt", (z_d_mean,z_d_std,z_m_mean,z_m_std), fmt='%0.2f',
           header='R1=Z(Daily Mean),R2=Z(Daily SD),R3=Z(Monthly Mean),R4=Z(Monthly SD)')
np.savetxt("ZScores_Depth_Multiyear_Site_"+str(stn)+".txt", (z_y_mean,z_y_std), fmt='%0.2f',header='R1=Z(Multiyear Mean),R2=Z(Multiyear SD)')
