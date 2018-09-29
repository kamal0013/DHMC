## This script calibrates the parameters of DHMC 
##...and generates multiple (say 1000) stochastic rainfall series.
## Change the filename of observed data and user inputs as required.
## Calibrated parameters and simulated rainfall can be saved in text files.

import numpy as np
import scipy.stats as st
import random
import itertools as itr

####=== Data and User Inputs ===####
stn = '66062' ##Site or Station ID

data = np.loadtxt("data_site_"+str(stn)+".txt") ## Observed daily rainfall timeseries
row, col = data.shape
print ('Rows-Columns in Data', row, col)

### User Inputs 
ths = 0.3  ##Threshold rainfall (mm)
syear = 1979 ##Start year
nyear = 30 ##Data period in number of years
ndec = int(nyear/10) ##Number of decades
nrun = 100 ##Number of replicates to be generated

### Want to save parameters and simulated rainfall or not
sav_decadal_MC_params = 'y' ##put 'y' to save decade-varied MC parameter values
sav_yearly_gamma_params = 'y' ##put 'y' to save yearly-varied gamma parameter values
sav_simulated_rain = 'y' ##put 'y' to save Simulated Rainfall

######====End of User Inputs=====##########



print ('Run DHMC for Station', stn)




##############======= Nothing below Should be Changed ===================##########
##############=======DHMC Calibration ==================##############
##==== Split into Decadal Data ====##
data1 = np.array_split(data, ndec, axis=0)

##=== Calibration of Decadal MC Parameters ===##
deca_dry = [[] for n in range(ndec)]
deca_wet = [[] for n in range(ndec)]

for decade in range(0,ndec):    
    dec_data = data1[decade]

    d2d = np.zeros(12)
    d2w = np.zeros(12)
    w2w = np.zeros(12)
    w2d = np.zeros(12)
    dry_dry = np.zeros(12)
    wet_wet = np.zeros(12)
    for i in range (0,len(dec_data)):
        for m in range(0,12):
            if dec_data[i,1] == m+1:            
                if dec_data[i-1,3] < ths:
                    if dec_data[i,3] < ths:            
                        d2d[m] += 1.0
                    else:
                        d2w[m] += 1.0
                if dec_data[i-1,3] >= ths:
                    if dec_data[i,3] < ths:           
                        w2d[m] += 1.0
                    else:
                        w2w[m] += 1.0

    for m in range(0,12):
        dry_dry[m] = d2d[m]/(d2d[m] + d2w[m])
        wet_wet[m] = w2w[m]/(w2d[m] + w2w[m])

    deca_dry[decade] = dry_dry
    deca_wet[decade] = wet_wet

dry = [(1- deca_wet[0][0])/(2- deca_dry[0][0] - deca_wet[0][0])]

####====Print and/or Save Decadal MC Parameters=====####
if sav_decadal_MC_params == 'y':
    np.savetxt("Param_MC_Dec_Dry_Site_"+str(stn)+".txt", (deca_dry), fmt='%0.3f')
    np.savetxt("Param_MC_Dec_Wet_Site_"+str(stn)+".txt", (deca_wet), fmt='%0.3f')


####====Calculate Gamma Parameters=======####
## Fitting Log-normal Parameters of Mean and Std of Wet-day Rainfall Depths ##
rain_samp = [[[] for i in range(nyear)] for i in range(12)]
samp_size = np.zeros((12,nyear))
mean_samp = np.zeros((12, nyear))
std_samp = np.zeros((12, nyear))

log_mean = np.zeros((12, nyear))
log_std = np.zeros((12, nyear))

for m in range(0,12):
    for y in range(0,nyear):
        for i in range (0,len(data)):          
            if data[i,1] == m+1 and data[i,0] == y+syear:
                if data[i,3] >= ths:
                    rain_samp[m][y].append(data[i,3])                    
                    samp_size[m][y] += 1
            
        if samp_size[m][y] > 2:
            mean_samp[m][y] = np.mean(rain_samp[m][y])        
            std_samp[m][y] = np.std(rain_samp[m][y])

        if mean_samp[m][y] > 0:
            log_mean[m][y]= np.log(mean_samp[m][y])
        if std_samp[m][y] > 0:
            log_std[m][y]= np.log(std_samp[m][y])

mu_mean = [np.mean(x) for x in log_mean]
mu_std = [np.mean(x) for x in log_std]
sigma_mean = [np.std(x) for x in log_mean]
sigma_std = [np.std(x) for x in log_std]

r = np.zeros(12)
for m in range(0,12):
    r[m] = np.corrcoef(log_mean[m],log_std[m])[1,0]

if sav_yearly_gamma_params == 'y':
    np.savetxt("Param_Gamma_Mean_Yearly_Site_"+str(stn)+".txt", mean_samp, fmt='%0.2f')
    np.savetxt("Param_Gamma_SD_Yearly_Site_"+str(stn)+".txt", std_samp, fmt='%0.2f')
    np.savetxt("Param_Gamma_Corr_Site_"+str(stn)+".txt", r, fmt='%0.2f')


#######==============DHMC Simulation===================#########
month_days = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
simrain = [[] for n in range(nrun)] ##List to store simulated rainfall

for n in range(0,(10+nrun)): ##Discard initial 10 runs
    mean1 = np.zeros(12)
    std1 = np.zeros(12)    
    beta1 = np.zeros(12)
    alpha1 = np.zeros(12)
    
    for decade in range(0,ndec):    
        for y in range(0,10):
            for m in range(0,12):
                mean1[m] = random.lognormvariate(mu_mean[m],sigma_mean[m])
                std1[m] = random.gauss((mu_std[m] + ((sigma_std[m] * r[m])*(mean1[m] - mu_mean[m])/sigma_mean[m])), np.sqrt(1 - pow(r[m],2))*(sigma_mean[m]))
                    
                beta1[m] = pow((std1[m]),2) / mean1[m]
                alpha1[m] = mean1[m] / beta1[m]              

########======Generate Rainfall Series========#########
                if (((y+decade*10)%4 == 0 and (y+decade*10)%100 != 0) or (y+decade*10)%400 == 0) and m==1: 
                    random_num = np.random.rand(29)
                    #print ('leapyear', y, m)
                else:
                    random_num = np.random.rand(month_days[m])
                        
                for day in range(0,len(random_num)):
                    if y==0 and m==0 and day==0: ##First day of the series
                        if random_num[0] <= dry: 
                            today_rain = 0                    
                            dry1 = True
                        else:                        
                            today_rain = max(ths, random.gammavariate(alpha1[0], beta1[0]))                            
                            dry1 = False

                    ##For rest of the days     
                    if dry1 is True: 
                        if random_num[day] <= deca_dry[decade][m]:
                            today_rain = 0                   
                            dry1 = True
                        else:                        
                            today_rain = max(ths, random.gammavariate(alpha1[m], beta1[m]))                            
                            dry1 = False
              
                    else:
                        if random_num[day] <= deca_wet[decade][m]:                        
                            today_rain = max(ths, random.gammavariate(alpha1[m], beta1[m]))                           
                            dry1 = False    
                        else:
                            today_rain = 0                   
                            dry1 = True

                    ##start appending simrain after first 10 runs
                    if n >= 10:
                        simrain[n-10].append(np.round(today_rain, decimals=2))

### Check Progress of Simulation        
    if n==0:
        print ('Simulation in Progress')
    if n==nrun/4:
        print ('25 Percent Completed')        
    if n==nrun/2:
        print ('50 Percent Completed')
    if n==3*nrun/4:
        print ('75 Percent Completed')
    if n==nrun:
        print ('100 Percent Completed')


simrain_array = np.asarray(simrain)
#print (simrain_array)
simrain_dated = np.vstack((data[:,0],data[:,1],data[:,2],simrain_array))
#print (simrain_dated)
print ('Columns and Rows in Simulated File', simrain_dated.shape)

if sav_simulated_rain == 'y':
    simrain_zip = list(zip(*simrain_dated))
    np.savetxt("DHMC_Simrain_Site_"+str(stn)+".txt", simrain_zip, fmt='%-8s', header = 'Year   Month     Day      SimRains')
