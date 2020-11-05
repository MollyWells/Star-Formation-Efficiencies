#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 12:47:13 2020

@author: mollywells
"""

#Modelling2.0

import math as m
import numpy as np
import matplotlib.pyplot as plt
import random as r
import csv
from scipy import stats
import numpy as np

# IMPORT INITIAL MASS FUNCTION 
# Assign stellar masses (range 0.1 to 120 Msun), logarithmically spaced 

sStellar_mass = []
x = 0.1

for n in range(5000):

    sStellar_mass.append(x)

    x = x*(10**0.003)

aStellar_mass = np.asarray(sStellar_mass)



# Assign stellar mass probabilities in accordance to Kroupa (2001) IMF
aC_Prob = np.genfromtxt('IMFC_Prob.csv', delimiter = ',', usecols = 2)


# IMPORT CLUMP MASS VALUES 
aClmp_M = np.genfromtxt('CMasses.csv', delimiter = ',', usecols = 2)

SFE_all = []
#START OF LOOP

runs = 1000

for k in range(len(aClmp_M)):

    

# Create empty strings for storing cluster properties, i.e Cluster 

#   luminosities, Star formation Efficiencies (SFEs), etc.

    sSFE = []   
    sClstr_L, sClstr_M, sLyFlx = [], [], []
    sMost_M, sMost_L, sMost_LyFlx = [], [], []        
    sTMC, sClstr_pop , sAvStar_M = [], [], []         
    Total_SFE, Total_Clstr_L, Total_AvStar_M = 0.0, 0.0, 0.0        
    Total_LyFlx, Total_Clstr_pop, Total_Most_M = 0.0, 0.0, 0.0

#==============================================================================

#===================== START of loop around each CLUSTER ======================

#==============================================================================



    for j in range(runs):

        

        # Define initial cluster properties

        Clstr_M, Clstr_L, LyFlx = 0.0, 0.0, 0.0     
        Prev_M, Prev_L = 0.0, 0.0                           
        Most_M, Most_L, Most_LyFlx, TMC = 0.1, 0.0, 0.0, 0.0

        Clstr_pop = 0
        

        

#==============================================================================

#======================= START of loop adding each STAR =======================
       
        
        print('running')
        counter = 0

        while (Clstr_L <= (Clstr_M**1.1849)):
            
            if counter == runs:
                break
            else: 
        
                Prev_M = Clstr_M
                Prev_L = Clstr_L   

            

            # Generate a star using random number generation and stellar mass 

            #   probabilities

                Random_number = r.random()         

                Index = np.argmax(np.where(aC_Prob < Random_number))     #maxvalue

                Star_M = aStellar_mass[Index]                        

            

            # Increase total cluster population

                Clstr_pop = Clstr_pop + 1

            

            # Calculate new cluster mass

                Clstr_M = Clstr_M + Star_M

            

            # Define constants for calculating cluster properties

                a = 3.94613-0.0242550*Star_M+0.000249373*Star_M**2-9.52930e-07*Star_M**3

                b = 39.0353+0.853175*Star_M-0.0255907*Star_M**2+0.000263813*Star_M**3

                c = 47.6503+0.0540930*Star_M-0.000517299*Star_M**2+1.80882e-06*Star_M**3

            

            # Calculate new cluster luminousity

                if (Star_M < 1):      

                        Clstr_L = Clstr_L+Star_M**4    

                else:                            

                        Clstr_L = Clstr_L+Star_M**(a) 


                try:
                    logClstr_L = m.log10(Clstr_L)
                except ValueError:
                    pass
                    
                logaClmp_M = m.log10(aClmp_M[k])

            

                if (Prev_L == 0):

                    logPrev_L = np.NaN

                else:    

                    logPrev_L = m.log10(Prev_L)
                
                  

                if (Prev_M > aClmp_M[k]*0.99):

                    break

                    Clstr_pop = Clstr_pop - 1

                    Clstr_M = Prev_M

                    Clstr_L = Prev_L 
                    
               
                counter += 1

        print ('Cluster', j, Clstr_L, Clstr_M)
        # Calculate the Star Formation Efficiency (SFE)

        SFE = (Clstr_M/aClmp_M[k])*100
                
        # Calculate Average Stellar Mass
        
        AvStar_M = Clstr_M/Clstr_pop 
        
        # Sum of cluster SFEs and luminosities etc. for calculation of average values 

        Total_SFE = Total_SFE + SFE  

        Total_Clstr_L = Total_Clstr_L + Clstr_L

        Total_AvStar_M = Total_AvStar_M + AvStar_M      
                
        Total_Clstr_pop = Total_Clstr_pop + Clstr_pop

        Total_Most_M = Total_Most_M + Most_M
            

        # Record cluster properties such as SFE and cluster luminosity, etc.

        sSFE.append(SFE)

        sClstr_L.append(Clstr_L)  
        
        sClstr_M.append(Clstr_M)      
            
        sMost_M.append(Most_M)    
    
        sMost_L.append(Most_L)      

        sClstr_pop.append(Clstr_pop)     
       
              
    SFE_all.append(sSFE)  
    

###############################################################################

####################### END of loop around each CLUSTER #######################

###############################################################################



    # Convert strings to arrays              
#
#    aSFE = np.asarray(sSFE)                   
#
#    aClstr_L = np.asarray(sClstr_L) 
#
#    aClstr_M = np.asarray(sClstr_M)   
#
#    aLyFlx = np.asarray(sLyFlx)         
#
#    aMost_M = np.asarray(sMost_M)              
#
#    aMost_L = np.asarray(sMost_L)           
#
#    aMost_LyFlx = np.asarray(sMost_LyFlx)   
#
#    aTMC = np.asarray(sTMC)               
#
#    aClstr_pop = np.asarray(sClstr_pop)
#
#    aAvStar_M = np.asarray(sAvStar_M)



                                                     

#---------------------------------------------------------------------  

# Calculate/store average values and edit arrays for plotting



    # Store average values for average lines on plots

#    AvL = Total_Clstr_L/runs
#
#    sAvL.append(AvL)
#
#    AvS = Total_SFE/runs
#
#    sAvS.append(AvS) 
#
#    AvA = Total_AvStar_M/runs
#
#    sAvA.append(AvA) 
#
#    AvN = Total_LyFlx/runs
#
#    sAvN.append(AvN) 
#
#    Avi = Total_Clstr_pop/runs
#
#    sAvi.append(Avi) 
#
#    AvMM = Total_Most_M/runs
#
#    sAvMM.append(AvMM)

    

    # Sort the arrays (size order)

#    sortaL = np.sort(aClstr_L)
#
#    sortaSFE = np.sort(aSFE)
#
#    sortaMM = np.sort(aMost_M)
#
#    sortCp = np.sort(aClstr_pop)
#
#    sortaAvM = np.sort(aAvStar_M)
#
#    sortaLy = np.sort(aLyFlx)
#
#
#
#    # Trim arrays to desired band size:
#
#    # Luminosity
#
#    midaLa = stats.trimboth(sortaL, 0.1)
#
#    midaLb = stats.trimboth(sortaL, 0.2)
#
#    midaLc = stats.trimboth(sortaL, 0.3)
#
#    midaLd = stats.trimboth(sortaL, 0.4)
#
#    # SFE
#
#    aSFEa = stats.trimboth(sortaSFE, 0.1)
#
#    aSFEb = stats.trimboth(sortaSFE, 0.2)
#
#    aSFEc = stats.trimboth(sortaSFE, 0.3)
#
#    aSFEd = stats.trimboth(sortaSFE, 0.4)   
#
#    # Cluster population
#
#    midCpa = stats.trimboth(sortCp, 0.1)
#
#    midCpb = stats.trimboth(sortCp, 0.2)
#
#    midCpc = stats.trimboth(sortCp, 0.3)
#
#    midCpd = stats.trimboth(sortCp, 0.4)
#
#    # Average stellar mass
#
#    midaAvMa = stats.trimboth(sortaAvM, 0.1)
#
#    midaAvMb = stats.trimboth(sortaAvM, 0.2)
#
#    midaAvMc = stats.trimboth(sortaAvM, 0.3)
#
#    midaAvMd = stats.trimboth(sortaAvM, 0.4)    
#
#    # Lyman Flux
#
#    midaLya = stats.trimboth(sortaLy, 0.1)
#
#    midaLyb = stats.trimboth(sortaLy, 0.2)
#
#    midaLyc = stats.trimboth(sortaLy, 0.3)
#
#    midaLyd = stats.trimboth(sortaLy, 0.4)
#
#    # Most Massive Star
#
#    midaMMa = stats.trimboth(sortaMM, 0.1)
#
#    midaMMb = stats.trimboth(sortaMM, 0.2)
#
#    midaMMc = stats.trimboth(sortaMM, 0.3)
#
#    midaMMd = stats.trimboth(sortaMM, 0.4)

    




###############################################################################

###############################################################################

##################### END of loop around each CLUMP MASS ######################

###############################################################################

###############################################################################

 #BOXPLOT SFE vs INITIAL CLUMP MASS
    
data1 = SFE_all[0] 
data2 = SFE_all[1] 
data3 = SFE_all[2] 
data4 = SFE_all[3] 
data5 = SFE_all[4] 
data6 = SFE_all[5] 
data7 = SFE_all[6] 
data8 = SFE_all[7] 
data9 = SFE_all[8] 
data10 = SFE_all[9] 
data11 = SFE_all[10] 
data12 = SFE_all[11]    
    
data = [data1,data2,data3,data4,data5,data6,data7,data8]

plt.boxplot(data, positions = [1,2,3,4,5,6,7,8], showfliers = False)
plt.xticks(np.array([1,2,3,4,5,6,7,8]), (125, 150, 2,3,4,5,6,7))
# show plot 
plt.show() 
    
#2D histogram  


#plt.cla()
x = [1,2,3,4,5,6,7,8]    
for xe, ye in zip(x, data):
    plt.hist2d([xe] * len(ye), ye, bins=(8)-0.5, cmap='Blues')
    
rwidth = 1    
plt.xticks(x)
plt.axes().set_xticklabels(['1', '2','3','4','5','6','7','8'])
plt.show()   
#plt.colorbar()
plt.show()    
    
    
    