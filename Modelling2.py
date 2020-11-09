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
data13 = SFE_all[12] 
data14 = SFE_all[13] 
data15 = SFE_all[14] 
data16 = SFE_all[15] 
data17 = SFE_all[16] 
data18 = SFE_all[17] 
data19 = SFE_all[18] 
data20 = SFE_all[19] 
data21 = SFE_all[20]    
data_all = np.concatenate([data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17,data18,data19,data20,data21])
data_all1 = np.concatenate([data1,data2,data3,data4,data5,data6])
data_all2 = np.concatenate([data7,data8,data9,data10,data11])
data_all3 = np.concatenate([data12,data13,data14,data15,data16])
data_all4 = np.concatenate([data17,data18,data19,data20,data21])



x1 = [2] * 1000
x2 = [2.1] * 1000
x3 = [2.2] * 1000
x4 = [2.3] * 1000
x5 = [2.4] * 1000
x6 = [2.5] * 1000
x7 = [2.6] * 1000
x8 = [2.7] * 1000
x9 = [2.8] * 1000
x10 = [2.9] * 1000
x11 = [3] * 1000
x12 = [3.1] * 1000
x13 = [3.2] * 1000
x14 = [3.3] * 1000
x15 = [3.4] * 1000
x16 = [3.5] * 1000
x17 = [3.6] * 1000
x18 = [3.7] * 1000
x19 = [3.8] * 1000
x20 = [3.9] * 1000
x21 = [4] * 1000

x_all = np.concatenate([x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21])
x_all1 = np.concatenate([x1,x2,x3,x4,x5,x6])
x_all2 = np.concatenate([x7,x8,x9,x10,x11])
x_all3 = np.concatenate([x12,x13,x14,x15,x16])
x_all4 = np.concatenate([x17,x18,x19,x20,x21])

data = [data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17,data18,data19,data20,data21]
plt.boxplot(data, positions = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21], showfliers = False)
plt.xticks(np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]), (2, 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4))
plt.xlabel('Log(Mass)')
plt.ylabel('SFE (%)')
plt.show() 



x = [2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4]    
xa = [2,2.1,2.2,2.3,2.4,2.5]
xb = [2.6,2.7,2.8,2.9,3]
xc = [3.1,3.2,3.3,3.4,3.5]
xd = [3.6,3.7,3.8,3.9,4]
#for xe, ye in zip(x, data):
    #plt.hist2d([xe] * len(ye), data, bins=[21,150] , cmap='Blues')
    
plt.hist2d(x_all,data_all, bins =[21,150], cmap = 'Blues' )  
plt.xticks(x)
plt.ylim(0,25)
#plt.axes().set_xticklabels(['1', '2','3','4','5','6','7','8'])
plt.colorbar()
plt.xlabel('Log(Mass)')
plt.ylabel('SFE (%)')
plt.show()   
    
plt.hist2d(x_all1,data_all1, bins =[6,150], cmap = 'Blues' )  
plt.xticks(xa)
plt.ylim(0.01,15)
#plt.axes().set_xticklabels(['1', '2','3','4','5','6','7','8'])
plt.colorbar()
plt.xlabel('Log(Mass)')
plt.ylabel('SFE (%)')
plt.show() 

plt.hist2d(x_all2,data_all2, bins =[5,150], cmap = 'Blues' )  
plt.xticks(xb)
plt.ylim(0.01,5)
#plt.axes().set_xticklabels(['1', '2','3','4','5','6','7','8'])
plt.colorbar()
plt.xlabel('Log(Mass)')
plt.ylabel('SFE (%)')
plt.show() 

plt.hist2d(x_all3,data_all3, bins =[5,200], cmap = 'Blues' )  
plt.xticks(xc)
plt.ylim(0.01,1)
#plt.axes().set_xticklabels(['1', '2','3','4','5','6','7','8'])
plt.colorbar()
plt.xlabel('Log(Mass)')
plt.ylabel('SFE (%)')
plt.show() 

plt.hist2d(x_all4,data_all4, bins =[5,200], cmap = 'Blues' )  
plt.xticks(xd)
plt.ylim(0.01,0.5)
#plt.axes().set_xticklabels(['1', '2','3','4','5','6','7','8'])
plt.colorbar()
plt.xlabel('Log(Mass)')
plt.ylabel('SFE (%)')
plt.show() 
    