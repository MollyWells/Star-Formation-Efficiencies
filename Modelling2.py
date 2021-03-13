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
import scipy as sp
from datetime import datetime

# IMPORT INITIAL MASS FUNCTION 
# Assign stellar masses (range 0.1 to 120 Msun), logarithmically spaced 

sStellar_mass = []
x = 0.1

for n in range(1026):

    sStellar_mass.append(x)

    x = x*(10**0.003)

aStellar_mass = np.asarray(sStellar_mass)

# Assign stellar mass probabilities in accordance to Kroupa (2001) IMF
aC_Prob = np.genfromtxt('IMFC_Prob.csv', delimiter = ',', usecols = 2)

# IMPORT CLUMP MASS VALUES 
aClmp_M = np.genfromtxt('CMasses.csv', delimiter = ',', usecols = 2)



SFE_all = []
#START OF LOOP
debug = True
Tol_L = 0.1
runs = 1000
for k in range(len(aClmp_M)):
    print(f'Processing Clump Mass {aClmp_M[k]:0.2f}' )
    

# Create empty strings for storing cluster properties, i.e Cluster 

#   luminosities, Star formation Efficiencies (SFEs), etc.

    sSFE = []   
    sClstr_L, sClstr_M = [], []         
    Total_SFE, Total_Clstr_L, Total_AvStar_M = 0.0, 0.0, 0.0        
   

#==============================================================================

#===================== START of loop around each CLUSTER ======================

#==============================================================================


    count_Ex_M = 0
    count_Ex_L = 0
    for j in range(runs):

        

        # Define initial cluster properties

        Clstr_M, Clstr_L, tot_massive_L, Star_L = 0.0, 0.0, 0.0, 0.0
        Prev_M, Prev_L = 0.0, 0.0                           
        Most_M, Most_L = 0.1, 0.0

        Clstr_pop = 0
        

        

#==============================================================================

#======================= START of loop adding each STAR =======================
       
        
        #print('running')
        counter = 0
        
        while True:
            
            #if counter == runs:
               # break
            #else: 
        
            Prev_M = Clstr_M
            Prev_L = Clstr_L   

            

            # Generate a star using random number generation and stellar mass probabilities

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

                    Star_L = Star_M**4
                    Clstr_L = Clstr_L + Star_L  
                    

            else:                            

                    Star_L = Star_M**(a)
                    Clstr_L = Clstr_L + Star_L
                    


            try:
                logClstr_L = m.log10(Clstr_L)
            except ValueError:
                pass
                
            logaClmp_M = m.log10(aClmp_M[k])


            if (Star_M > 10):

                    counter += 1
                    tot_massive_L = tot_massive_L + Star_L
                

            if (Prev_L == 0):

                logPrev_L = np.NaN

            else:    

                logPrev_L = m.log10(Prev_L)
            
            Ex_M = Prev_M > aClmp_M[k]*0.99

            #Ex_L = Clstr_L > (1+Tol_L)*21.64*(aClmp_M[k]**1.1849)
            #this one is fwhm upper limit
            #Ex_L = Clstr_L > (1+Tol_L)*33.27*(aClmp_M[k]**1.223)
            
            #Ex_L = Clstr_L > (1+Tol_L)*54.57*(aClmp_M[k]**0.9915)
            #Ex_L = Clstr_L > (1+Tol_L)*6.681*(aClmp_M[k]**1.0511)
            #this one is fwhm mean
            Ex_L = Clstr_L > (1+Tol_L)*12.8529*(aClmp_M[k]**1.0492)

            
            if (Ex_M and Ex_L): 
                
                Clstr_M, Clstr_L, tot_massive_L, Star_L = 0.0, 0.0, 0.0, 0.0
                Prev_M, Prev_L = 0.0, 0.0                           
                Most_M, Most_L = 0.1, 0.0

                Clstr_pop = 0

            if (Ex_M):

                Clstr_pop = Clstr_pop - 1

                Clstr_M = Prev_M

                Clstr_L = Prev_L 
                
                count_Ex_M += 1

                break

            if (Ex_L): 

                Clstr_pop = Clstr_pop - 1

                Clstr_M = Prev_M

                Clstr_L = Prev_L 
                
                count_Ex_L += 1

                break
            
            #counter += 1


        #if debug:
            #print ('Cluster', j, Clstr_L, Clstr_M)


        dateTimeObj = datetime.now()
        f = open('Clusterproperties.txt', 'a')
        f.write('Cluster' + ', ' + str(j) + ', ' + str(Clstr_L) + ',' + str(Clstr_M) + ', ' + str(Clstr_pop) + ', ' + str(aClmp_M[k]) + ', ' + str(dateTimeObj) + ', ' + str(counter) + ', ' + str(tot_massive_L) + "\n")


        # Calculate the Star Formation Efficiency (SFE)

        SFE = (Clstr_M/aClmp_M[k])*100
            

        # Record cluster SFE

        sSFE.append(SFE)
   
       
              
    SFE_all.append(sSFE)  
    
    print(count_Ex_L, count_Ex_M)

###############################################################################

####################### END of loop around each CLUSTER #######################

###############################################################################




                                                    
#BOXPLOT SFE vs INITIAL CLUMP MASS

#for n in SFE_all:
    #data_all.append(SFE_all[n])

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

print(np.std(data1))
print(np.std(data2))
print(np.std(data3))
print(np.std(data4))
print(np.std(data5))
print(np.std(data6))
print(np.std(data7))
print(np.std(data8))
print(np.std(data9))
print(np.std(data10))
print(np.std(data11))
print(np.std(data12))
print(np.std(data13))
print(np.std(data14))
print(np.std(data15))
print(np.std(data16))
print(np.std(data17))
print(np.std(data18))
print(np.std(data19))
print(np.std(data20))

print(np.mean(data1))
print(np.mean(data2))
print(np.mean(data3))
print(np.mean(data4))
print(np.mean(data5))
print(np.mean(data6))
print(np.mean(data7))
print(np.mean(data8))
print(np.mean(data9))
print(np.mean(data10))
print(np.mean(data11))
print(np.mean(data12))
print(np.mean(data13))
print(np.mean(data14))
print(np.mean(data15))
print(np.mean(data16))
print(np.mean(data17))
print(np.mean(data18))
print(np.mean(data19))
print(np.mean(data20))


print(np.std(data1)/np.mean(data1))
print(np.std(data2)/np.mean(data2))
print(np.std(data3)/np.mean(data3))
print(np.std(data4)/np.mean(data4))
print(np.std(data5)/np.mean(data5))
print(np.std(data6)/np.mean(data6))
print(np.std(data7)/np.mean(data7))
print(np.std(data8)/np.mean(data8))
print(np.std(data9)/np.mean(data9))
print(np.std(data10)/np.mean(data10))
print(np.std(data11)/np.mean(data11))
print(np.std(data12)/np.mean(data12))
print(np.std(data13)/np.mean(data13))
print(np.std(data14)/np.mean(data14))
print(np.std(data15)/np.mean(data15))
print(np.std(data16)/np.mean(data16))
print(np.std(data17)/np.mean(data17))
print(np.std(data18)/np.mean(data18))
print(np.std(data19)/np.mean(data19))
print(np.std(data20)/np.mean(data20))




#count = 1
#for n in range(2, 4, 0.1):
  #  x[count] = [n] * 1000
   # count+=1


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



xs = x_all
ys = data_all
myfit = sp.polyfit(xs,ys,2)
print(myfit[0],"x^2 + ", myfit[1], 'x + ', myfit[2])
fitted_ys =  myfit[0] * xs**2 + myfit[1] * xs + myfit[2]
#plt.boxplot(SFE_all, positions = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21], showfliers = False, showmeans = True)
#plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21], [2, 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4])
#plt.xlabel('Log_10(Mass)')
#plt.ylabel('SFE (%)')
#plt.show() 
#plt.savefig('boxplot.pdf', format = 'pdf', dpi = 1200)



#plt.plot(xs,ys,'k.')
plt.plot(xs,fitted_ys,'r-')
plt.xlabel('Log_10(Mass)')
plt.xlim
plt.ylabel('SFE (%)')
#plt.show()
plt.savefig('boxplot_line.pdf', format = 'pdf', dpi = 1200)




