#=== Importing libraries

import numpy as np

import datetime

import time
import timeit

from math import cos, sin, atan, pi, sqrt

####################################################################################################

number_of_periods = 8000
number_of_divisions_per_period = 8000
w = 5.61

####################################################################################################

inicio = timeit.default_timer()
####################################################################################################

data = np.genfromtxt("state_space.txt", skip_header = 8, unpack = True)

#==#
phi_poincare = data[0]
    
dotphi_poincare = data[1]

time = data[2]
    
#number_of_values = len(data[0])
    
f = open("poincare_orb.txt",'w')

i= datetime.datetime.now()
f.write("The simulation was run in" + str(i) + "\n" + "\n")

f.write("number_ of_periods = " + str(number_of_periods) + "\n")
f.write("number_of_divisions_per_period = " + str(number_of_divisions_per_period) + "\n")
f.write("w = " + str(w) + "\n\n")  
 
f.write( "     " + "        phi          " +"                 dotphi         " + "                 time          "+ "\n")

counter = 0   

print len(data[0])/number_of_divisions_per_period

for i in range(number_of_periods):
    for j in range(number_of_divisions_per_period):
        if j == 0:
           value = ("   %23.15E        %23.15E        %23.15E \n" %(data[0][counter], data[1][counter], data[2][counter]))
	   f.write(value)        
        counter += 1        




f.close()



####################################################################################################
fim = timeit.default_timer()
 
print ('duracao: %f' % (fim - inicio))
