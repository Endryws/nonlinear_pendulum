# -*- coding: utf-8 -*-
"""

Created on Fri Oct 27 09:59:53 2017

@author: endryws

This algoritim is based on Determining Lyapunov Exponents from a Time-Series
published by Wolf, Alan Swift, Jack B. Swinney, Harry L. Vastano, John A.
in Physica D (1985, volume 16, p 285 - 317)


"""
# ==================================================================== #
# =============== Determining the problem constants ================== #
# ==================================================================== #


A = 1.
B = 1.
C = 1.


# ==================================================================== #
# ============== Determinating the inicial conditions ================ #
# ==================================================================== #


x = [0,0,0]


# ==================================================================== #
# =============== Initializing the function vector =================== #
# ==================================================================== #


f = []


# ==================================================================== #
# ============= Initializing the Loretz's equation =================== #
# ==================================================================== #

f.append(A*(x[1] - x[0]))
f.append(B*x[0] - x[0]*x[2] - x[1])
f.append(x[0]*x[1] - C*x[2])

# ==================================================================== #
# ================== Linearizing the system ========================== #
# ==================================================================== #


for i in range(0,3):
    f.append()










