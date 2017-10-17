# =========================================================================== #
# ========================= Importing libraries ============================= #
# =========================================================================== #

import numpy as np

import datetime

import timeit

from math import cos, sin, atan, pi, sqrt

inicio = timeit.default_timer()
# =========================================================================== #
# ======================= Defining some constant ============================ #
# =========================================================================== #

# ------------------------- Geometric properties ---------------------------- #

a = 16.0
b = 6.0
D = 9.5
d = 4.8

ddotphi = 24.6

# -------------------------- Inertia properties --------------------------- #

m = 14.7

I = 1738.0

g = 981.0

# ------------------------- Measured properties --------------------------- #

zeta = 236.8
k = 2470.0

mu = 1272.0

# -------------------------- Another properties --------------------------- #

q = 10.**6.

w = 5.26738222216

w0 = sqrt((m*g*D)/I)

# ========================================================================= #
# ============ Defining the Fourth order Runge-Kutta function ============= #
# ========================================================================= #


def RK4(function, state_vector, dt, t):

    number_of_dimensions = len(state_vector)
    k1 = np.zeros(number_of_dimensions)
    k2 = k1
    k3 = k1
    k4 = k1
    u1 = np.zeros(number_of_dimensions)
    u2 = u1
    u3 = u1
    # K1
    k1 = function(state_vector, t)
    # K2
    t1 = t + dt/2.0
    u1 = state_vector + k1*dt/2.0
    k2 = function(u1, t1)
    # K3
    t2 = t + dt/2.0
    u2 = state_vector + k2*dt/2.0
    k3 = function(u2, t2)
    # K4
    t3 = t + dt
    u3 = state_vector + k3*dt
    k4 = function(u3, t3)
    state_vector = state_vector + dt*(k1 + 2.0*(k2 + k3) + k4)/6.0
    return state_vector

# ========================================================================== #
# =================== Defining nonlinear equation function================== #
# ========================================================================== #


def Pendulum(state_vector, actual_time):

    t = actual_time
    delta_f = sqrt(a**2 + b**2 - 2*a*b*cos(w*t)) - (a - b)
    dotphi = state_vector[1]
    ddotphi = (((-(k*(d**2)/(2*I)))*state_vector[0]) -
               (zeta/I)*state_vector[1] +
               (((k*d)/(2*I))*delta_f) -
               ((m*g*D*sin(state_vector[0]))/(2*I)) -
               ((2*mu)/(pi*I))*atan(q*state_vector[1]))
    state_vector_out = np.array([dotphi, ddotphi])
    return state_vector_out

# ========================================================================== #
# ====================== Defining the initial status ======================= #
# ========================================================================== #

w_bifurc = []
phi_bifurc = []
phi = 0.0
dotphi = 0.0
number_of_periods = 5000
number_of_divisions_per_period = 500
step = (2*pi)/(w*number_of_divisions_per_period)
state_vector = np.array([phi, dotphi])
out = state_vector
phi = []
dotphi = []
t = 0.0
trash = (number_of_periods/2)*number_of_divisions_per_period
nome = "state_space" + str(w) + ".txt"
f = open(nome, "w")

i = datetime.datetime.now()
f.write("The simulation was run in" + str(i) + "\n" + "\n")

f.write("number_ of_periods = " + str(number_of_periods) + "\n")
f.write(("number_of_divisions_per_period = " +
        str(number_of_divisions_per_period) + "\n"))
f.write("w = " + str(w) + "\n\n")
f.write(("     " + "        phi          " +
        "                 dotphi         " +
         "                 time          " + "\n"))
for i in range(1, number_of_periods):
    for j in range(1, number_of_divisions_per_period):
        print(timeit.default_timer() - inicio)
        t = t + step
        out = RK4(Pendulum, state_vector, step, t)
        state_vector = out
        if i*j >= trash:
            value = ("   %23.15E        %23.15E        %23.15E \n"
                     % (state_vector[0], state_vector[1], t))
            f.write(value)
fim = timeit.default_timer()
print('duracao: %f' % (fim - inicio))
