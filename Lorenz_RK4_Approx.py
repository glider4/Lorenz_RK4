# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 :12:03 2018
@author: mathemacode

RK4 system of 3 equations, all at once
Testing with Lorenz equations, 2D and 3D plots included

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

a = 0   # lower bound
b = 20 # upper bound
n = 20000  # number of steps
h = (b-a)/n  # step size
p = 0  # flag for while loop

# intialize lists for plot
t = np.zeros(int(n+1))
x1 = np.zeros(int(n+1))
x2 = np.zeros(int(n+1))
x3 = np.zeros(int(n+1))

# define initial values
t[0] = 0
x1[0] = 15
x2[0] = 15
x3[0] = 36

# initivalize slots for values of RK4
i = np.zeros(4)  # for x1
j = np.zeros(4)  # for x2
k = np.zeros(4)  # for x3

# derivatives (p for prime)
def x1p(t, x1, x2, x3):
    return 10 * (x2 - x1)

def x2p(t, x1, x2, x3):
    return x1 * (28 - x3) - x2

def x3p(t, x1, x2, x3):
    return (x1 * x2) - (8/3)*x3

# iterate n times
while p < n:
    
    ''' RK4 Calcuations, adapted into Numpy arrays '''

    i[0] = (h) * x1p(t[p], x1[p], x2[p], x3[p])
    j[0] = (h) * x2p(t[p], x1[p], x2[p], x3[p])
    k[0] = (h) * x3p(t[p], x1[p], x2[p], x3[p])
    
    i[1] = (h) * x1p(t[p] + (h/2), x1[p] + (1/2)*i[0], x2[p] + (1/2)*j[0], x3[p] + (1/2)*k[0])
    j[1] = (h) * x2p(t[p] + (h/2), x1[p] + (1/2)*i[0], x2[p] + (1/2)*j[0], x3[p] + (1/2)*k[0])
    k[1] = (h) * x3p(t[p] + (h/2), x1[p] + (1/2)*i[0], x2[p] + (1/2)*j[0], x3[p] + (1/2)*k[0])
    
    i[2] = (h) * x1p(t[p] + (h/2), x1[p] + (1/2)*i[1], x2[p] + (1/2)*j[1], x3[p] + (1/2)*k[1])
    j[2] = (h) * x2p(t[p] + (h/2), x1[p] + (1/2)*i[1], x2[p] + (1/2)*j[1], x3[p] + (1/2)*k[1])
    k[2] = (h) * x3p(t[p] + (h/2), x1[p] + (1/2)*i[1], x2[p] + (1/2)*j[1], x3[p] + (1/2)*k[1])
    
    i[3] = (h) * x1p(t[p] + (h), x1[p] + i[2], x2[p] + j[2], x3[p] + k[2])
    j[3] = (h) * x2p(t[p] + (h), x1[p] + i[2], x2[p] + j[2], x3[p] + k[2])
    k[3] = (h) * x3p(t[p] + (h), x1[p] + i[2], x2[p] + j[2], x3[p] + k[2])
    
    x1[p + 1] = x1[p] + (1/6) * ( i[0] + (2*i[1]) + (2*i[2]) + i[3] )
    x2[p + 1] = x2[p] + (1/6) * ( j[0] + (2*j[1]) + (2*j[2]) + j[3] )
    x3[p + 1] = x3[p] + (1/6) * ( k[0] + (2*k[1]) + (2*k[2]) + k[3] )
    
    t[p + 1] = t[p] + h
    
    p += 1  # advance flag variable


''' Plot Results of RK4 '''

# 2D plot to show consistency of estimations
plt.plot(t, x1, '-', label = "X1 RK4", linewidth = 2.0)
plt.plot(t, x2, '-', label = "X2 RK4", linewidth = 2.0)
plt.plot(t, x3, '-', label = "X3 RK4", linewidth = 2.0)

plt.xlim(a, b)
plt.title("RK4 Method Result: Lorenz Problem")
plt.xlabel("t")
plt.ylabel("Value")
plt.legend()
plt.show()

# 3D plot, because, Lorenz
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x1, x2, x3, linewidth = 1.0, color = "red")
plt.style.use('ggplot')
plt.show()
    

    
    
    
    
    
    
    
    
    



