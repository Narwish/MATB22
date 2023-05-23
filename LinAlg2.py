#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:42:24 2022

@author: Lukas Early, Björn Follin - Group 9
"""

from  numpy import *
from matplotlib.pyplot import *
from scipy import optimize

from sympy.solvers import solve
from sympy import Symbol
from mpl_toolkits import mplot3d


A = array([[1, 1, 2], [1, 2, 1], [2, 1, 1], [2, 2, 1]])
B = array([1, -1, 1, -1])

ATA = matmul(A.T, A)
ATB = matmul(A.T, B)

L = linalg.solve(ATA, ATB)

M = lambda L: linalg.norm(matmul(A,L)-B)

minimum = optimize.fmin(M, array([0, 0, 0]))


ATA_inverse = linalg.inv(ATA)
def B_veriable(a):
    B = array([1, a, 1, a])
    ATB = matmul(A.T, B)
    L = matmul(ATA, ATB)
    AL = matmul(A, L)
    
    return linalg.norm(AL - B)

def plotter (lower, upper, step):
    values = linspace(lower, upper, step)
    residuals = [B_veriable(a) for a in values]
    return [residuals,values]
    
figure()
plot(plotter(-10,10,100)[1],plotter(-10,10,100)[0])


#Task 2
def Mat(a0,b0,c0,stop):
    i = 1
    ai = a0+3*b0+2*c0
    bi = -3*a0 +4*b0 +3*c0
    ci =2*a0 +3*b0 +c0
    while i < stop:
        i += 1
        a = ai+3*bi+2*ci
        b = -3*ai +4*bi +3*ci
        c =2*ai +3*bi +ci
        ai = a
        bi = b
        ci = c
    return [ai,bi,ci]       
       
        
z = Mat(8,3,12,200)
z_bottom = (z[0]**2)+(z[1]**2)+(z[2]**2)
z_b = z_bottom**0.5
V = array([z[0]/z_b,z[1]/z_b,z[2]/z_b])
A = array([[1, 3, 2], [-3, 4, 3], [2, 3, 1]])
VTAV = matmul(V.T, matmul(A,V))

def epi_V(lim):
    ep=10
    i=0   
    while ep > k:
        i += 1
        znew = Mat(8,3,12,i)
        znew_bottom = (znew[0]**2)+(znew[1]**2)+(znew[2]**2)
        znew_b = znew_bottom**0.5
        Vnew = array([znew[0]/znew_b,znew[1]/znew_b,znew[2]/znew_b])
        z_final = [Vnew[0]-V[0],Vnew[1]-V[1],Vnew[2]-V[2]]
        z_f = (z_final[0]**2)+(z_final[1]**2)+(z_final[2]**2)
        ep = z_f**.5
    return (ep,i)

def epi_VTAV(lim):
    ep=10
    i=0   
    while ep > k:
        i += 1
        znew = Mat(8,3,12,i)
        znew_bottom = (znew[0]**2)+(znew[1]**2)+(znew[2]**2)
        znew_b = znew_bottom**0.5
        Vnew = array([znew[0]/znew_b,znew[1]/znew_b,znew[2]/znew_b])
        VTAVnew = matmul(Vnew.T, matmul(A,Vnew))
        ep = abs(VTAV-VTAVnew)
    return (ep,i)
runthrough = linspace(.1, e**-14, 1000)
figure()
for k in runthrough:
      plot(epi_V(k))
      plot(epi_VTAV(k))

#Task 3
x_list = list()
y_list = list()
for i in range (-10,11):
    x_list.append(i/10)
    y_list.append(i/10)

def function(x,y):
    z = Symbol ("z")
    a = 2*x**2
    b = y**2
    c = 10*x*y
    d = 4*x
    e = 10*y
    root = solve(((a)-(b)+(2*(z**2))-(c)-(d*z)+(e*z)-1),z)
    return [x, y, root]

Total = []   

for n in range (len(x_list)):
    Total.append(function(x_list[n], y_list[n]))

x_list = array(x_list)
y_list = array(y_list)
Total = array(Total)

figure()
ax = axes(projection='3d')
ax.plot_trisurf(x_list, y_list, Total)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
    
    
    



