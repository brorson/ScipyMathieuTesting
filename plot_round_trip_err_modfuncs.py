#!/usr/bin/env python

import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize



def modmc1(m, q, x):
    return scipy.special.mathieu_modcem1(m, q, x)
def modmc2(m, q, x):
    return scipy.special.mathieu_modcem2(m, q, x)

def modms1(m, q, x):
    return scipy.special.mathieu_modsem1(m, q, x)
def modms2(m, q, x):
    return scipy.special.mathieu_modsem2(m, q, x)

def mathieu_a(m, q):
    return scipy.special.mathieu_a(m, q)
def mathieu_b(m, q):
    return scipy.special.mathieu_b(m, q)

#==========================================================
# Round trip modmc1
NN = 100
v = np.linspace(2,5,NN)
h = 5e-6
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)

error = []
idx = 1
for q in qs:
    for m in range(MM):  
        fm3 = modmc1(m, q, v - 3 * h)[idx]
        fm2 = modmc1(m, q, v - 2 * h)[idx]
        fm1 = modmc1(m, q, v - h)[idx]
        fp1 = modmc1(m, q, v + h)[idx]
        fp2 = modmc1(m, q, v + 2 * h)[idx]
        fp3 = modmc1(m, q, v + 3 * h)[idx]

        ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60

        [y,yd] = modmc1(m,q,v)
        a = mathieu_a(m,q)

        r = ydd/(h) - (a - 2 * q * np.cosh(2*v)) * y
        stddev = np.std(r)
        l2norm = np.linalg.norm(y, ord=2) / np.sqrt(len(y))

        error.append( stddev/l2norm)
        if not np.isclose( l2norm , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {l2norm}')

            
M, Q = np.meshgrid(range(MM), qs)

error_array = np.array(error).reshape(M.shape)
levels = np.arange(-25, 36, 5)

plt.figure(1)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for modmc1")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render



#==========================================================
# Round trip modmc2
NN = 100
v = np.linspace(2,5,NN)
h = 5e-6
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)

error = []

for q in qs:
    for m in range(MM):    
        fm3 = modmc2(m, q, v - 3 * h)[1]
        fm2 = modmc2(m, q, v - 2 * h)[1]
        fm1 = modmc2(m, q, v - h)[1]
        fp1 = modmc2(m, q, v + h)[1]
        fp2 = modmc2(m, q, v + 2 * h)[1]
        fp3 = modmc2(m, q, v + 3 * h)[1]

        ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60

        [y,yd] = modmc2(m,q,v)
        a = mathieu_a(m,q)

        r = ydd/(h) - (a - 2 * q * np.cosh(2*v)) * y
        stddev = np.std(r)
        l2norm = np.linalg.norm(y, ord=2) / np.sqrt(len(y))

        error.append(stddev/l2norm)
        if not np.isclose( stddev/np.linalg.norm(y, ord=2), 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {l2norm}')

M, Q = np.meshgrid(range(MM), qs)

error_array = np.array(error).reshape(M.shape)
levels = np.arange(-25, 36, 5)

plt.figure(2)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for modmc2")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
import math
# Round trip modms1
NN = 100
v = np.linspace(2,5,NN)
h = 5e-6
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)

error = []

for q in qs:
    for m in range(MM):  
        fm3 = modms1(m, q, v - 3 * h)[1]
        fm2 = modms1(m, q, v - 2 * h)[1]
        fm1 = modms1(m, q, v - h)[1]
        fp1 = modms1(m, q, v + h)[1]
        fp2 = modms1(m, q, v + 2 * h)[1]
        fp3 = modms1(m, q, v + 3 * h)[1]

        ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60

        [y,yd] = modms1(m,q,v)
        a = mathieu_b(m,q)

        r = ydd/(h) - (a - 2 * q * np.cosh(2*v)) * y
        stddev = np.std(r)
        l2norm = np.linalg.norm(y, ord=2) 
        
        error.append(stddev / l2norm)
        if not np.isclose( stddev/ l2norm , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = {l2norm}')


M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-25, 36, 5)

plt.figure(3)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()

plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round rrip rel error for modms1")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


#==========================================================
# Round trip modms2
NN = 100
v = np.linspace(2,5,NN)
h = 5e-6
tol = 1e-3
MM = 50
qs = np.logspace(-4,4,10)
error = []

for q in qs:
    for m in range(MM):  
        fm3 = modms2(m, q, v - 3 * h)[1]
        fm2 = modms2(m, q, v - 2 * h)[1]
        fm1 = modms2(m, q, v - h)[1]
        fp1 = modms2(m, q, v + h)[1]
        fp2 = modms2(m, q, v + 2 * h)[1]
        fp3 = modms2(m, q, v + 3 * h)[1]

        ydd = -fm3/60 + 3*fm2/20 - 3*fm1/4 + 3*fp1/4 - 3*fp2/20 + fp3/60

        [y,yd] = modms2(m,q,v)
        a = mathieu_b(m,q)

        r = ydd/(h) - (a - 2 * q * np.cosh(2*v)) * y
        stddev = np.std(r)
        l2norm = np.linalg.norm(r, ord=2) 
        error.append( stddev / l2norm)

        if not np.isclose(stddev / l2norm , 0, atol=tol):
            print(f'm = {m}, q = {q}, error = { l2norm}')


M, Q = np.meshgrid(range(MM), qs)
error_array = np.array(error).reshape(M.shape)
levels = np.arange(-25, 36, 5)

plt.figure(4)
plt.contourf(M, np.log10(Q), np.log10(error_array), levels=levels, cmap='viridis')
plt.colorbar()
plt.xlabel("m")
plt.ylabel("log10(q)")
plt.title("Round trip rel error for modms2")
plt.draw()           # Draw the plt.figure
plt.pause(0.001)     # Give GUI time to render


# Do this to stop the program from closing all windows when
# exiting
input("Press Enter to close...") 
