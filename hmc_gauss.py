#!/usr/bin/env python3

# sample program of HMC 
#   Issaku Kanamori, 2024 Sep. 11
#   prepared for Lattice Theory Summer School 2024
#   https://akio-tomiya.github.io/latticeschool2024/
#      Licence: CC0

# 
#   distribtuion: exp(-x^2/2)
#   observable: x^2

import random
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# time to run
Ntraj=1000

# thermalization
Ntherm=100

# delta t = 1/Nstep  (Delta t =1)
Nstep=2

# initial configuration
x0=0

def Hamiltonian(p,x):
    return 0.5*p*p+0.5*x*x

def force(x):
    return -x

def MD_leapfrog(p, x, Nstep):
    Delta=1.0
    dt=Delta/Nstep

    # initial half step
    p+=force(x)*dt/2

    # leap frog
    for i in range(Nstep-1):
        x+=p*dt
        p+=force(x)*dt

    # last steps
    x+=p*dt
    p+=force(x)*dt/2

    return p,x

def calc_ave(vals, binsize=1):
    # make binned date
    binned=[]
    i=0
    tmp=0
    for val in vals:
        tmp+=val
        i+=1
        if(i==binsize):
            binned.append(tmp/binsize)
            i=0
            tmp=0

    num=len(binned)
    ave=sum(binned)/num
    err2=sum( (val-ave)**2 for val in binned)/num/(num-1)
    err=math.sqrt(err2)
    return num, ave, err

def show_histgram(conf):
    fig = plt.figure()
    plt.hist(conf, bins=20, density=True)
    points = np.linspace(-5,5,100)
    gauss=norm.pdf(points,0,1)
    plt.plot(points, gauss, color="r")
    plt.savefig("hmc_gauss.pdf")
    #plt.show()

print("# i x x^2 accpted dH exp(-dH)")

x=x0

obs=[]
configurations=[]

for i in range(Ntraj):
    x_try=x

    p=random.normalvariate(0,1.0) # initial momentum
    H_ini=Hamiltonian(p,x_try)    # initial Hamiltonian

    p, x_try = MD_leapfrog(p, x_try, Nstep)    # MD step

    H_fin=Hamiltonian(p,x_try)    # fintal Hamiltonian

    dH = H_fin -H_ini
    r=random.uniform(0,1)  # r in [0,1]
    accepted=0
    if(math.exp(-dH) > r):
        # new x is the trial one
        x=x_try
        accepted=1
    else:
        # new x is still x
        pass

    print(i, x, x*x, accepted, dH, math.exp(-dH))
    if(i>= Ntherm):
        obs.append(x*x)
        configurations.append(x)

n,ave,err=calc_ave(obs)
#n,ave,err=calc_ave(obs,10)  # average and error with binning

show_histgram(configurations)
print("# done")
print("# <x^2>: num, ave, err = ",n, ave, err)
