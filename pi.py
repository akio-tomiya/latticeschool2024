#!/usr/bin/env python3

# sample program of Monte Carlo integration 
#   Issaku Kanamori, 2024 Sep. 11
#   prepared for Lattice Theory Summer School 2024
#   https://akio-tomiya.github.io/latticeschool2024/
#      Licence: CC0

# estimates
#   4 \int_0^1 dx \sqrt(1-x^2)   (= \pi)
# by using uniform random number in [0,1]

import random
import math


# number of sampling
#   If N is too large, numerical rounding error causes incorrect sum.
#   It is a good exercise to rewrite the summation to overcome this issue
#   In the below, partial summation is used. 
N=10000000

N_show=100000

# range of integration
x_min=0
x_max=1

total=0
partial=0
psize=10000
for i in range(1,N+1):
    x=random.uniform(x_min, x_max)
    y=4*math.sqrt(1-x*x)
    partial+=y
    if(i % psize==0):
        total+=partial
        partial=0
    if(i % N_show == 0):
        average=(total+partial)/i
        print(i, average)
