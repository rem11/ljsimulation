#!/usr/bin/env python3

from random import random
from sys import argv

a = 20 # particles on rectangle side
s = 1 # distance between particles
N = a ** 2
T = 1000
dt_out = 3
dt = 0.0001
b = 1
r_m = 1
r_z = 3

d = 0.01 # coordinate dispersion

f = open(argv[1], 'w')
f.write(f'{N} {T} {dt_out} {dt} {b} {r_m} {r_z}\n')

for i in range(0, a):
    for j in range(0, a):
        x = - (a - 1) * s / 2 + s * i + d * (random() - 0.5)
        y = - (a - 1) * s / 2 + s * j + d * (random() - 0.5)
        f.write(f'{x} {y} 0 0\n')

f.close()
