#!/usr/bin/env python3

from os import mkdir, path
from sys import argv
import numpy as np
from imageio import imwrite
from shutil import rmtree
from subprocess import run

f = open(argv[1], 'r');

[N, M] = f.readline().split()

dirname = 'images'
rmtree(dirname, ignore_errors = True)
mkdir(dirname)

scale = 10
size = round(20 * scale * 1.5)

for i in range(0, int(M)):
    t = float(f.readline())
    frame = np.zeros([size, size], np.uint8)
    print(f'Writing frame {i} of {M}')
    for j in range(0, int(N)):
        [x, y, v_x, v_y] = f.readline().split()
        p_x = round(float(x) * scale + size / 2)
        p_y = round(float(y) * scale + size / 2)
        if (p_x >= 0 and p_x < size and p_y >= 0 and p_y < size):
            frame[p_x, p_y] = 255
        imwrite(path.join(dirname, f'{i}.png'), frame);

f.close()

print('Starting ffmpeg...')
run(['ffmpeg', '-framerate', '30', '-pattern_type', 'sequence', '-i', path.join(dirname, '%d.png'), '-c:v', 'libx264', '-pix_fmt', 'yuv420p', argv[2]])
