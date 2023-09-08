#!/usr/bin/python3
import os
import time
import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np

pmin = 0.16
pmax = 0.19

plt.ion()
fig, ax = plt.subplots(layout='constrained')


def read_data():
    sp.run(['scons'], check=True)
    data = {}
    for root, dirs, files in os.walk('result'):
        for fname in sorted(files):
            with open(f'{root}/{fname}') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    v = [float(word) for word in line.split()]
                    L = v[1]
                    p = v[2]
                    if p < pmin or p > pmax:
                        continue
                    if L not in data:
                        data[L] = []
                    data[L].append(v)

    return {key: np.array(value) for key, value in data.items()}

p = None
while True:
    data = read_data()
    ax.clear()
    for key, v in data.items():
        ax.errorbar(v[:, 2], v[:, 8]*key**(-1.75), yerr=v[:, 9]*key**(-1.75),
                    fmt='+', label=str(int(key)))
    ax.legend()
    plt.show()
    plt.pause(2)
