#!/usr/bin/python3
import os
import signal
import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np

pmin = 0.16
pmax = 0.19

plt.ion()
fig, ax = plt.subplots(layout='constrained')
terminating = False


def quit(*args):
    global terminating
    terminating = True
    fig.canvas.stop_event_loop()


def update_errorbar(errobj, x, y, yerr):
    ln, caps, bars = errobj
    barsy, = bars
    ln.set_data(x, y)
    barsy.set_segments([
        np.array([[x, yt], [x, yb]])
        for x, yt, yb in zip(x, y + yerr, y - yerr)
    ])


def read_data():
    sp.run(['scons', '-s'], check=True)
    data = {}
    for root, dirs, files in os.walk('result'):
        for fname in sorted(files):
            with open(f'{root}/{fname}') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    v = [float(word) for word in line.split()]
                    L = v[0]
                    p = v[1]
                    if p < pmin or p > pmax:
                        continue
                    if L not in data:
                        data[L] = []
                    data[L].append(v)

    return {key: np.array(value) for key, value in data.items()}


signal.signal(signal.SIGTERM, quit)
signal.signal(signal.SIGINT, quit)
fig.canvas.mpl_connect('close_event', quit)

graphs = {}
while not terminating:
    data = read_data()
    for key, v in data.items():
        x = v[:, 1]
        y = v[:, 3] * key**(-1.75)
        yerr = v[:, 4] * key**(-1.75)

        if key not in graphs:
            graphs[key] = ax.errorbar(x, y, yerr, fmt='+', label=str(int(key)))
        else:
            update_errorbar(graphs[key], x, y, yerr)

    ax.legend()
    fig.canvas.draw()
    fig.canvas.start_event_loop(2)
