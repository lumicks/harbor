#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def simulate_naive_1D_brownian_motion(duration=10, nSteps=10e4):
    print('Running simulation')
    time = np.linspace(0, duration, nSteps)
    position = np.random.standard_normal(size = nSteps)
    plt.plot(time, position)
    plt.show()

