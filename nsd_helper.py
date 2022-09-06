#!/usr/bin/env python
# -*- coding: utf-8 -*-
# developed & tested with Python 3.9

"""Little helper for NSD, e.g. noise & tone generators, nsd plot, engineering formatter"""

import numpy as np
import colorednoise as cn
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter

def nsd_noise(samples, noise_types=['white', 'brownian'], sps=50, corner=[0.1, 1e-9], seed=4):
    '''generate colored noise
    supported noise types are: white, white-pink, pink, brownian (1/f^(0; 0.5; 1; 2))
    sps & corner is not supported atm - for 2^20 samples, 50 sps the corner is 0.1Hz at 1nV
    '''
    sum_noise = [0] * samples
    
    # backwards compatibility, breaks with v2.1
    np.random.seed(seed) # same result for each run
    def noise_function(exponent, samples):
        try:
            return cn.powerlaw_psd_gaussian(exponent, samples, random_state=seed)
        except:
            return cn.powerlaw_psd_gaussian(exponent, samples)


    for noise_type in noise_types:
        match = False
        if noise_type == 'white':
            sum_noise += noise_function(0, samples)/199.449e6
            match = True
        if noise_type == 'white-pink':
            sum_noise += noise_function(0.5, samples)/0.55e9
            match = True
        if noise_type == 'pink':
            sum_noise += noise_function(1, samples)/0.85e9
            match = True
        if noise_type == 'brownian':
            sum_noise += noise_function(2, samples)/4e7
            match = True
        if match == False:
            raise ValueError(f'unsupported noise_type: {noise_type}')

    return sum_noise

def tone(rms, freq, sps, samples):
    '''generate a tone
    '''
    ur = []
    for i in range (0, samples):
        t = i / sps
        ur.append(rms * 2**0.5 * np.cos(2 * np.pi * freq * t))
    return np.array(ur)

def plot(values_list, title, filename=None):
    '''plot NSDs
    optional: provide a filename to save plot as .png
    '''
    fig, ax = plt.subplots(figsize=(15, 10))

    for values in iterable(values_list):
        plt.plot(values['nsd'][0], values['nsd'][1], '-', ms=1, lw=1, alpha=0.8, label=values['label'])

    plt.title(title, fontsize=13)
    plt.xlabel('frequency in Hz')
    plt.ylabel(r'noise in V/$\sqrt{Hz}$')

    plt.loglog()
    plt.grid(True, which="both")
    
    formatter = EngFormatter(sep="")
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().yaxis.set_minor_formatter(formatter)
    plt.legend()
    fig.tight_layout()                                             # Adjust spacings w.r.t. figsize

    #save plot
    plt.rcParams['savefig.facecolor']='white'                      # set background color for saving, standard is transparent
    if filename is not None:
        plt.savefig(f'{filename}.png')  # change name accordingly

    plt.show()

def eng(values):
    strings = []
    try:
        for value in iter(values):
            strings.append(EngFormatter(sep="", places=1).format_eng(value))
        return strings
    except TypeError:
        return EngFormatter(sep="", places=1).format_eng(values)

def iterable(object):
    '''make object iterable if not already
    '''
    try:
        for n in iter(object):
            return object
    except TypeError:
        return tuple(object)
    

