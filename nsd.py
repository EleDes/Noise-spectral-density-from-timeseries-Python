#!/usr/bin/env python
# -*- coding: utf-8 -*-
# developed & tested with Python 3.9

"""Generates an estimation of the noise amplitude spectral density (NSD) of a time series"""

import numpy as np
from scipy import signal, integrate
from scipy.optimize import curve_fit


def series_rms(values: tuple) -> float:
    """AC root mean square (ACRMS) of given values

    RMS without DC (mean removed)

    Parameters:
    -----------
    values : array_like

    Returns
    -------
    out : float
    """
    mean = values.mean()
    return ((values - mean)**2).mean()**0.5


def nsd_rms(nsd: tuple) -> float:
    """Root mean square (RMS) of given NSD

    Parameters:
    -----------
    nsd : [array_like, array_like]
        nsd[0]: frequencies, nsd[1]: NSD values

    Returns
    -------
    out : float
    """
    return integrate.trapezoid(y=nsd[1], x=nsd[0]**0.5)


def window_flattop(length: int):
    """Window flattop

    coeffs = [0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368]

    Parameters:
    -----------
    length : int
        number of samples to generate

    Returns
    -------
    out : [int, array]
        number of overlapping samples of the window (index 0) and array of the window values (index 1)
    """
    noverlap = int(length * 0.76)  # 0.76 is a guess
    return noverlap, signal.get_window("flattop", length, fftbins=False)


def window_HFT90D(length: int):
    """Window HFT90D - from https://holometer.fnal.gov/GH_FFT.pdf

    wj = 1 − 1.942604 cos(z) + 1.340318 cos(2z) − 0.440811 cos(3z) + 0.043097 cos(4z)

    Parameters:
    -----------
    length : int
        number of samples to generate

    Returns
    -------
    out : [ínt, array]
        number of overlapping samples of the window (index 0) and array of the window values (index 1)
    """
    noverlap = int(length * 0.76)
    window = np.empty(length)
    for i in range(length):
        z = (2.0 * np.pi * i) / length
        window[i] = (
            1
            - (1.942604 * np.cos(z))
            + (1.340318 * np.cos(2 * z))
            - (0.440811 * np.cos(3 * z))
            + (0.043097 * np.cos(4 * z))
        )

    return noverlap, window

# TODO overlapping
def smooth(ordered_nsd, nsd_bins=64, filter_function=np.mean):
    """Smooth the NSD evenly spaced in log space

    Parameters:
    -----------

    ordered_nsd : [array_like, array_like]
        ordered_nsd[0]: ascending ordered frequencies, ordered_nsd[1]: corresponding NSD values

    nsd_bins : int, optional
        number of NSD bins (points) to be calculated
        default: 64

    filter_function : function(array_like), optional
        returns: [array_like]
            filtered values
        default: np.mean

    Returns
    -------
    out : [array, array]
        The frequencies (index 0) and the corresponding smoothed NSD values (index 1)
    """
    geomspace = np.geomspace(ordered_nsd[0][0],ordered_nsd[0][-1],num=nsd_bins)
    frequencies = []
    nsd = []
    i1 = int(0)
    i2 = int(1)
    for freq in geomspace[1:]:
        while ordered_nsd[0][i2] < freq:
            i2 +=1
        freq_mean = ordered_nsd[0][i1:i2+1].mean()
        values_filtered = filter_function(ordered_nsd[1][i1:i2+1])
        frequencies.append(freq_mean)
        nsd.append(values_filtered)
        i1 = i2
    return (frequencies, nsd)

def fit_function(freq, slope, freq_exp, white):
    '''NSD fit function for 1/f^n + white noise

    f = (slope / freq**(freq_exp/2)) + white
    '''
    return ((slope / freq**(freq_exp/2)) + white)

def fit_function_loglog(freq, slope, freq_exp, white):
    '''NSD fit function for 1/f^n + white noise in the loglog space

    f = np.log10((10**slope / ((10**freq) ** ((10**freq_exp) / 2))) + 10**white)
    '''
    return np.log10((10**slope / ((10**freq) ** ((10**freq_exp) / 2))) + 10**white)

def fit(nsd, fit_function=fit_function):
    """Least squares curve fitting the NSD

       Gives worse fits than fitting in loglog space, see fit_loglog

    Parameters:
    -----------

    nsd : [array, array]
        NSD from get()

    fit_function : function, optional
        The function in loglog space to fit the nsd
        default: fit_function_loglog = np.log10((10**slope / ((10**freq) ** ((10**freq_exp) / 2)) + 10**white))

    Returns
    -------
    out : [array, array]
        The optimal values for the parameters (index 0) and the corresponding estimated covariance (index 1)
    """
    return curve_fit(fit_function, nsd[0], nsd[1])

def fit_loglog(nsd, fit_function=fit_function_loglog):
    """Least squares curve fitting the NSD in loglog space

       Gives better fits than fitting in linear space, see fit

    Parameters:
    -----------

    nsd : [array, array]
        NSD from get()

    fit_function : function, optional
        The function in loglog space to fit the nsd
        default: fit_function_loglog = np.log10((10**slope / ((10**freq) ** ((10**freq_exp) / 2)) + 10**white))

    Returns
    -------
    out : [array, array]
        The optimal values for the parameters (index 0) and the corresponding estimated covariance (index 1)
    """
    popt, pcov = curve_fit(fit_function, np.log10(nsd[0]), np.log10(nsd[1]))
    return (10**popt, pcov)

def get(
    ts_values: tuple,
    sample_frequency: float,
    nsd_bins: int=None,
    window_function=window_HFT90D,
    crop=np.s_[3:-1],
):
    """Estimation of the noise amplitude spectral density (NSD) of a time series.

    E.g. it is used in electronic component datasheets to show 1/f (1/f^n, n:real) corner and white noise part (flat) of voltage noise, but same can be applied to current noise.
    From the white noise part (flat) the rms value of a desired bandwidth can be calculated: rms = value * bandwidth^0.5

    Warning:
        If the acquisition time of the input signal (aperture) is lower than the time between samples of the output (ts_values), than the NSD is incorrect (e.g. AZ on ADCs/DMMs)!
        NSD using acquisition time (aperture) as sample_frequency only shows correct white noise part.
        NSD using time between samples as sample_frequency only shows correct 1/f (1/f^n, n:real) part, the white noise part is too high, being (Tsps/Tacq)^0.5, e.g. for AZ: (40ms/20ms)^0.5, so ~1.4
        Neither will give the correct 1/f (1/f^n, n:real) noise corner!

    Parameters:
    -----------

    ts_values : array_like
        Time series of the signal (amplitude)

    sample_frequency : int/float
        The sample frequency in Hz (SPS - Samples per second)

    nsd_bins : int, optional
        number of NSD bins (points) to be calculated
        needs to be lower than the ts_values count
        default: 1/4 count of ts_values

    window_function : function(length: int), optional
        returns: [int, array_like]
            first index is the optimum overlapping in samples for the window
            second index are the values of the window function
        default: window_HFT90D

    crop : slice, optional
        default: the last and the first 3 NSD values are dropped as they are not reliable

    Returns
    -------
    out : [array, array]
        The frequencies (index 0) and the corresponding NSD values (index 1)
    """
    if nsd_bins is None:
        nsd_bins = int(len(ts_values)/4)
    noverlap, window = window_function(nsd_bins)
    frequencies, psd = signal.welch(
        x=ts_values,
        fs=sample_frequency,
        window=window,
        # nperseg = len(window),
        noverlap=noverlap,
        # nfft = nperseg,
        # detrend = False,
        # return_onesided = True,
        # scaling = 'density',
        # axis = -1,
        # average = 'mean'
    )
    # crop & transform PSD to NSD
    return (frequencies[crop], psd[crop]**0.5)
