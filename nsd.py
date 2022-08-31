#!/usr/bin/env python
# -*- coding: utf-8 -*-
# developed & tested with Python 3.9

"""Generates an estimation of the noise amplitude spectral density (NSD) of a time series"""

from scipy import signal, integrate
import numpy as np


def series_rms(values: tuple) -> float:
    """AC root mean square (ACRMS) of given values

    RMS without DC (mean removed)

    Parameters:
    -----------
    values: array_like

    Returns
    -------
    out : float
    """
    mean = values.mean()
    return ((values - mean)**2).mean()**0.5


def nsd_rms(nsd: tuple) -> float:
    """Root mean square (RMS) of given nsd

    Parameters:
    -----------
    nsd: [array_like, array_like]
        nsd[0]: frequencies, nsd[1]: nsd values

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
    length: int
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
    length: int
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


def get(
    ts_values: tuple,
    sample_frequency: float,
    nsd_bins: int,
    window_function=window_HFT90D,
    crop=np.s_[3:-1],
):
    """Estimation of the noise amplitude spectral density (NSD) of a time series.

    E.g. it is used in electronic component datasheets to show 1/f (1/f^n, n:real) corner and white noise part (flat) of voltage noise, but same can be applied to current noise.
    From the white noise part (flat) the rms value of a desired bandwidth can be calculated: rms = value * bandwidth^0.5

    Warning:
        If the acquisition time (aperture) of the input signal is lower than the time between samples of the output (ts_values), then the nsd is not correct!
        Nsd using acquisition time (aperture) as sample_frequency only shows correct white noise part.
        Nsd using time between samples as sample_frequency only shows correct 1/f (1/f^n) part, the white noise part is too high, being (Tsps/Tacq)^0.5, e.g. for AZ: (40ms/20ms)^0.5, so ~1.4
        Neither will give the correct 1/f (1/f^n, n:real) noise corner!

    Parameters:
    -----------

    ts_values : array_like
        Time series of the signal (amplitude)

    sample_frequency : int/float
        The sample frequency in Hz (SPS - Samples per second)

    nsd_bins : int
        number of NSD bins (points) to be calculated
        needs to be lower than the ts_values count


    window_function : function(length: int), optional
        returns: [int, array_like]
            first index is the optimum overlapping in samples for the window
            second index are the values of the window function
        default: window_HFT90D

    crop : slice, optional
        default: the last and the first 3 nsd values are dropped as they are not reliable

    Returns
    -------
    out : [array, array]
        The frequencies (index 0) and the corresponding NSD values (index 1)

    """
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
    # crop & transform psd to nsd
    return (frequencies[crop], psd[crop]**0.5)
