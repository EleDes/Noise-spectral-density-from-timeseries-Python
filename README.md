# Noise spectral density (NSD) from time series in Python

This tool calculates the [estimated noise amplitude spectral density](https://en.wikipedia.org/wiki/Noise_spectral_density) of a time series.

From the white noise part the rms value of a desired bandwidth can be calculated: rms = value * bandwidth^0.5

E.g. it is used in electronic component datasheets to show 1/f (1/f^n, n:real) corner and white noise part (flat) of voltage, but same can be applied to current noise.


# NSD function get()
Based on [Welch's method](https://en.wikipedia.org/wiki/Welch%27s_method).


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



## Warning:
If the acquisition time of the input signal (aperture) is lower than the time between samples of the output (ts_values), than the NSD is incorrect (e.g. AZ on ADCs/DMMs)!

Nsd using acquisition time (aperture) as sample_frequency only shows correct white noise part.

Nsd using time between samples as sample_frequency only shows correct 1/f (1/f^n, n:real) part, the white noise part is too high, being (Tsps/Tacq)^0.5, e.g. for AZ: (40ms/20ms)^0.5, so ~1.4

Neither will give the correct 1/f (1/f^n, n:real) noise corner!

## Examples

[nsd_example](examples/nsd_example.py) generates white & brownian noise with a corner frequency of 0.1Hz/1nV:
![nsd_example](examples/white%20&%20brownian%20noise%20corner%200.1Hz%20x%201nV.png?raw=true "nsd_example")

[nsd_example_sps](examples/nsd_example_sps.py) shows how a difference between acquisition time and time between samples effects the NSD:
![nsd_example_sps](examples/white%20&%20white_pink%20(f^-0.5)%20noise%20-%20corner%200.1Hz%20x%201nV.png?raw=true "nsd_example_sps")

[nsd_example_csv](examples/nsd_example_csv.py) imports data from csv:
![nsd_example_csv](examples/K182%20low%20T-EMF%20short%203mV%2020ms%20no%20filter%202022-06-24.csv.png?raw=true "nsd_example_csv")

[nsd_example_smooth](examples/nsd_example_smooth.py) shows different filters:
![nsd_example_smooth](examples/white%20&%20brownian%20noise%20-%20filter%20comparison.png?raw=true "nsd_example_smooth")