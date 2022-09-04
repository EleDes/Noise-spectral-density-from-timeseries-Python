import numpy as np
import colorednoise as cn
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import nsd

sample_rate = 50 # in Hz
nsd_bins = 2**20 # number of nsd points

def get_noise():
    np.random.seed(4) # same result for each run

    samples = 2**22 # number of samples to generate 2**20 = 1Mio

    beta = 0 # the exponent - 0: gaussian, 1: pink (1/f), 2: brown (1/f**2), -1: blue, -2: violet
    white = cn.powerlaw_psd_gaussian(beta, samples)/2e8    #/2e8 for beta = 0, /5.5e8 for beta = 0.5, /8.5e8 for beta = 1, /4e7 for beta = 2 and 0.1Hzx1nV

    beta = 2 # the exponent - 0: gaussian, 1: pink (1/f), 2: brown (1/f**2), -1: blue, -2: violet
    brownian = cn.powerlaw_psd_gaussian(beta, samples)/4e7    #/2e8 for beta = 0, /5.5e8 for beta = 0.5, /8.5e8 for beta = 1, /4e7 for beta = 2 and 0.1Hzx1nV

    brownian_white = white + brownian
    return brownian_white

def plotnsd(freq, values, label):
    # Plot
    fig, ax = plt.subplots(figsize=(15, 10))
    plt.plot(freq, values, '-o', ms=0.5, lw=0.5, label=label)
    plt.xlabel('frequency in Hz')
    plt.ylabel(r'noise in V/$\sqrt{Hz}$')
    plt.loglog()
    plt.grid(True, which="both")
    formatter1 = EngFormatter(sep="")
    plt.gca().xaxis.set_major_formatter(formatter1)
    plt.gca().yaxis.set_major_formatter(formatter1)
    plt.gca().yaxis.set_minor_formatter(formatter1)
    plt.legend()

    #save plot
    fig.tight_layout()                                             # Adjust spacings w.r.t. figsize
    plt.rcParams['savefig.facecolor']='white'                      # set background color for saving, standard is transparent
    #plt.savefig(f'white & brownian noise corner 0.1Hz x 1nV.png')  # change name accordingly

    plt.show()

ts_noise = get_noise()

nsd_noise = nsd.get(ts_noise, sample_rate, nsd_bins)
    
series_rms=nsd.series_rms(ts_noise)
nsd_rms=nsd.nsd_rms(nsd_noise)
print(f'Series RMS: {series_rms:.2e}, NSD RMS: {nsd_rms:.2e}')
    
plotnsd(nsd_noise[0], nsd_noise[1], 'white and brownian noise - corner: 0.1Hz/1nV')
