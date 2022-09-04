import numpy as np
import colorednoise as cn
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import nsd

samples = 2**22 # number of samples to generate 2**20 = 1Mio - chart was generated with 2**26
sample_rate = 50 # in Hz
nsd_bins = 2**18 # number of nsd points - chart was generated with 2**22

def get_noise(samples):
    np.random.seed(4) # same result for each run

    beta = 0 # the exponent - 0: gaussian, 1: pink (1/f), 2: brown (1/f**2), -1: blue, -2: violet
    white = cn.powerlaw_psd_gaussian(beta, samples)/2e8    #/2e8 for beta = 0, /5.5e8 for beta = 0.5, /8.5e8 for beta = 1, /4e7 for beta = 2 and 0.1Hzx1nV

    beta = 0.5 # the exponent - 0: gaussian, 1: pink (1/f), 2: brown (1/f**2), -1: blue, -2: violet
    brownian = cn.powerlaw_psd_gaussian(beta, samples)/5.5e8    #/2e8 for beta = 0, /5.5e8 for beta = 0.5, /8.5e8 for beta = 1, /4e7 for beta = 2 and 0.1Hzx1nV

    brownian_white = white + brownian
    return brownian_white

def tones(rms, freq, fs, samples):
    ur = []
    for i in range (0, samples):
        t = i / fs
        ur.append(rms * 2**0.5 * np.cos(2 * np.pi * freq * t))
    return np.array(ur)

def plotnsd(values_list, title):
    # Plot
    fig, ax = plt.subplots(figsize=(15, 10))
    plt.title(title, fontsize=13)
    for values in values_list:
        plt.plot(values['nsd'][0], values['nsd'][1], '-', ms=1, lw=1, alpha=0.5, label=values['label'])
    #plt.plot(values[0], values[1], '-o', ms=0.5, lw=0.5, label=label)
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
    fig.tight_layout()                                                          # Adjust spacings w.r.t. figsize
    plt.rcParams['savefig.facecolor']='white'                                   # set background color for saving, standard is transparent
    #plt.savefig(f'white & white_pink (f^-0.5) noise - corner 0.1Hz x 1nV.png')  # change name accordingly

    plt.show()

ts_noise = get_noise(samples)

ts_noise += tones(6e-11, 1, sample_rate, samples) + tones(6e-11, 0.01, sample_rate, samples) 

nsd_noise = nsd.get(ts_noise, sample_rate, nsd_bins)

nth = 4
nsd_noise_nth = [nsd.get(ts_noise[::nth], sample_rate, int(nsd_bins/nth))]
nsd_noise_nth.append(nsd.get(ts_noise[::nth], sample_rate/nth, int(nsd_bins/nth)))
nsd_noise_nth.append([nsd_noise_nth[-2][0]/nth, nsd_noise_nth[-2][1]])
nsd_noise_nth.append([nsd_noise_nth[-2][0], nsd_noise_nth[-2][1]/nth**0.5])
    
series_rms=nsd.series_rms(ts_noise)
nsd_rms=nsd.nsd_rms(nsd_noise)
print(f'Series RMS: {series_rms:.2e}, NSD RMS: {nsd_rms:.2e}')
    
noises = [{'nsd': nsd_noise, 'label': f'target {sample_rate}SPS'}]
noises.append({'nsd': nsd_noise_nth[0], 'label': f'every {nth}th {sample_rate}SPS'})
noises.append({'nsd': nsd_noise_nth[1], 'label': f'every {nth}th {sample_rate/nth}SPS'})
noises.append({'nsd': nsd_noise_nth[2], 'label': f'every {nth}th {sample_rate}SPS "corrected"'})
noises.append({'nsd': nsd_noise_nth[3], 'label': f'every {nth}th {sample_rate/nth}SPS "corrected"'})

plotnsd(noises, r'White & white/pink (1/$\sqrt{f}$) noise - corner: 0.1Hz/1nV')

print('done')
