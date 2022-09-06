import numpy as np
import nsd
import nsd_helper as nh
from scipy import signal as sp


samples = 2**20                     # number of samples to generate 2**20 = 1Mio
sample_rate = 50                    # in Hz
#nsd_bins = 2**20                    # number of nsd points, optional, default: samples/4
noise_types = ['white', 'brownian'] # see nsd_helper nsd_noise

# generate white + brownian noise with two tones
ts_noise = nh.nsd_noise(samples, noise_types=noise_types)   # seed is optional, default: seed=4
ts_noise += nh.tone(1e-9, 1, sample_rate, samples) + nh.tone(1e-9, 0.01, sample_rate, samples) 

# get the NSD
nsd_noise = nsd.get(ts_noise, sample_rate)  #nsd_bins = number of nsd points, optional, default: samples/4

# compare ACRMS between time series and NSD, should be nearly equal
series_rms = nsd.series_rms(ts_noise)
nsd_rms = nsd.nsd_rms(nsd_noise)
print(f'Series RMS: {series_rms:.2e}, NSD RMS: {nsd_rms:.2e}')
    
# add NSD and NSD with different filters to dictionary
noises = [{'nsd': nsd_noise, 'label': f'target {sample_rate}SPS'}]

noises.append({'nsd': sp.savgol_filter(nsd_noise, window_length=64, polyorder=1, deriv=0, delta=1.0), 'label': f'target {sample_rate}SPS SavGol wl:64 poly:1'})

noises.append({'nsd': nsd.smooth(nsd_noise, 64, filter_function=np.median), 'label': f'target {sample_rate}SPS smoothed median 64pts'})
noises.append({'nsd': nsd.smooth(nsd_noise, 64), 'label': f'target {sample_rate}SPS smoothed mean 64pts'})

# fit
geomspace = np.geomspace(nsd_noise[0][0], nsd_noise[0][-1], num=64)

# linear space
popt, pcov = nsd.fit(nsd_noise)
perr = np.sqrt(np.diag(pcov))   # one standard deviation errors on the parameters
popt_eng = nh.eng(popt)
print (f'Parameter: {popt_eng}, one Std: {nh.eng(perr)}')
noises.append({'nsd': [geomspace, nsd.fit_function(geomspace, *popt)], 'label': f'target {sample_rate}SPS fit to a/f^b + c (least squares): a={popt_eng[0]}, b={popt_eng[1]}, c={popt_eng[2]}'})

# loglog space
popt1, pcov1 = nsd.fit_loglog(nsd_noise)
perr1 = np.sqrt(np.diag(pcov1))   # one standard deviation errors on the parameters
popt1_eng = nh.eng(popt1)
print (f'Parameter: {popt1_eng}, one Std: {nh.eng(perr1)}')
noises.append({'nsd': [geomspace, nsd.fit_function((geomspace), *popt1)], 'label': f'target {sample_rate}SPS fit to a/f^b + c (least squares in loglog space): a={popt1_eng[0]}, b={popt1_eng[1]}, c={popt1_eng[2]}'})

# plot
nh.plot(noises, f'NSD for {" & ".join(noise_types)} noise, {nh.eng(samples)} samples - filter comparison')
print('done')
