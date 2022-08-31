#!/usr/bin/env python
# -*- coding: utf-8 -*-
# developed & tested with Python 3.9

import unittest
import nsd
import numpy as np
import colorednoise as cn

class Test_test_NSD(unittest.TestCase):
	#np.random.seed(2)							# same result for each run

	def set_up(self):
		pass
	
	def test_gaussian_noise(self):
		beta = 0 # the exponent
		samples = 2**18 # number of samples to generate
		gaussian_noise = cn.powerlaw_psd_gaussian(beta, samples)
		gn_rms = NSD.series_rms(gaussian_noise)

		NSD.fs = 1
		nsd_gaussian_noise = NSD.nsd(gaussian_noise, NSD.fs)
		nsd_gn_rms = NSD.nsd_to_rms(nsd_gaussian_noise)

		#plotnsd(nsd_gaussian_noise[0], nsd_gaussian_noise[1])
		self.assertAlmostEqual(gn_rms, nsd_gn_rms, delta=0.01) # to 1% due to inaccuracies
		print(f'RMS Deltafactor fs {NSD.fs}Hz: {nsd_gn_rms/gn_rms-1}')

		NSD.fs = 1E3
		nsd_gaussian_noise = NSD.nsd(gaussian_noise, NSD.fs)
		nsd_gn_rms = NSD.nsd_to_rms(nsd_gaussian_noise)

		#plotnsd(nsd_gaussian_noise[0], nsd_gaussian_noise[1])
		self.assertAlmostEqual(gn_rms, nsd_gn_rms, delta=0.01) # to 1% due to inaccuracies
		print(f'RMS Deltafactor fs {NSD.fs}Hz: {nsd_gn_rms/gn_rms-1}')

	def test_tones(self):
		rms = 1
		frequency = 2**18
		samples = 2**20 # number of samples to generate
		tone = tones(rms, frequency, samples)
		
		#plot(tone)

		NSD.fs = 10
		NSD.nperseg = 2**12 # 2**12=4096
		NSD.noverlap = NSD.nperseg * 0.85
		NSD.window = NSD.window_dpss(2**13)
		#nsd_tones = NSD.nsd(tone, NSD.fs)
		nsd_tones = NSD.nsd_HFT90D(tone, NSD.fs, NSD.nperseg)
		#nsd_HFT90D(values, fs, n_fft)
		#nsd_index =	np.where(nsd_tones[0] == frequency)#[0][0]

		#get rms value of input
		tone_rms = NSD.series_rms(tone)
		
		#get rms value from NSD
		nsd_rms = NSD.nsd_to_rms(nsd_tones)

		#plot(nsd_tones[0])
		#plot(nsd_tones[1])
		plotnsd(nsd_tones[0], nsd_tones[1], 'one tone')
		#self.assertAlmostEqual(tone_rms, nsd_rms, delta=0.01) # to 1% due to inaccuracies
		print(f'RMS Deltafactor: {nsd_rms/rms-1}')

	def test_two_tones_GH(self):
		two_tones = two_tones_GH()
		
		#plot(tone)

		NSD.fs = 10000
		NSD.nperseg = 2**16 #4096
		NSD.noverlap = NSD.nperseg * 0.85
		NSD.window = NSD.window_dpss(NSD.nperseg)
		nsd_tones = NSD.nsd(two_tones, NSD.fs)
		#nsd_tones = NSD.nsd_HFT90D(two_tones, NSD.fs, NSD.nperseg)
		#nsd_index =	np.where(nsd_tones[0] == frequency)#[0][0]

		#get rms value of input
		tone_rms = NSD.series_rms(two_tones)
		
		#get rms value from NSD
		nsd_rms = NSD.nsd_to_rms(nsd_tones)

		#plot(nsd_tones[0])
		#plot(nsd_tones[1])
		plotnsd(nsd_tones[0], nsd_tones[1], 'two tones GH')
		#self.assertAlmostEqual(tone_rms, nsd_rms, delta=0.01) # to 1% due to inaccuracies
		print(f'RMS Deltafactor: {nsd_rms/tone_rms-1}')

def tones(rms, freq, fs):
	# sampling interval
	ts = 1.0/fs
	t = np.arange(0,1,ts)
	return rms*2**0.5*np.sin(2*np.pi*freq*t)

# from https://holometer.fnal.gov/GH_FFT.pdf#page=22
# Two tones with 1234Hz/2Vrms, 2500.2157Hz/0.707Vrms (1Vpk), noise floor ~4.08µV/Hz^0.5
def two_tones_GH():
	import math
	fs = 10000 # sampling frequency [Hz] */
	f1 = 1234 # first signal frequency [Hz] */
	amp1 = 2 * 2 ** 0.5 # 2 Vrms */
	f2 = 2500.2157 # second signal frequency [Hz] */
	amp2 = 1 # 0.707 Vrms */
	ulsb = 1E-3 # Value of 1 LSB in Volt */
	#int i
	#double t, u, ur
	ur = []
	for i in range (0, 1000000):
		t = i / fs
		u = amp1 * np.sin(2 * np.pi * f1 * t) + amp2 * np.sin(2 * np.pi * f2 * t)
		ur.append(math.floor(u / ulsb + 0.5) * ulsb) # Rounding */
		#printf ("%10.6f %8.5f\n", t, ur); # ASCII output */
		#fwrite (&ur, sizeof (double), 1, stdout); # alternative binary output */
	return np.array(ur)

def plotnsd(freq, values, label):
	import matplotlib.pyplot as plt
	from matplotlib.ticker import EngFormatter
	# Plot
	plt.subplots(figsize=(15, 10))
	plt.plot(freq, values, '-o', ms=0.5, lw=0.5, label=label)
	#plt.plot(freq, values, markersize=1)
	#plt.ylim(0.4)
	#plt.xlim([0.001, fs/2*0.99])
	plt.xlabel('frequency in Hz')
	plt.ylabel(r'Value noise (NSD) in ValueUnit/$\sqrt{Hz}$')
	#plt.semilogx()
	plt.loglog()
	plt.grid(True, which="both")
	formatter1 = EngFormatter(sep="")
	plt.gca().xaxis.set_major_formatter(formatter1)
	#plt.gca().xaxis.set_minor_formatter(formatter1)
	plt.gca().yaxis.set_major_formatter(formatter1)
	plt.gca().yaxis.set_minor_formatter(formatter1)
	plt.legend()
	plt.show()
	#plt.savefig('NSD of gaussian noise.png')

def plot(values):
	import matplotlib.pyplot as plt
	from matplotlib.ticker import EngFormatter
	# Plot
	plt.subplots(figsize=(15, 10))
	plt.plot(values, '-o', markersize=1, linewidth=0.5)
	#plt.plot(values, markersize=1, linewidth=0.5)
	#plt.ylim(0.4)
	#plt.xlim([0.001, fs/2*0.99])
	#plt.xlabel('frequency in Hz')
	plt.ylabel(r'Value')
	#plt.semilogx()
	#plt.loglog()
	plt.grid(True, which="both")
	formatter1 = EngFormatter(sep="")
	plt.gca().xaxis.set_major_formatter(formatter1)
	plt.gca().xaxis.set_minor_formatter(formatter1)
	plt.gca().yaxis.set_major_formatter(formatter1)
	plt.gca().yaxis.set_minor_formatter(formatter1)
	plt.legend()
	plt.show()
	#plt.savefig('NSD of gaussian noise.png')

if __name__ == '__main__':
    unittest.main()