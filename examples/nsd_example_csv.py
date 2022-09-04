import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import nsd
import pandas as pd

# change accordingly
data_path = 'C:\\'
data_dir = 'Temp\\'
data_file = 'OUR DATA.csv'

sample_rate = 50 # in Hz
nsd_bins = 2**18 # number of nsd points, needs to be <= datapoints

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
    fig.tight_layout()                              # Adjust spacings w.r.t. figsize
    plt.rcParams['savefig.facecolor']='white'       # set background color for saving, standard is transparent
    #plt.savefig(f'{data_file}.png')                 # change name accordingly

    plt.show()

df = pd.read_table (f'{data_path}{data_dir}{data_file}',    # our file
                    sep = ',',                              # column separator
                    header = 0,                             # none as we define the name ourself
                    skipinitialspace = True,                # convenient
                    skiprows = 0,                           # if there is a header, put the number of rows here
                    #nrows = 100,                           # only take n rows
                    usecols = [0],                          # which column has the values 1st = 0
                    names = ['value'],                      # define name of the column
                    encoding = 'latin1')                    # put the encoding type of the file here

nsd_values = nsd.get(df['value'], sample_rate, nsd_bins)
    
series_rms=nsd.series_rms(df['value'])
nsd_rms=nsd.nsd_rms(nsd_values)
print(f'Series RMS: {series_rms:.2e}, NSD RMS: {nsd_rms:.2e}')

plotnsd(nsd_values[0], nsd_values[1], data_file)    # nsd frequencies, nsd values, name of the series
