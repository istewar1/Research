import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime as dt
from itertools import groupby
from operator import itemgetter
import h5py
import os
import pandas as pd
import matplotlib
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
inPaths = ['/Volumes/Ian External HD/hdf5Files_12112018/Figures/']
length = np.array([])
for inPath in inPaths:
    dbFiles = [x for x in os.listdir(inPath) if 'Alarm_attributes.csv' in x] # location of csv
    for dbFile in dbFiles:
        print dbFile
        df = pd.read_csv(inPath+dbFile)
        l = df['posix']
        consecutive_list = []
        for k, g in groupby(enumerate(l), lambda (i, x): i - x):
            consecutive_list.append(map(itemgetter(1), g))
        consecutive_array = np.array(consecutive_list)
        l =[]
        for alarm_array in consecutive_array:
            if len(alarm_array) > 1:
                l.append(len(alarm_array))
        length = np.concatenate((length,np.array(l)))
print len(length)
#
fig, axs = plt.subplots(1, 1, tight_layout=True)
N, bins, patches = axs.hist(length, bins=np.arange(2,50+1,1),edgecolor='black', linewidth=1,density=True)
cmap = plt.cm.inferno
fracs = N / N.max()
norm = colors.Normalize(fracs.min(), fracs.max())
for thisfrac, thispatch in zip(fracs, patches):
    color = plt.cm.plasma(norm(thisfrac))
    thispatch.set_facecolor(color)
axs.set_xticks(np.arange(0,51,2))
plt.tick_params(axis='x', which='major', labelsize=8)
plt.yscale('linear')
plt.grid(alpha=0.5,which='both')
axs.set_ylabel('Normalized Frequency (%)')
axs.set_xlabel('Alarm Length (seconds)')
axs.set_axisbelow(True)
axs.yaxis.set_major_formatter(PercentFormatter(xmax=1))
#
cumalutive_sum = []; first = True; counter = 0
for i in N:
    if first:
        cumalutive_sum.append(i)
        first = False
    else:
        cumalutive_sum.append(i+cumalutive_sum[counter])
        counter+=1
axs2 = axs.twinx()
axs2.step(np.arange(2.5,50.5,1),cumalutive_sum,color='r',linestyle='--')
axs2.set_ylabel('Cumulative Percentage (%)', color='r')
axs2.tick_params('y', colors='r')
axs2.yaxis.set_major_formatter(PercentFormatter(xmax=1))
#
plt.savefig('/Volumes/Ian External HD/Node Data/figures/ORNL_histogram_alarm_length.png',dpi=600)
plt.show(False)
