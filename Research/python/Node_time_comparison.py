import numpy as np
import os
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime as dt
import time
from itertools import groupby
from operator import itemgetter
from scipy.optimize import curve_fit

directory_path = '/Volumes/Ian External HD/Node Data/'
csv_path = '/Volumes/Ian External HD/Node Data/'
nodes = [x for x in os.listdir(directory_path) if 'MUSE' in x]
outPath = directory_path + 'figures/'

# Y-min/Y-max for CPS plots
counts_range = {'MUSE01': [1200, 2100], 'MUSE04': [1200, 2500], \
                'MUSE06': [700, 1800], 'MUSE10': [1200, 2400], \
                'MUSE11': [1600, 2900], 'MUSE12': [1200, 2400]}

nodes = ['MUSE01','MUSE04','MUSE06','MUSE10','MUSE11','MUSE12']
Node_counts = {}
Node_times = {}
Node_spectra = {}
Node_datetimes = {}
Node_livetimes = {}
if False:
    for i in nodes:
        print 'Reading %s' % i
        inPath = directory_path + str(i) + '/'  # +'/hdf5Files/'
        dbFiles = [x for x in os.listdir(inPath) if '_ALG' not in x]
        dbFiles = [x for x in dbFiles if (
                    ('09-26T' in x) or ('09-27T' in x) or ('09-28T' in x) or ('09-29T' in x) or ('09-30T' in x) or (
                        '09-31T' in x) or ('10-01T' in x))]
        first_time_times = True;
        first_time_spectra = True;
        first_time_livetimes = True

        for dbFile in dbFiles:

            hdf5 = h5py.File(inPath + dbFile, 'r')
            print dbFile
            for j in hdf5.keys():

                if '2x4x16Spectra' in j:
                    spectra = np.array(hdf5.get(str(j)))
                    counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
                    if not first_time_spectra:
                        Node_counts[i] = np.concatenate([Node_counts[i], counts])
                        Node_spectra[i] = np.concatenate([Node_spectra[i], spectra])
                    else:
                        Node_counts[i] = counts
                        Node_spectra[i] = spectra
                        first_time_spectra = False

                if '2x4x16Times' in j:
                    if not first_time_times:
                        Node_times[i] = np.concatenate([Node_times[i], np.array(hdf5.get(str(j)))])
                    else:
                        Node_times[i] = np.array(hdf5.get(str(j)))
                        first_time_times = False

                if '2x4x16LiveTimes' in j:
                    if not first_time_livetimes:
                        Node_livetimes[i] = np.concatenate([Node_livetimes[i], np.array(hdf5.get(str(j)))])
                    else:
                        Node_livetimes[i] = np.array(hdf5.get(str(j)))
                        first_time_livetimes = False

            Node_datetimes[i] = np.array([dt.datetime.fromtimestamp(x) for x in Node_times[i]])
    for i in Node_counts:
        np.save(outPath+i+'_counts.npy',Node_counts[i])
        np.save(outPath+i+'_datetimes.npy',Node_datetimes[i])
        np.save(outPath+i+'_spectra.npy',Node_spectra[i])

if True:
    for i in nodes:
        Node_datetimes[i] = np.load(directory_path+'npy_arrays/'+i+'_datetimes.npy')
        Node_counts[i] = np.load(directory_path+'npy_arrays/'+i+'_counts.npy')
        Node_spectra[i] = np.load(directory_path+'npy_arrays/'+i+'_spectra.npy')

node_colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']

fig, ax = plt.subplots(figsize=(12,8), nrows=6, sharex=True)
for i in range(6):
    node = nodes[i]
    label = ['Node-1','Node-4','Node-6','Node-10','Node-11','Node-12']
    ax[i].plot( Node_counts[node], c=node_colors[i], linewidth=0.75, label=label[i])
    ax[i].set_ylim(counts_range[node][0], counts_range[node][1])
    ax[i].legend(loc='upper right',fancybox=False,shadow=True)
    ax[i].grid(alpha=0.5, linewidth=0.5)
    #ax[i].set_xlim(Node_datetimes[node][0], Node_datetimes[node][-1])
    ax[i].set_ylabel('CPS')


for ax in fig.axes:
    matplotlib.pyplot.sca(ax)
    plt.xticks(rotation=25)
plt.subplots_adjust(hspace=0.05)
#plt.savefig(outPath+'test.png',dpi=100)
plt.close('all')

'''
Creating spectra during alarm (hard-coded for specific known array indexes)
'''
energy = np.load('/Volumes/Ian External HD/Node Data/MUSE04/energy_pairs/eCal_2x4x16.npy')
Alarm_values = {'MUSE01':[190389,190394], \
                'MUSE04': [195977, 195985], \
                'MUSE06': [201422, 201433], \
                'MUSE10': [194634, 194641]}

fig,ax = plt.subplots()
for i in Alarm_values:
    indexes = Alarm_values[i]
    print indexes
    print i
    ax.plot(energy,np.sum(np.array(Node_spectra[i])[indexes[0]:indexes[1]],axis=0),label=str(i),c=node_colors[i])
ax.grid(alpha=0.5, linewidth=0.5)
ax.set_yscale('log')
ax.set_axisbelow(True)
ax.legend(loc='upper left',fancybox=False,shadow=True)


plt.show(True)
