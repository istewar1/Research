import numpy as np
import os
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime as dt


inPath = '/Volumes/IAN USB/Research/'
outPath = inPath + 'figures/'


Node_counts = {}
Node_times = {}
Node_spectra = {}
Node_datetimes = {}
Node_livetimes = {}

first_time_spectra = True;
first_time_livetimes = True
first_time_times = True

dbFiles = [x for x in os.listdir(inPath) if '.hdf5' in x]

for dbFile in dbFiles:
    print 'Reading %s' % dbFile

    hdf5 = h5py.File(inPath + dbFile, 'r')
    name = dbFile.split('_')[4]
    for j in hdf5.keys():
        i = 'MUSE1'
        if 'MUSE1_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Node_counts[i] = counts
            Node_spectra[i] = spectra
            first_time_spectra = False

        if 'MUSE1_Times' in j:
            Node_times[i] = np.array(hdf5.get(str(j)))
            first_time_times = False

        if 'MUSE1_LiveTimes' in j:
            Node_livetimes[i] = np.array(hdf5.get(str(j)))
            first_time_livetimes = False

        i = 'MUSE5'
        if 'MUSE5_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Node_counts[i] = counts
            Node_spectra[i] = spectra
            first_time_spectra = False

        if 'MUSE5_Times' in j:
            Node_times[i] = np.array(hdf5.get(str(j)))
            first_time_times = False

        if 'MUSE5_LiveTimes' in j:
            Node_livetimes[i] = np.array(hdf5.get(str(j)))
            first_time_livetimes = False

        i = 'MUSE6'
        if 'MUSE6_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Node_counts[i] = counts
            Node_spectra[i] = spectra
            first_time_spectra = False

        if 'MUSE6_Times' in j:
            Node_times[i] = np.array(hdf5.get(str(j)))
            first_time_times = False

        if 'MUSE6_LiveTimes' in j:
            Node_livetimes[i] = np.array(hdf5.get(str(j)))
            first_time_livetimes = False

        i = 'MUSE_OTHER'
        if 'MUSE_OTHER_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Node_counts[i] = counts
            Node_spectra[i] = spectra
            first_time_spectra = False

        if 'MUSE_OTHER_Times' in j:
            Node_times[i] = np.array(hdf5.get(str(j)))
            first_time_times = False

        if 'MUSE_OTHER_LiveTimes' in j:
            Node_livetimes[i] = np.array(hdf5.get(str(j)))
            first_time_livetimes = False

    for i in ['MUSE1','MUSE5','MUSE6','MUSE_OTHER']:
        Node_datetimes[i] = np.array([dt.datetime.fromtimestamp(x) for x in Node_times[i]])

    fig,ax = plt.subplots(figsize=(9,7),nrows=4,sharex=True)
    ax[0].plot(Node_datetimes['MUSE5'],Node_counts['MUSE5'],linewidth=1.5,label='Det. 1');ax[0].grid(alpha=0.5);ax[0].legend(loc='upper right',fancybox=True,shadow=True)
    ax[1].plot(Node_datetimes['MUSE6'],Node_counts['MUSE6'],'r',linewidth=1.5,label='Det. 2');ax[1].grid(alpha=0.5);ax[1].legend(loc='upper right',fancybox=True,shadow=True)
    ax[2].plot(Node_datetimes['MUSE_OTHER'],Node_counts['MUSE_OTHER'],'k',linewidth=1.5,label='Det. 3');ax[2].grid(alpha=0.5);ax[2].legend(loc='upper right',fancybox=True,shadow=True)
    ax[3].plot(Node_datetimes['MUSE1'],Node_counts['MUSE1'],'g',linewidth=1.5,label='Det. 4');ax[3].grid(alpha=0.5);ax[3].legend(loc='upper right',fancybox=True,shadow=True)
    ax[2].set_ylabel('Counts per 100-msec')
    ax[2] = plt.gca()
    ax[2].xaxis.set_major_formatter(md.DateFormatter('%H:%M:%S'))
    ax[3].set_xlabel('Time (Hour:Minute:Second)')
    title_string = name + '-source placed on Train Car'
    ax[0].set_title(title_string)
    for ax1 in fig.axes:
        matplotlib.pyplot.sca(ax1)
        plt.xticks(rotation=25)
    plt.subplots_adjust(hspace=0.05)
    plt.savefig(outPath+name+'.png',dpi=600)
    
plt.show(False)
