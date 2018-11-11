import numpy as np
import os
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime as dt

inPath = '/Volumes/IAN USB/Research/'
outPath = inPath + 'figures/'

dbFiles = ['archerair2012-7606_Train_4_Detectors_Bkgrd_2018-11-08T18.25.33.669.hdf5']
Bkgrd_counts = {}
Bkgrd_times = {}
Bkgrd_spectra = {}
Bkgrd_livetimes = {}

for dbFile in dbFiles:
    print 'Reading %s' % dbFile

    hdf5 = h5py.File(inPath + dbFile, 'r')
    name = dbFile.split('_')[4]
    for j in hdf5.keys():
        i = 'MUSE1'
        if 'MUSE1_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Bkgrd_counts[i] = counts
            Bkgrd_spectra[i] = spectra
            first_time_spectra = False

        if 'MUSE1_LiveTimes' in j:
            Bkgrd_livetimes[i] = np.array(hdf5.get(str(j)))
            first_time_livetimes = False

        i = 'MUSE5'
        if 'MUSE5_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Bkgrd_counts[i] = counts
            Bkgrd_spectra[i] = spectra
            first_time_spectra = False

        if 'MUSE5_LiveTimes' in j:
            Bkgrd_livetimes[i] = np.array(hdf5.get(str(j)))
            first_time_livetimes = False

        i = 'MUSE6'
        if 'MUSE6_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Bkgrd_counts[i] = counts
            Bkgrd_spectra[i] = spectra
            first_time_spectra = False

        if 'MUSE6_LiveTimes' in j:
            Bkgrd_livetimes[i] = np.array(hdf5.get(str(j)))
            first_time_livetimes = False

        i = 'MUSE_OTHER'
        if 'MUSE_OTHER_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Bkgrd_counts[i] = counts
            Bkgrd_spectra[i] = spectra

        if 'MUSE_OTHER_LiveTimes' in j:
            Bkgrd_livetimes[i] = np.array(hdf5.get(str(j)))
            first_time_livetimes = False
normed_bkgrd = {}
for i in Bkgrd_spectra:
    normed_bkgrd[i]=np.array(np.sum(Bkgrd_spectra[i],axis=0))/np.sum(Bkgrd_livetimes[i])

Node_counts = {}
Node_spectra = {}

dbFiles = ['archerair2012-7606_Train_4_Detectors_Cs137_x1-2018-11-08T18.42.45.359.hdf5']

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

        i = 'MUSE5'
        if 'MUSE5_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Node_counts[i] = counts
            Node_spectra[i] = spectra

        i = 'MUSE6'
        if 'MUSE6_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Node_counts[i] = counts
            Node_spectra[i] = spectra

        i = 'MUSE_OTHER'
        if 'MUSE_OTHER_Spectra' in j:
            spectra = np.array(hdf5.get(str(j)))
            counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
            Node_counts[i] = counts
            Node_spectra[i] = spectra

MUSE5_source = np.sum(Node_spectra['MUSE5'][142:162],axis=0)/float(162-142)
MUSE6_source = np.sum(Node_spectra['MUSE6'][176:190],axis=0)/float(190-176)
MUSE_OTHER_source = np.sum(Node_spectra['MUSE_OTHER'][206:219],axis=0)/float(219-206)
MUSE1_source = np.sum(Node_spectra['MUSE1'][240:258],axis=0)/float(258-240)

MUSE5_source[np.where(MUSE5_source==0)[0]] =  float('NaN')

energy = np.load('/Volumes/IAN USB/MUSE04_eCal_2x4x16.npy')

fig,ax = plt.subplots(2,2)
ax[0,0].plot(energy,MUSE5_source,'r--',linewidth=0.75,label='Det.1',zorder=2)
ax[0,0].plot(energy,normed_bkgrd['MUSE5'],'b--',linewidth=1.2,zorder=1)
ax[0,0].grid(alpha=0.5)
ax[0,0].set_yscale('log')
ax[0,0].legend(loc='upper right',fancybox=True,shadow=True)
ax[0,0].set_axisbelow(True)
ax[0,0].set_xlabel('Energy (keV)')

ax[1,0].plot(energy,MUSE6_source,'r--',linewidth=0.75,label='Det.2',zorder=2)
ax[1,0].plot(energy,normed_bkgrd['MUSE6'],'b--',linewidth=1.2,zorder=1)
ax[1,0].grid(alpha=0.5)
ax[1,0].set_yscale('log')
ax[1,0].legend(loc='upper right',fancybox=True,shadow=True)
ax[1,0].set_axisbelow(True)
ax[1,0].set_xlabel('Energy (keV)')

ax[0,1].plot(energy,MUSE_OTHER_source,'r--',linewidth=0.75,label='Det.3',zorder=2)
ax[0,1].plot(energy,normed_bkgrd['MUSE_OTHER'],'b--',linewidth=1.2,zorder=1)
ax[0,1].legend(loc='upper right',fancybox=True,shadow=True)
ax[0,1].grid(alpha=0.5)
ax[0,1].set_yscale('log')
ax[0,1].legend(loc='upper right',fancybox=True,shadow=True)
ax[0,1].set_axisbelow(True
ax[0,1].set_xlabel('Energy (keV)')

ax[1,1].plot(energy,MUSE1_source,'r--',linewidth=0.75,label='Det.4',zorder=2)
ax[1,1].plot(energy,normed_bkgrd['MUSE1'],'b--',linewidth=1.2,zorder=1)
ax[1,1].grid(alpha=0.5)
ax[1,1].set_yscale('log')
ax[1,1].legend(loc='upper right',fancybox=True,shadow=True)
ax[1,1].set_axisbelow(True)
ax[1,1].set_xlabel('Energy (keV)')
fig.tight_layout()

plt.show(True)
