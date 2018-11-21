import h5py
import numpy as np
import os
import datetime as dt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as md
from sklearn.decomposition import PCA

node = 'MUSE04'
inPath = '/Volumes/Ian External HD/Node Data/sqlitefiles/'+node+'/'

# list all database files
dbFiles = os.listdir(inPath)
dbFiles = [x for x in dbFiles if '.hdf5' in x]

times = []
spectra = []
datetimes = np.array([])
cps = np.array([])

for dbFile in dbFiles:
    print dbFile
    hdf5 = h5py.File(inPath + dbFile,'r')
    for j in hdf5.keys():
        if '3x3Spectra' in j:
            spectra.append(np.array(hdf5.get(str(j))))
        if '3x3Time' in j:
            times.append(np.array(hdf5.get(str(j))))
spectra = spectra[0]
times = times[0]
i0 = 0
for i in range(1200,len(spectra),600):
    print i
    fig,ax = plt.subplots(1,1)
    spectra1 = np.sum(np.array(spectra)[i0:i],axis=0)
    ax.plot(spectra1)
    title_string = '5min from %.1f to %.1f'%(times[i0],times[i])
    ax.legend()
    ax.grid(alpha=0.5)
    ax.set_yscale('log')
    i0=i+1
    save_title = node+'_%f_to_%f.png'%(times[i0],times[i])
    plt.savefig(inPath+'Figures/'+save_title,dpi=100)
    plt.close('all')
