import h5py
import numpy as np
import os
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as md
from sklearn.decomposition import PCA

inPath = '/Volumes/Ian External HD/Node Data/sqlitefiles/MUSE01/'

# list all database files
dbFiles = os.listdir(inPath)
dbFiles = [x for x in dbFiles if '.hdf5' in x]

times = []
#spectra = np.array([]).reshape(-1)
spectra = []
datetimes = np.array([])
cps = np.array([])

dbFiles = [x for x in dbFiles if ('T03' in x)or('T07' in x)]

for dbFile in dbFiles:
    print dbFile
    hdf5 = h5py.File(inPath + dbFile,'r')
    for j in hdf5.keys():
        if '_Spectra' in j:
            #spectra = np.concatenate((spectra,np.array(hdf5.get(str(j))).reshape(-1)))
            spectra.append(np.array(hdf5.get(str(j))))
            c = np.sum(np.array(hdf5.get(str(j))),axis=1)
            cps = np.concatenate((cps,c))
        if '_Times' in j:
            t = np.array(hdf5.get(str(j)))
            times.append(t)
            d = np.array([dt.datetime.fromtimestamp(x) for x in t])
            datetimes = np.concatenate((datetimes,d))

pca = PCA()
s = np.asarray(spectra)
pca.fit(spectra)
plt.figure()
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance');

fig,ax = plt.subplots(1,1)
cps = np.add.reduceat(cps, np.arange(0, len(cps), 10))
ax.plot(datetimes[::10],cps)
ax.grid(alpha=0.5)
ax = plt.gca()
ax.xaxis.set_major_formatter(md.DateFormatter('%m/%d %H:%M:%S'))
ax.set_xlabel('Time (month/day hour:minute)')
ax.set_ylabel('CPS')
for ax in fig.axes:
    plt.sca(ax)
    plt.xticks(rotation=25)
fig.tight_layout()
plt.show()