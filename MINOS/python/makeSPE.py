import numpy as np
import os
import matplotlib.pyplot as plt
from MINOS_analysis import database

dbFile = '/Volumes/IAN USB/MUSE04-2018-10-25T14.19.34.827.sqlite3'

detDB = database(dbFile)

dataTable = 'Det_2x4x16_data'
detDB.getColumn(dataTable, 'Spectrum__IntArray')
spectra_MUSE01 = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data])
detDB.getColumn(dataTable, 'Time')
Time_MUSE01 = np.array(detDB.data)[:, 0]
detDB.getColumn(dataTable, 'Live_Time')
Live_Time_MUSE01 = np.array(detDB.data)[:, 0]

times = [[19000,27500],[28500,31500],[32700,39200]]
for i in times:
    time0 = i[0]
    time1 = i[-1]
    spectra = np.sum(spectra_MUSE01[time0:time1],axis=0)
    fig,ax = plt.subplots()
    ax.plot(spectra)
    ax.set_yscale('log')
    ax.grid(alpha=0)
plt.show()