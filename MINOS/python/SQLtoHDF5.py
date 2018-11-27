import sqlite3
import numpy as np
import os
import h5py
from MINOS_analysis import database

outPath = '/Volumes/Ian External HD/Node Data/sqlitefiles/MUSE12/'
inPath = '/Volumes/Ian External HD/Node Data/sqlitefiles/MUSE12/'

# list all database files
dbFiles = os.listdir(inPath)
dbFiles = [x for x in dbFiles if 'sqlite3' in x]

for dbFile in dbFiles:

    f = h5py.File(outPath + dbFile.replace('.sqlite3', '.hdf5'), "w")

    # --------------
    # Open Det1 data
    # --------------
    try:
        dataTable = 'Det_3x3_data'
        detDB = database(inPath + dbFile)
        print dbFile
        detDB.getHeader(dataTable)
        header = np.array(detDB.header)
        detDB.getColumn(dataTable, 'Live_Time')
        liveTimes = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Time')
        times = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Spectrum__IntArray')
        spectra = [np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data]
        prefix = dbFile.split('-')[0]
        dset = f.create_dataset(prefix+'_Spectra', (len(spectra), len(spectra[0])), dtype='f', compression="gzip")
        dset[...] = spectra
        dset = f.create_dataset(prefix+'_LiveTimes', (len(liveTimes),), dtype='f', compression="gzip")
        dset[...] = liveTimes
        dset = f.create_dataset(prefix+'_Times', (len(times),), dtype='float64', compression="gzip")
        dset[...] = times

    except:
        pass

    f.close()