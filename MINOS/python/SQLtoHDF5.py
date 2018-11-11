import sqlite3
import numpy as np
import os
import h5py
from MINOS_analysis import database

outPath = '/Volumes/IAN USB/Research/'
inPath = '/Volumes/IAN USB/Research/'

# list all database files
dbFiles = os.listdir(inPath)
dbFiles = [x for x in dbFiles if 'sqlite3' in x]
dbFiles = ['archerair2012-7606_Train_4_Detectors_Ba133_10micro-2018-11-08T18.36.20.086.sqlite3','archerair2012-7606_Train_4_Detectors_Bkgrd-2018-11-08T18.25.33.669.sqlite3','archerair2012-7606_Train_4_Detectors_Co60-2018-11-08T18.48.32.419.sqlite3','archerair2012-7606_Train_4_Detectors_Cs137_x1-2018-11-08T18.42.45.359.sqlite3']

for dbFile in dbFiles:

    f = h5py.File(outPath + dbFile.replace('.sqlite3', '.hdf5'), "w")

    # --------------
    # Open Det1 data
    # --------------
    try:
        dataTable = 'Det_2x4x16_MUSE_1_data'
        detDB = database(inPath + dbFile)

        detDB.getHeader(dataTable)
        header = np.array(detDB.header)
        detDB.getColumn(dataTable, 'Live_Time')
        liveTimes = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Time')
        times = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Spectrum__IntArray')
        spectra = [np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data]
        dset = f.create_dataset('MUSE1_Spectra', (len(spectra), len(spectra[0])), dtype='f', compression="gzip")
        dset[...] = spectra
        dset = f.create_dataset('MUSE1_LiveTimes', (len(liveTimes),), dtype='f', compression="gzip")
        dset[...] = liveTimes
        dset = f.create_dataset('MUSE1_Times', (len(times),), dtype='float64', compression="gzip")
        dset[...] = times

        dataTable = 'Det_2x4x16_MUSE_5_data'
        detDB.getHeader(dataTable)
        header = np.array(detDB.header)
        detDB.getColumn(dataTable, 'Live_Time')
        liveTimes = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Time')
        times = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Spectrum__IntArray')
        spectra = [np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data]
        dset = f.create_dataset('MUSE5_Spectra', (len(spectra), len(spectra[0])), dtype='f', compression="gzip")
        dset[...] = spectra
        dset = f.create_dataset('MUSE5_LiveTimes', (len(liveTimes),), dtype='f', compression="gzip")
        dset[...] = liveTimes
        dset = f.create_dataset('MUSE5_Times', (len(times),), dtype='float64', compression="gzip")
        dset[...] = times

        dataTable = 'Det_2x4x16_MUSE_6_data'
        detDB.getHeader(dataTable)
        header = np.array(detDB.header)
        detDB.getColumn(dataTable, 'Live_Time')
        liveTimes = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Time')
        times = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Spectrum__IntArray')
        spectra = [np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data]
        dset = f.create_dataset('MUSE6_Spectra', (len(spectra), len(spectra[0])), dtype='f', compression="gzip")
        dset[...] = spectra
        dset = f.create_dataset('MUSE6_LiveTimes', (len(liveTimes),), dtype='f', compression="gzip")
        dset[...] = liveTimes
        dset = f.create_dataset('MUSE6_Times', (len(times),), dtype='float64', compression="gzip")
        dset[...] = times

        dataTable = 'Det_2x4x16_MUSE_other_data'
        detDB.getHeader(dataTable)
        header = np.array(detDB.header)
        detDB.getColumn(dataTable, 'Live_Time')
        liveTimes = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Time')
        times = np.array(detDB.data)[:, 0]
        detDB.getColumn(dataTable, 'Spectrum__IntArray')
        spectra = [np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data]
        dset = f.create_dataset('MUSE_OTHER_Spectra', (len(spectra), len(spectra[0])), dtype='f', compression="gzip")
        dset[...] = spectra
        dset = f.create_dataset('MUSE_OTHER_LiveTimes', (len(liveTimes),), dtype='f', compression="gzip")
        dset[...] = liveTimes
        dset = f.create_dataset('MUSE_OTHER_Times', (len(times),), dtype='float64', compression="gzip")
        dset[...] = times

    except:
        pass

    f.close()