from MINOS_analysis import database
import numpy as np
import matplotlib.pyplot as plt
import os

inPath_runestone = '/Users/i6o/MINOS/Rune_Runestone_v1.0.3_Testing/Runestone/'
inPath_rune = '/Users/i6o/MINOS/Rune_Runestone_v1.0.3_Testing/Runestone/'

dbFiles_rune = [x for x in os.listdir(inPath_rune) if '.sqlite3' in x]
dbFiles_runestone = [x for x in os.listdir(inPath_runestone) if '.sqlite3' in x]

spectra_size_comparison = []
total_spectra_comparison = []
cps_comparison = []

for dbFile in dbFiles_rune:
    print dbFile
    detDB_rune = database(inPath_rune + dbFile)
    dataTable = 'Osprey_data'
    detDB_rune.getColumn(dataTable, 'Spectrum__IntArray')
    spectra_rune = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB_rune.data])

    detDB_rune.getColumn(dataTable, 'Time')
    Time_rune = np.array(detDB_rune.data)[:, 0]
    detDB_rune.getColumn(dataTable, 'Live_Time')
    Live_Time_rune = np.array(detDB_rune.data)[:, 0]
    cps_rune = np.sum(spectra_rune,axis=1)
    total_spectra_rune = np.sum(spectra_rune,axis=0)

detDB_runestone = database(inPath_runestone + dbFile)
    dataTable = 'Osprey_data'
    detDB_runestone.getColumn(dataTable, 'Spectrum__IntArray')
    spectra_runestone = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB_runestone.data])
    detDB_runestone.getColumn(dataTable, 'Time')
    Time_runestone = np.array(detDB_runestone.data)[:, 0]
    detDB_runestone.getColumn(dataTable, 'Live_Time')
    Live_Time_runestone = np.array(detDB_runestone.data)[:, 0]
    cps_runestone = np.sum(spectra_runestone,axis=1)
    total_spectra_runestone = np.sum(spectra_runestone,axis=0)

    spectra_size_comparison.append(len(spectra_rune)-len(spectra_runestone))
    total_spectra_comparison.append(total_spectra_rune-total_spectra_runestone)
    cps_comparison.append(cps_rune-cps_runestone)

fig,ax = plt.subplots(nrows=3)
ax[0].plot(spectra_size_comparison,label='array size comparison');ax[0].grid(alpha=0.5)
ax[1].plot(total_spectra_comparison,label='spectra comparison');ax[1].grid(alpha=0.5)
ax[2].plot(cps_comparison,label='cps comparison');ax[2].grid(alpha=0.5)
plt.show()