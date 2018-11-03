import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as md
import os
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
from scipy.interpolate import UnivariateSpline
import sqlite3
from numba import jit

def generate_calibration_axis(bins1, energies1):
    """
    Perform a Akima spline interpolation between the given bins and energies, and a linear polynomial
    extrapolation beyond the known peaks, which is fit between the last two known peaks
    :param bins: Bins at which each energy is located
    :param energies: Energies being used for calibration
    :param display: set this to (characterizerinstance).NONE, (future)
    :return: the new x_axis energy data from the spline fit
    """
    bins = np.array(bins1)
    energies = np.array(energies1)
    bins = np.append(bins, 0)
    energies = np.append(energies, 0.0)
    # sort the bins and energies
    bins, energies = zip(*sorted(zip(bins, energies)))
    print("Using Bins: " + str(bins))
    print("Using Energies: " + str(energies))
    akima_fit = interpolate.Akima1DInterpolator(x=bins, y=energies)
    ext = np.polyfit([bins[-2], bins[-1]], [energies[-2], energies[-1]], 1)
    interpolated_energy_data = akima_fit(np.arange(0, bins[-1]))
    extrapolated_energy_data = np.polyval(ext, np.arange(bins[-1], 2048))
    energy_data = np.concatenate([interpolated_energy_data, extrapolated_energy_data])
    print(len(energy_data))
    return energy_data


global database

class database(object):  # database class

    def __init__(self, DBfilename):
        self.conn = sqlite3.connect(DBfilename)
        self.conn.row_factory = sqlite3.Row

    def getTables(self):
        c = self.conn.cursor()
        c.execute("select name from sqlite_master where type = 'table'")
        self.tables = c.fetchall()
        c.close()

    def getHeader(self, tblname):
        c = self.conn.cursor()
        c.execute("select rowid,* from %s" % (tblname))
        self.header = np.array(c.fetchone().keys())
        c.close()

    def getAllData(self, tblname):
        c = self.conn.cursor()
        c.execute("select rowid,* from %s order by Time asc" % (tblname))
        self.data = self.tableToArray(c)
        c.close()

    def getSomeData(self, tblname, label, t0, tf):
        c = self.conn.cursor()
        c.execute("select rowid,* from %s where %s>=%.1f and %s<%.1f" % (tblname, label, t0, label, tf))
        self.data = self.tableToArray(c)
        return

    def getSomeDataIf(self, tblname, label, t0, tf, x):
        c = self.conn.cursor()
        c.execute(
            "select rowid,* from %s where %s>=%.1f and %s<%.1f and Gate1=%d " % (tblname, label, t0, label, tf, x))
        self.data = self.tableToArray(c)
        return

    def getColumn(self, tbleName, column):
        c = self.conn.cursor()
        c.execute("select %s from %s" % (column, tbleName))
        self.data = c.fetchall()
        c.close()

    def getFirstAndLastRow(self, tblname):
        c = self.conn.cursor()
        c.execute("select rowid,* from %s order by rowid asc limit 1" % (tblname))
        firstRow = self.tableToArray(c)
        c.execute("select rowid,* from %s order by rowid desc limit 1" % (tblname))
        lastRow = self.tableToArray(c)
        return firstRow, lastRow

    def tableToArray(self, c):
        rows = c.fetchall()
        data = []
        for i in range(len(rows)):
            data.append(list(rows[i]))
        return data

    def closeConn(self):
        self.conn.close()
        self.data = []


# Gaussian + line
def peakFunc(x, a, b, c, mu, sigma):
    '''
    Peak fit function with Gaussian on sloped line to account
    for down scatter contribution from higher E photons.
    '''
    return (
        # line
            a + c * x +
            # photopeak
            b * np.exp(-(x - mu) ** 2 / (2.0 * sigma ** 2))
    )


def getIsotopeEnergies_large(isotope):
    '''
    Description: Provides energy and ROI for isotope.
    Param: isotope = string (e.g. 'Cs-137')
    '''
    peaks = [];
    peak_windows_lower = [];
    peak_windows_upper = []
    if isotope == "Cs-137":
        peaks.append(662.0)
        peak_windows_lower.append(570)
        peak_windows_upper.append(780)
    elif isotope == "Co-60":
        peaks.append(1173.0)
        peak_windows_lower.append(1050)
        peak_windows_upper.append(1240)
        peaks.append(1332.0)
        peak_windows_lower.append(1250)
        peak_windows_upper.append(1420)
    elif isotope == "Ba-133":
        peaks.append(80.0)
        peak_windows_lower.append(60)
        peak_windows_upper.append(100)
        peaks.append(356.0)
        peak_windows_lower.append(320)
        peak_windows_upper.append(450)
    elif isotope == "Co-57":
        peaks.append(122.0)
        peak_windows_lower.append(110)
        peak_windows_upper.append(135)
    elif isotope == "K-40":
        peaks.append(1461.0)
        peak_windows_lower.append(1280)
        peak_windows_upper.append(1600)
    elif isotope == "Background":
        peaks.append(1461.0)
        peak_windows_lower.append(1280)
        peak_windows_upper.append(1600)
    elif isotope == "Th-232":
        peaks.append(2614.0)
        peak_windows_lower.append(2430)
        peak_windows_upper.append(2800)
    elif isotope == "Eu-152":
        peaks.append(40.118)
        peak_windows_lower.append(25)
        peak_windows_upper.append(55)
        peaks.append(121.78)
        peak_windows_lower.append(100)
        peak_windows_upper.append(140)
        peaks.append(344.29)
        peak_windows_lower.append(300)
        peak_windows_upper.append(390)
        peaks.append(778.904)
        peak_windows_lower.append(710)
        peak_windows_upper.append(850)
        peaks.append(1408.006)
        peak_windows_lower.append(1300)
        peak_windows_upper.append(1550)
    elif isotope == "Am-241":
        peaks.append(59.5)
        peak_windows_lower.append(42)
        peak_windows_upper.append(75)
    elif isotope == "Ho-166":
        peaks.append(184.41)
        peak_windows_lower.append(160)
        peak_windows_upper.append(220)
        peaks.append(280.46)
        peak_windows_lower.append(250)
        peak_windows_upper.append(320)
        peaks.append(410.95)
        peak_windows_lower.append(375)
        peak_windows_upper.append(450)
        peaks.append(711.68)
        peak_windows_lower.append(680)
        peak_windows_upper.append(730)
        peaks.append(810.29)
        peak_windows_lower.append(785)
        peak_windows_upper.append(850)
    return peaks, peak_windows_lower, peak_windows_upper


inPath = '/Volumes/IAN USB/MINOS/Energy Calibration/C23_10.18.2018/'
energy = np.genfromtxt(inPath + 'Example_Energy_Pairs.json', delimiter=',')
energy = np.array(energy)
dbFiles = [x for x in os.listdir(inPath) if '.sqlite3' in x]
`

    cps_MUSE01 = np.sum(spectra_MUSE01,axis=1)
    dates_MUSE01 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_MUSE01])
    total_spectra_MUSE01 = np.sum(spectra_MUSE01,axis=0)

    dataTable = 'MUSE12_data'
    detDB.getColumn(dataTable, 'Spectrum__IntArray')
    spectra_MUSE12 = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data])
    detDB.getColumn(dataTable, 'Time')
    Time_MUSE12 = np.array(detDB.data)[:, 0]
    dateTime_MUSE12 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_MUSE12])
    detDB.getColumn(dataTable, 'Live_Time')
    Live_Time_MUSE12 = np.array(detDB.data)[:, 0]

    cps_MUSE12 = np.sum(spectra_MUSE12, axis=1)
    dates_MUSE12 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_MUSE12])
    total_spectra_MUSE12 = np.sum(spectra_MUSE12, axis=0)

    fig,ax = plt.subplots(nrows=2)
    ax[0].plot(cps_MUSE01,'b--',label='MUSE01');ax[0].grid(alpha=0.5)
    ax[1].plot(cps_MUSE12,'r--',label='MUSE12');ax[1].grid(alpha=0.5)

    fig,ax = plt.subplots(nrows=2)
    ax[0].plot(total_spectra_MUSE01,'k--',label='MUSE01')  ; ax[0].set_yscale('log')
    ax[1].plot(total_spectra_MUSE12,'r--',label='MUSE12') ; ax[1].set_yscale('log')

    plt.show()