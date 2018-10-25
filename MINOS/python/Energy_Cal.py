import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as md
import os
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
from scipy.interpolate import UnivariateSpline
import sqlite3


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
dbFiles = ['archerair2012-Characterization_C23_Cs137_Eu152-2018-10-18T15.58.32.839.sqlite3']

for dbFile in dbFiles:
    print dbFile
    detDB = database(inPath + dbFile)

    dataTable = 'MUSE01_data'
    detDB.getColumn(dataTable, 'Spectrum__IntArray')
    spectra_MUSE01 = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data])
    detDB.getColumn(dataTable, 'Time')
    Time_MUSE01 = np.array(detDB.data)[:, 0]
    dateTime_MUSE01 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_MUSE01])
    detDB.getColumn(dataTable, 'Live_Time')
    Live_Time_MUSE01 = np.array(detDB.data)[:, 0]

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

    '''
    Energy Calibration Section
    '''

    source_names = ['Th232', 'Cs-137', 'Eu-152']

    calibration_energies_pairs_MUSE01 = {}  # {energy_keV:channel,...}
    calibration_energies_pairs_MUSE12 = {}
    colors = ['g','r','k','y','m','c','orange']
    plotIt = True
    first = True
    count = 0
    for i in range(len(source_names)):
        name = source_names[i]
        peaks, peak_windows_lower, peak_windows_upper = getIsotopeEnergies_large(source_names[i])
        print '\t Fitting isotope:\t%s' % (source_names[i])
        for z in range(len(peaks)):

            j, k, l = (peaks[z], peak_windows_lower[z], peak_windows_upper[z])
            indexs = np.where((energy > k) & (energy < l))

            spec = total_spectra_MUSE01[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
            e = np.arange(1,2049)[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
            a = spec[-1]
            b = np.max(spec) - spec[0]
            c = -1.0 * np.abs(spec[0] - spec[-1]) / (energy[-1] - energy[0])
            mu = e[np.argmax(spec)]
            sigma = 0.9
            p0 = [a, b, c, mu, sigma]
            popt, pcov = curve_fit(peakFunc, e, spec, p0=p0)
            popt[4] = np.abs(popt[4])
            channel = (np.abs(e - e[np.argmax(peakFunc(e, *popt))])).argmin()
            if j not in calibration_energies_pairs_MUSE01.keys():
                calibration_energies_pairs_MUSE01[j] = [channel]
            else:
                calibration_energies_pairs_MUSE01[j].append(channel)

            spec = total_spectra_MUSE12[indexs[0][0] / 2:indexs[0][len(indexs[0]) - 1] / 2]
            channel_spec = np.arange(indexs[0][0] / 2, indexs[0][len(indexs[0]) - 1] / 2)
            a = spec[-1]
            b = np.max(spec) - spec[0]
            c = -1.0 * np.abs(spec[0] - spec[-1]) / (channel_spec[-1] - channel_spec[0])
            mu = channel_spec[np.argmax(spec)]
            sigma = 1.5
            p0 = [a, b, c, mu, sigma]
            popt, pcov = curve_fit(peakFunc, channel_spec, spec, p0=p0)
            popt[4] = np.abs(popt[4])
            channel = (np.abs(np.arange(0, len(total_spectra_MUSE12)) - channel_spec[
                np.argmax(peakFunc(channel_spec, *popt))])).argmin()
            if j not in calibration_energies_pairs_MUSE12.keys():
                calibration_energies_pairs_MUSE12[j] = [channel]
            else:
                calibration_energies_pairs_MUSE12[j].append(channel)

            if plotIt:
                '''
                Plots fitted data on top of spectra
                '''

                string = str(j) + ' ; ' +str(name)
                if first:
                    plt.figure()
                    first = False
                plt.plot(e, peakFunc(e, *popt),color=colors[count],linestyle='--', label=string, zorder=10)
            count+=1
    if plotIt:
        plt.plot(total_spectra_MUSE01, zorder=0)
        plt.yscale('log')
        plt.legend()
        plt.grid(True, which='both')
        plt.title('MUSE01')

        plotDiagnostics = False

    for i in range(len(source_names)):
        spec = total_spectra_MUSE12[indexs[0][0] / 2:indexs[0][len(indexs[0]) - 1] / 2]
        channel_spec = np.arange(indexs[0][0] / 2, indexs[0][len(indexs[0]) - 1] / 2)
        a = spec[-1]
        b = np.max(spec) - spec[0]
        c = -1.0 * np.abs(spec[0] - spec[-1]) / (channel_spec[-1] - channel_spec[0])
        mu = channel_spec[np.argmax(spec)]
        sigma = 1.5
        p0 = [a, b, c, mu, sigma]
        popt, pcov = curve_fit(peakFunc, channel_spec, spec, p0=p0)
        popt[4] = np.abs(popt[4])
        channel = (np.abs(np.arange(0, len(total_spectra_MUSE12)) - channel_spec[
            np.argmax(peakFunc(channel_spec, *popt))])).argmin()
        if j not in calibration_energies_pairs_MUSE12.keys():
            calibration_energies_pairs_MUSE12[j] = [channel]
        else:
            calibration_energies_pairs_MUSE12[j].append(channel)

        if plotIt:
            '''
            Plots fitted data on top of spectra
            '''

            string = str(j) + ' ; ' + str(name)
            if first:
                plt.figure()
                first = False
            plt.plot(e, peakFunc(e, *popt), color=colors[count], linestyle='--', label=string, zorder=10)
        count += 1
    if plotIt:
        plt.plot(total_spectra_MUSE12, zorder=0)
        plt.yscale('log')
        plt.legend()
        plt.title('MUSE12')
        plt.grid(True, which='both')

    '''
    Creating energy pairs
    '''
    energy_values_MUSE01 = [];
    energy_channel_MUSE01 = []
    for i in calibration_energies_pairs_MUSE01:
        if len(calibration_energies_pairs_MUSE01[i]) == 1:
            energy_values_MUSE01.append(i)
            energy_channel_MUSE01.append(int(calibration_energies_pairs_MUSE01[i][0]))
        else:
            energy_values_MUSE01.append(i)
            energy_channel_MUSE01.append(int(np.mean(calibration_energies_pairs_MUSE01[i])))
    energy_pairs_MUSE01 = generate_calibration_axis(energy_channel_MUSE01, energy_values_MUSE01)
    '''
    energy_values_MUSE12 = [];
    energy_channel_MUSE12 = []
    for i in calibration_energies_pairs_MUSE12:
        if len(calibration_energies_pairs_MUSE12[i]) == 1:
            energy_values_MUSE12.append(i)
            energy_channel_MUSE12.append(int(calibration_energies_pairs_MUSE12[i][0]))
        else:
            energy_values_MUSE12.append(i)
            energy_channel_MUSE12.append(int(np.mean(calibration_energies_pairs_MUSE12[i])))
    energy_pairs_MUSE12 = generate_calibration_axis(energy_channel_MUSE12, energy_values_MUSE12)
    '''
    '''
    us = UnivariateSpline(channels[0:-1],energies[0:-1],k=2,s=0.0,ext=0)
    interp = us(np.arange(0,channels[-2]))
    ext = np.polyfit(np.array([channels[-2],channels[-1]]),np.array([energies[-2],energies[-1]]),1)
    extrap = np.polyval(ext, np.arange(channels[-2], 2048))
    cs = np.concatenate([interp,extrap])
    x_range = np.arange(len(cs))
    fig,axes = plt.subplots(2,1)
    axes[0].plot(x_range, np.array(cs), label='Cubic Spline')
    axes[0].plot(channels,np.array(energies),'r*',label='Fitted Energies')
    axes[0].grid(alpha=0.5)
    axes[0].set_xlabel('Channel')
    axes[0].set_ylabel('Energy')
    axes[1].plot(x_range, np.array((x_range))/np.array(cs), label='Cubic Spline')
    axes[1].plot(channels,np.array(channels)/np.array(energies),'r*',label='Fitted Energies')
    axes[1].grid(alpha=0.5)
    axes[1].set_xlabel('Channel')
    axes[1].set_ylabel('Channel/Energy')
    fig.tight_layout()  
    '''

if (plotIt or plotDiagnostics):
    plt.show()