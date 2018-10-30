import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from MINOS_analysis import database
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

'''
def getIsotopeEnergies_large(isotope):
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
        peak_windows_upper.append(170)
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
'''
def getIsotopeEnergies_large(isotope):
    '''
    Description: Provides energy and ROI for isotope.
    Param: isotope = string (e.g. 'Cs-137')
    '''
    peaks = [];
    peak_windows_lower = [];
    peak_windows_upper = []
    if isotope == "Cs-137":
        peaks.append(32) #
        peak_windows_lower.append(22)
        peak_windows_upper.append(42)
        peaks.append(450.0) # 662
        peak_windows_lower.append(400)
        peak_windows_upper.append(500)
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
    elif isotope == "Th-232":
        peaks.append(1725.0) # 2614
        peak_windows_lower.append(1650)
        peak_windows_upper.append(1830)
    elif isotope == "Eu-152":
        peaks.append(88)
        peak_windows_lower.append(74)
        peak_windows_upper.append(100)
        peaks.append(240) # 344.2785
        peak_windows_lower.append(210)
        peak_windows_upper.append(275)
        peaks.append(650) # 964.079
        peak_windows_lower.append(615)
        peak_windows_upper.append(698)
        peaks.append(760) # 1112.074
        peak_windows_lower.append(699)
        peak_windows_upper.append(820)
        peaks.append(940) #1408
        peak_windows_lower.append(880)
        peak_windows_upper.append(1050)
    return peaks, peak_windows_lower, peak_windows_upper

energy_conversion = {'88':121.7817,'240':344.2785,'650':940.079,'760':1112.074,'940':1408.006,\
                     '1725.0':2614,\
                     '32':32,'450.0':661.7}

inPath = '/Volumes/IAN USB/WIND/GADRAS dat Creation/'
#energy = np.genfromtxt(inPath + 'Example_Energy_Pairs.json', delimiter=',')
#energy = np.array(energy)
dbFiles = [x for x in os.listdir(inPath) if '.sqlite3' in x]
dbFiles = ['archerair2012-WIND_EnergyCal_Eu152_Cs137-2018-10-30T15.49.20.526.sqlite3']

savePath = inPath
for dbFile in dbFiles:
    print dbFile
    detDB = database(inPath + dbFile)

    dataTable = 'Det1_BR_4in_data'
    detDB.getColumn(dataTable, 'Spectrum__IntArray')
    spectra_Det1 = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data])
    detDB.getColumn(dataTable, 'Time')
    Time_Det1 = np.array(detDB.data)[:, 0]
    dateTime_Det1 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_Det1])
    detDB.getColumn(dataTable, 'Live_Time')
    Live_Time_Det1 = np.array(detDB.data)[:, 0]

    cps_Det1 = np.sum(spectra_Det1,axis=1)
    total_spectra_Det1 = np.sum(spectra_Det1,axis=0)

    dataTable = 'Det3_BL_4in_data'
    detDB.getColumn(dataTable, 'Spectrum__IntArray')
    spectra_Det3 = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data])
    detDB.getColumn(dataTable, 'Time')
    Time_Det3 = np.array(detDB.data)[:, 0]
    dateTime_Det3 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_Det3])
    detDB.getColumn(dataTable, 'Live_Time')
    Live_Time_Det3 = np.array(detDB.data)[:, 0]

    cps_Det3 = np.sum(spectra_Det3, axis=1)
    total_spectra_Det3 = np.sum(spectra_Det3, axis=0)

    dataTable = 'Det5_FC_2in_data'
    detDB.getColumn(dataTable, 'Spectrum__IntArray')
    spectra_Det5 = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data])
    detDB.getColumn(dataTable, 'Time')
    Time_Det5 = np.array(detDB.data)[:, 0]
    detDB.getColumn(dataTable, 'Live_Time')
    Live_Time_Det5 = np.array(detDB.data)[:, 0]

    cps_Det5 = np.sum(spectra_Det5, axis=1)
    total_spectra_Det5 = np.sum(spectra_Det5, axis=0)

    '''
    Energy Calibration Section
    '''

    source_names = ['Th-232', 'Cs-137', 'Eu-152']

    calibration_energies_pairs_Det1 = {}  # {energy_keV:channel,...}
    calibration_energies_pairs_Det3 = {}
    calibration_energies_pairs_Det5 = {}  # {energy_keV:channel,...}
    colors = ['g','r','k','y','m','c','orange','purple','grey','darkred','plum','tan']
    plotIt = True
    first = True
    count = 0
    for i in range(len(source_names)):
        name = source_names[i]
        peaks, peak_windows_lower, peak_windows_upper = getIsotopeEnergies_large(source_names[i])
        print '\t Fitting isotope:\t%s' % (source_names[i])
        for z in range(len(peaks)):

            j, k, l = (peaks[z], peak_windows_lower[z], peak_windows_upper[z])
            #indexs = np.where((energy > k) & (energy < l))
            indexs = np.where((np.arange(0,2048) > k) & (np.arange(0,2048) < l))
            '''
            fig, ax = plt.subplots(nrows=2)
            ax[0].plot(total_spectra_MUSE01, 'b',linewidth=0.75, label='MUSE01');
            ax[0].plot(indexs[0],total_spectra_MUSE01[indexs[0]],'k--',linewidth=1)
            ax[0].set_yscale('log');
            ax[0].legend()
            ax[1].plot(total_spectra_MUSE12, 'r',linewidth=0.5, label='MUSE12');
            ax[1].plot(indexs[0],total_spectra_MUSE12[indexs[0]],'k--',linewidth=1)
            ax[1].set_yscale('log');
            ax[1].legend()
            plt.show()
            '''
            spec = total_spectra_Det1[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
            channel_spec = np.arange(1,2049)[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
            a = spec[-1]
            b = np.max(spec) - spec[0]
            c = -1.0 * np.abs(spec[0] - spec[-1]) / (channel_spec[-1] - channel_spec[0])
            mu = channel_spec[np.argmax(spec)]
            sigma = 5
            p0 = [a, b, c, mu, sigma]
            popt, pcov = curve_fit(peakFunc, channel_spec, spec, p0=p0)
            popt[4] = np.abs(popt[4])
            channel = (np.abs(channel_spec - channel_spec[np.argmax(peakFunc(channel_spec, *popt))])).argmin()+channel_spec[0]
            j = energy_conversion[str(j)]
            if j not in calibration_energies_pairs_Det1.keys():
                calibration_energies_pairs_Det1[j] = [channel]
            else:
                calibration_energies_pairs_Det1[j].append(channel)


            if plotIt:
                '''
                Plots fitted data on top of spectra
                '''

                string = str(j) + ' ; ' +str(name)
                if first:
                    plt.figure(1)
                    first = False
                plt.plot(indexs[0], peakFunc(indexs[0], *popt),color=colors[count],linestyle='--', label=string, zorder=10)
            count+=1
    if plotIt:
        plt.figure(1)
        plt.plot(total_spectra_Det1, zorder=0)
        plt.yscale('log')
        plt.legend()
        plt.title('MUSE01')
        plt.grid(True, which='both',alpha=0.5)
        #plt.savefig(savePath+'MUSE01_calibration.png',dpi=600)
        plotDiagnostics = False

    first = True
    count = 0
    for i in range(len(source_names)):
        name = source_names[i]
        peaks, peak_windows_lower, peak_windows_upper = getIsotopeEnergies_large(source_names[i])
        print '\t Fitting isotope:\t%s' % (source_names[i])
        for z in range(len(peaks)):

            j, k, l = (peaks[z], peak_windows_lower[z], peak_windows_upper[z])
            #indexs = np.where((energy > k) & (energy < l))
            indexs = np.where((np.arange(0,2048) > k) & (np.arange(0,2048) < l))

            spec = total_spectra_Det3[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
            channel_spec = np.arange(1, 2049)[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
            a = spec[-1]
            b = np.max(spec) - spec[0]
            c = -1.0 * np.abs(spec[0] - spec[-1]) / (channel_spec[-1] - channel_spec[0])
            mu = channel_spec[np.argmax(spec)]
            sigma = 1.5
            p0 = [a, b, c, mu, sigma]
            popt, pcov = curve_fit(peakFunc, channel_spec, spec, p0=p0)
            popt[4] = np.abs(popt[4])
            channel = (np.abs(np.arange(channel_spec[0],channel_spec[-1]) - channel_spec[
                np.argmax(peakFunc(channel_spec, *popt))])).argmin()+channel_spec[0]
            j = energy_conversion[str(j)]
            if j not in calibration_energies_pairs_Det3.keys():
                calibration_energies_pairs_Det3[j] = [channel]
            else:
                calibration_energies_pairs_Det3[j].append(channel)

            if plotIt:
                '''
                Plots fitted data on top of spectra
                '''

                string = str(j) + ' ; ' + str(name)
                if first:
                    plt.figure(2)
                    first = False
                plt.plot(channel_spec, peakFunc(channel_spec, *popt), color=colors[count], linestyle='--', label=string, zorder=10)
            count += 1
    if plotIt:
        plt.figure(2)
        plt.plot(total_spectra_Det3, zorder=0)
        plt.yscale('log')
        plt.legend()
        plt.title('Det3')
        plt.grid(True, which='both',alpha=0.5)
        #plt.savefig(savePath+'MUSE12_calibration.png',dpi=600)

    first = True
    count = 0
    for i in range(len(source_names)):
        name = source_names[i]
        peaks, peak_windows_lower, peak_windows_upper = getIsotopeEnergies_large(source_names[i])
        print '\t Fitting isotope:\t%s' % (source_names[i])
        for z in range(len(peaks)):

            j, k, l = (peaks[z], peak_windows_lower[z], peak_windows_upper[z])
            #indexs = np.where((energy > k) & (energy < l))
            indexs = np.where((np.arange(0,2048) > k) & (np.arange(0,2048) < l))

            spec = total_spectra_Det5[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
            channel_spec = np.arange(1, 2049)[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
            a = spec[-1]
            b = np.max(spec) - spec[0]
            c = -1.0 * np.abs(spec[0] - spec[-1]) / (channel_spec[-1] - channel_spec[0])
            mu = channel_spec[np.argmax(spec)]
            sigma = 1.5
            p0 = [a, b, c, mu, sigma]
            popt, pcov = curve_fit(peakFunc, channel_spec, spec, p0=p0)
            popt[4] = np.abs(popt[4])
            channel = (np.abs(np.arange(channel_spec[0],channel_spec[-1]) - channel_spec[
                np.argmax(peakFunc(channel_spec, *popt))])).argmin()+channel_spec[0]
            j = energy_conversion[str(j)]
            if j not in calibration_energies_pairs_Det5.keys():
                calibration_energies_pairs_Det5[j] = [channel]
            else:
                calibration_energies_pairs_Det5[j].append(channel)

            if plotIt:
                '''
                Plots fitted data on top of spectra
                '''

                string = str(j) + ' ; ' + str(name)
                if first:
                    plt.figure(3)
                    first = False
                plt.plot(channel_spec, peakFunc(channel_spec, *popt), color=colors[count], linestyle='--', label=string, zorder=10)
            count += 1
    if plotIt:
        plt.figure(3)
        plt.plot(total_spectra_Det5, zorder=0)
        plt.yscale('log')
        plt.legend()
        plt.title('Det3')
        plt.grid(True, which='both',alpha=0.5)
        #plt.savefig(savePath+'Det5_calibration.png',dpi=600)

if False:
    '''
    Creating energy pairs
    '''
    energy_MUSE01 = [];  energy_MUSE01.append(0.0)
    channel_MUSE01 = []; channel_MUSE01.append(0)
    energy_MUSE12 = [];  energy_MUSE12.append(0.0)
    channel_MUSE12 = []; channel_MUSE12.append(0)

    for i in calibration_energies_pairs_MUSE01:
        energy_MUSE01.append(float(i))
        channel_MUSE01.append(calibration_energies_pairs_MUSE01[i][0])
    for i in calibration_energies_pairs_MUSE12:
        energy_MUSE12.append(float(i))
        channel_MUSE12.append(calibration_energies_pairs_MUSE12[i][0])

    energy_MUSE01 = sorted(energy_MUSE01)
    energy_MUSE12 = sorted(energy_MUSE12)
    channel_MUSE01 = sorted(channel_MUSE01)
    channel_MUSE12 = sorted(channel_MUSE12)

    # Manually creating Channels and energies for spline for specific MUSE node
    channels = np.array(channel_MUSE12); channels = np.delete(channels,5)
    energies = np.array(energy_MUSE12) ; energies = np.delete(energies,5)

    us = UnivariateSpline(channels[0:-1],energies[0:-1],k=2,s=0.0,ext=0)
    interp = us(np.arange(0,channels[-2]))
    ext = np.polyfit(np.array([channels[-2],channels[-1]]),np.array([energies[-2],energies[-1]]),1)
    extrap = np.polyval(ext, np.arange(channels[-2], 2048))
    cs = np.concatenate([interp,extrap])
    np.save(savePath + 'eCal_MUSE12_updated.npy',cs)
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
    plt.savefig(savePath+'Calibration_Pairs_MUSE12.png',dpi=600)

if (plotIt or plotDiagnostics):
    plt.show(True)