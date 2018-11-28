
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
        peak_windows_lower.append(1680)
        peak_windows_upper.append(1880)
    elif isotope == "Eu-152":
        peaks.append(88)
        peak_windows_lower.append(74)
        peak_windows_upper.append(100)
        peaks.append(240) # 344.2785
        peak_windows_lower.append(210)
        peak_windows_upper.append(275)
        peaks.append(650) # 964.079
        peak_windows_lower.append(612)
        peak_windows_upper.append(705)
        peaks.append(760) # 1112.074
        peak_windows_lower.append(702)
        peak_windows_upper.append(820)
        peaks.append(940) #1408
        peak_windows_lower.append(880)
        peak_windows_upper.append(1050)
    return peaks, peak_windows_lower, peak_windows_upper

def getIsotopeEnergies_small(isotope):
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
        peaks.append(811.0) # 2614
        peak_windows_lower.append(750)
        peak_windows_upper.append(850)
    elif isotope == "Eu-152":
        peaks.append(88)
        peak_windows_lower.append(74)
        peak_windows_upper.append(100)
        peaks.append(240) # 344.2785
        peak_windows_lower.append(210)
        peak_windows_upper.append(275)
        peaks.append(650) # 964.079
        peak_windows_lower.append(612)
        peak_windows_upper.append(705)
        peaks.append(760) # 1112.074
        peak_windows_lower.append(702)
        peak_windows_upper.append(820)
        peaks.append(940) #1408
        peak_windows_lower.append(880)
        peak_windows_upper.append(1050)
    return peaks, peak_windows_lower, peak_windows_upper

energy_conversion = {'88':121.7817,'240':344.2785,'650':940.079,'760':1112.074,'940':1408.006,\
                     '1725.0':2614,\
                     '32':32,'450.0':661.7}

inPath = '/Users/i6o/Downloads/'
import h5py
dbFiles = ['MUSE01-2018-08-08T08.47.16.170.hdf5']
savePath = inPath
spectra = []
times = []
for dbFile in dbFiles:

    hdf5 = h5py.File(inPath + dbFile,'r')
    for j in hdf5.keys():
        if 'MUSE01_Spectra' in j:
            spectra.append(np.array(hdf5.get(str(j))))
        if 'MUSE01_Times' in j:
            times.append(np.array(hdf5.get(str(j))))
    background_spectrum = np.sum(spectra[0][10000:50000],axis=0)
    source_spectrum = np.sum(spectra[0][72000:79693],axis=0)
    plt.figure()
    plt.plot(source_spectrum,label='source');plt.yscale('log')
    plt.figure()
    plt.plot(background_spectrum,label='background');plt.yscale('log')
    plt.show()
    '''
    Energy Calibration Section
    '''

    source_names = ['Th-232', 'Eu-152']

    calibration_energies_pairs_Det1 = {}  # {energy_keV:channel,...}
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
            indexs = np.where((np.arange(0,1024) > k) & (np.arange(0,1024) < l))
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
        plt.title('Det1')
        plt.grid(True, which='both',alpha=0.5)
        #plt.savefig(savePath+'Det1_calibration.png',dpi=600)
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
            sigma = 5
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
        #plt.savefig(savePath+'Det3_calibration.png',dpi=600)

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
            sigma = 5
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
        plt.title('Det5')
        plt.grid(True, which='both',alpha=0.5)
        #plt.savefig(savePath+'Det5_calibration.png',dpi=600)

if False:
    '''
    Creating energy pairs
    '''
    energy_Det1 = [0.0];  #energy_Det1.append(0.0)
    channel_Det1 = [0]; #channel_Det1.append(0)
    energy_Det3 = [0.0];  #energy_Det3.append(0.0)
    channel_Det3 = [0]; #channel_Det3.append(0)
    energy_Det5 = [0.0];  #energy_Det5.append(0.0)
    channel_Det5 = [0]; #channel_Det5.append(0)

    for i in calibration_energies_pairs_Det1:
        energy_Det1.append(float(i))
        channel_Det1.append(calibration_energies_pairs_Det1[i][0])
    for i in calibration_energies_pairs_Det3:
        energy_Det3.append(float(i))
        channel_Det3.append(calibration_energies_pairs_Det3[i][0])
    for i in calibration_energies_pairs_Det5:
        energy_Det5.append(float(i))
        channel_Det5.append(calibration_energies_pairs_Det5[i][0])

    energy_Det1 = sorted(energy_Det1)
    energy_Det3 = sorted(energy_Det3)
    channel_Det1 = sorted(channel_Det1)
    channel_Det3 = sorted(channel_Det3)
    energy_Det5 = sorted(energy_Det5)
    channel_Det5 = sorted(channel_Det5)

    # Manually creating Channels and energies for spline for specific MUSE node
    channels = np.array(channel_Det5); channels = np.delete(channels,5)
    energies = np.array(energy_Det5) ; energies = np.delete(energies,5)
    print len(channels),len(energies)
    us = UnivariateSpline(channels[0:-1],energies[0:-1],k=2,s=0.0,ext=0)
    interp = us(np.arange(0,channels[-2]))
    ext = np.polyfit(np.array([channels[-2],channels[-1]]),np.array([energies[-2],energies[-1]]),1)
    extrap = np.polyval(ext, np.arange(channels[-2], 2048))
    cs = np.concatenate([interp,extrap])
    np.save(savePath + 'eCal_Det5.npy',cs)
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
    plt.savefig(savePath+'Calibration_Pairs_Det5.png',dpi=600)

if (plotIt or plotDiagnostics):
    plt.show(False)