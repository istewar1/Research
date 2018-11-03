# -*- coding: utf-8 -*-
'''---------------------------------------------------------
#Description:    Plotting HPGe spectrum (.spe) files within a
                    defined root directory.
---------------------------------------------------------'''
import matplotlib.pyplot as plt
import numpy as np
import struct
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.integrate import quad
import os

# plt.close('all')
'''---------------------------------------------------------
Define Functions
---------------------------------------------------------'''


def readSPE(dataFile):
    '''
    Reads spectrum (.spe) file and obtains values:
        date     : date of measurement
        liveTime : live time of measurement
        spec     : measurement spectrum
        cals     : energy calibration parameters
        e.g. spectrum,liveTime,a,b,c=readSPE(rootDir+dataFile)
            - This exports variables from dataFile in rootDir

    '''
    f = open(dataFile, 'r')
    data = f.read().split('\r\n')
    # get date
    index = data.index('$DATE_MEA:')
    date = data[index + 1]
    # get time
    index = data.index('$MEAS_TIM:')
    liveTime = float(data[index + 1].split(' ')[-1])
    # number of channels and spectrum
    index = data.index('$DATA:')
    numChans = int(data[index + 1].split(' ')[-1])
    spec = data[index + 2:index + 2 + numChans]
    spec = [int(x.replace(' ', '')) for x in spec]
    # calibration factors : a+bx+cx^2
    index = data.index('$MCA_CAL:')
    cals = [float(x) for x in data[index + 2].split(' ')[:-1]]
    a = cals[0]
    b = cals[1]
    c = cals[2]
    return spec, liveTime, a, b, c
    f.close()


def energyPerChannel(x, a, b, c=-1, d=-1):
    '''
    Calibrates spectrum to Energy per channel based on
    HPGe spectrum:
        x = channels (e.g. 0, 1, ..., 2048)
        a, b, c = fitting parameters provided in *.SPE
    '''
    if c == -1 and d == -1:
        return a + b * x
    elif c != -1 and d == -1:
        return a + b * x + c * x * x
    elif c != -1 and d != -1:
        return a + b * x + c * x * x + d * x * x * x


def gaussian(x, b, mu, sigma):
    '''
    Function to fit ROI peak with Gaussian-fit only
    '''
    return b * np.exp(-(x - mu) ** 2 / (2.0 * sigma ** 2))


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


'''----------------------------------------------------------
Main
----------------------------------------------------------'''
if __name__ == '__main__':
    # rootDir = '/Users/i6o/MINOS/REDC-HFIR_HPGe/7Nov2017/'
    rootDir = '/Users/i6o/Documents/REDC-HFIR-Characterization/grass/'
    dataFiles = os.listdir(rootDir)
    dataFiles = np.array([x for x in dataFiles if '.Spe' in x])

    plotSpec = False;  # plot total spectrum
    plotROI = [False, False, False, False]  # plot ROI spectrum (Pb-212, K-40, Tl-208)
    plotComp = False;  # plot intensities of each RIOs together
    plotKandTh = False;  # plot K-40 versus Thorium (Tl-208)
    pltKandCs = False

    count = 1

    Isotopes = ['Pb-212', 'K-40', 'Tl-208', 'Cs-137']
    E_min = [220, 1400, 2500, 655]
    E_max = [260, 1500, 2700, 667]
    Pb212ints = [];
    K40ints = [];
    Tl208ints = [];
    Cs137ints = []

    for dataFile in dataFiles:
        spectrum = [];
        energies = [];
        spectrum, liveTime, a, b, c = readSPE(rootDir + dataFile)
        energies = energyPerChannel(np.arange(len(spectrum)), a, b, c)

        # average count rate for spectra collected
        if count == 1:
            count_rate = [a / liveTime for a in spectrum]
        if count > 1:
            current_count_rate = [a / liveTime for a in spectrum]
            count_rate = [(a + b) / 2 for a, b in zip(current_count_rate, count_rate)]

        if plotSpec:
            plt.figure(count)
            plt.title(dataFile)
            plt.step(energies, spectrum)
            plt.xlabel('Energy (keV)')
            plt.ylabel(r'Counts')
            plt.grid(True, which="both")
            plt.show()

        # ------- Pb-212 ------- #
        eMin = E_min[0]
        eMax = E_max[0]
        indexs = np.where((energies > eMin) & (energies < eMax))
        spec = spectrum[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        energy = energies[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        # assigning linear line and gaussian values
        a = spec[-1]
        b = np.max(spec) - spec[0]
        c = -1.0 * np.abs(spec[0] - spec[-1]) / (eMax - eMin)
        mu = energy[np.argmax(spec)]
        sigma = 0.9
        p0 = [a, b, c, mu, sigma]
        popt, pcov = curve_fit(peakFunc, energy, spec, p0=p0)
        popt[4] = np.abs(popt[4])
        deltaX = np.average([energy[m] - energy[m - 1] for m in range(1, len(energy))])
        # integrate gaussian function from eMin to eMax for gaussian values in popt
        integral = quad(gaussian, eMin, eMax, args=(popt[1], popt[3], popt[4]))[0] / deltaX
        # Pb-212 intensity and standard deviation (#/sec)
        Pb212ints.append([integral / liveTime, np.sqrt(integral) / liveTime])
        Pb212_counts = [i[0] for i in Pb212ints]
        Pb212_er = [i[1] for i in Pb212ints]

        if plotROI[0]:
            plt.figure(count)
            # plotting raw data
            plt.step(energy, spec)
            # plotting fit function to data
            plt.plot(energy, peakFunc(energy, *popt), marker='', label=r'$^{212}$Pb: %.0f cnts' % (integral))
            plt.grid(True);
            plt.show()
            plt.title(dataFile + "\nIsotope:" + Isotopes[0])
            plt.xlabel('Energy (keV)')
            plt.ylabel(r'Counts')
            plt.grid(True, which="both")

            # ------- K-40 ------- #
        eMin = E_min[1]
        eMax = E_max[1]
        indexs = np.where((energies > eMin) & (energies < eMax))
        spec = spectrum[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        energy = energies[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        # assigning linear line and gaussian values
        a = spec[-1]
        b = np.max(spec) - spec[0]
        c = -1.0 * np.abs(spec[0] - spec[-1]) / (eMax - eMin)
        mu = energy[np.argmax(spec)]
        sigma = 0.9
        p0 = [a, b, c, mu, sigma]
        popt, pcov = curve_fit(peakFunc, energy, spec, p0=p0)
        popt[4] = np.abs(popt[4])
        deltaX = np.average([energy[m] - energy[m - 1] for m in range(1, len(energy))])
        # integrate gaussian function from eMin to eMax for gaussian values in popt
        # to find the intensity (counts) under peak
        integral = quad(gaussian, eMin, eMax, args=(popt[1], popt[3], popt[4]))[0] / deltaX
        # K-40 intensity and standard deviation (#/sec)
        K40ints.append([integral / liveTime, np.sqrt(integral) / liveTime])
        K40_counts = [i[0] for i in K40ints]
        K40_er = [i[1] for i in K40ints]

        if plotROI[1]:
            plt.figure(count)
            # plotting raw data
            plt.step(energy, spec)
            # plotting fit function to data
            plt.plot(energy, peakFunc(energy, *popt), marker='', label=r'$^{40}$K: %.0f cnts' % (integral))
            plt.grid(True);
            plt.show()
            plt.title(dataFile + "\nIsotope:" + Isotopes[0])
            plt.xlabel('Energy (keV)')
            plt.ylabel(r'Counts')
            plt.grid(True, which="both")

            # ------- Tl-208 ------- #
        eMin = E_min[2]
        eMax = E_max[2]
        indexs = np.where((energies > eMin) & (energies < eMax))
        spec = spectrum[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        energy = energies[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        # assigning linear line and gaussian values
        a = spec[-1]
        b = np.max(spec) - spec[0]
        c = -1.0 * np.abs(spec[0] - spec[-1]) / (eMax - eMin)
        mu = energy[np.argmax(spec)]
        sigma = 0.9
        p0 = [a, b, c, mu, sigma]
        popt, pcov = curve_fit(peakFunc, energy, spec, p0=p0)
        popt[4] = np.abs(popt[4])
        deltaX = np.average([energy[m] - energy[m - 1] for m in range(1, len(energy))])
        # integrate gaussian function from eMin to eMax for gaussian values in popt
        integral = quad(gaussian, eMin, eMax, args=(popt[1], popt[3], popt[4]))[0] / deltaX
        # Tl-208 intensity and standard deviation (#/sec)
        Tl208ints.append([integral / liveTime, np.sqrt(integral) / liveTime])
        Tl208_counts = [i[0] for i in Tl208ints]
        Tl208_er = [i[1] for i in Tl208ints]

        if plotROI[2]:
            plt.figure(count)
            # plotting raw data
            plt.step(energy, spec)
            # plotting fit function to data
            plt.plot(energy, peakFunc(energy, *popt), marker='', label=r'$^{208}$Tl: %.0f cnts' % (integral))
            plt.grid(True);
            plt.show()
            plt.title(dataFile + "\nIsotope:" + Isotopes[0])
            plt.xlabel('Energy (keV)')
            plt.ylabel(r'Counts')
            plt.grid(True, which="both")

            # ------- Cs-137 ------- #
        eMin = E_min[3]
        eMax = E_max[3]
        indexs = np.where((energies > eMin) & (energies < eMax))
        spec = spectrum[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        energy = energies[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        # assigning linear line and gaussian values
        a = spec[-1]
        b = np.max(spec) - spec[0]
        c = -1.0 * np.abs(spec[0] - spec[-1]) / (eMax - eMin)
        mu = energy[np.argmax(spec)]
        sigma = 0.9
        p0 = [a, b, c, mu, sigma]
        popt, pcov = curve_fit(peakFunc, energy, spec, p0=p0)
        popt[4] = np.abs(popt[4])
        deltaX = np.average([energy[m] - energy[m - 1] for m in range(1, len(energy))])
        # integrate gaussian function from eMin to eMax for gaussian values in popt
        integral = quad(gaussian, eMin, eMax, args=(popt[1], popt[3], popt[4]))[0] / deltaX
        # Cs-137 intensity and standard deviation (#/sec)
        Cs137ints.append([integral / liveTime, np.sqrt(integral) / liveTime])
        Cs137_counts = [i[0] for i in Cs137ints]
        Cs137_er = [i[1] for i in Cs137ints]

        if plotROI[3]:
            plt.figure(count)
            # plotting raw data
            plt.step(energy, spec)
            # plotting fit function to data
            plt.plot(energy, peakFunc(energy, *popt), marker='', label=r'$^{137}$Cs: %.0f cnts' % (integral))
            plt.grid(True);
            plt.show()
            plt.title(dataFile + "\nIsotope:" + Isotopes[0])
            plt.xlabel('Energy (keV)')
            plt.ylabel(r'Counts')
            plt.grid(True, which="both")

        count += 1
        print "\tRead in spectrum file:\t" + dataFile + "\t: COMPLETE"

    if plotComp:
        length = len(Pb212_counts)
        plt.figure()
        plt.scatter(np.full([1, length], 1), Pb212_counts, label='Pb-212 FTIG', color="b")
        plt.scatter(np.full([1, length], 2), K40_counts, label='K-40 FTIG', c="b")
        plt.scatter(np.full([1, length], 3), Tl208_counts, label='Tl-208 FTIG', c="b")
        plt.xticks([1, 2, 3], Isotopes);
        plt.grid(True);
        plt.title('Intensity of ROIs');
        plt.ylabel('Intensity (#/sec)');

    if plotKandTh:
        plt.figure(7)
        location = 'REDC Soil'
        plt.errorbar(K40_counts, Tl208_counts, xerr=K40_er, yerr=Tl208_er, ls='', marker='o', capsize=2, label=location)
        plt.xlabel('K-40 Intensity (#/sec)');
        plt.ylabel('Tl-208 Intensity (#/sec)');
        plt.title('K-40 vs Tl-208 Comparison');
        plt.grid(True);
        plt.show()
    if pltKandCs:
        plt.figure(1)
        location = 'FTIG Soil'
        plt.errorbar(K40_counts, Cs137_counts, xerr=K40_er, yerr=Cs137_er, ls='', marker='o', capsize=2, label=location)
        plt.xlabel('K-40 Intensity (#/sec)');
        plt.ylabel('Cs-137 Intensity (#/sec)');
        plt.title('K-40 vs Cs-137 Comparison');
        plt.grid(True);
        plt.show()
    plt.show()