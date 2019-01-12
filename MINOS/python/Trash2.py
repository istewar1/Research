import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import curve_fit
import h5py; import matplotlib.pyplot as plt

def quadratic(x,a,b,c):
    return a + b * x + c * x * x

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
        peak_windows_upper.append(520)
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
        #peaks.append(650) # 964.079
        #peak_windows_lower.append(612)
        #peak_windows_upper.append(705)
        peaks.append(760) # 1112.074
        peak_windows_lower.append(702)
        peak_windows_upper.append(820)
        peaks.append(940) #1408
        peak_windows_lower.append(880)
        peak_windows_upper.append(1050)
    return peaks, peak_windows_lower, peak_windows_upper

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

inPath = '/Volumes/IAN USB/WIND/GADRAS dat Creation/2018-11-21/'
filename = 'Eu152_Cs137_No_5770-2018-12-21T14.29.35.005.hdf5'
hdf5 = h5py.File(inPath+filename,'r')
spectra = np.sum(hdf5['WIND7_Spectra'],axis=0)

'''
Energy Calibration Section
'''
energy_conversion = {'88':121.7817,'240':344.2785,'650':940.079,'760':1112.074,'940':1408.006,\
                     '1725.0':2614,\
                     '32':32,'450.0':661.7}
source_names = ['Th-232', 'Cs-137', 'Eu-152']

calibration_energies_pairs_Det1 = {}  # {energy_keV:channel,...}
calibration_energies_pairs_Det3 = {}
calibration_energies_pairs_Det5 = {}  # {energy_keV:channel,...}
colors = ['g', 'r', 'k', 'y', 'm', 'c', 'orange', 'purple', 'grey', 'darkred', 'plum', 'tan']
plotIt = True
first = True
count = 0
for i in range(len(source_names)):
    name = source_names[i]
    peaks, peak_windows_lower, peak_windows_upper = getIsotopeEnergies_large(source_names[i])
    print '\t Fitting isotope:\t%s' % (source_names[i])
    for z in range(len(peaks)):

        j, k, l = (peaks[z], peak_windows_lower[z], peak_windows_upper[z])
        # indexs = np.where((energy > k) & (energy < l))
        indexs = np.where((np.arange(0, 2048) > k) & (np.arange(0, 2048) < l))
        spec = spectra[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        channel_spec = np.arange(1, 2049)[indexs[0][0]:indexs[0][len(indexs[0]) - 1]]
        a = spec[-1]
        b = np.max(spec) - spec[0]
        c = -1.0 * np.abs(spec[0] - spec[-1]) / (channel_spec[-1] - channel_spec[0])
        mu = channel_spec[np.argmax(spec)]
        sigma = 5
        p0 = [a, b, c, mu, sigma]
        popt, pcov = curve_fit(peakFunc, channel_spec, spec, p0=p0)
        popt[4] = np.abs(popt[4])
        channel = (np.abs(channel_spec - channel_spec[np.argmax(peakFunc(channel_spec, *popt))])).argmin() + \
                  channel_spec[0]
        j = energy_conversion[str(j)]
        if j not in calibration_energies_pairs_Det1.keys():
            calibration_energies_pairs_Det1[j] = [channel]
        else:
            calibration_energies_pairs_Det1[j].append(channel)

        if plotIt:
            '''
            Plots fitted data on top of spectra
            '''

            string = str(j) + ' ; ' + str(name)
            if first:
                plt.figure(1)
                first = False
            plt.plot(indexs[0], peakFunc(indexs[0], *popt), color=colors[count], linestyle='--', label=string,
                     zorder=10)
        count += 1
print calibration_energies_pairs_Det1
if plotIt:
    plt.figure(1)
    plt.plot(spectra, zorder=0)
    plt.yscale('log')
    plt.legend()
    plt.title('Det1')
    plt.grid(True, which='both', alpha=0.5)
    plt.show(False)
bins = []; energies = []
for j in calibration_energies_pairs_Det1:
    bins.append(calibration_energies_pairs_Det1[j][0])
    energies.append(j)
bins.append(0)
energies.append(0.0)
bins, energies = zip(*sorted(zip(bins, energies)))

spline_fit = interpolate.InterpolatedUnivariateSpline(x=bins,y=energies,k=2)
ext = np.polyfit([bins[-2], bins[-1]], [energies[-2], energies[-1]], 1)
interpolated_energy_data = spline_fit(np.arange(bins[0], bins[-1]))
extrapolated_energy_data = np.polyval(ext, np.arange(bins[-1], 2048))
energy = np.concatenate([interpolated_energy_data, extrapolated_energy_data])
popt, popc = curve_fit(quadratic, range(0, len(energy)), energy)
values = quadratic(np.arange(0, len(energy)), *popt)
plt.figure();plt.plot(energy,np.arange(0,len(energy))/energy,'--');plt.plot(energies,np.array(bins)/np.array(energies),'r*');
plt.figure();plt.plot(energy,spectra);plt.yscale('log')
plt.show()
print popt