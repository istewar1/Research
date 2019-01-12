import numpy as np
import h5py
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt

filename = '/Volumes/Ian External HD/RSL NOVArray/data/Aug1-Aug8_2018/Full Sensor Reports.2018-08-01 04_00_00+0000 - 2018-08-08 04_00_00+0000.4 rsl Sensors (4).hdf5'
hdf5 = h5py.File(filename, 'r')
MUSE01 = np.array(hdf5['digiBASE-RH 17178701/RadiationReading/spectrum/'])[30320:40000]
filename = '/Volumes/Ian External HD/Node Data/hdf5files/MUSE01/MUSE01-2018-08-01T00.35.06.865.hdf5'
MUSE01 = h5py.File(filename,'r')['2x4x16Spectra'][0:10000]
svd = TruncatedSVD(n_components=999)
SVD = svd.fit(MUSE01)
U = SVD.fit_transform(MUSE01)
Sigma = SVD.explained_variance_ratio_
VT = SVD.components_
var_ex = [];
first = True
for i in range(len(Sigma)):
    j = float(Sigma[i])*100.
    if first:
        var_ex.append(j)
        first = False
    else:
        var_ex.append(var_ex[i - 1] + j)
if True:
    fig, ax1 = plt.subplots()
    ax1.bar(np.arange(1, len(Sigma) + 1), np.array(Sigma)*100., label='Explained Variance')
    ax1.set_xlabel('Component No.')
    ax1.set_ylabel('Explained Variance')
    ax1.set_yscale('log')
    ax1.yaxis.label.set_color('blue')
    ax1.tick_params(axis='y', colors='blue')
    ax1.set_ylim(3E-4,max(Sigma)*100.)
    ax2 = ax1.twinx()
    ax2.plot(var_ex, 'r--', label='Cumulative Variance')
    ax2.set_ylabel('Cumulative Variance')
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
    # customizing legend box
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='center right', fancybox=False, shadow=True)
    #
    ax1.grid(alpha=0.5)
    ax1.set_axisbelow(True)
print np.where(np.array(var_ex)>95.)[0][0]
plt.show(True)
