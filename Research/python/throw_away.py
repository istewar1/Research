'''
Manually creating signal window for source moving through detector nodes
on Sept 28 2018

@ian_stewarts
'''

Node_spectra = {}
Node_counts = {}
import numpy as np
import matplotlib.pyplot as plt
nodes = ['MUSE01','MUSE04','MUSE06','MUSE10']
directory_path = '/Volumes/Ian External HD/Node Data/'
for i in nodes:
    print i
    Node_spectra[i] = np.load(directory_path+'npy_arrays/'+i+'_spectra.npy')
    Node_counts[i] = np.load(directory_path+'npy_arrays/'+i+'_counts.npy')

energy = 3*np.arange(0,1000)
Alarm_values = {'MUSE01':[190389,190394], \
                'MUSE04': [195977,195985], \
                'MUSE06': [201422,201433], \
                'MUSE10': [194634,194641]}
node_colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']

fig,ax = plt.subplots()
count=0
for i in Alarm_values:
    indexes = Alarm_values[i]
    spectrum = np.sum(Node_spectra[i][indexes[0]:indexes[1]],axis=0)/float(indexes[1]-indexes[0])
    print Node_counts[i][indexes[0]:indexes[1]]
    ax.plot(energy,spectrum,label=str(i),c=node_colors[count])
    count+=1
ax.grid(alpha=0.5, linewidth=0.5,which='both')
ax.set_yscale('linear')
ax.set_axisbelow(True)
ax.legend(loc='upper right',fancybox=False,shadow=True)
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Normalized Counts per bin')


plt.show(True)