import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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

def quadratic(x,a,b,c):
    return a + b * x + c * x * x

inPath = '/Volumes/IAN USB/WIND/GADRAS dat Creation/energy_pairs/'
e_file = 'eCal_Det5.npy'

energy_pairs = np.load(inPath+e_file)
popt,popc = curve_fit(quadratic,range(0,len(energy_pairs)),energy_pairs)
values = quadratic(np.array(range(1,2049)),*popt)
print popt
plt.figure();
plt.plot(values,'r--')
plt.plot(energy_pairs,'b')
plt.grid()
plt.show()