import numpy as np
from scipy.interpolate import UnivariateSpline
import numba

@numba.jit(nopython=True)
def MBC_smoother(data):
    '''

    Parameters
    ----------
    data

    Returns
    -------

    '''
    x = np.arange(0, len(data[0]))
    smooth1 = np.array([]); smooth_ratios = np.array([])
    for y in data:
        value = UnivariateSpline(x, y)(x)
        smooth1 = np.concatenate((smooth1,np.array(value)))
    ratios = np.array(data) / smooth1
    np.nan_to_num(0.0)
    ratios = np.nan_to_num(ratios)
    for y in ratios:
        value = UnivariateSpline(x, y)(x)
        smooth_ratios = np.concatenate((smooth_ratios,value))
    smooth2 = smooth_ratios * smooth1
    for i in range(len(smooth2)):
        smooth2[smooth2 < 0] = 0.0
    return np.array(smooth1),np.array(smooth2)