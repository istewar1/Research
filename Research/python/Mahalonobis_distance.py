def rebin_Archer(x1, y1, x2):
    """
    Rebin histogram values y1 from old bin edges x1 to new edges x2.

    Input
    -----
     * x1 : m+1 array of old bin edges.
     * y1 : m array of old histogram values. This is the total number in
              each bin, not an average.
     * x2 : n+1 array of new bin edges.

    Returns
    -------
     * y2 : n array of rebinned histogram values.

    The rebinning algorithm assumes that the counts in each old bin are
    uniformly distributed in that bin.

    Bins in x2 that are entirely outside the range of x1 are assigned 0.
    """

    x1 = np.asarray(x1)
    y1 = np.asarray(y1)
    x2 = np.asarray(x2)

    # allocating y2 vector
    m = y1.size
    n = x2.size - 1
    y2 = np.zeros(n, dtype=y1.dtype)

    i_place = np.searchsorted(x2, x1)

    # find out where x2 intersects with x1, this will determine which x2 bins
    # we need to consider
    start_pos = 0
    end_pos = n

    start_pos_test = np.where(i_place == 0)[0]
    if start_pos_test.size > 0:
        start_pos = start_pos_test[-1]

    end_pos_test = np.where(i_place == x1.size)[0]
    if end_pos_test.size > 0:
        end_pos = end_pos_test[0]

    # the first bin totally covers x1 range
    if (start_pos == end_pos - 1
            and i_place[start_pos] == 0
            and i_place[start_pos + 1] == x1.size):
        y2[start_pos] = y1[start_pos]
        return y2

    # print len(x1),len(y1),len(y2)
    # print len(i_place)
    # print i_place[980:len(i_place)]

    # first(0th) X1 bin
    ibin = 0
    i_place_lower = i_place[ibin]
    i_place_upper = i_place[ibin + 1]
    if i_place_lower == 0 and i_place_upper > 0:
        if i_place_upper == 1:
            y2_index = i_place_upper - 1
            y2[y2_index] += y1[ibin]
            # print 'First bin %d counts go to y2 bin-%d' % (y1[ibin],y2_index)
            # print y2
        if i_place_upper == 2:
            print 'first bin 2 overlaps'

    # redistributing X1 bins with [1,m-1] indeces
    for ibin in range(1, m - 1):
        if y1[ibin] == 0:
            continue
        i_place_lower = i_place[ibin]
        i_place_upper = i_place[ibin + 1]
        x1min = x1[ibin]
        x1max = x1[ibin + 1]
        # Stop at the end
        if i_place_lower >= n - 2:
            return y2

        # X1 bin fully inside X2 bin:
        if i_place_lower == i_place_upper:
            y2_index = i_place_upper - 1
            # print 'y2_index %d -- i_place_lower=%d -- ibin=%d' % (y2_index,i_place_lower, ibin)
            y2[y2_index] += y1[ibin]
            # print '%d-th bin %d counts go to y2 bin-%d: X1 overlaps with X2' % (ibin,y1[ibin],y2_index)
            # print y2
        # X1 bin overlaps w/ two X2 bins:
        if i_place_lower == i_place_upper - 1:
            # X2 value that "splits" the X1 bin
            x2_ov1 = x2[i_place_lower]
            # Roll the dice y1[ibin]-times
            for irand in range(0, int(y1[ibin])):
                probValue = np.random.uniform(x1min, x1max)
                # print 'rand-%d probV=%f x2_ov1=%f' %(irand,probValue,x2_ov1)
                if probValue < x2_ov1:
                    # send photon to lower bin :))
                    y2_index = i_place_upper - 2
                    y2[y2_index] += 1
                else:
                    y2_index = i_place_upper - 1
                    y2[y2_index] += 1
            # print '%d-th bin %d counts go to y2 bin-%d: X1 is completely in X2' % (ibin,y1[ibin],y2_index)
            # print y2
        # X1 bin overplaps w/ three X2 bins
        if i_place_lower == i_place_upper - 2:
            print '2 overlap bins'

    return y2

def MBC_smoother(data):
    from scipy.interpolate import UnivariateSpline
    x = np.arange(0, len(data[0]))
    smooth1 = np.array([UnivariateSpline(x,y)(x) for y in data])
    ratios  = np.array(data)/smooth1
    np.nan_to_num(0.0)
    ratios = np.nan_to_num(ratios)
    smooth_ratios = np.array([UnivariateSpline(x,y)(x) for y in ratios])
    smooth2 = smooth_ratios*smooth1
    return np.array(smooth1),np.array(smooth2)

from scipy.interpolate import UnivariateSpline
import numpy as np
import h5py
import matplotlib.pyplot as plt

filename = '/Users/i6o/DHS_TSI/data/train_set/archerair2012-7606_Train_4_Detectors_Co60_2018-11-08T18.48.32.419.hdf5'
spectra = h5py.File(filename, 'r')['MUSE6_Spectra']
spectra = np.array([rebin_Archer(np.arange(0, 2048), x, np.arange(0, 1001)) for x in spectra])
smooth_spectra1,smooth_spectra2 = MBC_smoother(spectra)
spectra = np.array([np.add.reduceat(i, np.arange(0, len(i), 10)) for i in spectra])
CPS = np.sum(spectra,axis=1)

def kalman_filter(data,window):
    # smoothing data
    x = np.arange(0, len(data[0]))
    smooth1 = np.array([UnivariateSpline(x,y)(x) for y in data])
    ratios  = np.array(data)/smooth1
    np.nan_to_num(0.0)
    ratios = np.nan_to_num(ratios)
    smooth_ratios = np.array([UnivariateSpline(x,y)(x) for y in ratios])
    data_smooth = smooth_ratios*smooth1
    #
    # estimating Q parameter in Kalman Filter
    #
    Q_outside = (1/(window-2))
    Q_0 = 0
    for i in range(2,window):
        Q_0+=(data_smooth[i]-data_smooth[i-1])*(data_smooth[i]-data_smooth[i-1]).T
    Q = Q_outside*Q_0
    #
    # estimating R parameter in Kalman Filter
    #
    R_outside = (1/(window-1))
    R_0 = 0
    for i in range(1,window):
        R_0+=(data[i]-data[i-1])*(data[i]-data[i-1]).T
    R = R_outside*R_0
    #
    # initiating b_hat,Sigma,x_hat,S_hat
    #
    
