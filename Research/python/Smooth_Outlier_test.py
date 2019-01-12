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

if __name__ == "__main__":
    #
    import h5py
    from scipy.spatial import distance
    from scipy.interpolate import UnivariateSpline
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from MBC_smoother import MBC_smoother
    #
    filename = '/Users/i6o/DHS_TSI/data/train_set/archerair2012-7606_Train_4_Detectors_Co60_2018-11-08T18.48.32.419.hdf5'
    MUSE01 = h5py.File(filename,'r')['MUSE6_Spectra']
    MUSE01 = np.array([rebin_Archer(np.arange(0,2048),x,np.arange(0,1000)) for x in MUSE01])
    MUSE01 = np.array([np.add.reduceat(i, np.arange(0, len(i), 10)) for i in MUSE01])
    CPS = np.sum(MUSE01,axis=1)
    #
    x = np.arange(0, len(MUSE01[0]))
    smooth1 = np.array([UnivariateSpline(x,y)(x) for y in MUSE01])
    ratios  = np.array(MUSE01)/smooth1
    np.nan_to_num(0.0)
    ratios = np.nan_to_num(ratios)
    smooth_ratios = np.array([UnivariateSpline(x,y)(x) for y in ratios])
    smooth2 = smooth_ratios*smooth1
    #
    count = 0; D_regular = []; D_smooth = []; value = 100
    for i in range(len(MUSE01)):
        if count > value:
            # Mahalanobis Distance - MUSE01
            mean = np.mean((MUSE01[i-value:i]),axis=0)
            cov  = np.linalg.inv(np.cov(MUSE01[i-value:i]))
            D_regular.append(distance.mahalanobis(MUSE01[i],mean,cov))
            # Mahalanobis Distance - Smoother
            mean = np.mean((smooth2[i-value:i]),axis=0)
            cov  = np.linalg.inv(np.cov(smooth2[i-value:i]))
            D_smooth.append(distance.mahalanobis(smooth2[i],mean,cov))
        count+=1
    zeros = np.zeros(value)
    D_regular = np.concatenate((zeros,np.array(D_regular)))
    D_smooth = np.concatenate((zeros, np.array(D_smooth)))
    #
    fig, ax1 = plt.subplots()
    ax1.plot(D_regular,'bo',label='Mahalnobis - REGULAR',zorder=10,alpha=0.25)
    ax1.set_xlabel('Entry Value')
    ax1.set_ylabel('Mahalnobis Distance')
    ax1.legend(fancybox=False,shadow=True)
    ax1.set_yscale('linear')
    ax1.yaxis.label.set_color('blue')
    ax1.tick_params(axis='y', colors='blue')
    ax2 = ax1.twinx()
    ax2.plot(CPS, 'r', label='CPS')
    ax2.set_ylabel('CPs')
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
    #
    fig, ax1 = plt.subplots()
    ax1.plot(D_smooth,'bo',label='Mahalnobis - SMOOTH',zorder=10,alpha=0.25)
    ax1.set_xlabel('Entry Value')
    ax1.set_ylabel('Mahalnobis Distance')
    ax1.set_yscale('linear')
    ax1.legend(fancybox=False,shadow=True)
    ax1.yaxis.label.set_color('blue')
    ax1.tick_params(axis='y', colors='blue')
    ax2 = ax1.twinx()
    ax2.plot(CPS, 'r', label='CPS')
    ax2.set_ylabel('CPs')
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
    plt.show(False)
    #
    fig,ax = plt.subplots(5,1,sharex=True)
    for i in range(5):
        ax[i].plot(MUSE01[i],label='MUSE%i'%i)
        ax[i].plot(smooth2[i],label='Smooth2-%i'%i)
        ax[i].plot(smooth1[i],label='Smooth1-%i'%i)
        ax[i].set_yscale('linear')
        ax[i].legend()
        ax[i].grid(alpha=0.5); ax[i].set_axisbelow(True)
    #
    x = np.array([]); y=np.array([]); first = True; count = 0; j = 0
    for i in range(len(MUSE01)):
        j+=1
        if j == 21:
            print i
            j = 0
        if first:
            value_counter = 0
            for element in MUSE01[count]:
                if element > 0.:
                    element = int(element)
                    x = np.concatenate((x,np.ones(element)+value_counter))
                    y = np.concatenate((y, np.zeros(element)))
                value_counter+=1
            first = False
        else:
            value_counter = 0
            for element in MUSE01[count]:
                if element >0.:
                    element = int(element)
                    x = np.concatenate((x, np.ones(element) + value_counter))
                    y = np.concatenate((y,np.zeros(element)+count))
                value_counter += 1
        count += 1
    bins_x = np.arange(0,len(MUSE01[0]))
    bins_y = np.arange(0,i)
    plt.figure()
    plt.hist2d(x,y,bins=(bins_x,bins_y),norm=LogNorm())
    plt.colorbar()
    plt.savefig('/Users/i6o/Desktop/TRAIN_unsmoothed_spectra.png',dpi=400)
    #
    x = np.array([]); y=np.array([]); first = True; count = 0; j = 0
    for i in range(len(smooth2)):
        count = int(count)
        j+=1
        if j == 21:
            print i
            j = 0
        if first:
            value_counter = 0
            for element in smooth2[count]:
                if element > 0.:
                    element = int(element)
                    x = np.concatenate((x,np.ones(element)+value_counter))
                    y = np.concatenate((y,np.zeros(element)))
                value_counter+=1
            first = False
        else:
            value_counter = 0
            for element in smooth2[count]:
                if element > 0.:
                    element = int(element)
                    x = np.concatenate((x,np.ones(element)+value_counter))
                    y = np.concatenate((y,np.zeros(element)+count))
                value_counter += 1
        count += 1
    bins_x = np.arange(0,len(smooth2[0]))
    bins_y = np.arange(0,i)
    plt.figure()
    plt.hist2d(x,y,bins=(bins_x,bins_y),norm=LogNorm())
    plt.colorbar()
    plt.savefig('/Users/i6o/Desktop/TRAIN_smoothed_spectra.png',dpi=400)
    #
    plt.show()