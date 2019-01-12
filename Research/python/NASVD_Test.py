import pandas as pd

def nasvd(spectra_data):
    from sklearn.decomposition import TruncatedSVD
    import numpy as np
    '''

    Parameters
    ----------
    spectra_data

    Returns
    -------

    '''
    summed_spectra = np.sum(spectra_data, axis=1) # summed spectra (counts in spectra)
    norm_spectra = [np.array(x)/y for x,y in zip(spectra_data,summed_spectra)] # normalized spectra
    #import matplotlib.pyplot as plt;plt.figure();plt.plot(spectra_data[0],'b');plt.plot(norm_spectra[0],'r');plt.yscale('log');plt.show()
    total_summed_spectra = np.sum(spectra_data,axis=0) # normalized unit spectra
    unit_norm_spectra = np.array(norm_spectra) / np.array(total_summed_spectra)
    denom = np.array([np.array(x)*y for x,y in zip(unit_norm_spectra,summed_spectra)])
    na_spectra = spectra_data / np.sqrt(denom)
    na_spectra = np.array([np.array(np.nan_to_num(x)) for x in na_spectra])
    #
    # SVD of noise-adjusted spectra
    #
    svd_tr = TruncatedSVD(n_components=99)
    return svd_tr.fit(na_spectra)


if __name__ == "__main__":
    import h5py; import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from scipy.linalg import svd
    filename = '/Volumes/Ian External HD/Node Data/hdf5files/MUSE01/MUSE01-2018-08-01T00.35.06.865.hdf5'
    MUSE01 = h5py.File(filename,'r')['2x4x16Spectra'][0:10000]
    MUSE01 = np.array([np.add.reduceat(x, np.arange(0, 1000, 10)) for x in MUSE01])
    #MUSE01 = np.add.reduceat(MUSE01, np.arange(0, len(MUSE01), 5),axis=0)
    SVD = nasvd(MUSE01)
    U = SVD.fit_transform(MUSE01)
    Sigma = SVD.explained_variance_ratio_
    VT = SVD.components_

    #
    '''
    A = MUSE01
    U,s,VT = svd(A)
    # create m x n Sigma matrix
    Sigma = np.zeros((A.shape[0], A.shape[1]))
    # populate Sigma with n x n diagonal matrix
    Sigma[:A.shape[1], :A.shape[1]] = np.diag(s)
    A = U.dot(Sigma.dot(VT))
    '''
    var_ex = [] ; first = True
    for i in range(len(Sigma)):
        j = Sigma[i]
        if first:
            var_ex.append(j)
            first = False
        else:
            var_ex.append(var_ex[i-1]+j)
    #
    makeSCREE = True
    if makeSCREE:
        fig,ax1 = plt.subplots()
        ax1.bar(np.arange(1,len(Sigma)+1),Sigma,label='Explained Variance')
        ax1.set_xlabel('Component No.')
        ax1.set_ylabel('Explained Variance')
        ax1.set_yscale('linear')
        ax1.yaxis.label.set_color('blue')
        ax1.tick_params(axis='y', colors='blue')
        ax2 = ax1.twinx()
        ax2.plot(var_ex,'r--',label='Cumulative Variance')
        ax2.set_ylabel('Cumulative Variance')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
        # customizing legend box
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='center right',fancybox=False,shadow=True)

        ax1.grid(alpha=0.5)
        ax1.set_axisbelow(True)
        plt.show()
        plt.savefig('/Users/i6o/DHS_TSI/figures/Scree_plot_LOG_ORNL_data.png',dpi=600)
    #
    fig,ax = plt.subplots(10,1,sharex=True)
    for i in range(10):
        ax[i].plot(MUSE01[i],label=str(i))
        ax[i].plot(B[i],'r--',label=str(i))
        ax[i].legend()
        ax[i].set_yscale('linear')
        #ax[i].set_ylim(0.1,150)
        ax[i].grid(alpha=0.5,which='both')
    plt.show(True)
    print [sum(x-y) for x,y in zip(MUSE01[0:10],A[0:10])]
    #
    makePCcomparison = False
    if makePCcomparison:
        fig,ax = plt.subplots(nrows=2,ncols=3)
        ax[0,0].plot(VT[0],label='PC1');ax[0,0].legend();ax[0,0].grid(which='both');ax[0,0].set_axisbelow(True)
        ax[0,1].plot(VT[1], label='PC2');ax[0,1].legend();ax[0,1].grid(which='both');ax[0,1].set_axisbelow(True)
        ax[0,2].plot(VT[2], label='PC3');ax[0,2].legend();ax[0,2].grid(which='both');ax[0,2].set_axisbelow(True)
        ax[1,0].plot(VT[3], label='PC4');ax[1,0].legend();ax[1,0].grid(which='both');ax[1,0].set_axisbelow(True)
        ax[1,1].plot(VT[4], label='PC5');ax[1,1].legend();ax[1,1].grid(which='both');ax[1,1].set_axisbelow(True)
        ax[1,2].plot(VT[5], label='PC6');ax[1,2].legend();ax[1,2].grid(which='both');ax[1,2].set_axisbelow(True)
        fig,ax = plt.subplots(nrows=2,ncols=3)
        ax[0,0].plot(U[0],label='PC1');ax[0,0].legend();ax[0,0].grid(which='both');ax[0,0].set_axisbelow(True)
        ax[0,1].plot(U[1], label='PC2');ax[0,1].legend();ax[0,1].grid(which='both');ax[0,1].set_axisbelow(True)
        ax[0,2].plot(U[2], label='PC3');ax[0,2].legend();ax[0,2].grid(which='both');ax[0,2].set_axisbelow(True)
        ax[1,0].plot(U[3], label='PC4');ax[1,0].legend();ax[1,0].grid(which='both');ax[1,0].set_axisbelow(True)
        ax[1,1].plot(U[4], label='PC5');ax[1,1].legend();ax[1,1].grid(which='both');ax[1,1].set_axisbelow(True)
        ax[1,2].plot(U[5], label='PC6');ax[1,2].legend();ax[1,2].grid(which='both');ax[1,2].set_axisbelow(True)
        fig, ax = plt.subplots(nrows=2, ncols=3)
        ax[0,0].plot(U[0],U[1],'o',label='PC1 vs PC2',markersize=0.75,mfc='none', markeredgecolor='r');ax[0,0].grid(which='both');ax[0,0].set_axisbelow(True);ax[0,0].set_xlabel('PC1');ax[0,0].set_ylabel('PC2')
        ax[0,1].plot(U[0],U[2],'o',label='PC1 vs PC3',markersize=0.75,mfc='none', markeredgecolor='r');ax[0,1].grid(which='both');ax[0,1].set_axisbelow(True);ax[0,1].set_xlabel('PC1');ax[0,1].set_ylabel('PC3')
        ax[0,2].plot(U[0],U[3],'o',label='PC1 vs PC4',markersize=0.75,mfc='none', markeredgecolor='r');ax[0,2].grid(which='both');ax[0,2].set_axisbelow(True);ax[0,2].set_xlabel('PC1');ax[0,2].set_ylabel('PC4')
        ax[1,0].plot(U[0],U[4],'o',label='PC1 vs PC5',markersize=0.75,mfc='none', markeredgecolor='r');ax[1,0].grid(which='both');ax[1,0].set_axisbelow(True);ax[1,0].set_xlabel('PC1');ax[0,0].set_ylabel('PC5')
        ax[1,1].plot(U[0],U[5],'o',label='PC1 vs PC6',markersize=0.75,mfc='none', markeredgecolor='r');ax[1,1].grid(which='both');ax[1,1].set_axisbelow(True);ax[1,1].set_xlabel('PC1');ax[1,1].set_ylabel('PC6')
        ax[1,2].plot(U[0],U[6],'o',label='PC1 vs PC7',markersize=0.75,mfc='none', markeredgecolor='r');ax[1,2].grid(which='both');ax[1,2].set_axisbelow(True);ax[1,2].set_xlabel('PC1');ax[1,2].set_ylabel('PC7')
        plt.show()
    #
    # Making 2D histograms
    #
    plt.figure()
    x = np.array([]); y=np.array([]); first = True; count = 0; j = 0
    for i in range(100):
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
    bins_x = np.arange(0,1000)
    bins_y = np.arange(0,i)
    plt.hist2d(x,y,bins=(bins_x,bins_y),norm=LogNorm())
    plt.colorbar()
    plt.savefig('/Users/i6o/Desktop/unsmoothed_spectra.png',dpi=400)
    plt.close('all')
    #
    U0 = np.zeros(1000)
    for array in VT[0:10]:
        U0+=array
    x = np.array([]); y=np.array([]); first = True; count = 0; j = 0
    for i in range(100):
        count = int(count)
        j+=1
        if j == 21:
            print i
            j = 0
        if first:
            value_counter = 0
            for element in U0[count]:
                if element > 0.:
                    element = int(element)
                    x = np.concatenate((x,np.ones(element)+value_counter))
                    y = np.concatenate((y, np.zeros(element)))
                value_counter+=1
            first = False
        else:
            value_counter = 0
            for element in U0[count]:
                if element > 0.:
                    element = int(element)
                    x = np.concatenate((x, np.ones(element) + value_counter))
                    y = np.concatenate((y,np.zeros(element)+count))
                value_counter += 1
        count += 1
    bins_x = np.arange(0,1000)
    bins_y = np.arange(0,i)
    plt.figure()
    plt.hist2d(x,y,bins=(bins_x,bins_y),norm=LogNorm())
    plt.colorbar()
    plt.savefig('/Users/i6o/Desktop/NASVD_unsmooted_spectra.png',dpi=400)
    #
    plt.show(False)


