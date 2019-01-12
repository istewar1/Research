import numpy as np
from scipy.spatial import distance
from scipy.linalg import svd
def pca_spectra(spectral_data,threshold,plotter=False)
    '''
    
    Parameters
    ----------
    spectral_data
    threshold
    plotter

    Returns
    -------

    '''
    spectra_background = spectral_data[0:len(spectral_data[0])]
    U, s, VT = svd(spectra_background)
    Sigma = np.zeros((spectra_background.shape[0], spectra_background.shape[1]))
    Sigma[:background.shape[1], :background.shape[1]] = np.diag(s)
    Y = U.dot(Sigma.dot(VT))  # original data into PC space