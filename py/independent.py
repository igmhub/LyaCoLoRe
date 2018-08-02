import numpy as np
from astropy.io import fits

lya = 1215.67

#Function to generate random Gaussian skewers with a given standard deviation.
def get_gaussian_skewers(generator,N_cells,sigma_G=1.0,N_skewers=1):

    if N_cells*N_skewers == 1:
        size = 1
    else:
        size=(N_skewers,N_cells)

    """
    gaussian_skewers = fast_prng.normal(size=size,scale=sigma_G)
    """

    gaussian_skewers = generator.normal(size=(N_skewers,N_cells),scale=sigma_G)

    return gaussian_skewers

#Function to generate random Gaussian fields at a given redshift.
#From lya_mock_functions
def get_gaussian_fields(generator,N_cells,z=0.0,dv_kms=10.0,N_skewers=1,white_noise=True,n=0.7,k1=0.001,A0=58.6):

    # number of Fourier modes
    NF = int(N_cells/2+1)

    # get frequencies (wavenumbers in units of s/km)
    k_kms = np.fft.rfftfreq(N_cells)*2*np.pi/dv_kms

    # get power evaluated at each k_kms
    P_kms = power_kms(z,k_kms,dv_kms,white_noise=white_noise,n=n,k1=k1,A0=A0)

    # generate random Fourier modes
    modes = np.empty([N_skewers,NF], dtype=complex)
    modes[:].real = np.reshape(generator.normal(size=N_skewers*NF),[N_skewers,NF])
    modes[:].imag = np.reshape(generator.normal(size=N_skewers*NF),[N_skewers,NF])

    # normalize to desired power (and enforce real for i=0, i=NF-1)
    modes[:,0] = modes[:,0].real * np.sqrt(P_kms[0])
    modes[:,-1] = modes[:,-1].real * np.sqrt(P_kms[-1])
    modes[:,1:-1] *= np.sqrt(0.5*P_kms[1:-1])

    # inverse FFT to get (normalized) delta field
    delta = np.fft.irfft(modes,n=N_cells) * np.sqrt(N_cells/dv_kms)

    return delta

#Function to return a gaussian P1D in k.
#From lya_mock_functions
def power_amplitude(z,A0=58.6):
    """Add redshift evolution to the Gaussian power spectrum."""
    return A0*pow((1+z)/4.0,-2.82)

#Function to return a gaussian P1D in k.
#From lya_mock_functions
def power_kms(z_c,k_kms,dv_kms,white_noise,n=0.7,k1=0.001,A0=58.6):
    """Return Gaussian P1D at different wavenumbers k_kms (in s/km), fixed z_c.

      Other arguments:
        dv_kms: if non-zero, will multiply power by top-hat kernel of this width
        white_noise: if set to True, will use constant power of 100 km/s
    """
    if white_noise: return np.ones_like(k_kms)*100.0
    # power used to make mocks in from McDonald et al. (2006)
    A = power_amplitude(z_c,A0=A0)
    #k1 = 0.001
    #n = 0.7
    R1 = 5.0
    # compute term without smoothing
    P = A * (1.0+pow(0.01/k1,n)) / (1.0+pow(k_kms/k1,n))
    # smooth with Gaussian and top hat
    kdv = np.fmax(k_kms*dv_kms,0.000001)
    P *= np.exp(-pow(k_kms*R1,2)) * pow(np.sin(kdv/2)/(kdv/2),2)
    return P
