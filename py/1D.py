import numpy as np
from astropy.io import fits

#Function to generate random Gaussian skewers with a given standard deviation.
def get_gaussian_skewers(generator,N_cells,sigma_G=1.0,N_skewers=1):

    if N_cells*N_skewers == 1:
        size = 1
    else:
        size=(N_skewers,N_cells)

    gaussian_skewers = fast_prng.normal(size=size,scale=sigma_G)

    """
    gaussian_skewers = np.random.normal(size=(N_skewers,N_cells),scale=sigma_G)
    """
    return gaussian_skewers

#Function to generate random Gaussian fields at a given redshift.
#From lya_mock_functions
def get_gaussian_fields(z,N_cells,dv_kms=10.0,N_skewers=1,new_seed=None,white_noise=True):
    """Generate N_skewers Gaussian fields at redshift z_c.

      If new_seed is set, it will reset random generator with it."""

    random_state = np.random.RandomState(new_seed)

    # number of Fourier modes
    NF = int(N_cells/2+1)

    # get frequencies (wavenumbers in units of s/km)
    k_kms = np.fft.rfftfreq(N_cells)*2*np.pi/dv_kms

    # get power evaluated at each k_kms
    P_kms = power_kms(z,k_kms,dv_kms,white_noise=white_noise)

    # generate random Fourier modes
    modes = np.empty([N_skewers,NF], dtype=complex)
    modes[:].real = np.reshape(random_state.normal(size=N_skewers*NF),[N_skewers,NF])
    modes[:].imag = np.reshape(random_state.normal(size=N_skewers*NF),[N_skewers,NF])

    # normalize to desired power (and enforce real for i=0, i=NF-1)
    modes[:,0] = modes[:,0].real * np.sqrt(P_kms[0])
    modes[:,-1] = modes[:,-1].real * np.sqrt(P_kms[-1])
    modes[:,1:-1] *= np.sqrt(0.5*P_kms[1:-1])

    # inverse FFT to get (normalized) delta field
    delta = np.fft.irfft(modes,n=N_cells) * np.sqrt(N_cells/dv_kms)

    return delta

#Function to return a gaussian P1D in k.
#From lya_mock_functions
def power_kms(z_c,k_kms,dv_kms,white_noise):
    """Return Gaussian P1D at different wavenumbers k_kms (in s/km), fixed z_c.

      Other arguments:
        dv_kms: if non-zero, will multiply power by top-hat kernel of this width
        white_noise: if set to True, will use constant power of 100 km/s
    """
    if white_noise: return np.ones_like(k_kms)*100.0
    # power used to make mocks in from McDonald et al. (2006)
    A = power_amplitude(z_c)
    k1 = 0.001
    n = 0.7
    R1 = 5.0
    # compute term without smoothing
    P = A * (1.0+pow(0.01/k1,n)) / (1.0+pow(k_kms/k1,n))
    # smooth with Gaussian and top hat
    kdv = np.fmax(k_kms*dv_kms,0.000001)
    P *= np.exp(-pow(k_kms*R1,2)) * pow(np.sin(kdv/2)/(kdv/2),2)
    return P

#Function to integrate under the 1D power spectrum to return the value of sigma_dF at a given redshift.
def get_sigma_dF_P1D(z,l_hMpc=0.25,Om=0.3):
    #Choose log spaced values of k
    k_hMpc_max = 100.0/l_hMpc
    k_hMpc = np.logspace(-5,np.log10(k_hMpc_max),10**5)

    # TODO: generalise the conversion in here
    # need to go from Mpc/h to km/s, using dv / dX = H(z) / (1+z)
    # we will define H(z) = 100 h E(z)
    # with E(z) = sqrt(Omega_m(1+z)^3 + Omega_L), and assume flat universe
    E_z = np.sqrt(Om*(1+z)**3 + (1-Om))
    dkms_dhMpc = 100. * E_z / (1+z)

    # transform input wavenumbers to s/km
    k_kms = k_hMpc / dkms_dhMpc

    # get power in units of km/s
    pk_kms = P1D_z_kms_PD2013(z,k_kms)

    # transform to h/Mpc
    pk_hMpc = pk_kms / dkms_dhMpc

    # compute Fourier transform of Top-Hat filter of size l_hMpc
    W_hMpc = np.sinc((k_hMpc*l_hMpc)/(2*np.pi))
    sigma_dF = np.sqrt((1/np.pi)*np.trapz((W_hMpc**2)*pk_hMpc,k_hMpc))

    return sigma_dF
