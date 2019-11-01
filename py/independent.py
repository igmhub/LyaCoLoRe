import numpy as np
import time

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
def get_gaussian_fields(generator,N_cells,z=0.0,dv_kms=10.0,N_skewers=1,white_noise=False,n=0.7,k1=0.001,A0=58.6,R_kms=25.0,norm=True,remove_P1D_data=None,remove_P1D_amp=None):
    #print(generator,N_cells,z,dv_kms,N_skewers,white_noise,n,k1,A0)

    times = []
    start = time.time(); times += [start]
    # number of Fourier modes
    NF = int(N_cells/2+1)

    # get frequencies (wavenumbers in units of s/km)
    k_kms = np.fft.rfftfreq(N_cells)*2*np.pi/dv_kms

    # get power evaluated at each k_kms
    P_kms = power_kms(z,k_kms,dv_kms,white_noise=white_noise,n=n,k1=k1,A0=A0,R_kms=R_kms,smooth=True,norm=norm)
    #P_kms = alternative_power_kms(z,k_kms,dv_kms,A0=A0,k0=k1,E1=n,E2=-0.1,R1=R1,smooth=True)

    if remove_P1D_data is not None:
        k_rem = remove_P1D_data['k']
        Pk_rem = remove_P1D_data['Pk']
        Pk_rem = np.interp(k_kms,k_rem,Pk_rem)
        P_kms -= Pk_rem/remove_P1D_amp

    times += [time.time()]
    # generate random Fourier modes
    modes = np.empty([N_skewers,NF], dtype=complex)
    modes[:].real = np.reshape(generator.normal(size=N_skewers*NF),[N_skewers,NF])
    modes[:].imag = np.reshape(generator.normal(size=N_skewers*NF),[N_skewers,NF])
    times += [time.time()]

    #print('rand numbers size',modes.shape)
    #print('start of rand numbers',modes[0,:10])
    # normalize to desired power (and enforce real for i=0, i=NF-1)
    modes[:,0] = modes[:,0].real * np.sqrt(P_kms[0])
    modes[:,-1] = modes[:,-1].real * np.sqrt(P_kms[-1])
    modes[:,1:-1] *= np.sqrt(0.5*P_kms[1:-1])
    times += [time.time()]
    # inverse FFT to get (normalized) delta field
    delta = np.fft.irfft(modes,n=N_cells) * np.sqrt(N_cells/dv_kms)

    #check
    #pk_rows = np.fft.rfft(delta,axis=1) / np.sqrt(N_cells/dv_kms)
    #pk_rows = np.abs(pk_rows)**2
    #pk_measured = np.average(pk_rows,axis=0)
    #print(k_kms)
    #print('Measured',pk_measured)
    #print('Input   ',pk_measured)

    times += [time.time()]
    #print('sigma of Pk added (no smoothing)',np.sqrt((1/np.pi)*np.trapz(power_kms(z,k_kms,dv_kms,white_noise=white_noise,n=n,k1=k1,A0=A0,smooth=False,norm=True),k_kms)))
    #print('sigma of Pk added (smoothing)   ',np.sqrt((1/np.pi)*np.trapz(power_kms(z,k_kms,dv_kms,white_noise=white_noise,n=n,k1=k1,A0=A0,smooth=True,norm=True),k_kms)))
    #print('sigma measured inside ggf       ',np.sqrt((1/np.pi)*np.trapz(pk_measured,k_kms)))
    #print('std measured inside ggf         ',np.std(delta))
    #print(' ')

    #print(P_kms)
    #print(pk_measured)
    #print(np.array(times)-start)
    return delta

#Function to return a gaussian P1D in k.
#From lya_mock_functions
def power_amplitude(z,A0=58.6):
    """Add redshift evolution to the Gaussian power spectrum."""
    return A0*pow((1+z)/4.0,-2.82)

#Function to return a gaussian P1D in k.
#From lya_mock_functions
def power_kms(z_c,k_kms,dv_kms,white_noise=False,n=0.7,k1=0.001,A0=58.6,R_kms=25.0,smooth=True,norm=False):
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
    #R1 = 5.0
    # compute term without smoothing
    P = A * (1.0+pow(0.01/k1,n)) / (1.0+pow(k_kms/k1,n))
    if smooth:
        # smooth with Gaussian and top hat
        kdv = np.fmax(k_kms*dv_kms,0.000001)
        P *= np.exp(-pow(k_kms*R_kms,2)) * pow(np.sin(kdv/2)/(kdv/2),2)
    if norm:
        # normalise to have sigma of 1
        sigma2 = (1/np.pi) * np.trapz(P,k_kms)
        P /= sigma2
    return P

def alternative_power_kms(z_c,k_kms,dv_kms,white_noise=False,A0=58.6,k0=0.009,E1=-0.55,E2=-0.1,R_kms=25.0,smooth=True,norm=False):

    A = power_amplitude(z_c,A0=A0)
    P = np.zeros(k_kms.shape)
    cells = k_kms>0
    P[cells] = A * (k_kms[cells]/k0) ** (E1 + E2*np.log(k_kms[cells]/k0))
    #Flatten the low_k:
    k_min = k_kms[np.argmax(P)]
    k_kms = np.fmax(k_kms,k_min)
    P = A * (k_kms/k0) ** (E1 + E2*np.log(k_kms/k0))
    if smooth:
        # smooth with Gaussian and top hat
        kdv = np.fmax(k_kms*dv_kms,0.000001)
        P *= np.exp(-pow(k_kms*R_kms,2)) * pow(np.sin(kdv/2)/(kdv/2),2)
    if norm:
        # normalise to have sigma of 1
        sigma2 = (1/np.pi) * np.trapz(P,k_kms)
        P /= sigma2

    return P

def get_sigma_G(z_c,k_kms,dv_kms,white_noise=False,n=0.7,k1=0.001,A0=58.6):

    Pk_kms = power_kms(z_c,k_kms,dv_kms,white_noise=False,n=0.7,k1=0.001,A0=58.6)

    return sigma_G
