import numpy as np

# code to make mock Lya spectra following McDonald et al. (2006)
# copied from c++ code in Cosmology/LNLyaF

# adapted and added to by James
#TO DO:
#make the mask maker into a function of its own?
#put histogram maker into a function of its own
#pass density stats the whole hdulist rather than just the density skewers?

#Produce a binned average density skewer
def density_stats(z,z_qso,density_skewers):
    """Determine the average density, density variance and delta variance over a number of redshift bins"""

    N_cells = density_skewers.shape[1]
    density_squared_skewers = density_skewers**2
    N_qso = len(z_qso)

    #Make a mask to remove irrelevant data (densities for z>z_qso)
    mask = np.ones(density_skewers.shape)
    max_pixel_qso = [0.]*N_qso
    for j in range(N_qso):
        max_pixel_qso[j] = (np.argmax(z>z_qso[j]))%N_cells
        mask[j,max_pixel_qso[j]+1:]=np.zeros(1,(mask.shape[1]-max_pixel_qso[j]-1))

    #Calculate the average density over all skewers
    mean_density = np.zeros(N_cells)
    ignore = []
    for j in range(N_cells):
        if sum(mask[:,j])!=0:
            mean_density[j] = np.average(density_skewers[:,j],weights=mask[:,j])
        else:
            ignore += [j]

    #Calculate the average density squared over all skewers
    mean_density_squared = np.zeros(N_cells)
    for j in range(N_cells):
        if sum(mask[:,j])!=0:
            mean_density_squared[j] = np.average(density_squared_skewers[:,j],weights=mask[:,j])

    #Set up the binned data structure
    N_bins = 100
    binned_z = np.zeros(N_bins)
    binned_z_location = np.zeros(N_bins)
    binned_mean_density = np.zeros(N_bins)
    binned_mean_density_alt = np.zeros(N_bins)
    binned_density_var = np.zeros(N_bins)
    binned_delta_var = np.zeros(N_bins)

    if N_cells<N_bins:
        #ERROR, not enough cells for the N_bins specified
        exit('Not enough cells for the N_bins specified.')

    #Construct the conversion matrix: an N_bins by N_cells matrix detailing the contribution of each cell to each bin
    #conversion[i,j] denotes the contribution of the jth cell to the ith bin
    conversion = np.zeros((N_bins,N_cells))
    ratio = N_cells/N_bins
    for i in range(N_bins):
        unassigned_row_contribution=ratio
        for j in range(N_cells):
            unassigned_col_contribution=1-sum(conversion[:,j])
            if unassigned_row_contribution>0:
                conversion[i,j] = min(unassigned_row_contribution,unassigned_col_contribution,1)
                unassigned_row_contribution -= conversion[i,j]
                unassigned_col_contribution -= conversion[i,j]
            else:
                break
        #Construct output vectors
        binned_z_location[i] = (i+0.5)*ratio + 0.5

        binned_mean_density[i] = sum(conversion[i,:]*mean_density)/ratio
        binned_density_var[i] = sum(conversion[i,:]*mean_density_squared)/ratio - (binned_mean_density[i])**2
        binned_delta_var[i] = binned_density_var[i]/((binned_mean_density[i])**2)

    #Calculate the binned z values by linear interpolation
    binned_z = np.interp(binned_z_location,list(range(N_cells)),z)

    for j in range(N_bins):
        conversion_mask = mask*conversion[j,:]
        binned_mean_density_alt[j] = sum(sum(density_skewers*conversion_mask))/sum(sum(conversion_mask))

    return binned_z, binned_mean_density, binned_density_var, binned_delta_var

def power_amplitude(z):
    """Add redshift evolution to the Gaussian power spectrum."""
    return 58.6*pow((1+z)/4.0,-2.82)

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

def get_tau(z,density):
    """transform lognormal density to optical depth, at each z"""
    # add redshift evolution to mean optical depth
    A = 0.374*pow((1+z)/4.0,5.10)
    return A*density

class MockMaker(object):
    """Class to generate 1D mock Lyman alpha skewers."""

    # central redshift, sets center of skewer and pivot point in z-evolution
    z_c=1.9

    def __init__(self, N2=15, dv_kms=10.0, seed=666, white_noise=False):
        """Construct object to make 1D Lyman alpha mocks.

          Optional arguments:
            N2: the number of cells in the skewer will be 2^N2
            dv_kms: cell width (in km/s)
            seed: starting seed for the random number generator
            white_noise: use constant power instead of realistic P1D. """
        self.N = np.power(2,N2)
        self.dv_kms = dv_kms
        # setup random number generator using seed
        self.gen = np.random.RandomState(seed)
        self.white_noise = white_noise

    def get_density(self,var_delta,z,delta):
        """Transform Gaussian field delta to lognormal density, at each z."""
        tau_pl=2.0
        # relative amplitude
        rel_amp = power_amplitude(z)/power_amplitude(self.z_c)
        return np.exp(tau_pl*(delta*np.sqrt(rel_amp)-0.5*var_delta*rel_amp))

    def get_redshifts(self):
        """Get redshifts for each cell in the array (centered at z_c)."""
        N = self.N
        L_kms = N * self.dv_kms
        c_kms = 2.998e5
        if (L_kms > 4 * c_kms):
            print('Array is too long, approximations break down.')
            raise SystemExit
        # get indices
        i = range(N)
        z = (1+self.z_c)*pow(1-(i-N/2+1)*self.dv_kms/2.0/c_kms,-2)-1
        return z

    def get_gaussian_fields(self,Ns=1,new_seed=None):
        """Generate Ns Gaussian fields at redshift z_c.

          If new_seed is set, it will reset random generator with it."""
        if new_seed:
            self.gen = np.random.RandomState(new_seed)
        # length of array
        N = self.N
        # number of Fourier modes
        NF=int(N/2+1)
        # get frequencies (wavenumbers in units of s/km)
        k_kms = np.fft.rfftfreq(N)*2*np.pi/self.dv_kms
        # get power evaluated at each k
        P_kms = power_kms(self.z_c,k_kms,self.dv_kms,self.white_noise)
        # compute also expected variance, will be used in lognormal transform
        dk_kms = 2*np.pi/(N*self.dv_kms)
        var_delta=np.sum(P_kms)*dk_kms/np.pi
        # Nyquist frecuency is counted twice in variance, and it should not be
        var_delta *= NF/(NF+1)

        # generate random Fourier modes
        modes = np.empty([Ns,NF], dtype=complex)
        modes[:].real = np.reshape(self.gen.normal(size=Ns*NF),[Ns,NF])
        modes[:].imag = np.reshape(self.gen.normal(size=Ns*NF),[Ns,NF])
        # normalize to desired power (and enforce real for i=0, i=NF-1)
        modes[:,0] = modes[:,0].real * np.sqrt(P_kms[0])
        modes[:,-1] = modes[:,-1].real * np.sqrt(P_kms[-1])
        modes[:,1:-1] *= np.sqrt(0.5*P_kms[1:-1])
        # inverse FFT to get (normalized) delta field
        delta = np.fft.irfft(modes) * np.sqrt(N/self.dv_kms)

        return delta, var_delta

    def get_lya_skewers(self,Ns=10,new_seed=None):
        """Return Ns Lyman alpha skewers (wavelength, flux).

          If new_seed is set, it will reset random generator with it."""
        if new_seed:
            self.gen = np.random.RandomState(new_seed)
        # get redshift for all cells in the skewer
        z = self.get_redshifts()
        delta, var_delta = self.get_gaussian_fields(Ns)
        #var_delta = np.var(delta)
        density = self.get_density(var_delta,z,delta)
        tau = get_tau(z,density)
        flux = np.exp(-tau)
        wave = 1215.67*(1+z)
        return wave, flux
