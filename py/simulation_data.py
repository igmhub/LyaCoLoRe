import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
import time

from . import utils, read_files, bias, convert, RSD, DLA, independent, absorber, absorber_data, stats

lya = utils.lya_rest

#Function to create a SimulationData object given a specific pixel, information about the complete simulation, and the filenames.
def make_gaussian_pixel_object(pixel,base_filename,input_format,MOCKID_lookup,lambda_min=0,IVAR_cutoff=lya):

    #Determine which file numbers we need to look at for the current pixel.
    relevant_keys = [key for key in MOCKID_lookup.keys() if key[1]==pixel]
    files_included = 0

    #For each relevant file, extract the data and aggregate over all files into a 'combined' object.
    for key in relevant_keys:
        #Get the MOCKIDs of the relevant quasars: those that are located in the current pixel, stored in the current file.
        file_number = key[0]
        relevant_MOCKIDs = MOCKID_lookup[key]
        N_relevant_qso = len(relevant_MOCKIDs)

        #If there are some relevant quasars, open the data file and make it into a SimulationData object.
        #We use SimulationData.get_reduced_data to avoid loading all of the file's data into the object.
        if N_relevant_qso > 0:
            filename = base_filename+'{}.fits'.format(file_number)
            working = SimulationData.get_gaussian_skewers_object(filename,file_number,input_format,MOCKIDs=relevant_MOCKIDs,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff)

        #Combine the data from the working file with that from the files already looked at.
        if files_included > 0:
            combined = SimulationData.combine_files(combined,working,gaussian_only=True)
            files_included += 1
        else:
            combined = working
            files_included += 1

    pixel_object = combined

    return pixel_object

#Definition of a generic SimulationData class, from which it is easy to save in new formats.
class SimulationData:
    #Initialisation function.
    def __init__(self,N_qso,N_cells,SIGMA_G,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP):

        self.N_qso = N_qso
        self.N_cells = N_cells
        self.SIGMA_G = SIGMA_G

        # catalog information
        self.TYPE = TYPE
        self.RA = RA
        self.DEC = DEC
        self.Z_QSO = Z_QSO
        self.DZ_RSD = DZ_RSD
        self.MOCKID = MOCKID
        self.PLATE = PLATE
        self.MJD = MJD
        self.FIBER = FIBER

        # skewer information used by all absorbers
        self.GAUSSIAN_DELTA_rows = GAUSSIAN_DELTA_rows
        self.DENSITY_DELTA_rows = DENSITY_DELTA_rows
        self.VEL_rows = VEL_rows

        # used in picca files to mask outside Lya region
        self.LOGLAM_MAP = LOGLAM_MAP
        self.IVAR_rows = IVAR_rows

        # coordinates for the skewer cells (might get rid of some of these)
        self.R = R
        self.Z = Z
        self.D = D
        self.V = V

        # these will store the absorbing fields (Lya, Lyb, metals...)
        self.lya_absorber = absorber.AbsorberData(name='Lya',rest_wave=lya,flux_transform_m=1.0)
        self.lyb_absorber = None
        self.metals = None

        # these will store the DLA (if asked for)
        self.DLA_table = None

        # these will store the RSD weights (if asked for)
        self.RSD_weights = None

        return

    def setup_Lyb_absorber(self):

        self.lyb_absorber = absorber_data.get_lyb_absorber()

        return

    def setup_metal_absorbers(self):

        # get a dictionary with multiple absorbers, one for each metal line
        self.metals = absorber_data.get_metal_dict()

        return

    #Method to extract reduced data from an input file of a given format, with a given list of MOCKIDs.
    @classmethod
    def get_gaussian_skewers_object(cls,filename,file_number,input_format,MOCKIDs=None,lambda_min=0,IVAR_cutoff=lya,SIGMA_G=None):

        h = fits.open(filename)

        #Get data about the catalog and cosmology.
        h_MOCKID = read_files.get_MOCKID(h,input_format,file_number)
        h_R, h_Z, h_D, h_V = read_files.get_COSMO(h,input_format)
        h_lya_lambdas = read_files.get_lya_lambdas(h,input_format)
        h_lya_lambdas_edges = utils.get_edges(h_lya_lambdas)

        #Work out which rows in the hdulist we are interested in.
        if MOCKIDs is not None:
            rows = []
            s = set(MOCKIDs)
            for i, qso in enumerate(h_MOCKID):
                if qso in s:
                    rows += [i]
        else:
            rows = list(range(h_MOCKID.shape[0]))

        #Calculate the first_relevant_cell.
        first_relevant_cell = np.searchsorted(h_lya_lambdas_edges[1:],lambda_min)

        if input_format == 'physical_colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE'][rows]
            RA = h[1].data['RA'][rows]
            DEC = h[1].data['DEC'][rows]
            Z_QSO = h[1].data['Z_COSMO'][rows]
            DZ_RSD = h[1].data['DZ_RSD'][rows]
            DENSITY_DELTA_rows = h[2].data[rows,:][:,first_relevant_cell:]
            VEL_rows = h[3].data[rows,:][:,first_relevant_cell:]
            Z = h[4].data['Z'][first_relevant_cell:]
            R = h[4].data['R'][first_relevant_cell:]
            D = h[4].data['D'][first_relevant_cell:]
            V = h[4].data['V'][first_relevant_cell:]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = Z.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive the MOCKID and LOGLAM_MAP.
            MOCKID = read_files.get_MOCKID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Calculate the Gaussian skewers.
            GAUSSIAN_DELTA_rows = convert.lognormal_delta_to_gaussian(DENSITY_DELTA_rows,SIGMA_G,D)

            #Make binary IVAR_rows for picca.
            IVAR_rows = utils.make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

            #Insert placeholder values for remaining variables.
            DENSITY_DELTA_rows = None
            PLATE = MOCKID
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)

        elif input_format == 'gaussian_colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE'][rows]
            RA = h[1].data['RA'][rows]
            DEC = h[1].data['DEC'][rows]
            Z_QSO = h[1].data['Z_COSMO'][rows]
            DZ_RSD = h[1].data['DZ_RSD'][rows]
            GAUSSIAN_DELTA_rows = h[2].data[rows,:][:,first_relevant_cell:]
            VEL_rows = h[3].data[rows,:][:,first_relevant_cell:]
            Z = h[4].data['Z'][first_relevant_cell:]
            R = h[4].data['R'][first_relevant_cell:]
            D = h[4].data['D'][first_relevant_cell:]
            V = h[4].data['V'][first_relevant_cell:]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = Z.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive the MOCKID and LOGLAM_MAP.
            if MOCKIDs != None:
                MOCKID = MOCKIDs
            else:
                MOCKID = read_files.get_MOCKID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Make binary IVAR_rows for picca.
            IVAR_rows = utils.make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

            #Insert placeholder values for remaining variables.
            DENSITY_DELTA_rows = None
            PLATE = MOCKID
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)
        else:
            print('Input format not recognised: current options are "colore" and "picca".')
            print('Please choose one of these options and try again.')

        h.close()

        return cls(N_qso,N_cells,SIGMA_G,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    #Function to trim skewers according to a minimum value of lambda. QSOs with no relevant cells are removed.
    def trim_skewers(self,lambda_min=0.,min_catalog_z=0.,extra_cells=0,lambda_max=None,whole_lambda_range=False,remove_irrelevant_QSOs=False):

        #Find the first relevant cell, and the last one if desired.
        #We make sure to include all frequencies within (lambda_min,lambda_max).
        lambdas = 10**(self.LOGLAM_MAP)
        lambda_edges = utils.get_edges(lambdas)
        first_relevant_cell = np.searchsorted(lambda_edges[1:],lambda_min)
        if lambda_max:
            last_relevant_cell = np.searchsorted(lambda_edges[1:],lambda_max) - 1
        else:
            last_relevant_cell = -1 % self.N_cells

        #Calculate the actual values of lambda min and max.
        actual_lambda_min = lambda_edges[:-1][first_relevant_cell]
        actual_lambda_max = lambda_edges[:-1][last_relevant_cell]

        #If we want to keep any extra_cells, we subtract from the first_relevant_cell.
        #If we cannot add enough extra cells, then we just set the first relevant cell to 0.
        if first_relevant_cell>extra_cells:
            first_relevant_cell -= extra_cells
        else:
            first_relevant_cell = 0

        #Determine which QSOs have any relevant cells to keep.
        relevant_QSOs = (self.Z_QSO>min_catalog_z)
        if remove_irrelevant_QSOs:
            relevant_QSOs *= (np.sum(self.IVAR_rows,axis=1)>0)

        #If we want the entirety of the lambda range to be relevant (i.e. with IVAR=1), we must remove skewers that do not have this
        if whole_lambda_range:
            relevant_QSOs *= (self.IVAR_rows[:,first_relevant_cell] == 1) * (self.IVAR_rows[:,last_relevant_cell] == 1)

        #Remove QSOs no longer needed.
        self.N_qso = np.sum(relevant_QSOs)

        self.TYPE = self.TYPE[relevant_QSOs]
        self.RA = self.RA[relevant_QSOs]
        self.DEC = self.DEC[relevant_QSOs]
        self.Z_QSO = self.Z_QSO[relevant_QSOs]
        self.DZ_RSD = self.DZ_RSD[relevant_QSOs]
        self.MOCKID = self.MOCKID[relevant_QSOs]
        self.PLATE = self.PLATE[relevant_QSOs]
        self.MJD = self.MJD[relevant_QSOs]
        self.FIBER = self.FIBER[relevant_QSOs]

        self.GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[relevant_QSOs,:]
        if self.DENSITY_DELTA_rows is not None:
            self.DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[relevant_QSOs,:]
        self.VEL_rows = self.VEL_rows[relevant_QSOs,:]
        self.IVAR_rows = self.IVAR_rows[relevant_QSOs,:]
        if self.lya_absorber.tau_computed():
            self.lya_absorber.tau = self.lya_absorber.tau[relevant_QSOs,:]
        if self.lya_absorber.RSDs_applied:
            self.lya_absorber.tau_noRSD = self.lya_absorber.tau_noRSD[relevant_QSOs,:]
        if self.lyb_absorber:
            if self.lyb_absorber.tau_computed():
                self.lyb_absorber.tau = self.lyb_absorber.tau[relevant_QSOs,:]
            if self.lyb_absorber.RSDs_applied:
                self.lyb_absorber.tau_noRSD = self.lyb_absorber.tau_noRSD[relevant_QSOs,:]
        if self.metals:
            for metal in iter(self.metals.values()):
                if metal.tau_computed():
                    metal.tau = metal.tau[relevant_QSOs,:]
                if metal.RSDs_applied:
                    metal.tau_noRSD = metal.tau_noRSD[relevant_QSOs,:]

        #Now trim the skewers of the remaining QSOs.
        self.N_cells = last_relevant_cell - first_relevant_cell + 1

        self.GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[:,first_relevant_cell:last_relevant_cell + 1]
        if self.DENSITY_DELTA_rows is not None:
            self.DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[:,first_relevant_cell:last_relevant_cell + 1]
        self.VEL_rows = self.VEL_rows[:,first_relevant_cell:last_relevant_cell + 1]
        self.IVAR_rows = self.IVAR_rows[:,first_relevant_cell:last_relevant_cell + 1]
        if self.lya_absorber.tau_computed():
            self.lya_absorber.tau = self.lya_absorber.tau[:,first_relevant_cell:last_relevant_cell + 1]
        if self.lya_absorber.RSDs_applied:
            self.lya_absorber.tau_noRSD = self.lya_absorber.tau_noRSD[:,first_relevant_cell:last_relevant_cell + 1]
        if self.lyb_absorber:
            if self.lyb_absorber.tau_computed():
                self.lyb_absorber.tau = self.lyb_absorber.tau[:,first_relevant_cell:last_relevant_cell + 1]
            if self.lyb_absorber.RSDs_applied:
                self.lyb_absorber.tau_noRSD = self.lyb_absorber.tau_noRSD[:,first_relevant_cell:last_relevant_cell + 1]
        if self.metals:
            for metal in iter(self.metals.values()):
                if metal.tau_computed():
                    metal.tau = metal.tau[:,first_relevant_cell:last_relevant_cell + 1]
                if metal.RSDs_applied:
                    metal.tau_noRSD = metal.tau_noRSD[:,first_relevant_cell:last_relevant_cell + 1]

        self.R = self.R[first_relevant_cell:last_relevant_cell + 1]
        self.Z = self.Z[first_relevant_cell:last_relevant_cell + 1]
        self.D = self.D[first_relevant_cell:last_relevant_cell + 1]
        self.V = self.V[first_relevant_cell:last_relevant_cell + 1]
        self.LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:last_relevant_cell + 1]
        if not isinstance(self.SIGMA_G,float):
            self.SIGMA_G = self.SIGMA_G[first_relevant_cell:last_relevant_cell + 1]
        if hasattr(self,'sample_SIGMA_G'):
            self.sample_SIGMA_G = self.sample_SIGMA_G[first_relevant_cell:last_relevant_cell + 1]

        #Modify the RSD weights to remove QSOs and cut off cells simultaneously
        if self.RSD_weights:
            trimmed_RSD_weights = {}
            k = 0
            for i in self.RSD_weights.keys():
                if relevant_QSOs[i]:
                    #Extract the matrix from the old dictionary.
                    weights = self.RSD_weights[i]

                    #Trim in both dimensions.
                    weights = weights[first_relevant_cell:last_relevant_cell + 1,:]
                    weights = weights[:,first_relevant_cell:last_relevant_cell + 1]

                    #Add the new weights to a new dictionary.
                    trimmed_RSD_weights[k] = weights
                    k += 1

            #Add the new dictionary to the object.
            self.RSD_weights = trimmed_RSD_weights

        #Remove DLAs that are no longer relevant, either because their QSO has
        #been removed, or they are outside the wavelength range.
        if self.DLA_table is not None:
            DLA_lambdas = lya*(1+self.DLA_table['Z_DLA_NO_RSD'])
            relevant_DLAs = [id for id in range(self.DLA_table['MOCKID'].shape[0]) if self.DLA_table['MOCKID'][id] in self.MOCKID and DLA_lambdas[id]>actual_lambda_min and DLA_lambdas[id]<actual_lambda_max]
            self.DLA_table = self.DLA_table[relevant_DLAs]

        return

    # Function to add transformation to the object formally.
    def add_transformation(self,transformation):
        self.transformation = transformation
        return

    # Function to scale the velocities in the skewers up.
    def scale_velocities(self,use_transformation=True,a_v=None):

        if use_transformation and a_v != None:
            print('WARNING: asked to use transformation and specified a_v.\n -> Using transformation.')
            try:
                self.VEL_rows *= self.transformation.a_v
            except AttributeError:
                print('WARNING: No tranformation found.\n -> Using specified value of a_v.')
                self.VEL_rows *= a_v

        elif use_transformation:
            try:
                self.VEL_rows *= self.transformation.a_v
            except AttributeError:
                print('WARNING: No tranformation found.\n -> Using a_v=1.')

        elif a_v == None:
            print('WARNING: Not asked to use tranformation and no value of a_v specified.\n -> Using a_v=1.')

        return

    #Function to add small scale gaussian fluctuations.
    def add_small_scale_gaussian_fluctuations(self,cell_size,generator,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya,use_transformation=True,n=None,k1=None,A0=None,R_kms=None):

        if use_transformation and n != None:
            print('WARNING: asked to use transformation and specified parameter values.\n -> Using transformation.')
            try:
                n = self.transformation.n
                k1 = self.transformation.k1
                R_kms = self.transformation.R_kms
                A0 = 58.6
            except AttributeError:
                print('WARNING: No tranformation found.\n -> Using specified values of parameters.')
                self.VEL_rows *= a_v

        elif use_transformation:
            try:
                n = self.transformation.n
                k1 = self.transformation.k1
                R_kms = self.transformation.R_kms
                A0 = 58.6
            except AttributeError:
                print('WARNING: No tranformation found.\n -> Using default parameter values.')
                n = 0.7
                k1 = 0.001
                R_kms = 25.0
                A0 = 58.6

        elif n == None:
            print('WARNING: Not asked to use tranformation and no parameter values specified.\n -> Using default parameter values.')
            n = 0.7
            k1 = 0.001
            R_kms = 25.0
            A0 = 58.6

        #Define the new R grid. Ensure that we include the entire range of R.
        old_R = self.R
        old_R_edges = utils.get_edges(old_R)
        R_edge_min = old_R_edges[0]
        R_edge_max = old_R_edges[-1]
        new_N_cells = (R_edge_max - R_edge_min) // cell_size + 1
        R_edge_max = R_edge_min + cell_size * new_N_cells
        new_R_edges = np.linspace(R_edge_min,R_edge_max,new_N_cells+1)
        new_R = utils.get_centres(new_R_edges)
        new_N_cells = new_R.shape[0]

        #Get the nearest grid points.
        NGPs = interp1d(old_R,list(range(self.N_cells)),kind='nearest',fill_value=(0,self.N_cells-1),bounds_error=False)(new_R).astype('int')

        #Redefine the necessary variables (N_cells, Z, D etc)
        self.N_cells = new_N_cells
        self.R = new_R

        #Compute the new Z, D and V by interpolating.
        old_Z = self.Z
        old_D = self.D
        old_V = self.V
        old_Z_edges = interp1d(old_R,old_Z,fill_value='extrapolate')(old_R_edges)
        old_D_edges = interp1d(old_R,old_D,fill_value='extrapolate')(old_R_edges)
        old_V_edges = interp1d(old_R,old_V,fill_value='extrapolate')(old_R_edges)
        new_Z_edges = interp1d(old_R_edges,old_Z_edges,fill_value='extrapolate')(new_R_edges)
        new_D_edges = interp1d(old_R_edges,old_D_edges,fill_value='extrapolate')(new_R_edges)
        new_V_edges = interp1d(old_R_edges,old_V_edges,fill_value='extrapolate')(new_R_edges)
        self.Z = interp1d(new_R_edges,new_Z_edges,fill_value='extrapolate')(new_R)
        self.D = interp1d(new_R_edges,new_D_edges,fill_value='extrapolate')(new_R)
        self.V = interp1d(new_R_edges,new_V_edges,fill_value='extrapolate')(new_R)
        self.LOGLAM_MAP = np.log10(lya*(1+self.Z))

        #Expand the skewers.
        expanded_GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[:,NGPs]
        expanded_VEL_rows = self.VEL_rows[:,NGPs]
        expanded_IVAR_rows = self.IVAR_rows[:,NGPs]

        #Get the extra sigma_G values from the transformation.
        extra_sigma_G = self.transformation.get_seps(self.Z)

        #Generate extra variance, either white noise or correlated.
        dkms_dhMpc = utils.get_dkms_dhMpc(0.)
        dv_kms = cell_size * dkms_dhMpc
        extra_var = independent.get_gaussian_fields(generator,self.N_cells,dv_kms=dv_kms,N_skewers=self.N_qso,white_noise=white_noise,n=n,k1=k1,A0=A0,R_kms=R_kms,norm=True)
        extra_var *= extra_sigma_G

        #Add the extra fluctuations to the expanded rows.
        expanded_GAUSSIAN_DELTA_rows += extra_var

        #Mask beyond the QSOs.
        new_Z_ledges = new_Z_edges[:-1]
        LOGLAM_ledges = np.log10(lya*(1+new_Z_ledges))
        lya_lr_mask = utils.make_IVAR_rows(lya,self.Z_QSO,LOGLAM_ledges)
        expanded_GAUSSIAN_DELTA_rows *= lya_lr_mask
        expanded_VEL_rows *= lya_lr_mask
        expanded_IVAR_rows *= lya_lr_mask

        #Assign the new skewers.
        self.GAUSSIAN_DELTA_rows = expanded_GAUSSIAN_DELTA_rows
        self.VEL_rows = expanded_VEL_rows
        self.IVAR_rows = expanded_IVAR_rows
        self.SIGMA_G = np.sqrt(extra_sigma_G**2 + (self.SIGMA_G)**2)

        return

    def return_cosmology(self):

        dtype = [('R', 'f4'), ('Z', 'f4'), ('D', 'f4'), ('V', 'f4')]
        cosmology = np.array(list(zip(self.R,self.Z,self.D,self.V)),dtype=dtype)

        return cosmology

    #Function to add physical skewers to the object via a lognormal transformation.
    def compute_physical_skewers(self,density_type='lognormal'):

        self.DENSITY_DELTA_rows = convert.gaussian_to_lognormal_delta(self.GAUSSIAN_DELTA_rows,self.SIGMA_G,self.D)

        return

    #Function to add tau skewers to an absorber using FGPA.
    def compute_tau_skewers(self,absorber):

        tau0 = self.transformation.get_tau0(self.Z)
        texp = self.transformation.get_texp(self.Z)

        # scale optical depth for this particular absorber (=1 for Lya)
        absorber_tau0 = tau0*absorber.flux_transform_m
        absorber.tau = convert.density_to_tau(self.DENSITY_DELTA_rows+1,absorber_tau0,texp)

        #Set tau to 0 beyond the quasars.
        R_edges = utils.get_edges(self.R)
        Z_edges = interp1d(self.R,self.Z,fill_value='extrapolate')(R_edges)
        LOGLAM_edges = np.log10(lya*(1+Z_edges))
        lya_lr_mask = utils.make_IVAR_rows(lya,self.Z_QSO,LOGLAM_edges[:-1])
        absorber.tau *= lya_lr_mask

        return

    #Function to compute tau for all absorbers.
    def compute_all_tau_skewers(self):

        # for each absorber, compute its optical depth skewers
        self.compute_tau_skewers(self.lya_absorber)

        # optical depth for Ly_b
        if self.lyb_absorber is not None:
            self.compute_tau_skewers(self.lyb_absorber)

        # loop over metals in dictionary
        if self.metals is not None:
            for metal in iter(self.metals.values()):
                self.compute_tau_skewers(metal)

        return

    #Get the weights for going into redshift space.
    def compute_RSD_weights(self,thermal=False,d=0.0,z_r0=2.5):

        density = 1. + self.DENSITY_DELTA_rows
        RSD_weights = RSD.get_weights(density,self.VEL_rows,self.Z,self.R,self.Z_QSO,thermal=thermal,d=d,z_r0=z_r0)
        self.RSD_weights = RSD_weights

        return

    #Get the weights dictionary required to make measurements of b_eta.
    def get_bias_eta_RSD_weights(self,z_values,d=0.,z_width=0.2,thermal=False,lambda_buffer=None):

        bias_eta_weights = bias.get_bias_eta_weights(self,z_values,d=d,z_width=z_width,include_thermal_effects=thermal,lambda_buffer=lambda_buffer)

        return bias_eta_weights

    #Function to add RSDs from the velocity skewers, with an option to include thermal effects too.
    def add_RSDs(self,absorber,weights=None,thermal=False,d=0.0,z_r0=2.5):

        density = 1 + self.DENSITY_DELTA_rows
        if weights:
            self.RSD_weights = weights
        elif not self.RSD_weights:
            self.compute_RSD_weights(thermal=thermal,d=d,z_r0=z_r0)

        new_tau = RSD.add_skewer_RSDs(absorber.tau,density,self.VEL_rows,self.Z,self.R,self.Z_QSO,thermal=thermal,weights=self.RSD_weights,d=d,z_r0=z_r0)
        tau_noRSD = absorber.tau

        #Overwrite the tau skewers and set a flag to True.
        absorber.tau = new_tau
        absorber.RSDs_applied = True
        absorber.tau_noRSD = tau_noRSD

        return

    #Function to add RSDs for all absorbers.
    def add_all_RSDs(self,weights=None,thermal=False,d=0.0,z_r0=2.5):

        # for each absorber, add RSDs
        self.add_RSDs(self.lya_absorber,weights=weights,thermal=thermal,d=d,z_r0=z_r0)

        # RSD for Ly-b
        if self.lyb_absorber is not None:
            self.add_RSDs(self.lyb_absorber,weights=weights,thermal=thermal,d=d,z_r0=z_r0)

        # loop over metals in dictionary
        if self.metals is not None:
            for metal in iter(self.metals.values()):
                self.add_RSDs(metal,thermal=thermal,weights=weights,d=d,z_r0=z_r0)

        return

    #Function to measure mean flux.
    def get_mean_quantity(self,quantity,z_value=None,z_width=None,single_value=True,power=1,all_absorbers=True):

        if quantity == 'gaussian':
            skewer_rows = self.GAUSSIAN_DELTA_rows ** power
        elif quantity == 'density':
            skewer_rows = (self.DENSITY_DELTA_rows + 1) ** power
        elif quantity == 'tau':
            skewer_rows = self.lya_absorber.tau
            if all_absorbers:
                if self.lyb_absorber is not None:
                    #Shift the skewers according to absorber rest wavelength.
                    lyb_lam = (10**self.LOGLAM_MAP)*self.lyb_absorber.rest_wave/lya
                    lyb_skewers = interp1d(lyb_lam,self.lyb_absorber.tau,axis=1,fill_value=(0.,0.),bounds_error=False)(10**self.LOGLAM_MAP)
                    #Add tau contribution and rest wavelength to header.
                    skewer_rows += lyb_skewers
                if self.metals is not None:
                    for metal in iter(self.metals.values()):
                        #Shift the skewers according to absorber rest wavelength.
                        metal_lam = (10**self.LOGLAM_MAP)*metal.rest_wave/lya
                        metal_skewers = interp1d(metal_lam,metal.tau,axis=1,fill_value=(0.,0.),bounds_error=False)(10**self.LOGLAM_MAP)
                        #Add tau contribution and rest wavelength to header.
                        skewer_rows += metal_skewers
            skewer_rows = skewer_rows ** power

        elif quantity == 'flux':
            skewer_rows = self.lya_absorber.transmission()
            if all_absorbers:
                if self.lyb_absorber is not None:
                    #Shift the skewers according to absorber rest wavelength.
                    lyb_lam = (10**self.LOGLAM_MAP)*self.lyb_absorber.rest_wave/lya
                    lyb_skewers = interp1d(lyb_lam,self.lyb_absorber.transmission(),axis=1,fill_value=(1.,1.),bounds_error=False)(10**self.LOGLAM_MAP)
                    #Add tau contribution and rest wavelength to header.
                    skewer_rows *= lyb_skewers
                if self.metals is not None:
                    for metal in iter(self.metals.values()):
                        #Shift the skewers according to absorber rest wavelength.
                        metal_lam = (10**self.LOGLAM_MAP)*metal.rest_wave/lya
                        metal_skewers = interp1d(metal_lam,metal.transmission(),axis=1,fill_value=(1.,1.),bounds_error=False)(10**self.LOGLAM_MAP)
                        #Add tau contribution and rest wavelength to header.
                        skewer_rows *= metal_skewers
            skewer_rows = skewer_rows ** power

        elif quantity == 'FlnF':
            #Use that ln(F)=-tau so FlnF = -F*tau
            # TODO: metals in this option?
            skewer_rows = (-self.lya_absorber.transmission() * self.lya_absorber.tau) ** power
        elif quantity == 'FlnFlnF':
            #Use that ln(F)=-tau so FlnFlnF = F*tau**2
            # TODO: metals in this option?
            skewer_rows = (self.lya_absorber.transmission() * (self.lya_absorber.tau)**2) ** power

        #If no z value, then compute the mean as a function of redshift.
        if not z_value:
            mean = np.zeros(self.N_cells)
            cells = np.sum(self.IVAR_rows,axis=0)>0
            mean[cells] = np.average(skewer_rows[:,cells],axis=0,weights=self.IVAR_rows[:,cells])

        #Else if there's no width, compute the mean of the cells neighbouring the z_value.
        elif not z_width:
            j_value_upper = np.searchsorted(self.Z,z_value)

            j_value_lower = max(j_value_upper - 1,0)
            relevant_rows = [i for i in range(self.N_qso) if np.sum(self.IVAR_rows[j_value_lower,j_value_upper]) == 2]

            if j_value_lower > -1:
                weight_upper = (z_value - self.Z[j_value_lower])/(self.Z[j_value_upper] - self.Z[j_value_lower])
                weight_lower = (self.Z[j_value_upper] - z_value)/(self.Z[j_value_upper] - self.Z[j_value_lower])

            else:
                weight_upper = 1
                weight_lower = 0

            weights = np.ones((self.N_qso,2))
            weights[:,0] *= weight_lower
            weights[:,1] *= weight_upper

            mean = np.average(skewer_rows[relevant_rows,j_value_lower:j_value_upper+1],weights=weights)

        #Else, compute the mean of the chunk of width z_width centred on z_value.
        else:
            j_value_upper = np.searchsorted(self.Z,z_value + z_width/2.) - 1
            j_value_lower = np.max([0,np.searchsorted(self.Z,z_value - z_width/2.)])

            if single_value:
                mean = np.average(skewer_rows[:,j_value_lower:j_value_upper+1],weights=self.IVAR_rows[:,j_value_lower:j_value_upper+1])
            else:
                mean = np.average(skewer_rows[:,j_value_lower:j_value_upper+1],weights=self.IVAR_rows[:,j_value_lower:j_value_upper+1],axis=0)

        return mean

    #Function to measure pdf.
    def get_pdf_quantity(self,quantity,z_value=None,z_width=None,single_value=True,bins=100,power=1):

        if type(bins) != float:
            N_bins = bins.shape[0]

        if quantity == 'gaussian':
            skewer_rows = self.GAUSSIAN_DELTA_rows ** power
        elif quantity == 'density':
            skewer_rows = (self.DENSITY_DELTA_rows + 1) ** power
        elif quantity == 'tau':
            skewer_rows = self.lya_absorber.tau ** power
        elif quantity == 'flux':
            skewer_rows = self.lya_absorber.transmission() ** power
        elif quantity == 'FlnF':
            #Use that ln(F)=-tau so FlnF = -F*tau
            skewer_rows = (-self.lya_absorber.transmission() * self.lya_absorber.tau) ** power
        elif quantity == 'FlnFlnF':
            #Use that ln(F)=-tau so FlnFlnF = F*tau**2
            skewer_rows = (self.lya_absorber.transmission() * (self.lya_absorber.tau)**2) ** power

        #If no z value, then compute the pdf as a function of redshift.
        if not z_value:
            hist = np.zeros(N_bins,self.N_cells)
            edges = np.zeros(N_bins+1,self.N_cells)
            for i in range(self.N_cells):
                hist_i,edges_i = np.histogram(skewer_rows[:,i],bins=bins,weights=self.IVAR_rows[:,i],density=True)
                hist[:,i] = hist_i
                edges[:,i] = edges_i

        #Else if there's no width, compute the pef of the cells neighbouring the z_value.
        elif not z_width:
            j_value_upper = np.searchsorted(self.Z,z_value)
            j_value_lower = max(j_value_upper - 1,0)
            relevant_rows = [i for i in range(self.N_qso) if np.sum(self.IVAR_rows[j_value_lower,j_value_upper]) == 2]

            if j_value_lower > -1:
                weight_upper = (z_value - self.Z[j_value_lower])/(self.Z[j_value_upper] - self.Z[j_value_lower])
                weight_lower = (self.Z[j_value_upper] - z_value)/(self.Z[j_value_upper] - self.Z[j_value_lower])

            else:
                weight_upper = 1
                weight_lower = 0

            weights = np.ones((self.N_qso,2))
            weights[:,0] *= weight_lower
            weights[:,1] *= weight_upper

            hist,edges = np.histogram(skewer_rows[relevant_rows,j_value_lower,j_value_upper+1],bins=bins,weights=weights,density=True)

        #Else, compute the mean of the chunk of width z_width centred on z_value.
        else:
            j_value_upper = np.searchsorted(self.Z,z_value + z_width/2.) - 1
            j_value_lower = np.max([0,np.searchsorted(self.Z,z_value - z_width/2.)])
            if single_value:
                hist,edges = np.histogram(skewer_rows[:,j_value_lower:j_value_upper+1],bins=bins,weights=self.IVAR_rows[:,j_value_lower:j_value_upper+1],density=True)
            else:
                hist = np.zeros(N_bins,j_value_upper+1-j_value_lower)
                edges = np.zeros(N_bins+1,j_value_upper+1-j_value_lower)
                for i in range(j_value_upper+1-j_value_lower):
                    hist_i,edges_i = np.histogram(skewer_rows[:,i],bins=bins,weights=self.IVAR_rows[:,i],density=True)
                    hist[:,i] = hist_i
                    edges[:,i] = edges_i

        return hist,edges

    #Function to measure sigma dF.
    def get_sigma_dF(self,absorber,z_value=None,z_width=None,mean_F=None):

        if not mean_F:
            mean_F = self.get_mean_quantity('flux',z_value=z_value,z_width=z_width)

        F = absorber.transmission()
        dF = F/mean_F

        if not z_value:
            # TODO: there's no weighting in here?
            sigma_dF = np.std(dF,axis=0)

        elif not z_width:
            j_value_upper = np.searchsorted(self.Z,z_value)

            j_value_lower = np.max([j_value_upper - 1,0])
            relevant_rows = [i for i in range(self.N_qso) if np.sum(self.IVAR_rows[j_value_lower,j_value_upper]) == 2]

            sigma_dF = np.std(dF[relevant_rows,j_value_lower:j_value_upper+1])

        else:
            j_value_upper = np.searchsorted(self.Z,z_value + z_width/2.)
            j_value_lower = np.max([0,np.searchsorted(self.Z,z_value - z_width/2.) - 1])

            sigma_dF = np.std(dF[:,j_value_lower:j_value_upper+1])

        #print('gsd return',sigma_dF)
        return sigma_dF

    #Method to combine data from two objects into one.
    # TODO: add something to check that we can just take values from 1 of the objects
    @classmethod
    def combine_files(cls,object_A,object_B,gaussian_only=False):

        N_qso = object_A.N_qso + object_B.N_qso

        """
        something to check N_cells is the same in both files
        """

        N_cells = object_A.N_cells
        SIGMA_G = object_A.SIGMA_G

        TYPE = np.concatenate((object_A.TYPE,object_B.TYPE),axis=0)
        RA = np.concatenate((object_A.RA,object_B.RA),axis=0)
        DEC = np.concatenate((object_A.DEC,object_B.DEC),axis=0)
        Z_QSO = np.concatenate((object_A.Z_QSO,object_B.Z_QSO),axis=0)
        DZ_RSD = np.concatenate((object_A.DZ_RSD,object_B.DZ_RSD),axis=0)
        MOCKID = np.concatenate((object_A.MOCKID,object_B.MOCKID),axis=0)
        PLATE = np.concatenate((object_A.PLATE,object_B.PLATE),axis=0)
        MJD = np.concatenate((object_A.MJD,object_B.MJD),axis=0)
        FIBER = np.concatenate((object_A.FIBER,object_B.FIBER),axis=0)

        GAUSSIAN_DELTA_rows = np.concatenate((object_A.GAUSSIAN_DELTA_rows,object_B.GAUSSIAN_DELTA_rows),axis=0)
        VEL_rows = np.concatenate((object_A.VEL_rows,object_B.VEL_rows),axis=0)
        IVAR_rows = np.concatenate((object_A.IVAR_rows,object_B.IVAR_rows),axis=0)
        if gaussian_only:
            DENSITY_DELTA_rows = None
        else:
            DENSITY_DELTA_rows = np.concatenate((object_A.DENSITY_DELTA_rows,object_B.DENSITY_DELTA_rows),axis=0)

        """
        Something to check this is ok?
        """

        Z = object_A.Z
        LOGLAM_MAP = object_A.LOGLAM_MAP
        R = object_A.R
        D = object_A.D
        V = object_A.V

        return cls(N_qso,N_cells,SIGMA_G,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    #Function to save in the colore format.
    def save_as_colore(self,quantity,filename,header,overwrite=False,cell_size=None,compress=True):
        t = time.time()

        #Organise the catalog data into a colore-format array.
        colore_1_data = []
        for i in range(self.N_qso):
            colore_1_data += [(self.TYPE[i],self.RA[i],self.DEC[i],self.Z_QSO[i],self.DZ_RSD[i],self.MOCKID[i])]
        dtype = [('TYPE', 'f4'), ('RA', 'f4'), ('DEC', 'f4'), ('Z_COSMO', 'f4'), ('DZ_RSD', 'f4'), ('MOCKID', int)]
        colore_1 = np.array(colore_1_data,dtype=dtype)

        #Choose the right skewers according to input quantity.
        if quantity == 'gaussian':
            colore_2 = self.GAUSSIAN_DELTA_rows
        elif quantity == 'density':
            colore_2 = self.DENSITY_DELTA_rows + 1
        elif quantity == 'tau':
            colore_2 == self.lya_absorber.tau
        elif quantity == 'flux':
            colore_2 = self.lya_absorber.transmission()
        colore_2 = colore_2.astype('float32')

        #Add the velocity skewers and cosmology data
        colore_3 = self.VEL_rows.astype('float32')
        colore_4_data = []
        for i in range(self.N_cells):
            colore_4_data += [(self.R[i],self.Z[i],self.D[i],self.V[i])]
        dtype = [('R', 'f4'), ('Z', 'f4'), ('D', 'f4'), ('V', 'f4')]
        colore_4 = np.array(colore_4_data,dtype=dtype)

        #Construct HDUs from the data arrays.
        prihdr = fits.Header()
        prihdu = fits.PrimaryHDU(header=prihdr)
        cols_CATALOG = fits.ColDefs(colore_1)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')
        hdu_GAUSSIAN = fits.ImageHDU(data=colore_2,header=header,name=quantity.upper())
        hdu_VEL = fits.ImageHDU(data=colore_3,header=header,name='VELOCITY')
        cols_COSMO = fits.ColDefs(colore_4)
        hdu_COSMO = fits.BinTableHDU.from_columns(cols_COSMO,header=header,name='COSMO')

        #Combine the HDUs into an HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([prihdu, hdu_CATALOG, hdu_GAUSSIAN, hdu_VEL, hdu_COSMO])
        hdulist.writeto(filename,overwrite=overwrite)
        hdulist.close()

        #print('--> saving takes {:2.3f}s'.format(time.time()-t))
        t = time.time()
        #Compress the file if desired.
        if compress:
            utils.compress_file(filename)
        #print('--> compressing takes {:2.3f}s'.format(time.time()-t))

        return

    #Function to save in the picca format.
    def save_as_picca_delta(self,quantity,filename,header,mean_data=None,overwrite=False,min_number_cells=2,cell_size=None,notnorm=False,add_QSO_RSDs=True,compress=True,all_absorbers=False):

        t = time.time()
        lya_lambdas = 10**self.LOGLAM_MAP

        #If not normalising:
        if notnorm:
            if quantity == 'gaussian':
                skewer_rows = self.GAUSSIAN_DELTA_rows
            elif quantity == 'density':
                skewer_rows = self.DENSITY_DELTA_rows
            elif quantity == 'tau':
                skewer_rows = self.lya_absorber.tau
                #If desired, add all metal absorbers.
                if all_absorbers:
                    if self.lyb_absorber is not None:
                        #Shift the skewers according to absorber rest wavelength.
                        lyb_lam = (10**self.LOGLAM_MAP)*self.lyb_absorber.rest_wave/lya
                        lyb_skewers = interp1d(lyb_lam,self.lyb_absorber.tau,axis=1,fill_value=(0.,0.),bounds_error=False)(10**self.LOGLAM_MAP)
                        #Add tau contribution and rest wavelength to header.
                        skewer_rows += lyb_skewers
                        header[self.lyb_absorber.HDU_name] = self.lyb_absorber.rest_wave
                    if self.metals is not None:
                        for metal in iter(self.metals.values()):
                            #Shift the skewers according to absorber rest wavelength.
                            metal_lam = (10**self.LOGLAM_MAP)*metal.rest_wave/lya
                            metal_skewers = interp1d(metal_lam,metal.tau,axis=1,fill_value=(0.,0.),bounds_error=False)(10**self.LOGLAM_MAP)
                            #Add tau contribution and rest wavelength to header.
                            skewer_rows += metal_skewers
                            header[metal.HDU_name] = metal.rest_wave
            elif quantity == 'flux':
                skewer_rows = self.lya_absorber.transmission()
                #If desired, add all metal absorbers.
                if all_absorbers:
                    if self.lyb_absorber is not None:
                        #Shift the skewers according to absorber rest wavelength.
                        lyb_lam = (10**self.LOGLAM_MAP)*self.lyb_absorber.rest_wave/lya
                        lyb_skewers = interp1d(lyb_lam,self.lyb_absorber.transmission(),axis=1,fill_value=(1.,1.),bounds_error=False)(10**self.LOGLAM_MAP)
                        #Add tau contribution and rest wavelength to header.
                        skewer_rows *= lyb_skewers
                        header[self.lyb_absorber.HDU_name] = self.lyb_absorber.rest_wave
                    if self.metals is not None:
                        for metal in iter(self.metals.values()):
                            #Shift the skewers according to absorber rest wavelength.
                            metal_lam = (10**self.LOGLAM_MAP)*metal.rest_wave/lya
                            metal_skewers = interp1d(metal_lam,metal.transmission(),axis=1,fill_value=(1.,1.),bounds_error=False)(10**self.LOGLAM_MAP)
                            #Add tau contribution and rest wavelength to header.
                            skewer_rows *= metal_skewers
                            header[metal.HDU_name] = metal.rest_wave

        #If normalising:
        else:
            if quantity == 'gaussian':
                skewer_rows = self.GAUSSIAN_DELTA_rows
            elif quantity == 'density':
                skewer_rows = self.DENSITY_DELTA_rows
            elif quantity == 'tau':
                skewer_rows = np.zeros(self.lya_absorber.tau.shape)
                cells = np.sum(self.IVAR_rows,axis=0)>0
                #print('min mean tau:',np.min(mean[cells]))
                skewer_rows[:,cells] = self.lya_absorber.tau[:,cells]
                #If desired, add all metal absorbers.
                if all_absorbers:
                    if self.lyb_absorber is not None:
                        #Shift the skewers according to absorber rest wavelength.
                        lyb_lam = (10**self.LOGLAM_MAP)*self.lyb_absorber.rest_wave/lya
                        lyb_skewers = interp1d(lyb_lam,self.lyb_absorber.tau,axis=1,fill_value=(0.,0.),bounds_error=False)(10**self.LOGLAM_MAP)
                        #Add tau contribution and rest wavelength to header.
                        skewer_rows += lyb_skewers
                        header[self.lyb_absorber.HDU_name] = self.lyb_absorber.rest_wave
                    if self.metals is not None:
                        for metal in iter(self.metals.values()):
                            #Shift the skewers according to absorber rest wavelength.
                            metal_lam = (10**self.LOGLAM_MAP)*metal.rest_wave/lya
                            metal_skewers = interp1d(metal_lam,metal.tau,axis=1,fill_value=(0.,0.),bounds_error=False)(10**self.LOGLAM_MAP)
                            #Add tau contribution and rest wavelength to header.
                            skewer_rows += metal_skewers
                            header[metal.HDU_name] = metal.rest_wave
                #Get mean with redshift.
                if mean_data is None:
                    mean = self.get_mean_quantity('tau',all_absorbers=all_absorbers)
                else:
                    mean = np.interp(self.Z,mean_data['z'],mean_data['mean'])
                #Normalise to get deltas
                skewer_rows[:,cells] = skewer_rows[:,cells]/mean[cells] - 1
            elif quantity == 'flux':
                skewer_rows = np.zeros(self.lya_absorber.transmission().shape)
                cells = np.sum(self.IVAR_rows,axis=0)>0
                #print('min mean flux:',np.min(mean[cells]))
                skewer_rows[:,cells] = self.lya_absorber.transmission()[:,cells]
                #If desired, add all metal absorbers.
                if all_absorbers:
                    if self.lyb_absorber is not None:
                        #Shift the skewers according to absorber rest wavelength.
                        lyb_lam = (10**self.LOGLAM_MAP)*self.lyb_absorber.rest_wave/lya
                        lyb_skewers = interp1d(lyb_lam,self.lyb_absorber.transmission(),axis=1,fill_value=(1.,1.),bounds_error=False)(10**self.LOGLAM_MAP)
                        #Add tau contribution and rest wavelength to header.
                        skewer_rows *= lyb_skewers
                        header[self.lyb_absorber.HDU_name] = self.lyb_absorber.rest_wave
                    if self.metals is not None:
                        for metal in iter(self.metals.values()):
                            #Shift the skewers according to absorber rest wavelength.
                            metal_lam = (10**self.LOGLAM_MAP)*metal.rest_wave/lya
                            metal_skewers = interp1d(metal_lam,metal.transmission(),axis=1,fill_value=(1.,1.),bounds_error=False)(10**self.LOGLAM_MAP)
                            #Add tau contribution and rest wavelength to header.
                            skewer_rows *= metal_skewers
                            header[metal.HDU_name] = metal.rest_wave
                #Get mean with redshift.
                if mean_data is None:
                    mean = self.get_mean_quantity('flux',all_absorbers=all_absorbers)
                else:
                    mean = np.interp(self.Z,mean_data['z'],mean_data['mean'])
                #Normalise to get deltas
                skewer_rows[:,cells] = skewer_rows[:,cells]/mean[cells] - 1

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        relevant_QSOs = []
        for i in range(self.N_qso):
            if np.sum(self.IVAR_rows[i,:]) >= min_number_cells:
                relevant_QSOs += [i]
        #relevant_QSOs = np.sum(self.IVAR_rows,axis=1) >= min_number_cells
        non_rel_QSOs = np.sum(self.IVAR_rows,axis=1) < min_number_cells

        #Trim data according to the relevant cells and QSOs.
        relevant_skewer_rows = skewer_rows[relevant_QSOs,:]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,:]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[:]

        #If desired, add in QSO RSDs.
        if add_QSO_RSDs:
            Z_QSO = self.Z_QSO + self.DZ_RSD
        else:
            Z_QSO = self.Z_QSO

        #Organise the data into picca-format arrays.
        picca_0 = relevant_skewer_rows.T.astype('float32')
        picca_1 = relevant_IVAR_rows.T
        picca_2 = relevant_LOGLAM_MAP.astype('float32')
        picca_3_data = list(zip(self.RA[relevant_QSOs],self.DEC[relevant_QSOs],Z_QSO[relevant_QSOs],self.PLATE[relevant_QSOs],self.MJD[relevant_QSOs],self.FIBER[relevant_QSOs],self.MOCKID[relevant_QSOs]))
        dtype = [('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('PLATE', int), ('MJD', 'f4'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Add cell size to the header (average)
        dr_hMpc = (self.R[-1] - self.R[0])/(self.N_cells - 1)
        header['dr_hMpc'] = dr_hMpc

        #Make the data into suitable HDUs.
        hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        hdu_CATALOG = fits.BinTableHDU(picca_3,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        #print('--> organising takes {:2.3f}s'.format(time.time()-t))
        t = time.time()
        hdulist.writeto(filename,overwrite=overwrite)
        #print('--> writing takes {:2.3f}s'.format(time.time()-t))
        t = time.time()
        hdulist.close()

        #Compress the file if desired.
        if compress:
            utils.compress_file(filename)
        #print('--> compressing takes {:2.3f}s'.format(time.time()-t))

        return

    #Compute transmission for a particular absorber, on a particular grid
    def compute_grid_transmission(self,absorber,wave_grid):

        #Get data from the absorber.
        F_skewer = absorber.transmission()
        rest_wave = absorber.rest_wave
        wave_skewer = rest_wave*(1+self.Z)

        #Create the F_grid.
        N_los = F_skewer.shape[0]
        N_w = wave_grid.shape[0]
        F_grid = np.empty([N_los,N_w])

        #Interpolate the skewers.
        for i in range(N_los):
            F_grid[i,] = np.interp(wave_grid,wave_skewer,F_skewer[i],left=1.0,right=1.0)

        return F_grid

    #Function to save data as a transmission file.
    def save_as_transmission(self,filename,header,overwrite=False,wave_min=3550.,wave_max=6500.,wave_step=0.2,fmt='final',add_QSO_RSDs=True,compress=True):

        t = time.time()

        # define common wavelength grid to be written in files (in Angstroms)
        wave_grid = np.arange(wave_min,wave_max,wave_step).astype('float32')

        # compute Lyman alpha transmission on grid of wavelengths
        F_grid_Lya = self.compute_grid_transmission(self.lya_absorber,wave_grid).astype('float32')

        #construct quasar catalog HDU
        if add_QSO_RSDs:
            Z_QSO = self.Z_QSO + self.DZ_RSD
        else:
            Z_QSO = self.Z_QSO
        catalog_data = list(zip(self.RA,self.DEC,Z_QSO,self.Z_QSO,self.MOCKID))
        dtype = [('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_noRSD', 'f4'), ('MOCKID', int)]
        catalog_data = np.array(catalog_data,dtype=dtype)

        #Construct HDUs from the data arrays.
        prihdr = fits.Header()
        prihdu = fits.PrimaryHDU(header=prihdr)
        cols_METADATA = fits.ColDefs(catalog_data)
        hdu_METADATA = fits.BinTableHDU.from_columns(cols_METADATA,header=header,name='METADATA')
        hdu_WAVELENGTH = fits.ImageHDU(data=wave_grid,header=header,name='WAVELENGTH')

        #Combine the HDUs into an HDUlist (including DLAs and metals, if they have been computed)
        list_hdu = [prihdu, hdu_METADATA, hdu_WAVELENGTH]

        #Set up the absorber HDUs according to the input format 'fmt'.
        if fmt=='single_HDU':

            #Transmission of all absorbers.
            abs_header = header.copy()
            abs_header['LYA'] = self.lya_absorber.rest_wave
            F_grid = F_grid_Lya
            if self.lyb_absorber is not None:
                abs_header[self.lyb_absorber.HDU_name] = self.lyb_absorber.rest_wave
                F_grid *= self.compute_grid_transmission(self.lyb_absorber,wave_grid).astype('float32')
            if self.metals is not None:
                for metal in iter(self.metals.values()):
                    abs_header[metal.HDU_name] = metal.rest_wave
                    F_grid *= self.compute_grid_transmission(metal,wave_grid).astype('float32')
            hdu_F = fits.ImageHDU(data=F_grid,header=abs_header,name='F')
            list_hdu += [hdu_F]

        elif fmt == 'final':

            #Gives transmission of Lya only
            lya_header = header
            lya_header['LYA'] = self.lya_absorber.rest_wave
            list_hdu += [fits.ImageHDU(data=F_grid_Lya,header=lya_header,name='F_LYA')]

            # compute Lyman beta transmission on grid of wavelengths
            if self.lyb_absorber is not None:
                F_grid_Lyb = self.compute_grid_transmission(self.lyb_absorber,wave_grid).astype('float32')
                HDU_name = 'F_'+self.lyb_absorber.HDU_name
                lyb_header = header.copy()
                lyb_header[self.lyb_absorber.HDU_name] = self.lyb_absorber.rest_wave
                list_hdu += [fits.ImageHDU(data=F_grid_Lyb,header=lyb_header,name=HDU_name)]

            #Add an HDU for each metal computed.
            if self.metals is not None:
                F_grid_all_metals = np.ones_like(F_grid_Lya)
                met_header = header.copy()
                # compute metals' transmission on grid of wavelengths
                for metal in iter(self.metals.values()):
                    F_grid_all_metals *= self.compute_grid_transmission(metal,wave_grid).astype('float32')
                    met_header[metal.HDU_name] = metal.rest_wave
                HDU_name = 'F_METALS'
                list_hdu += [fits.ImageHDU(data=F_grid_all_metals,header=met_header,name=HDU_name)]

        elif fmt == 'develop':

            #Gives transmission of Lya only
            lya_header = header
            lya_header['LYA'] = self.lya_absorber.rest_wave
            list_hdu += [fits.ImageHDU(data=F_grid_Lya,header=lya_header,name='F_LYA')]

            # compute Lyman beta transmission on grid of wavelengths
            if self.lyb_absorber is not None:
                F_grid_Lyb = self.compute_grid_transmission(self.lyb_absorber,wave_grid).astype('float32')
                HDU_name = 'F_'+self.lyb_absorber.HDU_name
                lyb_header = header.copy()
                lyb_header[self.lyb_absorber.HDU_name] = self.lyb_absorber.rest_wave
                list_hdu += [fits.ImageHDU(data=F_grid_Lyb,header=lyb_header,name=HDU_name)]

            #Add an HDU for each metal computed.
            if self.metals is not None:
                # compute metals' transmission on grid of wavelengths
                for metal in iter(self.metals.values()):
                    F_grid_metal = self.compute_grid_transmission(metal,wave_grid).astype('float32')
                    HDU_name = 'F_'+metal.HDU_name
                    met_header = header.copy()
                    met_header[metal.HDU_name] = metal.rest_wave
                    list_hdu += [fits.ImageHDU(data=F_grid_metal,header=met_header,name=HDU_name)]

        # add table of DLAs
        if self.DLA_table is not None:
            hdu_DLAs = fits.hdu.BinTableHDU(self.DLA_table,header=header,name='DLA')
            list_hdu.append(hdu_DLAs)

        #Save as a new file. Close the HDUlist.
        hdulist = fits.HDUList(list_hdu)
        hdulist.writeto(filename,overwrite=overwrite)
        hdulist.close()
        #print('--> saving takes {:2.3f}s'.format(time.time()-t))
        t = time.time()

        #Compress the file if desired.
        if compress:
            utils.compress_file(filename)
        #print('--> compressing takes {:2.3f}s'.format(time.time()-t))

        return

    #Function to calculate the mean and variance of the different quantities as a function of Z.
    def get_means(self,all_absorbers=True):

        #For each cell, determine the number of skewers for which it is relevant.
        N_skewers = np.sum(self.IVAR_rows,axis=0)

        #Calculate the mean in each cell of the gaussian delta and its square.
        GM = self.get_mean_quantity('gaussian',all_absorbers=all_absorbers)
        GSM = self.get_mean_quantity('gaussian',power=2.,all_absorbers=all_absorbers)

        #Calculate the mean in each cell of the density delta and its square.
        DM = self.get_mean_quantity('density',all_absorbers=all_absorbers)
        DSM = self.get_mean_quantity('density',power=2.,all_absorbers=all_absorbers)

        #Calculate the mean in each cell of the tau and its square.
        TM = self.get_mean_quantity('tau',all_absorbers=all_absorbers)
        TSM = self.get_mean_quantity('tau',power=2.,all_absorbers=all_absorbers)

        #Calculate the mean in each cell of the flux and its square.
        FM = self.get_mean_quantity('flux',all_absorbers=all_absorbers)
        FSM = self.get_mean_quantity('flux',power=2.,all_absorbers=all_absorbers)

        #Calculate the mean in each cell of the flux delta and its square.
        #FDB = np.average(relevant_delta_F,weights=relevant_IVAR+small,axis=0)*relevant_cells
        #FDSB = np.average(relevant_delta_F**2,weights=relevant_IVAR+small,axis=0)*relevant_cells

        #Stitch together the means into a binary table.
        dtype = [('Z', 'f4'), ('N', 'int'), ('GAUSSIAN', 'f4'), ('GAUSSIAN_SQUARED', 'f4'), ('DENSITY', 'f4'), ('DENSITY_SQUARED', 'f4'), ('TAU', 'f4'), ('TAU_SQUARED', 'f4'), ('F', 'f4'), ('F_SQUARED', 'f4')]
        #, ('F_DELTA', 'f4'), ('F_DELTA_SQUARED', 'f4')]
        means = np.array(list(zip(self.Z,N_skewers,GM,GSM,DM,DSM,TM,TSM,FM,FSM)),dtype=dtype)
        #,FDB,FDSB

        return means

    #Function to save the means as a function of z.
    def save_statistics(self,filepath,overwrite=False,compress=True,all_absorbers=True):

        means = self.get_means(all_absorbers=all_absorbers)
        statistics = stats.means_to_statistics(means)
        stats.write_statistics(filepath,statistics,overwrite=overwrite,compress=compress)

        return statistics

    #Function to add DLAs to a set of skewers.
    def add_DLA_table(self,seed,dla_bias=2.0,evol='b_const',method='global'):

        #If extrapolate_z_down is set to a value below the skewer, then we extrapolate down to that value.
        #Otherwise, we start placing DLAs at the start of the skewer.
        extrapolate_z_down = None
        DLA_table = DLA.get_DLA_table(self,dla_bias=dla_bias,extrapolate_z_down=extrapolate_z_down,seed=seed,method=method)
        self.DLA_table = DLA_table

        return


    ####
    """
    Obsolete functions

    #Function to save data as a Gaussian colore file.
    def save_as_gaussian_colore(self,filename,header,overwrite=False):

        #Organise the data into colore-format arrays.
        colore_1_data = []
        for i in range(self.N_qso):
            colore_1_data += [(self.TYPE[i],self.RA[i],self.DEC[i],self.Z_QSO[i],self.DZ_RSD[i],self.MOCKID[i])]

        dtype = [('TYPE', 'f8'), ('RA', 'f8'), ('DEC', 'f8'), ('Z_COSMO', 'f8'), ('DZ_RSD', 'f8'), ('MOCKID', int)]
        colore_1 = np.array(colore_1_data,dtype=dtype)
        colore_2 = self.GAUSSIAN_DELTA_rows
        colore_3 = self.VEL_rows

        colore_4_data = []
        for i in range(self.N_cells):
            colore_4_data += [(self.R[i],self.Z[i],self.D[i],self.V[i])]

        dtype = [('R', 'f8'), ('Z', 'f8'), ('D', 'f8'), ('V', 'f8')]
        colore_4 = np.array(colore_4_data,dtype=dtype)

        #Construct HDUs from the data arrays.
        prihdr = fits.Header()
        prihdu = fits.PrimaryHDU(header=prihdr)
        cols_CATALOG = fits.ColDefs(colore_1)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')
        hdu_GAUSSIAN = fits.ImageHDU(data=colore_2,header=header,name='GAUSSIAN_DELTA')
        hdu_VEL = fits.ImageHDU(data=colore_3,header=header,name='VELOCITY')
        cols_COSMO = fits.ColDefs(colore_4)
        hdu_COSMO = fits.BinTableHDU.from_columns(cols_COSMO,header=header,name='COSMO')

        #Combine the HDUs into an HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([prihdu, hdu_CATALOG, hdu_GAUSSIAN, hdu_VEL, hdu_COSMO])
        hdulist.writeto(filename,overwrite=overwrite)
        hdulist.close

        return

    #Function to save data as a picca density file.
    def save_as_picca_gaussian(self,filename,header,overwrite=False,zero_mean_delta=False,min_number_cells=2,mean_DELTA=None):

        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        relevant_QSOs = []
        for i in range(self.N_qso):
            if np.sum(self.IVAR_rows[i,:]) >= min_number_cells:
                relevant_QSOs += [i]

        #Trim data according to the relevant cells and QSOs.
        relevant_GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[relevant_QSOs,:]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,:]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[:]

        #If desired, enforce that the Delta rows have zero mean.
        if zero_mean_delta == True:
            relevant_GAUSSIAN_DELTA_rows = utils.normalise_deltas(relevant_GAUSSIAN_DELTA_rows,mean_DELTA)

        #Organise the data into picca-format arrays.
        picca_0 = relevant_GAUSSIAN_DELTA_rows.T
        picca_1 = relevant_IVAR_rows.T
        picca_2 = relevant_LOGLAM_MAP

        picca_3_data = []
        for i in range(self.N_qso):
            if i in relevant_QSOs:
                picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

        dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('PLATE', int), ('MJD', 'f8'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Make the data into suitable HDUs.
        hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        cols_CATALOG = fits.ColDefs(picca_3)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        hdulist.writeto(filename,overwrite=overwrite)
        hdulist.close()

        return

    #Function to save data as a Lognormal colore file.
    def save_as_physical_colore(self,filename,header):

        #Organise the data into colore-format arrays.
        colore_1_data = []
        for i in range(self.N_qso):
            colore_1_data += [(self.TYPE[i],self.RA[i],self.DEC[i],self.Z_QSO[i],self.DZ_RSD[i],self.MOCKID[i])]

        dtype = [('TYPE', 'f8'), ('RA', 'f8'), ('DEC', 'f8'), ('Z_COSMO', 'f8'), ('DZ_RSD', 'f8'), ('MOCKID', int)]
        colore_1 = np.array(colore_1_data,dtype=dtype)

        colore_2 = self.DENSITY_DELTA_rows
        colore_3 = self.VEL_rows

        colore_4_data = []
        for i in range(self.N_cells):
            colore_4_data += [(self.R[i],self.Z[i],self.D[i],self.V[i])]

        dtype = [('R', 'f8'), ('Z', 'f8'), ('D', 'f8'), ('V', 'f8')]
        colore_4 = np.array(colore_4_data,dtype=dtype)

        #Construct HDUs from the data arrays.
        prihdr = fits.Header()
        prihdu = fits.PrimaryHDU(header=prihdr)
        cols_CATALOG = fits.ColDefs(colore_1)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')
        hdu_DELTA = fits.ImageHDU(data=colore_2,header=header,name='PHYSICAL_DELTA')
        hdu_VEL = fits.ImageHDU(data=colore_3,header=header,name='VELOCITY')
        cols_COSMO = fits.ColDefs(colore_4)
        hdu_COSMO = fits.BinTableHDU.from_columns(cols_COSMO,header=header,name='COSMO')

        #Combine the HDUs into an HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([prihdu, hdu_CATALOG, hdu_DELTA, hdu_VEL, hdu_COSMO])
        hdulist.writeto(filename)
        hdulist.close

        return

    #Function to save data as a picca density file.
    def save_as_picca_density(self,filename,header,zero_mean_delta=False,min_number_cells=2,mean_DELTA=None):

        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        relevant_QSOs = []
        for i in range(self.N_qso):
            if np.sum(self.IVAR_rows[i,:]) >= min_number_cells:
                relevant_QSOs += [i]

        #Trim data according to the relevant cells and QSOs.
        relevant_DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[relevant_QSOs,:]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,:]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[:]

        #If desired, enforce that the Delta rows have zero mean.
        if zero_mean_delta == True:
            relevant_DENSITY_DELTA_rows = utils.normalise_deltas(relevant_DENSITY_DELTA_rows,mean_DELTA)

        #Organise the data into picca-format arrays.
        picca_0 = relevant_DENSITY_DELTA_rows.T
        picca_1 = relevant_IVAR_rows.T
        picca_2 = relevant_LOGLAM_MAP

        picca_3_data = []
        for i in range(self.N_qso):
            if i in relevant_QSOs:
                picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

        dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('PLATE', int), ('MJD', 'f8'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Make the data into suitable HDUs.
        hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        cols_CATALOG = fits.ColDefs(picca_3)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        hdulist.writeto(filename)
        hdulist.close()

        return

    #Function to save data as a picca density file.
    def save_as_picca_tau(self,absorber,filename,header,overwrite=False,min_number_cells=2):

        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        relevant_QSOs = []
        for i in range(self.N_qso):
            if np.sum(self.IVAR_rows[i,:]) >= min_number_cells:
                relevant_QSOs += [i]

        #Trim data according to the relevant cells and QSOs.
        relevant_TAU_rows = absorber.tau[relevant_QSOs,:]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,:]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[:]

        #Organise the data into picca-format arrays.
        picca_0 = relevant_TAU_rows.T
        picca_1 = relevant_IVAR_rows.T
        picca_2 = relevant_LOGLAM_MAP

        picca_3_data = []
        for i in range(self.N_qso):
            if i in relevant_QSOs:
                picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

        dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('PLATE', int), ('MJD', 'f8'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Make the data into suitable HDUs.
        hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        cols_CATALOG = fits.ColDefs(picca_3)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        hdulist.writeto(filename,overwrite=overwrite)
        hdulist.close()

        return

    #Function to save data as a picca flux file.
    def save_as_picca_flux(self,filename,header,min_number_cells=2,mean_F_data=None,rebin_size_hMpc=None):

        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        relevant_QSOs = []
        for i in range(self.N_qso):
            if np.sum(self.IVAR_rows[i,:]) >= min_number_cells:
                relevant_QSOs += [i]

        # get Lya transmission
        F = self.lya_absorber.transmission()

        #Trim data according to the relevant cells and QSOs.
        relevant_F = F[relevant_QSOs,:]
        relevant_IVAR = self.IVAR_rows[relevant_QSOs,:]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[:]
        relevant_Z = self.Z[:]

        #Calculate mean F as a function of z for the relevant cells, then delta_F.
        try:
            mean_F_z_values = mean_F_data[:,0]
            mean_F = mean_F_data[:,1]
            relevant_mean_F = np.interp(relevant_Z,mean_F_z_values,mean_F)
        except ValueError:
            #This is done with a 'hack' to avoid problems with weights summing to zero.
            small = 1.0e-10
            relevant_mean_F = np.average(relevant_F,weights=relevant_IVAR+small,axis=0)

        relevant_delta_F = ((relevant_F)/relevant_mean_F - 1)*relevant_IVAR

        if rebin_size_hMpc:
            grid_start = self.R[0]
            rebin_map = np.array((self.R - self.R[0])//rebin_size_hMpc,dtype='int')
            new_relevant_delta_F = np.zeros((self.N_qso,max(rebin_map)))
            new_relevant_IVAR = np.zeros((self.N_qso,max(rebin_map)))
            new_Z = np.zeros(max(rebin_map))
            for j in range(max(rebin_map)):
                j_lo = np.argmax(rebin_map==j)
                j_hi = np.argmin(rebin_map==j)
                new_relevant_delta_F[:,j] = np.average(relevant_delta_F[:,j_lo:j_hi],axis=1)
                new_relevant_IVAR[:,j] = (relevant_IVAR[:,j_lo:j_hi] == j_hi - j_lo)
                new_Z[j] = np.average(self.Z[j_lo:j_hi])
            new_relevant_LOGLAM_MAP = np.log(lya*(1+new_Z))

            picca_0 = new_relevant_delta_F.T
            picca_1 = new_relevant_IVAR.T
            picca_2 = new_relevant_LOGLAM_MAP

        else:
            #Organise the data into picca-format arrays.
            picca_0 = relevant_delta_F.T
            picca_1 = relevant_IVAR.T
            picca_2 = relevant_LOGLAM_MAP

        picca_3_data = []
        for i in range(self.N_qso):
            if i in relevant_QSOs:
                picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

        dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('PLATE', int), ('MJD', 'f8'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Make the data into suitable HDUs.
        hdu_F = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        cols_CATALOG = fits.ColDefs(picca_3)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_F, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        hdulist.writeto(filename)
        hdulist.close()

        return

    # TODO: Do we really want this? Does it make much sense?
    #Function to save data as a picca velocity file.
    def save_as_picca_velocity(self,filename,header,zero_mean_delta=False,min_number_cells=2,overwrite=False):

        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        relevant_QSOs = []
        for i in range(self.N_qso):
            if np.sum(self.IVAR_rows[i,:]) >= min_number_cells:
                relevant_QSOs += [i]

        #Trim data according to the relevant cells and QSOs.
        relevant_VEL = self.VEL_rows[relevant_QSOs,:]
        relevant_IVAR = self.IVAR_rows[relevant_QSOs,:]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[:]

        #Organise the data into picca-format arrays.
        picca_0 = relevant_VEL.T
        picca_1 = relevant_IVAR.T
        picca_2 = relevant_LOGLAM_MAP

        picca_3_data = []
        for i in range(self.N_qso):
            if i in relevant_QSOs:
                picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

        dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('PLATE', int), ('MJD', 'f8'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Make the data into suitable HDUs.
        hdu_VEL = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        cols_CATALOG = fits.ColDefs(picca_3)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_VEL, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        hdulist.writeto(filename,overwrite=overwrite)
        hdulist.close()

        return

    """
