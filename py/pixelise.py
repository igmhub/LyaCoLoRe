import numpy as np
from astropy.io import fits
import time

import general
import input
import convert
import RSD
import DLA
import independent
import absorber
import metals

lya = 1215.67

#Function to create a SimulationData object given a specific pixel, information about the complete simulation, and the location/filenames of data files.
def make_gaussian_pixel_object(pixel,original_file_location,original_filename_structure,input_format,MOCKID_lookup,lambda_min=0,IVAR_cutoff=lya):

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
            filename = original_file_location + '/' + original_filename_structure.format(file_number)
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
    def __init__(self,N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP,A):

        self.N_qso = N_qso
        self.N_cells = N_cells
        self.SIGMA_G = SIGMA_G
        self.ALPHA = ALPHA

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
        self.A = A

        #self.absorbers = []
        #self.absorbers.append(absorber.AbsorberData(name='Lya',rest_wave=1215.67,flux_transform_m=1.0))
        #self.absorbers.append(absorber.AbsorberData(name='Lyb',rest_wave=1024.0,flux_transform_m=0.1))
        #self.absorbers.append(absorber.AbsorberData(name='SiII',rest_wave=1205.0,flux_transform_m=0.05))

        self.lya_absorber=absorber.AbsorberData(name='Lya',rest_wave=1215.67,flux_transform_m=1.0)
        self.lyb_absorber=None
        self.metals=None

        return

    def setup_Lyb_absorber(self):

        self.lyb_absorber=absorber.AbsorberData(name='Lyb',rest_wave=1025.72,flux_transform_m=0.1)

        return

    def setup_metal_absorbers(self):

        # get a dictionary with multiple absorbers, one for each metal line 
        self.metals=metals.get_metal_dict()

        return

    #Method to extract reduced data from an input file of a given format, with a given list of MOCKIDs.
    @classmethod
    def get_gaussian_skewers_object(cls,filename,file_number,input_format,MOCKIDs=None,lambda_min=0,IVAR_cutoff=lya,SIGMA_G=None):

        start = time.time()
        times = [0.]

        lya = 1215.67

        h = fits.open(filename)

        h_MOCKID = input.get_MOCKID(h,input_format,file_number)
        h_R, h_Z, h_D, h_V = input.get_COSMO(h,input_format)
        h_lya_lambdas = input.get_lya_lambdas(h,input_format)

        times += [time.time()-start-np.sum(times[:-1])]

        if MOCKIDs != None:
            #Work out which rows in the hdulist we are interested in.
            rows = []
            s = set(MOCKIDs)

            for i, qso in enumerate(h_MOCKID):
                if qso in s:
                    rows += [i]
        else:
            rows = list(range(h_MOCKID.shape[0]))

        #Calculate the first_relevant_cell.
        first_relevant_cell = np.searchsorted(h_lya_lambdas,lambda_min)
        actual_lambda_min = h_lya_lambdas[first_relevant_cell]

        times += [time.time()-start-np.sum(times[:-1])]

        if input_format == 'physical_colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE'][rows]
            RA = h[1].data['RA'][rows]
            DEC = h[1].data['DEC'][rows]
            Z_QSO = h[1].data['Z_COSMO'][rows]
            DZ_RSD = h[1].data['DZ_RSD'][rows]

            DENSITY_DELTA_rows = h[2].data[rows,first_relevant_cell:]

            VEL_rows = h[3].data[rows,first_relevant_cell:]

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
            MOCKID = input.get_MOCKID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Calculate the Gaussian skewers.
            GAUSSIAN_DELTA_rows = convert.lognormal_delta_to_gaussian(DENSITY_DELTA_rows,SIGMA_G,D)

            #Set the remaining variables to None
            DENSITY_DELTA_rows = None
            A = None
            ALPHA = None

            #Insert placeholder values for remaining variables.
            PLATE = MOCKID
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)

            IVAR_rows = general.make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

        elif input_format == 'gaussian_colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE'][rows]
            RA = h[1].data['RA'][rows]
            DEC = h[1].data['DEC'][rows]
            Z_QSO = h[1].data['Z_COSMO'][rows]
            DZ_RSD = h[1].data['DZ_RSD'][rows]

            GAUSSIAN_DELTA_rows = h[2].data[rows,first_relevant_cell:]

            VEL_rows = h[3].data[rows,first_relevant_cell:]

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
                MOCKID = input.get_MOCKID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Set the remaining variables to None
            DENSITY_DELTA_rows = None
            A = None
            ALPHA = None

            #Insert placeholder values for remaining variables.
            PLATE = MOCKID
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)

            IVAR_rows = general.make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

        elif input_format == 'picca_density':

            #Extract data from the HDUlist.
            DENSITY_DELTA_rows = h[0].data.T[rows,first_relevant_cell:]

            IVAR_rows = h[1].data.T[rows,first_relevant_cell:]

            LOGLAM_MAP = h[2].data[first_relevant_cell:]

            RA = h[3].data['RA'][rows]
            DEC = h[3].data['DEC'][rows]
            Z_QSO = h[3].data['Z'][rows]
            PLATE = h[3].data['PLATE'][rows]
            MJD = h[3].data['MJD'][rows]
            FIBER = h[3].data['FIBER'][rows]
            MOCKID = h[3].data['THING_ID'][rows]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = LOGLAM_MAP.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive Z and transmitted flux fraction.
            Z = (10**LOGLAM_MAP)/lya - 1

            #Calculate the Gaussian skewers.
            GAUSSIAN_DELTA_rows = convert.lognormal_delta_to_gaussian(DENSITY_DELTA_rows,SIGMA_G,D)

            #Set the remaining variables to None
            DENSITY_DELTA_rows = None
            A = None
            ALPHA = None

            """
            Can we calculate DZ_RSD,R,D,V?
            """

            #Insert placeholder variables for remaining variables.
            TYPE = np.zeros(RA.shape[0])
            R = np.zeros(Z.shape[0])
            D = np.zeros(Z.shape[0])
            V = np.zeros(Z.shape[0])
            DZ_RSD = np.zeros(RA.shape[0])
            VEL_rows = np.zeros(DENSITY_DELTA_rows.shape)

        else:
            print('Input format not recognised: current options are "colore" and "picca".')
            print('Please choose one of these options and try again.')

        times += [time.time()-start-np.sum(times[:-1])]

        h.close()

        #print('{:3.0%} {:3.0%} {:3.0%} {:3.0%}'.format(times[0]/np.sum(times),times[1]/np.sum(times),times[2]/np.sum(times),times[3]/np.sum(times)))

        return cls(N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP,A)

    #Function to trim skewers according to a minimum value of lambda. QSOs with no relevant cells are removed.
    def trim_skewers(self,lambda_min,min_catalog_z,extra_cells=0,lambda_max=None,whole_lambda_range=False):

        lambdas = 10**(self.LOGLAM_MAP)
        first_relevant_cell = np.searchsorted(lambdas,lambda_min)
        if lambda_max:
            last_relevant_cell = np.searchsorted(lambdas,lambda_max) - 1
        else:
            last_relevant_cell = -1 % self.N_cells

        #If we want to keep any extra_cells, we subtract from the first_relevant_cell.
        first_relevant_cell -= extra_cells

        #Determine which QSOs have any relevant cells to keep.
        """
        relevant_QSOs = []
        for i in range(self.N_qso):
            lambda_QSO = lya*(1 + self.Z_QSO[i])
            if self.IVAR_rows[i,first_relevant_cell] > 0:
                relevant_QSOs += [i]
        """
        relevant_QSOs = (self.Z_QSO>min_catalog_z)
        #If we want the entirety of the lambda range to be relevant (i.e. with IVAR=1), we must remove skewers that do not have this
        if whole_lambda_range:
            relevant_QSOs *= (self.IVAR_rows[:,first_relevant_cell] == 1) * (self.IVAR_rows[:,last_relevant_cell + 1] == 1)

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

        #Now trim the skewers of the remaining QSOs.
        self.N_cells = last_relevant_cell - first_relevant_cell + 1

        self.GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[:,first_relevant_cell:last_relevant_cell + 1]
        if self.DENSITY_DELTA_rows is not None:
            self.DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[:,first_relevant_cell:last_relevant_cell + 1]
        self.VEL_rows = self.VEL_rows[:,first_relevant_cell:last_relevant_cell + 1]
        self.IVAR_rows = self.IVAR_rows[:,first_relevant_cell:last_relevant_cell + 1]
        if self.lya_absorber.tau_computed():
            self.lya_absorber.tau = self.lya_absorber.tau[:,first_relevant_cell:last_relevant_cell + 1]

        self.R = self.R[first_relevant_cell:last_relevant_cell + 1]
        self.Z = self.Z[first_relevant_cell:last_relevant_cell + 1]
        self.D = self.D[first_relevant_cell:last_relevant_cell + 1]
        self.V = self.V[first_relevant_cell:last_relevant_cell + 1]
        self.LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:last_relevant_cell + 1]

        return

    #Function to add small scale gaussian fluctuations.
    def add_small_scale_gaussian_fluctuations(self,cell_size,sigma_G_z_values,extra_sigma_G_values,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya,n=0.7,k1=0.001):

        # TODO: Is NGP really the way to go?

        #Add small scale fluctuations
        old_R = self.R
        Rmax = np.max(old_R)
        Rmin = np.min(old_R)
        new_R = np.arange(Rmin,Rmax,cell_size)
        new_N_cells = new_R.shape[0]

        NGPs = general.get_NGPs(old_R,new_R).astype(int)
        expanded_GAUSSIAN_DELTA_rows = np.zeros((self.N_qso,new_N_cells))

        for i in range(self.N_qso):
            expanded_GAUSSIAN_DELTA_rows[i,:] = self.GAUSSIAN_DELTA_rows[i,NGPs]

        #Redefine the necessary variables (N_cells, Z, D etc)
        self.N_cells = new_N_cells
        self.R = new_R

        # TODO: Ideally would want to recompute these rather than interpolating?
        self.Z = np.interp(new_R,old_R,self.Z)
        self.D = np.interp(new_R,old_R,self.D)
        self.V = np.interp(new_R,old_R,self.V)
        self.LOGLAM_MAP = np.log10(lya*(1+self.Z))

        # TODO: What to do with this?
        self.VEL_rows = self.VEL_rows[:,NGPs]

        #Make new IVAR rows.
        self.IVAR_rows = general.make_IVAR_rows(IVAR_cutoff,self.Z_QSO,self.LOGLAM_MAP)

        #For each skewer, determine the last relevant cell
        first_relevant_cells = np.zeros(self.N_qso)
        last_relevant_cells = np.zeros(self.N_qso)
        for i in range(self.N_qso):
            first_relevant_cell = np.searchsorted(10**(self.LOGLAM_MAP),lambda_min)
            # it is not clear whether to cut at Z_QSO or Z_QSO + DZ_RSD
            last_relevant_cell = np.searchsorted(self.Z,self.Z_QSO[i]) - 1

            #Clip the gaussian skewers so that they are zero after the quasar.
            #This avoids effects from NGP interpolation).
            expanded_GAUSSIAN_DELTA_rows[i,last_relevant_cell + 1:] = 0

            first_relevant_cells[i] = first_relevant_cell
            last_relevant_cells[i] = last_relevant_cell

        extra_var = np.zeros(expanded_GAUSSIAN_DELTA_rows.shape)
        extra_sigma_G = np.interp(self.Z,sigma_G_z_values,extra_sigma_G_values)


        # TODO: dv is not constant at the moment - how to deal with this
        #Generate extra variance, either white noise or correlated.
        dkms_dhMpc = general.get_dkms_dhMpc(0.)
        dv_kms = cell_size * dkms_dhMpc
        extra_var = independent.get_gaussian_fields(generator,self.N_cells,dv_kms=dv_kms,N_skewers=self.N_qso,white_noise=white_noise,n=n,k1=k1)

        #Normalise the extra variance to have unit variance
        k_kms = np.fft.rfftfreq(self.N_cells)*2*np.pi/dv_kms
        mean_P = np.average(independent.power_kms(0.,k_kms,dv_kms,white_noise))
        extra_var /= np.sqrt(mean_P/dv_kms)

        extra_var *= extra_sigma_G

        mask = general.make_IVAR_rows(lya,self.Z_QSO,self.LOGLAM_MAP)
        extra_var *= mask



        """
        # TODO: Improve this
        #1 by 1? Or just N_skewers=N_qso?
        extra_variance = get_gaussian_fields(z,self.N_cells,sigma_G=extra_sigma_G,dv_kms=10.0,N_skewers=self.N_qso,white_noise=white_noise)
        final_GAUSSIAN_DELTA_rows = expanded_GAUSSIAN_DELTA_rows + amplitude*extra_variance
        """

        """
        for j in range(self.N_cells):
            relevant_QSOs = [i for i in range(self.N_qso) if first_relevant_cells[i]<=j and last_relevant_cells[i]>=j]
            extra_var[relevant_QSOs,j] = generator.normal(scale=extra_sigma_G[j],size=len(relevant_QSOs))
        """


        expanded_GAUSSIAN_DELTA_rows += amplitude*extra_var

        """
        for i in range(self.N_qso):
            first_relevant_cell = first_relevant_cells[i].astype('int32')
            last_relevant_cell = last_relevant_cells[i].astype('int32')

            #Number of cells needed is either the dist between the first and last relevant cells, or 0
            N_cells_needed = np.max([(last_relevant_cell - first_relevant_cell).astype('int32'),0])

            extra_var = np.zeros(N_cells_needed)


            #Pass the generator to get_gaussian_skewers, along with the required sigma, and the size, and the seed
            seed = self.MOCKID[i]
            extra_var = get_gaussian_skewers(generator,N_cells_needed,extra_sigma_G[first_relevant_cell:last_relevant_cell],new_seed=seed)
            #Generate a skewer of the right size with the given seed into 'extra_var'
            #Add on extra var *amplitude


            for j in range(first_relevant_cell,last_relevant_cell):
                extra_sigma_G_cell = extra_sigma_G[j]
                extra_var[j-first_relevant_cell] = get_gaussian_skewers(1,extra_sigma_G_cell)

            if last_relevant_cell >= 0:

                expanded_GAUSSIAN_DELTA_rows[i,first_relevant_cell:last_relevant_cell] += amplitude*extra_var
            """
        self.GAUSSIAN_DELTA_rows = expanded_GAUSSIAN_DELTA_rows
        self.SIGMA_G = np.sqrt(extra_sigma_G**2 + (self.SIGMA_G)**2)

        dtype = [('R', 'f8'), ('Z', 'f8'), ('D', 'f8'), ('V', 'f8')]
        new_cosmology = np.array(list(zip(self.R,self.Z,self.D,self.V)),dtype=dtype)

        return new_cosmology

    #Function to add physical skewers to the object via a lognormal transformation.
    def compute_physical_skewers(self,density_type='lognormal'):

        self.DENSITY_DELTA_rows = convert.gaussian_to_lognormal_delta(self.GAUSSIAN_DELTA_rows,self.SIGMA_G,self.D)

        return

    #Function to add physical skewers to the object via a lognormal transformation.
    def compute_tau_skewers(self,absorber,alpha,beta):

        # scale optical depth for this particular absorber (=1 for Lya)
        absorber_alpha = alpha*absorber.flux_transform_m
        absorber.tau = convert.density_to_tau(self.DENSITY_DELTA_rows+1,alpha,beta)

        #Set tau to 0 beyond the quasars.
        for i in range(self.N_qso):
            last_relevant_cell = np.searchsorted(self.Z,self.Z_QSO[i]) - 1
            absorber.tau[i,last_relevant_cell+1:] = 0

        return

    def compute_all_tau_skewers(self,alpha,beta):

        # for each absorber, compute its optical depth skewers
        self.compute_tau_skewers(self.lya_absorber,alpha,beta)

        # optical depth for Ly_b
        if self.lyb_absorber is not None:
            self.compute_tau_skewers(self.lyb_absorber,alpha,beta)
    
        # loop over metals in dictionary
        if self.metals is not None:
            for metal in iter(self.metals.values()):
                self.compute_tau_skewers(metal,alpha,beta)

        return


    #Function to add thermal RSDs from the velocity skewers.
    def add_RSDs(self,absorber,alpha,beta,thermal=False):

        density = 1 + self.DENSITY_DELTA_rows
        new_tau = RSD.add_skewer_RSDs(absorber.tau,density,self.VEL_rows,self.Z,self.R,thermal=thermal)

        #Overwrite the tau skewers and set a flag to True.
        absorber.tau = new_tau

        return


    def add_all_RSDs(self,alpha,beta,thermal=False):

        # for each absorber, add RSDs
        self.add_RSDs(self.lya_absorber,alpha,beta,thermal)

        # RSD for Ly-b
        if self.lyb_absorber is not None:
            self.add_RSDs(self.lyb_absorber,alpha,beta,thermal)

        # loop over metals in dictionary
        if self.metals is not None:
            for metal in iter(self.metals.values()):
                self.add_RSDs(metal,alpha,beta,thermal)

        return


    #Function to measure mean flux.
    def get_mean_flux(self,absorber,z_value=None,z_width=None):

        F = absorber.transmission()
        if not z_value:
            mean_F = np.average(F,axis=0)

        elif not z_width:
            j_value_upper = np.searchsorted(self.Z,z_value)
            j_value_lower = j_value_upper - 1

            if j_value_lower > -1:
                weight_upper = (z_value - self.Z[j_value_lower])/(self.Z[j_value_upper] - self.Z[j_value_lower])
                weight_lower = (self.Z[j_value_upper] - z_value)/(self.Z[j_value_upper] - self.Z[j_value_lower])

            else:
                weight_upper = 1
                weight_lower = 0

            weights = np.ones((self.N_qso,2))
            weights[:,0] *= weight_lower
            weights[:,1] *= weight_upper

            mean_F = np.average(F[:,j_value_lower:j_value_upper+1],weights=weights)

        else:
            j_value_upper = np.searchsorted(self.Z,z_value + z_width/2.)
            j_value_lower = np.max(0,np.searchsorted(self.Z,z_value - z_width/2.) - 1)
            mean_F = np.average(F[j_value_lower:j_value_upper+1])
            #print(self.N_qso)
            #print(j_value_lower,j_value_upper)
        return mean_F

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
        ALPHA = object_A.ALPHA

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
        A = object_A.A

        return cls(N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP,A)

    #Function to save data as a Gaussian colore file.
    def save_as_gaussian_colore(self,location,filename,header,overwrite=False):

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
        hdulist.writeto(location+filename,overwrite=overwrite)
        hdulist.close

        return

    #Function to save data as a picca density file.
    def save_as_picca_gaussian(self,location,filename,header,overwrite=False,zero_mean_delta=False,min_number_cells=2,mean_DELTA=None):

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
            relevant_GAUSSIAN_DELTA_rows = general.normalise_deltas(relevant_GAUSSIAN_DELTA_rows,mean_DELTA)

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
        hdulist.writeto(location+filename,overwrite=overwrite)
        hdulist.close()

        return

    #Function to save data as a Lognormal colore file.
    def save_as_physical_colore(self,location,filename,header):

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
        hdulist.writeto(location+filename)
        hdulist.close

        return

    #Function to save data as a picca density file.
    def save_as_picca_density(self,location,filename,header,zero_mean_delta=False,min_number_cells=2,mean_DELTA=None):

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
            relevant_DENSITY_DELTA_rows = general.normalise_deltas(relevant_DENSITY_DELTA_rows,mean_DELTA)

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
        hdulist.writeto(location+filename)
        hdulist.close()

        return

    #Compute transmission for a particular absorber, on a particular grid
    def compute_grid_transmission(self,absorber,wave_grid):
        #Get transmission on each cell, from tau stored in absorber 
        F_skewer = absorber.transmission()
        #Get rest-frame wavelength for this particular absorber
        rest_wave = absorber.rest_wave
        #Get wavelength on each original cell, for this particular absorber
        wave_skewer = rest_wave*(1+self.Z)

        # interpolate F into the common grid
        N_los = F_skewer.shape[0]
        N_w = wave_grid.shape[0]
        F_grid = np.empty([N_los,N_w])
        for i in range(N_los):
            F_grid[i,] = np.interp(wave_grid,wave_skewer,F_skewer[i])

        return F_grid

    #Function to save data as a transmission file.
    def save_as_transmission(self,location,filename,header):
        
        # define common wavelength grid to be written in files (in Angstroms)
        wave_min=3550
        wave_max=6500
        wave_step=0.1
        wave_grid=np.arange(wave_min,wave_max,wave_step)

        # now we should loop over the different absorbers, combine them and 
        # write them in HDUs. I suggest to have two HDU:
        # - TRANSMISSION will contain both Lya and Lyb
        # - METALS will contain all metal absorption

        # compute Lyman alpha transmission on grid of wavelengths
        F_grid_Lya = self.compute_grid_transmission(self.lya_absorber,wave_grid)
        
        # compute Lyman beta transmission on grid of wavelengths
        if self.lyb_absorber is None:
            F_grid_Lyb = np.ones_like(F_grid_Lya)
        else:
            F_grid_Lyb = self.compute_grid_transmission(self.lyb_absorber,wave_grid)

        # compute metals' transmission on grid of wavelengths
        F_grid_met = np.ones_like(F_grid_Lya)
        for metal in iter(self.metals.values()):
            F_i=self.compute_grid_transmission(metal,wave_grid) 
            F_grid_met *= F_i

        F_grid_met

        # construct quasar catalog HDU
        Z_RSD = self.Z_QSO + self.DZ_RSD
        catalog_data = list(zip(self.RA,self.DEC,Z_RSD,self.Z_QSO,self.MOCKID))
        dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('Z_noRSD', 'f8'), ('MOCKID', int)]
        catalog_data = np.array(catalog_data,dtype=dtype)

        #Construct HDUs from the data arrays.
        prihdr = fits.Header()
        prihdu = fits.PrimaryHDU(header=prihdr)
        cols_METADATA = fits.ColDefs(catalog_data)
        hdu_METADATA = fits.BinTableHDU.from_columns(cols_METADATA,header=header,name='METADATA')
        hdu_WAVELENGTH = fits.ImageHDU(data=wave_grid,header=header,name='WAVELENGTH')
              #Gives transmission with the different species
        #hdu_TRANSMISSION = fits.ImageHDU(data=F_grid_Lya*F_grid_Lyb,header=header,name='TRANSMISSION')
        #hdu_TRANSMISSION = fits.ImageHDU(data=F_grid_Lya,header=header,name='TRANSMISSION')
        hdu_TRANSMISSION = fits.ImageHDU(data=F_grid_met,header=header,name='TRANSMISSION')

        #Combine the HDUs into an HDUlist (including DLAs, if they have been computed)
        if hasattr(self,'DLA_table') == True:
            hdu_DLAs = fits.hdu.BinTableHDU(data=self.DLA_table,header=header,name='DLA')
            hdulist = fits.HDUList([prihdu, hdu_METADATA, hdu_WAVELENGTH, hdu_TRANSMISSION, hdu_DLAs])
        else:
            hdulist = fits.HDUList([prihdu, hdu_METADATA, hdu_WAVELENGTH, hdu_TRANSMISSION])

        #Save as a new file. Close the HDUlist.
        hdulist.writeto(location+filename)
        hdulist.close()

        return

    #Function to save data as a picca flux file.
    def save_as_picca_flux(self,location,filename,header,min_number_cells = 2,mean_F_data=None):

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
        hdulist.writeto(location+filename)
        hdulist.close()

        return

    #Function to save data as a picca velocity file.
    def save_as_picca_velocity(self,location,filename,header,zero_mean_delta=False,min_number_cells=2,overwrite=False):

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
        hdulist.writeto(location+filename,overwrite=overwrite)
        hdulist.close()

        return

    #Function to save the mean and variance of the different quantities as a function of Z.
    def get_means(self,lambda_min=0.0):

        #Determine the relevant cells and QSOs.
        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the first cell which corresponds to a lya_line at wavelength > lambda_min
        first_relevant_cell = general.get_first_relevant_index(lambda_min,lya_lambdas)

        #Determine the furthest cell which is still relevant: i.e. the one in which at least one QSO has non-zero value of IVAR.
        furthest_QSO_index = np.argmax(self.Z_QSO)
        #last_relevant_cell = (np.argmax(self.IVAR_rows[furthest_QSO_index,:]==0) - 1) % self.N_cells
        last_relevant_cell = self.N_cells - 1

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        relevant_QSOs = [i for i in range(self.N_qso) if self.IVAR_rows[i,first_relevant_cell] == 1]

        #Trim data according to the relevant cells and QSOs.
        relevant_DENSITY_DELTA = self.DENSITY_DELTA_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_GAUSSIAN_DELTA = self.GAUSSIAN_DELTA_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        F = self.lya_absorber.transmission()
        relevant_F = F[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_IVAR = self.IVAR_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:last_relevant_cell+1]

        #For each cell, determine the number of skewers for which it is relevant.
        N_relevant_skewers = np.sum(relevant_IVAR,axis=0)
        relevant_cells = N_relevant_skewers>0

        #Calculate delta_F from F.
        #Introduce a small 'hack' in order to get around the problem of having cells with no skewers contributing to them.
        # TODO: find a neater way to deal with this
        small = 1.0e-10
        relevant_mean_F = np.average(relevant_F,weights=relevant_IVAR+small,axis=0)
        relevant_delta_F = ((relevant_F)/relevant_mean_F - 1)*relevant_IVAR

        #Calculate the mean in each cell of the gaussian delta and its square.
        GDB = np.average(relevant_GAUSSIAN_DELTA,weights=relevant_IVAR+small,axis=0)*relevant_cells
        GDSB = np.average(relevant_GAUSSIAN_DELTA**2,weights=relevant_IVAR+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the density delta and its square.
        DDB = np.average(relevant_DENSITY_DELTA,weights=relevant_IVAR+small,axis=0)*relevant_cells
        DDSB = np.average(relevant_DENSITY_DELTA**2,weights=relevant_IVAR+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the flux and its square.
        FB = np.average(relevant_F,weights=relevant_IVAR+small,axis=0)*relevant_cells
        FSB = np.average(relevant_F**2,weights=relevant_IVAR+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the flux delta and its square.
        FDB = np.average(relevant_delta_F,weights=relevant_IVAR+small,axis=0)*relevant_cells
        FDSB = np.average(relevant_delta_F**2,weights=relevant_IVAR+small,axis=0)*relevant_cells

        #Stitch together the means into a binary table.
        dtype = [('N', 'f4'),('GAUSSIAN_DELTA', 'f4'), ('GAUSSIAN_DELTA_SQUARED', 'f4'), ('DENSITY_DELTA', 'f4'), ('DENSITY_DELTA_SQUARED', 'f4')
                , ('F', 'f4'), ('F_SQUARED', 'f4'), ('F_DELTA', 'f4'), ('F_DELTA_SQUARED', 'f4')]
        means = np.array(list(zip(N_relevant_skewers,GDB,GDSB,DDB,DDSB,FB,FSB,FDB,FDSB)),dtype=dtype)

        return means

    #Function to add DLAs to a set of skewers.
    def add_DLA_table(self):

        dla_bias = 2.0
        #If extrapolate_z_down is set to a value below the skewer, then we extrapolate down to that value.
        #Otherwise, we start placing DLAs at the start of the skewer.
        extrapolate_z_down = None
        DLA.add_DLA_table_to_object(self,dla_bias=dla_bias,extrapolate_z_down=extrapolate_z_down)

        return
