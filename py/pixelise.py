import numpy as np
from astropy.io import fits

import general
import input
import convert
import RSD
import DLA

lya = 1215.67

#Function to create a 'simulation_data' object given a specific pixel, information about the complete simulation, and the location/filenames of data files.
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

        #If there are some relevant quasars, open the data file and make it into a simulation_data object.
        #We use simulation_data.get_reduced_data to avoid loading all of the file's data into the object.
        if N_relevant_qso > 0:
            filename = original_file_location + '/' + original_filename_structure.format(file_number)
            working = simulation_data.get_gaussian_skewers_object(filename,file_number,input_format,MOCKIDs=relevant_MOCKIDs,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff)

        #Combine the data from the working file with that from the files already looked at.
        if files_included > 0:
            combined = simulation_data.combine_files(combined,working,gaussian_only=True)
            files_included += 1
        else:
            combined = working
            files_included += 1

    pixel_object = combined

    return pixel_object

#Definition of a generic 'simulation_data' class, from which it is easy to save in new formats.
class simulation_data:
    #Initialisation function.
    def __init__(self,N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,F_rows,R,Z,D,V,LOGLAM_MAP,A):

        self.N_qso = N_qso
        self.N_cells = N_cells
        self.SIGMA_G = SIGMA_G
        self.ALPHA = ALPHA

        self.TYPE = TYPE
        self.RA = RA
        self.DEC = DEC
        self.Z_QSO = Z_QSO
        self.DZ_RSD = DZ_RSD
        self.MOCKID = MOCKID
        self.PLATE = PLATE
        self.MJD = MJD
        self.FIBER = FIBER

        self.GAUSSIAN_DELTA_rows = GAUSSIAN_DELTA_rows
        self.DENSITY_DELTA_rows = DENSITY_DELTA_rows
        self.VEL_rows = VEL_rows
        self.IVAR_rows = IVAR_rows
        self.F_rows = F_rows

        self.R = R
        self.Z = Z
        self.D = D
        self.V = V
        self.LOGLAM_MAP = LOGLAM_MAP
        self.A = A

        self.linear_skewer_RSDs_added = False
        self.thermal_skewer_RSDs_added = False

        self.density_computed = False
        self.tau_computed = False
        self.flux_computed = False

        return

    #Method to extract reduced data from an input file of a given format, with a given list of MOCKIDs.
    @classmethod
    def get_gaussian_skewers_object(cls,filename,file_number,input_format,MOCKIDs=None,lambda_min=0,IVAR_cutoff=lya,SIGMA_G=None):

        lya = 1215.67

        h = fits.open(filename)

        h_MOCKID = input.get_MOCKID(h,input_format,file_number)
        h_R, h_Z, h_D, h_V = input.get_COSMO(h,input_format)
        h_lya_lambdas = input.get_lya_lambdas(h,input_format)

        if MOCKIDs != None:
            #Work out which rows in the hdulist we are interested in.
            rows = ['']*len(MOCKIDs)
            s = set(MOCKIDs)
            j = 0
            for i, qso in enumerate(h_MOCKID):
                if qso in s:
                    rows[j] = i
                    j = j+1
        else:
            rows = list(range(h_MOCKID.shape[0]))

        #Calculate the first_relevant_cell.
        first_relevant_cell = np.searchsorted(h_lya_lambdas,lambda_min)
        actual_lambda_min = h_lya_lambdas[first_relevant_cell]

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
            F_rows = None

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
            F_rows = None

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
            F_rows = None

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

        h.close()

        return cls(N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,F_rows,R,Z,D,V,LOGLAM_MAP,A)

    #Function to trim skewers according to a minimum value of lambda. QSOs with no relevant cells are removed.
    def trim_skewers(self,lambda_min,min_catalog_z,extra_cells=0):

        lambdas = 10**(self.LOGLAM_MAP)
        first_relevant_cell = np.searchsorted(lambdas,lambda_min)

        #Determine which QSOs have any relevant cells to keep.
        """
        relevant_QSOs = []
        for i in range(self.N_qso):
            lambda_QSO = lya*(1 + self.Z_QSO[i])
            if self.IVAR_rows[i,first_relevant_cell] > 0:
                relevant_QSOs += [i]
        """
        relevant_QSOs = self.Z_QSO>min_catalog_z

        #If we want to keep any extra_cells, we subtract from the first_relevant_cell.
        first_relevant_cell -= extra_cells

        #Remove QSOs no longer needed.
        self.N_qso = len(relevant_QSOs)

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
        if self.density_computed == True:
            self.DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[relevant_QSOs,:]
        self.VEL_rows = self.VEL_rows[relevant_QSOs,:]
        self.IVAR_rows = self.IVAR_rows[relevant_QSOs,:]
        if self.tau_computed == True:
            self.TAU_rows = self.TAU_rows[relevant_QSOs,:]
        if self.flux_computed == True:
            self.F_rows = self.F_rows[relevant_QSOs,:]

        #Now trim the skewers of the remaining QSOs.
        self.N_cells -= first_relevant_cell

        self.GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[:,first_relevant_cell:]
        if self.density_computed == True:
            self.DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[:,first_relevant_cell:]
        self.VEL_rows = self.VEL_rows[:,first_relevant_cell:]
        self.IVAR_rows = self.IVAR_rows[:,first_relevant_cell:]
        if self.tau_computed == True:
            self.TAU_rows = self.TAU_rows[:,first_relevant_cell:]
        if self.flux_computed == True:
            self.F_rows = self.F_rows[:,first_relevant_cell:]

        self.R = self.R[first_relevant_cell:]
        self.Z = self.Z[first_relevant_cell:]
        self.D = self.D[first_relevant_cell:]
        self.V = self.V[first_relevant_cell:]
        self.LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:]

        return

    #Function to add small scale gaussian fluctuations.
    def add_small_scale_gaussian_fluctuations(self,cell_size,sigma_G_z_values,extra_sigma_G_values,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya):

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
            if self.linear_skewer_RSDs_added == True:
                last_relevant_cell = np.searchsorted(self.Z,self.Z_QSO[i]+self.DZ_RSD[i]) - 1
            else:
                last_relevant_cell = np.searchsorted(self.Z,self.Z_QSO[i]) - 1

            #Clip the gaussian skewers so that they are zero after the quasar.
            #This avoids effects from NGP interpolation).
            expanded_GAUSSIAN_DELTA_rows[i,last_relevant_cell + 1:] = 0

            first_relevant_cells[i] = first_relevant_cell
            last_relevant_cells[i] = last_relevant_cell

        extra_var = independent.get_gaussian_fields(generator,self.N_cells,dv_kms=,N_skewers=self.N_qso,white_noise=white_noise)

        """
        # TODO: Improve this
        #1 by 1? Or just N_skewers=N_qso?
        extra_variance = get_gaussian_fields(z,self.N_cells,sigma_G=extra_sigma_G,dv_kms=10.0,N_skewers=self.N_qso,white_noise=white_noise)
        final_GAUSSIAN_DELTA_rows = expanded_GAUSSIAN_DELTA_rows + amplitude*extra_variance
        """

        extra_sigma_G = np.interp(self.Z,sigma_G_z_values,extra_sigma_G_values)

        extra_var = np.zeros(expanded_GAUSSIAN_DELTA_rows.shape)

        for j in range(self.N_cells):
            relevant_QSOs = [i for i in range(self.N_qso) if first_relevant_cells[i]<=j and last_relevant_cells[i]>=j]
            extra_var[relevant_QSOs,j] = generator.normal(scale=extra_sigma_G[j],size=len(relevant_QSOs))

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
        self.density_computed = True

        return

    #Function to add physical skewers to the object via a lognormal transformation.
    def compute_tau_skewers(self,alpha,beta):

        self.TAU_rows = convert.density_to_tau(self.DENSITY_DELTA_rows+1,alpha,beta)
        self.tau_computed = True

        return

    #Function to add flux skewers to the object.
    def compute_flux_skewers(self):

        #self.TAU_rows = get_tau(self.Z,self.DENSITY_DELTA_rows+1,alpha,beta)
        self.F_rows = np.exp(-self.TAU_rows)
        #self.F_rows = density_to_flux(self.DENSITY_DELTA_rows+1,alpha,beta)

        #Set the skewers to 1 beyond the quasars.
        for i in range(self.N_qso):
            if self.linear_skewer_RSDs_added == True:
                last_relevant_cell = np.searchsorted(self.Z,self.Z_QSO[i]+self.DZ_RSD[i]) - 1
            else:
                last_relevant_cell = np.searchsorted(self.Z,self.Z_QSO[i]) - 1
            self.F_rows[i,last_relevant_cell+1:] = 1

        self.flux_computed = True

        return

    ## TODO: remove this, now defunct
    #Function to add linear RSDs from the velocity skewers.
    def add_linear_RSDs(self,alpha,beta):

        #add RSDs to these physical density rows
        new_TAU_rows = RSD.add_linear_skewer_RSDs(self.TAU_rows,self.VEL_rows,self.Z)

        ## TODO: find a neater way to do this
        #For the moment, we add a very small value onto the tau skewers, to avoid problems in the inverse lognormal transformation
        #In future, when we don't care about the gaussian skewers, we can get rid of this
        moodified_new_TAU_rows = new_TAU_rows + (new_TAU_rows==0)*1.0e-10

        #convert the new tau rows back to physical density
        new_density_rows = convert.tau_to_density(moodified_new_TAU_rows,alpha,beta)
        new_density_delta_rows = new_density_rows - 1

        #convert the new physical density rows back to gaussian
        new_gaussian_rows = convert.lognormal_delta_to_gaussian(new_density_delta_rows,self.SIGMA_G,self.D)

        #Make a mask where the physical skewers are zero.
        #mask = (new_density_rows != 0)
        #self.IVAR_rows *= mask

        #Overwrite the physical and tau skewers and set a flag to True.
        self.TAU_rows = new_TAU_rows
        self.DENSITY_DELTA_rows = new_density_delta_rows
        self.GAUSSIAN_DELTA_rows = new_gaussian_rows
        self.linear_skewer_RSDs_added = True

        return

    #Function to add thermal RSDs from the velocity skewers.
    def add_RSDs(self,alpha,beta,thermal=False):

        initial_density_rows = 1 + self.DENSITY_DELTA_rows
        new_TAU_rows = RSD.add_skewer_RSDs(self.TAU_rows,initial_density_rows,self.VEL_rows,self.Z,self.R,thermal=thermal)

        ## TODO: find a neater way to do this
        #For the moment, we add a very small value onto the tau skewers, to avoid problems in the inverse lognormal transformation
        #In future, when we don't care about the gaussian skewers, we can get rid of this
        moodified_new_TAU_rows = new_TAU_rows + (new_TAU_rows==0)*1.0e-10

        #convert the new tau rows back to physical density
        new_density_rows = convert.tau_to_density(moodified_new_TAU_rows,alpha,beta)
        new_density_delta_rows = new_density_rows - 1

        #convert the new physical density rows back to gaussian
        new_gaussian_rows = convert.lognormal_delta_to_gaussian(new_density_delta_rows,self.SIGMA_G,self.D)

        #Make a mask where the physical skewers are zero.
        #mask = (new_density_rows != 0)
        #self.IVAR_rows *= mask

        #Overwrite the physical and tau skewers and set a flag to True.
        self.TAU_rows = new_TAU_rows
        self.DENSITY_DELTA_rows = new_density_delta_rows
        self.GAUSSIAN_DELTA_rows = new_gaussian_rows
        self.thermal_skewer_RSDs_added = True

        return

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

        if not gaussian_only:
            GAUSSIAN_DELTA_rows = np.concatenate((object_A.GAUSSIAN_DELTA_rows,object_B.GAUSSIAN_DELTA_rows),axis=0)
            DENSITY_DELTA_rows = np.concatenate((object_A.DENSITY_DELTA_rows,object_B.DENSITY_DELTA_rows),axis=0)
            VEL_rows = np.concatenate((object_A.VEL_rows,object_B.VEL_rows),axis=0)
            IVAR_rows = np.concatenate((object_A.IVAR_rows,object_B.IVAR_rows),axis=0)
            F_rows = np.concatenate((object_A.F_rows,object_B.F_rows),axis=0)
        else:
            GAUSSIAN_DELTA_rows = np.concatenate((object_A.GAUSSIAN_DELTA_rows,object_B.GAUSSIAN_DELTA_rows),axis=0)
            DENSITY_DELTA_rows = None
            VEL_rows = np.concatenate((object_A.VEL_rows,object_B.VEL_rows),axis=0)
            IVAR_rows = np.concatenate((object_A.IVAR_rows,object_B.IVAR_rows),axis=0)
            F_rows = None

        """
        Something to check this is ok?
        """

        Z = object_A.Z
        LOGLAM_MAP = object_A.LOGLAM_MAP
        R = object_A.R
        D = object_A.D
        V = object_A.V
        A = object_A.A

        return cls(N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,F_rows,R,Z,D,V,LOGLAM_MAP,A)

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

    #Function to save data as a transmission file.
    def save_as_transmission(self,location,filename,header):
        lya_lambdas = 10**self.LOGLAM_MAP

        Z_RSD = self.Z_QSO + self.DZ_RSD

        transmission_1_data = list(zip(self.RA,self.DEC,Z_RSD,self.Z_QSO,self.MOCKID))

        dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('Z_noRSD', 'f8'), ('MOCKID', int)]
        transmission_1 = np.array(transmission_1_data,dtype=dtype)

        transmission_2 = 10**(self.LOGLAM_MAP)
        transmission_3 = self.F_rows

        #Construct HDUs from the data arrays.
        prihdr = fits.Header()
        prihdu = fits.PrimaryHDU(header=prihdr)
        cols_METADATA = fits.ColDefs(transmission_1)
        hdu_METADATA = fits.BinTableHDU.from_columns(cols_METADATA,header=header,name='METADATA')
        hdu_WAVELENGTH = fits.ImageHDU(data=transmission_2,header=header,name='WAVELENGTH')
        hdu_TRANSMISSION = fits.ImageHDU(data=transmission_3,header=header,name='TRANSMISSION')

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

        #Trim data according to the relevant cells and QSOs.
        relevant_F_rows = self.F_rows[relevant_QSOs,:]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,:]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[:]
        relevant_Z = self.Z[:]

        #Calculate mean F as a function of z for the relevant cells, then F_DELTA_rows.
        try:
            mean_F_z_values = mean_F_data[:,0]
            mean_F = mean_F_data[:,1]
            relevant_F_BAR = np.interp(relevant_Z,mean_F_z_values,mean_F)
        except ValueError:
            #This is done with a 'hack' to avoid problems with weights summing to zero.
            small = 1.0e-10
            relevant_F_BAR = np.average(relevant_F_rows,weights=relevant_IVAR_rows+small,axis=0)

        relevant_F_DELTA_rows = ((relevant_F_rows)/relevant_F_BAR - 1)*relevant_IVAR_rows

        #Organise the data into picca-format arrays.
        picca_0 = relevant_F_DELTA_rows.T
        picca_1 = relevant_IVAR_rows.T
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
        relevant_VEL_rows = self.VEL_rows[relevant_QSOs,:]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,:]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[:]

        #Organise the data into picca-format arrays.
        picca_0 = relevant_VEL_rows.T
        picca_1 = relevant_IVAR_rows.T
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
        relevant_DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_F_rows = self.F_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:last_relevant_cell+1]

        #For each cell, determine the number of skewers for which it is relevant.
        N_relevant_skewers = np.sum(relevant_IVAR_rows,axis=0)
        relevant_cells = N_relevant_skewers>0

        #Calculate F_DELTA_rows from F_rows.
        #Introduce a small 'hack' in order to get around the problem of having cells with no skewers contributing to them.
        # TODO: find a neater way to deal with this
        small = 1.0e-10
        relevant_F_BAR = np.average(relevant_F_rows,weights=relevant_IVAR_rows+small,axis=0)
        relevant_F_DELTA_rows = ((relevant_F_rows)/relevant_F_BAR - 1)*relevant_IVAR_rows

        #Calculate the mean in each cell of the gaussian delta and its square.
        GDB = np.average(relevant_GAUSSIAN_DELTA_rows,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells
        GDSB = np.average(relevant_GAUSSIAN_DELTA_rows**2,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the density delta and its square.
        DDB = np.average(relevant_DENSITY_DELTA_rows,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells
        DDSB = np.average(relevant_DENSITY_DELTA_rows**2,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the flux and its square.
        FB = np.average(relevant_F_rows,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells
        FSB = np.average(relevant_F_rows**2,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the flux delta and its square.
        FDB = np.average(relevant_F_DELTA_rows,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells
        FDSB = np.average(relevant_F_DELTA_rows**2,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells

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
