import numpy as np
import healpy as hp

class Footprint:

    def __init__(self,footprint_name,nside=16,pixels=None):

        #See if we can use desimodel. This is preferable as it will be the most
        #up-do-date footprint.
        try:
            from desimodel.footprint import tiles2pix, is_point_in_desi
            from desimodel.io import load_tiles
            self.desimodel_installed = True
        except ModuleNotFoundError:
            print('WARN: desimodel is not installed; footprint pixel data will be read from file.')
            self.desimodel_installed = False

        self.footprint_name = footprint_name
        self.nside = nside
        self.pixels = pixels

        return

    def __call__(self,RA,DEC):

        #If we have desimodel and want to replicate the footprint precisely, use
        #function "is_point_in_desi".
        if self.desimodel_installed and self.footprint_name=='desi':
            tiles = load_tiles()
            w = is_point_in_desi(tiles,RA,DEC)

        #If not, but we still want to filter...
        elif self.footprint_name in ['desi','desi_pixel','desi_pixel_plus']:

            #If desimodel is installed, then we use "tiles2pix" to determine which
            #pixels to include.
            if self.desimodel_installed:
                from desimodel.footprint import tiles2pix
                if footprint=='desi_pixel':
                    valid_pixels = tiles2pix(self.nside)
                elif footprint=='desi_pixel_plus':
                    valid_pixels = tiles2pix(self.nside)
                    valid_pixels = add_pixel_neighbours(valid_pixels)

            #Otherwise, we load pixel lists from file. Note: using desimodel is
            #preferable to loading from file as the footprint could change, and
            #desimodel will be more up to date than the lists in this case.
            else:
                if self.footprint_name=='desi':
                    print('WARN: desimodel is not installed; footprint pixel data will be read from file.')
                    valid_pixels = np.loadtxt('input_files/pixel_footprints/DESI_pixels.txt',dtype=int)
                elif self.footprint_name=='desi_pixel':
                    valid_pixels = np.loadtxt('input_files/pixel_footprints/DESI_pixels.txt',dtype=int)
                elif self.footprint_name=='desi_pixel_plus':
                    valid_pixels = np.loadtxt('input_files/pixel_footprints/DESI_pixels_plus.txt',dtype=int)

            #With a list of valid pixels, we now can make a filter.
            theta = (np.pi/180.0)*(90.0-DEC)
            phi = (np.pi/180.0)*RA
            pix = hp.ang2pix(self.nside,theta,phi,nest=True)
            w = np.in1d(pix,valid_pixels)

        #Else if we don't want to filter at all, set the filter to "None".
        elif self.footprint_name=='full_sky':
            w = None

        else:
            if self.footprint_name is not None:
                print('Footprint {} not recognised; no filter applied.'.format(self.footprint(name)))
            w = None

        return w
