import matplotlib.pyplot as plt
import numpy as np
import pyfits
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from PyAstronomy import pyaGui
from PyAstronomy import pyasl
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord

def spec_readspec_multiext(file,header='NO OUTPUT HEADER',ext_num=0):
    wave=[]
    flux=[]
    hdr=[]
    split = file.split('.')
    exten = split[-1]
    if (exten == 'fits') or (exten == 'fit'):
        hdu = pyfits.open(file)
        hdr = hdu[ext_num].header
        if 'crpix1' in hdu[ext_num].header:
            flux = hdu[ext_num].data[0] 
            wave = readlambda_multiext(flux,hdu,ext_num=ext_num)

        else:
            print('!!!	Wavelength keyword not found in FITS HEADER 	!!!')
            return
    else:
        print("Not yet supported!")
        return
    hdu.close	
    if header == 'NO OUTPUT HEADER':
        return wave,flux
    else:
        return wave,flux,hdr


def readlambda_multiext(spec, hdu_sp,ext_num=0):
    
    crpix1 = hdu_sp[ext_num].header['crpix1']
    #/ value of ref pixel
    crval1 = hdu_sp[ext_num].header['crval1']
    #/ delta per pixel
    if 'cd1_1' in hdu_sp[ext_num].header:
        cd1_1 = hdu_sp[ext_num].header['cd1_1']
    #cd1_1 is sometimes called cdelt1.
    else:
        cd1_1 = hdu_sp[ext_num].header['cdelt1']
    if cd1_1 == 0:
        print("NOT WORKING")
        return
#    print "CR Pix",crpix1,"CR VAL1",crval1,"CD 1 1",cd1_1

    n_lambda = len(spec)
    wave = np.zeros(n_lambda)
    for l  in xrange(n_lambda):
        wave[l] = (l+1.0-crpix1)*cd1_1+crval1
        #Use pixel number starting with 0 if no lambda information is found.
	if (np.min(wave)+np.max(wave) == 0.0):
		print('No lambda information found: used pixel number starting with 0')
		for l  in xrange(n_lambda):
                    wave[l] = l	
    return wave

def readlambda_multiext_cube(cube, hdu_cube,ext_num=0):
    # Here you look for the 3rd conponent as that is the wl 
    # direction
    crpix3 = hdu_cube[ext_num].header['crpix3']
    #/ value of ref pixel
    crval3 = hdu_cube[ext_num].header['crval3']
    #/ delta per pixel
    if 'cd3_3' in hdu_cube[ext_num].header:
        cd3_3 = hdu_cube[ext_num].header['cd3_3']
    #cd1_ is sometimes called cdelt3.
    #In the case of KMOS cube they are both the
    #same. 
        print crpix3,crval3,cd3_3
    else:
        cd3_3 = hdu_cube[ext_num].header['cdelt3']
    if cd3_3 == 0:
        print("NOT WORKING")
        return

    n_lambda = len(cube)
    wave = np.zeros(n_lambda)
    for l  in xrange(n_lambda):
        wave[l] = (l+1.0-crpix3)*cd3_3+crval3
        #Use pixel number starting with 0 if no lambda information is found.
	if (np.min(wave)+np.max(wave) == 0.0):
		print('No lambda information found: used pixel number starting with 0')
		for l  in xrange(n_lambda):
			wave[l] = l	
   
    return wave



class SpecCube:

    def __init__(self,file_name,ext_num):
        self.file_name = file_name
        self.ext_num = ext_num
        self.wl,self.flux,self.hdr = spec_readspec_multiext(self.file_name,"Return Header",self.ext_num)

    def get_sn_spectrum(self,spectrum,wl_range):
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_range[0] < wl < wl_range[1]] 
        return np.mean(spectrum[wl_indices])/np.std(spectrum[wl_indices])

    def normalise_spec(self,wl_range):
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_range[0] < wl < wl_range[1]] 
	self.norm_flux = self.flux/np.mean(self.flux[wl_indices])

    def return_MJD(self):
        hdu_list = pyfits.open(self.file_name)
        return hdu_list[self.ext_num].header['MJD']

    def ew(self,wl_win,wl_cont = [0,0,0,0]):
        #Selecting the region of the line
        # and finding the mean flux across the line. (equivalent to integrating across the line)
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_win[0] <= wl <= wl_win[1]]
        mean_flux_line = np.mean(self.norm_flux[wl_indices])


        # Where specific cont windows are given
        # Find estimate of the continuum either side of the line
    	# The mean of this gives you an estimate of what the continuum should be over the line
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_cont[0] <= wl <= wl_cont[1]]
        mean_flux_cont1 = np.mean(self.norm_flux[wl_indices])
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_cont[2] <= wl <= wl_cont[3]]
        mean_flux_cont2 = np.mean(self.norm_flux[wl_indices])
        
        mean_flux_cont = (mean_flux_cont1 + mean_flux_cont2)/2
        EW_tmp = (mean_flux_line - mean_flux_cont)/mean_flux_cont
        EW = EW_tmp*(wl_win[1]-wl_win[0])   # Multiply by the size of the line window "integrated" over
        print "Equivalent Width = %.4f " % EW
        
        return EW


    def sigma_clip_spectra(self,sigma,window_size,wl_range):
        #This is not working very well. 
        print "This one does not work very well"

	len_spectra = len(self.wl)
        
        self.clean_flux = self.flux*0.0

	for element in range(window_size,len_spectra-window_size):
	     #perform sigma clipping
             flux_ele = self.flux[element]
             mean_flux_in_win = np.mean([flux for index,flux in enumerate(self.flux) if element-window_size < index < element+window_size])
             std_flux_in_win = np.std([flux for index,flux in enumerate(self.flux) if element-window_size < index < element+window_size])
 
	     if flux_ele > sigma*std_flux_in_win:
                 self.clean_flux[element] = mean_flux_in_win #0.0
             else:
                 self.clean_flux[element] = flux_ele   

        wl_indices = [index for index,wl in enumerate(self.wl) if wl_range[0] < wl < wl_range[1]]
        self.norm_clean_flux = self.clean_flux/np.mean(self.clean_flux[wl_indices])


    
    def bary_shift_correction(self):

        MJD_float = self.return_MJD()
    
        hdu_list = pyfits.open(self.file_name)
        RA = hdu_list[self.ext_num].header['CAT-RA']
        DEC = hdu_list[self.ext_num].header['CAT-DEC']
        coord_sys = hdu_list[self.ext_num].header['RADECSYS']
        if coord_sys == 'FK5':
            coord_sys = 'fk5'
        elif coord_sys == 'FK4':
            coord_sys = 'fk4'
    
        EPOCH = hdu_list[self.ext_num].header['CAT-EPOC']

        coordinates = SkyCoord(RA+' '+DEC,unit=(u.hourangle,u.deg),frame=coord_sys)
        coord_string = coordinates.to_string('decimal')
        RA_decimal = float(coord_string.split(' ')[0])
        DEC_decimal = float(coord_string.split(' ')[1])


        JD = MJD_float + 2400000.5

        #heli,bary = pyasl.baryvel(JD,deq=EPOCH)
        #print "Heliocentric vel",heli
        #print "Barycentric vel",bary


        vh,vb = pyasl.baryCorr(JD,RA_decimal,DEC_decimal,deq=EPOCH)
        #print "Barycentric velocity of earth towards star",vb
        #print "Heliocentric velocity of earth towards star",vh

        heli_corr = vh*u.kilometer/u.s
        bary_corr = vb*u.kilometer/u.s

        wl_units = hdu_list[self.ext_num].header['CUNIT1']
        print "WL units",wl_units
        if wl_units == 'Angstroms':
            wl_in_units = self.wl*u.AA
        else:
            print "Setting wl units to nm"
            wl_in_units = self.wl*u.nm

        print "Mean Correction",(np.mean(wl_in_units.cgs)*(bary_corr.cgs/const.c.cgs)).to(u.AA)

        wl_corr = wl_in_units.cgs*(1.*u.dimensionless_unscaled + bary_corr.cgs/const.c.cgs)
        if wl_units == 'Angstroms':
            self.wl_corr = wl_corr.to(u.AA)
        else:
            self.wl_corr = wl_corr.to(u.nm)

        hdu_list.close()
        return (np.mean(wl_in_units.cgs)*(bary_corr.cgs/const.c.cgs)).to(u.AA)
 


class Datacube:

    def __init__(self,file_name,ext_num):
        self.file_name = file_name
        self.ext_num = ext_num
	self.get_cube()

    def get_cube(self):
        self.hdu_list = pyfits.open(self.file_name)
        #Extension specific to FRODOspec data 
        raw_data_cube = pyfits.getdata(self.file_name,self.ext_num)
        
        # set nan's in data cube to zero
        print "nan's --> zero"
        self.cube = np.nan_to_num(raw_data_cube)
        
    def get_sum_cube(self):
        self.cube_sum = np.sum(self.cube,0)
        
    def get_cube_wl(self):
        self.wl = readlambda_multiext_cube(self.cube,self.hdu_list,self.ext_num)
	
    def cube_slice(self,wl_begin,wl_end):
        #Selecting a wavelength slice of the cube
        self.get_cube_wl()
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_begin < wl <wl_end] 
        min_wl_index = np.amin(wl_indices)
        max_wl_index = np.amax(wl_indices)

        #return the sum of the resulting slice
        return np.sum(self.cube[min_wl_index:max_wl_index,:,:],0)

    def get_sn_spectrum(self,spectrum,wl_range):
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_range[0] < wl < wl_range[1]] 
        return np.mean(spectrum[wl_indices])/np.std(spectrum[wl_indices])

    def normalise_spec(self,wl_range):
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_range[0] < wl < wl_range[1]]
        self.norm_flux = self.sum_spectrum/np.mean(self.sum_spectrum[wl_indices])

    def return_obs_date(self):
        hdu_list = pyfits.open(self.file_name)
        return hdu_list[self.ext_num].header['DATE-OBS']

    def return_MJD(self):
        hdu_list = pyfits.open(self.file_name)
        return hdu_list[self.ext_num].header['MJD']



    def extract_spec_sum_pixels(self,cube_num_pix,sn_wl_range):
        ## Extracting spectra through adding succesive
        ## pixels. Pixels are chosen according to the max
        ## value in the cube
        s2n = []
        self.get_cube_wl()
        #So that sum_cube is not overwritten
        pixel_map = 1*self.cube_sum
        for i in range(0,cube_num_pix):
            #Find max pixel
            maxINDEX = np.where(pixel_map == np.max(pixel_map))
            if i == 0:
                sum_spectrum = self.cube[:,maxINDEX[0],maxINDEX[1]]
            else:
                sum_spectrum += self.cube[:,maxINDEX[0],maxINDEX[1]]

            s2n.append(self.get_sn_spectrum(sum_spectrum,sn_wl_range))
            pixel_map[maxINDEX[0],maxINDEX[1]] = 1 
            # Set to 1 so that it is different from NAN

#        print "sum spectrum",sum_spectrum.shape

        for i in range(0,len(pixel_map[:,0])):
                for j in range(0,len(pixel_map[0,:])):
                    if pixel_map[i,j] != 1:
                        pixel_map[i,j] = 0

        self.sum_spectrum = sum_spectrum
        self.extract_map = pixel_map

        return sum_spectrum,s2n,pixel_map

    def find_sky_pixels(self):
       
        #These are the pixels not used in the
        #spectrum extraction
        sky_map = np.subtract(self.extract_map,1)
        sky_map = np.multiply(sky_map,-1)
        
        sky_cube = np.multiply(self.cube,0)
    
        # Multiply the data cube by this sky map
        for i in range(0,len(sky_map[:,0])):
                for j in range(0,len(sky_map[0,:])):
                    sky_cube[:,i,j] = sky_map[i,j]*self.cube[:,i,j]


        sky_cube_sum = np.sum(sky_cube,0)

        #Identify level for sky signal
        median_av_sky = np.median(sky_cube_sum)
            #Find max pixel
        INDEX = np.where(sky_cube_sum <= median_av_sky)
    

        sky_spectrum = sky_cube[:,INDEX[0],INDEX[1]]
        sky_spectrum_sum = np.sum(sky_spectrum,1)
    
        return sky_spectrum_sum


    def ew(self,wl_win,wl_cont = [0,0,0,0]):
        #Selecting the region of the line
        # and finding the mean flux across the line. (equivalent to integrating across the line)
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_win[0] <= wl <= wl_win[1]]
        mean_flux_line = np.mean(self.norm_flux[wl_indices])


        # Where specific cont windows are given
        # Find estimate of the continuum either side of the line
    	# The mean of this gives you an estimate of what the continuum should be over the line
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_cont[0] <= wl <= wl_cont[1]]
        mean_flux_cont1 = np.mean(self.norm_flux[wl_indices])
        wl_indices = [index for index,wl in enumerate(self.wl) if wl_cont[2] <= wl <= wl_cont[3]]
        mean_flux_cont2 = np.mean(self.norm_flux[wl_indices])
        
        mean_flux_cont = (mean_flux_cont1 + mean_flux_cont2)/2
        EW_tmp = (mean_flux_line - mean_flux_cont)/mean_flux_cont
        EW = EW_tmp*(wl_win[1]-wl_win[0])   # Multiply by the size of the line window "integrated" over
        print "Equivalent Width = %.4f " % EW
        
        return EW




    def check_EW_parameters(self):  
        pp = pyaGui.Picker()
        pp.a.plot(self.wl,self.norm_flux)
        points = pp.pick()
        wl_win = [points[0][0],points[1][0]]
        wlc11,wlc12,wlc21,wlc22 = points[2][0],points[3][0],points[4][0],points[5][0]
        wl_cont_win = [wlc11,wlc12,wlc21,wlc22]

        plt.plot(self.wl,self.norm_flux)
        plt.plot([wl_win[0],wl_win[0]],[0.9*np.min(self.norm_flux),1.1*np.max(self.norm_flux)])
        plt.plot([wl_win[1],wl_win[1]],[0.9*np.min(self.norm_flux),1.1*np.max(self.norm_flux)])
        plt.plot([wlc11,wlc12],[1,1])
        plt.plot([wlc21,wlc22],[1,1])
        plt.xlim([wlc11,wlc22])
        plt.ylim([0.9*np.min(self.norm_flux),1.1*np.max(self.norm_flux)])
        plt.show()

        EQW=self.ew(wl_win,wl_cont_win)

        return wl_win,wl_cont_win

	


	

