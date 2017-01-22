import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
##Plotting
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
#Data class
from Frodospec_class import *
## Interactive EW measurement
from PyAstronomy import pyaGui
from PyAstronomy import pyasl

from scipy.interpolate import interp1d

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
  
rcParams.update({'lines.linewidth':1.5})
font = {'family' : 'normal','weight' : 'bold','size'   : 22}
rcParams.update({'font.size':12})


C_kms = 299792458e-3
C_ms = 299792458
C_As = 299792458e10



time_stamp = ['07.02.2013','12.13.2013','27.07.2013','20.08.2013','09.09.2013','30.09.2013','21.10.2013','15.11.2013','18.12.2013','14.01.2014','08.03.2014','24.03.2014','05.04.2014','10.09.2014','06.01.2015','12.01.2015','17.01.2015','27.01.2015','26.02.2015']
 

def convert_wl2vel(wl_array,central_wl):
    tmp1 = np.subtract(wl_array,central_wl)
    tmp2 = np.divide(tmp1,central_wl)
    velocity_array =  np.multiply(tmp2,C_kms)
    return velocity_array

def vel_at_wl(central_wl,wl_of_interest):
    tmp1 = np.subtract(wl_of_interest,central_wl)
    tmp2 = np.divide(tmp1,central_wl)
    velocity =  np.multiply(tmp2,C_kms)
    return velocity

#WL over which to normalise the spectra
NORM_WL = [6600,6700]

O3_LINE_WL = [6294.0,6306.0] # O3 [6285,6310]
O3_LINE_centre = [6300.0] # O3

Halpha_LINE_WL = [6540.0,6580.0]  # Halpha [6520,6600]
Halpha_LINE_centre = [6562.81]

He_LINE_WL = [6676.0,6698.0]  # He [6665,6690]
He_LINE_centre = [6678.2]


#############################################################
# Choose line
#######################################
line = 'O3'
if line == 'he':
    LINE_WL = He_LINE_WL
    LINE_centre = He_LINE_centre
    spectral_offset = 0.5
    plot_title ="HeI 6678.2 $\AA$"
    annotate_shift = 0.1
elif line == 'halpha':
    LINE_WL = Halpha_LINE_WL
    LINE_centre = Halpha_LINE_centre
    plot_title ="H$\alpha$ 6562.81 $\AA$"
    annotate_shift = 0.5
elif line == 'O3':
    LINE_WL = O3_LINE_WL
    LINE_centre = O3_LINE_centre
    spectral_offset = 1.2
    plot_title ="[OI] 6300 $\AA$"
    annotate_shift = 0.2



#Number of different files to test EW parameters in
list_of_lists = ['JL13BO1_extraction.txt','JL13A08_extraction.txt','JL14AO4_extraction.txt','JL14B02_extraction.txt']


#This loop is to read in the list of files and find the dates of observations. 
#These are then stored in an array, along with the file names
obs_day_mjd = [] # Initiasing arrays. NOTE: here the obs_mjd is rounded to the nearest full number
file_name_list = []
for epoch in range(0,len(list_of_lists)):
    fits_list  = np.genfromtxt(list_of_lists[epoch],dtype='S40')
    
    fits_info  = np.genfromtxt(list_of_lists[epoch],dtype=[('file_name','S40'),('DATE-OBS','S40'),('SN','f4'),('no_pix','i2')])

    for file_loop in range(0,len(fits_info)):

        file_name = fits_info['file_name'][file_loop]
        data = Datacube(file_name,2)
        obs_day_mjd.append(np.around(data.return_MJD())) #
        file_name_list.append(file_name)
        


first_day = np.min(obs_day_mjd) # Set the first day of observations

set_obs_day_mjd = set(obs_day_mjd) #Find the number of unique days observations
unique_obs_day_mjd = list(set_obs_day_mjd)
unique_obs_day_mjd.sort()
print "Number of unique MJDs to nearest in",len(unique_obs_day_mjd)

#Initialise sum spectrum array. 
data.get_cube()
data.get_cube_wl()
median_spectra = np.ndarray(shape=(len(data.wl),len(unique_obs_day_mjd)), dtype=float, order='F')

# To check which spectra have the wrong length wl arrays
wrong_spec = []
wrong_spec_date = []


###################
#setting up plot
plt.figure(num=None, figsize=(8, 16), dpi=80, facecolor='w', edgecolor='k')


for date_loop in range(0,len(unique_obs_day_mjd)):

    
    date = unique_obs_day_mjd[date_loop]
    day = date - first_day #Day starts 0 at first observation
    #Find the files that were observed on this day.
    print "Date Loop", date_loop,"date",date,"day",day
    date_indexes = [index for index,x in enumerate(obs_day_mjd) if x == date]
    print date_indexes
    source_files = [file_name_list[index] for index in date_indexes]

    spectra4date = [] # Initilise loop to store all spectra for that date

    for source_loop in range(0,len(source_files)):
        source = source_files[source_loop]
        pipe_spec = SpecCube(source,5)
        pipe_spec.normalise_spec(NORM_WL)
        shift = pipe_spec.bary_shift_correction()

        if len(pipe_spec.wl) < 3000:
            print "The Wl array is too short?, ",source," ",date
            print "Length of wl array",len(pipe_spec.wl)
            print "Length of flux array",len(pipe_spec.flux)
            print "Date of short spectrum",date
            wrong_spec.append(source)
            wrong_spec_date.append(date)

        f = interp1d(pipe_spec.wl_corr.value,pipe_spec.norm_flux,bounds_error=False, fill_value=0.0)
        resampled_flux = f(pipe_spec.wl) 
        print len(pipe_spec.wl)
       
        spectra4date.append(resampled_flux)



    #Find the median spectrum of all the osberved spectra for that day
    median_spectra = spectra4date[0]*0

    for spectral_ele in range(0,len(spectra4date[1])):

        median_spectra[spectral_ele] = np.median([spectra4date[spectra][spectral_ele] for spectra in range(0,len(spectra4date))])  

    velocity_array = convert_wl2vel(pipe_spec.wl,LINE_centre)

    offset_all_spectra = np.add(spectra4date,spectral_offset*date_loop)
    plt.fill_between(velocity_array,offset_all_spectra[0],offset_all_spectra[1],facecolor='red',alpha=0.1)
    plt.fill_between(velocity_array,offset_all_spectra[1],offset_all_spectra[2],facecolor='red',alpha=0.1)
    plt.fill_between(velocity_array,offset_all_spectra[2],offset_all_spectra[0],facecolor='red',alpha=0.1)

    offset = median_spectra + spectral_offset*date_loop
    mean_offset = np.median(offset)+annotate_shift
    x_pos = LINE_centre[0]-10



    plt.plot(velocity_array,offset,label=int(day),color='red')
    plt.annotate(time_stamp[date_loop],xy=(vel_at_wl(LINE_centre,x_pos),mean_offset))


plt.ylim(0,25)
print 'line centre',LINE_centre,"LINE_WL[0]",LINE_WL[0],"LINE_WL[1]",LINE_WL[1]

plt.xlim(vel_at_wl(LINE_centre,LINE_WL[0]),vel_at_wl(LINE_centre,LINE_WL[1]))

#plt.legend(loc=4)
plt.xlim(-600,600)
plt.xlabel("Velocity [km/s]")
plt.ylabel("Normalised Flux + Offset")
plt.title(plot_title)
plt.show()


