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

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

rcParams.update({'lines.linewidth':1.5})
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

rcParams.update({'font.size':12})

spec_quality = np.genfromtxt('quality_of_spectra.txt',dtype=[('Date','float'),('Flag','int'),('file_name','S80')])
matching_flag = []
matching_flag_MJD = []


he_ew_file = "pipe_heI_ew.lst"
he_EW_array = np.genfromtxt(he_ew_file,dtype=[('date',int),('EW',float),('minEW',float),('maxEW',float),('NoEW',int)],delimiter='&')



for date in he_EW_array['date']:
    ## Flag index 
    flag_index = [index for index,day in enumerate(np.rint(spec_quality['Date'])) if day == date]
    matching_flag.append(spec_quality['Flag'][flag_index])
    matching_flag_MJD.append(spec_quality['Date'][flag_index])


he_yerr_lower = np.subtract(he_EW_array['EW'],he_EW_array['minEW'])
he_yerr_upper = np.subtract(he_EW_array['maxEW'],he_EW_array['EW'])



fig, ax1 = plt.subplots()
ax1.errorbar(he_EW_array['date'],he_EW_array['EW'],yerr=[he_yerr_lower,he_yerr_upper],fmt='o')
ax1.yaxis.label.set_color('blue')
#ax1.plot(he_EW_array['date'],he_EW_array['EW'],'b')

flag_indices = [index for index,flag in enumerate(matching_flag) if np.sum(flag) >= 2]
#for each_index in flag_indices:
#    ax1.plot(Ha_EW_array['date'][each_index],Ha_EW_array['EW'][each_index],'.r')


ax1.set_xlabel('Date')


ax1.set_ylabel("He EW [$\AA$]")

new_labels=['01.2013','0.4.2013','07.2013','11.2013','02.2014','05.2014','08.2014','12.2014','03.2015']
new_labels=['08.01.2013','18.04.2013','27.07.2013','04.11.2013','12.02.2014','23.05.2014','31.08.2014','09.12.2014','03.03.2015']

ax1.set_xticklabels(new_labels, rotation=40, ha='right')



plt.show()





