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


ha_ew_file = "pipe_ha_ew.lst"
Ha_EW_array = np.genfromtxt(ha_ew_file,dtype=[('date',int),('EW',float),('minEW',float),('maxEW',float),('NoEW',int)])


oI_ew_file = "pipe_oI_ew.lst"
oI_EW_array = np.genfromtxt(oI_ew_file,dtype=[('date',int),('EW',float),('minEW',float),('maxEW',float),('NoEW',int)])


for date in Ha_EW_array['date']:
    ## Flag index 
    flag_index = [index for index,day in enumerate(np.rint(spec_quality['Date'])) if day == date]
    matching_flag.append(spec_quality['Flag'][flag_index])
    matching_flag_MJD.append(spec_quality['Date'][flag_index])


Ha_yerr_lower = np.subtract(Ha_EW_array['EW'],Ha_EW_array['minEW'])
Ha_yerr_upper = np.subtract(Ha_EW_array['maxEW'],Ha_EW_array['EW'])

oI_yerr_lower = np.subtract(oI_EW_array['EW'],oI_EW_array['minEW'])
oI_yerr_upper = np.subtract(oI_EW_array['maxEW'],oI_EW_array['EW'])


fig, ax1 = plt.subplots()
ax1.errorbar(Ha_EW_array['date'],Ha_EW_array['EW'],yerr=[Ha_yerr_lower,Ha_yerr_upper],fmt='o')
ax1.yaxis.label.set_color('blue')
ax1.plot(Ha_EW_array['date'],Ha_EW_array['EW'],'--b')

flag_indices = [index for index,flag in enumerate(matching_flag) if np.sum(flag) >= 2]
#for each_index in flag_indices:
#    ax1.plot(Ha_EW_array['date'][each_index],Ha_EW_array['EW'][each_index],'.r')

ax2 = ax1.twinx()
ax2.yaxis.label.set_color('green')
ax2.errorbar(oI_EW_array['date'],oI_EW_array['EW'],yerr=[oI_yerr_lower,oI_yerr_upper],color='g',fmt='^')
#for each_index in flag_indices:
#    ax2.plot(oI_EW_array['date'][each_index],oI_EW_array['EW'][each_index],'.r')
ax2.plot(oI_EW_array['date'],oI_EW_array['EW'],'--g')



ax1.set_xlabel('Date')

ax2.set_ylabel('[OI] EW [$\AA$]')

ax1.set_ylabel("Ha EW [$\AA$]")

new_labels=['01.2013','0.4.2013','07.2013','11.2013','02.2014','05.2014','08.2014','12.2014','03.2015']
new_labels=['08.01.2013','18.04.2013','27.07.2013','04.11.2013','12.02.2014','23.05.2014','31.08.2014','09.12.2014','03.03.2015']

ax1.set_xticklabels(new_labels, rotation=40, ha='right')



plt.show()





