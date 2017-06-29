import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
##Plotting
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
#Data class
from Frodospec_class import *

def plot_cube(some_cube,string):
    plt.pcolormesh(some_cube,cmap='PuBuGn')
    plt.colorbar()
    plt.title(string)


def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def multiplot(some_cube,some_spec,some_s2n,extract_map,extra_spec,line_wl,file_name):
    
    # some cube is frodo data cube class, where some spec comes from
    # some sn is from this frodo cube class as well
    # extract map is the map used to extract spectrum
    # extra spec is some extrac spec to overplot
    # Here line_wl is the window over which to plot line

    fig = plt.figure()


    gs = gridspec.GridSpec(2, 3) 
    gs.update(hspace=0.05,wspace=0.05)

    ax_cubesum = plt.subplot(gs[0,0])
    ax_extract = plt.subplot(gs[0,1])
    ax_sn = plt.subplot(gs[0,2])
    ax_spec = plt.subplot(gs[-1,:])

    rescaled_sum_cube = some_cube.cube_sum/100000

    im = ax_cubesum.pcolormesh(rescaled_sum_cube,cmap='PuBuGn')
    ax_cubesum.set_xlim(0,11)
    ax_cubesum.set_ylim(0,11)
    ax_cubesum.set_xticklabels([])
    ax_cubesum.set_yticklabels([])
    ax_cubesum.set_title("Cube Sum/100000")
    plt.colorbar(im,ax=ax_cubesum)



    im = ax_extract.pcolormesh(extract_map,cmap='Greys')
    ax_extract.set_xlim(0,11)
    ax_extract.set_ylim(0,11)
    ax_extract.set_xticklabels([])
    ax_extract.set_yticklabels([])
    ax_extract.set_title("Extraction Map")

    ax_sn.plot(range(0,len(some_s2n)),some_s2n)
    ax_sn.yaxis.tick_right()
    ax_sn.xaxis.tick_top()
    ax_sn.annotate('S/N Ratio',xy=(0.91,0.6),xycoords='axes fraction',rotation='vertical')
    ax_sn.annotate('No. Pixels',xy=(0.3,0.05),xycoords='axes fraction')


    ax_spec.plot(some_cube.wl,some_spec,'r')
    ax_spec.plot(some_cube.wl,extra_spec,'-b')
    ax_spec.set_xlim(min(some_cube.wl),max(some_cube.wl))
    ax_spec.set_xlabel("Wavelength [Angstrom]")
    ax_spec.set_ylabel("Counts")
 
    ax_spec.annotate('Extracted Spec.',xy=(0.05,0.6),xycoords='axes fraction',color='red')
    ax_spec.annotate('Pipeline Spec.',xy=(0.05,0.5),xycoords='axes fraction',color='blue')
   
    if np.max(extra_spec) > np.max(some_spec):
        ax_spec.set_ylim([0.9*np.min(extra_spec),1.1*np.max(extra_spec)])
    else:
        ax_spec.set_ylim([0.9*np.min(some_spec),1.1*np.max(some_spec)])

    line_inset = add_subplot_axes(ax_spec,[0.5,0.3,0.35,0.65])
    line_inset.plot(some_cube.wl,some_spec,'r')
    line_inset.plot(some_cube.wl,extra_spec,'b')

    line_inset.set_xlim(line_wl)
    wl_indices = [index for index,wl in enumerate(some_cube.wl) if line_wl[0] <= wl <= line_wl[1]]
    if np.max(extra_spec[wl_indices]) > np.max(some_spec[wl_indices]):
        line_inset.set_ylim([0.9*np.min(extra_spec[wl_indices]),1.1*np.max(extra_spec[wl_indices])])
    else:
        line_inset.set_ylim([0.9*np.min(some_spec[wl_indices]),1.1*np.max(some_spec[wl_indices])])

#    fig.savefig('plots/'+file_name+'.png',bbox_inches='tight')

#    plt.show()
    plt.close()
    

date = []
quality_flag = []
file_name_array = []

#Number of different epochs of observations
list_of_lists = ['J13BO1.lst','JL13A08.lst','JL14AO4.lst','JL14B02.lst']

#WL over which to normalise the spectra
NORM_WL = [6600,6700]
LINE_WL = [6520,6600]

for epoch in range(0,len(list_of_lists)):
    fits_list  = np.genfromtxt(list_of_lists[epoch],dtype='S40')

    for file_loop in range(0,len(fits_list)):
        print "File Name",fits_list[file_loop]
    
        file_name = fits_list[file_loop]

        data = Datacube(file_name,2)
        data.get_cube()
        data.get_sum_cube()

        spectrum,signal_noise,extract_pix_map = data.extract_spec_sum_pixels(50,NORM_WL)
        data.normalise_spec(NORM_WL)
    
        sky_spectrum = data.find_sky_pixels()

        pipe_spec = SpecCube(file_name,5)
    
        pipe_spec.normalise_spec(NORM_WL)

        multiplot(data,data.norm_flux,signal_noise,extract_pix_map,pipe_spec.norm_flux,LINE_WL,file_name)

	date.append(data.return_MJD())
#        quality = raw_input("Well extracted (1) Yes (2) No : ")
#        quality_flag.append(quality)
        file_name_array.append(file_name)


quality = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
print 'Length of quality flag',len(quality)
print 'length of date',len(date)

data_array = np.array(zip(date,quality,file_name_array))
print data_array
np.savetxt('quality_of_spectra.txt',data_array,fmt=['%s']+['%s']+['%s'],delimiter='    ')

