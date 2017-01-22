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

    plt.show()
#    plt.close()
    
pipe_EW = []
extracted_EW = []
obs_mjd = []


#Number of different files to test EW parameters in
list_of_lists = ['JL13BO1_extraction.txt','JL13A08_extraction.txt','JL14AO4_extraction.txt','JL14B02_extraction.txt']

#WL over which to normalise the spectra
NORM_WL = [6600,6700]
LINE_WL = [6520,6600]


for epoch in range(0,len(list_of_lists)):
    fits_list  = np.genfromtxt(list_of_lists[epoch],dtype='S40')
    
    fits_info  = np.genfromtxt(list_of_lists[epoch],dtype=[('file_name','S40'),('DATE-OBS','S40'),('SN','f4'),('no_pix','i2')])
#    print type(fits_info)

    for file_loop in range(0,len(fits_info)):
#        print fits_info[file_loop]

        file_name = fits_info['file_name'][file_loop]
        pix4extraction = fits_info['no_pix'][file_loop]

        #Make the Frodosepc class, 
        #Read in the Data cube for Frodospec
        #Make sum of Data cube along wl direction
        data = Datacube(file_name,2)
        data.get_cube()
        data.get_sum_cube()

        #Append MJD to array
        obs_mjd.append(data.return_MJD())

        #Extract the spectrum and normalise
        spectrum,signal_noise,extract_pix_map = data.extract_spec_sum_pixels(pix4extraction,NORM_WL)
        data.normalise_spec(NORM_WL)

        #find the sky spectrum
        sky_spectrum = data.find_sky_pixels()

        # extract pipeline spectrum for comparison
        pipe_spec = SpecCube(file_name,5)
        pipe_spec.normalise_spec(NORM_WL)


        #Average EW parameter
        line_win =[6536.95,6599.36]
        cont_win = [6522.50,6534.89,6601.94177612,6632.2609158]


        pipe_EW.append(pipe_spec.ew(line_win,cont_win))
        extracted_EW.append(data.ew(line_win,cont_win))

        plt.plot(pipe_spec.wl,pipe_spec.norm_flux)
        plt.plot(data.wl,data.norm_flux)
        plt.plot

        plt.plot([line_win[0],line_win[0]],[0.9*np.min(pipe_spec.norm_flux),1.1*np.max(pipe_spec.norm_flux)])
        plt.plot([line_win[1],line_win[1]],[0.9*np.min(pipe_spec.norm_flux),1.1*np.max(pipe_spec.norm_flux)])
        plt.plot([cont_win[0],cont_win[1]],[1,1])
        plt.plot([cont_win[2],cont_win[3]],[1,1])
        plt.xlim([cont_win[0],cont_win[3]])
        plt.ylim([0.9*np.min(pipe_spec.norm_flux),1.1*np.max(pipe_spec.norm_flux)])

        plt.savefig('plots/'+file_name+'_EW.png',bbox_inches='tight')
        plt.close()

plt.scatter(obs_mjd,pipe_EW,marker='.',label='Pipe EW',s=80)
plt.scatter(obs_mjd,extracted_EW,marker='^',label='Extracted EW',s=35)
plt.legend()
plt.xlabel('MJD')
plt.ylabel('EW $\AA$')
plt.show()


plt.hist(obs_mjd)
plt.show()


#EW and errors

spec_quality = np.genfromtxt('quality_of_spectra.txt',dtype=[('Date','float'),('Flag','int'),('file_name','S80')])
print spec_quality
matching_flag = []
matching_flag_MJD = []

print 'length of spec quality',
pipe_file_out = "pipe_ha_ew.lst"

open(pipe_file_out, 'w').close()

round_obs_mjd = np.rint(obs_mjd)
set_round_obs_mjd = list(set(round_obs_mjd))
set_round_obs_mjd.sort()
mean_EW_list = []
min_EW_list = []
max_EW_list = []
EW_date_list = []
for i in range(0,len(set_round_obs_mjd)):
    date = set_round_obs_mjd[i]
    print "Date",date

    indices = [index for index,day in enumerate(round_obs_mjd) if day == date]
        
    EWs = []
    for index in indices:
        EWs.append(pipe_EW[index])

    print EWs
    max_EW = np.max(EWs)
    min_EW = np.min(EWs)
    mean_EW = np.mean(EWs)
    len_EW_list = len(EWs)
    print EWs

    #Write results to file
    with open(pipe_file_out, "a") as myfile:
#        myfile.write(str(int(date))+' &   '+str(round(mean_EW,2))+' &   '+str(round(min_EW,2))+' &  '+str(round(max_EW,2))+'  &  '+str(int(len_EW_list))+'\n')
        myfile.write(str(int(date))+'    '+str(round(mean_EW,2))+'    '+str(round(min_EW,2))+'    '+str(round(max_EW,2))+'     '+str(int(len_EW_list))+'\n')
    mean_EW_list.append(round(mean_EW,2))
    max_EW_list.append(round(max_EW,2))
    min_EW_list.append(round(min_EW,2))
    EW_date_list.append(round(date,2))


    ## Flag index 
    flag_index = [index for index,day in enumerate(np.rint(spec_quality['Date'])) if day == date]
    matching_flag.append(spec_quality['Flag'][flag_index])
    matching_flag_MJD.append(spec_quality['Date'][flag_index])


yerr_lower = np.subtract(mean_EW_list,min_EW_list)
yerr_upper = np.subtract(max_EW_list,mean_EW_list)

plt.errorbar(EW_date_list,mean_EW_list,yerr=[yerr_lower,yerr_upper],fmt='o')

flag_indices = [index for index,flag in enumerate(matching_flag) if np.sum(flag) >= 2]
for each_index in flag_indices:
    plt.plot(EW_date_list[each_index],mean_EW_list[each_index],'.r')
plt.ylabel("EW [$\AA$]")
plt.xlabel("Time [MJD]")
plt.title('Ha  EW')
plt.show()

