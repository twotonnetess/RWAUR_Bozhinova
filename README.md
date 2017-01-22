# RWAUR_Bozhinova
===================



Introduction
-------------

This scripts were written to run on FRODOSpec. Further information on the instrument can be found here http://telescope.livjm.ac.uk/TelInst/Inst/FRODOspec/

One problem we encountered with this data is that it is a robotic telescope, so for some observations the star was not centered in the field of view and in some others they star was not in the field of view at all. 

This first step in working with this data was to pick those data that were of usable quality and also test for signal loses in the cases that the star was not centered in the field of view. 

In order to do this, I wrote my own reduction script (all of the functions are within Frodospec_class.py)), details of this script and what it entailed can be found in **spectrum_report.pdf.** 

Scripts
------------
The following scripts are a selection of those used in the analysis of the FRODOSpec observations of RWAUR. 
- *Frodospec_class.py* is the class written to reduce the IFU cubes, extract the spectra and analyse the resulting spectra. 

- *measure_ew.py* measures the equivalent width of the Halpha emission in the reduced spectra. 

- *measre_oII_ew.py* measures the equivalent width of the OI emisison line at 6300 angstrom. 

- *make_mean_profile.py* This script plots the mean line profile of either the Halpha, OI or HeI lines. 

- *make_mean_profile_Ha_split.py* This script plots the mean profile of the Halpha line, but in this case it is split into a mean profile before and after the dimming event. 

- *plot_spectra.py* This script plots a number of things beginning with a flattened IFU cube showing the peak of the emission, the extraction map, the S/N ratio in the spectrum versus number of pixels usedin the extraction (the turn over was used to determine the most suitable extraction map). Finally the extracted spectrum is plotted along side the extracted sky spectrum. Finally a note is recorded whether it was a well extracted spectrum or not. These plots can be found in *plots/* and are split into folders based on the observing date, JL13A08, JL13BO1, JL14AO4 and JL14B02.

 
