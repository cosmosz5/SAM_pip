# SAM_pip

#To reduce the SAM data, as part of CASSINI, we developed a software called SAMPip. It is included in this repository. 
The software consist in a series of Python scripts to perform the reduction. 
We included a test.py script to load a simulated data as example for the reduction. 
Just clone the repository and type on the Terminal

python test.py

This will run SAMPip to create raw observables form the comissioning data on ABDor and its calibrators. The data filenames are loaded using txt files where
the filenames are defined for SCI or CAL. They can be loaded secuentialy into the test.py script. In case you want to change the filter, just uptadate the central
wavelength and the corresponding bandwith. 

The previous command will go through the different data sets included and it will produce a series of .fits files for quality check purposes 
(to view the .fits files, users can use the DS9 software, together with the final OIFITS files. The output files are the following ones:

1. CENTERED_*.fits files > These are the centered input files. SAMPip does a fine adjustment of the interferogram centroid in the middle of the pixel grid (this option is, at the moment deactivated, but the files are produced).
2. MODEL_interferogram_*.fits > These are the models of the interferograms produced by the code. They should be quite similar to the input data files. We can use the MODEL files to inspect (quickly) visually how good is our fringe modeling (i.e., visibility extraction) 
3. MODEL_interferogram_windowed_*.fits > These are the models of the interferograms produced by the code but they are windowed. Only the valid pixels for the model extraction have values different from zero. 
4. MODEL_residuals+*.fits > These .fits files show the residuals between our fringe model and the input data.
5. MTF.fits > This .fits file includes the Mutual Trsfer Funcion of the interferogram
6. PSF.fits > This .fits file includes the model of the interferogram produced solely by the geometry of the non-redundant mask. It is the IDEAL PSF
7. cube_bl.fits > This .fits file contains the cube of the different interference pattern produce by each pair of pin-holes in the non-redundant mask. 
8. hexa.fits > This file contains the diffraction pattern of the pin-hole geometry
9. hexagon.fits > This file contains the geometry of the pin-hole mask
10. WINDDAT_*.fits > These files are the windowed input data sets used for the model extraction
11. SIM_DATA_uncalib.*.fits > These are the raw OFITIS files produced from the data. They are in a standard OIFITS format and they include all the interferometric observables extracted with SAMPip. 

To calibrate the data use the script called calib_sam_obo.py
As input you need txt files to include the raw OITIFS files created in the previous step using test.py
This calib_sam_obo.py will create the calibrated observables and it wil store them on OIFITS files with the prefix CALIB_*

You can combine these CALIB_*.fits files into a single OIFITS file to merge an epoch of observation. For this, you can use the script called oi_combine.py
Again, you need a txt file with the CALIB_*.fits files that you want to merge and to define an output filename. 

