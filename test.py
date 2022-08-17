####### Import the necessary modules ############
import sim_SAM
import pdb
import calibrate_SAM
import os
import pylab
pylab.ion()

#######  The following parameters are necessary to run the SAM pipelie ##################
mask_filename = '7holes_jwst_mask_corr.txt' #A text file with the mask geometry
wave = 3.82650365041922e-06 #3.828e-06 #4.817e-6 #3.828e-06 #4.2934e-6 #Central wavelength of the observations (meters)
bandwidth = 0.205e-06 #1.93756431390019e-07 #0.205e-06 #0.298e-06 #0.205e-06 #0.202e-6 #Bandwidth of the observations (meters)
hole_size = 0.8 #in meters
imsize = 80 # in pixels
px_scale = 65.6 #in mas
hole_geom = 'HEXAGON' #HEXAGON for the JWST
inst = 'JWST' 
arrname = 'DATA' ## This could be DATA or SIM
rotation_angle = 0.0 ### In case the mask is not propely aligned with the position indicated in the manual
oversample = 1.0 ## Number of times that you want to oversample the data

#print('Data 380')
#data_filename = 'data_380.txt'
#source = 'AB-DOR'
#sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
#           arrname, rotation_angle, oversample)

#print('Cal1 380')
#data_filename = 'cal_hd37093_380.txt'
#source = 'HD-37093'
#sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
#           arrname, rotation_angle, oversample)

#print('Cal2 380')
#data_filename = 'cal_hd36085_380.txt'
#source = 'HD-36085'
#sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
#           arrname, rotation_angle, oversample)


wave = 4.285-6 #3.828e-06 #4.817e-6 #3.828e-06 #4.2934e-6 #Central wavelength of the observations (meters)
bandwidth = 0.202e-6 #1.93756431390019e-07 #0.205e-06 #0.298e-06 #0.205e-06 #0.202e-6 #Bandwidth of the observations (meters)

print('Data 430')
data_filename = 'data_430.txt'
source = 'AB-DOR'
sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)

print('Cal1 430')
data_filename = 'cal_hd37093_430.txt'
source = 'HD-37093'
sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)

print('Cal2 430')
data_filename = 'cal_hd36085_430.txt'
source = 'HD-36085'
sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)


wave = 4.817e-6 #3.828e-06 #4.817e-6 #3.828e-06 #4.2934e-6 #Central wavelength of the observations (meters)
bandwidth = 0.298e-06 #1.93756431390019e-07 #0.205e-06 #0.298e-06 #0.205e-06 #0.202e-6 #Bandwidth of the observations (meters)

print('Data 480')
data_filename = 'data_480.txt'
source = 'AB-DOR'
sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)

print('Cal1 480')
data_filename = 'cal_hd37093_480.txt'
source = 'HD-37093'
sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)

print('Cal2 480')
data_filename = 'cal_hd36085_480.txt'
source = 'HD-36085'
sim_SAM.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)
