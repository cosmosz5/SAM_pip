####### Import the necessary modules ############
import pdb
import calibrate_SAM
import os

## For calibration
sc_files = 'data_380_oifits.txt'
cal_files = 'cal_hd37093_380_oifits.txt'
delta = 1000.0  ### Hours
calibrate_SAM.calibrate_SAM(sc_files, cal_files, delta)
