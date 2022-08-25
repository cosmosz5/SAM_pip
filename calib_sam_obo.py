####### Import the necessary modules ############
import sim_SAM
import pdb
import calibrate_SAM_obo
import readcol
import os

input_sci_path = ''
input_cal_path = ''
calibrator_filename = 'cal_hd37093_480_oifits.txt'
science_filename = 'data_480_oifits.txt'
cal_name = 'HD37093'

[sc_files] = readcol.readcol(input_sci_path+science_filename, twod=False)
[cal_files] = readcol.readcol(input_cal_path+calibrator_filename, twod=False)
for i in range(len(sc_files)):

    sc_dat = sc_files[i]
    cal_dat = cal_files[i]
    delta = 1000.0  ### Hours
    calibrate_SAM_obo.calibrate_SAM(input_sci_path, input_cal_path, sc_dat, cal_dat, cal_name, delta)
    print('done')
