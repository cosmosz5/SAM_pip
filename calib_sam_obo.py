####### Import the necessary modules ############
import sim_SAM
import pdb
import calibrate_SAM_obo
import readcol
import os

input_sci_path = ''
input_cal_path = ''
calibrator_filename = 'cal_hd36085_380_oifits.txt' ### Calibrator filenames
science_filename = 'data_380_oifits.txt' ## Science filenames
cal_name = 'HD36085'

[sc_files] = readcol.readcol(input_sci_path+science_filename, twod=False)
[cal_files] = readcol.readcol(input_cal_path+calibrator_filename, twod=False)
for i in range(len(sc_files)):

    sc_dat = sc_files[i]
    cal_dat = cal_files[1]
    delta = 1000.0  ### The alrogithm takes into account the calibrators over a given period of time. I deactivate this using a large number to include all the calibrators needed.
    calibrate_SAM_obo.calibrate_SAM(input_sci_path, input_cal_path, sc_dat, cal_dat, cal_name, delta)
    print('done')
