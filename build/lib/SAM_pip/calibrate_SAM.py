
def calibrate_SAM(calibrator_filename, science_filename, delta=2.0):
    import numpy as np
    import astropy.io.fits as pyfits
    from .readcol import readcol as readcol
    import pandas as pd
    from . import js_oifits as oifits

    [sc_files] = readcol(calibrator_filename, twod=False)
    [cal_files] = readcol(science_filename, twod=False)

    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    for i in range(len(sc_files)):
        sc_data = pyfits.open(sc_files[i])
        sc_v = sc_data['OI_VIS'].data['VISAMP']
        sc_v_err = sc_data['OI_VIS'].data['VISAMPERR']
        sc_vphi = sc_data['OI_VIS'].data['VISPHI']
        sc_vphi_err = sc_data['OI_VIS'].data['VISPHIERR']
        sc_v2 = sc_data['OI_VIS2'].data['VIS2DATA']
        sc_v2_err = sc_data['OI_VIS2'].data['VIS2ERR']
        sc_t3phi = sc_data['OI_T3'].data['T3PHI']
        sc_t3phi_err = sc_data['OI_T3'].data['T3PHIERR']
        sc_mjd = sc_data['OI_VIS2'].data['MJD'][0]

        ### Define the arrays for the calibrators: ####
        cal_v = np.zeros([len(cal_files), len(sc_v)])
        cal_v_err = np.zeros([len(cal_files), len(sc_v)])
        cal_vphi = np.zeros([len(cal_files), len(sc_v)])
        cal_vphi_err = np.zeros([len(cal_files), len(sc_v)])
        cal_v2 = np.zeros([len(cal_files), len(sc_v2)])
        cal_v2_err = np.zeros([len(cal_files), len(sc_v2)])
        cal_t3phi = np.zeros([len(cal_files), len(sc_t3phi)])
        cal_t3phi_err = np.zeros([len(cal_files), len(sc_t3phi)])
        cal_mjd = np.zeros([len(cal_files)])

        for j in range(len(cal_files)):
            cal_data = pyfits.open(cal_files[j])
            cal_v[j, :] = cal_data['OI_VIS'].data['VISAMP']
            cal_v_err[j, :] = cal_data['OI_VIS'].data['VISAMPERR']
            cal_vphi[j, :] = cal_data['OI_VIS'].data['VISPHI']
            cal_vphi_err[j, :] = cal_data['OI_VIS'].data['VISPHIERR']
            cal_v2[j, :] = cal_data['OI_VIS2'].data['VIS2DATA']
            cal_v2_err[j, :] = cal_data['OI_VIS2'].data['VIS2ERR']
            cal_t3phi[j, :] = cal_data['OI_T3'].data['T3PHI']
            cal_t3phi_err[j, :] = cal_data['OI_T3'].data['T3PHIERR']
            cal_mjd[j] = cal_data['OI_VIS2'].data['MJD'][0]

        diff_MJD = (sc_mjd - cal_mjd) * 24.0  # Convert to hours
        [ind_MJD] = np.where(np.abs(diff_MJD) <= delta)

        aver_cal_vis = np.zeros([sc_v.shape[0]])
        aver_cal_vis_err = np.zeros([sc_v.shape[0]])
        aver_cal_visphi = np.zeros([sc_v.shape[0]])
        aver_cal_visphi_err = np.zeros([sc_v.shape[0]])
        aver_cal_vis2 = np.zeros([sc_v2.shape[0]])
        aver_cal_vis2_err = np.zeros([sc_v2.shape[0]])
        aver_cal_t3phi = np.zeros([sc_t3phi.shape[0]])
        aver_cal_t3phi_err = np.zeros([sc_t3phi.shape[0]])

        ##### To calibrate V2 #########
        for n in range(cal_v2.shape[1]):
            pf = pd.DataFrame({'CAL_V2_mod': cal_v2[ind_MJD, n]})
            bootstrap = pd.DataFrame(
                {'mean_CAL_V2': [pf.sample(int(len(cal_v2[ind_MJD, n]) / 3), replace=True).CAL_V2_mod.mean() for j in
                                 range(1000)]})
            aver_cal_vis2[n] = bootstrap.mean_CAL_V2.mean()
            aver_cal_vis2_err[n] = bootstrap.mean_CAL_V2.std()
        calib_v2 = sc_v2 / aver_cal_vis2  #### CALIBRATED VIS2
        calib_v2_err = calib_v2 * (
            np.sqrt((sc_v2_err / sc_v2) ** 2 + (aver_cal_vis2_err / aver_cal_vis2) ** 2))  ### CALIBRATED VIS2 ERROR

        ##### To calibrate VIS and PHASES #########
        for n in range(cal_v.shape[1]):
            pf = pd.DataFrame({'CAL_V_mod': cal_v[ind_MJD, n]})
            bootstrap = pd.DataFrame(
                {'mean_CAL_V': [pf.sample(int(len(cal_v[ind_MJD, n]) / 3), replace=True).CAL_V_mod.mean() for j in
                                range(1000)]})
            aver_cal_vis[n] = bootstrap.mean_CAL_V.mean()
            aver_cal_vis_err[n] = bootstrap.mean_CAL_V.std()

            pf2 = pd.DataFrame({'CAL_VISPHI_mod': cal_vphi[ind_MJD, n]})
            bootstrap2 = pd.DataFrame(
                {'mean_CAL_VISPHI': [pf2.sample(int(len(cal_vphi[ind_MJD, n]) / 3), replace=True).CAL_VISPHI_mod.mean()
                                     for j in
                                     range(1000)]})
            aver_cal_visphi[n] = bootstrap2.mean_CAL_VISPHI.mean()
            aver_cal_visphi_err[n] = bootstrap2.mean_CAL_VISPHI.std()
        calib_v = sc_v / aver_cal_vis  #### CALIBRATED VIS
        calib_v_err = calib_v * (
            np.sqrt((sc_v_err / sc_v) ** 2 + (aver_cal_vis_err / aver_cal_vis) ** 2))  ### CALIBRATED VIS ERROR
        calib_visphi = sc_vphi - aver_cal_visphi  #### CALIBRATED  VISPHI
        calib_visphi_err = np.sqrt((sc_vphi_err) ** 2 + (aver_cal_visphi_err) ** 2)  ### CALIBRATED VISPHIERR

        ##### To calibrate closupre phases #######3
        for n in range(cal_t3phi.shape[1]):
            pf = pd.DataFrame({'CAL_T3PHI_mod': cal_t3phi[ind_MJD, n]})
            bootstrap = pd.DataFrame(
                {'mean_CAL_T3PHI': [pf.sample(int(len(cal_t3phi[ind_MJD, n]) / 3), replace=True).CAL_T3PHI_mod.mean()
                                    for j in
                                    range(1000)]})
            aver_cal_t3phi[n] = bootstrap.mean_CAL_T3PHI.mean()
            aver_cal_t3phi_err[n] = bootstrap.mean_CAL_T3PHI.std()
        calib_t3phi = sc_t3phi - aver_cal_t3phi  #### CALIBRATED  T3PHI
        calib_t3phi_err = np.sqrt((sc_t3phi_err) ** 2 + (aver_cal_t3phi_err) ** 2)  ### CALIBRATED T3PHIERR

        ##### TO SAVE THE CALIBRATED DATA INTO THE OIFITS FILE #########
        oi_file = oifits.open(sc_files[i])
        vis_file = oi_file.vis
        vis2_file = oi_file.vis2
        t3phi_file = oi_file.t3

        calib_oifile = oifits.oifits()
        calib_oifile.array = oi_file.array
        calib_oifile.target = oi_file.target
        calib_oifile.wavelength = oi_file.wavelength

        oi_file.vis = oifits.OI_VIS(1, vis_file.dateobs, vis_file.arrname, vis_file.insname, vis_file.target_id,
                                    vis_file.time, vis_file.mjd, \
                                    vis_file.int_time, calib_v, calib_v_err, calib_visphi, calib_visphi_err,
                                    vis_file.ucoord, \
                                    vis_file.vcoord, vis_file.sta_index, vis_file.flag)

        oi_file.vis2 = oifits.OI_VIS2(1, vis2_file.dateobs, vis2_file.arrname, vis2_file.insname, vis2_file.target_id,
                                      vis2_file.time, \
                                      vis2_file.mjd, vis2_file.int_time, calib_v2, calib_v2_err, vis2_file.ucoord, \
                                      vis2_file.vcoord, vis2_file.sta_index, vis2_file.flag)

        oi_file.t3 = oifits.OI_T3(1, t3phi_file.dateobs, t3phi_file.arrname, t3phi_file.insname, t3phi_file.target_id,
                                  t3phi_file.time, \
                                  t3phi_file.mjd, t3phi_file.int_time, t3phi_file.t3amp, t3phi_file.t3amperr,
                                  calib_t3phi, \
                                  calib_t3phi_err, t3phi_file.u1coord, t3phi_file.v1coord, t3phi_file.u2coord, \
                                  t3phi_file.v2coord, t3phi_file.sta_index, t3phi_file.flag)

        oi_file.write('CALIB_' + sc_files[i])

    return 0




