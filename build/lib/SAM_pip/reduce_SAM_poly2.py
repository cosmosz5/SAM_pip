def weighted_avg_and_std(values, weights):
    import numpy as np
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))



def redSAM(data_file, PSF_im, uwindow, bl_x, bl_y, wave, hole_size, imsize, px_scale, nholes, nbl, ncp, \
           baselines, closure_phases, x, y, air, cube_bl, cube_bl_sin, xcoord, ycoord, source, bandwidth, instrument, \
           oversample, arrname = 'SIM'):
    import numpy as np
    import astropy.io.fits as pyfits
    import matplotlib.pyplot as plt
    from . import sam_tools
    import pdb
    import pandas as pd
    from . import js_oifits as oifits
    from astropy.time import Time
    from datetime import datetime
    from scipy.stats import trim_mean
    import scipy.ndimage as sciim

    data_cube = pyfits.open(data_file)
    data_header = data_cube[0].header
    data_cube = data_cube[0].data
    sz = np.shape(data_cube)
    data_cube = sciim.zoom(data_cube, (1, oversample, oversample), order=3)

    sz = np.shape(data_cube)
    V2_mod = np.zeros([sz[0],nbl])
    PHASE_mod = np.zeros([sz[0],nbl])
    windowed_data = np.zeros([sz[0], sz[1], sz[2]])
    MODEL_IM = np.zeros([sz[0],imsize, imsize])

    for i in range(sz[0]):

        windowed_PSF = PSF_im * uwindow
        windowed_data[i,:,:] = data_cube[i,:,:] * uwindow
        ind = np.where(uwindow != 0.0)
        BB = windowed_data[i,ind[0], ind[1]]

        #### Define the matrix for the P2VM:
        P2VM = np.zeros([BB.shape[0], 2*nbl+2])

        P2VM[:, 0] = air[ind] * nholes
        P2VM[:, -1] = 1.0

        f=0; m = 0
        for k in range(2*nbl):
            if (k %2 == 0):
                P2VM[:, k+1] = air[ind] * 2. *cube_bl[f, ind[0], ind[1]]
                f = f + 1
            elif (k%2 !=0):
                P2VM[:, k + 1] = air[ind] * 2. *cube_bl_sin[m, ind[0], ind[1]]
                m = m + 1

        P2VMt = np.linalg.pinv(P2VM)

        MOD = np.dot(P2VMt,BB)
        MOD_fringes = np.dot(P2VM,MOD)
        MODEL_IM[i,ind[0], ind[1]] = MOD_fringes

        for ll in range(nbl):
            V2_mod[i, ll] = (MOD[ll*2 +1]**2 + MOD[ll*2+2]**2) / MOD[0]**2
            PHASE_mod[i, ll] = np.arctan2(MOD[ll*2+2], MOD[ll*2+1])


    if arrname == 'SIM':
        mean_parang = 0.0
    else:
        ##### Transform the baselines according to the rotation of the parallactic angle:
        mean_parang = (data_header['HIERARCH ESO TEL PARANG START'] + data_header['HIERARCH ESO TEL PARANG END']) / 2.0

    bl_xp = bl_x * np.cos(np.deg2rad(mean_parang)) - bl_y * np.sin(np.deg2rad(mean_parang))
    bl_yp = bl_x * np.sin(np.deg2rad(mean_parang)) + bl_y * np.cos(np.deg2rad(mean_parang))
    bl_x1_cp, bl_y1_cp, bl_x2_cp, bl_y2_cp, t3amp_mod, t3phi_mod = sam_tools.compute_closure_phases(nbl, ncp, baselines, closure_phases, \
                                                                               bl_xp, bl_yp, np.sqrt(V2_mod), \
                                                                               np.rad2deg(PHASE_mod))

    ### Here we can do some statistics to get the errorbars of the observables:
    ### We used a bootstraping method to estimate the mean and standard deviation of the observables
    VIS_mod = np.sqrt(V2_mod)

    #fig1, (ax1,ax2) = plt.subplots(1, 2)
    VIS_aver = np.zeros([VIS_mod.shape[1]])
    VIS_aver_err = np.zeros([VIS_mod.shape[1]])
    PHASE_aver  = np.zeros([VIS_mod.shape[1]])
    PHASE_aver_err = np.zeros([VIS_mod.shape[1]])
    V2_aver = np.zeros([V2_mod.shape[1]])
    V2_aver_err = np.zeros([V2_mod.shape[1]])
    CP_aver = np.zeros([t3phi_mod.shape[1]])
    CP_aver_err = np.zeros([t3phi_mod.shape[1]])
    T3AMP_aver = np.zeros([t3phi_mod.shape[1]])
    T3AMP_aver_err = np.zeros([t3phi_mod.shape[1]])

    if sz[0] <= 3:
        VIS_aver = np.mean(VIS_mod, axis=0)
        VIS_aver_err = np.std(VIS_mod, axis=0)
        PHASE_aver = np.mean(PHASE_mod, axis=0)
        PHASE_aver_err = np.std(PHASE_mod, axis=0)
        V2_aver = np.mean(V2_mod, axis=0)
        V2_aver_err = np.std(V2_mod, axis=0)
        CP_aver = np.mean(t3phi_mod, axis=0)
        CP_aver_err = np.std(t3phi_mod, axis=0)
        T3AMP_aver = np.mean(t3amp_mod, axis=0)
        T3AMP_aver_err = np.std(t3amp_mod, axis=0)
    else:
        for i in range(VIS_mod.shape[1]):
            pf = pd.DataFrame({'VIS_mod': VIS_mod[:, i]})
            bootstrap = pd.DataFrame({'mean_VIS_aver': [
                pf.sample(int(len(VIS_mod[:, i]) / 3), replace=True).VIS_mod.mean() for j in range(1000)]})
            # hist = bootstrap.mean_VIS_aver.hist(bins=3)
            VIS_aver[i] = bootstrap.mean_VIS_aver.mean()
            VIS_aver_err[i] = bootstrap.mean_VIS_aver.std()

            pf_phi = pd.DataFrame({'PHASE_mod': PHASE_mod[:, i]})
            bootstrap_phi = pd.DataFrame({'mean_PHASE_aver': [
                pf_phi.sample(int(len(PHASE_mod[:, i]) / 3), replace=True).PHASE_mod.mean() for j in range(1000)]})
            PHASE_aver[i] = bootstrap_phi.mean_PHASE_aver.mean()
            PHASE_aver_err[i] = bootstrap_phi.mean_PHASE_aver.std()

        for i in range(V2_mod.shape[1]):
            pf = pd.DataFrame({'V2_mod': V2_mod[:, i]})
            bootstrap = pd.DataFrame({'mean_V2_aver': [pf.sample(int(len(V2_mod[:, i]) / 3), replace=True).V2_mod.mean()
                                                       for j in range(1000)]})
            V2_aver[i] = bootstrap.mean_V2_aver.mean()
            V2_aver_err[i] = bootstrap.mean_V2_aver.std()

        for i in range(t3phi_mod.shape[1]):
            pf = pd.DataFrame({'t3phi_mod': t3phi_mod[:, i]})
            bootstrap = pd.DataFrame({'mean_t3phi_aver': [
                pf.sample(int(len(t3phi_mod[:, i]) / 3), replace=True).t3phi_mod.mean() for j in range(1000)]})
            CP_aver[i] = bootstrap.mean_t3phi_aver.mean()
            CP_aver_err[i] = bootstrap.mean_t3phi_aver.std()

            pf_t3amp = pd.DataFrame({'t3amp_mod': t3amp_mod[:, i]})
            bootstrap_t3amp = pd.DataFrame({'mean_t3amp_aver': [
                pf_t3amp.sample(int(len(t3amp_mod[:, i]) / 3), replace=True).t3amp_mod.mean() for j in range(1000)]})
            T3AMP_aver[i] = bootstrap_t3amp.mean_t3amp_aver.mean()
            T3AMP_aver_err[i] = bootstrap_t3amp.mean_t3amp_aver.std()


    ###########################################
    ####### Plot the visibilities and closure phases versus spatial frequency

    sam_tools.plot_v2_cp(source, instrument, nbl, ncp, wave, bl_xp, bl_yp, bl_x1_cp, bl_y1_cp, bl_x2_cp, bl_y2_cp,\
                         V2_aver, V2_aver_err, CP_aver, CP_aver_err)

    ####### Save the OIFITS file ##############

    oi_file = oifits.oifits()
    if arrname == 'SIM':
        oi_arrname = instrument+'_'+arrname
        oi_target = source+'_'+arrname
        oi_ra = 0.0
        oi_dec = 0.0
        oi_pmra = 0.0
        oi_pmdec = 0.0
        oi_time = np.zeros([nbl]) + 2000.0
        tjd = Time(datetime.now().strftime('%Y-%m-%d'))
        oi_dateobs = datetime.now().strftime('%Y-%m-%d')
        oi_mjd = np.zeros([nbl]) + tjd.mjd
        oi_inttime = np.zeros([nbl]) + 10.0
        oi_visflag = np.zeros([nbl], dtype='i1')
        oi_v2flag = np.zeros([nbl], dtype='i1')
        oi_targetid = np.zeros([nbl]) + 1
        oi_t3targetid = np.zeros([ncp]) + 1
        oi_t3flag = np.zeros([ncp], dtype='i1')
        oi_t3inttime = np.zeros([ncp]) + 10.0
        oi_t3time = np.zeros([ncp]) + 2000.0
        oi_t3mjd = np.zeros([ncp]) + tjd.mjd
        catg = 'SIM_DATA'
        equinox = 2000.0
        radvel = 0.0
        parallax = 0.0

    else:
        oi_arrname = instrument + '_' + arrname
        oi_target = source + '_' + arrname
        oi_ra = data_header['HIERARCH ESO TEL TARG ALPHA']
        oi_dec = data_header['HIERARCH ESO TEL TARG DELTA']
        oi_pmra = data_header['HIERARCH ESO TEL TARG PMA']
        oi_pmdec = data_header['HIERARCH ESO TEL TARG PMD']
        oi_time = np.zeros([nbl]) + data_header['UTC']
        oi_dateobs = data_header['DATE-OBS']
        oi_datenow = datetime.now().strftime('%Y-%m-%d')
        oi_mjd = np.zeros([nbl]) + data_header['MJD-OBS']
        if instrument == 'NACO':
            oi_inttime = np.zeros([nbl]) + data_header['HIERARCH ESO DET DIT']
        else:
            oi_inttime = np.zeros([nbl]) + data_header['HIERARCH ESO DET SEQ1 DIT']
        oi_visflag = np.zeros([nbl], dtype='i1')
        oi_v2flag = np.zeros([nbl], dtype='i1')
        oi_targetid = np.zeros([nbl]) + 1
        oi_t3targetid = np.zeros([ncp]) + 1
        oi_t3flag = np.zeros([ncp], dtype='i1')
        if instrument == 'NACO':
            oi_t3inttime = np.zeros([ncp]) + data_header['HIERARCH ESO DET DIT']
        else:
            oi_t3inttime = np.zeros([ncp]) + data_header['HIERARCH ESO DET SEQ1 DIT']
        oi_t3time = np.zeros([ncp]) + data_header['EQUINOX']
        oi_t3mjd = np.zeros([ncp]) + data_header['MJD-OBS']
        catg = data_header['HIERARCH ESO DPR CATG']
        equinox = data_header['HIERARCH ESO TEL TARG EQUINOX']
        radvel = data_header['HIERARCH ESO TEL TARG RADVEL']
        parallax = data_header['HIERARCH ESO TEL TARG PARALLAX']

    oi_telname = np.zeros([nholes], dtype='S16')
    oi_staname = np.zeros([nholes], dtype='S16')
    oi_staindex = np.zeros([nholes], dtype='>i2')
    oi_sta_coord = np.zeros([nholes, 3], dtype='>f8')
    oi_size = np.zeros([nholes], dtype='f4')
    for i in range(nholes):
        oi_telname[i] = 'Hole' + str(i + 1)
        oi_staname[i] = 'P' + str(i + 1)
        oi_staindex[i] = i + 1
        oi_size[i] = hole_size
        oi_sta_coord[i, 0] = xcoord[i, 0]
        oi_sta_coord[i, 1] = ycoord[i, 0]
        oi_sta_coord[i, 2] = 0

    oi_file.array = oifits.OI_ARRAY(1, oi_arrname, 'GEOCENTRIC', 0.,0.,0.,oi_telname, oi_staname, oi_staindex, oi_size, oi_sta_coord)

    oi_file.target = oifits.OI_TARGET(1, np.array([1]), np.array([oi_target]), np.array([oi_ra]), np.array([oi_dec]), \
                                      np.array([equinox]), np.array([0.]), np.array([0.]), np.array([radvel]), \
                                      np.array(['UNKNOWN']), np.array(['OPTICAL']), np.array([oi_pmra]), np.array([oi_pmdec]), \
                                      np.array([0.]), np.array([0.]), np.array([parallax]), \
                                      np.array([0.]), np.array(['UNKNOWN']))
    oi_file.wavelength = oifits.OI_WAVELENGTH(1,instrument, np.array([wave]), np.array([bandwidth]))

    oi_file.vis = oifits.OI_VIS(1, oi_dateobs, oi_arrname, instrument, oi_targetid, oi_time, oi_mjd, oi_inttime, VIS_aver, \
                                VIS_aver_err, PHASE_aver, PHASE_aver_err, bl_xp, bl_yp, baselines+1, oi_visflag)

    oi_file.vis2 = oifits.OI_VIS2(1, oi_dateobs, oi_arrname, instrument, oi_targetid, oi_time, oi_mjd, oi_inttime, V2_aver,\
                                  V2_aver_err, bl_xp, bl_yp, baselines+1, oi_v2flag)

    oi_file.t3 = oifits.OI_T3(1, oi_dateobs, oi_arrname, instrument, oi_t3targetid, oi_t3time, oi_t3mjd, oi_t3inttime, \
                              T3AMP_aver, T3AMP_aver_err, CP_aver, CP_aver_err, bl_x1_cp, bl_y1_cp, bl_x2_cp, bl_y2_cp,\
                              closure_phases+1, oi_t3flag)

    oi_file.write(catg+'_uncalib_'+data_file[:-5]+'.oifits')
    #MODEL_IM = sciim.zoom(MODEL_IM, (1, 1/oversample, 1/oversample), order=3)
    pyfits.writeto('MODEL_interferogram_'+data_file, MODEL_IM, overwrite=True)
    #windowed_data = sciim.zoom(windowed_data, (1 / oversample, 1 / oversample), order=3)
    pyfits.writeto('WINDDAT_'+data_file, windowed_data, overwrite=True)
    return
