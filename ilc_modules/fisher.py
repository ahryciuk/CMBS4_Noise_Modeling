

    ###############################################Input Parameters#############################################

def fisher(paramfile, which_spectra, 
           nl_dic, cl_residual, 
           freq_bands, ell, 
           fsky, cwd, 
           fwhm_bands, 
           f3db_beams=None, 
           pct_unc_lpf=None):
    
    import numpy as np, scipy as sc, sys, os
    import ilc, misc
    import json
    os.chdir(cwd + '/ilc_modules')
    from Systematics import systematics
    import tools
    param_dict = misc.fn_get_param_dict('params.ini')
    os.chdir(cwd)
    
    """
    #CMB fisher forecasting
    
    pass in an nl_dic to add systematics to then generate a new cl_residual to 
    propogate thru fisher calculation
    """  

    assert which_spectra in ('unlensed_scalar', 'lensed_scalar'),"Input either unlensed_scalar or lensed_scalar"
    
    #two parameter files
    
    
    #paramfile
    #which_spectra choose from ['lensed_scalar', 'unlensed_scalar']
    #nl_dic foreground cleaned
    #fsky
    #cambfolder?
    #cambfname?
    #cl_dic?
    
    #change files to relative to be cwd
    
    #excel inputs different sheets?
    #-prior sheet
    #-calc parameter sheet? like min/max ell?
    #-desired param is one value?
    

    #nl_dic = get_nl_from_cl_residuals(cl_residual, ell)
    bl_dic = get_bl_dic(freq_bands, fwhm_bands, ell)
    
    #BUG: the f3db_beams and pct_unc_lpf is used to turn on beam_syst but must be none. When it is passed as None here an error is thrown
    #check if inputs are strings
    if type(freq_bands[0]) is str:
        freq_bands = json.loads(freq_bands)
    if type(fwhm_bands[0]) is str:
        fwhm_bands = json.loads(fwhm_bands)
    if type(f3db_beams[0]) is str:
        f3db_beams = json.loads(f3db_beams)
    if type(pct_unc_lpf) is str:
        pct_unc_lpf = float(pct_unc_lpf)  
    if type(fsky) is str:
        fsky = float(fsky)
        
    #get experiment specs
    print('\n\tget experiment specs and nl')
    pspectra_to_use = ['TT', 'EE', 'TE']
    
    #pass ell mins and maxes from nl_dic
    fr = freq_bands[0]
    
    ####FIX#########
    #min_l_temp = 10
    #max_l_temp = 5000
    #min_l_pol = 10
    #max_l_pol = 5000
    ################
    min_l_temp = ell[0]
    max_l_temp = ell[-1]
    min_l_pol = ell[0]
    max_l_pol = ell[-1]    
    
    #####FIX#######
    #more inputs
    fix_params = [] #['mnu', 'ws'] but curretnly nothing to fix as we only have a 6+1(neff) LCDM model
    prior_dic = {'tau':0.007}#02}#02} #Planck tau prior
    desired_param = 'Neff' #'neff' #desired parameter for which we are computing the constraints. set to None if you want to analyse the full fisher matrix
    include_gal = 0
    gal_mask = 3 #only valid if galaxy is included
    #fsky = 0.77
    ##############

    ############################################################################################################


    ############################################################################################################
    #folders containing inputs
    camb_folder = '/ilc_modules/data/CMB_spectra_derivatives_for_code_comparison/'
    draft_results_folder = '/ilc_modules/data/DRAFT_results_20200601/s4like_mask/TT-EE-TE/baseline/'
    ############################################################################################################
    ############################################################################################################
    ###########################################Power Spectral Generation########################################
    
    camb_results = calculate_cmb(lmax = max_l_temp-50)
    cl_dic = get_power_spectra(camb_results,which_spectra, ell)
    
    #cl_dic = {}
    #cl_dic['TT'] = cl_camb[:,0]
    #cl_dic['EE'] = cl_camb[:,1]
    #cl_dic['TE'] = cl_camb[:,3]
    
    ############################################################################################################
    ############################################################################################################
    ############################################################################################################
    #get fiducual LCDM power spectra computed using CAMB
    print('\n\tget/read fiducual LCDM power spectra computed using CAMB')
    #camb_fname = cwd + '%s/cmb_spectra_%s_Srini.txt' %(camb_folder, which_spectra.replace('_scalar',''))

    #cl_camb = np.loadtxt(camb_fname)
    #els = cl_camb[:,0]
    #cl_dic = {}
    #cl_dic['TT'] = cl_camb[:,1]
    #cl_dic['EE'] = cl_camb[:,2]
    #cl_dic['TE'] = cl_camb[:,3]

    ############################################################################################################
    #read derivatives
    print('\n\tget/read derivatives')
    #camb_deriv_fname = cwd + '%s/cmb_spectra_derivs_%s_Srini.npy' %(camb_folder, which_spectra.replace('_scalar',''))
    #cl_deriv_dic_tmp = np.load(camb_deriv_fname, allow_pickle = 1).item()
    #cl_deriv_dic = {}
    #param_names = []
    #for p in sorted( cl_deriv_dic_tmp ):
    #    if p == 'ell': continue
    #    cl_deriv_dic[p]={}
    #    cl_deriv_dic[p]['TT'] = cl_deriv_dic_tmp[p][0]
    #    cl_deriv_dic[p]['EE'] = cl_deriv_dic_tmp[p][1]
    #    cl_deriv_dic[p]['TE'] = cl_deriv_dic_tmp[p][2]
    #    param_names.append( p )
        
    
    #FIXME: Param step size dict should be part of the input to this module    
    param_names=['As','mnu','Neff','ns','ombh2','ombh2','tau','thetastar']
    param_step_size_dic = {'As':0.1e-9,'mnu':0.02,'Neff':0.08,'ns':0.01,'ombh2':0.0008,'omch2':0.003,'tau':0.02,'thetastar':0.00005}    
    cl_deriv_dic = get_cmb_derivs_dic(ell, which_spectra, param_step_size_dic)    
    
    for param in param_names:
        if param == 'As':
            for TP in ['TT','EE','TE']:
                cl_deriv_dic[param][TP] = 1e-9*cl_deriv_dic[param][TP]
    ############################################################################################################
    #time constant transfer function error
    
    
    #Now want such that we can quantify the bias in the shifted peak due to this time constant error#
    #need to run with and without the error in the gaussian beam and subtract the peak bandpower values
    #in the measured spectra#
    
    syst = None
    if fwhm_bands is not None and f3db_beams is not None and pct_unc_lpf is not None:
        
        #initialize a systematics object
        syst = systematics(fsky, freq_bands, fwhm_bands, f3db_beams, ell)

        #calculate the systematic for all frequencies
        beam_syst = get_beam_systematics(bl_dic, fsky, freq_bands, fwhm_bands, f3db_beams, ell, pct_unc_lpf)

        #the beam systematics should be uncorrelated across frequency bands so add in quadrature
        pct_diff_avg = np.zeros(len(ell))
        for freq in beam_syst.keys():
            pct_diff_avg += beam_syst[freq]**2.
        pct_diff = np.sqrt(pct_diff_avg)

        #output power spectra will either be enhanced or suppressed due to uncertainty in detector time constant
        #BUG: The two arrays need to be the same length. Only multiply up to a point?
        for TP in cl_dic.keys():
            cl_dic[TP] = cl_dic[TP] * (1. + pct_diff)
            for param in param_names:
                cl_deriv_dic[param][TP] = cl_deriv_dic[param][TP] * (1. + pct_diff)
            
        #new cl_dic with shifted peak due to pct_diff? so can subtract later
    
    ############################################################################################################
    
    #add systematics to nl_dic here
    
    #TParr = ['T','P']
    
    ####FIX########
    #TParr = ['TT','EE','TE']
    #for TP in TParr:
    #    if TP == 'TT':
    #        nl_dic[TP] = nl_dic['T'][(fr,fr)]
    #   elif TP == 'EE':
    #        nl_dic[TP] = nl_dic['P'][(fr,fr)]
    #    else:
    #        nl_dic[TP] = np.zeros(len(nl_dic['T'][(fr,fr)]))
    ###############
    
    #parameters FIX:
    #readout
    NEI_list = [1,1,1,1,1,1]
    S_responsivity = [1e6,2e6,1e6,2e6,1e6,2e6]
    ell_knee_readout = [200,200,200,200,200,200]
    alpha_ind_readout = [1.5,1.5,1.5,1.5,1.5,1.5]

    #detector
    ell_knee_det = [150,150,150,200,200,200]
    alpha_ind_det = [2.2,2.2,2.2,2.2,2.2,2.2]
    white_noise_det = [5e-5,3e-5,4e-5,3e-5,5e-5,3e-5]

    #atmosphere
    ell_knee_atm = [150,150,150,150,150,150]
    alpha_ind_atm = [2.,2.1,2.1,2.1,2.1,2.1]
    white_noise_atm = [5e-4,2e-4,3e-4,2e-4,5e-4,3e-4]
    
    if syst is not None:
        readout_syst = get_readout_systematics(syst, ell, freq_bands, NEI_list, S_responsivity, ell_knee_readout, alpha_ind_readout)
        atm_syst = get_atm_systematics(syst, ell, freq_bands, alpha_ind_atm, ell_knee_atm, white_noise_atm)
        det_syst = get_det_systematics(syst, ell, freq_bands, alpha_ind_det, ell_knee_det, white_noise_det)
    
         #add systematics to noise spectrum
        nl_dic = add_systematics(syst, nl_dic, freq_bands, beam_syst=beam_syst, readout_syst=readout_syst, atm_syst=atm_syst, det_syst=det_syst)
            
         #calculate new cl_residual and proceed with fisher calc
   
       
    
    
    #up to here the nl_dic is normal nl spectra and needs to be converted to cl_residual
    
    #BUG: Only need this if adding systematics, the previous nl_dic is represented by the cl_residuals
    if nl_dic is not None:
        cl_residual = recalculate_cl_residual(freq_bands, ell, nl_dic, bl_dic)
    
    #convert nl_dic to foreground cleaned
    nl_dic = get_nl_from_cl_residuals(cl_residual, ell)
    
    
    ############################################################################################################
    
    #TODO: Build inclusion of galactic model
    
    #get nl
    #print('\n\tget nl')
    #if not include_gal:
    #    nlfile = cwd + '%s/S4_ilc_galaxy0_27-39-93-145-225-278_TT-EE-TE_lf2-mf12-hf5_AZ.npy' %(draft_results_folder)
    #else:
    #    nlfile = cwd + '%s/S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf2-mf12-hf5_galmask%s_AZ.npy' %(draft_results_folder, gal_mask)
    #nl_dic, fsky = tools.get_nldic(nlfile, ell)
    ############################################################################################################
    #get delta_cl
    print('\n\tget delta Cl')
    delta_cl_dic = tools.fn_delta_Cl(ell, cl_dic, nl_dic, fsky)
    ############################################################################################################
    #get Fisher / COV matrices
    print('\n\tget fisher')
    F_mat = tools.get_fisher_mat(ell, cl_deriv_dic, delta_cl_dic, param_names, pspectra_to_use = pspectra_to_use,\
                min_l_temp = min_l_temp, max_l_temp = max_l_temp, min_l_pol = min_l_pol, max_l_pol = max_l_pol)
    #print(F_mat)
    ############################################################################################################
    
    #Below is where Fisher matrix is refined (parameters are fixed)
    
    #fix params
    print('\n\tfixing paramaters, if need be')
    F_mat, param_names = tools.fn_fix_params(F_mat, param_names, fix_params)
    param_names = np.asarray(param_names)
    #print(F_mat)

    
    ############################################################################################################
    #add prior
    print('\n\tadding prior')
    F_mat = tools.fn_add_prior(F_mat, param_names, prior_dic)
    ############################################################################################################
    #get cov matrix now
    print('\n\tget covariance matrix')
    Cov_mat = sc.linalg.pinv2(F_mat) #made sure that COV_mat_l * Cinv_l ~= I
    ############################################################################################################
    #extract sigma(neff)
    if desired_param is not None:
        print('\n\textract sigma(neff)')
        pind = np.where(param_names == desired_param)[0][0]
        pcntr1, pcntr2 = pind, pind
        cov_inds_to_extract = [(pcntr1, pcntr1), (pcntr1, pcntr2), (pcntr2, pcntr1), (pcntr2, pcntr2)]
        cov_extract = np.asarray( [Cov_mat[ii] for ii in cov_inds_to_extract] ).reshape((2,2))
        sigma = cov_extract[0,0]**0.5
        sigma = round(sigma, 4)
        opline = '\n\t\t\simga(neff) = %.4f using observables = %s; fsky = %s; power spectra = %s\n' %(sigma, str(pspectra_to_use), fsky, which_spectra)
        print(opline)
    #sys.exit()
    
    else:
        for i in range(len(Cov_mat)):
            print(Cov_mat[i,i])
    
    return
############################################################################################################


def add_systematics(syst, nl_dic, freq_bands, beam_syst=None, readout_syst=None, atm_syst=None, det_syst=None):
    import numpy as np
    #from Systematics import systematics
    
    #initialize a systematics object
    #syst = systematics(fsky, freq_bands, fwhm_bands, f3db_beams, ell)
    
    TParr = ['TT','EE']
    #freq_bands
    
    if beam_syst is not None:
        
        for TP in TParr:
            for freq1 in freq_bands:
                for freq2 in freq_bands:
                    
                    if freq1 == freq2:
                        nl_dic[TP][(freq1,freq2)] = np.sqrt(nl_dic[TP][(freq1,freq2)]**2. + beam_syst[freq1]**2.)
                    
                    else:
                        nl_dic[TP][(freq1,freq2)] = np.sqrt(nl_dic[TP][(freq1,freq2)]**2. + beam_syst[freq1]**2. + beam_syst[freq2]**2.)
    
    if readout_syst is not None:
        
        for TP in TParr:
            for freq1 in freq_bands:
                for freq2 in freq_bands:
                    
                    if freq1 == freq2:
                        nl_dic[TP][(freq1,freq2)] = np.sqrt(nl_dic[TP][(freq1,freq2)]**2. + readout_syst[freq1]**2.)
                    
                    else:
                        nl_dic[TP][(freq1,freq2)] = np.sqrt(nl_dic[TP][(freq1,freq2)]**2. + readout_syst[freq1]**2. + readout_syst[freq2]**2.)
        
    if atm_syst is not None:
        
        for TP in TParr:
            for freq1 in freq_bands:
                for freq2 in freq_bands:
                    
                    if freq1 == freq2:
                        nl_dic[TP][(freq1,freq2)] = np.sqrt(nl_dic[TP][(freq1,freq2)]**2. + atm_syst[freq1]**2.)
                    
                    else:
                        nl_dic[TP][(freq1,freq2)] = np.sqrt(nl_dic[TP][(freq1,freq2)]**2. + atm_syst[freq1]**2. + atm_syst[freq2]**2.)
    
    if det_syst is not None:
        
        for TP in TParr:
            for freq1 in freq_bands:
                for freq2 in freq_bands:
                    
                    if freq1 == freq2:
                        nl_dic[TP][(freq1,freq2)] = np.sqrt(nl_dic[TP][(freq1,freq2)]**2. + det_syst[freq1]**2.)
                    
                    else:
                        nl_dic[TP][(freq1,freq2)] = np.sqrt(nl_dic[TP][(freq1,freq2)]**2. + det_syst[freq1]**2. + det_syst[freq2]**2.)
    
    return nl_dic

def get_beam_systematics(bl_dic, fsky, freq_bands, fwhm_bands, f3db_beams, ell, pct_unc_lpf):
    
    from Systematics import systematics
    
    #initialize a systematics object
    syst = systematics(fsky, freq_bands, fwhm_bands, f3db_beams, ell)
    
    pct_diff = syst.get_beam_unc(freq_bands, fwhm_bands, f3db_beams, pct_unc_lpf)
    
    beam_syst = {}
    for freq in freq_bands:
        beam_syst[freq] = bl_dic[freq]*pct_diff[freq]
    
    return beam_syst

def get_readout_systematics(syst, ell, freq_bands, NEI_list, S_responsivity, ell_knee_list, alpha_ind_list):
    
    #from Systematics import systematics
    
    #initialize a systematics object
    #syst = systematics(fsky, freq_bands, fwhm_bands, f3db_beams, ell)
    
    readout_syst = {}
    for freq in freq_bands:
        i = freq_bands.index(freq)
        readout_syst[freq] = syst.oneovereff_readout(ell, NEI_list[i], S_responsivity[i], ell_knee_list[i], alpha_ind_list[i])
    
    
    return readout_syst


def get_atm_systematics(syst, ell, freq_bands, alpha_ind_list, ell_knee_list, white_noise_level_list):
    
    #from Systematics import systematics
    
    """
    Generates atmosphere 1/f systematics spectra.
    
    Arguments:
    -------------------------------------------------------
    syst[]: Initialized object of Systematics class found in
            Systematics.py
    ell[]: array of angular scale values to generate 
            beams over
    freq_bands[]: list of band centers in GHz
    alpha_ind_list[]: list of 1/f indices for each frequency
            band
    ell_knee_list[]: list of ell knee values for each frequency
            band
    white_noise_level_list[]: list of white noise levels for
            each frequency band
            
    Output:
    -------------------------------------------------------
    atm_syst: dictionary of atmosphere systematic spectra organized
               by band centers
    """
    
    #initialize a systematics object
    #syst = systematics(fsky, freq_bands, fwhm_bands, f3db_beams, ell)
    
    atm_syst = {}
    for freq in freq_bands:
        i = freq_bands.index(freq)
        atm_syst[freq] = syst.oneovereff_atm(ell, alpha_ind_list[i], ell_knee_list[i], white_noise_level_list[i])
    
    
    return atm_syst

def get_det_systematics(syst, ell, freq_bands, alpha_ind_list, ell_knee_list, white_noise_level_list):
    
    #from Systematics import systematics
    
    """
    Generates detector 1/f systematics spectra.
    
    Arguments:
    -------------------------------------------------------
    syst[]: Initialized object of Systematics class found in
            Systematics.py
    ell[]: array of angular scale values to generate 
            beams over
    freq_bands[]: list of band centers in GHz
    alpha_ind_list[]: list of 1/f indices for each frequency
            band
    ell_knee_list[]: list of ell knee values for each frequency
            band
    white_noise_level_list[]: list of white noise levels for
            each frequency band
            
    Output:
    -------------------------------------------------------
    det_syst: dictionary of detector systematic spectra organized
               by band centers
    """
    
    #initialize a systematics object
    #syst = systematics(fsky, freq_bands, fwhm_bands, f3db_beams, ell)
    
    det_syst = {}
    for freq in freq_bands:
            i = freq_bands.index(freq)
            det_syst[freq] = syst.oneovereff_det(ell, alpha_ind_list[i], ell_knee_list[i], white_noise_level_list[i])
    
    return det_syst

def get_nl_from_cl_residuals(cl_residual, ell):
    import numpy as np
    
    """
    Extracts noise spectra in terms of TT, EE, and TE
    keys from foreground cleaned noise spectra cl_residual
    
    Arguments:
    -------------------------------------------------------
    cl_residual[]: dictionary of foreground cleaned noise
                    spectra
    ell[]: array of angular scale values
            
    Output:
    -------------------------------------------------------
    nl_dic[]: dictionary of noise spectra arrays organized
               by keys TT, EE, and TE
    """
    
    el_nl = ell
    
    if 'T' in cl_residual:
        nl_TT, nl_EE = cl_residual['T'], cl_residual['P']
        nl_TT = np.interp(ell, el_nl, nl_TT)
        nl_EE = np.interp(ell, el_nl, nl_EE)
        nl_TE = None
    else:
        nl_TT, nl_EE = cl_residual['TT'], cl_residual['EE']
        if 'TE' in cl_residual:
            nl_TE = cl_residual['TE']
        else:
            nl_TE = None
        el_nl = np.arange(len(nl_TT))
        nl_TT = np.interp(ell, el_nl, nl_TT)
        nl_EE = np.interp(ell, el_nl, nl_EE)
        if nl_TE is not None:
            nl_TE = np.interp(ell, el_nl, nl_TE)
        else:
            nl_TE = np.zeros(len(ell))

    nldic = {}
    nldic['TT'] = nl_TT
    nldic['EE'] = nl_EE
    nldic['TE'] = nl_TE
    
    
    return nldic

def get_bl_dic(freq_bands, beam_sizes, ell):
    
    import numpy as np
    
    """
    Creates dictionary of beam arrays as a function of ell.
    Returns a dictionary of Gaussian beams specified by a 
    full-width half max beam size in arcmin and band centers
    
    Arguments:
    -------------------------------------------------------
    freq_bands[]: array or list of band centers in GHz
    beam_sizes[]: array or list of full-width half max
                    beam sizes in arcmin
    ell[]: array of angular scale values to generate 
            beams over
            
    Output:
    -------------------------------------------------------
    bl_dic[]: dictionary of gaussian beam arrays organized
               by band centers
    """
    
    sigma = beam_sizes / np.sqrt(8. * np.log(2)) /60. * np.pi/180.
    
    bl_dic = {}
    for freq in freq_bands:
        i=freq_bands.index(freq)
        bl_dic[freq] = np.exp(ell * (ell + 1) * sigma[i]**2.)
    
    return bl_dic

def calculate_cmb(H0=None, ombh2=0.0222, omch2=0.1197, mnu=0.06, omk=0, 
                    tau=0.06, As=2.196e-9, ns=0.9655, r=0, lmax=5000, max_l_tensor=5000,
                    lens_potential_accuracy = 1, thetastar=0.010409,Alens=1, 
                    Neff = 3.046):
    import numpy as np
    import camb
    import os
    from camb import model, initialpower
    
    # initialize camb parameter, no tensor here
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk, 
                       tau=tau, nnu=Neff, thetastar=thetastar, Alens=Alens)
    pars.num_nu_massless = Neff
    #no tensor yet, set r=0 here, calculate inflation B mode separately
    pars.InitPower.set_params(As=As, ns=ns, r=0)
    pars.set_for_lmax(lmax, lens_potential_accuracy=lens_potential_accuracy)
    pars.max_l_tensor=max_l_tensor
    #calculate results for these parameters
    results = camb.get_results(pars)
    powers =results.get_cmb_power_spectra(pars, CMB_unit='muK', raw_cl=True)
    # return the parameter, results, and power
    return {'pars': pars, 'results': results, 'power': powers}
   
def get_power_spectra(camb_results, which_spectra, ell):
    import numpy as np
    
    cl_dic = {}
    #index = np.where(ell==lmax)[0][0]
    #ell = ell[:index]
    #get last accurate index of ell array
    #last_acc_index = np.where(ell == lmax-1)[0][0]
    #last_acc_index = lmax
    #cls = camb_results['power'][which_spectra][:last_acc_index,:]
    
    #FIXME: the ell vector thru all this should really get sorted out.
    #right now the lmin up to know is 2 with delta_ell = 1
    cls = camb_results['power'][which_spectra][2:,:]
    
    #dl_to_cl = np.zeros(len(ell)+2)
    #dl_to_cl[0] = 0
    #dl_to_cl[1] = np.pi
    #dl_to_cl[2:] = 2*np.pi / (ell * (ell + 1) )
    #lmax = ell[-1]
    
    #cl_dic['BB'] = cls[:,2]
    
    #NOTE: Starts from ell = 2
    #cl_dic['TT'] = cls[:lmax,0]# * dl_to_cl
    #cl_dic['EE'] = cls[:lmax,1]# * dl_to_cl
    #cl_dic['TE'] = cls[:lmax,3]# * dl_to_cl
    cl_dic['TT'] = cls[:,0]# * dl_to_cl
    cl_dic['EE'] = cls[:,1]# * dl_to_cl
    cl_dic['TE'] = cls[:,3]# * dl_to_cl
    
    
    return cl_dic

def get_cmb_derivs(ell, which_spectra, param, param_value, param_step_size):
    import numpy as np
    
    assert param in ('As', 'mnu', 'Neff', 'ns', 'ombh2', 'omch2', 'tau', 'thetastar'),"The parameter passed is not in the LCDM Model"
    
    cl_camb_deriv = {}
    camb_results_plus = 0
    camb_results_plus = 0
    lmax = ell[-1]
    
    if param is 'As':
        camb_results_plus = calculate_cmb(As = param_value+param_step_size, lmax=lmax)
        camb_results_minus = calculate_cmb(As = param_value-param_step_size, lmax=lmax)
    elif param is 'mnu':
        camb_results_plus = calculate_cmb(mnu = param_value+param_step_size, lmax=lmax)
        camb_results_minus = calculate_cmb(mnu = param_value-param_step_size, lmax=lmax)
    elif param is 'Neff':
        camb_results_plus = calculate_cmb(Neff = param_value+param_step_size, lmax=lmax)
        camb_results_minus = calculate_cmb(Neff = param_value-param_step_size, lmax=lmax)
    elif param is 'ns':
        camb_results_plus = calculate_cmb(ns = param_value+param_step_size, lmax=lmax)
        camb_results_minus = calculate_cmb(ns = param_value-param_step_size, lmax=lmax)
    elif param is 'ombh2':
        camb_results_plus = calculate_cmb(ombh2 = param_value+param_step_size, lmax=lmax)
        camb_results_minus = calculate_cmb(ombh2 = param_value-param_step_size, lmax=lmax)
    elif param is 'omch2':
        camb_results_plus = calculate_cmb(omch2 = param_value+param_step_size, lmax=lmax)
        camb_results_minus = calculate_cmb(omch2 = param_value-param_step_size, lmax=lmax)
    elif param is 'tau':
        camb_results_plus = calculate_cmb(tau = param_value+param_step_size, lmax=lmax)
        camb_results_minus = calculate_cmb(tau = param_value-param_step_size, lmax=lmax)
    elif param is 'thetastar':
        camb_results_plus = calculate_cmb(thetastar = param_value+param_step_size, lmax=lmax)
        camb_results_minus = calculate_cmb(thetastar = param_value-param_step_size, lmax=lmax)

    power_spectra_plus = get_power_spectra(camb_results_plus, which_spectra, ell)
    power_spectra_minus = get_power_spectra(camb_results_minus, which_spectra, ell)

    for TP in power_spectra_plus.keys():
        cl_camb_deriv[TP] = (power_spectra_plus[TP] - power_spectra_minus[TP]) / (2. * param_step_size)
    
    return cl_camb_deriv

def get_cmb_derivs_dic(ell, which_spectra, param_step_size_dic):
    
    #lcdm_params = ['As', 'mnu', 'neff', 'ns', 'ombh2', 'omch2', 'tau', 'thetastar']
    
    #FIXME: Fiducial values should be an input? at least for Neff?
    #fiducial values
    param_value_dic = {'As':2.196e-9,'mnu':0.06,'Neff':3.046,'ns':0.9655,'ombh2':0.0222,'omch2':0.1197,'tau':0.06,'thetastar':0.010409}
    
    cl_camb_deriv_dic = {}
    for param in param_value_dic.keys():
        cl_camb_deriv_dic[param] = get_cmb_derivs(ell, which_spectra, param, param_value_dic[param], param_step_size_dic[param])
    
    return cl_camb_deriv_dic

def recalculate_cl_residual(freq_bands, ell, nl_dic, bl_dic):
    
    import misc, ilc
    
    """
    Recalculates the foreground cleaned noise levels after
    adding various systematic spectra
    
    Arguments:
    -------------------------------------------------------
    freq_bands[]: array or list of band centers in GHz
    ell[]: array of angular scale values to generate 
            beams over
    nl_dic[]: dictionary of noise spectra with added systematics
            organized by T and P keys and frequency auto and 
            cross spectra
    bl_dic[]: dictionary of beam arrays organized by frequency
            bands
            
    Output:
    -------------------------------------------------------
    cl_residual[]: dictionary of foreground cleaned noise
            spectra organized by TT, EE, and TE keys
    """
    
    #take in nl, cl, param_dict from ilc, ell, freq_band
    
    paramfile = 'param.ini'
    param_dict = misc.fn_get_param_dict(paramfile)
       
    el, cl_dic = ilc.get_analytic_covariance(param_dict, freq_bands, nl_dic=nl_dic, \
                                             bl_dic=bl_dic)
    
    cl_residual = misc.residual_power(param_dict, freq_bands, ell, cl_dic)    
    
    
    return cl_residual
