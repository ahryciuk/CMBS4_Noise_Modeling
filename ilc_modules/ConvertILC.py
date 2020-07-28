import sys, os
import pandas as pd

#import ilc modules


def calc_cls(nl_dic, bl_dic, ell, band_centers, cwd):
    
    os.chdir(cwd + '/ilc_modules')
    import ilc, flatsky, misc, foregrounds as fg
    os.chdir(cwd)
    
    ############
    which_spec_arr = ['TT','EE']
    final_comp = 'cmb'
    freqcalib_fac = None
    ############
    
    #params
    paramfile = 'params.ini'
    # read and store param dict
    os.chdir(cwd + '/ilc_modules')
    param_dict = misc.fn_get_param_dict(paramfile)
    
    #get the CMB, noise, and foreground covariance
    try:
        ignore_fg = param_dict['ignore_fg']
    except:
        ignore_fg = []

    ignore_fg.append(final_comp.lower()) #the required component need not go into the covariance matrix.
    print('Ignoring ' + str(ignore_fg))

    #freqarr = [145]
    #param_dict['which_gal_mask'] = 0

    #Generate Cl spectrum
    cl_dic = {}
    include_gal = param_dict['include_gal']
    for which_spec in which_spec_arr:
        if which_spec == 'TT':
            el, cl_dic[which_spec] = ilc.get_analytic_covariance(param_dict, band_centers, \
                    nl_dic = nl_dic['T'], ignore_fg = ignore_fg, include_gal = include_gal, bl_dic = bl_dic)
        else:
            el, cl_dic[which_spec] = ilc.get_analytic_covariance\
                        (param_dict, band_centers, nl_dic = nl_dic['P'], ignore_fg = ignore_fg, which_spec = which_spec, \
                        pol_frac_per_cent_dust = param_dict['pol_frac_per_cent_dust'], \
                        pol_frac_per_cent_radio = param_dict['pol_frac_per_cent_radio'], \
                        pol_frac_per_cent_tsz = param_dict['pol_frac_per_cent_tsz'], \
                        pol_frac_per_cent_ksz = param_dict['pol_frac_per_cent_ksz'], \
                        include_gal = include_gal, bl_dic = bl_dic)


    #ilc weights
    lmin = ell[0]
    weights_dic = {}
    for which_spec in which_spec_arr:
        weights_dic[which_spec] = ilc.get_multipole_weightsarr(final_comp, band_centers, el, cl_dic[which_spec], lmin, 
                                                               freqcalib_fac)#, ignore_fg)


    #calculate residuals and separate
    #get the residual power now
    cl_residual = {}
    for which_spec in which_spec_arr:
        cl_residual[which_spec] = ilc.residual_power(param_dict, band_centers, el, cl_dic[which_spec], 
                                                                              final_comp = final_comp, 
                                                                              freqcalib_fac = freqcalib_fac, 
                                                                              return_weights = 0)

    return cl_dic, weights_dic, cl_residual

def convert_params_excel(exp_dir, cwd):
    
    params_file = 'params.ini'
    
    os.chdir(cwd + '/ilc_modules')
    params = pd.read_csv('params.ini',comment='#',sep='=',names=['Parameter','Value'])
    
    os.chdir(exp_dir)
    params.to_excel('comp_sep_params.xlsx')    
    
    return 



def append_params(exp_dir, cwd):
    
    os.chdir(exp_dir)
    append_params = pd.read_excel('comp_sep_params.xlsx')
    
    os.chdir(cwd + '/ilc_modules')
    append_params.to_csv('params.ini',sep='=',columns=['Parameter','Value'],index=False)
   
    print('Done!')
    
    return





