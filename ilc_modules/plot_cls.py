def plot_cls(param_dict, freqarr, el, cwd, cl_residual, weights_dic, nl_dic):
    import matplotlib.pyplot as plt
    from pylab import subplots_adjust,subplot2grid,plot
    import os
    import numpy as np

    
    #params
    #paramfile = 'params.ini'
    # read and store param dict
    #os.chdir(cwd + '/ilc_modules')
    #import misc, flatsky, ilc, foregrounds as fg
    #param_dict = misc.fn_get_param_dict(paramfile)

    ##Constants##
    include_gal = param_dict['include_gal']
    Tcmb = param_dict['T_cmb']
    which_spec_arr = ['TT','EE']
    remove_atm = param_dict['remove_atm']
    corr_noise = 1
    #############

    freq0, lmax = param_dict['freq0'], param_dict['lmax']
    if include_gal:
        foregrounds_to_plot = ['kSZ', 'tSZ', 'DG-Po', 'DG-Cl', 'RG', 'dust', 'sync']
        pol_foregrounds_to_plot = ['dust', 'sync']
    else:
        foregrounds_to_plot = ['kSZ', 'tSZ', 'DG-Po', 'DG-Cl', 'RG']
        pol_foregrounds_to_plot = []

    #CAMB output for plotting
    #camb_file = param_dict['Dlfile_len']
    #camb_file = '%s/%s' %(param_dict['data_folder'], param_dict['Dlfile_len'])
    camb_file = '/home/ahryciuk/CMBS4_Noise_Modeling/ilc_modules/output_planck_r_0.0_2015_cosmo_lensedCls.dat'
    el_camb = np.loadtxt(camb_file, usecols = [0])
    dl_camb = np.loadtxt(camb_file, usecols = [1,2,3,4])

    cl_camb = ( Tcmb**2. * dl_camb * 2 * np.pi ) / ( el_camb[:,None] * (el_camb[:,None] + 1) )
    cl_camb *= 1e12
    cl_TT, cl_EE, cl_BB, cl_TE = cl_camb.T


    #clf(); 
    fsval = 8
    tr, tc = 6, len(which_spec_arr)
    subplots_adjust(wspace=0.1, hspace = 0.1)
    xmin, xmax = 20, 7000
    #first plot weights
    colordic = {27:'indigo', 39:'blue', 93: 'royalblue', 145: 'darkgreen', 225: 'goldenrod', 280: 'darkred'}
    rspan, cspan = 2, 1
    curr_row = 0
    for cntr, which_spec in enumerate( which_spec_arr ):
        ax = subplot2grid((tr,tc), (curr_row, cntr), rowspan = rspan, colspan = cspan, xscale = 'log')#, yscale = 'log')
        for frqcntr, freq in enumerate( freqarr ):
            plot(weights_dic[which_spec][frqcntr], color = colordic[freq], label = r'%s' %(freq), lw = 0.5)
        plot(np.sum(weights_dic[which_spec], axis = 0), 'k--', label = r'Sum', lw = 0.5)
        plt.axhline(lw=0.3);
        #xlabel(r'Multipole $\ell$');
        plt.setp(ax.get_xticklabels(which = 'both'), visible=False)
        if cntr == 0:
            plt.ylabel(r'Weight $W_{\ell}$')
            plt.legend(loc = 3, fontsize = 5, ncol = 4, handlelength = 2., handletextpad = 0.1)
        else:
            plt.setp(ax.get_yticklabels(which = 'both'), visible=False)
        plt.ylim(-3., 3.);
        plt.xlim(xmin, xmax);
        for label in ax.get_xticklabels(): label.set_fontsize(fsval)
        for label in ax.get_yticklabels(): label.set_fontsize(fsval)        

        plt.title(r'%s' %(which_spec))#, fontsize = 10)

    curr_row = rspan
    rspan = tr - rspan
    for cntr, which_spec in enumerate( which_spec_arr ):
        #ax = subplot(1,2,cntr+1, xscale = 'log', yscale = 'log')
        ax = subplot2grid((tr,tc), (curr_row, cntr), rowspan = rspan, colspan = cspan, xscale = 'log', yscale = 'log')
        plot(el, cl_residual[which_spec], 'black', lw = 2., label = r'Residual')
        if which_spec == 'TT':
            plot(el_camb, cl_TT, 'gray', lw = 1., label = r'TT')
            '''
            cl_fg = np.zeros(len(el))
            for curr_fg in foregrounds_to_plot:
                if curr_fg == 'dust':
                    el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'dust', 145, 145, 'TT', bl_dic = bl_dic, el = el)
                elif curr_fg == 'sync':
                    el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'sync', 145, 145, 'TT', bl_dic = bl_dic, el = el)
                else:
                    el_, cl_curr_fg = fg.get_foreground_power_spt(curr_fg, freq1 = freq0, lmax = lmax)
                #plot(el, cl_curr_fg, lw = 0.5, ls = '--', label = r'150: %s' %(curr_fg), alpha = 0.4)
                cl_fg += cl_curr_fg
            plot(el, cl_fg, lw = 5., ls = '--', label = r'150: All foregrounds', alpha = 1.)
            '''
        elif which_spec == 'EE':
            plot(el_camb, cl_EE, 'gray', lw = 0.5)#, label = r'EE')
            plot(el_camb, cl_TE, 'gray', ls = '-', lw = 0.5)#, label = r'TE')        
            plot(el_camb, abs( cl_TE ), 'navy', ls = '--', lw = 0.5) 
            '''
            for curr_fg in pol_foregrounds_to_plot:
                if curr_fg == 'dust':
                    el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'dust', 145, 145, 'EE', bl_dic = bl_dic, el = el)
                elif curr_fg == 'sync':
                    el, cl_curr_fg = fg.get_cl_galactic(param_dict, 'sync', 145, 145, 'EE', bl_dic = bl_dic, el = el)
                plot(el, cl_curr_fg, lw = 0.5, ls = '--', label = r'150: %s' %(curr_fg), alpha = 0.4)
            '''
        #for freq in freqarr:
        #    plot(el, cl_dic[which_spec][(freq,freq)], color = colordic[freq], lw = 0.5, ls = '-', label = r'%s' %(freq), alpha = 1.)        
        for freq in freqarr:
            if which_spec == 'TT':
                nl = nl_dic['T'][(freq, freq)]
            elif which_spec == 'EE':
                nl = nl_dic['P'][(freq, freq)]
            elif which_spec == 'TE':
                nl = nl_dic['T'][(freq, freq)] * 0.
            plot(el, nl, color = colordic[freq], lw = 0.5, ls = '--', label = r'Noise: %s' %(freq))#, alpha = 0.5)
        #legend(loc=3, fancybox=1, ncol = 4, fontsize = 6);

        plt.xlim(xmin, xmax);
        plt.ylim(1e-8,1e6);
        plt.xlabel(r'Multipole $\ell$')
        if cntr == 0: 
            plt.ylabel(r'$C_{\ell}$')
            plt.legend(loc = 1, fontsize = 6, ncol = 2, handlelength = 2., handletextpad = 0.1)
        else:
            plt.setp(ax.get_yticklabels(which = 'both'), visible=False)
        for label in ax.get_xticklabels(): label.set_fontsize(fsval)
        for label in ax.get_yticklabels(): label.set_fontsize(fsval)

    #tit = 'Galaxy = %s; Mask = %s; Bands = %s' %(include_gal, param_dict['which_gal_mask'], str(freqarr))
    if remove_atm:
        tit = 'Galaxy = %s; Mask = %s; Bands = %s; no 1/f' %(include_gal, param_dict['which_gal_mask'], str(freqarr))
    else:
        if include_gal:
            tit = 'Galaxy = %s; Mask = %s; Bands = %s' %(include_gal, param_dict['which_gal_mask'], str(freqarr))    
        else:
            tit = 'Galaxy = %s; Mask = N/A; Bands = %s' %(include_gal, str(freqarr))    
    if not corr_noise:
        tit = '%s; No corr. noise' %(tit)
    plt.suptitle(r'%s' %tit, x = 0.53, y = 1.)
    #savefig(plname)
    plt.show()