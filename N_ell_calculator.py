#from __future__ import print_function
import numpy as np
import os as os


####################################################################
####################################################################
### LAT & SAT CALCULATOR ###
####################################################################
####################################################################

#Inputs: array of band centers, array of beams, sensitivities from calcBolos sim, f_knees, C's, and alphas, 
#survey time, f_sky, retention after observing cuts = 0.2, non-uniformity parameter = 0.85, ell_pivot

#number of tubes per camera?


def S4_noise(band_centers=[27.,39.,93.,145.,225.,280.], beam_sizes=[7.4,5.1,2.2,1.4,1.0,0.9], Sens=[48.,24.,5.4,6.7,15.,36.], ell_max=1e4, f_knees=[700.,700.,700.,700.,700.,700.], Cs=[200,7.7,1800,12000,68000,124000], alpha_temp=[-3.5,-3.5,-3.5,-3.5,-3.5,-3.5], survey_time=5., f_sky=0.4, ret_after_obs_cuts=0.2, non_uniformity_param=0.85, ell_pivot=[1000.,1000.,1000.,1000.,1000.,1000.],delta_ell=1,alpha_pol=[-0.4,-0.4,-0.4,-0.4,-0.4,-0.4],NTubes=[1,1,1,1,1,1],model_num=1):
    
    print("band centers: ", band_centers, "[GHz]")
    print("beam sizes: "  , beam_sizes, "[arcmin]")
   
    #NTubes = []
    #sect = len(band_centers)/3
    #for i in range(len(band_centers)):
    #    if i < len(band_centers)/3:
    #        NTubes.append(NTubes_LF)
    #    elif i >= len(band_centers)/3 and i < 2*len(band_centers)/3:
    #        NTubes.append(NTubes_MF)
    #    else:
    #        NTubes.append(NTubes_UHF)
    
    # sensitivity in uK*sqrt(s)
    # set noise to irrelevantly high value when NTubes=0

    #input error checking
    assert(f_sky > 0 and f_sky <= 1.)
    assert(model_num in [0,1])
    assert(len(band_centers) == len(beam_sizes))
    assert(len(beam_sizes) == len(Sens))
    assert(len(Sens) == len(f_knees))
    assert(len(f_knees) == len(Cs))
    assert(len(Cs) == len(alpha_temp))
    assert(len(alpha_temp) == len(ell_pivot))
    assert(len(ell_pivot) == len(alpha_pol))
    assert(len(alpha_pol) == len(NTubes))
    
    ####################################################################
    ## calculate the survey area and time
    #survey_time = 5. #years
    t = survey_time * 365.25 * 24. * 3600.    ## convert years to seconds
    t = t * ret_after_obs_cuts   ## retention after observing efficiency and cuts
    t = t * non_uniformity_param  ## a kludge for the noise non-uniformity of the map edges
    A_SR = 4. * np.pi * f_sky  ## sky areas in Steradians
    A_deg =  A_SR * (180/np.pi)**2  ## sky area in square degrees
    A_arcmin = A_deg * 3600.
    print("sky area: ", A_deg, "degrees^2")

    ####################################################################
    ## make the ell array for the output noise curves
    ell = np.arange(2,ell_max,delta_ell)

    ####################################################################
    ###   CALCULATE N(ell) for Temperature
    ## calculate the experimental weight
    Weights=[]
    Map_white_noise_levels=[]
    AN_T=[]
    for i in range(len(Sens)):
        Weights.append(Sens[i]/np.sqrt(t))
        Map_white_noise_levels.append(Sens[i] * np.sqrt(A_arcmin) / np.sqrt(t))
        #AN_T.append(Cs[i] * (ell/ell_pivot[i])**alpha_temp[i] * A_SR / t / (2.*NTubes[i]) )
        
        #Atacama model
        if not model_num:
            AN_T.append(Cs[i] * (ell/ell_pivot[i])**alpha_temp[i] * A_SR / t / (2.*NTubes[i]) )
        
        #Changed normalization to line up with comp. sep. code---Ask Jeff maybe reflects SO not SP
        #SP model
        else:
            AN_T.append(Sens[i]**2. * (A_SR) /(t) * (ell/ell_pivot[i])**alpha_temp[i] / (2.*NTubes[i]) )


#    Map_white_noise_levels= np.array([MN_T_27,MN_T_39,MN_T_93,MN_T_145,MN_T_225,MN_T_280])
    print("white noise level: ",Map_white_noise_levels ,"[uK-arcmin]")

    ## calculate the atmospheric contribution for T (based on Matthew's model)
    # the 2*NTube factor comes from Matthew H.'s email on 1-25-18
    # handle cases where there are zero tubes of some kind
    
   

    ## calculate N(ell)
    ## include the impact of the beam
    LA_beams = beam_sizes / np.sqrt(8. * np.log(2)) /60. * np.pi/180.
                                ## lac beams as a sigma expressed in radians
        
    # include cross-frequency correlations in the atmosphere
    # Matthew H.: the most well-motivated model for atmospheric correlation between bands is as follows:
    #   - the atmospheric noise is 100% correlated between bands in a single optics tube.
    #   - the atmospheric noise is not correlated, at least for ell>400, between adjacent optics tubes.
    # use correlation coefficient of r=0.9 within each dichroic pair and 0 otherwise
    r_atm = 0.9
    
    #if in same optics tube then gets r_atm correlation, otherwise make array of zeros#
    
#    for i in range(0,len(AN_T)-1):
#        AN_T_cross.append( (AN_T[i]*AN_T[i+1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[i]**2. + LA_beams[i+1]**2.)) )
#        if i == len(AN_T)-2:
#            AN_T_cross.append( (AN_T[0]*AN_T[-1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[0]**2. + LA_beams[-1]**2.)) )
            
    #AN_T_cross.append( r_atm * np.sqrt(AN_T[0]*AN_T[1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[0]**2. + LA_beams[1]**2.)) )
    #AN_T_cross.append( r_atm * np.sqrt(AN_T[2]*AN_T[3])*np.exp( (ell*(ell+1)/2.)*(LA_beams[2]**2. + LA_beams[3]**2.)) )
    #AN_T_cross.append( r_atm * np.sqrt(AN_T[4]*AN_T[5])*np.exp( (ell*(ell+1)/2.)*(LA_beams[4]**2. + LA_beams[5]**2.)) )
   
    N_ell_T_LA = {}
    
    #allow margin for small perturbation in central frequencies
    #corr_freq = [27., 93., 225.]
    LF1,LF2,LF3,MF1,MF2,UHF1,UHF2 = 0,0,0,0,0,0,0
    for freq in band_centers:
        if freq < 24 and freq > 18:
            LF1 = freq
        elif freq < 30 and freq > 25:
            LF2 = freq
        elif freq < 42 and freq > 35:
            LF3 = freq
        elif freq < 96 and freq > 84:
            MF1 = freq
        elif freq < 156 and freq > 144:
            MF2 = freq
        elif freq < 226 and freq > 214:
            UHF1 = freq
        elif freq < 286 and freq > 274:
            UHF2 = freq
            
    corr_freq = {LF1:LF1,LF2:LF3,LF3:LF2,MF1:MF2,MF2:MF1,UHF1:UHF2,UHF2:UHF1}
            
    for freq1 in band_centers:
        for freq2 in band_centers:
            i = band_centers.index(freq1)
            j = band_centers.index(freq2)
            if freq1 == freq2:
                N_ell_T_LA[(freq1,freq2)] = (Weights[i]**2.*A_SR + AN_T[i]) * np.exp( ell*(ell+1) * LA_beams[i]**2. )
      
            else:
                if freq2 is corr_freq[freq1]:
                    N_ell_T_LA[(freq1,freq2)] = r_atm * np.sqrt(AN_T[i]*AN_T[j]) * np.exp( (ell*(ell+1))*(LA_beams[i]**2. + LA_beams[j]**2.) )
                else:
                    N_ell_T_LA[(freq1,freq2)] = np.zeros(len(ell))
                    
#            elif j == i+1:
#                if freq1 in corr_freq:
#                    N_ell_T_LA[(freq1,freq2)] = r_atm * np.sqrt(AN_T[i]*AN_T[j]) * np.exp( (ell*(ell+1))*(LA_beams[i]**2. + #LA_beams[j]**2.) )
#                else:
#                    N_ell_T_LA[(freq1,freq2)] = np.zeros(len(ell))
#            else:
#                N_ell_T_LA[(freq1,freq2)] = np.zeros(len(ell))
                

    
    #for i in range(len(N_ell_T)):               
    #    N_ell_T[i] *= np.exp( ell*(ell+1) * LA_beams[i]**2.)
        

    ## make an array of noise curves for T
    # include cross-correlations due to atmospheric noise
    #N_ell_T_LA = np.array([N_ell_T, N_ell_T_cross])

    ####################################################################
    ###   CALCULATE N(ell) for Polarization
     ## calculate the astmospheric contribution for P
    
    try:
        AN_P=[]
        N_ell_P=[]
        N_ell_P_atm=[]
        #white_noise_level_P=[ for i in range(len(Weights)), (Weights[i] * np.sqrt(2))**2. * A_SR for j in range(len(ell)) ]

        white_noise_level_P = np.full( (len(f_knees) , len(ell)), 0.)

        #for i in range(len(f_knees)):
        #    for j in range(len(ell)):
        #        white_noise_level_P[i][j] = (Weights[i] * np.sqrt(2))**2. * A_SR
                #Map_white_noise_levels.append(white_noise_level_P[i][j])
                #white_noise_level_P[i][j] = white_noise[i]

       
        for i in range(len(f_knees)):
            white_noise_level_P[i] = (Weights[i] * np.sqrt(2))**2. * A_SR
    
        #Atacama Model (Atmosphere limited regime):
        if model_num == 0:
            for i in range(len(f_knees)):
                #AN_P.append( (ell / f_knees[i] )**alpha_pol[i] + 1.)
                AN_P.append( (ell / f_knees[i] )**alpha_pol[i] )

            for i in range(len(f_knees)):
                #N_ell_P.append( (Weights[i] * np.sqrt(2))**2. * A_SR * AN_P[i] )
                N_ell_P.append( white_noise_level_P[i] * (AN_P[i]+1) )
                #N_ell_P_atm.append( (Weights[i] * np.sqrt(2))**2. * A_SR * AN_P[i] )
                N_ell_P_atm.append( white_noise_level_P[i] * AN_P[i] )

            for i in range(len(N_ell_P)):
                N_ell_P[i] *= np.exp( ell*(ell+1) * LA_beams[i]**2 )




        #South Pole Atmosphere Model (Receiver limited regime):
        elif model_num == 1:
            for i in range(len(f_knees)):
                #N_ell_P_atm.append( (Weights[i] * np.sqrt(2))**2. * A_SR * (ell/f_knees[i])**alpha_pol )
                #N_ell_P_atm.append( (7e-9 * np.sqrt(2))**2. * A_SR * (ell/f_knees[i])**alpha_pol )
                AN_P.append( (ell/f_knees[i])**alpha_pol[i] )
                N_ell_P_atm.append(white_noise_level_P[i] * (ell/f_knees[i])**alpha_pol[i] )

            for i in range(len(f_knees)):
                N_ell_P.append( white_noise_level_P[i]*np.sqrt(AN_P[i]**2. + 1.) )
                #N_ell_P.append(N_ell_P_atm[i]**2. + white_noise_level_P[i]**2. )
                
            for i in range(len(N_ell_P)):
                N_ell_P[i] *= np.exp( ell*(ell+1) * LA_beams[i]**2. )
            
    except IndexError:
        print('The lengths of the list parameters are not consistent')
            
    N_ell_P_LA = {}
    for freq1 in band_centers:
        for freq2 in band_centers:
            i = band_centers.index(freq1)
            j = band_centers.index(freq2)
            if freq1 == freq2:
                N_ell_P_LA[(freq1,freq2)] = N_ell_P[i]
            else:
                if freq2 is corr_freq[freq1]:
                    N_ell_P_LA[(freq1,freq2)] = r_atm * np.sqrt( N_ell_P[i] * N_ell_P[j] )
                else:
                    N_ell_P_LA[(freq1,freq2)] = np.zeros(len(ell))
#            elif j == i+1:
#                if freq1 in corr_freq:
                    #N_ell_P_LA[(freq1,freq2)] = r_atm * np.sqrt(N_ell_P_atm[i] * N_ell_P_atm[j])*np.exp( (ell*(ell+1)/2.) * (LA_beams[i]**2. + LA_beams[j]**2.) )
                    #N_ell_P_LA[(freq1,freq2)] = r_atm * np.sqrt((Weights[i] * np.sqrt(2))**2. * A_SR *AN_P[i] * (Weights[i] * np.sqrt(2))**2. * A_SR * AN_P[j])*np.exp( (ell*(ell+1)) * (LA_beams[i]**2. + LA_beams[j]**2.) )
#                    N_ell_P_LA[(freq1,freq2)] = r_atm * np.sqrt( N_ell_P[i] * N_ell_P[j] )
#                else:
#                    N_ell_P_LA[(freq1,freq2)] = np.zeros(len(ell))
#            else:
#                N_ell_P_LA[(freq1,freq2)] = np.zeros(len(ell))
         
  
    return(ell, N_ell_T_LA, N_ell_P_LA, Map_white_noise_levels,corr_freq)



def package_nl_dic(band_centers, ell, 
                   N_ell_T_LA, N_ell_P_LA):
    
    #Spectra Keys
    TParr = ['T','P']
    
    #Power Spectra Keys
    which_spec_arr = ['TT','EE']
    
    #Frequency Bands
    specs_dic = {}
    
    #initialize N_ell dictionary
    nl_dic = {}
    for TP in TParr:
        nl_dic[TP] = {}
        for freq1 in band_centers:
            for freq2 in band_centers:
                if freq1 == freq2:
                    if TP == 'T':
                        nl_dic[TP][(freq1,freq2)] = N_ell_T_LA[0][band_centers.index(freq1)]
                    elif TP == 'P':
                        nl_dic[TP][(freq1,freq2)] = N_ell_P_LA[0][band_centers.index(freq1)]
                else:
                    if TP == 'T':
                        nl_dic[TP][(freq1,freq2)] = N_ell_T_LA[1][band_centers.index(freq2)]
                    elif TP == 'P':
                        nl_dic[TP][(freq1,freq2)] = N_ell_P_LA[1][band_centers.index(freq2)]


                        
    return nl_dic

def package_beam_dic( white_noise_levels_T, white_noise_levels_P,
                     band_centers=[27.,39.,93.,145.,225.,278.], 
                     beam_sizes=[7.4,5.1,2.2,1.4,1.0,0.9], 
                     ellmax=1e4 ):
    import misc
    
    #need freqarr beam arr noise arr elkneearr alphakneearr
    #specs_dic = {}
    #for freq in band_centers:
    #    specs_dic.update({freq: [beam_sizes[band_centers.index(freq)], white_noise_levels_T[band_centers.index(freq)],
    #                             elknee_T[band_centers.index(freq)], alpha_temp[band_centers.index(freq)], 
    #                             white_noise_levels_P[band_centers.index(freq)], elknee_P[band_centers.index(freq)], 
    #                             alpha_pol[band_centers.index(freq)]]})
  
    
    #freqarr = sorted( specs_dic.keys() )
    #make beam noise dictionary
    #collect beam and noise into a dic; elknee and alpha into a dic
    
    beam_noise_dic = {}
    #elknee_dic = {}
    TParr = ['T','P']
    for TP in TParr:
        beam_noise_dic[TP] = {}
        #elknee_dic[TP] = {} 
        if TP == 'T':
            for freq in band_centers:
                beam_noise_dic[TP][freq] = [beam_sizes[band_centers.index(freq)],white_noise_levels_T[band_centers.index(freq)]]
            #freqarr, beamarr, noisearr, elkneearr, alphakneearr = freqarr, beamarr, noisearr_T, elkneearr_T, alphakneearr_T
        elif TP == 'P':
            for freq in band_centers:
                beam_noise_dic[TP][freq] = [beam_sizes[band_centers.index(freq)], white_noise_levels_P[band_centers.index(freq)]]
                
    
    #initialize beam dictionary
    bl_dic = misc.get_beam_dic(band_centers, beam_noise_dic['T'], ellmax)
       
    return bl_dic


def gen_cl_dic(param_dict, band_centers, nl_dic = None, 
               bl_dic = None, ignore_fg = [], which_spec = 'TT', 
               pol_frac_per_cent_dust = 0.02, pol_frac_per_cent_radio = 0.03, 
               pol_frac_per_cent_tsz = 0., pol_frac_per_cent_ksz = 0., include_gal = 0):
    
    #things from param_dict
    #freq0
    #fg_model
    #spec_index_dg_po,spec_index_dg_clus
    #Tcib
    
    #for get_cl_galactic
    #cl_gal_dic_dust_fname, cl_gal_dic_sync_fname
    #cl_gal_folder
    #cl_gal_dic_sync_fname_forced
    
    
    return cl


