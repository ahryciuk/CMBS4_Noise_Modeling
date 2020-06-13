from __future__ import print_function
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


def S4_noise(band_centers=[27.,39.,93.,145.,225.,280.], beam_sizes=[7.4,5.1,2.2,1.4,1.0,0.9], Sens=[48.,24.,5.4,6.7,15.,36.], ell_max=1e4, f_knees=[700.,700.,700.,700.,700.,700.], Cs=[200,7.7,1800,12000,68000,124000], alpha_temp=-3.5, survey_time=5., f_sky=0.4, ret_after_obs_cuts=0.2, non_uniformity_param=0.85, ell_pivot=1000.,delta_ell=5,alpha_pol=0.4,NTubes_LF=1, NTubes_MF=4,NTubes_UHF=2,model_num=0):
    
    print("band centers: ", band_centers, "[GHz]")
    print("beam sizes: "  , beam_sizes, "[arcmin]")
   
    NTubes = []
    sect = len(band_centers)/3
    for i in range(len(band_centers)):
        if i < len(band_centers)/3:
            NTubes.append(NTubes_LF)
        elif i >= len(band_centers)/3 and i < 2*len(band_centers)/3:
            NTubes.append(NTubes_MF)
        else:
            NTubes.append(NTubes_UHF)
    
    # sensitivity in uK*sqrt(s)
    # set noise to irrelevantly high value when NTubes=0

    assert(f_sky > 0 and f_sky <= 1.)
    
    ####################################################################
    ## calculate the survey area and time
    #survey_time = 5. #years
    t = survey_time * 365.25 * 24. * 3600.    ## convert years to seconds
    t = t * 0.2   ## retention after observing efficiency and cuts
    t = t * 0.85  ## a kludge for the noise non-uniformity of the map edges
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
#    W_T_27  = S_LA_27[sensitivity_mode]  / np.sqrt(t)
#    W_T_39  = S_LA_39[sensitivity_mode]  / np.sqrt(t)
#    W_T_93  = S_LA_93[sensitivity_mode]  / np.sqrt(t)
#    W_T_145 = S_LA_145[sensitivity_mode] / np.sqrt(t)
#    W_T_225 = S_LA_225[sensitivity_mode] / np.sqrt(t)
#    W_T_280 = S_LA_280[sensitivity_mode] / np.sqrt(t)
    Weights=[]
    Map_white_noise_levels=[]
    AN_T=[]
    for i in range(len(Sens)):
        Weights.append(Sens[i]/np.sqrt(t))
        Map_white_noise_levels.append(Sens[i] * np.sqrt(A_arcmin) / np.sqrt(t))
        AN_T.append(Cs[i] * (ell/ell_pivot)**alpha_temp * A_SR / t / (2.*NTubes[i]))

    ## calculate the map noise level (white) for the survey in uK_arcmin for temperature
#    MN_T_27  = W_T_27  * np.sqrt(A_arcmin)
#    MN_T_39  = W_T_39  * np.sqrt(A_arcmin)
#    MN_T_93  = W_T_93  * np.sqrt(A_arcmin)
#    MN_T_145 = W_T_145 * np.sqrt(A_arcmin)
#    MN_T_225 = W_T_225 * np.sqrt(A_arcmin)
#    MN_T_280 = W_T_280 * np.sqrt(A_arcmin)
#    Map_white_noise_levels= np.array([MN_T_27,MN_T_39,MN_T_93,MN_T_145,MN_T_225,MN_T_280])
    print("white noise level: ",Map_white_noise_levels ,"[uK-arcmin]")

    ## calculate the atmospheric contribution for T (based on Matthew's model)
    # the 2*NTube factor comes from Matthew H.'s email on 1-25-18
    # handle cases where there are zero tubes of some kind
    
   
#    AN_T_27  = C_27  * (ell/ell_pivot)**alpha_temp * A_SR / t / (2.*NTubes_LF)
#    AN_T_39  = C_39  * (ell/ell_pivot)**alpha_temp * A_SR / t / (2.*NTubes_LF)

#    AN_T_93  = C_93  * (ell/ell_pivot)**alpha_temp * A_SR / t / (2.*NTubes_MF)
#    AN_T_145 = C_145 * (ell/ell_pivot)**alpha_temp * A_SR / t / (2.*NTubes_MF)

#    AN_T_225 = C_225 * (ell/ell_pivot)**alpha_temp * A_SR / t / (2.*NTubes_UHF)
#    AN_T_280 = C_280 * (ell/ell_pivot)**alpha_temp * A_SR / t / (2.*NTubes_UHF)
        
    ## calculate N(ell)
    N_ell_T=[]
    for i in range(len(Weights)):
        N_ell_T.append(Weights[i]**2. * A_SR + AN_T[i])                
                    
    ## include the impact of the beam
    LA_beams = beam_sizes / np.sqrt(8. * np.log(2)) /60. * np.pi/180.
                                ## lac beams as a sigma expressed in radians
        
    # include cross-frequency correlations in the atmosphere
    # Matthew H.: the most well-motivated model for atmospheric correlation between bands is as follows:
    #   - the atmospheric noise is 100% correlated between bands in a single optics tube.
    #   - the atmospheric noise is not correlated, at least for ell>400, between adjacent optics tubes.
    # use correlation coefficient of r=0.9 within each dichroic pair and 0 otherwise
    r_atm = 0.9
    AN_T_cross=[]               
    for i in range(0,len(AN_T)-1):
        AN_T_cross.append( (AN_T[i]*AN_T[i+1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[i]**2. + LA_beams[i+1]**2.)) )
        if i == len(AN_T)-2:
            AN_T_cross.append( (AN_T[0]*AN_T[-1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[0]**2. + LA_beams[-1]**2.)) )
        
    N_ell_T_cross = AN_T_cross
#    AN_T_27x39 = r_atm * np.sqrt(AN_T_27 * AN_T_39)
#    AN_T_93x145 = r_atm * np.sqrt(AN_T_93 * AN_T_145)
#    AN_T_225x280 = r_atm * np.sqrt(AN_T_225 * AN_T_280)

    
                               
                               
#    N_ell_T_27   = (W_T_27**2. * A_SR) + AN_T_27
#    N_ell_T_39   = (W_T_39**2. * A_SR) + AN_T_39
#    N_ell_T_93   = (W_T_93**2. * A_SR) + AN_T_93
#    N_ell_T_145  = (W_T_145**2. * A_SR) + AN_T_145
#    N_ell_T_225  = (W_T_225**2. * A_SR) + AN_T_225
#    N_ell_T_280  = (W_T_280**2. * A_SR) + AN_T_280
    # include cross-correlations due to atmospheric noise
#    N_ell_T_27x39 = AN_T_27x39
#    N_ell_T_93x145 = AN_T_93x145
#    N_ell_T_225x280 = AN_T_225x280

    
    for i in range(len(N_ell_T)):               
        N_ell_T[i] *= np.exp( ell*(ell+1)* LA_beams[i]**2.)
        
#    for i in range(len(N_ell_T_cross)):
#        N_ell_T_cross[i] *= np.exp( (ell*(ell+1)/2.) )
                            
#    N_ell_T_27  *= np.exp( ell*(ell+1)* LA_beams[0]**2. )
#    N_ell_T_39  *= np.exp( ell*(ell+1)* LA_beams[1]**2. )
#    N_ell_T_93  *= np.exp( ell*(ell+1)* LA_beams[2]**2. )
#    N_ell_T_145 *= np.exp( ell*(ell+1)* LA_beams[3]**2. )
#    N_ell_T_225 *= np.exp( ell*(ell+1)* LA_beams[4]**2. )
#    N_ell_T_280 *= np.exp( ell*(ell+1)* LA_beams[5]**2. )
#    N_ell_T_27x39 *= np.exp( (ell*(ell+1)/2.) * (LA_beams[0]**2. + LA_beams[1]**2.) )
#    N_ell_T_93x145 *= np.exp( (ell*(ell+1)/2.) * (LA_beams[2]**2. + LA_beams[3]**2.) )
#    N_ell_T_225x280 *= np.exp( (ell*(ell+1)/2.) * (LA_beams[4]**2. + LA_beams[5]**2.) )

    ## make an array of noise curves for T
    # include cross-correlations due to atmospheric noise
    N_ell_T_LA = np.array([N_ell_T, N_ell_T_cross])

    ####################################################################
    ###   CALCULATE N(ell) for Polarization
     ## calculate the astmospheric contribution for P
    
    try:
        AN_P=[]
        N_ell_P=[]
        N_ell_P_atm=[]
        N_ell_P_cross=[]
        #white_noise_level_P=[ for i in range(len(Weights)), (Weights[i] * np.sqrt(2))**2. * A_SR for j in range(len(ell)) ]

        white_noise_level_P = np.full( (len(f_knees) , len(ell)), 0.)

        #hard coded sanity check
        white_noise = [1e-5,1e-5,1.5e-5,1.5e-5,2.5e-4,3e-4]
        for i in range(len(f_knees)):
            for j in range(len(ell)):
                white_noise_level_P[i][j] = (Weights[i] * np.sqrt(2))**2. * A_SR
                #white_noise_level_P[i][j] = white_noise[i]

        #Atacama Model (Atmosphere limited regime):
        if model_num == 0:
            for i in range(len(f_knees)):
                AN_P.append( (ell / f_knees[i] )**alpha_pol + 1.)

            for i in range(len(f_knees)):
                N_ell_P.append( (Weights[i] * np.sqrt(2))**2. * A_SR * AN_P[i] )
                N_ell_P_atm.append( (Weights[i] * np.sqrt(2))**2. * A_SR * (ell / f_knees[i] )**alpha_pol )

            for i in range(0,len(N_ell_P)-1):
                N_ell_P_cross.append( r_atm * np.sqrt(N_ell_P_atm[i] * N_ell_P_atm[i+1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[i]**2. + LA_beams[i+1]**2.)) )
                if i == len(N_ell_P)-2:
                        N_ell_P_cross.append(r_atm * np.sqrt(N_ell_P_atm[0] * N_ell_P_atm[-1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[0]**2. + LA_beams[-1]**2.)) )

            for i in range(len(N_ell_P)):
                N_ell_P[i] *= np.exp( ell*(ell+1) * LA_beams[i]**2 )




        #South Pole Atmosphere Model (Receiver limited regime):
        elif model_num == 1:
            for i in range(len(f_knees)):
                #N_ell_P_atm.append( (Weights[i] * np.sqrt(2))**2. * A_SR * (ell/f_knees[i])**alpha_pol )
                #N_ell_P_atm.append( (7e-9 * np.sqrt(2))**2. * A_SR * (ell/f_knees[i])**alpha_pol )
                N_ell_P_atm.append( (ell/f_knees[i])**alpha_pol )

            for i in range(len(f_knees)):
                N_ell_P.append( white_noise_level_P[i]*np.sqrt(N_ell_P_atm[i]**2. + 1.) )
                #N_ell_P.append(N_ell_P_atm[i]**2. + white_noise_level_P[i]**2. )


            for i in range(0,len(N_ell_P)-1):
                N_ell_P_cross.append( r_atm * np.sqrt(N_ell_P_atm[i] * N_ell_P_atm[i+1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[i]**2. + LA_beams[i+1]**2.)) )
                if i == len(N_ell_P)-2:
                    N_ell_P_cross.append( r_atm * np.sqrt(N_ell_P_atm[0] * N_ell_P_atm[-1])*np.exp( (ell*(ell+1)/2.)*(LA_beams[0]**2. + LA_beams[-1]**2.)) ) 

            for i in range(len(N_ell_P)):
                N_ell_P[i] *= np.exp( ell*(ell+1) * LA_beams[i]**2. )
            
    except IndexError:
        print('The lengths of the list parameters are not consistent')
        
        
        
    
#    AN_P_27  = (ell / f_knee_pol_LA_27 )**alpha_pol + 1. 
#    AN_P_39  = (ell / f_knee_pol_LA_39 )**alpha_pol + 1.
#    AN_P_93  = (ell / f_knee_pol_LA_93 )**alpha_pol + 1.  
#    AN_P_145 = (ell / f_knee_pol_LA_145)**alpha_pol + 1.  
#    AN_P_225 = (ell / f_knee_pol_LA_225)**alpha_pol + 1.  
#    AN_P_280 = (ell / f_knee_pol_LA_280)**alpha_pol + 1.

    ## calculate N(ell)
#    N_ell_P_27   = (W_T_27  * np.sqrt(2))**2. * A_SR * AN_P_27
#    N_ell_P_39   = (W_T_39  * np.sqrt(2))**2. * A_SR * AN_P_39
#    N_ell_P_93   = (W_T_93  * np.sqrt(2))**2. * A_SR * AN_P_93
#    N_ell_P_145  = (W_T_145 * np.sqrt(2))**2. * A_SR * AN_P_145
#    N_ell_P_225  = (W_T_225 * np.sqrt(2))**2. * A_SR * AN_P_225
#    N_ell_P_280  = (W_T_280 * np.sqrt(2))**2. * A_SR * AN_P_280
                    
    # include cross-correlations due to atmospheric noise
    # different approach than for T -- need to subtract off the white noise part to get the purely atmospheric part
#    N_ell_P_27_atm = (W_T_27  * np.sqrt(2))**2. * A_SR * (ell / f_knee_pol_LA_27 )**alpha_pol
#    N_ell_P_39_atm = (W_T_39  * np.sqrt(2))**2. * A_SR * (ell / f_knee_pol_LA_39 )**alpha_pol
#    N_ell_P_93_atm = (W_T_93  * np.sqrt(2))**2. * A_SR * (ell / f_knee_pol_LA_93 )**alpha_pol
#    N_ell_P_145_atm = (W_T_145  * np.sqrt(2))**2. * A_SR * (ell / f_knee_pol_LA_145 )**alpha_pol
#    N_ell_P_225_atm = (W_T_225  * np.sqrt(2))**2. * A_SR * (ell / f_knee_pol_LA_225 )**alpha_pol
#    N_ell_P_280_atm = (W_T_280  * np.sqrt(2))**2. * A_SR * (ell / f_knee_pol_LA_280 )**alpha_pol
#    N_ell_P_27x39 = r_atm * np.sqrt(N_ell_P_27_atm * N_ell_P_39_atm)
#    N_ell_P_93x145 = r_atm * np.sqrt(N_ell_P_93_atm * N_ell_P_145_atm)
#    N_ell_P_225x280 = r_atm * np.sqrt(N_ell_P_225_atm * N_ell_P_280_atm)
    
  
                    
#    for i in range(len(N_ell_P_cross)):
#        N_ell_P_cross[i] *= np.exp( (ell*(ell+1)/2.) )
                    
    ## include the imapct of the beam
#    N_ell_P_27  *= np.exp( ell*(ell+1)* LA_beams[0]**2 )
#    N_ell_P_39  *= np.exp( ell*(ell+1)* LA_beams[1]**2 )
#    N_ell_P_93  *= np.exp( ell*(ell+1)* LA_beams[2]**2 )
#    N_ell_P_145 *= np.exp( ell*(ell+1)* LA_beams[3]**2 )
#    N_ell_P_225 *= np.exp( ell*(ell+1)* LA_beams[4]**2 )
#    N_ell_P_280 *= np.exp( ell*(ell+1)* LA_beams[5]**2 )
#    N_ell_P_27x39 *= np.exp( (ell*(ell+1)/2.) * (LA_beams[0]**2. + LA_beams[1]**2.) )
#    N_ell_P_93x145 *= np.exp( (ell*(ell+1)/2.) * (LA_beams[2]**2. + LA_beams[3]**2.) )
#    N_ell_P_225x280 *= np.exp( (ell*(ell+1)/2.) * (LA_beams[4]**2. + LA_beams[5]**2.) )

    ## make an array of noise curves for P
    N_ell_P_LA = np.array([N_ell_P, N_ell_P_cross])
         
        
        
        
        
        
   
  
    '''                    
    #Plotting:
    colors = ['b','r','g','m','k','y']
    fig, (plt1, plt2) = plt.subplots(2)
    for i in range(len(band_centers)):
        plt1.loglog(ell,N_ell_T_LA[0][i], label=str(band_centers[i]) + ' GHz', color=colors[i], ls='-', lw=2.)
        #plt.loglog(ell,N_ell_V3_T_white[i], color=colors[i], ls='-', lw=0.5) #white noise
        
             
    # include correlated atmospheric noise across frequencies
    plt1.loglog(ell, N_ell_T_LA[1][0], label=r'$27 \times 39$ GHz atm.', color='orange', lw=1.5)
    plt1.loglog(ell, N_ell_T_LA[1][1], label=r'$93 \times 145$ GHz atm.', color='fuchsia', lw=1.5)
    plt1.loglog(ell, N_ell_T_LA[1][2], label=r'$225 \times 280$ GHz atm.', color='springgreen', lw=1.5)
    plt1.set_title('$N(\ell$) Temperature', fontsize=18)
    plt1.set_ylabel('$N(\ell$) [$\mu$K${}^2$]', fontsize=16)
    plt1.set_xlabel('$\ell$', fontsize=16)
    plt1.set_ylim(5e-7,1)
    plt1.set_xlim(100,10000)
    plt1.legend(loc='lower left', ncol=2, fontsize=8)
    plt1.grid()
    #plt.savefig('V3_calc_mode'+str(mode)+'_fsky'+str(fsky)+'_defaultdist_noise_LAT_T.pdf')
    #plt.close()
                   
    ## plot the polarization noise curves
    
    for i in range(len(band_centers)):
        plt2.loglog(ell,N_ell_P_LA[0][i], label=str(band_centers[i])+' GHz (V3)', color=colors[i], ls='-', lw=2.)
        #plt.loglog(ell,N_ell_V3_T_white[i], color=colors[i], ls='-', lw=0.5) #white noise
        i+=1
    # include correlated atmospheric noise across frequencies
    plt2.loglog(ell, N_ell_P_LA[1][0], label=r'$27 \times 39$ GHz atm.', color='orange', lw=1.5)
    plt2.loglog(ell, N_ell_P_LA[1][1], label=r'$93 \times 145$ GHz atm.', color='fuchsia', lw=1.5)
    plt2.loglog(ell, N_ell_P_LA[1][2], label=r'$225 \times 280$ GHz atm.', color='springgreen', lw=1.5)
    plt2.set_title(r"$N(\ell$) Polarization", fontsize=18)
    plt2.set_ylabel(r"$N(\ell$) [$\mu$K${}^2$]", fontsize=16)
    plt2.set_xlabel(r"$\ell$", fontsize=16)
    plt2.set_ylim(5e-7,1)
    plt2.set_xlim(100,10000)
    plt2.legend(loc='upper left', ncol=2, fontsize=9)
    plt2.grid()
    #plt.savefig('V3_calc_mode'+str(mode)+'_fsky'+str(fsky)+'_defaultdist_noise_LAT_P.pdf')
    #plt2.close()        
    #fig.show()        
    '''
    ####################################################################
    return(ell, N_ell_T_LA, N_ell_P_LA, Map_white_noise_levels)

#S4_noise()


def package_nl_dic(band_centers, ell, N_ell_T_LA, N_ell_P_LA):
    
    #Spectra Keys
    TParr = ['T','P']
    
    #Power Spectra Keys
    which_spec_arr = ['TT','EE']
    
    #Frequency Bands
    freqarr = [93,145,225,280]
    specs_dic = {}
    
    #initialize N_ell dictionary
    nl_dic = {}
    for TP in TParr:
        nl_dic[TP] = {}
        for freq1 in freqarr:
            for freq2 in freqarr:
                if freq1 == freq2:
                    if TP == 'T':
                        nl_dic[TP][(freq1,freq2)] = N_ell_T_LA[0][freqarr.index(freq1)]
                    elif TP == 'P':
                        nl_dic[TP][(freq1,freq2)] = N_ell_P_LA[0][freqarr.index(freq1)]
                else:
                    if TP == 'T':
                        nl_dic[TP][(freq1,freq2)] = N_ell_T_LA[1][freqarr.index(freq1)]
                    elif TP == 'P':
                        nl_dic[TP][(freq1,freq2)] = N_ell_P_LA[1][freqarr.index(freq1)]


                        
    return nl_dic


