

class systematics():
    
    
    
    def __init__(self, fsky, freq_bands, fwhm_bands, f3db_beams, ell):
        
        self.fsky = fsky    #experiment sky fraction
        self.freq_bands = freq_bands    #array of frequency bands
        self.fwhm_bands = fwhm_bands    #array of beam sizes
        self.f3db_beams = f3db_beams    #array of detector 3db values
        self.ell = ell    #array of ell values for simulation
        
        return
    
    
    
    def get_beam_unc(self, freq_bands, fwhm_bands, f3db_beams, pct_unc_lpf):
        
        #percent difference 
        pct_diff = {}  
        #iterate over band frequencies
        for i in range(len(fwhm_bands)):
            
            #convert pct_unc_beam to f3db
            f3db = [f3db_beams[i] , (f3db_beams[i] - pct_unc_lpf * f3db_beams[i])]
        
            #calculate get_tc_uncertainty for each band
            pct_diff[freq_bands[i]] = self.get_tc_uncertainty(self.ell, self.fsky, self.fwhm_bands[i], f3db)
            
                
        #put in dictionary structure like nl_dic
        
        #return dictionary        
        
        return pct_diff
    
    def transfer_function_error(cl_dic, pct_diff):
        
        assert (len(cl_dic) == len(pct_diff)),"Dimensions of Cls and Beam are not the same!"
        
        for TP in cl_dic.keys():            
            cl_dic[TP] = cl_dic[TP] * (1. + pct_diff)
        
        return cl_dic
    
    def get_tc_uncertainty(self, l, fsky, fwhm, f3db):
        
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.optimize import curve_fit
        
        #parameters
        el= 45. #scan elevation in degrees
        fscan= 1.#scan speed of the telescope deg/s
        #fsky=fscan*np.cos(el*np.pi/180.) #deg/s on sky scan speed
        
        
        #is fsky correct?
        
        
        #frequency space Gaussian beam
        sigma=fwhm/(60.*2.*fsky*np.sqrt(2.*np.log(2)))
        t_range=100000.
        t_points=10000000
        t_freq=np.linspace(-t_range,t_range,t_points)
        sample_rate=(t_points/(2.*t_range))
        gauss_2=np.exp(-2.*np.pi**2*t_freq**2*sigma**2)

        pshift=np.zeros(len(f3db))
        FWHM_new=np.zeros(len(f3db))
        for ii in range(len(f3db)):

                #Make Lowpass Filter (and apply to negative side--this is just the scan direction)
                h=(1.-1j*(-t_freq/f3db[ii]))/(1.+np.square(-t_freq/f3db[ii]))
                filtered_2=np.copy(gauss_2)*h
                filtered=np.copy(filtered_2)

                #Take FFT of Convoluted Beam
                gffted=np.abs(np.fft.fft(gauss_2))
                ffted=np.fft.fft(filtered)
                rffted=np.absolute(ffted)
                #make position array of the same length
                pos1=np.linspace(0,len(rffted)-1,len(rffted))
                #normalize to one and scale by sample rate and convert into degrees
                pos2=(pos1/len(rffted))*sample_rate*fsky

                #put the plot back together
                #First section of fft
                fftr1=rffted[0:int(len(rffted)/2)]
                #print pointing offset
                pshift[ii]=pos2[np.argmax(fftr1)]*60.
                #print f3db[ii], 'Hz', pshift[ii], 'arcmin'
                #Find x value where value is half
                x=np.copy(pos2[0:int(len(pos2)/2)])
                yf=np.copy(fftr1)
                y=yf/np.max(yf)
                xp=np.fliplr([x])[0]
                yp=np.fliplr([y])[0]
                half=np.interp([0.5],yp,xp)
                FWHM_new[ii]=(half[0]*60.-pshift[ii])*2.
                #print "Beam FWHM: ",FWHM_new[ii], "arcmin"
                #print "Change in FWHM: ",FWHM_new[ii]-fwhm, "arcmin"

        # calculate error in pointing from time constant uncertainty
        p_10=pshift[1]-pshift[0]
        #print "Pointing error: ", p_10*60., "arcsec"
        #print "Pointing error: ", p_10*60., "arcsec"
        #Take units in arcminutes
        fwhm_shift=np.zeros(len(f3db))
        
        for jj in range(len(f3db)):
            #print '%0.2f Hz:'%f3db[jj]
            angle=np.linspace(-500,500,100000)
            poff=pshift[jj]
            sigma2=FWHM_new[jj]/(2.*np.sqrt(2.*np.log(2)))
            a=1./(sigma2*np.sqrt(2*np.pi))
            top1=-0.5*((angle-poff)**2/sigma2**2)
            top2=-0.5*((angle+poff)**2/sigma2**2)
            #These are the two gaussians
            gauss_a=a*np.exp(top1)
            gauss_b=a*np.exp(top2)

            #Make the combined beam
            gauss=gauss_a+gauss_b

            #Find FWHM and thus the best guess of sigma (assume centered around zero)
            aa=np.copy(angle)
            yyf=np.copy(gauss)
            yy=yyf/np.max(yyf)
            half=np.interp([0.5],yy[0:int(len(angle)/2)],aa[0:int(len(angle)/2)])
            fwhm_fit=np.absolute(2*half[0])
            #print 'Guess: ',fwhm_fit
            sigfit=fwhm_fit/(2.*np.sqrt(2.*np.log(2)))

            popt,pcov = curve_fit(self.gaus,angle,gauss,p0=[1,sigfit])
            fwhm_shift[jj]=popt[1]*(2.*np.sqrt(2.*np.log(2)))
            #print 'FWHM Fit: ',fwhm_shift[jj], "arcmin"
            #print 'Change in FWHM: ',popt[1]*(2.*np.sqrt(2.*np.log(2)))-fwhm
            #take these and look at their Gaussian window functions
            
        #l=np.linspace(0,50000,50001)
        sigma_nt=(fwhm/60.)*(np.pi/180.)/np.sqrt(8*np.log(2))
        #Gaussian window function, no time constant
        bl=np.exp(-l*(l+1)*sigma_nt**2)
        #with base time constant
        sigma_tc=(fwhm_shift[0]/60.)*(np.pi/180.)/np.sqrt(8*np.log(2))
        blm_orig=np.exp(-l*(l+1)*sigma_tc**2)
        #base time constant + err
        sigma_mod=(fwhm_shift[1]/60.)*(np.pi/180.)/np.sqrt(8*np.log(2))
        blm=np.exp(-l*(l+1)*sigma_mod**2)
        #pct_diff
        pct_diff_tc=(blm_orig-blm_orig)/bl
        pct_diff=(blm-blm_orig)/bl
        
        
        return pct_diff
  
    
    
    #Fit Gauss to a Gaussian
    def gaus(self,x,a,sig):
        import numpy as np
        return a*np.exp(-(x)**2/(2.*sig**2))
    
    
    def oneovereff(self, ell, A, ell_knee, alpha_ind):
        import numpy as np
        
        #oneovereff_spec = np.zeros(len(ell))
        oneovereff_spec = A*(ell_knee/ell)**alpha_ind      

        return oneovereff_spec
    
    def oneovereff_readout(self, ell, NEI, S_responsivity, ell_knee, alpha_ind):
        import numpy as np
        
        #convert to white noise level
        white_noise = NEI/S_responsivity
        
        #call generic oneovereff function
        readout = self.oneovereff(ell, white_noise, ell_knee, alpha_ind)
       
        
        #normalize to white noise level and add in quadrature
        readout_spec = np.sqrt(readout**2. + white_noise**2.)
        
        return readout_spec
    
    
    def oneovereff_atm(self, ell, alpha_ind, ell_knee, white_noise_level):
        import numpy as np
        
        #call generc oneovereff function
        atm = self.oneovereff(ell, white_noise_level, ell_knee, alpha_ind)
        
        
        #normalize to white noise level
        atm_spec = np.sqrt(atm**2. + white_noise_level**2.)
        
        return atm_spec
    
    
    def oneovereff_det(self, ell, alpha_ind, ell_knee, white_noise_level):
        import numpy as np
        
        #call generic oneovereff function
        det = self.oneovereff(ell, white_noise_level, ell_knee, alpha_ind)
        
        #normalize to white noise level
        det_spec = np.sqrt(det**2. + white_noise_level**2.)
        
        
        return det_spec
    
    def generate_oneovereff_spectra(self, A_list, ell_knee_list, alpha_ind_list):
        
        oneovereff_spectra = {}
        for freq in self.freq_bands:
            i = self.freq_band.index(freq)
            oneovereff_spectra[freq] = self.oneovereff(self.ell, A_list[i], ell_knee_list[i], alpha_ind_list[i])
            
            
        #return all spectra?
        return oneovereff_spectra
    
    def scale_nl(self, nl_dic, scale_dic):
        
        TParr = ['T','P']
        
        for TP in TParr:
            for freq1 in self.freq_bands:
                for freq2 in self.freq_bands:
                    nl_dic[TP][(freq1,freq2)] = scale_dic[TP][(freq1,freq2)] * nl_dic[TP][(freq1,freq2)]
        
        return nl_dic
    
    def inv_var_weight_nl(self, freq_bands, nl_dic, which_spectra):
        
        nl_inv_var_dic = {}
        for freq1 in freq_bands:
            nl_inv_var_dic[which_spectra] += 1. / nl_dic[(freq1,freq1)]
            
        nl_inv_var_dic[which_spectra] = 1. / nl_inv_var_dic[which_spectra]
        
        return nl_inv_var_dic
    
    
    