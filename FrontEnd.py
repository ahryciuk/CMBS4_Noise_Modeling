

class FrontEnd():
    
    def __init__(self, exp_dir, cwd, bolocalc_wd):
        
        self.exp_dir = exp_dir
        self.cwd = cwd
        self.bolocalc_wd = bolocalc_wd
        
        return

    def gen_interface(self, exp_dir):

        import ipywidgets as ipw
        from ipywidgets import interact, interact_manual, Layout, Button, Box, HBox, VBox
        import IPython.display as display
        from BoloCalcConverters import TeleCamNames, AppendFiles, runCalcBolos
        from N_ell_calculator import S4_noise
        import numpy as np
        #from ilc_modules import flatsky, misc, ilc
        #from ilc_modules import foregrounds as fg
        import pandas as pd
        import json

        #Run to Generate Simulation Interface



        #Import simulation from Bolocalc to get value of status bar?

        #Instantiate command classes
        run_sim = runCalcBolos.runCalcBolos()
        ap = AppendFiles.AppendFiles()


        ################################################################################# 

        #Title of this block of interface
        title = ipw.HTML("<b><font size=6><font color='green'>Simulations!</font></b>")

        #The button to run a simulation
        run_bolocalc_layout = Layout(width='400px',height='100px')
        run_bolocalc = ipw.Button(description = 'Run calcBolos', layout=run_bolocalc_layout, button_style='primary')

        #Run Simulation button click definition:
        #if a simulation has not been run---button click = run simulation
        #if SimulationCheck.SimulationCheck.(exp_dir) == True:
        #    run_bolocalc.on_click(runCalcBolos.runCalcBolos(exp_dir))
        #else:
        #    print('You have previously ran a simulation')
        #    print('The output files will be overwritten')
        #    print('If you still want to proceed click run calcBolos')
        #    run_bolocalc.on_click(runCalcBolos.runCalcBolos(exp_dir))

        #Text box display of current experimental directory
        try:
            display_exp_path = ipw.Text(description = 'Current Experiment', value=exp_dir)

        except NameError:
            display_exp_path = ipw.Text(description='Current Experiment',value='Null')
            print('Did you enter a path to your experiment?')

        #Progress bar for viewing pleasure
        sim_progress = ipw.IntProgress(value=0,
                                       min=0,
                                       max=100,
                                       step=0.5,
                                       description='Simulation Progress',
                                       bar_style='success',
                                       orientation='horizontal')

        #sim_progress.value = Simulation._bar_len

        #################################################################################

        #Additional buttons
        add_box_layout = Layout(display = 'flex',
                           flex_flow = 'column',
                            justify_content='space-between',
                            align_items = 'stretch',
                           border = 'solid',
                           width = '300px',
                            height = '75px')
        append_input_button = ipw.Button(description='Convert Input Files to Excel',layout=add_box_layout)
        append_input_button.style.button_color='green'
        convert_output_button = ipw.Button(description='Convert Output files to Excel',layout=add_box_layout)
        convert_output_button.style.button_color='orange'
        save_output_button = ipw.Button(description='Overwrite the Input Files',layout=add_box_layout)
        save_output_button.style.button_color='lightblue'
        plots_button = ipw.Button(description='Save Plot Parameters From Excel',layout=add_box_layout)
        plots_button.style.button_color='red'

        #################################################################################

        #save to path and name
        save_layout = Layout(width = '200px',height='75px')
        save_in_button = ipw.Button(description='Save Input Files',layout=save_layout, button_style='info')
        save_out_button = ipw.Button(description='Save Output Files',layout=save_layout, button_style='info')
        save_buttons = VBox([save_in_button,save_out_button])

        save_in_to = ipw.Text(description='Save Input to:',placeholder='Path to save to')
        save_in_name = ipw.Text(description='Input Directory Name:')
        save_out_to = ipw.Text(description='Save Output to:',placeholder='Path to save to')
        save_out_name = ipw.Text(description='Output Directory Name:')
        save_to = VBox([save_in_to, save_in_name, save_out_to, save_out_name])

        save_int = HBox([save_buttons, save_to])
        save_int.layout.margin = '0 0 0 275px'

        add_buttons = HBox([VBox([append_input_button,convert_output_button]),VBox([save_output_button,plots_button])])
        add_buttons.layout.margin = '0 0 0 175px'

        #################################################################################

        #Button Click Execution Commands:
        save_in_button.on_click(self.savefiles_in)
        save_out_button.on_click(self.savefiles_out)
        append_input_button.on_click(self.appfiles)
        save_output_button.on_click(self.appinput)
        convert_output_button.on_click(self.convout)
        run_bolocalc.on_click(self.runsim)
        plots_button.on_click(self.saveplotparams)

        #################################################################################

        #Formatting and Display
        display_col = VBox([display_exp_path,sim_progress])

        run_sim_row = HBox([run_bolocalc,display_col])
        run_sim_row.layout.margin = '25px 0 0 150px'
        total_layout=Layout(border='solid 15px blue')
        total = VBox([title, run_sim_row, add_buttons, save_int],layout=total_layout,)
        total.layout.margin = '0px 0px 0px 0px'



        return total



    ################################################################################# 

    #Button Click Definitions:
    def appfiles(self, a):
        from BoloCalcConverters import TeleCamNames, AppendFiles, runCalcBolos
        import IPython.display as display
        display.clear_output()
        ap = AppendFiles.AppendFiles()

        try:
            ap.InputConvert(self.exp_dir)
        except NameError:
            print('Please Enter an Experimental Path')
        return
    def appinput(self, a):
        from BoloCalcConverters import TeleCamNames, AppendFiles, runCalcBolos
        import IPython.display as display
        display.clear_output()
        ap = AppendFiles.AppendFiles()

        try:
            ap.AppendInputs(self.exp_dir)
        except NameError:
            print('Please Enter an Experimental Path')
        return
    def convout(self, a):
        from BoloCalcConverters import TeleCamNames, AppendFiles, runCalcBolos
        import IPython.display as display
        display.clear_output()
        ap = AppendFiles.AppendFiles()

        try:
            ap.ConvertOutputFiles(self.exp_dir)
        except NameError:
            print('Please Enter an Experimental Path')
        return
    def savefiles_in(self, a):
        from BoloCalcConverters import TeleCamNames, AppendFiles, runCalcBolos
        import IPython.display as display
        display.clear_output()
        ap = AppendFiles.AppendFiles()

        try:
            a=0
            ap.SaveFiles(self.exp_dir, save_in_to.value, save_in_name.value, a)
        except NameError:
            print('Please Enter an Experimental Path')
        return
    def savefiles_out(self, a):
        from BoloCalcConverters import TeleCamNames, AppendFiles, runCalcBolos
        import IPython.display as display
        display.clear_output()
        
        ap = AppendFiles.AppendFiles()

        try:
            a=1
            ap.SaveFiles(self.exp_dir, save_out_to.value, save_out_name.value, a)
        except NameError:
            print('Please Enter an Experimental Path')
        return
    def runsim(self, a):
        from BoloCalcConverters import TeleCamNames, AppendFiles, runCalcBolos
        import IPython.display as display
        display.clear_output()
        
        run_sim = runCalcBolos.runCalcBolos()
        try:
            run_sim.runSim(self.exp_dir,self.bolocalc_wd)
        except NameError:
            print('Enter a path to your experiment.')
        return
    def saveplotparams(self, a):
        from BoloCalcConverters import TeleCamNames, AppendFiles, runCalcBolos
        import pandas as pd
        import os
        #ap = AppendFiles.AppendFiles()

        import json
        import IPython.display as display
        display.clear_output()
        #try:
        print('Saving plot parameters. Shift-Enter the next cell to display')
        os.chdir(self.exp_dir)# + '/' + 'OutputExcelFiles')
        plot_params_set = pd.read_excel('OutputExcelFiles.xlsx',sheet_name='N_ell_Plotting_Parameters')

        global plot_params_names
        global plot_params_values
        plot_params_names=[]
        plot_params_values=[]
        for i in range(len(plot_params_set)):
            plot_params_names.append(plot_params_set.iloc[i,0])
            plot_params_values.append(plot_params_set.iloc[i,1])

        for i in range(len(plot_params_values)):
            if type(plot_params_values[i]) is str:
                plot_params_values[i] = json.loads(plot_params_values[i])
        print('Done!')
        #except NameError:
        #    print('Enter a path to your experiment.')
        return
    
    def N_ell_calculator(self, a):
        
        import matplotlib.pyplot as plt
        import matplotlib as pltt
        import subprocess
        import IPython.display as display
        from N_ell_calculator import S4_noise
        
        display.clear_output()

        global ell, N_ell_T_LA, N_ell_P_LA, Map_white_noise_levels
        ell, N_ell_T_LA, N_ell_P_LA, Map_white_noise_levels, corr_freq = S4_noise(band_centers=plot_params_values[0],
                                                                  beam_sizes=plot_params_values[1],Sens=plot_params_values[15],
                                                                  f_knees=plot_params_values[2],Cs=plot_params_values[3],
                                                                  alpha_temp=plot_params_values[4],
                                                                   survey_time=plot_params_values[5],f_sky=plot_params_values[6],
                                                                   ret_after_obs_cuts=plot_params_values[7],
                                                                   non_uniformity_param=plot_params_values[8],
                                                                   ell_max=plot_params_values[9],ell_pivot=plot_params_values[10],
                                                                   delta_ell=plot_params_values[11],
                                                                   alpha_pol=plot_params_values[12],NTubes=plot_params_values[13],
                                                                      model_num=plot_params_values[14])


        global band_centers
        band_centers = plot_params_values[0]

        colors = ['b','r','g','m','k','y']

        fig, (plt1,plt2) = plt.subplots(2,1,figsize=(10,8))
        fig.tight_layout(pad=3.0)

        #TParr = ['T','P']
        for freq1 in band_centers:
            for freq2 in band_centers:

                if freq1 == freq2:
                    i = band_centers.index(freq1)
                    plt1.loglog(ell, N_ell_T_LA[(freq1,freq2)],label=str(band_centers[i])+' GHz', color=colors[i], ls='-', lw=2.)
                else:
                    i = band_centers.index(freq1)
                    j = band_centers.index(freq2)
                    if freq2 is corr_freq[freq1] and freq2 > freq1:
                        plt1.loglog(ell, N_ell_T_LA[(freq1,freq2)],label=str(band_centers[i])+
                                'x'+str(band_centers[j]) + ' GHz', color=colors[i],lw=1.5)

        #for i in range(0,len(band_centers)):
        #    plt1.loglog(ell,N_ell_T_LA[0][i], label=str(band_centers[i]) + ' GHz', color=colors[i], ls='-', lw=2.)
            #plt.loglog(ell,N_ell_V3_T_white[i], color=colors[i], ls='-', lw=0.5) #white noise


        # include correlated atmospheric noise across frequencies
        #plt1.loglog(ell, N_ell_T_LA[1][0], label=r'$27 \times 39$ GHz atm.', color='orange', lw=1.5)
        #plt1.loglog(ell, N_ell_T_LA[1][1], label=r'$93 \times 145$ GHz atm.', color='fuchsia', lw=1.5)
        #plt1.loglog(ell, N_ell_T_LA[1][2], label=r'$225 \times 280$ GHz atm.', color='springgreen', lw=1.5)
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

        #for i in range(len(band_centers)):
        #   plt2.loglog(ell,N_ell_P_LA[0][i], label=str(band_centers[i])+' GHz (V3)', color=colors[i], ls='-', lw=2.)
            #plt.loglog(ell,N_ell_V3_T_white[i], color=colors[i], ls='-', lw=0.5) #white noise
        #    i+=1
        # include correlated atmospheric noise across frequencies

        for freq1 in band_centers:
            for freq2 in band_centers:

                if freq1 == freq2:
                    i = band_centers.index(freq1)
                    plt2.loglog(ell, N_ell_P_LA[(freq1,freq2)],label=str(band_centers[i])+' GHz', color=colors[i], ls='-', lw=2.)

                else:
                    i = band_centers.index(freq1)
                    j = band_centers.index(freq2)
                    if freq2 is corr_freq[freq1] and freq2 > freq1:
                        plt2.loglog(ell, N_ell_P_LA[(freq1,freq2)],label=str(band_centers[i])+
                                'x'+str(band_centers[j]) + ' GHz', color=colors[i],lw=1.5)


        #plt2.loglog(ell, N_ell_P_LA[1][0], label=r'$27 \times 39$ GHz atm.', color='orange', lw=1.5)
        #plt2.loglog(ell, N_ell_P_LA[1][1], label=r'$93 \times 145$ GHz atm.', color='fuchsia', lw=1.5)
        #plt2.loglog(ell, N_ell_P_LA[1][2], label=r'$225 \times 280$ GHz atm.', color='springgreen', lw=1.5)
        plt2.set_title(r"$N(\ell$) Polarization", fontsize=18)
        plt2.set_ylabel(r"$N(\ell$) [$\mu$K${}^2$]", fontsize=16)
        plt2.set_xlabel(r"$\ell$", fontsize=16)
        plt2.set_ylim(5e-7,1)
        plt2.set_xlim(100,10000)
        plt2.legend(loc='upper left', ncol=2, fontsize=9)
        plt2.grid()

        return
    
    def gen_nl_interface(self):

        import ipywidgets as ipw
        from ipywidgets import interact, interact_manual, Layout, Button, Box, HBox, VBox
        import IPython.display as display
        
        #############Text Box inputs to calculation################

        try:
            params = [0 for i in range(len(plot_params_names))]
            for i in range(len(plot_params_names)):
                params[i] = ipw.Text(description=plot_params_names[i],value=str(plot_params_values[i]))
        except NameError:
            params = [0 for i in range(18)]
            for i in range(len(params)):
                params[i] = ipw.Text(description='')
            print('Save your plot parameters with the red button in the previous section')
            print('Then try Shift-Enter on this cell again.')

        ####################Display###############################

        params_box1=HBox(params[0:3])
        params_box2=HBox(params[3:6])
        params_box3=HBox(params[6:9])
        params_box4=HBox(params[9:12])
        params_box5=HBox(params[12:15])
        params_box6=HBox(params[15:])
        params_box = VBox([params_box1,params_box2,params_box3,params_box4,params_box5,params_box6])

        plot_layout = Layout(width='600px',height='100px')
        plot_button = ipw.Button(description = 'Plot N(ell)', layout=plot_layout, button_style='danger')
        plot_button.layout.margin = '0px 0px 0px 175px'
        plot_button.on_click(self.N_ell_calculator)

        ##################Run Interact button#####################
        plots = VBox([params_box,plot_button])

        return plots
    
    def gen_cl_int(self):
        
        import os
        import ipywidgets as ipw
        from ipywidgets import interact, interact_manual, Layout, Button, Box, HBox, VBox
        import IPython.display as display
        os.chdir(self.cwd + '/ilc_modules')
        import misc, ilc, flatsky, foregrounds as fg
        os.chdir(self.cwd)
        #Need to pass Nl's to make dictionary, make bl dictionary, ell array
        from N_ell_calculator import package_beam_dic, package_nl_dic
        




        ###############Interface###############################
        #Interface
        #params
        
        #os.chdir(self.cwd + '/ilc_modules')
        

        #Layout of buttons
        calc_cl_layout = Layout(display = 'flex',
                           flex_flow = 'column',
                            justify_content='space-between',
                            align_items = 'stretch',
                           border = 'solid',
                           width = '600px',
                            height = '100px')
        add_box_layout = Layout(display = 'flex',
                           flex_flow = 'column',
                            justify_content='space-between',
                            align_items = 'stretch',
                           border = 'solid',
                           width = '300px',
                            height = '75px')
        convert_parameters = ipw.Button(description='Convert Comp Sep Params to Excel',layout=add_box_layout)
        append_parameters = ipw.Button(description='Overwrite Comp Sep Text File',layout=add_box_layout)
        append_parameters.style.button_color='lightblue'
        convert_parameters.style.button_color='violet'

        calc_cl = ipw.Button(description='Calculate Cl Component Spectra',layout=calc_cl_layout,button_style='primary')
        calc_cl.style.button_color='green'

        #Define clicking
        convert_parameters.on_click(self.convert_comp_sep)
        append_parameters.on_click(self.app_comp_sep)
        calc_cl.on_click(self.calc_cl_spec)


        param_buttons = HBox([convert_parameters,append_parameters])

        total_cls = VBox([param_buttons,calc_cl])
        total_cls.layout.margin = '0px 0px 0px 175px'

        return total_cls
        ###############################################
    
    def convert_comp_sep(self,a):
        import IPython.display as display
        display.clear_output()
        convert_params_excel(self.exp_dir,self.cwd)
        return
    def app_comp_sep(self,a):
        import IPython.display as display
        display.clear_output()
        append_params(self.exp_dir,self.cwd)
        return
    def calc_cl_spec(self,a):
        
        import IPython.display as display
        display.clear_output()
        
        paramfile = 'params.ini'
        
        import os
        os.chdir(self.cwd + '/ilc_modules')
        import misc, ilc, flatsky, foregrounds as fg
        import plot_cls
        from ConvertILC import convert_params_excel, append_params, calc_cls
        param_dict = misc.fn_get_param_dict(paramfile)
        os.chdir(self.cwd)
        import matplotlib
        
        #try:
        band_centers = plot_params_values[0]
        sensitivities = plot_params_values[15]
        beam_sizes = plot_params_values[1]
        ell_max = plot_params_values[9]
        nl_dic = {'T':N_ell_T_LA,'P':N_ell_P_LA}
        #nl_dic['T'][(145,93)] = nl_dic['T'][(93,145)]
        #nl_dic['T'][(280,225)] = nl_dic['T'][(225,280)]
        #nl_dic['P'][(145,93)] = nl_dic['P'][(93,145)]
        #nl_dic['P'][(280,225)] = nl_dic['P'][(225,280)]
        #bl_dic = package_beam_dic(white_noise_levels_T=sensitivities,
        #                          white_noise_levels_P=sensitivities,
        #                          band_centers = band_centers,
        #                         beam_sizes=beam_sizes,ellmax=ell_max)
        bl_dic = {}
        for freq in band_centers:
            i = band_centers.index(freq)
            bl_dic[freq] = misc.get_bl(beam_sizes[i],ell)

        global cl_residual
        el, cl_dic, weights_dic, cl_residual = calc_cls(nl_dic, bl_dic, ell, band_centers, self.cwd)



        plot_cls.plot_cls(param_dict, band_centers, el, self.cwd, cl_residual,weights_dic, nl_dic)
        #except NameError:
        #    print('Did you calculate the N_ell?')


        return