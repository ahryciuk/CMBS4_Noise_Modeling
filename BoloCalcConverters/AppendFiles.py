#All functions related to appending files

##############################Append input files######################################

class AppendFiles():

    def __init__(self):
        
        
        return

    ##ConfigReader##
    def ConfigReader(self, file_name_input, file_name_output, writer, file_path_output, sep_p=None):

        import os
        import pandas as pd
        
        '''
        Reads in a text file to parse data in a dataframe. This dataframe can be written to
        it's own excel file, or appended to an existing excel file
        
        file_name_input: name of file assuming already in the same directory
        file_name_output: name for the excel file output to be written to some directory
        writer: the writer object that specifies the excel file or sheet to write to
        file_path_output: the path to where the excel file will be saved
        '''
        
        
        #store input as a dataframe
        if sep_p is None:
            inputs = pd.read_csv(file_name_input, sep = '|', comment = '#',error_bad_lines=False)
        else:
            inputs = pd.read_csv(file_name_input, sep = sep_p, comment = '#',error_bad_lines=False)
            
        os.chdir(file_path_output)
        inputs.to_excel(writer, sheet_name = file_name_output)


        return
    
    ##Convert Inputs Parameters to Excel##
    def InputConvert(self, exp_dir):

       
        print('Converting Inputs to Excel......')

        import os,shutil
        import pandas as pd
        import subprocess, sys
        from BoloCalcConverters import TeleCamNames

        '''
        This module walks through all the config directories in the passed
        experimental directory and converts them to a master spreadsheet.
        The master spreadsheet can be modified by the user and mapped back
        to append the original config text files in the experimental directory 
        with the AppendInputs module below.

        exp_dir: the experimental directory passed by the interface.
        '''

        telescope_names, camera_names = TeleCamNames.TeleCamNames(exp_dir)

        #try:
        #    os.mkdir(exp_dir + '/' +  'InputExcelParameters')
        #except FileExistsError:
        #    print('Overwriting lastest Input directory')
        #    shutil.rmtree(exp_dir + '/' + 'InputExcelParameters')
        #    os.mkdir(exp_dir+'/'+'InputExcelParameters')

        file_path_output = exp_dir# + '/' + 'InputExcelParameters'
        InputExcel = 'InputExcelParameters.xlsx'

        os.chdir(exp_dir + '/' + 'config')

        with pd.ExcelWriter(InputExcel, engine = 'xlsxwriter') as writer:

            #plotting parameter inputs (hardcoded from LAT model)
            params = ['band_centers','beam_sizes','f_knees','Cs','alpha_temp',
               'survey_time','f_sky','ret_after_obs_cuts','non_uniformity_param','ell_max','ell_pivot',
                      'delta_ell','alpha_pol','NTubes','model_num']

            #Hard coded from SO LAT model
            default_values = {'LAT Model Parameters':[[27.,39.,93.,145.,225.,280.],[7.4,5.1,2.2,1.4,1.0,0.9], [700.,700.,700.,700.,700.,700.],[200,7.7,1800,12000,68000,124000],[-3.5,-3.5,-3.5,-3.5,-3.5,-3.5], 5., 0.4, 0.2,0.85,1e4,[1000.,1000.,1000.,1000.,1000.,1000.],1,[-0.4,-0.4,-0.4,-0.4,-0.4,-0.4],[1,1,1,1,1,1],1]}

            #store the previous LAT model parameters in a dataframe  and write to the master sheet
            plot_params = pd.DataFrame(default_values,index=params)
            os.chdir(exp_dir)# + '/' + 'InputExcelParameters')
            plot_params.to_excel(writer, sheet_name = 'Atmosphere_Model_Parameters')

            os.chdir(exp_dir + '/' + 'config')
            self.ConfigReader('foregrounds.txt', 'foregrounds', writer, file_path_output)
       
            #Telescope Directories#
            for i in range(len(telescope_names)):
                for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + 'config'):
                    for k in range(len(fileList)):
                        if os.path.exists(exp_dir + '/' + telescope_names[i] + '/config/Dist/' + fileList[k]):
                            pass
                        else:
                            if fileList[k] == 'telescope.txt':
                                os.chdir(exp_dir + '/' + telescope_names[i] + '/' + 'config')
                                self.ConfigReader(fileList[k], fileList[k].rstrip('.txt') + '_' + telescope_names[i], writer, file_path_output)
                
                #Camera Directories# 
                for j in range(len(camera_names[i])):
                    for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/' + 'config'):
                        for k in range(len(fileList)):
                            if os.path.exists(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/config/Dist/Detectors/' + fileList[k]):
                                pass
                            elif os.path.exists(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/config/Dist/Optics/' + fileList[k]):
                                pass
                            elif os.path.exists(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/config/Bands/Detectors/' + fileList[k]):
                                pass
                            elif os.path.exists(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/config/Bands/Optics/' + fileList[k]):
                                pass
                            else:
                                if fileList[k] == 'optics.txt' or fileList[k] == 'camera.txt' or fileList[k] == 'channels.txt':
                                    os.chdir(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/' + 'config')
                                    self.ConfigReader(fileList[k], fileList[k].rstrip('.txt') + '_' + telescope_names[i] + '_' + camera_names[i][j], writer, file_path_output)
                            
            writer.save()
            print('Done!')
            print('\n')
                
     
        return

    def AppendConfigFiles(self, file_name, sheet_name, output_path, output_file_name):
        import pandas as pd
        import os
    
        #input file to be converted to excel
        inputs = pd.read_excel(file_name, sheet_name)
        inputs.drop(inputs.columns[[0]],axis=1,inplace=True)

        #Change path to whatever the specified output is
        os.chdir(output_path)
        if not os.path.isfile(output_file_name):
            inputs.to_csv(output_file_name, sep = '|', index = False)
        else:
            inputs.to_csv(output_file_name, sep = '|', mode = 'w', index = False)
       
        return

    

    ##Function to walk through directory##
    def AppendInputs(self, exp_dir):
        
        import os
        import pandas as pd
        from BoloCalcConverters import TeleCamNames
        
        '''
        This module takes the InputExcelFiles master spreadsheet, which has been
        appended for a new BoloCalc calculation and maps those back to the original
        config files that define the simulation parameters. From there, a new BoloCalc
        simulation can be run. 
        
        exp_dir: the experimental directory that contains the InputExcelFiles directory
        '''
        
        #import the names of the telescopes and cameras from the given directory
        telescope_names, camera_names = TeleCamNames.TeleCamNames(exp_dir)
        
        try:
            print('Appending the Input Config Text Files')

            #InputExcelParameters Path
            input_path = exp_dir# + '/' + 'InputExcelParameters'
            
            #Change the directory path to where the inputs have been appended
            os.chdir(input_path)
            file_name = 'InputExcelParameters.xlsx' 
            
            #Append the foregrounds text file 
            self.AppendConfigFiles(file_name, 'foregrounds', 'config', 'foregrounds.txt')

            #Loop through all the telescopes and camera directories and append each text file from the corresponding sheet
            for i in range(len(telescope_names)):
                for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + 'config'):
                    for k in range(len(fileList)):
                        if fileList[k] == 'telescope.txt':
                            output_file_name = fileList[k]
                            output_path = exp_dir + '/' + telescope_names[i] + '/' + 'config'
                            sheet_name = fileList[k].rstrip('.txt') + '_' + telescope_names[i]
                            os.chdir(input_path)
                            self.AppendConfigFiles(file_name, sheet_name, output_path, output_file_name)
                for j in range(len(camera_names[i])):
                    for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/' + 'config'):
                        for k in range(len(fileList)):
                            if fileList[k] == 'optics.txt' or fileList[k] == 'camera.txt' or fileList[k] == 'channels.txt':
                                output_file_name = fileList[k]
                                output_path = exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/' + 'config'
                                sheet_name = fileList[k].rstrip('.txt') + '_' + telescope_names[i] + '_' + camera_names[i][j]
                                os.chdir(input_path)
                                self.AppendConfigFiles(file_name, sheet_name, output_path, output_file_name)

            print('Done!')
            print('\n')
            
        #If the user has yet to convert the inputs to excel and append them    
        #except (FileExistsError, FileNotFoundError):
        except NameError:
            print('Have you appended the inputs in excel?')
            print('\n')
            
        return

#############################Saves Files to some directory#########################################
  
    def SaveFiles(self, exp_dir, save_path, save_name, in_or_out):
        import shutil, os
        
        '''
        This module allows the user to save the InputExcelFiles folder to another place
        on their machine in the case where multiple simulation want to be run and the 
        inputs/outputs are documented for later use. The InputExcelFiles directory
        is then copied into the path specified by save_path under the name specified
        by save_name.
        
        exp_dir: the experimental directory containing the OutputExcelFiles directory
        save_path: path on this machine to copy and save the OutputExcelFiles directory
        save_name: name for the OutputExcelFiles directory to be under
        '''
        
        try:
            if not os.path.exists(save_path + '/' + save_name):
                print('Saving your file to ' + save_path + ' as ' + save_name)
                
                if in_or_out == 0:
                    os.chdir(exp_dir)
                    shutil.copytree('InputExcelParameters', save_path + '/' + save_name)
                elif in_or_out == 1:
                    os.chdir(exp_dir)
                    shutil.copytree('OutputExcelFiles', save_path + '/' + save_name)
                else:
                    print('Something went wrong saving your file')
            else:
                print('This file name already exists in this directory...Care to try another?')
        
        except FileExistsError:
            print('It appears the name you entered already exists...Care to try another?')
        except FileNotFoundError:
            print('Either the directory does not exist of you have yet to enter what the file name is')
   
        return

    def SaveFile_Out(self, exp_dir, save_path, save_name):
        import shutil, os
        
        '''
        This module allows the user to save the OutputExcelFiles folder to another place
        on their machine in the case where multiple simulation want to be run and the 
        inputs/outputs are documented for later use. The OutputExcelFiles directory
        is then copied into the path specified by save_path under the name specified
        by save_name.
        
        exp_dir: the experimental directory containing the OutputExcelFiles directory
        save_path: path on this machine to copy and save the OutputExcelFiles directory
        save_name: name for the OutputExcelFiles directory to be under
        '''
        
        try:
            if not os.path.exists(save_path + '/' + save_name):
                print('Saving your file to ' + save_path + ' as ' + save_name)
            
                os.chdir(exp_dir + '/' + 'OutputExcelFiles')
                shutil.copy('OutputExcelFiles.xlsx', save_path + '/' + save_name)
            else:
                print('This file name already exists in this directory...Care to try another?')
        
        except FileExistsError:
            print('It appears the name you entered already exists...Care to try another?')
        except FileNotFoundError:
            print('That directory does not seem to exist...Care to try another?')
        
        return
    
    
################################Convert all the Output Files functions#########################################

    
    def WritetoSheet(self, inputs, file_name_out, file_path_out, writer):
        import pandas as pd
        import os
        
        '''
        This module takes an input dataframe (from conversion modules below)
        and writes the input to a master sheet specified by writer. 
        
        inputs: the input dataframe from conversion modules
        file_name_out: name of the sheet
        file_path_out: path to excel sheet
        writer: master spreadsheet to be appended with input dataframe
        '''
        
        os.chdir(file_path_out)
        if type(inputs)==list:
            for i in range(len(inputs)):
                inputs[i].to_excel(writer, sheet_name=file_name_out + '_' + str(i+1))
        else:
            inputs.to_excel(writer, sheet_name=file_name_out)        
        
        return
    
    def ConvertOutputFiles(self, exp_dir):
        
        from BoloCalcConverters import TeleCamNames
        import os
        import pandas as pd
        
        '''
        This module will take in an experimental directory that has calculated the outputs
        from a BoloCalc simulation. The code will walk through the experimental directory
        that is passed by the interface and use the conversion modules below to convert all
        the output files to write all the data to a spreadsheet directory as well as a master
        spreadsheet containing all the output files named OutputExcelFiles.xlsx. In addition,
        another speadsheet is created with modifiable parameters that will be mapped to later
        instrument models within the interface of BoloCalc. The user can walk through all the
        output files and pick what parameters are interesting for later models. 
        
        exp_dir: the experimental directory containing the standard BoloCalc simulation files.
        It is assumed that the directory has run through a BoloCalc simulation calculation.
        '''
        
        #This calls the module TeleCamNames, which walks through the directory and stores
        #lists containing the telescope and camera names
        telescope_names, camera_names = TeleCamNames.TeleCamNames(exp_dir)
        
        #error handling for if the files have already been stored
      
        print('Converting Outputs to Excel. If a warning pops up do not fear.')

        #loop thru camera files to extract output files

        #master spreadsheet name
        OutputExcel = 'OutputExcelFiles.xlsx'

        #instantiating variables for loop
        path_in = 0
        path_out = 0

        #get sensitivities from LAT output files
        os.chdir(exp_dir + '/' + telescope_names[0])

        f = open('sensitivity.txt','r')
        f_read = f.readlines()

        col_header = f_read[0].split('|')
        for i in range(len(col_header)):
            col_header[i] = col_header[i].strip(' ')


        sens_frame = pd.read_csv('sensitivity.txt',sep='|',comment='-',skiprows=[0,1,2,-1],names=col_header)

        sens_list = list(sens_frame['Array NET_CMB'])
        for i in range(len(sens_list)):
            sens_list[i] = sens_list[i].rstrip(' +/')
            sens_list[i] = sens_list[i].lstrip(' ')
            sens_list[i] = float(sens_list[i])
        sens_list.pop(-1)

        #make directory to store spreadsheets within the experimental directory
        #try:
        #    os.mkdir(exp_dir + '/' + 'OutputExcelFiles')
        #except FileExistsError:
        #    print('Overwriting latest Output Directory')
        #    shutil.rmtree(exp_dir + '/' + 'OutputExcelFiles')
        #    os.mkdir(exp_dir + '/' + 'OutputExcelFiles')

        #spreadsheet writer for master spreadsheet
        with pd.ExcelWriter(OutputExcel, engine = 'xlsxwriter') as writer:

            ####################
            try:
                #read in excel sheet with parameters
                os.chdir(exp_dir)# + '/' + 'InputExcelParameters')
                input_model_frame = pd.read_excel('InputExcelParameters.xlsx',sheet_name='Atmosphere_Model_Parameters')
                params = list(input_model_frame.iloc[:,0])
                param_values = list(input_model_frame['LAT Model Parameters'])



                #add sensitivities
                params.append('Sensitivities')
                param_values.append(sens_list)

                param_values = {'LAT Model Parameters':param_values}
                #sens_row = pd.DataFrame(data={'LAT Model Parameters':[[48.,24.,5.4,6.7,15.,36.]]},index=['Sensitivities'])
                #input_model_frame = input_model_frame.append(sens_row,ignore_index=False)
                #frames = [input_model_frame,sens_row]

                #plot_params = pd.concat(frames)
                plot_params = pd.DataFrame(param_values,index=params)
            ####################
            except NameError:
                #plotting parameter inputs (hard coded from LAT model)
                params = ['band_centers','beam_sizes','Sensitivities','f_knees','Cs','alpha_temp',
                   'survey_time','f_sky','ret_after_obs_cuts','non_uniformity_param','ell_max','ell_pivot',
                          'delta_ell','alpha_pol','NTubes','model_num']

                #Hard coded from LAT model
                default_values = {'LAT Model Parameters':[[27.,39.,93.,145.,225.,280.],[7.4,5.1,2.2,1.4,1.0,0.9],[48.,24.,5.4,6.7,15.,36.], [700.,700.,700.,700.,700.,700.],[200,7.7,1800,12000,68000,124000],[-3.5,-3.5,-3.5,-3.5,-3.5,-3.5], 5., 0.4, 0.2,0.85,1e4,[1000.,1000.,1000.,1000.,1000.,1000.],1,[-0.4,-0.4,-0.4,-0.4,-0.4,-0.4],[1,1,1,1,1,1],1]}

                #store the previous LAT model parameters in a dataframe  and write to the master sheet
                plot_params = pd.DataFrame(default_values,index=params)

            os.chdir(exp_dir)# + '/' + 'OutputExcelFiles')
            plot_params.to_excel(writer, sheet_name = 'N_ell_Plotting_Parameters')

            #Loop through all the telescope and camera directories to convert output files
            for i in range(len(telescope_names)):
                #os.mkdir(exp_dir + '/' + 'OutputExcelFiles'+ '/' + telescope_names[i])
                for j in range(len(camera_names[i])):
                    #Directory to store filter by camera
                    #os.mkdir(exp_dir + '/' + 'OutputExcelFiles' + '/' + telescope_names[i] + '/' + camera_names[i][j])
                    for dirName, subdirList, filelist in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]):
                        for k in range(len(filelist)):
                            if filelist[k] == 'output.txt':
                                path_in = exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                name_in = 'output.txt'
                                name_out = 'output' + '_' + telescope_names[i] + '_' + camera_names[i][j] + '.xls'
                                path_out = exp_dir# + '/' + 'OutputExcelFiles' + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                data = self.ConvertOuttoExcel(path_in, path_out, name_in, name_out)
                                path_out = exp_dir# + '/' + 'OutputExcelFiles'
                                self.WritetoSheet(data, name_out.rstrip('.xls'), path_out, writer)

                            elif filelist[k] == 'optical_power.txt':
                                path_in = exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                name_in = 'optical_power.txt'
                                name_out = 'optical_power' + '_' + telescope_names[i] + '_' + camera_names[i][j] + '.xls'
                                path_out = exp_dir# +  '/' + 'OutputExcelFiles' + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                data = self.ConvertOptPowtoExcel(path_in, path_out, name_in, name_out)
                                path_out = exp_dir# + '/' + 'OutputExcelFiles'
                                self.WritetoSheet(data, name_out.rstrip('.xls'), path_out, writer)

            #Loop through and convert Sensitivity output files
            for dirName, subdirList, filelist in os.walk(exp_dir):
                for k in range(len(filelist)):
                    if filelist[k] == 'sensitivity.txt':
                        path_in = exp_dir
                        path_out = exp_dir# + '/' + 'OutputExcelFiles'
                        name_in = 'sensitivity.txt'
                        name_out = 'sensitivity.xls'
                        data = self.ConvertSensitivitytoExcel(path_in, path_out, name_in, name_out, 0)
                        self.WritetoSheet(data, name_in.rstrip('.txt'), path_out, writer)

            for i in range(len(telescope_names)):
                for dirName, subdirList, filelist in os.walk(exp_dir + '/' + telescope_names[i]):
                    for k in range(len(filelist)):
                        if filelist[k] == 'sensitivity.txt':
                            path_in = exp_dir + '/' + telescope_names[i]
                            path_out = exp_dir# +  '/' + 'OutputExcelFiles' + '/' + telescope_names[i]
                            name_in = 'sensitivity.txt'
                            name_out = 'sensitivity' + '_' + telescope_names[i] + '.xls'
                            data = self.ConvertSensitivitytoExcel(path_in, path_out, name_in, name_out, 0)
                            path_out = exp_dir# + '/' + 'OutputExcelFiles'
                            self.WritetoSheet(data, name_out.rstrip('.xls'), path_out, writer)
                for j in range(len(camera_names[i])):
                    for dirName, subdirList, filelist in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]):
                        for k in range(len(filelist)):
                            if filelist[k] == 'sensitivity.txt':
                                path_in = exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                path_out = exp_dir# + '/' + 'OutputExcelFiles' + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                name_in = 'sensitivity.txt'
                                name_out = 'sensitivity' + '_' + telescope_names[i] + '_' + camera_names[i][j] + '.xls'
                                data = self.ConvertSensitivitytoExcel(path_in, path_out, name_in, name_out, 1)
                                path_out = exp_dir# + '/' + 'OutputExcelFiles'
                                self.WritetoSheet(data, name_out.rstrip('.xls'), path_out, writer)




        print('Done!')
        print('\n')
        
        return

    def ConvertOuttoExcel(self, file_path_input,file_path_output,file_name_input,file_name_output):
        import os
        import pandas as pd
        
        '''
        The module to convert Monte Carlo simulation outputs to output.txt to a spreadsheet. Not much 
        parsing or formattig options need to be taken into account for these outputs since the formatting
        is more or less self expalanatory. The output spreadsheet contains a | deliminator to separate the
        simulated frequency bands in the output files. The dataframe is returned so it can be written
        to the master spreadsheet.
        
        file_path_input: the path to the directory containing output.txt
        file_path_output: the path to the directory to save the converted spreadsheet
        file_name_input: the name of the input file to be converted (output.txt)
        file_name_output: the name the converted excel file will be under (output.xls)
        '''
        
        #Change directory to whereever the input file is
        os.chdir(file_path_input)

        #for ex: if file_name_input = 'sensitivity.txt' 
            #use comment = '#' to skip these rows
        data = pd.read_csv(file_name_input,delim_whitespace=True)

        #Save file in desired directory as desired format
        os.chdir(file_path_output)
        
        #If want to write to own excel sheet
        #data.to_excel(file_name_output)

        return data


    def ConvertOptPowtoExcel(self, file_path_input, file_path_output, file_name_input, file_name_output):
        import os
        import pandas as pd
        import decimal as Decimal
        
        '''
        This module will take an optical_power.txt input and convert it to an excel spreadsheet. The module
        returns the parsed dataframe so that it can be written to a master spreadsheet. The commented out
        code is the additional format option that is not necessarily useful. It parses the columns with the
        number +/- (number, number) format and appends the dataframe with a parsed version of the dataframe
        with each number in separate columns. 
        
        file_path_input: the path to the directory containing optical_power.txt
        file_path_output: the path to the directory to save the converted spreadsheet
        file_name_input: the name of the input file to be converted (optical_power.txt)
        file_name_output: the name the converted excel file will be under (optical_power.xls)
        '''
        
        
        #Change directory to where the input file is              
        os.chdir(file_path_input)
        f = open(file_name_input,'r')

        #Store lines of string from input file
        f_opt_pow = f.readlines()

        #For optical_power.txt the column header is the second row of text for all cameras
        col_header = f_opt_pow[1].split('|')
        col_header = col_header[1:-1]
        for i in range(0,len(col_header)):
            col_header[i] = col_header[i].strip(' ')

        #Initialize empty list of string data
        data = [0 for i in range(len(f_opt_pow))]

        #Get the number of parameters in the file for each camera
        num_params = 0
        for i in range(1,len(f_opt_pow)):
            if f_opt_pow[i][0] == '|':
                num_params += 1
            elif f_opt_pow[i][0] == '*':
                break

        num_params = num_params - 2
        #Fill the empty list with text lines deliminated by '|'
        for i in range(0,len(f_opt_pow)):
            if f_opt_pow[i][0] == '|':
                data.append((f_opt_pow[i].split('|')))

        #Remove the extraneous zeros and strings from initialization
        try:
            while True:
                data.remove(0)

        except ValueError:
            pass
        for i in range(0,len(data)):
                    try: 
                        while True:
                            data[i].remove('')
                            data[i].remove('\n')
                    except ValueError:
                        pass
        for i in range(0,len(data)):
            for j in range(0,len(data[i])):
                data[i][j] = data[i][j].strip(' ')

        #eliminate superfluous column headers
        num_cams = int(len(data)/(num_params+2))
        for i in range(num_cams):
            data.pop(i*num_params)
            data.pop(i*num_params)

        #Store filtered strings in a dataframe for easier manipulation
        d_cam = pd.DataFrame(data,columns = col_header)

        #The next couple blocks parses the strings in the format '0.00 +/- (0.00 , 0.00)'
        #Separating so the individual numbers can be used for data analysis
        '''
        new = d_cam[col_header[1]].str.split()

        d_cam['Power from Sky'] = 1
        d_cam['UpLim in Sky Power'] = 1
        d_cam['LowLim in Sky Power'] = 1
        for i in range(0,len(new)):
            d_cam.loc[i,['Power from Sky']] = Decimal.Decimal(new[i][0].strip(' '))
            d_cam.loc[i,['UpLim in Sky Power']] = Decimal.Decimal(new[i][2].strip(' ('))
            d_cam.loc[i,['LowLim in Sky Power']] = Decimal.Decimal(new[i][4].strip(') '))


        new = d_cam[col_header[2]].str.split()

        d_cam['Power to Detect'] = 1
        d_cam['UpLim in Detected Power'] = 1
        d_cam['LowLim in Detected Power'] = 1
        for i in range(0,len(new)):
            new[i][2] = new[i][2].split(',')
            d_cam.loc[i,['Power to Detect']] = Decimal.Decimal(new[i][0].strip(' '))
            d_cam.loc[i,['UpLim in Detected Power']] = Decimal.Decimal(new[i][2][0].strip(' ('))
            d_cam.loc[i,['LowLim in Detected Power']] = Decimal.Decimal(new[i][2][1].strip(') '))

        new = d_cam[col_header[3]].str.split()
        d_cam['Cumulative Eff'] = 1
        d_cam['UpLim in Cumulative Eff'] = 1
        d_cam['LowLim in Cumulative Eff'] = 1
        for i in range(0,len(new)):
            new[i][2] = new[i][2].split(',')
            d_cam.loc[i,['Cumulative Eff']] = Decimal.Decimal(new[i][0].strip(' '))
            d_cam.loc[i,['UpLim in Cumulative Eff']] = Decimal.Decimal(new[i][2][0].strip(' ('))
            d_cam.loc[i,['LowLim in Cumulative Eff']] = Decimal.Decimal(new[i][2][1].strip(') '))                                    

        reordered_columns = ['Element', 'Power from Sky', 'UpLim in Sky Power', 'LowLim in Sky Power', 'Power to Detect', 'UpLim in Detected Power', 'LowLim in Detected Power', 'Cumulative Eff', 'UpLim in Cumulative Eff', 'LowLim in Cumulative Eff']     
        d_cam = d_cam[reordered_columns] 
        '''

        
        #Write each camera to a new sheet
        d_cam_sep = [0 for i in range(num_cams)] 
        for i in range(0,num_cams):
            d_cam_sep[i] = d_cam.loc[i*num_params:(i+1)*num_params-1,d_cam.columns]

        #Navigate directory to where the converted file is to be stored    
        os.chdir(file_path_output)

        #Write converted dataframe to excel with each camera in a separate sheet.
        #If want to write to own excel file
        #with pd.ExcelWriter(file_name_output) as writer:
        #    for i in range(0,len(d_cam_sep)):
        #        d_cam_sep[i].to_excel(writer,sheet_name='Camera ' + str(i+1))

        #Close the text file reader
        f.close() 

        return d_cam_sep


    def ConvertSensitivitytoExcel(self, file_path_input, file_path_output, file_name_input, file_name_output,choose):
        import os
        import pandas as pd
        import decimal as Decimal
        
        '''
        This module will walk through the multiple formats of the different sensitivity files BoloCalc will output.
        With choose equal to 0, the code walks thru the compilation sensitivity files in the experiment directory.
        With choose being 1, the code will walk thru the individual camera sensitivity files as the deliminator
        for the text file is not homogenous and requires a little more massaging. The module returns the dataframe
        that has been created so it can be written to the master spreadsheet.
        
        file_path_input: the path to the directory containing the sensitivity file
        file_path_output: the path to where the sensitivity file should be saved to 
        file_name_input: the name of the file to be converted to excel (sensitivity.txt)
        file_name_output: the name of the excel sheet will be saved under (sensitivity.xls)
        choose: choosing which format of sensitivity file is converted to excel
        '''
        
        
        #if want to convert experiment sensitivity files
        if choose == 0:
            os.chdir(file_path_input)
            f = open(file_name_input,'r')

            f_sens = f.readlines()

            col_header = f_sens[0].split('|')
            for i in range(0,len(col_header)):
                col_header[i] = col_header[i].strip(' \n')

            data = pd.read_csv(file_name_input,sep='|',skiprows=lambda x: x % 2 !=0,names=col_header)
            data = data.drop(data.index[0:2])
            
            
            #additional format option to parse the number +/- (number,number) format and separate
            #into different columns
            '''
            new = data[col_header[2]].str.split()

            data['Array NET_CMB'] = 1
            data['UpLim in Array NET_CMB'] = 1
            data['LowLim in Array NET_CMB'] = 1
            for i in range(2,len(new)+2):
                data.loc[i,['Array NET_CMB']] = Decimal.Decimal(new[i][0].strip(' '))
                data.loc[i,['UpLim in Array NET_CMB']] = Decimal.Decimal(new[i][2].strip(' ('))
                data.loc[i,['LowLim in Array NET_CMB']] = Decimal.Decimal(new[i][4].strip(') '))


            new = data[col_header[3]].str.split()

            data['Array NET_RJ'] = 1
            data['UpLim in Array NET_RJ'] = 1
            data['LowLim in Array NET_RJ'] = 1
            for i in range(2,len(new)+2):
                data.loc[i,['Array NET_RJ']] = Decimal.Decimal(new[i][0].strip(' '))
                data.loc[i,['UpLim in Array NET_RJ']] = Decimal.Decimal(new[i][2].strip(' ('))
                data.loc[i,['LowLim in Array NET_RJ']] = Decimal.Decimal(new[i][4].strip(') '))

            new = data[col_header[4]].str.split()

            data['CMB Map Depth'] = 1
            data['UpLim in CMB Map Depth'] = 1
            data['LowLim in CMB Map Depth'] = 1
            for i in range(2,len(new)+2):
                data.loc[i,['CMB Map Depth']] = Decimal.Decimal(new[i][0].strip(' '))
                data.loc[i,['UpLim in CMB Map Depth']] = Decimal.Decimal(new[i][2].strip(' ('))
                data.loc[i,['LowLim in CMB Map Depth']] = Decimal.Decimal(new[i][4].strip(') '))

            new = data[col_header[5]].str.split()

            data['RJ Map Depth'] = 1
            data['UpLim in RJ Map Depth'] = 1
            data['LowLim in RJ Map Depth'] = 1
            for i in range(2,len(new)+2):
                data.loc[i,['RJ Map Depth']] = Decimal.Decimal(new[i][0].strip(' '))
                data.loc[i,['UpLim in RJ Map Depth']] = Decimal.Decimal(new[i][2].strip(' ('))
                data.loc[i,['LowLim in RJ Map Depth']] = Decimal.Decimal(new[i][4].strip(') '))

            #reorder columns
                reordered_header = ['Array NET_RJ', 'UpLim in Array NET_RJ', 'LowLim in Array NET_RJ', 'CMB Map Depth', 'UpLim in CMB Map Depth', 'LowLim in CMB Map Depth', 'RJ Map Depth', 'UpLim in RJ Map Depth', 'LowLim in RJ Map Depth']

            data = data[reordered_header]
            '''
            
            #change directory looking at to where to save the converted excel file
            os.chdir(file_path_output)
            
            #If want to write to own excel file
            #data.to_excel(file_name_output)
           
        #if want to convert camera sensitivity files
        if choose == 1:
            os.chdir(file_path_input)
            
            f = open(file_name_input,'r')

            f_sens = f.readlines()

            col_header = f_sens[0].split('|')
            for i in range(0,len(col_header)):
                col_header[i] = col_header[i].strip(' \n')

            data = pd.read_csv(file_name_input,sep='|',skiprows=lambda x: x % 2 !=0,names=col_header)
            data = data.drop(data.index[0:2])
            
            #shift third row values to correct columns (deliminator not homogenous)
            last_index=-1
            data.iloc[last_index]['Array NET_CMB'] = data.iloc[last_index]['Optical Power']
            data.iloc[last_index]['Optical Power'] = ''
            data.iloc[last_index]['Array NET_RJ'] = data.iloc[last_index]['Telescope Temp']
            data.iloc[last_index]['Telescope Temp'] = ''
            data.iloc[last_index]['CMB Map Depth'] = data.iloc[last_index]['Photon NEP']
            data.iloc[last_index]['Photon NEP'] = ''
            data.iloc[last_index]['RJ Map Depth'] = data.iloc[last_index]['Bolometer NEP']
            data.iloc[last_index]['Bolometer NEP'] = ''
            
            os.chdir(file_path_output)
            
            #If want to write to own excel file
            #data.to_excel(file_name_output)

        return data
    
    
    
    
    
    
    


    ##Function to ##