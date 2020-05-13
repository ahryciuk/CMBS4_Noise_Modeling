#All functions related to appending files

##############################Append input files######################################

class AppendFiles():

    def __init__(self):
        
        
        return

    def printsome(self):
        
        print('Works')
        
        return
    
    ##ConfigReader##
    def ConfigReader(self, file_name_input, file_name_output, writer, file_path_output):

        import os
        import pandas as pd
  
        inputs = pd.read_csv(file_name_input, sep = '|', comment = '#')
        os.chdir(file_path_output)
        inputs.to_excel(writer, sheet_name = file_name_output)


        return
    
    ##Convert Inputs Parameters to Excel##
    def InputConvert(self, exp_dir):

        try:
            print('Converting Inputs to Excel......')

            import os
            import pandas as pd
            import subprocess, sys
            from BoloCalcConverters import TeleCamNames
            


            telescope_names, camera_names = TeleCamNames.TeleCamNames(exp_dir)

            os.mkdir(exp_dir + '/' +  'InputExcelParameters')
            file_path_output = exp_dir + '/' + 'InputExcelParameters'
            InputExcel = 'InputExcelParameters.xlsx'

            os.chdir(exp_dir + '/' + 'config')
            

            #Try with .rstrip
            with pd.ExcelWriter(InputExcel, engine = 'xlsxwriter') as writer:
                os.chdir(exp_dir + '/' + 'config')
                self.ConfigReader('foregrounds.txt', 'foregrounds', writer, file_path_output)
                for i in range(len(telescope_names)):
                    for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + 'config'):
                        for k in range(len(fileList)):
                            os.chdir(exp_dir + '/' + telescope_names[i] + '/' + 'config')
                            self.ConfigReader(fileList[k], fileList[k].rstrip('.txt') + '_' + telescope_names[i], writer, file_path_output)
                    for j in range(len(camera_names[i])):
                        for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/' + 'config'):
                            for k in range(len(fileList)):
                                os.chdir(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/' + 'config')
                                self.ConfigReader(fileList[k], fileList[k].rstrip('.txt') + '_' + telescope_names[i] + '_' + camera_names[i][j], writer, file_path_output)
                writer.save()
                print('Done! Opening file!')
                print('\n')
                
        except FileExistsError:
            print('You have already converted the input files to excel')
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
        
        #import the names of the telescopes and cameras from the given directory
        telescope_names, camera_names = TeleCamNames.TeleCamNames(exp_dir)
        
        try:
            print('Appending the Input Config Text Files')

            #InputExcelParameters Path
            input_path = exp_dir + '/' + 'InputExcelParameters'
            
            #Change the directory path to where the inputs have been appended
            os.chdir(input_path)
            file_name = 'InputExcelParameters.xlsx' 
            
            #Append the foregrounds text file 
            self.AppendConfigFiles(file_name, 'foregrounds', '/mnt/c/Users/12622/Desktop/CMB_Research/S4/Noise_Modelling/V3_baseline_2/V3_baseline_2/config', 'foregrounds.txt')

            #Loop through all the telescopes and camera directories and append each text file from the corresponding sheet
            for i in range(len(telescope_names)):
                for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + 'config'):
                    for k in range(len(fileList)):
                        output_file_name = fileList[k]
                        output_path = exp_dir + '/' + telescope_names[i] + '/' + 'config'
                        sheet_name = fileList[k].rstrip('.txt') + '_' + telescope_names[i]
                        os.chdir(input_path)
                        self.AppendConfigFiles(file_name, sheet_name, output_path, output_file_name)
                for j in range(len(camera_names[i])):
                    for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/' + 'config'):
                        for k in range(len(fileList)):
                            output_file_name = fileList[k]
                            output_path = exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/' + 'config'
                            sheet_name = fileList[k].rstrip('.txt') + '_' + telescope_names[i] + '_' + camera_names[i][j]
                            os.chdir(input_path)
                            self.AppendConfigFiles(file_name, sheet_name, output_path, output_file_name)

            print('Done!')
            print('\n')
            
        #If the user has yet to convert the inputs to excel and append them    
        except (FileExistsError, FileNotFoundError):
            print('Have you appended the inputs in excel?')
            print('\n')
            
        return

#############################Saves Files to some directory#########################################
  
    def SaveFiles(self, exp_dir, save_path, save_name, in_or_out):
        import shutil, os
        
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
        
        telescope_names, camera_names = TeleCamNames.TeleCamNames(exp_dir)
        
        try:
            print('Converting Outputs to Excel. If a warning pops up do not fear.')

            #loop thru camera files to extract output files
            OutputExcel = 'OutputExcelFiles.xlsx'
            path_in = 0
            path_out = 0
            os.mkdir(exp_dir + '/' + 'OutputExcelFiles')
            
            with pd.ExcelWriter(OutputExcel, engine = 'xlsxwriter') as writer:
                #plotting parameter inputs
                params = ['band_centers','beam_sizes','Sensitivities','f_knees','Cs','alpha_temp',
                   'survey_time','f_sky','ret_after_obs_cuts','non_uniformity_param','ell_max','ell_pivot',
                          'delta_ell','alpha_pol','NTubes_LF','NTubes_MF','NTubes_UHF']
                
                default_values = {'Plotting Parameters':[[27.,39.,93.,145.,225.,280.],[7.4,5.1,2.2,1.4,1.0,0.9],[48.,24.,5.4,6.7,15.,36.],
                                                         [700.,700.,700.,700.,700.,700.],[200,7.7,1800,12000,68000,124000],
                                                         -3.5,5.,0.4,0.2,0.85,1e4,1000.,5,0.4,1,4,2]}
                plot_params = pd.DataFrame(default_values,index=params)
                os.chdir(exp_dir + '/' + 'OutputExcelFiles')
                plot_params.to_excel(writer, sheet_name = 'N_ell_Plotting_Parameters')
                
                for i in range(len(telescope_names)):
                    os.mkdir(exp_dir + '/' + 'OutputExcelFiles'+ '/' + telescope_names[i])
                    for j in range(len(camera_names[i])):
                        #Directory to store filter by camera
                        os.mkdir(exp_dir + '/' + 'OutputExcelFiles' + '/' + telescope_names[i] + '/' + camera_names[i][j])
                        for dirName, subdirList, filelist in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]):
                            for k in range(len(filelist)):
                                if filelist[k] == 'output.txt':
                                    path_in = exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                    name_in = 'output.txt'
                                    name_out = 'output' + '_' + camera_names[i][j] + '.xls'
                                    path_out = exp_dir + '/' + 'OutputExcelFiles' + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                    data = self.ConvertOuttoExcel(path_in, path_out, name_in, name_out)
                                    path_out = exp_dir + '/' + 'OutputExcelFiles'
                                    self.WritetoSheet(data, name_out.rstrip('.xls'), path_out, writer)
                                    
                                elif filelist[k] == 'optical_power.txt':
                                    path_in = exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                    name_in = 'optical_power.txt'
                                    name_out = 'optical_power' + '_' + camera_names[i][j] + '.xls'
                                    path_out = exp_dir +  '/' + 'OutputExcelFiles' + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                    data = self.ConvertOptPowtoExcel(path_in, path_out, name_in, name_out)
                                    path_out = exp_dir + '/' + 'OutputExcelFiles'
                                    self.WritetoSheet(data, name_out.rstrip('.xls'), path_out, writer)
                                    
                #Convert Sensitivity output files
                for dirName, subdirList, filelist in os.walk(exp_dir):
                    for k in range(len(filelist)):
                        if filelist[k] == 'sensitivity.txt':
                            path_in = exp_dir
                            path_out = exp_dir + '/' + 'OutputExcelFiles'
                            name_in = 'sensitivity.txt'
                            name_out = 'sensitivity.xls'
                            data = self.ConvertSensitivitytoExcel(path_in, path_out, name_in, name_out, 0)
                            self.WritetoSheet(data, name_in.rstrip('.txt'), path_out, writer)
                            
                for i in range(len(telescope_names)):
                    for dirName, subdirList, filelist in os.walk(exp_dir + '/' + telescope_names[i]):
                        for k in range(len(filelist)):
                            if filelist[k] == 'sensitivity.txt':
                                path_in = exp_dir + '/' + telescope_names[i]
                                path_out = exp_dir +  '/' + 'OutputExcelFiles' + '/' + telescope_names[i]
                                name_in = 'sensitivity.txt'
                                name_out = 'sensitivity' + '_' + telescope_names[i] + '.xls'
                                data = self.ConvertSensitivitytoExcel(path_in, path_out, name_in, name_out, 0)
                                path_out = exp_dir + '/' + 'OutputExcelFiles'
                                self.WritetoSheet(data, name_out.rstrip('.xls'), path_out, writer)
                    for j in range(len(camera_names[i])):
                        for dirName, subdirList, filelist in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]):
                            for k in range(len(filelist)):
                                if filelist[k] == 'sensitivity.txt':
                                    path_in = exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                    path_out = exp_dir + '/' + 'OutputExcelFiles' + '/' + telescope_names[i] + '/' + camera_names[i][j]
                                    name_in = 'sensitivity.txt'
                                    name_out = 'sensitivity' + '_' + telescope_names[i] + '_' + camera_names[i][j] + '.xls'
                                    data = self.ConvertSensitivitytoExcel(path_in, path_out, name_in, name_out, 1)
                                    path_out = exp_dir + '/' + 'OutputExcelFiles'
                                    self.WritetoSheet(data, name_out.rstrip('.xls'), path_out, writer)
                                    
            
                
                
            print('Done!')
            print('\n')
        except FileExistsError:
            print('You have already converted the outputs to excel')
            print('\n')
        return

    def ConvertOuttoExcel(self, file_path_input,file_path_output,file_name_input,file_name_output):
        import os
        import pandas as pd
        
        #Change directory to whereever the input file is
        os.chdir(file_path_input)

        #for ex: if file_name_input = 'sensitivity.txt' 
            #use comment = '#' to skip these rows
        data = pd.read_csv(file_name_input,delim_whitespace=True)

        #Save file in desired directory as desired format
        os.chdir(file_path_output)
        data.to_excel(file_name_output)

        return data


    def ConvertOptPowtoExcel(self, file_path_input, file_path_output, file_name_input, file_name_output):
        import os
        import pandas as pd
        import decimal as Decimal
        
        
        #Change directory to where the input file is              
        os.chdir(file_path_input)
        f = open(file_name_input,'r')

        #Store lines of string from input file
        f_opt_pow = f.readlines()

        #For optical_power.txt the column header is the second row of text for all cameras
        col_header = f_opt_pow[1].split('|')
        col_header = col_header[1:5]
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

        #Write each camera to a new sheet
        d_cam_sep = [0 for i in range(num_cams)] 
        for i in range(0,num_cams):
            d_cam_sep[i] = d_cam.loc[i*num_params:(i+1)*num_params-1,d_cam.columns]

        #Navigate directory to where the converted file is to be stored    
        os.chdir(file_path_output)

        #Write converted dataframe to excel with each camera in a separate sheet.
        with pd.ExcelWriter(file_name_output) as writer:
            for i in range(0,len(d_cam_sep)):
                d_cam_sep[i].to_excel(writer,sheet_name='Camera ' + str(i+1))

        #Close the text file reader
        f.close() 

        return d_cam_sep


    def ConvertSensitivitytoExcel(self, file_path_input, file_path_output, file_name_input, file_name_output,choose):
        import os
        import pandas as pd
        import decimal as Decimal
        
        
        if choose == 0:
            os.chdir(file_path_input)
            f = open(file_name_input,'r')

            f_sens = f.readlines()

            col_header = f_sens[0].split('|')
            for i in range(0,len(col_header)):
                col_header[i] = col_header[i].strip(' \n')

            data = pd.read_csv(file_name_input,sep='|',skiprows=lambda x: x % 2 !=0,names=col_header)
            data = data.drop(data.index[0:2])



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

            os.chdir(file_path_output)
            data.to_excel(file_name_output)
            
        if choose == 1:
            os.chdir(file_path_input)
            
            f = open(file_name_input,'r')

            f_sens = f.readlines()

            col_header = f_sens[0].split('|')
            for i in range(0,len(col_header)):
                col_header[i] = col_header[i].strip(' \n')

            data = pd.read_csv(file_name_input,sep='|',skiprows=lambda x: x % 2 !=0,names=col_header)
            data = data.drop(data.index[0:2])
            
            os.chdir(file_path_output)
            data.to_excel(file_name_output)

        return data
    
    
    
    
    
    
    


    ##Function to ##