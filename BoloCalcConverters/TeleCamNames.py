#Store telescope and camera names:
def TeleCamNames(exp_dir):
    
    import os
    
    '''
    This module will walk through the experimental directory passed by
    the interface so that all the telescope names and camera names that
    are passed can be stored in a list and outputed so other modules
    have the ability to loop through and convert/append files
    
    exp_dir: the experimental directory contianing the typical BoloCalc structure
    '''
    
    subdir = [0]
    #List subdirectories in exp_dir
    for dirName, subdirList, filelist in os.walk(exp_dir):
        subdir.append(subdirList)

    subdir.pop(0)

    #instantiate list of telescope names
    telescope_names = [0]
    
    #loop through and store telescope names
    for i in range(len(subdir[0])):
        if subdir[0][i] != 'config' and subdir[0][i] != 'paramVary' and subdir[0][i] != 'OutputExcelFiles' and subdir[0][i] != 'InputExcelParameters':
            telescope_names.append(subdir[0][i])

    #drop the index used for instantiation
    telescope_names.pop(0)

    #instantiate list of camera names
    camera_names = [[0] for i in range(len(telescope_names))]
    subdir = [[[0,0]] for i in range(len(telescope_names))]
    
    #loop through and store camera names
    for i in range(0,len(telescope_names)):
        for dirName, subdirList, filelist in os.walk(exp_dir + '/' + telescope_names[i]):
            subdir[i].append(subdirList)

    for i in range(len(telescope_names)):
        for j in range(len(subdir[i][1])):
            if subdir[i][1][j] != 'config':
                camera_names[i].append(subdir[i][1][j])

    #get rid of the extraneous initialization zeros
    for i in range(len(camera_names)):
        camera_names[i].pop(0)
        
    return telescope_names, camera_names