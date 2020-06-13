#This runs calcBolos.py



class runCalcBolos():
    import sys,os
    from sys import platform
    
    def __init__(self):
        
        
        
        return
    
    def runSim(self, exp_dir, cwd):
        import os
        from sys import platform
        
        '''
        This module executes calcBolos.py from the command line from the
        experimental directory passed by the interface. The typical calcBolos.py
        outputs will be printed to the command line rather than the Jupyter interface.
        The if statements will hopefully modulate based on the operating system that 
        the user's machine is using (since these are bash commands). 
        
        
        exp_dir: the experimental directory to be simulated
        '''
        
        print('Running Simulation')
        print('Sit back........')
        print('.................relax')
        print('.......................and enjoy the show')
        print('\n')

#############May change depending on what operating system is being used!##################
        
        if platform == 'linux' or platform == 'linux2':
            os.chdir(cwd)
            os.system('python3 calcBolos.py ' + exp_dir)
        elif platform == 'darwin':
            os.chdir(cwd)
            os.system('python3 calcBolos.py ' + exp_dir)
        elif platform == 'win32':
            os.chdir(cwd)
            os.system('python3.exe calcBolos.py ' + exp_dir)
        else:
            print('Hmm I am not familiar with this operating system')
         
        print('Done!')
        return

    def runCalcBolos_check(self, exp_dir):

        print('My records indicate you have previously ran a simulation')
        print('The output files will be overwritten')
        print('Are you sure you would like to proceed?')



        return


    def SimulationCheck(self, exp_dir):

        sim_check = False

        subdirs = [0]
        files = [0]
        for dirName, subdirList, fileList in os.walk(exp_dir):
            subdirs.append(subdirList)
            files.append(fileList)

        if 'OutputExcelFiles' in subdirs or 'sensitivity.txt' in files:
            sim_check = True
        else:
            sim_check = False


        return sim_check
    
    def SimulationProgress(self):
        
        import sys, os
        from src import simulation
        
        #sim = simulation.Simulation()
        #Find instance of simulation
        #see how to link the instance of the running simulation to this one
        #maybe stick this into the runSim function
        #export this to the interface so it can be displayed by the progress bar
        
        sim_prog = 10
        
        
        
        return sim_prog