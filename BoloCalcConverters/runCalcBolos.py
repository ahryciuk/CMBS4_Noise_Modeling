#This runs calcBolos.py



class runCalcBolos():
    import sys,os
    from sys import platform
    
    def __init__(self):
        
        
        
        return
    
    def runSim(self, exp_dir):
        import os
        
        print('Running Simulation')
        print('Sit back........')
        print('.................relax')
        print('.......................and enjoy the show')
        print('\n')

#############May change depending on what operating system is being used!##################
        
        if platform == 'linux' or platform == 'linux2':
            !python3 calcBolos.py $exp_dir
        elif platform == 'darwin':
            !python3 calcBolos.py $exp_dir
        elif platform == 'win32':
            !python3.exe calcBolos.py $exp_dir
          
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