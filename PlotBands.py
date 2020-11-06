

def PlotBand(bandfile):
    import numpy as np
    import os
    import pandas as pd
    
    bandpass = pd.read_csv(bandfile,sep='\t')

    freqs = bandpass.iloc[:,0]
    trans = bandpass.iloc[:,1]
    import matplotlib.pyplot as plt
    plt.plot(freqs,trans)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Transmission')
    plt.show()
    
    return


def PlotAllBands(exp_dir, cwd):
    import pandas as pd
    from BoloCalcConverters import TeleCamNames
    import matplotlib.pyplot as plt
    import os
    
    #exp_dir = '/mnt/c/Users/12622/Desktop/CMB_Research/S4/Noise_Modelling/BoloCalc_Simulations/V5/S4_baseline/CD/Real_Bands'


    #import the names of the telescopes and cameras from the given directory
    telescope_names, camera_names = TeleCamNames.TeleCamNames(exp_dir)

    #walk thru config directories
    plot_index = 1
    try:
        for i in range(len(telescope_names)):
            for j in range(len(camera_names[i])):
                for dirName, subdirList, fileList in os.walk(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/config/Bands/Detectors'):
                    os.chdir(exp_dir + '/' + telescope_names[i] + '/' + camera_names[i][j] + '/config/Bands/Detectors')
                    for k in range(len(fileList)):
                        data = pd.read_csv(fileList[k],sep='\t')
                        name = fileList[k].rstrip('.txt')

                        freqs = data.iloc[:,0]
                        trans = data.iloc[:,1]

                        plt.subplot(2,3,plot_index)
                        plt.plot(freqs,trans,color='g')
                        plt.title(name)
                        plt.xlabel('Frequency (GHz)')
                        plt.ylabel('Transmission')
                        plot_index += 1
                        if name == 'UHF_2':
                            break

                    else:
                        continue
                    break
                else:
                    continue
                break
            else:
                continue
            break
    except FileNotFound:
        print('Real Bandpass files not found!')

                    #bandpass[fileList[k]] = {}

                    #bandpass[fileList[k]]['freqs'] = data.iloc[:,0]
                    #bandpass[fileList[k]]['trans'] = data.iloc[:,1]

    #plt.show()
    
    os.chdir(cwd)
    return


def calc_band_center(bandfile):
    
    #integrate freq*bandpass / integrate bandpass
    
    #print expectation value
    
    return