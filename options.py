import numpy as np
import os
c = 2.9979e10

class model_variables:
    def __init__(self):

        self.epsilonLog = -1. #Radiative efficiency, logarithmic. Only applied when fixed_epsilon == True and radiativeLosses == True
        self.epsiloneLog = np.log10(0.2)#-0.33#np.log10(2.5e-1) #Fraction of energy in the electrons
        self.epsilonpLog = -1 #Above, for protons
        self.epsilone3Log = -3
        self.epsilonp3Log = -2
        self.E0log = np.log10(1e55)    #Initial isotropic energy
        self.nCMLog = -0 # Circumstellar density, constant medium
        self.A0Log = 33 # Circumstellar density, windlike medium
        self.s = 0.0 #CBM power law slope
        self.R_ISM_log = 18.  # Radius of transition to interstellar matter. Logarithmic
        self.Gamma0log = np.log10(10.) #Initial gamma
        self.eBlog = -3#np.log10(3.5e-4) #Magnetic microphysical constant in Forward Shock (region 2)
        self.eB3log = -4   #Magnetic microphysical constant in Reverse Shock (region 3)
        self.p = 2.1  #Synchrotron emission power law
        self.logt0 = -1. #Time between R=0 and prompt time   (logarithmic)
        self.theta0 = np.pi / 180. * 5       #Opening angle
        self.alpha =  np.pi / 180. * 0          #Off-axis angle
        self.z = 1.62         #Redshift
        self.tprompt_log = 3.         #Duration time of prompt phase. Needed for the reverse shock. Input time is observer frame, the module returns the progenitor frame time. Note from 29/1 -14
        self.pRS = 2.05          #Ratio between fixed reverse shock and forward shock microphysical ceofficients. Assumes coeffs of forward shock to set coeffs of reverse shock. If you want to use this, set fixedRSFSratio = TRUE in beginning of mainFit.py.
        self.const_names = []  ### Empty array. Will be filled by calling constant_names() method

        ### Saving attribute names of self class
        for attr , value in self.__dict__.items():
            if not attr == 'const_names':
                self.const_names.append(attr)


        #return epsilonLog,epsiloneLog,epsilonpLog,epsilone3Log,epsilonp3Log,E0,nCMLog,A0Log,s,R_ISM,Gamma0log,eB,eB3,p,logt0,theta0,alpha,tN,logrphoto,Tprim0,N0,tprompt_log,pRS,R0_Clog,N_Clog,t_outflow_log,Ctheta,GammaClog,z,WV

    ### Method to change values of objects in ModVar class
    def new_value(self, index , new_value):
        
        try:
        #if True:
            for i in range(len(index)):### If index is not an array, this will fail and go into except
                exec('self.%s = new_value[i]'%self.const_names[i])
        except: ### input index is not an array
            raw_input()
            #exec('self.%s = new_value'%self.const_names[i])

    def echo_value(self, index):
        exec('out_value = self.%s'%self.const_names[index])
        return out_value


class userOptions:
    
########################################
#             User options             #
########################################

    def __init__(self):
#Input data
        self.inputData = '%s/Data_files/'%(os.path.expanduser('~'))  #Where the input data-files are stored. Store the files in a folder under said folder e.g. named 990510fit for GRB 990510. Then specify the burst date to GRBlabel
        self.fileExtention = 'dat'                                         #File extentions for input data files
        self.breakName = 'lc.dat'                                            #Type what comes after the frequency in the file name
        self.inputTimeFormat = 'auto'                                     #What format is input time in? 'day', 'hour', 'min' or 'sec'; 'auto' discriminates between formats seconds and days
        self.inputFreqUnit = 'Hz'                                        #In what unit the frequency is stated in the filename. 'cm','Hz','GHz','MHz'
        self.GRBlabel = '990510'                                           #Name of the GRB. If you want manual input, set GRBinput to True

        #Options
        self.allowDataRemoval = False       #True: allow the programme to ask if the user wants to remove previously saved data
    
        self.thermalComp = False             #Include a thermal component?
        self.photosphere = False             #Use the photosphere as the origin of the photons? If False, the rim of the cocoon will be used
        self.radiativeLosses = False          # Radiative losses in dynamics calculations?

        self.useEATS = True                 #Equal arrival time surface integrator on or off? True: on
        #surfaceRings = 15              #If not using EATS, how many rings should the shock front surface be integrated over?
        self.plotOutput = False              #If running the fitter, and you want to plot the output. Plot: True

        self.gridStart,self.gridEnd,self.gridStep = 12,24,2000          #Start of grid (log10), end of grid (log10), and number of grid points
        self.printProcess = True           #Print the values of the variable values as the fitter is running. Prints when fitter hits a new lowest chi^2
        self.allowPrint = True             #Allow the programme to return prints? For using in cluster, False is recommended to avoid obscene logs. Crucial messages will be printed if printCrucial == True, no matter this value
        self.printCrucial = True            #Print crucial messages

        #Plot options
        self.daysOrSec = 'd'                #Plot output in days or seconds? 'd': days,   'h': hours,   'm': minutes   's': seconds.
        self.fluxAxis = 'Jy'               #Choose flux units for the flux out-put. 'mJy' = milli-Janskys   ;   'mJy' = Janskys
        #                      e_rad e_e2  e_p2  e_e3  e_p3  E0    n_CM   A0    s   R_ISM Gam0  e_B2  e_B3   p2         theta0 alpha                            tprompt p3  
        self.preferredPlotScale = ['lin','log','lin','log','lin','log','log','lin','lin','log','log','log','log','lin','None','deg','deg','None','None','None','None','lin','lin','log','log','log','deg','log','lin','lin'] #Determines what x-scale the probability plots should have. Parameter preferredScale in mainFit.py determines scale used in fitter

        #Dynamics options
        self.reverseShock = True            #Include a reverse shock component?
        self.exponential_outflow = True        #Ejecta density distribution exponentially (True) or constant (False)?
        self.opticalDepth = True           #Take synchrotron self-absorption optical depth into account? True: yes
        self.fixedRSFSratio = False          #Fixed ratio between the microphysical coefficients esilone,epsilonp,epsilonB and p of the reverse shock and the forward shock?
        self.fixed_epsilon = True           #Use a fix radiative loss parameter epsilon_rad (see Johansson et al. 2016)?

        #Sanity check options
        self.printStats = True               #Print the characteristics of the jet evolution

        #Fit options
        self.runOption = 'fit'               #Run option. Fitting routine: 'fit'; Plot 3D likelihood surface of two variables: 'surf'; Print lightcurves of values from constants.txt: 'LC'
        self.save_params = True             #Print parameters to file?
        self.chi2_type = 'log'              #Evaluate chi2 logarithmically ('log') or lineraly ('lin')?
        self.livePoints = 1000               #Number of live points

        #Mock observations options


        self.mockDistribution = 'log'           #Distribution of the mock data points in the temporal grid. For manual time points, enter a float numpy array, e.g. np.array([4.3,7.2e3,4.1e5]). Options: 'log': logarithmic; 'lin': linear
        self.numberOfMock = 50                 #Number of mock data points. This option is only activated if mockDistribution is set to an automatic option
        self.mockIntervalTime = [3600,20*86400]##       #Interval of the temporal grid in mock distribution ([lowerValue,upperValue])
        self.mockIntervalFreq = [1.e10,1.e20]  #See above, for frequency (used when createMock == True and mockDim == 'E' in the runOption == 'LC' mode)
        self.createMock = False                 #Creates mock observations.
        self.gaussianError = 0.000001               #Gaussian error (1-sigma) in mock distribution
        self.gaussianOffset = True             #Activates a gaussian offset. This is the only option (for now)
        self.offsetType = 'log'                #Offset type of mock data. lin - linear offset, log - logarithmic offset
        self.useData = True                    #True: Uses observations for LC production. False: Uses input data
        self.plotComponents = True              #Plot thermal, reverse shock and forward shock components? Yes, plot them: True; No, only plot total lightcurve: False
        self.plotMock = False                   #Plots mock observation. If creating mock observation, this is automatically over-run to TRUE. Recommendation is to keep this FALSE at all times.
        self.mockDim = 'T'                      #The dimensions to create mock observations in. 'T': time; 'E': energy
        #freqGrid = [5.e14]
        #    freqGrid = [5.e14, 2.e18]
        #freqGrid = [4.56e14, 2.4e17]       #R-band and 1 keV
        #freqGrid = [2.4e18]
        #freqGrid = [1.4e9, 4.8e9, 8.4e9, 22.5e9, 4.3e10, 1.e11, 3.e11] #Radio
        #freqGrid = [2.4e18]  #x-rays
        #freqGrid = [ 6.7e14, 8.2e14, 5.4e14, 4.6e14, 3.7e14, 2.5e14, 1.4e14, 4.8e14] #Optical
        #freqGrid = [6.7e14, 8.2e14, 5.4e14, 4.6e14, 3.7e14, 2.5e14, 1.4e14, 4.8e14,1.e16,2.4e18] #xrays-UV-optical
        #freqGrid = [1.e16, 2.4e18] #UV xrays
        #freqGrid = np.concatenate([np.linspace(8.2e14,8.7e14,10) , np.linspace(2.4e17,2.4e18,10)])
        #freqGrid = [1.4e9, 4.8e9, 8.4e9, 22.5e9, 4.3e10, 1.e11, 3.e11, 6.7e14, 8.2e14, 5.4e14, 4.6e14, 3.7e14, 2.5e14, 1.4e14, 4.8e14, 2.4e18]   # Will be over-run if (runOption == 'LC') and (createMock == True)
        self.freqGrid = np.logspace(9,18,50)
        #freqGrid = [1.4e9, 2.4e18]   # Will be over-run if (runOption == 'LC') and (createMock == True)
        self.tgrid = [1.3e1,1.e8]        #Only used to create mock observations resolved in energies
    
        #Parameter constraint options

    
    
        #Output data
        self.figTypes = 'eps'                                              #File extentions for saving figures
    

        #Constants

        ############################################
        #                Intervals                 #
        ############################################

        #outputCounter = 0

        self.paramLimits = np.array([
            [-.8,0.]              #epsilon
            ,[-6.,-.1]              #epsilon e (logarithmic)
            ,[-8.,0.]              #epsilon p (logarithmic)
            ,[-6.,-.1]              #epsilon e RS (logarithmic)
            ,[-8.,0.]              #epsilon p RS (logarithmic)
            ,[48.,58.]            #log10( E0 ) (isotropic)
            ,[-4.,3.]             #Constant medium number density n (logarithmic)
            ,[32.,43.]              #Windlike medium constant A0
            ,[0.,2.1]            #s
            ,[16,21]         #R_ISM (logarithmic)
            ,[1.,4.]         #Gamma0 (logarithmic)
            ,[-8.,-1.]           #epsilon B (FS, logarithmic)
            ,[-8.,-1.]           #epsilon B (RS, logarithmic) 
            ,[2.,2.8]            #Electron energy distribution slope p for the forward shock
            ,[-5.,5.]             #t0 (logarithmic)
            ,[.0001,90*np.pi/180.]  #theta0 (half opening angle)
            ,[0.,90*np.pi/180]            #alpha (off-axis angle)
            ,[-2,7]            #Duration of outflow after trigger time (progenitor frame). Only affects the reverse shock scenario. Logarithmic
            ,[2.,3.]        #Electron energy distribution slope p for reverse shock
            ,[1e-2,10]     #Redshift
        ])
        ##########################################
        #        Variables for fitting           #
        ##########################################

                  #Order:       epsil,eps_E,eps_P, eE_RS,eP_RS, E0  n/A0   s    R_ISM Gamma;epsB; eB_RS   p    t0  theta0 alpha prompt, pRS   z 
        if self.reverseShock:
            self.parametrar = [False,True ,False,True ,False,True ,True ,False,True ,True ,True ,True ,True ,False,True ,True ,True ,True ,False] #False = Constant. 

        else:           
            self.parametrar = [False,True ,False,False,False,True ,True ,False,True ,True ,True ,False,True ,False,True ,True ,False,False,False] #False = Constant. 


