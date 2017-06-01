def modelFunc(R,ModVar,UseOp,PlotDetails,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,numberOfPoints,ndims,Plot_Exceptions,plot_SED=False):
    import numpy as np
    from dynamics import Dynamics

    from cosmocalc import cosmocalc
    if UseOp.runOption=='LC': import time
    from EATS_func import eats_function
    from EATS_func import BL_constants
    from EATS_func import alphanu_func
    from EATS_func import kappa_constants
    import warnings
    import os
    from radiation_modules import radiation_function
    from useful_modules import cgs_constants
    from radiation_modules import rad_var
    from radiation_modules import weights
    from radiation_modules import flux_allocation
    from radiation_modules import self_absorption
    if UseOp.runOption != 'fit':
        from matplotlib import pyplot as plt

    #Natural constants

    NatCon = cgs_constants()
    Kappas = kappa_constants(ModVar.p)

    # Radiation constants
    RadCon = BL_constants(ModVar.p)
    if UseOp.reverseShock:
        RadConRS = BL_constants(ModVar.pRS)
        Kappas_RS = kappa_constants(ModVar.pRS)

    D = cosmocalc(ModVar.z,H0=67.8,WM=0.308)['DL_cm']    #Distance to burst from observer in cm. Values gathered from http://adsabs.harvard.edu/abs/2015arXiv150201589P


    ModVar.eB  = 10**ModVar.eBlog
    ModVar.epsilone = 10**ModVar.epsiloneLog
    ModVar.epsilon = 10**ModVar.epsilonLog
    ModVar.R_ISM = 10**ModVar.R_ISM_log
    if ModVar.s == 0: 
        ModVar.A0 = 10**ModVar.nCMLog
    else: 
        ModVar.A0 = 10 ** ModVar.A0Log
    #ModVar.t0 = 10**ModVar.logt0
    ModVar.tprompt = 10**ModVar.tprompt_log
    ModVar.Gamma0 = 10**ModVar.Gamma0log
    
    if UseOp.fixedRSFSratio:
        ModVar.eB3 = np.copy(ModVar.eB)
        ModVar.epsilone3 = np.copy(ModVar.epsilone)
        ModVar.pRS = np.copy(ModVar.p)
    else:
        ModVar.eB3 = 10**ModVar.eB3log
        ModVar.epsilone3 = 10**ModVar.epsilone3Log
    
    selfAbs = BL_constants(ModVar.p)
    if UseOp.reverseShock:
        selfAbsRS = BL_constants(ModVar.pRS)


    ModVar.E0 = 10**ModVar.E0log * (1 - np.cos(ModVar.theta0)) #E0log is the isotropic energy, here it is corrected for the opening angle ModVar.theta0
    ModVar.M0 = ModVar.E0*NatCon.c**-2 / ModVar.Gamma0
    ModVar.M03 = 0.#2*pi*R[0]**3*(1-np.cos(ModVar.theta0)) * n * NatCon.mp /3              #M0/1000
    twoPi = 2*np.pi
    #ModVar.A0 = n * (NatCon.mp + NatCon.me) * 10**(profile_cutoff*s)
    mmean = (NatCon.me+NatCon.mp)/2.   #Constant
    

    
    #Get lightcurve from model. Exctract points corresponding to data. Compare (get chi**2) and return the log-likelihood

    #Find the last time that we want to calculate model for. 
    tobsEnd = np.max(tdata)# + t0)
    if UseOp.runOption=='fit': tobsRedUpper = tobsEnd
    else:
        tobsRedLower = .01
        tobsRedUpper = 2.16e12 #250 days
        tobsEnd = np.copy(tobsRedUpper)

    if (UseOp.runOption=='LC'):
        dynStartTime = time.time()




    ###############################
    ### Running dynamics module ###
    ###############################

    ### Tolerance in adaptive stepsize routine
    tol = 1e-3
    sensible_value = False
    while not sensible_value:
        if True:
        #try:
            #if UseOp.reverseShock: 
            Dyn = Dynamics(R[0],ModVar,UseOp,NatCon,tobsRedUpper*2,tol)

            sensible_value = True
            
        else:
        
        #except NameError as in_error:
            if str(in_error) in ['stepsize','bisect','gamma_min > gamma_max']:
                print in_error
                return 1e20,None,None
            else:
                raise NameError(in_error)
            """
            print tol
            print 'lowering tol'
            if tol < 1e-8:
                print '\n\n----------------------------\n\nOoops it crashed!\n\n'
                print ModVar.tprompt
                return float('-inf'),None,None

            tol /= 10
            """
                




    if (UseOp.runOption=='LC'): print "Dynamics module time use: %f"%(time.time()-dynStartTime)


    cosTheta = np.cos(Dyn.theta)


    chi2 = 0.
    timeGrid = 100   #How tight the lightcurve grid should be when runOption=='LC'
    timeGridSigma = 400



    if (UseOp.runOption == 'LC'):
        if UseOp.createMock:
            lightcurve = np.zeros(np.shape(tdata))
            tobsGrid = []
        elif plot_SED:
            lightcurve = np.zeros(np.shape(tdata))
            tobsGrid = np.copy(tdata)
        else:
            lightcurve = np.zeros([iterationLength,timeGrid])
            tobsGrid = np.zeros([iterationLength,timeGrid])

    if UseOp.runOption == 'one-sigma': #Evaluating lightcurves to plot one-sigma
        lightcurve = np.zeros([iterationLength,timeGridSigma])
        tobsGrid = np.zeros([iterationLength,timeGridSigma])


    if (UseOp.runOption=='LC'): dynTimeStart = time.time()


    #############################
    #### Spectrum generation ####
    #############################
    
    ## Setting spectrum generation constants and vectors

    distance_factor = (1+ModVar.z)/(2*D**2)
    
    #gamma_min = (p-2)/(p-1)*(1+NatCon.mp/NatCon.me*ModVar.epsilone*(Dyn.Gamma-1))
    Dyn.gamma_min[np.where(Dyn.gamma_min<1)] = 1.

    if UseOp.reverseShock:
        tobs_RS_cutoff = Dyn.tobs[Dyn.RS_elements_upper - 1]
        Rad = rad_var(Dyn , ModVar , UseOp , NatCon , RadCon , RadConRS)
        RS_in_EATS = True
    else:
        Rad = rad_var(Dyn , ModVar , UseOp , NatCon , RadCon)

    
    if UseOp.reverseShock: 
        #gamma_min_RS_out = (ModVar.pRS-2)/(ModVar.pRS-1)*(1+NatCon.mp/NatCon.me*ModVar.epsilone3*Dyn.gamma43_minus_one)
        Dyn.gamma_min_RS[np.where(Dyn.gamma_min_RS<1.)] = 1.


    #upperRimLim = Dyn.theta + ModVar.alpha
    #upperRimLim[np.where((Dyn.theta+ModVar.alpha)>np.pi/2)] = np.pi/2
    #tobsRim = (1+z) * (Dyn.tburst - Dyn.R * np.cos(upperRimLim) / c)     # Observing time at the rim for each radial point. Will use this in setting EATS grid. Note from 9/12 -13

    if UseOp.runOption=='LC' and not UseOp.createMock:
        Flux = flux_allocation(UseOp,iterationLength,Plot_Exceptions,timeGrid,freq,timeGrid)

    for nuIte in range(iterationLength):     #Loop over all input frequencies or time steps

        freqArr = np.array([freq[nuIte]])
        onePzFreq = (1+ModVar.z) * freqArr
        tobsRed = tdata[nuIte,:numberOfEmpties[nuIte]] #+ t0

        noChi2 = ((UseOp.runOption == 'LC') & (UseOp.createMock)) | (UseOp.runOption == 'one-sigma') | plot_SED  #This is true if we want to produce mock observations, and don't want to read in data
            
        if not noChi2:
            if (UseOp.runOption=='LC') and (not UseOp.createMock):  #Creating an equally spaced temporal grid to make smoother plots
                tdataLC = tobsRed
                tobsRed = np.logspace(np.log10(tobsRedLower),np.log10(tobsRedUpper),timeGrid)
                tobsGrid[nuIte] = tobsRed
            Fdata = FdataInput[nuIte,:numberOfEmpties[nuIte]]
            errorbar = errorbarInput[nuIte,:numberOfEmpties[nuIte]]
        elif UseOp.runOption == 'one-sigma': #If plotting the one-sigma range
            tobsRed = np.logspace(np.log10(tobsRedLower),np.log10(tobsRedUpper),timeGridSigma)
            tobsGrid[nuIte] = tobsRed
        #if useEATS:
        if True:
        #Equal Arrival Time Surface (EATS) integrator
            
            

            #Allocating space
            
#            if mockDim == 'T':
            
            EATSsteps = len(tobsRed)
            
            F = np.zeros(EATSsteps)
            
            if UseOp.reverseShock: PRSprim = np.zeros(EATSsteps)
            if UseOp.thermalComp: 
                thermal_component = np.zeros(EATSsteps)



            for rimI in range(EATSsteps):  #If we want a frequency resolution, len(tobsRed) = 1 . Note that len(tobsRed) must be a vector still, but with one element. Test: Run a file with only one point! 
                #Finding index for the nearest point behind seaked rim radius.
                
                #tobsRimNow = tobsRed[rimI]     #Tells the program what observer's time we want to find emission for
                #indRim = np.argmin(np.abs(tobsRimNow-tobsRim))    #Find what index the point at the rim has with the observer's time we are looking for.
                #indRim -= (tobsRim[indRim] > tobsRimNow) * (indRim != 0)          #Making sure found index is behind seaked point
                

                """
                    #Weights for the rim interpolation
                weightRim,weight1Rim,weight2Rim = np.log10(Dyn.tobs[indRim+1] / Dyn.tobs[indRim])   ,   np.log10(Dyn.tobs[indRim+1] / tobsRimNow)   ,   np.log10(tobsRimNow / Dyn.tobs[indRim])

                thetaRimPre = np.copy(Dyn.theta[indRim]) 
                lowerLimit_low = Dyn.theta[indRim] - ModVar.alpha
                upperLimit_low = Dyn.theta[indRim] + ModVar.alpha
                if upperLimit_low > np.pi/2: 
                    upperLimit_low = np.pi/2

                thetaRimPre2 = np.copy(Dyn.theta[indRim+1])
                lowerLimit_high = Dyn.theta[indRim+1] - ModVar.alpha
                upperLimit_high = Dyn.theta[indRim+1] + ModVar.alpha
                if upperLimit_high > np.pi/2:
                    upperLimit_high = np.pi/2

                if weightRim < 1e-2 or indRim == 0:
                    #thetaRim = np.copy(thetaRimPre)
                    lowerLimit = np.copy(lowerLimit_low)
                    upperLimit = np.copy(upperLimit_low)
                else:
                    #thetaRim = (thetaRimPre * weight1Rim + thetaRimPre2 * weight2Rim) / weightRim
                    lowerLimit = (lowerLimit_low * weight1Rim + lowerLimit_high * weight2Rim) / weightRim
                    upperLimit = (upperLimit_low * weight1Rim + upperLimit_high * weight2Rim) / weightRim
    
                """

                ### New approach. Instead of creating a grid that we interpolate the dynamical values on, we base the EATSurface on the dynamical points.

                ### Finding first index behind the centre of the EATS
                
                
                first_index = np.argmin(np.abs(Dyn.tobs - tobsRed[rimI])) ### Index to the point right behind the foremost point on the EATS
                if Dyn.tobs[first_index] > tobsRed[rimI]: ### If this point is in fact ahead of the time point we are on, we step back one step in order to be behind it
                    if first_index != 0:
                        first_index -= 1
                    else:
                        print 'Something wrong in EATS calculation!'

                ### Finding the index just behind the point at the very rim of the jet
                total_angle = Dyn.theta + ModVar.alpha
                total_angle[np.where(total_angle > np.pi/2)] = np.pi/2
                last_index = np.argmin(np.abs((1+ModVar.z)*(Dyn.tburst - Dyn.R*np.cos(total_angle)/NatCon.c) - tobsRed[rimI]))
                            
                if ((1+ModVar.z)*(Dyn.tburst[last_index] - Dyn.R[last_index]*np.cos(total_angle[last_index])/NatCon.c)) > tobsRed[rimI]:
                    last_index -= 1
                if last_index <= 0:
                    raise NameError("Data point is too early! Please decrease lower limit of radius R array in options.py and run again!")

                tobs_behind = (1+ModVar.z)*(Dyn.tburst[last_index] - Dyn.R[last_index]*np.cos(total_angle[last_index])/NatCon.c)
                tobs_before = (1+ModVar.z)*(Dyn.tburst[last_index+1] - Dyn.R[last_index+1]*np.cos(total_angle[last_index+1])/NatCon.c)
                        

                

                ### Now we want to create an array with obs times of all points on the EATS. Then we integrate using trapzoidal rule
                intermid_ind = np.arange(last_index+1,first_index+1) ### intermediate indeces, ranging from last_index+1 to first_index (all indeces inside the EATS)
                
      

                ### Angle Phi is defined from setting
                ### tobs = (1+z)*(tburst-R*cos(Phi)/c)
                ### This yields Phi = arccos(c/R*(tburst-tobs/(1+z)))

                Phi_factor = NatCon.c/Dyn.R[intermid_ind]*(Dyn.tburst[intermid_ind] - tobsRed[rimI]/(1+ModVar.z))

                nonzero_Phi_factor = np.sum((Phi_factor > 1) | (Phi_factor < -1))

                if nonzero_Phi_factor != 0: ### Jet has expanded to 90 degrees
                    print 'EATS crossed pi/2!'
                    intermid_ind = intermid_ind[np.where(Phi_factor < 1)]


                EATSrings = len(intermid_ind) + 2 ### Number of points in Dyn class inside EATSurface

                ### Same for RS
                if UseOp.reverseShock:
                    where_RS = np.where(Dyn.tobs[intermid_ind] <= tobs_RS_cutoff)[0] ### Finds what rings on the EATS has an RS counter part.
                    intermid_ind_RS = intermid_ind[where_RS] ### Finds what indeces has an RS counter part.


                    ### Check if we hit RS cutoff
                    try:
                        if np.max(intermid_ind_RS) == (Dyn.RS_elements_upper-1):
                            intermid_ind_RS = intermid_ind_RS[:-1]
                            where_RS = where_RS[:-1]

                        EATSrings_RS = len(intermid_ind_RS) + 2                      

                        ### innermost and edge elements of RS
                        first_index_RS = intermid_ind_RS[-1]
                        print 'first_index_RS =',first_index_RS
                    
                        last_index_RS = intermid_ind_RS[0] - 1


                    except: ### len(intermid_ind_RS) = 0
                        RS_in_EATS = False
                        print 'no RS'
                        print where_RS
                        pass



                ### Weights for interpolating the front point and the edge point
                if UseOp.reverseShock:
                    InterWeights = weights(Dyn , UseOp , Rad , ModVar , NatCon , tobsRed[rimI] , tobs_behind , tobs_before , first_index , last_index , onePzFreq , first_index_RS , last_index_RS)
                else:
                    InterWeights = weights(Dyn , UseOp , Rad , ModVar , NatCon , tobsRed[rimI] , tobs_behind , tobs_before , first_index , last_index , onePzFreq)    


                ### Setting array containing angle from LoS to EATS rings
                Phi = np.zeros(EATSrings)
                #if len(Phi) == 2: ### If there are no Dyn points inside the EATS, we have to interpolate the edges
                                       
                Phi[1:-1] = np.arccos(NatCon.c/Dyn.R[intermid_ind]*(Dyn.tburst[intermid_ind] - tobsRed[rimI]/(1+ModVar.z)))

                Phi[0] = InterWeights.Phi_edge
                #Phi[-1] = #NatCon.c/InterWeights.R_front*(InterWeights.tburst_front - tobsRed[rimI] / (1+ModVar.z))
                Phi[-1] = 0. ### Per definition




                """
                if np.sum(Phi < 0) != 0:
                    print Phi
                if len(intermid_ind) == 0:
                    print len(intermid_ind)
                    print nuIte
                    raw_input(Phi)
                """

                phiInter = np.ones(EATSrings) * 2 * np.pi
                if ModVar.alpha != 0: #Off-axis
                    partialRingsInd = np.where((Phi[1:-1] > Dyn.theta[intermid_ind] - ModVar.alpha) & (Phi[1:-1] < Dyn.theta[intermid_ind] + ModVar.alpha))   #Rings crossing the rim. Happens only when ModVar.alpha != 0
                    partialRingsInd_edge = (Phi[0] > InterWeights.theta_edge - ModVar.alpha) & (Phi[0] < InterWeights.theta_edge + ModVar.alpha)

                    offAxisFoV = (Dyn.theta[intermid_ind[partialRingsInd]]**2 - ModVar.alpha**2 - Phi[1:-1][partialRingsInd]**2) / (2*ModVar.alpha*Phi[1:-1][partialRingsInd])
                    if partialRingsInd_edge:
                        offAxisFoV_edge = (InterWeights.theta_edge**2 - ModVar.alpha**2 - Phi[0]**2) / (2*ModVar.alpha*Phi[0])
                        phiInter[0] = 2*np.pi - 2*np.arccos(offAxisFoV_edge)
                    ### phiInter at the centre is always 2*pi if theta is larger than alpha, and vice verse if alpha is larger than theta (orphan burst)
                    if InterWeights.theta_edge < ModVar.alpha: ### Orphan
                        phiInter[-1] = 0.
                    offAxisFoV[np.where(offAxisFoV<-1)] = -1.
                    phiInter[partialRingsInd] = 2*np.pi - 2*np.arccos(offAxisFoV)



                
                
                if np.isnan(Phi[-1]):
                    print len(intermid_ind)
                    """
                    print '---'
                    argmin = np.argmin(np.abs(tobsRed[rimI] - Dyn.tobs))
                    print tobsRed[rimI]
                    print Dyn.tobs[argmin]
                    print Dyn.tobs[argmin-1]
                    print Dyn.tobs[argmin+1]
                    print (1+ModVar.z)*(Dyn.tburst[argmin]-Dyn.R[argmin]*np.cos(Phi[0])/NatCon.c)
                    
                    print intermid_ind
                    print NatCon.c/R_front*(tburst_front - tobsRed[rimI]/(1+ModVar.z))
                    print 'tobsEnd =',Dyn.tobs[-1]/86400
                """


                nuPrim = np.zeros(EATSrings)
                nuPrim[1:-1] = onePzFreq * Dyn.Gamma[intermid_ind] * (1-Dyn.beta[intermid_ind] * np.cos(Phi[1:-1]))
                nuPrim[0] = InterWeights.nuPrim_edge

                nuPrim[-1] = InterWeights.nuPrim_front
                
                PprimTemp = np.zeros(EATSrings)
                PprimTemp[1:-1] = radiation_function(Dyn , Rad , UseOp , ModVar , nuPrim , Phi[1:-1] , intermid_ind , Kappas, False , True)
                PprimTemp[0] , PprimTemp[-1] = radiation_function(Dyn , Rad , UseOp , ModVar , nuPrim , Phi , intermid_ind , Kappas, False , False , InterWeights , last_index , first_index) ### Interpolating edge points

                """
                if (tobsRed[rimI] > 86400*84) and (freqArr > 1e13):
                #if np.count_nonzero(PprimTemp) != EATSrings:
                    print 'tobs =',tobsRed[rimI]/86400
                    print Phi
                    print PprimTemp
                    print freqArr
                    file_name = 'plot_dir/%s.txt'%(tobsRed[rimI]/86400)
                    #np.savetxt(file_name , [Phi , PprimTemp])
                    plt.plot(Phi,PprimTemp)
                    plt.yscale('log')
                    plt.show()
                
                if len(np.where(PprimTemp<=0)[0]) > 0:
                    plt.plot(Phi,PprimTemp)
                    plt.show()
                """
                if UseOp.opticalDepth:
                    tauFS = self_absorption(Dyn , ModVar , selfAbs , Rad , NatCon , InterWeights , nuPrim , intermid_ind , False)

                    tau_factor = np.ones(EATSrings)
                    high_tauFS = np.where(tauFS > 1e-2)
                    medium_tauFS = np.where((tauFS<=1e-2) & (tauFS > 1e-8))
                    tau_factor[high_tauFS] = (1-np.exp(-tauFS[high_tauFS]))/tauFS[high_tauFS]
                    tau_factor[medium_tauFS] = (tauFS[medium_tauFS] - tauFS[medium_tauFS]**2/2 + tauFS[medium_tauFS]**4/4 - tauFS[medium_tauFS]**6/6)/tauFS[medium_tauFS]
                    """
                    if (tobsRed[rimI] > 80*86400) and (freqArr < 1e10):
                        print tauFS
                        print high_tauFS
                        print medium_tauFS
                        print len(high_tauFS[0]) + len(medium_tauFS[0])
                        print EATSrings
                        plt.plot(Phi , tau_factor)
                        plt.yscale('log')
                        plt.show()
                    """
                    if np.count_nonzero(tau_factor) != EATSrings:
                        raw_input('hold it!')
                    ### any lower tau will give tau factor 1
                    PprimTemp *= tau_factor
                if UseOp.reverseShock and RS_in_EATS:
                    print 'len(intermid_ind) =',len(intermid_ind)
                    print 'len(intermid_ind_RS) =',len(intermid_ind_RS)
                    print len(nuPrim[where_RS])
                    PRSprimTemp = np.zeros(EATSrings_RS)
                    PRSprimTemp[1:-1] = radiation_function(Dyn , Rad , UseOp , ModVar , nuPrim[where_RS] , Phi[where_RS] , intermid_ind_RS , Kappas_RS , True , True )
                    PRSprimTemp[0] , PRSprimTemp[-1] = radiation_function(Dyn , Rad , UseOp , ModVar , nuPrim[where_RS] , Phi[where_RS] , intermid_ind_RS , Kappas_RS , True , False , InterWeights , last_index_RS , first_index_RS)

                    

                    if UseOp.opticalDepth:
                        tauRS_component = self_absorption(Dyn , ModVar , selfAbsRS , Rad , NatCon , InterWeights , nuPrim[where_RS] , intermid_ind_RS , True)
                        print tauRS_component.shape
                        print tauFS.shape
                        print np.max(where_RS)
                        print np.shape(tauFS[:where_RS[-1]+3])
                        tauRS = tauFS[:where_RS[-1]+3] + tauRS_component
                        PRSprimTemp *= (1-np.exp(-tauRS))/tauRS

                    #PRSprimTemp[where_RS] = radiation_function(Dyn , numRS , nucRS , nuPrim , beta , Phi , pRS , PRSmaxF , PRSmaxS)
                    #PRSprimTemp[where_angleInd_RS] , rho3primBehind , rho3primForward , thicknessRS_behind , thicknessRS_forward = interpolation_function(Dyn.R , Dyn.Gamma , Dyn.rho4 , Dyn.M3 , numRS , nucRS , nuPrimBehind[where_angleInd_RS] , nuPrimForward[where_angleInd_RS] , Dyn.theta , beta , angleInd_RS , cosAng[where_angleInd_RS] , pRS , PmaxF_RS , PmaxS_RS , weight[where_angleInd_RS] , weight1[where_angleInd_RS] , weight2[where_angleInd_RS])


                ### Now we 



                """
                    #Angular grid
                xAngInt = np.linspace(0,upperLimit,surfaceRings+1)    #The shell is divided into surfaceRing number of rings, with surfaceRings+1 number of borders. xAnd has the central angle of each ring segment, while xAngInt has the border angle. 
                xAng = (xAngInt[:-1] + xAngInt[1:]) / 2
                

                ### Use interpolating when calculating forward shock?
                F_interpolation = True

                if xAng[-1] > 0.5:
                    cosAngm1 = np.cos(xAng) - 1
                else:
                    cosAngm1 = -xAng**2/2 + xAng**4/24 - xAng**6/720 ### Taylor expansion
                cosAng = np.cos(xAng)
                ### Looping along the Equal Arrival Time Surface

                #
                angleInd[0] = np.argmin(np.abs((1+z)*(tburst - Dyn.R*cosAngm1[0]/c) - tobsRed[rimI])) + 1 #Index of the time at the surface point in the observer's LoS
                #if tobs_to_EATS_diff 
                

                while tobs[angleInd[0]] == tobs[angleInd[0]+1]: ### In case of duplicate entries
                    angleInd[0] += 1
                if angleInd[0] == 0:
                    PprimTemp = 0.
                    if UseOp.reverseShock:
                        PRSprimTemp = 0.
                    continue
                where_angleInd = np.ones(surfaceRings,dtype=bool)  ### If this point is before t=0

                ### Finds the elements to the EATS


                
                for indAng in range(0,surfaceRings):
                    if indAng > 0:
                        angleInd[indAng] = np.copy(angleInd[indAng-1])
                    if useEATS:
                        EATS_time_shift = (Dyn.R[angleInd[indAng]] * cosAngm1[indAng])/c
                    else:
                        EATS_time_shift = Dyn.R[angleInd[indAng]] / c
                    while tobsRed[rimI] < ((1+z) * (tburst[angleInd[indAng]] - EATS_time_shift)):
                        if angleInd[indAng] == 0:  ### If a too early point is sought
                            where_angleInd[indAng] = False
                            break
                        angleInd[indAng] -= 1  #Finding the index corresponding to the shell behind the currently evaluated point on the EATS

                    if tobsRed[rimI] > ((1+z) * (tburst[angleInd[indAng]+1] - (Dyn.R[angleInd[indAng]+1] * cosAngm1[indAng])/c)):
                        print 'tobsFront behind!!!!'
                if angleInd[indAng] < 0:
                    print 'angleInd < 0'
                    raise SystemExit(0)

                intermid_ind = angleInd[where_angleInd]
                
                ### Finds what elements on the EATsurface has a counterpart in the RS
                if UseOp.reverseShock:
                    where_angleInd_RS = np.where((intermid_ind >= RS_elements_lower) & ( intermid_ind < (RS_elements_upper-1)))
                    angleInd_RS = np.copy(intermid_ind[where_angleInd_RS])   ### These elements should be used when calculating radiation from RS

                """
                
                
                #angle_integ = (np.cos(xAngInt[:-1]) - np.cos(xAngInt[1:])) * phiInter   #Angular integration segments
                if not Plot_Exceptions.RS_only and not Plot_Exceptions.FS_only:
                    if UseOp.reverseShock and RS_in_EATS:
                        PprimTot = np.copy(PprimTemp)
                        PprimTot[:where_RS[-1]+3] += PRSprimTemp
                    else:
                        PprimTot = PprimTemp
                elif Plot_Exceptions.FS_only or not UseOp.reverseShock:
                    PprimTot = PprimTemp
                elif Plot_Exceptions.RS_only:
                    PprimTot = PRSprimTemp
                """
                if (tobsRed[rimI] > 84*86400 ) and (freqArr < 1e10):
                    print tobsRed[rimI] / 86400
                    print freqArr
                    print phiInter[0]
                    plt.plot(Phi,phiInter)
                    plt.show()
                """
                F[rimI] = np.trapz(PprimTot * phiInter ,  np.cos(Phi)) * distance_factor
                ### Negative sign is because integration is reversed on the x-axis (angle axis)


                if UseOp.runOption == 'LC' and not UseOp.createMock:
                    if not Plot_Exceptions.RS_only:
                        Flux.FFS[nuIte,rimI] = np.trapz(PprimTemp * phiInter , Phi) * distance_factor
                    if UseOp.reverseShock and not Plot_Exceptions.FS_only and RS_in_EATS:
                        Flux.FRS[nuIte,rimI] = np.trapz(PRSprimTemp * phiInter[:where_RS[-1]+3] , Phi[:where_RS[-1]+3]) * distance_factor
                    Flux.Ftotal[nuIte,rimI] = np.copy(F[rimI])

        if UseOp.plotComponents and UseOp.runOption == 'LC':
            if not Plot_Exceptions.RS_only:
                plt.plot(tobsRed/PlotDetails.scalePlotTime[UseOp.daysOrSec] , Flux.FFS[nuIte] * PlotDetails.scaleFluxAxis[UseOp.fluxAxis] , '%s--'%PlotDetails.colourCycle[nuIte])
            if UseOp.reverseShock:
                if not Plot_Exceptions.FS_only:
                    plt.plot(tobsRed/PlotDetails.scalePlotTime[UseOp.daysOrSec] , Flux.FRS[nuIte] * PlotDetails.scaleFluxAxis[UseOp.fluxAxis] , '%s:'%PlotDetails.colourCycle[nuIte])


        if UseOp.thermalComp:
            Fthermal = (1+z)/(2*D**2)*thermal_component
            if UseOp.plotComponents and (UseOp.runOption == 'LC'):
                colourCycle = ['b','g','r','c','m','y','k']     #Cycle of colours in plotting. Matplotlib standard cycle                                                
                scalePlotTime = {'d': 86400. , 'h' : 3600. , 'm' : 60. , 's' : 1. }
                scaleFluxAxis = {'mJy' : 1.e3 , 'Jy' : 1. }
                if not UseOp.reverseShock and not thermal_only: plt.plot(tobsRed/scalePlotTime[UseOp.daysOrSec],F * scaleFluxAxis[UseOp.fluxAxis],'%s--'%colourCycle[nuIte%len(colourCycle)])
                plt.plot(tobsRed/scalePlotTime[UseOp.daysOrSec],Fthermal * scaleFluxAxis[UseOp.fluxAxis],'%s:'%colourCycle[nuIte%len(colourCycle)])
            if (not Plot_Exceptions.RS_only) and (not Plot_Exceptions.FS_only):
                if thermal_only: F = Fthermal
                else:  F += Fthermal

        if noChi2 == False: ### noChi2 is true if data is not loaded
            if UseOp.runOption=='LC':
                for fInd2 in range(numberOfEmpties[nuIte]):
                    #testF[fInd2] = F[np.argmin(np.abs(tdataLC[fInd2]-tobsRed))]
                    ### Interpolating produced points onto data points
                    middleInd = np.argmin(np.abs(tdataLC[fInd2]-tobsRed))
                    behindInd = middleInd - (tdataLC[fInd2] < tobsRed[middleInd]) ### Index directly behind the model point
 
                    Fweight1 , Fweight2 , Fweight = np.log10(tdataLC[fInd2]) - np.log10(tobsRed[behindInd]) , np.log10(tobsRed[behindInd+1]) - np.log10(tdataLC[fInd2]) , np.log10(tobsRed[behindInd+1]) - np.log10(tobsRed[behindInd])
                    Finter = (F[behindInd] * Fweight2 + F[behindInd + 1] * Fweight1) / Fweight

                    if UseOp.chi2_type == 'lin': chi2 += ((Fdata[fInd2] - Finter) / errorbar[fInd2])**2
                    elif UseOp.chi2_type == 'log': chi2 += (np.log10(Fdata[fInd2]/Finter)) ** 2 / np.log10((Fdata[fInd2]+errorbar[fInd2])/Fdata[fInd2])**2
                    else: 
                        print 'Bad chi2 type %s. Now exiting'%UseOp.chi2_type
                        raise SystemExit(0)
                
            else:
                #if (np.sum(F < 0) == 0): # Avoiding negativ fluxes
                
                if UseOp.chi2_type == 'lin': 
                    F[np.where(F<0)] = 0.
                    chi2 += np.sum(((Fdata - F) / errorbar)**2)
                elif UseOp.chi2_type == 'log': 
                    F[np.where(F<=0)] = 1e-30
                    chi2 += np.sum(np.log10(Fdata/F)**2 / np.log10((Fdata+errorbar)/Fdata)**2)
                    if np.count_nonzero(F) != len(F):
                        raw_input('F has %d nonzeros but length %d'%(np.count_nonzero(F),len(F)))
                    if np.isnan(chi2) or np.isinf(chi2):
                        raw_input(F)
                else: 
                    print 'Bad chi2 type %s. Now exiting'%UseOp.chi2_type
                    raise SystemExit(0)
                if chi2 <= 0: 
                    return float('inf'),None,None

                #else:
                #    print 'Bad flux output. Returning chi2 = \'inf\''
                #    return float('inf'),float('nan'),float('nan') #If we get a negativ flux, return 'NaN'
        
        if (UseOp.runOption == 'LC') or (UseOp.runOption == 'one-sigma'): 
#            print "Loop time = %f"%loopTimerTotal
            lightcurve[nuIte] = np.copy(F)#np.concatenate([F ,   [-1]*(len(lightcurve[nuIte]) - len(F))])

    if (UseOp.runOption == 'LC') or (UseOp.runOption == 'one-sigma'): 

        print "Synchrotron time use: %f s"%(time.time() - dynTimeStart)

        ### Saving flux
        if not os.path.isdir('Flux'):
            os.mkdir('Flux/')
            print 'Created directory Flux'
        if not Plot_Exceptions.RS_only:
            if os.path.isfile('Flux/FFS.txt'):
                os.system('rm Flux/FFS.txt')
            np.savetxt('Flux/FFS.txt' , Flux.FFS)
        if UseOp.reverseShock and not Plot_Exceptions.FS_only:
            if os.path.isfile('Flux/FRS.txt'):
                os.system('rm Flux/FRS.txt')
            np.savetxt('Flux/FRS.txt' , Flux.FRS)
        if os.path.isfile('Flux/Ftotal.txt'):
            os.system('rm Flux/Ftotal.txt')
        np.savetxt('Flux/Ftotal.txt' , Flux.Ftotal)
        np.savetxt('Flux/tobs.txt' , tobsRed)
        if UseOp.allowPrint & (noChi2 == False): print "chi2 = %s\nReduced chi2 = %s"%(chi2,chi2 / (numberOfPoints-ndims) )
        #print lightcurve

        
                #outputCounter = 0
            #else: outputCounter += 1
        fovAngle = 1 / Dyn.Gamma   ### Field of view angle. Not used in EATS integration
        startJB_index = np.argmin(np.abs(Dyn.theta-ModVar.alpha-fovAngle))

        if startJB_index >= len(Dyn.tobs)-1: startJetBreak = - Dyn.tobs[startJB_index] ### If jetbreak starts at really late times
        else:
            startJB_index -= ((Dyn.theta[startJB_index] - ModVar.alpha - fovAngle[startJB_index]) < 0)  ### Making sure the field of view is still just a little bit smaller than the rim
            startJB_weight1 = Dyn.theta[startJB_index] - ModVar.alpha - fovAngle[startJB_index]
            startJB_weight2 = fovAngle[startJB_index+1] - Dyn.theta[startJB_index+1] + ModVar.alpha
            print startJB_weight1
            print startJB_weight2
        
        endJB_index = np.argmin(np.abs(Dyn.theta+ModVar.alpha-fovAngle))
        if endJB_index >= len(Dyn.tobs)-1: endJetBreak = - Dyn.tobs[endJB_index]  ### if jetbreak starts at really late times
        else:    

            endJB_index -= ((Dyn.theta[endJB_index] + ModVar.alpha) < fovAngle[endJB_index])  ### Making sure the field of view is still just a little bit smaller than the last crossing of the rim
            endJB_weight1 = Dyn.theta[endJB_index] + ModVar.alpha - fovAngle[endJB_index]
            endJB_weight2 = fovAngle[endJB_index+1] - Dyn.theta[endJB_index+1] - ModVar.alpha


        if startJB_index < len(Dyn.tobs)-1: startJetBreak = (Dyn.tobs[startJB_index] * startJB_weight2 + Dyn.tobs[startJB_index+1] * startJB_weight1) / (startJB_weight1 + startJB_weight2)
        if endJB_index < len(Dyn.tobs)-1: endJetBreak = (Dyn.tobs[endJB_index] * endJB_weight2 + Dyn.tobs[endJB_index+1] * endJB_weight1) / (endJB_weight1 + endJB_weight2)

        startJB_text = '%s %f'%('='*(startJetBreak > 0) + '>'*(startJetBreak < 0) , startJetBreak*(((startJetBreak>0)*2) - 1) / 86400)
        endJB_text = '%s %f'%('='*(endJetBreak > 0) + '>'*(endJetBreak < 0) , endJetBreak*(((endJetBreak>0)*2) - 1) / 86400)

        print "Field of view started crossing the rim at tobs %s days and covered the entire rim at tobs %s days."%(startJB_text,endJB_text)

            

        return lightcurve , startJetBreak , endJetBreak , tobsGrid
    elif UseOp.runOption == 'one-sigma': 
        return lightcurve , 0. , 0.  ,  tobsGrid
    else: 
        return chi2 , None , None
