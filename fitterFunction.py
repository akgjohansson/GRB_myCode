def modelFunc(R,ModVar,UseOp,PlotDetails,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,numberOfPoints,ndims,Plot_Exceptions,plot_SED=False):
    import numpy as np
    from dynamics import Dynamics
    from matplotlib import pyplot as plt
    from cosmocalc import cosmocalc
    if UseOp.runOption=='LC': import time
    from EATS_func import eats_function
    from EATS_func import BL_constants
    from EATS_func import alphanu_func
    from EATS_func import kappa_constants
    import warnings
    from radiation_modules import radiation_function
    from useful_modules import cgs_constants
    from radiation_modules import rad_var
    from radiation_modules import weights
    from radiation_modules import flux_allocation
    from radiation_modules import self_absorption
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
        tobsRedLower = 1.
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
    else:
        Rad = rad_var(Dyn , ModVar , UseOp , NatCon , RadCon)

    
    if UseOp.reverseShock: 
        #gamma_min_RS_out = (ModVar.pRS-2)/(ModVar.pRS-1)*(1+NatCon.mp/NatCon.me*ModVar.epsilone3*Dyn.gamma43_minus_one)
        Dyn.gamma_min_RS[np.where(Dyn.gamma_min_RS<1.)] = 1.


    #upperRimLim = Dyn.theta + ModVar.alpha
    #upperRimLim[np.where((Dyn.theta+ModVar.alpha)>np.pi/2)] = np.pi/2
    #tobsRim = (1+z) * (Dyn.tburst - Dyn.R * np.cos(upperRimLim) / c)     # Observing time at the rim for each radial point. Will use this in setting EATS grid. Note from 9/12 -13

    if UseOp.runOption=='LC' and not UseOp.createMock:
        Flux = flux_allocation(UseOp,iterationLength,Plot_Exceptions,timeGrid,freq)

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
                    print 'last_index is smaller or equal to 0!!!'

                tobs_behind = (1+ModVar.z)*Dyn.tburst[last_index] - Dyn.R[last_index]*np.cos(total_angle[last_index])/NatCon.c
                tobs_before = (1+ModVar.z)*Dyn.tburst[last_index+1] - Dyn.R[last_index+1]*np.cos(total_angle[last_index+1])/NatCon.c
                        

                

                ### Now we want to create an array with obs times of all points on the EATS. Then we integrate using trapzoidal rule
                intermid_ind = np.arange(last_index+1,first_index+1) ### intermediate indeces, ranging from last_index+1 to first_index (all indeces inside the EATS)
                
                ### Weights for interpolating the front point and the edge point
                InterWeights = weights(Dyn , UseOp , Rad , ModVar , NatCon , tobsRed[rimI] , tobs_behind , tobs_before ,  , first_index , last_index , onePzFreq , Phi)
                                       

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
                    where_RS = np.where(Dyn.tobs[intermid_ind] <= tobs_RS_cutoff) ### Finds what rings on the EATS has an RS counter part.
                    intermid_ind_RS = intermid_ind[where_RS] ### Finds what indeces has an RS counter part.
                    EATSrings_RS = len(intermid_ind_RS) + 2

                ### Setting array containing angle from LoS to EATS rings
                Phi = np.zeros(EATSrings)
                #if len(Phi) == 2: ### If there are no Dyn points inside the EATS, we have to interpolate the edges
                                       
                Phi[1:-1] = np.arccos(NatCon.c/Dyn.R[intermid_ind]*(Dyn.tburst[intermid_ind] - tobsRed[rimI]/(1+ModVar.z)))

                Phi[0] = InterWeights.Phi_edge
                Phi[-1] = NatCon.c/InterWeights.R_front*(InterWeights.tburst_front - tobsRed[rimI] / (1+ModVar.z))
                #Phi[-1] = 0.

[0])



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
                nuPrim[0] = onePzFreq * InterWeights.Gamma_edge * (1-InterWeights.beta_edge * np.cos(Phi[0]))

                nuPrim[-1] = onePzFreq * InterWeights.Gamma_front * (1-InterWeights.beta_front * np.cos(Phi[-1]))

                PprimTemp = np.zeros(EATSrings)
                PprimTemp[1:-1] = radiation_function(Dyn , Rad , UseOp , ModVar , nuPrim , Phi[1:-1] , intermid_ind , Kappas, False , True)
                PprimTemp[0] , PprimTemp[-1] = radiation_function(Dyn , Rad , UseOp , ModVar , nuPrim , Phi , intermid_ind , Kappas, False , False , InterWeights , last_index , first_index)

                if UseOp.opticalDepth:
                    tauFS = np.zeros(EATSrings)
                    tauFS[1:-1] = self_absorption(Dyn , ModVar , selfAbs , Rad , NatCon , InterWeights , nuPrim , intermid_ind , False)
                    
                    
                    PprimTemp *= (1-np.exp(-tauFS))/tauFS

                if UseOp.reverseShock:

                    where_angleInd_RS = np.where((intermid_ind >= RS_elements_lower) & ( intermid_ind < (RS_elements_upper-1)))
                    angleInd_RS = np.copy(intermid_ind[where_angleInd_RS])   ### These elements should be used when calculating radiation from RS

                    PRSprimTemp = np.zeros(EATSrings)
                    PRSprimTemp[where_RS] = radiation_function(Dyn , Rad , UseOp , ModVar , nuPrim[where_RS] , Phi[where_RS] , intermid_ind_RS , Kappas_RS , True)

                    if opticalDepth:
                        tauRS = np.zeros(EATSrings_RS)
                        tauRS[1:-1] = tauFS + self_absorption(Dyn , ModVar , selfAbsRS , Rad , NatCon , InterWeights , nuPrim[where_RS] , intermid_ind_RS , True)
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
                """
                ### Finds what elements on the EATsurface has a counterpart in the RS
                if UseOp.reverseShock:
                    where_angleInd_RS = np.where((intermid_ind >= RS_elements_lower) & ( intermid_ind < (RS_elements_upper-1)))
                    angleInd_RS = np.copy(intermid_ind[where_angleInd_RS])   ### These elements should be used when calculating radiation from RS

                
                

                """
                if useEATS:
                    tobsBehind = (1+z) * (Dyn.tburst[intermid_ind] - (Dyn.R[intermid_ind] * cosAngm1 + Dyn.R[intermid_ind]) / c)
                else:
                    tobsBehind = (1+z) * (Dyn.tburst[intermid_ind] - Dyn.R[intermid_ind] / c)
                


                phiInter = (theta[angleInd]+ModVar.alpha <= lowerLimit) * 2*np.pi
                if theta[first_index] + ModVar.alpha <= lowerLimit
                
                if ModVar.alpha != 0: #Off-axis
                    partialRingsInd = np.where((xAng > Dyn.theta[intermid_ind] - ModVar.alpha) & (xAng < Dyn.theta[intermid_ind] + ModVar.alpha))   #Rings crossing the rim. Happens only when ModVar.alpha != 0
                    offAxisFoV = (Dyn.theta[angleInd[partialRingsInd]]**2 - ModVar.alpha**2 - xAng[partialRingsInd]**2) / (2*ModVar.alpha*xAng[partialRingsInd])
                    offAxisFoV[np.where(offAxisFoV<-1)] = -1.
                    if np.sum(offAxisFoV > 1) != 0: offAxisFoV[np.where(offAxisFoV > 1)] = 1.
                    phiInter[partialRingsInd] = 2*np.pi - 2*np.arccos(offAxisFoV)


                #oneDzFreq = freqArr / (1+z)
                cosAng2 = cosAng**2
                if F_interpolation:
                    nuPrimForward = onePzFreq * Dyn.Gamma[intermid_ind+1] * (1-beta[intermid_ind+1] * cosAng)

                    nuPrimBehind = onePzFreq * Dyn.Gamma[intermid_ind] * (1-beta[intermid_ind] * cosAng)


                    if useEATS:
                        tobsFront = (1+z) * (Dyn.tburst[intermid_ind+1] - (Dyn.R[intermid_ind+1] * cosAngm1 + Dyn.R[intermid_ind+1]) / c)
                    else:
                        tobsFront = (1+z) * (Dyn.tburst[intermid_ind+1] - Dyn.R[intermid_ind+1] / c)

                    weight,weight1,weight2 = np.log10( tobsFront / tobsBehind )  ,  np.log10( tobsFront / tobsRimNow )  ,  np.log10( tobsRimNow / tobsBehind )
                    

                    if np.count_nonzero(weight) != len(weight): ### Grid tight enough that tobsFront and tobsBehind are the same. We go about this by choosing the weight to be totally on the first point
                        where_weight = np.where(weight == 0)
                        weight[where_weight] = 1.
                        weight1[where_weight] = 1.
                        weight2[where_weight] = 0.
                    #PprimTemp = (Pbehind ** (1-fweight) * Pforward ** (fweight))
                    #PprimTemp = (Pbehind * weight1 + Pforward * weight2)  /  weight
                    
                    PprimTemp , rhoPrimBehind , rhoPrimForward , thickBehind , thickForward = radiation_function(Dyn , num , nuc , nuPrim , beta , Phi , p , PmaxF , PmaxS)




                if UseOp.reverseShock: 

                    PRSprimTemp = np.zeros(surfaceRings)

                    PRSprimTemp[where_angleInd_RS] , rho3primBehind , rho3primForward , thicknessRS_behind , thicknessRS_forward = interpolation_function(Dyn.R , Dyn.Gamma , Dyn.rho4 , Dyn.M3 , numRS , nucRS , nuPrimBehind[where_angleInd_RS] , nuPrimForward[where_angleInd_RS] , Dyn.theta , beta , angleInd_RS , cosAng[where_angleInd_RS] , pRS , PmaxF_RS , PmaxS_RS , weight[where_angleInd_RS] , weight1[where_angleInd_RS] , weight2[where_angleInd_RS])

                """
                """
                if UseOp.opticalDepth:

                    tau = self_absorption(num,nuc,alpha0F,alpha0S,nuPrim)


                    numCentral = (num[intermid_ind] * weight1 + num[intermid_ind+1] * weight2) / weight
                    nucCentral = (nuc[intermid_ind] * weight1 + nuc[intermid_ind+1] * weight2) / weight
                                        
                    where_slow_cooling = np.where(numCentral <= nucCentral)
                    where_fast_cooling = np.where(numCentral > nucCentral)
                    
                    alpha0Fforward = alpha0Ffactor * rhoPrimForward[where_fast_cooling] / Dyn.B[intermid_ind+1][where_fast_cooling]  *  Dyn.gammac[intermid_ind+1][where_fast_cooling]**(-5)
                    alpha0Fbehind = alpha0Ffactor * rhoPrimBehind[where_fast_cooling] / Dyn.B[intermid_ind][where_fast_cooling]  *  Dyn.gammac[intermid_ind][where_fast_cooling]**(-5)

                    alpha0Sforward =  alpha0Sfactor * rhoPrimForward[where_slow_cooling]  * Dyn.gamma_min[intermid_ind+1][where_slow_cooling]**(-5) /  Dyn.B[intermid_ind+1][where_slow_cooling]
                    alpha0Sbehind =  alpha0Sfactor * rhoPrimBehind[where_slow_cooling]  * Dyn.gamma_min[intermid_ind][where_slow_cooling]**(-5) /  Dyn.B[intermid_ind][where_slow_cooling]

                    nuPrimCentral = (nuPrimBehind * weight1 + nuPrimForward * weight2) / weight
                    alpha0Fcentral = (alpha0Fbehind * weight1[where_fast_cooling] + alpha0Fforward * weight2[where_fast_cooling]) / weight[where_fast_cooling]
                    alpha0Scentral = (alpha0Sbehind * weight1[where_slow_cooling] + alpha0Sforward * weight2[where_slow_cooling]) / weight[where_slow_cooling]
                    

                    alphanu = np.zeros(len(numCentral))
                    ### Fast cooling
                    alphanu[where_fast_cooling] = alpha0Fcentral * alphanu_func(nuPrimCentral[where_fast_cooling] , numCentral[where_fast_cooling] , nucCentral[where_fast_cooling], True , p)
                    #alphanu[where_fast_cooling] = alpha0Fcentral * ( (nuPrimCentral[where_fast_cooling]<=nucCentral[where_fast_cooling]) * (nuPrimCentral[where_fast_cooling]/nucCentral[where_fast_cooling])**(-5/3.) + ((nucCentral[where_fast_cooling]<nuPrimCentral[where_fast_cooling])&(nuPrimCentral[where_fast_cooling]<numCentral[where_fast_cooling])) * (nuPrimCentral[where_fast_cooling]/nucCentral[where_fast_cooling])**(-3) + (numCentral[where_fast_cooling]<=nuPrimCentral[where_fast_cooling]) * (numCentral[where_fast_cooling]/nucCentral[where_fast_cooling])**-3*(nuPrimCentral[where_fast_cooling]/numCentral[where_fast_cooling])**(-(p+5)/2.) )
                    ### Slow cooling

                    alphanu[where_slow_cooling] = alpha0Sforward * alphanu_func(nuPrimCentral[where_slow_cooling] , numCentral[where_slow_cooling] , nucCentral[where_slow_cooling] , False , p)
                    #alphanu[where_slow_cooling] = alpha0Sforward * ( (nuPrimCentral[where_slow_cooling]<=numCentral[where_slow_cooling])*(nuPrimCentral[where_slow_cooling]/numCentral[where_slow_cooling])**(-5/3.) + ((numCentral[where_slow_cooling]<nuPrimCentral[where_slow_cooling])&(nuPrimCentral[where_slow_cooling]<nucCentral[where_slow_cooling])) * (nuPrimCentral[where_slow_cooling]/numCentral[where_slow_cooling])**(-(p+4)/2.) + (nucCentral[where_slow_cooling]<=nuPrimCentral[where_slow_cooling]) * (nucCentral[where_slow_cooling]/numCentral[where_slow_cooling])**(-(p+4)/2.) * (nuPrimCentral[where_slow_cooling]/nucCentral[where_slow_cooling])**(-(p+5)/2.) )
                    tau = alphanu * (thickBehind * weight1 + thickForward * weight2) / 2 / weight     #tau is always negative

                    tau[np.where(tau<0)] = 0.

                        
                    where_low_tau = np.where((tau > 1e-10) & (tau<1e-2))
                    where_high_tau = np.where(tau > 1e-2)
                    where_no_tau = np.where(tau < 1e-10)

                    tau_factor = np.zeros(len(tau))
                    tau_factor[where_low_tau] = (tau[where_low_tau] - tau[where_low_tau]**2/2 + tau[where_low_tau]**3/6) / tau[where_low_tau]
                    tau_factor[where_high_tau] = (1-np.exp(-tau[where_high_tau]))/tau[where_high_tau]
                    tau_factor[where_no_tau] = 1.

                    PprimTemp *= tau_factor

                    
                    if UseOp.reverseShock:

                        

                        numRSCentral = (numRS[intermid_ind_RS] * weight1[where_intermid_ind_RS] + numRS[intermid_ind_RS] * weight2[where_intermid_ind_RS]) / weight[where_intermid_ind_RS]
                        nucRSCentral = (nucRS[intermid_ind_RS] * weight1[where_intermid_ind_RS] + nucRS[intermid_ind_RS] * weight2[where_intermid_ind_RS]) / weight[where_intermid_ind_RS]

                        alphanuRS = np.zeros(len(intermid_ind_RS))

                        ### Fast cooling
                        where_fast_cooling_RS = np.where(numRSCentral>nucRSCentral)
                        alpha0FRSforward = alpha0FRSfactor * rho3primForward[where_fast_cooling_RS] / BRS[intermid_ind_RS+1][where_fast_cooling_RS]  *  gammacRS[intermid_ind_RS+1][where_fast_cooling_RS]**(-5)
                        alpha0FRSbehind = alpha0FRSfactor * rho3primBehind[where_fast_cooling_RS] / BRS[intermid_ind_RS][where_fast_cooling_RS]  *  gammacRS[intermid_ind_RS][where_fast_cooling_RS]**(-5)

                        alpha0FRScentral = (alpha0FRSbehind * weight1[where_intermid_ind_RS][where_fast_cooling_RS] + alpha0FRSforward * weight2[where_intermid_ind_RS][where_fast_cooling_RS]) / weight[where_intermid_ind_RS][where_fast_cooling_RS]

                        alphanuRS[where_fast_cooling_RS] = alpha0FRScentral * ( (nuPrimCentral[where_intermid_ind_RS][where_fast_cooling_RS]<=nucRSCentral[where_fast_cooling_RS]) * (nuPrimCentral[where_intermid_ind_RS][where_fast_cooling_RS]/nucRSCentral[where_fast_cooling_RS])**(-5/3.) + ((nucRSCentral[where_fast_cooling_RS]<nuPrimCentral[where_intermid_ind_RS][where_fast_cooling_RS])&(nuPrimCentral[where_intermid_ind_RS][where_fast_cooling_RS]<numRSCentral[where_fast_cooling_RS])) * (nuPrimCentral[where_intermid_ind_RS][where_fast_cooling_RS]/nucRSCentral[where_fast_cooling_RS])**(-3) + (numRSCentral[where_fast_cooling_RS]<=nuPrimCentral[where_intermid_ind_RS][where_fast_cooling_RS]) * (numRSCentral[where_fast_cooling_RS]/nucRSCentral[where_fast_cooling_RS])**-3*(nuPrimCentral[where_intermid_ind_RS][where_fast_cooling_RS]/numRSCentral[where_fast_cooling_RS])**(-(ModVar.pRS+5)/2.) )

                        ### Slow cooling
                        where_slow_cooling_RS = np.where(numRSCentral>=nucRSCentral)
                        
                        alpha0SRSforward =  alpha0SRSfactor * rho3primForward[where_slow_cooling_RS]  * gamma_min_RS[intermid_ind_RS+1][where_slow_cooling_RS]**(-5) /  BRS[intermid_ind_RS+1][where_slow_cooling_RS]
                        alpha0SRSbehind =  alpha0SRSfactor * rho3primBehind[where_slow_cooling_RS]  * gamma_min_RS[intermid_ind_RS][where_slow_cooling_RS]**(-5) /  BRS[intermid_ind_RS][where_slow_cooling_RS]
                        alpha0SRScentral = (alpha0SRSbehind * weight1[where_intermid_ind_RS][where_slow_cooling_RS] + alpha0SRSforward * weight2[where_intermid_ind_RS][where_slow_cooling_RS]) / weight[where_intermid_ind_RS][where_slow_cooling_RS]

                        alphanuRS[where_slow_cooling_RS] = alpha0SRScentral * ( (nuPrimCentral[where_intermid_ind_RS][where_slow_cooling_RS]<=numRSCentral[where_slow_cooling_RS])*(nuPrimCentral[where_intermid_ind_RS][where_slow_cooling_RS]/numRSCentral[where_slow_cooling_RS])**(-5/3.) + ((numRSCentral[where_slow_cooling_RS]<nuPrimCentral[where_intermid_ind_RS][where_slow_cooling_RS])&(nuPrimCentral[where_intermid_ind_RS][where_slow_cooling_RS]<nucRSCentral[where_slow_cooling_RS])) * (nuPrimCentral[where_intermid_ind_RS][where_slow_cooling_RS]/numRSCentral[where_slow_cooling_RS])**(-(ModVar.pRS+4)/2.) + (nucRSCentral[where_slow_cooling_RS]<=nuPrimCentral[where_intermid_ind_RS][where_slow_cooling_RS]) * (nucRSCentral[where_slow_cooling_RS]/numRSCentral[where_slow_cooling_RS])**(-(ModVar.pRS+4)/2.) * (nuPrimCentral[where_intermid_ind_RS][where_slow_cooling_RS]/nucRSCentral[where_slow_cooling_RS])**(-(ModVar.pRS+5)/2.) )

                        

                        #alphanuRS = alpha0FRScentral * (numRSCentral>nucRSCentral) * ( (nuPrimCentral[where_intermid_ind_RS]<=nucRSCentral) * (nuPrimCentral[where_intermid_ind_RS]/nucRSCentral)**(-5/3.) + ((nucRSCentral<nuPrimCentral[where_intermid_ind_RS])&(nuPrimCentral[where_intermid_ind_RS]<numRSCentral)) * (nuPrimCentral[where_intermid_ind_RS]/nucRSCentral)**(-3) + (numRSCentral<=nuPrimCentral[where_intermid_ind_RS]) * (numRSCentral/nucRSCentral)**-3*(nuPrimCentral[where_intermid_ind_RS]/numRSCentral)**(-(ModVar.pRS+5)/2.) )    +    alpha0SRSforward * (numRSCentral<=nucRSCentral) * ( (nuPrimCentral[where_intermid_ind_RS]<=numRSCentral)*(nuPrimCentral[where_intermid_ind_RS]/numRSCentral)**(-5/3.) + ((numRSCentral<nuPrimCentral[where_intermid_ind_RS])&(nuPrimCentral[where_intermid_ind_RS]<nucRSCentral)) * (nuPrimCentral[where_intermid_ind_RS]/numRSCentral)**(-(ModVar.pRS+4)/2.) + (nucRSCentral<=nuPrimCentral[where_intermid_ind_RS]) * (nucRSCentral/numRSCentral)**(-(ModVar.pRS+4)/2.) * (nuPrimCentral[where_intermid_ind_RS]/nucRSCentral)**(-(ModVar.pRS+5)/2.) )

                        tauRS = alphanuRS * (thicknessRS_behind * weight1[where_intermid_ind_RS] + thicknessRS_forward * weight2[where_intermid_ind_RS]) / 2 / weight[where_intermid_ind_RS]     +    2*tau[where_intermid_ind_RS]#tau is always negative

                        tauRS[np.where(tauRS<0)] = 0.

                        where_low_tau_RS = np.where((tauRS > 1e-10) & (tauRS<1e-2))
                        where_high_tau_RS = np.where(tauRS>1e-2)

                        PRSprimTemp[where_intermid_ind_RS][where_low_tau_RS] *= (tauRS[where_low_tau_RS]**2/2 + tauRS[where_low_tau_RS]**3/6) / tauRS[where_low_tau_RS]
                        PRSprimTemp[where_intermid_ind_RS][where_high_tau_RS] *= (1-np.exp(-tauRS[where_high_tau_RS])) / tauRS[where_high_tau_RS]
                        
                        where_low_tau = np.where((tau > 1e-10) & (tau<1e-2))
                        where_high_tau = np.where(tau > 1e-2)

                        PprimTemp[where_low_tau] *= (tau[where_low_tau] - tau[where_low_tau]**2/2 + tau[where_low_tau]**3/6) / tau[where_low_tau]
                        PprimTemp[where_high_tau] *= (1-np.exp(-tau[where_high_tau]))/tau[where_high_tau]
                """
                
                #angle_integ = (np.cos(xAngInt[:-1]) - np.cos(xAngInt[1:])) * phiInter   #Angular integration segments
                if not Plot_Exceptions.RS_only and not Plot_Exceptions.FS_only:
                    if UseOp.reverseShock:
                        PprimTot = PprimTemp + PRSprimTemp
                    else:
                        PprimTot = PprimTemp
                elif Plot_Exceptions.FS_only or not UseOp.reverseShock:
                    PprimTot = PprimTemp
                elif Plot_Exceptions.RS_only:
                    PprimTot = PRSprimTemp

                F[rimI] = np.trapz(PprimTot * phiInter  ,  Phi) * distance_factor

                if UseOp.runOption == 'LC' and not UseOp.createMock:
                    if not Plot_Exceptions.FS_only and not Plot_Exceptions.RS_only: ### Plot as usual
                        exec('Flux.FFS_%d[rimI] = np.trapz(PprimTemp * phiInter , Phi) * distance_factor'%nuIte)
                        if UseOp.reverseShock:
                            exec('Flux.FRS_%d[rimI] = np.trapz(PRSprimTemp * phiInter , Phi) * distance_factor'%nuIte)
                    elif Plot_Exceptions.FS_only:
                        exec('Flux.FFS_%d[rimI] = np.trapz(PprimTemp * phiInter , Phi) * distance_factor'%nuIte)
                    elif Plot_Exceptions.RS_only:
                        exec('Flux.FRS_%d[rimI] = np.trapz(PRSprimTemp * phiInter , Phi) * distance_factor'%nuIte)
                    else:
                        raise NameError('This exception should not be entered')

    

        if UseOp.plotComponents and UseOp.runOption == 'LC':
            if not Plot_Exceptions.RS_only:
                exec('plt.plot(tobsRed/PlotDetails.scalePlotTime[UseOp.daysOrSec] , Flux.FFS_%d * PlotDetails.scaleFluxAxis[UseOp.fluxAxis] , \'%%s--\'%%PlotDetails.colourCycle[nuIte])'%nuIte)
            if UseOp.reverseShock:
                if not Plot_Exceptions.FS_only:
                    exec('plt.plot(tobsRed/PlotDetails.scalePlotTime[UseOp.daysOrSec] , Flux.FRS_%d * PlotDetails.scaleFluxAxis[UseOp.fluxAxis] , \'%s:\'%PlotDetails.colourCycle[nuIte])'%nuIte)


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
