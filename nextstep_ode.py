def nextstep_func(R , params , ModVar, UseOp , gamma_c_w ,  gamma_c_w_RS=None , shutOff=None):
    import numpy as np
    from scipy import optimize
    from gammamin_fzero import minimize_gammamin


    """
    0 - tburst
    1 - tcomoving
    2 - Gamma
    3 - Eint2
    4 - Eint3
    5 - theta
    6 - Erad2
    7 - Erad3
    8 - Esh2
    9 - Esh3
    10 - Ead2
    11 - Ead3
    12 - M2
    13 - M3
    14 - deltaR4
    """

    adiabaticLosses = True
 
    c,pi,mp,me,kB = 2.9979e10,np.pi, 1.6726e-24, 9.1094e-28 , 1.38065e-16
    qe,mppme = 4.803204e-10,mp+me
    sigma_B = 5.6704e-5  ### Stephan Boltzann constant in cgs                  
    mup,mue = 1.,me / mp
    sigmaT = 6.6524e-25

    """
    ### Booleans
    UseOp.radiativeLosses = booleans[0]
    if UseOp.radiativeLosses:
        remix_radloss = booleans[1]
        continuous_radloss = booleans[2]
        fixed_epsilon = booleans[3]
    reverseShock = booleans[4]
    if reverseShock:
        UseOp.exponential_outflow = booleans[5]    
    """

    ### Setting variables from params array
    tburst = params[0]
    tcomoving = params[1]
    Gamma = params[2]
    
    if Gamma < 1 or Gamma > ModVar.Gamma0:
        return [-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10]
    
    Eint2 = params[3]
    Eint3 = params[4]
    if Eint3 < 0:
        Eint3 = 0.
    """
    if Eint2 < 0 or Eint3 < 0:
        if Eint2 < 0:
            print 'Eint2 < 0'
        if Eint3 < 0:
            print 'Eint3 < 0'
        print 'R =',R
        #return [-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10,-1e10]
        #raw_input()
        #return [float('-inf')]*15
    """
    theta = params[5]
    M2 = params[12]
    M3 = params[13]
    deltaR4 = params[14]


    ###################################################
    ### Setting variables from analytic expressions ###
    ###################################################

    ### Lorentz factors
    #one_over_beta = 1 + 0.5/Gamma**2 + 3/8/Gamma**4 + 5/16/Gamma**6 + 35/128/Gamma**8
    #one_over_beta0 = 1 + 0.5/ModVar.Gamma0**2 + 3/8/ModVar.Gamma0**4 + 5/16/ModVar.Gamma0**6 + 35/128/ModVar.Gamma0**8
    beta = np.sqrt(1-Gamma**-2)
    beta0 = np.sqrt(1-ModVar.Gamma0**-2)
    
    dGammaEffdGamma = 4/3. +1/Gamma**2/3 + 2/Gamma**3/3  ### Note from 2016-06-21
    gammaAdi = (4+1/Gamma) / 3.
    GammaEff = (gammaAdi*Gamma**2 - gammaAdi+1)/Gamma
    if UseOp.reverseShock:

        #if (Gamma > 0.9*ModVar.Gamma0) and (ModVar.Gamma0 > 10):
            #beta43 = (ModVar.Gamma0**2/2+ModVar.Gamma0**2/8/Gamma**2 - Gamma**2/2-Gamma**2/8/ModVar.Gamma0**2) / (ModVar.Gamma0**2/2 + ModVar.Gamma0**2/8/Gamma**2 + Gamma**2/2 + Gamma**2/8/ModVar.Gamma0**2 - 0.25)
            #gamma43_minus_one = 0.5*beta43**2 + 3/8.*beta43**4 + 5/16.*beta43**6
        if Gamma > 0.999995*ModVar.Gamma0:
            gamma43_minus_one = 0.
        else:
            gamma43_minus_one = Gamma*ModVar.Gamma0 * (1/ModVar.Gamma0**2 + 1/Gamma**2 - 1/Gamma**2/ModVar.Gamma0**2) / (1+beta*beta0) - 1
        #else:
            #gamma43_minus_one = Gamma*ModVar.Gamma0*(1-beta*beta0) - 1
        gammaAdi3 = (4+1/(gamma43_minus_one+1)) / 3.
        dgamma43dGamma = 0.5/ModVar.Gamma0 + 0.5/Gamma**2*(0.5/ModVar.Gamma0 - ModVar.Gamma0 - 3*ModVar.Gamma0/8/Gamma**2) ### Note from 2016-06-21. Using asymptotic approximation of beta
        GammaEff3 = (gammaAdi3*Gamma**2 - gammaAdi3+1)/Gamma
        dGammaEff3dGamma = gammaAdi3*(1+Gamma**-2) - dgamma43dGamma/3/(gamma43_minus_one+1)**2*(Gamma-1/Gamma) - Gamma**-2  ### Note from 2016-06-21


    ### Angles
    if theta > 1e-3:
        one_minus_costheta = 1 - np.cos(theta)
    else:
        one_minus_costheta = theta**2/2 - theta**4/24 + theta**6/120
    if ModVar.theta0 > 1e-3:
        one_minus_costheta0 = 1 - np.cos(ModVar.theta0)
    else:
        one_minus_costheta0 = ModVar.theta0**2/2 - ModVar.theta0**4/24 + ModVar.theta0**6/120
        
    if theta < pi/2:
        vperp = np.sqrt( gammaAdi*(gammaAdi-1)*(Gamma-1) / (1 + gammaAdi*(Gamma-1)) ) * c
        dthetadR =vperp/R/Gamma/beta/c 
    else:
        dthetadR = 0.


    ### Densities and volumes
    if ModVar.s == 0: ### CM
        rho = ModVar.A0 / ModVar.M0 * mppme
        dlnrho1dR = 0.
    else:
        if R < ModVar.R_ISM:
            rho = ModVar.A0 / ModVar.M0 * R**(-ModVar.s) * mppme
            dlnrho1dR = -ModVar.s/R
        else:
            rho = ModVar.A0 / ModVar.M0 * ModVar.R_ISM**(-ModVar.s) * mppme
            dlnrho1dR = 0.
        
    rho2 = 4 * Gamma * rho

    dM2dR = 2*pi*R**(2.-ModVar.s)*ModVar.A0 / ModVar.M0*(mp+me)*one_minus_costheta
    
    V2 = M2 / rho2
    if Eint2 > 0:
        B = np.sqrt(8*pi*ModVar.eB)*np.sqrt(Eint2)/np.sqrt(V2) * c * np.sqrt(ModVar.M0)
    else:
        B = 0.
    if UseOp.reverseShock:

        if UseOp.exponential_outflow and (Gamma < ModVar.Gamma0):

            ### ddeltaR4dR
            """
            if Gamma > 0.8*ModVar.Gamma0:
                #deltaR4 = (one_over_beta0**-2*c**2*tburst**2 - R**2) / (c*tburst/one_over_beta0 + R)
            else:
                deltaR4 = c*tburst/one_over_beta0 - R
            """
            alpha_of = ModVar.tprompt*beta0*c
            xi04 = 1 / alpha_of
            ddeltaR4dR = (1/beta**4 - 1/beta0**4) / (1/beta + 1/beta0) / (1/beta**2+1/beta0**2) *beta0
            #ddeltaR4dR = (1/Gamma - 1/ModVar.Gamma0) * (1/Gamma + 1/ModVar.Gamma0) / (beta * (beta+beta0))

            rho4_scale_factor = -deltaR4/alpha_of
            rho4_factor = xi04 / 2/np.pi/R**2/one_minus_costheta0
            ### Injection rate. See Johansson et al. (2016) for details
            #if rho4_scale_factor < -0.1: ### To ensure accuracy for exponents in the vicinity of zero

            #if rho4_scale_factor < 100:

            dM3dR = 1/alpha_of * np.exp(rho4_scale_factor) 


            dM3dR *= ddeltaR4dR
            rho4 = rho4_factor*np.exp(rho4_scale_factor)
                #M3 = ModVar.M0 * (1-np.exp(rho4_scale_factor))
            #else:
            #dM3dR = ModVar.M0/alpha_of * ddeltaR4dR * (1 + rho4_scale_factor + rho4_scale_factor**2/2 + rho4_scale_factor**3/6 + rho4_scale_factor**4/24 + rho4_scale_factor**5/120)
            #rho4 = rho4_factor*(1 + rho4_scale_factor + rho4_scale_factor**2/2 + rho4_scale_factor**3/6 + rho4_scale_factor**4/24 + rho4_scale_factor**5/120)
                #M3 = -ModVar.M0 * (rho4_scale_factor + rho4_scale_factor**2/2 + rho4_scale_factor**3/6 + rho4_scale_factor**4/24 + rho4_scale_factor**5/120)

            rho3prim = 4*Gamma*rho4
            """
            if rho3prim < 1e-70 or (M3 > (R**3 * one_minus_costheta)*rho3prim):
                shutOff = True
            else:
            """
            if not shutOff:
                V3 = M3 / rho3prim

                if M3 > 0 and Eint3 > 0:
                    BRS = np.sqrt(8*pi*ModVar.eB3)*np.sqrt(Eint3)/np.sqrt(V3) * c * np.sqrt(ModVar.M0)
                else:
                    BRS = 0.
                dlnrho4dR = -ddeltaR4dR /alpha_of - 2/R
            else:
                dM3dR = 0.
                rho4 = 0.
                BRS = 0.
                dlnrho4dR = 0.
                
        elif not UseOp.exponential_outflow: ### Constant outflow
            if ( M3 <= 1 ):  

                bb0m1 = 0.5/Gamma**2 + 0.5/ModVar.Gamma0**2 + Gamma**(-4)/8 + ModVar.Gamma0**(-4)/8 - 0.25/(ModVar.Gamma0**2)/(Gamma**2)
                dM3dR = 1 * (beta0**2 - beta**2)/(beta0 + beta) / beta0 / ModVar.tprompt / beta / c / bb0m1
                rho4 = 1 / (2*pi*R**2*ModVar.tprompt*one_minus_costheta*beta*c)
                dlnrho4dR = -2/R
            else:

                dM3dR = 0.
                rho4 = 0.
        else:
            dM3dR = 0.
            dlnrho4dR = 0.
            rho4 = 0.
            #M3 = 0.
            ddeltaR4dR = 0.
            """
            if Gamma == ModVar.Gamma0:
                shutOff = False
            else:
                shutOff = True
            """ 
    else:
        dM3dR = 0.
        dlnrho4dR = 0.

        ddeltaR4dR = 0.




    ### Mass derivatives

    

    ### dGammadR
    f_2 = GammaEff*(gammaAdi-1)*Eint2/Gamma
    h_2 = GammaEff*(gammaAdi-1)*Eint2*(dM2dR/M2-dlnrho1dR)

    if UseOp.reverseShock:
        fh_factor3 = GammaEff3*(gammaAdi3-1)*Eint3
        f_3 = fh_factor3 / (gamma43_minus_one+1) * dgamma43dGamma
        if Eint3 != 0:
            #try:
            h_3 = fh_factor3 * (dM3dR/M3 - dlnrho4dR)
            #except:
            #print 'M3 =',M3
            #    print 'deltaR4dR =',deltaR4
            #    raw_input(Gamma/ModVar.Gamma0)
        else:
            h_3 = 0.
        """
        print 'Gamma-ModVar.Gamma0+GammaEff3*gamma43_minus_one = %s'%(Gamma-ModVar.Gamma0+GammaEff3*gamma43_minus_one)
        print (Gamma-1)*(GammaEff+1)*dM2dR*c**2
        print (Gamma-ModVar.Gamma0+GammaEff3*gamma43_minus_one)*dM3dR*c**2
        print 'h_2 =',h_2
        print 'h_3 =',h_3
        print ((M2+M3)*c**2 + Eint2*dGammaEffdGamma + Eint3*dGammaEff3dGamma + f_2 + f_3)
        """
        dGammadR = -(  (Gamma-1)*(GammaEff+1)*dM2dR + (Gamma-ModVar.Gamma0+GammaEff3*gamma43_minus_one)*dM3dR - h_2 - h_3   ) / ((M2+M3) + Eint2*dGammaEffdGamma + Eint3*dGammaEff3dGamma + f_2 + f_3)
    else:
        dGammadR = -((Gamma-1)*(GammaEff+1)*dM2dR - h_2) / ((1+M2) + Eint2*dGammaEffdGamma + f_2)
    
    if dGammadR > 0:
        if Gamma > 0.95*ModVar.Gamma0:
            dGammadR = 0.


    
    ################
    ### Energies ###
    ################
    
    ### Shocked energy
    dEsh2dR = (Gamma-1)*dM2dR
    ### Expansion energy
    
    dlnV2dR = dM2dR/M2 - dlnrho1dR - dGammadR/Gamma
    if adiabaticLosses:
        dEad2dR = -(gammaAdi-1)*Eint2*dlnV2dR
    else:
        dEad2dR = 0.
    ### Radiative losses
    if UseOp.radiativeLosses and not fixed_epsilon and Eint2 > 0:
    
        dErad2dR = 0.

        #gamma_c_w = 6*pi*me*c/sigmaT/B**2/tcomoving
        #if gamma_c_w < 1:
            #gamma_c_w = 1.
            
        gamma_max = (6*pi*qe/sigmaT/B)**.5
        #gamma_min = (p-2)/(p-1)*(ModVar.epsilone/mue*(Gamma-1)+1)


        ### Finding gamma_min
        #gamma_min = minimize_gammamin(gamma_c_w , gamma_max , p , ModVar.epsilone , Gamma-1 , mue)
        """
        gmin_func = lambda gmin: gmin_fzero(gmin , gamma_c_w , gamma_max , p , ModVar.epsilone , Gamma-1 , mue , False)
        try:
            gamma_min = optimize.newton(gmin_func , gamma_min)
        except:

            res = optimize.fsolve(gmin_func,gamma_min)  
            if len(res) > 1:  ### Testing for multiple solutions. Hasn't occured yet
                print 'lenght of res > 1 !!!!!!!!!!!!!!'
            else:
                gamma_min = res[0]

        """
        if gamma_min_inj > 0.9995*gamma_max:
            N0_inj = 0.
        else:
            #if gamma_c_w > gamma_min_inj:
            N0_inj = dM2dR * ModVar.epsilone / me * (Gamma-1) / ((gamma_max**(2-p)-gamma_min_inj**(2-p))/(2-p) - (gamma_max**(1-p)-gamma_min_inj**(1-p))/(1-p))
            #else: ### Fast cooling
            #    N0_inj = dM2dR / mp / ((gamma_min_inj**(-p)-gamma_max**(-p))/p + gamma_min_inj**(1-p)*(gamma_c_w**(-1)-gamma_min_inj**(-1)))
        """
        if gamma_c_w > gamma_min:
            N0_inj = dM2dR /mp / ((gamma_c_w**(-p-1)-gamma_min**(-p+1)) / (1-p) - (gamma_c_w*gamma_max**(-p) - gamma_c_w**(-p+1))/p)
        else:
            #N0_inj = dM2dR /mp / ((gamma_max**(1-p)-gamma_min**(1-p))/(1-p) + gamma_min**(2-p)/gamma_c_w - gamma_min**(1-p))
            N0_inj = dM2dR /mp / ((gamma_max**(1-p)-gamma_min**(1-p))/(1-p))
        """
            
        if gamma_c_w > gamma_max:
            gamma_c_w = gamma_max

        if remix_radloss and gamma_c_w > 0:
            if gamma_c_w < gamma_max: ### No losses when gamma_c_w > gamma_max
                #dErad2dR += N0_inj*me * (gamma_min_inj**(2-p)/(p-2)  -  (gamma_c_w*gamma_max**(1-p)-gamma_c_w**(2-p))/(1-p)  -  (gamma_c_w**(2-p)-gamma_min_inj**(2-p))/(2-p))
                if gamma_c_w > gamma_min_inj: ### Slow cooling
                    dErad2dR += N0_inj * me * ((gamma_max**(2-ModVar.p)-gamma_c_w**(2-ModVar.p))/(2-ModVar.p) - (gamma_max**(1-ModVar.p)-gamma_c_w**(1-ModVar.p))/(1-ModVar.p) - gamma_c_w*((gamma_max**(1-ModVar.p)-gamma_c_w**(1-ModVar.p))/(1-ModVar.p) + (gamma_max**-ModVar.p-gamma_c_w**-ModVar.p)/ModVar.p))
                else:


                    dErad2dR += N0_inj * me * ((gamma_max**(2-ModVar.p)-gamma_min_inj**(2-ModVar.p))/(2-ModVar.p) - (gamma_max**(1-ModVar.p)-gamma_min_inj**(1-ModVar.p))/(1-ModVar.p) - gamma_min_inj*((gamma_max**(1-ModVar.p)-gamma_min_inj**(1-ModVar.p))/(1-ModVar.p) + (gamma_max**-ModVar.p-gamma_min_inj**-ModVar.p)/ModVar.p))
                    #dErad2dR += N0_inj * me * ((gamma_max**(2-p)-gamma_min_inj**(2-p))/(2-p) - 2*(gamma_max**(1-p)-gamma_min_inj**(1-p))/(1-p) - gamma_min_inj**(1-p)*(np.log(gamma_min/gamma_c_w) + gamma_min - gamma_c_w) + (gamma_min**(-p)-gamma_max**(-p))/p)
                    #dErad2dR += N0_inj * me * ((gamma_max**(2-p)-gamma_min**(2-p))/(2-p) - (gamma_max**(1-p)-gamma_min**(1-p))/(1-p)   - (gamma_min**(1-p)*(np.log(gamma_min/gamma_c_w) + gamma_min - gamma_c_w) + (gamma_max**(1-p)-gamma_min**(1-p))/(1-p) + (gamma_max**(-p)-gamma_min**(-p))/p))


        if continuous_radloss and gamma_c_w > 0:
            if gamma_c_w > gamma_min: ### Slow cooling
                N0_total = M2/mp*((gamma_c_w**(1-ModVar.p)-gamma_min**(1-ModVar.p))/(1-ModVar.p) - gamma_c_w*(gamma_max**(-ModVar.p)-gamma_c_w**(-ModVar.p))/ModVar.p)
                #N0_total = M2/mp*(p-1)*gamma_min**(p-1)
                #N0_total = Eint2 / c**2 /  me * (p-2) * gamma_min**(p-2) * ModVar.epsilone
                if gamma_c_w > gamma_max:
                    dErad2dR += N0_total*sigmaT*B**2/6/pi * (gamma_max**(3-ModVar.p)-gamma_min**(3-ModVar.p))/(3-ModVar.p) / beta/c**2/Gamma
                else:
                    dErad2dR += N0_total*sigmaT*B**2/6/pi * ((gamma_c_w**(3-ModVar.p)-gamma_min**(3-ModVar.p))/(3-ModVar.p)  +  (gamma_c_w*gamma_max**(2-ModVar.p) - gamma_c_w**(3-ModVar.p))/(2-ModVar.p)) / beta/c**2/Gamma
            else: ### Fast cooling
                N0_total = M2 / mp * ModVar.p / (gamma_min**(-ModVar.p) - gamma_max**(-ModVar.p))
                #N0_total = M2 / mp * (gamma_min**(1-p)*(gamma_c_w**-1 - gamma_min**-1) + gamma_min**(-p)/p)
                #N0_total = M2 / mp * gamma_min**(-p)/p
                #dErad2dR += N0_total * sigmaT*B**2/6/pi * (gamma_min**(1-p)*(gamma_min-gamma_c_w) + (gamma_max**(2-p)-gamma_min**(2-p))/(2-p)) / beta/c**2/Gamma
                dErad2dR += N0_total * sigmaT*B**2/6/pi * (gamma_max**(2-ModVar.p)-gamma_min**(2-ModVar.p))/(2-ModVar.p) / beta/c**2/Gamma
        if UseOp.reverseShock and dM3dR > 0 and not shutOff:   ### Radiative losses for RS
            dErad3dR = 0.

            if BRS > 1e-70:
                
                gamma_max_RS = (6*pi*qe/sigmaT/BRS)**.5
                #gamma_min_RS = (ModVar.pRS-2)/(ModVar.pRS-1)*(ModVar.epsilone3/mue*gamma43_minus_one+1)
                #gamma_c_w_RS = 6*pi*me*c/sigmaT/BRS**2/tcomoving
                #raw_input('gamma_c_w_RS = %s'%gamma_c_w_RS)

                #if gamma_c_w_RS < 1:
                    #gamma_c_w_RS = 1.
                if gamma_c_w_RS > gamma_max_RS:
                    gamma_c_w_RS = gamma_max_RS

                #gamma_min_RS = minimize_gammamin(gamma_c_w_RS , gamma_max_RS , ModVar.pRS , ModVar.epsilone3 , gamma43_minus_one , mue)

                """
                gmin_func_RS = lambda gmin: gmin_fzero(gmin , gamma_c_w_RS , gamma_max_RS , ModVar.pRS , ModVar.epsilone3 , gamma43_minus_one , mue , False)

                res = optimize.fsolve(gmin_func_RS,gamma_min_RS)

                gamma_min_RS = res[0]
                
                except:

                    res = optimize.fsolve(gmin_func_RS , gamma_min_RS)
                    
                    if len(res) > 1:
                        print 'lenght of res (RS) > 1 !!!!!!!!!!!!!'
                    else:
                        gamma_min_RS = res[0]
                """

                
                if gamma_min_RS < 1:
                    gamma_min_RS = 1.
                
                #if gamma_c_w_RS > gamma_min_RS: ### Slow cooling
                    #N0_inj_RS = dM3dR / mp / ((gamma_c_w_RS**(-ModVar.pRS+1)-gamma_min_RS**(-ModVar.pRS+1)) / (1-ModVar.pRS) - (gamma_c_w_RS*gamma_max_RS**(-ModVar.pRS) - gamma_c_w_RS**(-ModVar.pRS+1))/ModVar.pRS)
                    #N0_inj_RS = dM3dR / mp * (1-ModVar.pRS) / (gamma_max_RS**(1-ModVar.pRS) - gamma_min_RS**(1-ModVar.pRS))
                if gamma_min_RS_inj > 0.9995*gamma_max_RS:
                    N0_inj_RS = 0.
                else:
                    N0_inj_RS = dM3dR * ModVar.epsilone3 / me * gamma43_minus_one / ((gamma_max_RS**(2-ModVar.pRS)-gamma_min_RS_inj**(2-ModVar.pRS))/(2-ModVar.pRS) - (gamma_max_RS**(1-ModVar.pRS)-gamma_min_RS_inj**(1-ModVar.pRS))/(1-ModVar.pRS))
                """
                else: ### Fast cooling
                    try:
                        N0_inj_RS = dM3dR /mp / ((gamma_max_RS**(1-ModVar.pRS)-gamma_min_RS**(1-ModVar.pRS))/(1-ModVar.pRS))# + gamma_min_RS**(2-ModVar.pRS)/gamma_c_w_RS - gamma_min_RS**(1-ModVar.pRS))
                    except:
                        print 'crashed'
                        print Eint3
                        print BRS
                        print dM3dR
                        print gamma_max_RS
                        print gamma_min_RS
                        print ModVar.pRS
                        print gamma_c_w_RS
                """
                if remix_radloss and not shutOff:
                    if gamma_c_w_RS < gamma_max_RS:  ### Only losses when gamma_c_w_RS is smaller than gamma_max_RS
                        if gamma_c_w_RS > gamma_min_RS_inj:
                            dErad3dR += N0_inj_RS * me * ((gamma_max_RS**(2-ModVar.pRS)-gamma_c_w_RS**(2-ModVar.pRS))/(2-ModVar.pRS) - (gamma_max_RS**(1-ModVar.pRS)-gamma_c_w_RS**(1-ModVar.pRS))/(1-ModVar.pRS) - gamma_c_w_RS*((gamma_max_RS**(1-ModVar.pRS)-gamma_c_w_RS**(1-ModVar.pRS))/(1-ModVar.pRS) + (gamma_max_RS**-ModVar.pRS-gamma_c_w_RS**-ModVar.pRS)/ModVar.pRS))
                            #dErad3dR += N0_inj_RS*me * ((gamma_min_RS_inj**(2-ModVar.pRS))/(ModVar.pRS-2)  -  (gamma_c_w_RS*gamma_max_RS**(1-ModVar.pRS)-gamma_c_w_RS**(2-ModVar.pRS))/(1-ModVar.pRS)  -  (gamma_c_w_RS**(2-ModVar.pRS)-gamma_min_RS_inj**(2-ModVar.pRS))/(2-ModVar.pRS))
                        else:
                            #dErad3dR += N0_inj_RS * me * ((gamma_max_RS**(1-ModVar.pRS)-gamma_min_RS_inj**(1-ModVar.pRS))/(1-ModVar.pRS) + (gamma_max_RS**-ModVar.pRS - gamma_min_RS_inj**-ModVar.pRS)/p)


                            dErad3dR += N0_inj_RS * me * ((gamma_max_RS**(2-ModVar.pRS)-gamma_min_RS_inj**(2-ModVar.pRS))/(2-ModVar.pRS) - (gamma_max_RS**(1-ModVar.pRS)-gamma_min_RS_inj**(1-ModVar.pRS))/(1-ModVar.pRS) - gamma_min_RS_inj*((gamma_max_RS**(1-ModVar.pRS)-gamma_min_RS_inj**(1-ModVar.pRS))/(1-ModVar.pRS) + (gamma_max_RS**-ModVar.pRS-gamma_min_RS_inj**-ModVar.pRS)/ModVar.pRS))


                            #dErad3dR += N0_inj_RS * me * ((gamma_max_RS**(2-ModVar.pRS)-gamma_min_RS_inj**(2-ModVar.pRS))/(2-ModVar.pRS)  -  2*(gamma_max_RS**(1-ModVar.pRS)-gamma_min_RS**(1-ModVar.pRS))/(1-ModVar.pRS) - gamma_min_RS**(1-ModVar.pRS)*(np.log(gamma_min_RS/gamma_c_w_RS) + gamma_min_RS - gamma_c_w_RS) + (gamma_min_RS**(-ModVar.pRS)-gamma_max_RS**(-ModVar.pRS))/ModVar.pRS)
                            #dErad3dR += N0_inj_RS * me * ((gamma_max_RS**(2-ModVar.pRS)-gamma_min_RS**(2-ModVar.pRS))/(2-ModVar.pRS) - (gamma_max_RS**(1-ModVar.pRS)-gamma_min_RS**(1-ModVar.pRS))/(1-ModVar.pRS)   - (gamma_min_RS**(1-ModVar.pRS)*(np.log(gamma_min_RS/gamma_c_w_RS) + gamma_min_RS - gamma_c_w_RS) + (gamma_max_RS**(1-ModVar.pRS)-gamma_min_RS**(1-ModVar.pRS))/(1-ModVar.pRS) + (gamma_max_RS**(-ModVar.pRS)-gamma_min_RS**(-ModVar.pRS))/ModVar.pRS))

                if continuous_radloss:
                    ### Continuous radiation from the entire shocked region. We need to calculate N0 separately, since N0_inj is for the newly shocked particles
                    if gamma_c_w_RS > gamma_min_RS: ### Slow cooling
                        N0_total_RS = M3/mp/((gamma_c_w_RS**(1-ModVar.pRS)-gamma_min_RS**(1-ModVar.pRS))/(1-ModVar.pRS) - gamma_c_w_RS*(gamma_max_RS**(-ModVar.pRS)-gamma_c_w_RS**(-ModVar.pRS))/ModVar.pRS)
                        #N0_total_RS = M3/mp*(p-1)*gamma_min_RS**(p-1)
                        if gamma_c_w_RS > gamma_max_RS:
                            dErad3dR += N0_total_RS*sigmaT*BRS**2/6/pi * (gamma_max_RS**(3-ModVar.pRS)-gamma_min_RS**(3-ModVar.pRS))/(3-ModVar.pRS) / beta/c**2/Gamma
                        else:
                            dErad3dR += N0_total_RS*sigmaT*BRS**2/6/pi * ((gamma_c_w_RS**(3-ModVar.pRS)-gamma_min_RS**(3-ModVar.pRS))/(3-ModVar.pRS)  +  (gamma_c_w_RS*gamma_max_RS**(2-ModVar.pRS) - gamma_c_w_RS**(3-ModVar.pRS))/(2-ModVar.pRS)) /beta/c**2/Gamma
                        
                    else: ### Fast cooling
                        N0_total_RS = M3 / mp * ModVar.pRS / (gamma_min_RS**(-ModVar.pRS)-gamma_max_RS**(-ModVar.pRS))
                        #N0_total_RS = M3 / mp * (gamma_min_RS**(1-ModVar.pRS)*(gamma_c_w_RS**-1 - gamma_min_RS**-1) + gamma_min_RS**(-ModVar.pRS)/ModVar.pRS)
                        #N0_total_RS = M3 / mp * (gamma_min_RS**(-ModVar.pRS)/ModVar.pRS)
                        #dErad3dR += N0_total_RS * sigmaT*BRS**2/6/pi * (gamma_min_RS**(1-ModVar.pRS)*(gamma_min_RS-gamma_c_w_RS) + (gamma_max_RS**(2-ModVar.pRS)-gamma_min_RS**(2-ModVar.pRS))/(2-ModVar.pRS)) / beta/c**2/Gamma
                        dErad3dR += N0_total_RS * sigmaT*BRS**2/6/pi * ((gamma_max_RS**(2-ModVar.pRS)-gamma_min_RS**(2-ModVar.pRS))/(2-ModVar.pRS)) / beta/c**2/Gamma
                """
                print 'M3             =',M3
                print 'gamma_min_RS   =',gamma_min_RS
                print 'gamma_c_w_RS   =',gamma_c_w_RS
                #print 'N0_total_RS    =',N0_total_RS
                print 'N0_inj_RS      =',N0_inj_RS
                print 'gamma_max_RS   =',gamma_max_RS
                print 'BRS            =',BRS
                print 'dErad3dR       =',dErad3dR
                print 'dM3dR          =',dM3dR
                raw_input(Gamma)
                """
        else:
            dErad3dR = 0.
    elif UseOp.radiativeLosses and fixed_epsilon:
        dErad2dR = microParams[6] * dEsh2dR * ModVar.epsilone
        if UseOp.reverseShock:
            dErad3dR = microParams[6] * gamma43_minus_one * dM3dR * ModVar.epsilone3
        else:
            dErad3dR = 0.
    else:
        dErad2dR = 0.
        dErad3dR = 0.
    if UseOp.reverseShock and M3>0:  ### rho4 becomes 0 when RS injecta is cut-off
        ### Shocked energy
        dEsh3dR = gamma43_minus_one*dM3dR

        ### Expansion energy
        #dlnGamma43dR = dGammadR * dgamma43dGamma / (gamma43_minus_one+1)
        dlnV3dR = dM3dR / M3 - dlnrho4dR - dGammadR / Gamma
        """
        print '---'
        print 'Gamma/ModVar.Gamma0 = %s'%(Gamma/ModVar.Gamma0)
        print 'dM3dR / M3   = %s'%(dM3dR / M3)
        print 'dM3dR / ddeltaR4dR = %s'%(dM3dR / ddeltaR4dR)
        print 'e(rho4_scale_factor) = %s'%(np.exp(rho4_scale_factor))
        print 'dlnrho4dR    = %s'%(-dlnrho4dR)
        print 'dM3dR/M3 - dlnrho4dR = %s'%(dM3dR/M3 - dlnrho4dR)
        print 'ddeltaR4dR / alpha_of = %s'%(ddeltaR4dR / alpha_of)
        print '2/R = %s'%(2/R)
        print 'dlnGammadR = %s'%(-dGammadR/Gamma)
        print 'dlnV3dR =',dlnV3dR
        raw_input('Eint3 = %s'%Eint3)
        """
        if adiabaticLosses:
            dEad3dR = -(gammaAdi3-1)*Eint3*dlnV3dR#dgamma43dGamma - h_3 / GammaEff3
        else:
            dEad3dR = 0.
    else:
        dEsh3dR = 0.
        dEad3dR = 0.

    dEint2dR = dEsh2dR + dEad2dR - dErad2dR
    dEint3dR = dEsh3dR + dEad3dR - dErad3dR
    """
    if dErad3dR > dEsh3dR:
        print 'dErad3dR > dEsh3dR'
    if dErad2dR > dEsh2dR:
        print 'dErad2dR > dEsh2dR'
    
    #print 'R =',R
    #if dErad3dR > 0:
    print '---------------------'
    print 'R =',R
    print 'dthetadR =',dthetadR
    print 'dM2dR =',dM2dR
    print 'dM3dR =',dM3dR
    #    print 'dGamma43dGamma =',dgamma43dGamma
    print 'ddeltaR4dR =',ddeltaR4dR    
    print 'dEint2dR =',dEint2dR
    print 'dEint3dR =',dEint3dR
    print 'dEad3dR  =',dEad3dR
    #print 'dEad3dR / Ead3  = %s'%(dEad3dR/Ead3)
    print 'dEsh3dR  =',dEsh3dR
    print 'dErad3dR =',dErad3dR
    try:
        print 'epsilonr = %s'%(dErad3dR / dEsh3dR)
    except:None
    #    print 'gamma43-1 =',gamma43_minus_one
    raw_input('dGammadR = %s'%dGammadR)

    if Eint2 < 0 or Eint3 < 0:
        print 'dEsh3dR  =',dEsh3dR
        print 'dErad3dR =',dErad3dR
        print 'dEint3dR =',dEint3dR
        print 'dGammadR =',dGammadR
        print 'deltaR4  =',deltaR4
        print 'Gamma/ModVar.Gamma0 = %s'%(Gamma/ModVar.Gamma0)
        print 'Gamma43 - 1 =',gamma43_minus_one
        if R > 1e13:
            if Eint2 < 0:
                print 'B =',B
            else:
                print 'BRS =',BRS
                print 'rho4 =',rho4
                print 'shutOff =',shutOff
    """            
    #raw_input(np.array([1/beta/c  ,  1/beta/Gamma/c  ,  dGammadR  , dEint2dR , dEint3dR , dthetadR , dErad2dR , dErad3dR , dEsh2dR , dEsh3dR , dEad2dR , dEad3dR , dM2dR , dM3dR , ddeltaR4dR]))

    return np.array([1/beta/c  ,  1/beta/Gamma/c  ,  dGammadR  , dEint2dR , dEint3dR , dthetadR , dErad2dR , dErad3dR , dEsh2dR , dEsh3dR , dEad2dR , dEad3dR , dM2dR , dM3dR , ddeltaR4dR])

#    else:
#        return np.array([1*one_over_beta/c  ,  1*one_over_beta/Gamma/c  ,  dGammadR  , dEint2dR , dEint3dR , dthetadR , dErad2dR , dErad3dR ,  dEsh2dR , dEad2dR , dM2dR , ddeltaR4dR])
