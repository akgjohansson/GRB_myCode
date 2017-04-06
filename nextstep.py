def next_step(full_step , Gamma , R , M2 , M0 , deltaR , tburst , Gamma0 , beta0 , useSpread , reverseShock , radiativeLosses , fixed_epsilon , gammaAdi_minus_one , theta, theta0 , R_ISM , s , A0 , rho , p , epsilone , eB , Eint2 , gamma_c_w , dGammadR , GammaEff , dGammaEffdGamma, print_error,shutOff=None  , tprompt=None , xi04=None , M3=None , gamma43_minus_one=None , gammaAdi3_minus_one=None , delta_tco=None , pRS=None , epsilone3=None , eB3=None , Eint3=None , gamma_c_w_RS=None , GammaEff3=None , dGammaEff3dGamma=None , exponential_outflow=None , dgamma43dGamma=None):
    import numpy as np
    from scipy import optimize
    from matplotlib import pyplot as plt
    from test_epsilon_rad_fit import test_epsrad
    from gammamin_fzero import gmin_fzero
    c,pi,mp,me,kB = 2.9979e10,np.pi, 1.6726e-24, 9.1094e-28 , 1.38065e-16
    qe = 4.803204e-10
    sigma_B = 5.6704e-5  ### Stephan Boltzann constant in cgs                  

    print_error = True
 

    mup,mue = 1.,me / mp
    sigmaT = 6.6524e-25
    beta = (1 - Gamma**(-2))**0.5
    one_over_beta = 1 + 0.5/Gamma**2 + 3/8/Gamma**4 + 5/16/Gamma**6 + 35/128/Gamma**8
    one_over_beta0 = 1 + 0.5/Gamma0**2 + 3/8/Gamma0**4 + 5/16/Gamma0**6 + 35/128/Gamma0**8

    ### Boolean variables for testing of radiative loss calculation
    cooltime_integ = False
    continuous_radloss = True
    remix_radloss = True

    if theta < 0.1:
        one_minus_costheta = theta**2/2 - theta**4/24 + theta**6/120
    else:
        one_minus_costheta = 1 - np.cos(theta)

    dM2dR = 2*pi*R**(2.-s)*A0*(mp+me)*one_minus_costheta 
    dM2 = 2*pi*A0*(mp+me)*((R+deltaR)**(3-s) - R**(3-s))*one_minus_costheta / (3-s)
    M2_next = M2 + dM2

    if useSpread:
        vperp = ( (gammaAdi_minus_one+1)*gammaAdi_minus_one*(Gamma-1) / (1 + (gammaAdi_minus_one+1)*(Gamma-1)) )**(1/2.) * c
        dthetadR = 2*vperp/(2*R+deltaR)/Gamma*one_over_beta/c               ###(theta-theta[i-1])/(R-R[i-1]))
        if theta >= np.pi / 2: 
            theta_next = np.pi/2
        else:
            theta_next = theta + dthetadR * deltaR    ###2*vperp / (R+R[i+1]) * (tcomoving[i+1]-tcomoving[i])    #Spreading jet
    else:
        dthetadR = np.append(dthetadR,0.)
        theta_next = theta

    #Injection of mass. Note from 11/2 2014
    if reverseShock and not shutOff:
        shutOff_next = False
        if exponential_outflow:
            deltaR4 = tburst*beta0*c - R

            if deltaR4 < 0: ### Jet is accererating, and is not accreting material from outflow
                print 'deltaR4 < 0'
                print Gamma
                print Gamma0
                print tburst
                print beta0
                print c
                print deltaR4
                raw_input(R)
                dM3dR = 0.
                rho4 = 0.
                print 'deltaR4 < 0'
            else:
                alpha_of = tprompt*beta0*c
                                    
                    
                rho4_scale_factor = -deltaR4/alpha_of
                rho4_factor = xi04 / 2/np.pi/R**2/one_minus_costheta

                ### Injection rate. See Johansson et al. (2016) for details
                if rho4_scale_factor < -0.1:
                    beta0dbeta = beta0*one_over_beta
                    if beta0dbeta > 1.1:
                        dM3dR = M0/alpha_of * np.exp(rho4_scale_factor) * (beta0dbeta - 1)
                    else:
                        dM3dR = M0/alpha_of * np.exp(rho4_scale_factor) * (one_over_beta**2 - one_over_beta0**2)/(one_over_beta + one_over_beta0)
                else:
                    dM3dR = M0/alpha_of * (one_over_beta**2 - one_over_beta0**2)/(one_over_beta + one_over_beta0) * (1 + rho4_scale_factor + rho4_scale_factor**2/2 + rho4_scale_factor**3/6 + rho4_scale_factor**4/24 + rho4_scale_factor**5/120)




                M3_next = M0 * (1 - np.e**(rho4_scale_factor))
                
                dM3 = dM3dR * deltaR
                rho4 = rho4_factor*np.exp(rho4_scale_factor)
                if rho4 == 0: 
                    print 'rho4 == 0!!! epsilone = ',epsilone
                    print rho4_scale_factor
                    print rho4_factor_log
                    print 'R =',R
                    print 'xi04 =',xi04
                if rho4 < 0: print 'rho4 < 0!',rho4

                """
                if (rho4_factor_log + np.e**rho4_scale_factor < -70) or M3 >= M0: ### If this is False, the RS feed is as close to finished as it will get

                    shutOff_next = True
                """

                    
                        
        else:

            if ( M3 <= M0 ):  
                bb0m1 = 0.5/Gamma**2 + 0.5/Gamma0**2 + Gamma**(-4)/8 + Gamma0**(-4)/8 - 0.25/(Gamma0**2)/(Gamma**2)
                dM3dR = M0 * (beta0 - beta) * one_over_beta0 / tprompt / beta / c / bb0m1
                dM3 = dM3dR * deltaR
                M3_next = M3 + dM3
                if dM3 < 0:
                    print 'dM3 < 0'
                    print 'beta/betai = %s'%(beta0*one_over_beta)
                    print 'delta tco = %s'%delta_tco
                    raise RuntimeError()

                        
                
                rho4 = M0 / (2*pi*R**2*tprompt*one_minus_costheta*beta*c)
            
                
            if M3_next > M0:

                if shutOff == False: 
                    
                    shutOff_next = True

                dM3 = 0.

        rho3prim = 4 * rho4 * Gamma

        #M3_next = M3 + dM3            

    elif reverseShock:
        
        shutOff_next = True
        M3_next = M3
        dM3dR = 0.
        dM3 = 0.
        rho4 = 0.
        dgamma43dGamma_next = 0.
        gammaAdi3_minus_one_next = 4/3.


    #######################################
    ### Radiative efficiency estimation ###
    #######################################


    rho2prim = 4*Gamma*rho
    V2 = M2 / rho2prim
    if reverseShock and not shutOff: 
        V3 = M3 / rho3prim
        if V3 > (R**3 * one_minus_costheta):
            shutOff_next = True
    
    B = np.sqrt(8*pi*eB*Eint2/V2)
    if np.isnan(B):
        print 'Eint2 =',Eint2
        print 'V2 =',V2
        raw_input()
        raise RuntimeError()
    ### New radiative losses estimation. See draft from 2015-10-14 or later
    gamma_max = (6*pi*qe/sigmaT/B)**.5
    gamma_min = (p-2)/(p-1)*(epsilone/mue*(Gamma-1)+1)
    if gamma_c_w > 0:
        """
        if gmin_fzero(gamma_c_w , gamma_c_w , gamma_max , p , epsilone , Gamma-1 , mue , False) > 0: ### Slow cooling
            gmin_func = lambda gmin: gmin_fzero(gmin , gamma_c_w , gamma_max , p , epsilone , Gamma-1 , mue , False)
            gamma_min = optimize.newton(gmin_func,gamma_min)
        else:
            gmin_func = lambda gmin: gmin_fzero(gmin , gamma_c_w , gamma_max , p , epsilone , Gamma-1 , mue , True)
            gamma_min = optimize.newton(gmin_func,gamma_c_w+1)
        """
        gmin_func = lambda gmin: gmin_fzero(gmin , gamma_c_w , gamma_max , p , epsilone , Gamma-1 , mue , False)
        try:
            gamma_min = optimize.newton(gmin_func , gamma_min)
        except:

            res = optimize.fsolve(gmin_func,gamma_min)  
            if len(res) > 1:
                print 'lenght of res > 1 !!!!!!!!!!!!!!'
                raw_input(res)
            else:
                gamma_min = res[0]

        
    #gamma_min = (p-2)/(p-1)*(Eint2/M2/mue*epsilone/c**2 + 1)



    ### estimating radiative losses. This only applies if the radiative losses are active, and if epsilon_rad is not used
    if radiativeLosses and not fixed_epsilon:
        
        
        dE2sh_el = (Gamma-1) * dM2 * c**2 * epsilone
        if gamma_c_w > gamma_min:
            N0_inj = dM2 /mp / ((gamma_c_w**(-p-1)-gamma_min**(-p+1)) / (1-p) - (gamma_c_w*gamma_max**(-p) - gamma_c_w**(-p+1))/p)
        else:
            N0_inj = dM2 /mp / ((gamma_max**(1-p)-gamma_min**(1-p))/(1-p) + gamma_min**(2-p)/gamma_c_w - gamma_min**(1-p))
        

        
        
        
        gamma_c_mean = np.sqrt(gamma_min * gamma_max)
        if gamma_c_mean > gamma_max:
            gamma_c_mean = gamma_max

        dErad2 = 0.
        if cooltime_integ:
        
            if gamma_min < gamma_max and dE2sh_el > 0.:

                dErad2 += N0_inj*me*c**2*((gamma_c_mean**(2-p)-gamma_max**(2-p))/(3-p)/(p-2) - (gamma_min**(3-p)*(gamma_c_mean**-1 - gamma_max**-1))/(3-p) + (gamma_c_mean**(2-p) - gamma_max**(2-p))/(2-p)**2   -  gamma_max**(2-p)*np.log(gamma_c_mean/gamma_max) / (2-p) )


                
                

        if remix_radloss and gamma_c_w > 0:
            if gamma_min > gamma_max:
                dErad2 += dE2sh_el
            else:
                if gamma_c_w < gamma_max: ### No losses when gamma_c_w > gamma_max
                    dErad2 += N0_inj*me*c**2 * (gamma_min**(2-p)/(p-2)  -  (gamma_c_w*gamma_max**(1-p)-gamma_c_w**(2-p))/(1-p)  -  (gamma_c_w**(2-p)-gamma_min**(2-p))/(2-p))

        if continuous_radloss and gamma_c_w > 0:
            if gamma_c_w > gamma_min: ### Slow cooling
                N0_total = M2/mp*(p-1)*gamma_min**(p-1)
                #N0_total = Eint2 / c**2 /  me * (p-2) * gamma_min**(p-2) * epsilone
                if gamma_c_w > gamma_max:
                    dErad2 += N0_total*sigmaT*B**2*c/6/pi * (gamma_max**(3-p)-gamma_min**(3-p))/(3-p) * deltaR*one_over_beta/c/Gamma
                else:
                    dErad2 += N0_total*sigmaT*B**2*c/6/pi * ((gamma_c_w**(3-p)-gamma_min**(3-p))/(3-p)  +  (gamma_c_w*gamma_max**(2-p) - gamma_c_w**(3-p))/(2-p)) * deltaR*one_over_beta/c/Gamma
            else: ### Fast cooling
                N0_total_RS = M2 / mp * (gamma_min**(1-p)*(gamma_c_w**-1 - gamma_min**-1) + gamma_min**(-p)/p)
                dErad3 += N0_total * sigmaT*c*BRS**2/6/pi * (gamma_min**(1-p)*(gamma_min-gamma_c_w) + (gamma_max**(2-p)-gamma_min**(2-p))/(2-p)) * deltaR*one_over_beta/c/Gamma



        if dE2sh_el == 0:
            epsilon_rad = 0.
            dErad2 = 0.
        else:
            if dErad2 < 0:
                if print_error:
                    print 'dErad2 < 0 !!!'
                raise RuntimeError()
            epsilon_rad = dErad2 / dE2sh_el

#            if epsilon_rad > 1:
#                print epsilon_rad
                

        if np.isnan(epsilon_rad):
            print np.log10(dE2sh_el)

        if np.isnan(dErad2):
            print 'dErad2 isnan!'
            print Gamma
            print Gamma0
            print gamma_c_mean
            print gamma_max
            print gamma_min
        if gamma_c_mean > gamma_max:
            print 'gamma_c_mean larger than gamma_max!'
            print gamma_c_mean
            print gamma_max


        

 

    elif fixed_epsilon and radiativeLosses:
        dErad2 = epsilon_rad * epsilone * (Gamma-1) * dM2 * c**2
    else:
        dErad2 = 0.
        epsilon_rad = 0.


    if reverseShock and not shutOff:

        if V3>0:
            term1 = np.sqrt(8*pi*eB3)
            term2 = np.sqrt(Eint3)
            sqrt_one_over_V3 = np.sqrt(rho3prim / M3)
            BRS = term1 * term2 * sqrt_one_over_V3
        else:
            BRS = 0.
        dE3sh_el = epsilone3*gamma43_minus_one * dM3 * c**2

        if BRS < 1e-70:
            dErad3 = 0.
            epsilon_radRS = 0.
            gamma_min_RS = 1.
            gamma_max_RS = float('inf')
        else:

            gamma_max_RS = (6*pi*qe/sigmaT/BRS)**.5
            gamma_min_RS = (pRS-2)/(pRS-1)*(epsilone3/mue*gamma43_minus_one+1)


            if gamma_c_w_RS > 0:
                gmin_func = lambda gmin: gmin_fzero(gmin , gamma_c_w_RS , gamma_max_RS , pRS , epsilone3 , gamma43_minus_one , mue , False)

                #### Test ####

                if np.isnan(gamma_c_w_RS):
                    print 'isnan gammac!!!'
                    print np.isnan(gamma_max_RS)
                if np.isinf(gamma_max_RS):
                    print 'isinf gamma_max'
                if np.isnan(gamma_max_RS):
                    print 'isnan gammamax!!!'
                    
                
                #### End test ####

                try:
                    gamma_min_RS = optimize.newton(gmin_func,gamma_min_RS)
                except:

                    res = optimize.fsolve(gmin_func , gamma_min_RS)

                    if len(res) > 1:
                        print 'lenght of res (RS) > 1 !!!!!!!!!!!!!'
                        raw_input(res)
                        print 'gamma_min_RS =',gamma_min_RS
                        print 'gamma_c_w_RS =',gamma_c_w_RS
                        print 'gamma_max_RS =',gamma_max_RS
                        x_plot = np.linspace(1,1e6,1e6)
                        y_plot = np.zeros(1e6)
                        for ip,i_xyplot in enumerate(x_plot):
                            y_plot[ip] = gmin_fzero(i_xyplot , gamma_c_w_RS , gamma_max_RS , pRS , epsilone3 , gamma43_minus_one , mue , False)
                        plt.plot(x_plot , y_plot)
                        plt.show()
                    else:
                        gamma_min_RS = res[0]


            if np.isnan(gamma_min_RS):
                print 'gamma_min_RS isnan !!!'
                print 'gamma_c =',gamma_c_w_RS
                print 'gamma_max =',gamma_max_RS
            
            if gamma_min_RS < 1:
                gamma_min_RS = 1.

            #gamma_min_RS = (pRS-2)/(pRS-1)*(Eint3/M3*epsilone3/mue/c**2 + 1)

                    
                

            
            if gamma_c_w_RS > gamma_min_RS: ### Slow cooling
                if gamma_c_w_RS > gamma_max_RS:
                    N0_inj_RS = dM3 / mp * (pRS-1) * gamma_min_RS**(pRS-1)
                else:
                    N0_inj_RS = dM3 / mp / ((gamma_c_w_RS**(-pRS+1)-gamma_min_RS**(-pRS+1)) / (1-pRS) - (gamma_c_w_RS*gamma_max_RS**(-pRS) - gamma_c_w_RS**(-pRS+1))/pRS)
            else: ### Fast cooling
                N0_inj_RS = dM3 /mp / ((gamma_max_RS**(1-pRS)-gamma_min_RS**(1-pRS))/(1-pRS) + gamma_min_RS**(2-pRS)/gamma_c_w_RS - gamma_min_RS**(1-pRS))

            gamma_c_mean_RS = np.sqrt(gamma_max_RS*gamma_min_RS)
            dErad3 = 0.
            if cooltime_integ:

                if gamma_min_RS < gamma_max_RS: 
                    dErad3 += N0_inj_RS*me*c**2*((gamma_c_mean_RS**(2-pRS)-gamma_max_RS**(2-pRS))/(3-pRS)/(pRS-2) - (gamma_min_RS**(3-p)*(gamma_c_mean_RS**-1 - gamma_max_RS**-1))/(3-pRS) + (gamma_c_mean_RS**(2-pRS) - gamma_max_RS**(2-pRS))/(2-pRS)**2  - gamma_max_RS**(2-p)*np.log(gamma_c_mean_RS/gamma_max_RS)/(2-p)  )

                if dErad3 < 0:
                    dErad3 = 0.
                
            if remix_radloss and gamma_c_w_RS > 0:
                if gamma_min_RS > gamma_max_RS:
                    dErad3 += dE3sh_el
                    print 'gamma_min_RS > gamma_max_RS !!!!!!!!!!!!!!!!'
                else:
                    if gamma_c_w_RS < gamma_max_RS:  ### Only losses when gamma_c_w_RS is smaller than gamma_max_RS
                        dErad3 += N0_inj_RS*me*c**2 * ((gamma_min_RS**(2-pRS))/(pRS-2)  -  (gamma_c_w_RS*gamma_max_RS**(1-pRS)-gamma_c_w_RS**(2-pRS))/(1-pRS)  -  (gamma_c_w_RS**(2-pRS)-gamma_min_RS**(2-pRS))/(2-pRS))
            if dErad3 < 0:
                print 'remix'

            if dErad3 > gamma43_minus_one*dM3*c**2:
                print 'remix larger than shocked energy!'

            if continuous_radloss and gamma_c_w_RS > 0:
                ### Continuous radiation from the entire shocked region. We need to calculate N0 separately, since N0_inj is for the newly shocked particles
                if gamma_min_RS < gamma_max_RS:
                    if gamma_c_w_RS > gamma_min_RS: ### Slow cooling
                        if dM3 > 0:
                            N0_total_RS = M3/mp*(p-1)*gamma_min_RS**(p-1)
                            #N0_total_RS = M3 / me * (pRS-2) * gamma_min_RS**(pRS-2) * (gamma43 - 1)
                            if gamma_c_w_RS > gamma_max_RS:
                                dErad3 += N0_total_RS*sigmaT*BRS**2*c/6/pi * (gamma_max_RS**(3-pRS)-gamma_min_RS**(3-pRS))/(3-pRS) * deltaR*one_over_beta/c/Gamma
                            else:
                                dErad3 += N0_total_RS*sigmaT*BRS**2*c/6/pi * ((gamma_c_w_RS**(3-pRS)-gamma_min_RS**(3-pRS))/(3-pRS)  +  (gamma_c_w_RS*gamma_max_RS**(2-pRS) - gamma_c_w_RS**(3-pRS))/(2-pRS)) * deltaR*one_over_beta/c/Gamma
                        
                        else:
                            print 'dM3 = 0 !!!!!!!!!!!!!!!!!!!!!!!!'
                    else: ### Fast cooling

                        N0_total_RS = M3 / mp * (gamma_min_RS**(1-pRS)*(gamma_c_w_RS**-1 - gamma_min_RS**-1) + gamma_min_RS**(-pRS)/pRS)
                        dErad3 += N0_total_RS * sigmaT*c*BRS**2/6/pi * (gamma_min_RS**(1-pRS)*(gamma_min_RS-gamma_c_w_RS) + (gamma_max_RS**(2-pRS)-gamma_min_RS**(2-pRS))/(2-pRS)) * deltaR*one_over_beta/c/Gamma
                        
            if dErad3 < 0:
                raw_input('dErad3 < 0 !!!!!!!!!!!!!!')
            



        if dM3 == 0 or dE3sh_el == 0: 
            epsilon_radRS = 0.
        else:




            if gamma43_minus_one < 0: 
                epsilon_radRS = 0
                if print_error:
                    print 'gamma43 < 1 !'
                raise RuntimeError()
            else:
                epsilon_radRS = dErad3 / dE3sh_el
#                if epsilon_radRS > 1:
#                    epsilon_radRS = 1.
#                    dErad3 = epsilone3*(gamma43-1) * dM3 * c**2




            if np.isinf(epsilon_radRS):
                print i
                print 'dErad3 = %s'%dErad3
                print 'dEsh3 = %s'%dE3sh_el

                print gamma43_minus_one


    elif reverseShock:
        dErad3 = 0.
        epsilon_radRS = 0.



    if R > R_ISM: 
        dlnrho1dR = 0.
    else:
        dlnrho1dR = -s/R

                

    f_2 = GammaEff*gammaAdi_minus_one*Eint2/Gamma
    h_2 = GammaEff*gammaAdi_minus_one*Eint2*(dM2dR/M2-dlnrho1dR)


    ### Energy of region 3 is still decreasing through adiabatic expansion

    if reverseShock: 

        #Approach: The region 2 (Nava 2013) is the forward shock, i.e. same as in the case without the reverse shock. Now we are adding terms for the reverse shocks, in the region labeled 3.
        
        if M3>0 and not shutOff:
            ######################################
            ### Reverse Shock Energy Evolution ###
            ######################################
            

            if exponential_outflow: 
                dlnrho4dR = -beta0*(one_over_beta**2-one_over_beta0**2)/(one_over_beta+one_over_beta0) /alpha_of - 2/R
#                dlnrho4dR = -2/R - (beta0 / beta - 1) / alpha_of

            else: ### Constant outflow
                dlnrho4dR = -2/R


            fh_factor3 = GammaEff3*gammaAdi3_minus_one*Eint3
            f_3 = fh_factor3 / (gamma43_minus_one+1) * dgamma43dGamma
            h_3 = fh_factor3 * (dM3dR/M3 - dlnrho4dR)
            dEsh3 = gamma43_minus_one*dM3*c**2


            
        else: 
            dEad3 = 0.
            dEsh3 = 0.
            f_3 = 0.
            h_3 = 0.
            dlnV3dR = 0.
            delta_Eint3 = 0.
            dgamma43dGamma = 0.   
            gamma_min_RS = 1.


            ###
            ### Temporary
            ###

            dlnrho4dR = 0.



        
            
            

                    

        dGammadR_next = -(  (Gamma-1)*(GammaEff+1)*dM2dR*c**2 + (Gamma-Gamma0+GammaEff3*gamma43_minus_one)*dM3dR*c**2 - h_2 - h_3   ) / ((M2+M3)*c**2 + Eint2*dGammaEffdGamma + Eint3*dGammaEff3dGamma + f_2 + f_3)



    else:  #Without reverse shock
                

        dGammadR_next = -((Gamma-1)*(GammaEff+1)*dM2dR*c**2 - h_2) / ((M0+M2)*c**2 + Eint2*dGammaEffdGamma + f_2)

    
    Gamma_next = Gamma + dGammadR_next * deltaR

    ######################################
    ### Forward Shock Energy Evolution ###
    ######################################
            
    dlnV2dR = dM2dR/M2 - dlnrho1dR - dGammadR_next/Gamma
    dEad2dR = -gammaAdi_minus_one*Eint2*dlnV2dR

    dEsh2 = (Gamma-1)*dM2*c**2
    delta_Eint2 = dEad2dR*deltaR + dEsh2 - dErad2   
    Eint2_next = Eint2 + delta_Eint2

    if reverseShock:
        if M3>0 and not shutOff:
            dlnGamma43dR = dGammadR_next * dgamma43dGamma / (gamma43_minus_one+1)
            dlnV3dR = dM3dR / M3 - dlnrho4dR - dlnGamma43dR
            
            dEad3 = -gammaAdi3_minus_one*Eint3*dlnV3dR * deltaR
            if dEad3 > 0:
                print 'dEad3 > 0 !!!'
                print 'R =',R
                raw_input('dEad2 = %s'%(dEad2dR*deltaR))

            delta_Eint3 = dEsh3 + dEad3 - dErad3

    if Eint2_next < 0:
        if print_error: 
            print 'Eint2_next < 0'
        raise RuntimeError()

    Eint3_next = Eint3 + delta_Eint3
    if Eint3_next <= 0.:
        ### Early in the evolution, sudden expansion of the shocked region may give a negative internal energy. In this case, we set it to 0. If later in evolution, we ask adaptive step-size module to take a smaller step
             
        if Eint3 > 0:
            if print_error:
                print 'Eint3_next < 0'
            raise RuntimeError()
        else:
            Eint3_next = 0.


    ##################
    ### Fail-safes ###
    ##################
    
    ### Gamma may increase a very small fraction during the early evolution. This causes trouble in the RS description, so Gamma is forced not to exceed Gamma0
    if Gamma_next > Gamma0: 
        if print_error:
            print 'Constraining Gamma to Gamma0'
#        Gamma_next = Gamma0

    if Gamma_next < 1: 
        if print_error:
            print 'gamma<1'
            print 'step = ',deltaR
            print 'R =',R
            print Eint2
            print Eint3
            print Gamma_next
        raise RuntimeError('Gamma < 1')

    #if Gamma_next / Gamma < 0.95 or Gamma_next / Gamma > 1.05: ### To large change, need to reduce stepsize
        #if print_error:
            #print 'gamma changing fast'
        #raise RuntimeError()


    ### Returning Gamma when we are calculating midsteps

    if not full_step: ### Boolean vaiable full_step is True if full return is desired. Otherwise we return only Gamma
        if reverseShock:
            return Gamma_next , M2_next , M3_next , Eint2_next , Eint3_next
        else:
            return Gamma_next , M2_next , Eint2_next


    #############################################
    ### Evolving other variables to next step ###
    #############################################
       
    gammaAdi_minus_one_next = (1+1/Gamma_next) / 3.

    GammaEff_next = (gammaAdi_minus_one_next*Gamma_next**2+Gamma_next**2 - gammaAdi_minus_one_next+2)/Gamma_next

    dGammaEffdGamma_next = 4/3. +1/Gamma**2/3 + 2/Gamma**3/3  ### Note from 2016-06-21

    beta_next = (1-1/Gamma_next**2)**(1/2.)
    if reverseShock and not shutOff:
        ### Using Taylor expansion to calculate (1 - beta * beta0) to avoid numerical errors
        bbm1 = 0.5/Gamma_next**2 + 0.5/Gamma0**2 + Gamma_next**(-4)/8 + Gamma0**(-4)/8 - 0.25/(Gamma0**2)/(Gamma_next**2)
        G0Gi = Gamma0*Gamma_next
            
        if Gamma_next >= Gamma0:
            gamma43_minus_one_next = 0.
        else:
            if Gamma > 0.9*Gamma0:
                ### Using McLaurin expansion to calculate gamma43 - 1
                beta43 = (Gamma0**2/2+Gamma0**2/8/Gamma**2 - Gamma**2/2-Gamma**2/8/Gamma0**2) / (Gamma0**2/2 + Gamma0**2/8/Gamma**2 + Gamma**2/2 + Gamma**2/8/Gamma0**2 - 0.25)
                gamma43_minus_one_next = 0.5*beta43**2 + 3/8.*beta43**4 + 5/16.*beta43**6

                """
                ### Error estimation ###
                Gamma_return_func = lambda Gamma: Gamma/2/Gamma0 + Gamma/8/Gamma0**3 + Gamma0/2/Gamma - 0.25/Gamma0/Gamma + Gamma0/8/Gamma**3 - gamma43_minus_one_next - 1
                print 'gamma43-1  =  ',gamma43_minus_one_next
                print 'Gamma_next =  ',Gamma_next
                print 'Gamma0 =      ',Gamma0
                print 'beta43 =      ',beta43
                raw_input('Gamma_return = %s'%(optimize.newton(Gamma_return_func , Gamma0)))
                """
            else:
                gamma43_minus_one_next = G0Gi * bbm1 - 1

        if gamma43_minus_one_next < 0:
            gamma43_minus_one_next = 0.0
            if print_error:
                print 'gamma43 is smaller than 1'
            raise RuntimeError()

        gammaAdi3_minus_one_next = (1+1/(gamma43_minus_one_next+1)) / 3.    #Is it supposed to be the same as for the FS?

        GammaEff3_next = (gammaAdi3_minus_one_next*Gamma_next**2+Gamma_next**2 - gammaAdi3_minus_one_next+2)/Gamma_next
        dgamma43dGamma_next = 0.5/Gamma0 + 0.5/Gamma_next**2*(0.5/Gamma0 - Gamma0 - 3*Gamma0/8/Gamma_next**2) ### Note from 2016-06-21. Using asymptotic approximation of beta
        dGammaEff3dGamma_next = (gammaAdi3_minus_one_next+1)*(1+Gamma_next**-2) - dgamma43dGamma_next/3/(gamma43_minus_one_next+1)**2*(Gamma_next-1/Gamma_next) - Gamma_next**-2  ### Note from 2016-06-21
    else:
        gamma43_minus_one_next = 0.

        dGammaEff3dGamma_next = 0.
        GammaEff3_next = 0.
        BRS = 0.
        shutOff_next = True


              
        


    if reverseShock: 

        return Gamma_next , theta_next , gammaAdi_minus_one_next , gamma43_minus_one_next , beta_next , dM2 , dM3 , M2_next , M3_next , dGammaEffdGamma_next , dGammaEff3dGamma_next , GammaEff_next , GammaEff3_next , Eint2_next , Eint3_next , B , BRS , rho4 , dGammadR_next  , shutOff_next , epsilon_rad , epsilon_radRS , gamma_min , gamma_min_RS , dgamma43dGamma_next  , gammaAdi3_minus_one_next 

    else: 
        
        return Gamma_next , theta_next , gammaAdi_minus_one_next , beta_next , dM2 , M2_next , dGammaEffdGamma_next , GammaEff_next , Eint2_next , B , dGammadR_next , epsilon_rad , gamma_min

