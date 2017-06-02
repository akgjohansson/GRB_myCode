class Dynamics:
    def __init__(self,R0,ModVar,UseOp,NatCon,tobsEnd,tol): 
        import numpy as np
        self.Gamma,self.tburst,self.tobs,self.tcomoving,self.theta,self.m,self.M3,self.Eint2,self.Eint3,self.B,self.BRS,self.rho,self.rho4,self.R,self.gamma43_minus_one,self.gamma_min,self.gamma_min_RS,self.gammac,self.gammacRS,self.RS_elements_upper = dyn(R0,ModVar,UseOp,tobsEnd,tol)

        self.nprim = 4 * self.Gamma * self.rho / NatCon.mp
        #self.gammaAd = (4 + 1/self.Gamma) / 3
        self.beta = np.sqrt(1-self.Gamma**-2)
        self.thickness_FS = self.m/ (8*np.pi*(1.-np.cos(self.theta)) * self.Gamma**2*self.rho*self.R**2)
        if UseOp.reverseShock:
            self.rho3prim = 4 * self.Gamma[:self.RS_elements_upper] * self.rho4
            self.nprim3 = self.rho3prim / NatCon.mp
            self.thickness_RS = self.M3[:self.RS_elements_upper]/ (8*np.pi*(1.-np.cos(self.theta[:self.RS_elements_upper])) * self.Gamma[:self.RS_elements_upper]**2*self.rho4[:self.RS_elements_upper]*self.R[:self.RS_elements_upper]**2)

#This module considers the Nava et al. (2013) paper. 
def dyn(R0,ModVar,UseOp,tobsEnd,tol):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import ode
    from scipy import optimize
    from gammamin_fzero import minimize_gammamin
    from nextstep_ode import nextstep_func
    if UseOp.save_params:
        import os
    if tol < 1e-5:
        print_error = True
    else:
        print_error = False

    c,pi,mp,me,kB = 2.9979e10,np.pi, 1.6726e-24, 9.1094e-28 , 1.38065e-16
    qe = 4.803204e-10
    sigma_B = 5.6704e-5  ### Stephan Boltzann constant in cgs                   
    rad_const = 4 * sigma_B / c   #### Radiation constant                 
#Initial conditions
    #ModVar.M0 = ModVar.E0*c**-2
    if UseOp.reverseShock: grandNan = [float('nan')] * 19
    else: grandNan = [float('nan')] * 12

    useSpread = True   #If you want to turn spread off, turn it off on this line!
    
    

    

    UseOp.remix_radloss = True
    UseOp.continuous_radloss = True

    #################
    ### Constants ###
    #################

    mup,mue = 1.,me / mp
    sigmaT = 6.6524e-25
    numeric_Eint = True   ### numeric or analytical approach to solve Eint?



    ##############################
    ### Runge-Kutta 853 method ###
    ##############################

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

    ### Allocating space and setting initial values
    Rstop = 1e23
    nsteps = int(np.ceil(100*np.log10(Rstop/R0)))

    beta0 = 1 - 0.5/ModVar.Gamma0**2 - 1/8./ModVar.Gamma0**4
    beta0c = beta0*c
    
    M20 = 2*pi*R0**3*(1-np.cos(ModVar.theta0)) * ModVar.A0 * (mp+me) * R0**-ModVar.s /3
    tburst0,tcomoving0 = R0/beta0c  ,  R0/beta0c/ModVar.Gamma0
    tobs,rho,B = np.zeros(nsteps),np.zeros(nsteps),np.zeros(nsteps)
    Eint20 = M20*(ModVar.Gamma0-1)*c**2
    gamma_c_w = np.zeros(nsteps)
    gamma_c_w[0] = 6*me*c/(ModVar.A0*R0**(-ModVar.s)*(mp+me))/sigmaT/8/tcomoving0/ModVar.eB/Eint20*M20             ### Mass weighted cooling Lorentz factor

    if UseOp.reverseShock:
        gamma_c_w_RS = np.ones(nsteps)
        gamma_c_w_RS *= float('inf')

        rho4,BRS,gamma43_minus_one = np.zeros(nsteps),np.zeros(nsteps),np.zeros(nsteps)

        alpha_of = ModVar.tprompt*beta0c
        shutOff = False
    if ModVar.theta0 < 1e-3:
        one_minus_costheta0 = ModVar.theta0**2/2 - ModVar.theta0**4/24 + ModVar.theta0**6/120
    else:
        one_minus_costheta0 = 1 - np.cos(ModVar.theta0)

    initial_values = [tburst0,tcomoving0,ModVar.Gamma0,(ModVar.Gamma0-1)*M20/ModVar.M0,0.,ModVar.theta0,0.,0.,0.,0.,0.,0.,M20/ModVar.M0,0.,0.]
    #initial_values = [tburst0,tcomoving0,ModVar.Gamma0,(ModVar.Gamma0-1)*M20/ModVar.M0,0.,ModVar.theta0,0.,0.,0.,0.,0.,0.,M20/ModVar.M0,0.,0.] ### Logarithmic
    """
    ### An initial very small step, to initialize reverse shock quantities
    infmal_R0 = R0*1e-9

    initial_values += nextstep_func(R0 , initial_values , s , ModVar.Gamma0  , z , ModVar.theta0 , ModVar.tprompt , ModVar.M0 , ModVar.A0/ModVar.M0 , ModVar.R_ISM , [UseOp.radiativeLosses,remix_radloss,continuous_radloss,fixed_epsilon,UseOp.reverseShock,UseOp.exponential_outflow] , [ModVar.epsilone,ModVar.eB,p,ModVar.epsilone3,ModVar.eB3,pRS],float('inf'),float('inf')) * infmal_R0
    

    R0 = R0 + infmal_R0
    """
    odeinstance = ode(nextstep_func)
    #firststep = R0/10000
    odeinstance.set_integrator("dop853", rtol=1e-5, nsteps=100000, first_step=R0 * 1e-9)
    if UseOp.reverseShock:
        odeinstance.set_f_params(ModVar , UseOp , float('inf') , float('inf') , False)
    else:
        odeinstance.set_f_params(ModVar , UseOp , float('inf'))

    odeinstance.set_initial_value(initial_values, R0)

    
    

    R = np.zeros(nsteps)
    out = np.zeros([nsteps,len(initial_values)])
    for i in range(nsteps):

        ###################
        ### Integrating ###
        ###################

        nextR = R0*(Rstop/R0)**(float(i+1)/float(nsteps))

        
        #try:
        out[i] = odeinstance.integrate(nextR)

        """
        except:
            while firststep > 1:
                if i == 0:
                    firststep /= 10
                    print 'reducing first step to %s'%(firststep)
                    try:
                        odeinstance.set_integrator("dop853", rtol=1e-3, nsteps=100000, first_step=firststep)
                        out[i] = odeinstance.integrate(nextR)
                        break
                    except:
                         firststep /= 10
                   
                else:
                    print '\n\n\n!!!\n\n\n'
                    raise NameError('i is not 0 and nextstep_ode.py still crashed')
        """
        R[i] = np.copy(odeinstance.t)
        tobs[i] = (1+ModVar.z) * (out[i,0] - R[i] / c)
        beta = np.sqrt(1-out[i,2]**-2)


        ######################################
        ### Magnetic fields and densities  ###
        ######################################

        rho1_fac = ModVar.A0 * (mp+me)
        if ModVar.s == 0: ### CM
            rho[i] = ModVar.A0 * (mp+me)
        else:
            if R[i] < ModVar.R_ISM:
                rho[i] = ModVar.A0 * R[i]**(-ModVar.s) * (mp+me)
            else:
                rho[i] = ModVar.A0 * R[i]**(-ModVar.s) * (mp+me) 
    
        rho2 = 4*rho1_fac * out[i,2] * R[i]**(-ModVar.s)        
        B_fac = np.sqrt(8*pi*ModVar.eB)
        B[i] = B_fac * np.sqrt(out[i,3]) / np.sqrt(out[i,12]) * np.sqrt(rho2) * c  ### Factor c is because Eint is on the form E/ModVar.M0/c**2 and M on the form M/ModVar.M0
        gamma_max = (6*pi*qe/sigmaT/B[i])**.5
        if UseOp.reverseShock:
            if out[i,4] < 0:   ### Negative energy in the RS
                out[i,4] = 0.
            if beta < 0.99:
                gamma43_minus_one[i] = out[i,2] * ModVar.Gamma0 * (1-beta*beta0)
            else:
                gamma43_minus_one[i] = out[i,2]*ModVar.Gamma0 * (1/ModVar.Gamma0**2 + 1/out[i,2]**2 - 1/out[i,2]**2/ModVar.Gamma0**2) / (1+beta*beta0) - 1
            if out[i,13]>0:
                rho4_fac_1 = ModVar.M0 / 2 / alpha_of / np.pi / one_minus_costheta0
                rho4_fac = rho4_fac_1  / R[i]**2
                rho4[i] = rho4_fac * np.exp(-out[i,14]/alpha_of) 
                if rho4[i] < 1e-60:
                    rho4[i] = 0. ### Avoiding numerical overflows
                    BRS[i] = 0.
                    shutOff = True
                    odeinstance.set_f_params(ModVar , UseOp , gamma_c_w[i] , gamma_c_w_RS[i], shutOff)

                else:
                    rho3 = 4*out[i,2]*rho4[i]
                    BRS_fac = np.sqrt(8*pi*ModVar.eB3)
                    BRS[i] = BRS_fac * np.sqrt(out[i,4]) / np.sqrt(out[i,13]) * np.sqrt(rho3) * c ### Factor c is because Eint is on the form E/ModVar.M0/c**2 and M on the form M/ModVar.M0
                    if np.isnan(BRS[i]):
                        print 'E3 =',out[i,4]
                        print 'M3 =',out[i,13]
                        print 'rho3 =',rho3
                    if BRS[i] == 0:
                        BRS[i] = 1e-70
                    gamma_max_RS = np.sqrt(6*pi*qe/sigmaT/BRS[i])
                    

                            
        #############################
        ### Calculating gamma_c_w ###
        #############################
        if i > 0:
            dM2 = (out[1:i+1,12]-out[:i,12])
            gamma_c_w_fac = 6*pi*me*c/sigmaT 
            rho2_array = 4 * rho1_fac * out[:i,2] * R[:i] ** (-ModVar.s)
            B_array = B_fac * np.sqrt(out[:i,3]) / np.sqrt(out[:i,12]) * np.sqrt(rho2_array) * c
            gamma_c_w_array = gamma_c_w_fac / B_array**2 / out[:i,1]
            gamma_c_w[i] = np.sum(dM2 * gamma_c_w_array) / (out[i-1,12]-M20/ModVar.M0)
            if gamma_c_w[i] > gamma_max:
                gamma_c_w[i] = np.copy(gamma_max)
            if UseOp.reverseShock and BRS[i] > 0:
                dM3 = (out[1:i+1,13] - out[:i,13])
                rho3_array = 4 * out[:i,2] * rho4_fac_1 * np.exp(-out[:i,14]/alpha_of) * R[:i]**-2
                BRS_array = BRS[:i]#BRS_fac * np.sqrt(out[:i,4]) / np.sqrt(out[:i,13]) * np.sqrt(rho3_array) * c

                
                gamma_c_w_RS_array = gamma_c_w_fac / BRS_array**2 / out[:i,1]
                gamma_c_w_RS[i] = np.sum(dM3 * gamma_c_w_RS_array) / out[i-1,13]
                """
                if gamma_c_w_RS[i] > gamma_max_RS:
                    gamma_c_w_RS[i] = np.copy(gamma_max_RS)
                """

 
        ############################
        ### Calulating gamma_min ###
        ############################

        """
            
        #gamma_min[i] = (p-2)/(p-1)*(ModVar.epsilone/mue*(out[i,2]-1)+1)
        

        gamma_min[i] = minimize_gammamin(gamma_c_w[i] , gamma_max , p , ModVar.epsilone , out[i,2] - 1,mue)
        gamma_min_inj = minimize_gammamin(gamma_max , gamma_max , p ,ModVar.epsilone , out[i,2] - 1,mue)
        
        
        if UseOp.reverseShock:
            
            if (out[i,13]>0) and (BRS[i]>0) and (gamma_c_w_RS[i] > 1) and (i > 0): 
                gamma_min_RS[i] = minimize_gammamin(gamma_c_w_RS[i] , gamma_max_RS , pRS , ModVar.epsilone3 , gamma43_minus_one[i] , mue)
                gamma_min_RS_inj = minimize_gammamin(gamma_max_RS , gamma_max_RS , pRS , ModVar.epsilone3 , gamma43_minus_one[i] , mue)
            else:
                gamma_min_RS[i] = (pRS-2)/(pRS-1)*(ModVar.epsilone3/mue*gamma43_minus_one[i]+1)
                gamma_min_RS_inj = np.copy(gamma_min_RS[i])
            odeinstance.set_f_params(s , ModVar.Gamma0  , z , ModVar.theta0 , ModVar.tprompt , ModVar.M0 , ModVar.A0/ModVar.M0 , ModVar.R_ISM , [UseOp.radiativeLosses,remix_radloss,continuous_radloss,fixed_epsilon,UseOp.reverseShock,UseOp.exponential_outflow] , [ModVar.epsilone,ModVar.eB,p,ModVar.epsilone3,ModVar.eB3,pRS] , gamma_c_w[i] , gamma_min_inj , gamma_min[i] , gamma_c_w_RS[i], gamma_min_RS_inj , gamma_min_RS[i] , shutOff)
        else:
            odeinstance.set_f_params(s , ModVar.Gamma0  , z , ModVar.theta0 , ModVar.tprompt , ModVar.M0 , ModVar.A0/ModVar.M0 , ModVar.R_ISM , [UseOp.radiativeLosses,remix_radloss,continuous_radloss,fixed_epsilon,UseOp.reverseShock,UseOp.exponential_outflow] , [ModVar.epsilone,ModVar.eB,p,ModVar.epsilone3,ModVar.eB3,pRS] , gamma_c_w[i] , gamma_min_inj , gamma_min[i])
        if ( not odeinstance.successful()):
            raise NameError('stepsize')
        if ( tobs[i-3] > tobsEnd):#*365.35*10 ):
            break #Stop at 10 years
        """

    out[:,3] *= ModVar.M0*c**2
    out[:,4] *= ModVar.M0*c**2
    out[:,6] *= ModVar.M0*c**2
    out[:,7] *= ModVar.M0*c**2
    out[:,8] *= ModVar.M0*c**2
    out[:,9] *= ModVar.M0*c**2
    out[:,10] *= ModVar.M0*c**2
    out[:,11] *= ModVar.M0*c**2
    out[:,12] *= ModVar.M0
    out[:,13] *= ModVar.M0
    
    """
    plt.plot(R[:i+1] , out[:i+1,2])
    plt.ylabel(r'$\Gamma$')
    plt.loglog()
    plt.show()

    plt.plot(R[:i] , (out[1:i+1,13]-out[:i,13]) / (R[1:i+1]-R[:i]))
    plt.loglog()
    plt.ylabel(r'$\frac{dM_3}{dR}$')
    plt.show()

    plt.plot(R[:i+1] , out[:i+1,12])
    plt.ylabel(r'$M_2$')
    plt.loglog()
    plt.show()

    plt.plot(R[:i+1] , out[:i+1,13])
    plt.ylabel(r'$M_3$')
    plt.loglog()
    plt.show()
    plt.plot(R[:i+1] , out[:i+1,3])
    plt.ylabel(r'$E_{\rm int,2}$')
    plt.loglog()
    plt.show()
    plt.plot(R[:i+1] , out[:i+1,4])
    plt.ylabel(r'$E_{\rm int,3}$')
    plt.loglog()
    plt.show()
    plt.plot(R[:i+1] , out[:i+1,14])
    plt.ylabel(r'$\Delta R_4$')
    plt.loglog()
    plt.show()
    """
    ### Calculating gamma_min ###
    gamma_min = (ModVar.p-2)/(ModVar.p-1)*(1+mp/me*ModVar.epsilone*(out[:i+1,2]-1))
    gamma_min[np.where(gamma_min<1)] = 1.
    if UseOp.reverseShock:
        ### Calculating gamma_min_RS ###
        gamma_min_RS = (ModVar.pRS-2)/(ModVar.pRS-1)*(1+mp/me*ModVar.epsilone3*gamma43_minus_one[:i+1])

    if UseOp.save_params:

        os.system('rm Parameters/*')
        np.savetxt('Parameters/tobs.txt',tobs[:i+1])
        np.savetxt('Parameters/tburst.txt',out[:i+1,0])
        np.savetxt('Parameters/tcomoving.txt',out[:i+1,1])
        np.savetxt('Parameters/Gamma.txt',out[:i+1,2])
        np.savetxt('Parameters/R.txt',R[:i+1])
        #np.savetxt('Parameters/ModVar.epsilon_rad.txt',ModVar.epsilon_rad[:i+1])
        np.savetxt('Parameters/dM2dR.txt',(out[1:i+1,12]-out[:i,12])/(R[1:i+1] - R[:i]))
        np.savetxt('Parameters/M2.txt',out[:i+1,12])
        np.savetxt('Parameters/Eint2.txt',out[:i+1,3])
        np.savetxt('Parameters/Esh2.txt',out[:i+1,8])
        np.savetxt('Parameters/Ead2.txt',out[:i+1,10])
        np.savetxt('Parameters/Erad2.txt',out[:i+1,6])
        np.savetxt('Parameters/gamma_c_w.txt',gamma_c_w[:i+1])
        np.savetxt('Parameters/gamma_min.txt',gamma_min)
        np.savetxt('Parameters/B.txt' , B[:i+1])
        np.savetxt('Parameters/theta.txt' , out[:i+1,5])
        np.savetxt('Parameters/rho.txt' , rho[:i+1])
        if UseOp.reverseShock:

            np.savetxt('Parameters/dM3dR.txt',(out[1:i+1,13]-out[:i,13])/(R[1:i+1] - R[:i]))
            np.savetxt('Parameters/M3.txt',out[:i+1,13])
            #np.savetxt('Parameters/epsilon_rad_RS.txt',ModVar.epsilon_rad_RS[:i+1])
            np.savetxt('Parameters/Eint3.txt',out[:i+1,4])
            np.savetxt('Parameters/Esh3.txt',out[:i+1,9])
            np.savetxt('Parameters/Ead3.txt',out[:i+1,11])
            np.savetxt('Parameters/Erad3.txt',out[:i+1,7])
            np.savetxt('Parameters/Gamma43_minus_one.txt',gamma43_minus_one[:i+1])
            np.savetxt('Parameters/gamma_c_w_RS.txt',gamma_c_w_RS[:i+1])
            np.savetxt('Parameters/gamma_min_RS.txt',gamma_min_RS)
            np.savetxt('Parameters/BRS.txt',BRS[:i+1])
            np.savetxt('Parameters/rho4.txt',rho4[:i+1])


    if UseOp.reverseShock:
        ### RS_elements_upper reduces RS arrays to where there is an RS component
        RS_elements_upper = np.count_nonzero(rho4)

        return out[:i+1,2],out[:i+1,0],tobs[:i+1],out[:i+1,1],out[:i+1,5],out[:i+1,12],out[:RS_elements_upper,13],out[:i+1,3],out[:RS_elements_upper,4],B[:i+1],BRS[:RS_elements_upper],rho[:i+1],rho4[:RS_elements_upper],R[:i+1],gamma43_minus_one[:RS_elements_upper],gamma_min,gamma_min_RS,gamma_c_w[:i+1],gamma_c_w_RS[:RS_elements_upper] , RS_elements_upper
    else:
        return out[:i+1,2],out[:i+1,0],tobs[:i+1],out[:i+1,1],out[:i+1,5],out[:i+1,12],None,out[:i+1,3],None,B[:i+1],None,rho[:i+1],None,R[:i+1],None,gamma_min,None,gamma_c_w[:i+1],None,None
