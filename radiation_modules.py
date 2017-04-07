class rad_var:
    def __init__(self,Dyn,ModVar,UseOp,NatCon,RadCon,RadConRS=None):
        import numpy as np
        self.PmaxF =  RadCon.phipF*2.234*NatCon.qe**3*Dyn.nprim*Dyn.B/NatCon.me/NatCon.c**2  *  Dyn.Gamma
        self.PmaxS = RadCon.phipS*11.17*(ModVar.p-1)*NatCon.qe**3*Dyn.nprim*Dyn.B/(3*ModVar.p-1)/NatCon.me/(NatCon.c**2)  *  Dyn.Gamma


        gamToNuFactor = NatCon.qe * Dyn.B / 2 / np.pi / NatCon.me / NatCon.c
        self.num = Dyn.gamma_min**2 * gamToNuFactor
        self.nuc = Dyn.gammac**2 * gamToNuFactor
        
        if UseOp.reverseShock:
            self.PmaxF_RS = RadConRS.phipF*2.234*NatCon.qe**3*Dyn.nprim3*Dyn.BRS/NatCon.me/NatCon.c**2   *  Dyn.Gamma[:Dyn.RS_elements_upper]   ### The last Gamma is to transform from comoving frame to observer's frame
            self.PmaxS_RS = RadConRS.phipS*11.17*(ModVar.pRS-1)*NatCon.qe**3*Dyn.nprim3*Dyn.BRS/(3*ModVar.pRS-1)/NatCon.me/(NatCon.c**2)  *  Dyn.Gamma[:Dyn.RS_elements_upper]


            gamToNuFactorRS = NatCon.qe * Dyn.BRS / 2 / np.pi / NatCon.c / NatCon.me
            where_nucRS_finite = np.isfinite(Dyn.gammacRS)
            where_nucRS_inf = np.isinf(Dyn.gammacRS)
            self.nucRS = np.zeros(Dyn.RS_elements_upper)
            self.nucRS[where_nucRS_finite] = Dyn.gammacRS[where_nucRS_finite]**2 * gamToNuFactorRS[where_nucRS_finite]
            self.nucRS[where_nucRS_inf] = float('inf')
            self.numRS = Dyn.gamma_min_RS[:Dyn.RS_elements_upper]**2 * gamToNuFactorRS
            self.nucRS[where_nucRS_finite] = Dyn.gammacRS[where_nucRS_finite]**2 * gamToNuFactorRS[where_nucRS_finite]
        

def radiation_function(Dyn , Rad , UseOp , ModVar , nu , Phi , elements , kappas , RS , intermediate,InterWeights=None , last_index=None , first_index=None):
    from EATS_func import eats_function
    import numpy as np

    if intermediate:
        if RS:
            where_slow_cooling = np.where(Rad.numRS[elements] <= Rad.nucRS[elements])[0]

            rhoprim = 4 * Dyn.rho4[elements] * Dyn.Gamma[elements]

            where_fast_cooling = np.where(Rad.numRS[elements] > Rad.nucRS[elements])[0]

            P_fac = 1e23 * Dyn.R[elements]**2 * Dyn.thickness_RS[elements] / (Dyn.Gamma[elements] ** 3 * (1-Dyn.beta[elements]*np.cos(Phi))**3)

        else:
            where_slow_cooling = np.where(Rad.num[elements] <= Rad.nuc[elements])[0]

            rhoprim = 4 * Dyn.rho[elements] * Dyn.Gamma[elements]

            where_fast_cooling = np.where(Rad.num[elements] > Rad.nuc[elements])[0]


            P_fac = 1e23 * Dyn.R[elements]**2 * Dyn.thickness_FS[elements] / (Dyn.Gamma[elements] ** 3 * (1-Dyn.beta[elements]*np.cos(Phi))**3)
        

        return P_fac * eats_function(ModVar , UseOp , Rad , nu[1:-1] , elements , where_slow_cooling , where_fast_cooling , kappas)
    else:
        regions = ['edge' , 'front']
        out = np.zeros(2)
        for i , region in enumerate(regions):
            if RS:
                num = InterWeights.interpolator(Rad.numRS[last_index],Rad.numRS[last_index+1],region)
                nuc = InterWeights.interpolator(Rad.nucRS[last_index],Rad.nucRS[last_index+1],region)
                slow_cooling = numRS <= nucRS
                Gamma = InterWeights.interpolator(Dyn.Gamma[last_index],Dyn.Gamma[last_index+1],region)
                rho = InterWeights.interpolator(Dyn.rho4[last_index],Dyn.rho4[last_index+1],region)
                rhoprim = 4*InterWeights.interpolator(Dyn.Gamma[last_index]*Dyn.rho4[last_index],Dyn.Gamma[last_index+1]*Dyn.rho4[last_index+1],region)

                fast_cooling = numRS > nucRS
                thickness = InterWeights.interpolator(Dyn.M3[last_index]/ ((1.-np.cos(Dyn.theta[last_index])) * Dyn.Gamma[last_index]**2*Dyn.rho4[last_index]*Dyn.R[last_index]**2)  ,  Dyn.M3[last_index+1]/ ((1.-np.cos(Dyn.theta[last_index+1])) * Dyn.Gamma[last_index+1]**2*Dyn.rho4[last_index+1]*Dyn.R[last_index+1]**2) , region) / 8/np.pi
            
            else:

                num = InterWeights.interpolator(Rad.num[last_index],Rad.num[last_index+1],region)
                nuc = InterWeights.interpolator(Rad.nuc[last_index],Rad.nuc[last_index+1],region)
                slow_cooling = num <= nuc
                Gamma = InterWeights.interpolator(Dyn.Gamma[last_index],Dyn.Gamma[last_index+1],region)
                rho = InterWeights.interpolator(Dyn.rho[last_index],Dyn.rho[last_index+1],region)
                rhoprim = 4*InterWeights.interpolator(Dyn.Gamma[last_index]*Dyn.rho[last_index],Dyn.Gamma[last_index+1]*Dyn.rho[last_index+1],region)

                fast_cooling = num > nuc
                thickness = InterWeights.interpolator(Dyn.m[last_index]/ ((1.-np.cos(Dyn.theta[last_index])) * Dyn.Gamma[last_index]**2*Dyn.rho[last_index]*Dyn.R[last_index]**2)  ,  Dyn.m[last_index+1]/ ((1.-np.cos(Dyn.theta[last_index+1])) * Dyn.Gamma[last_index+1]**2*Dyn.rho[last_index+1]*Dyn.R[last_index+1]**2) , region) / 8/np.pi
            if region == 'edge':
                P_b_fac = 1e23 * thickness * InterWeights.interpolator(Dyn.R[last_index]**2 / Dyn.Gamma[last_index]**3 / (1-Dyn.beta[last_index]*np.cos(Phi[0]))**3  ,  Dyn.R[last_index+1]**2 / Dyn.Gamma[last_index+1]**3 / (1-Dyn.beta[last_index+1]*np.cos(Phi[0]))**3  , region)
                if slow_cooling:
                    out[i] =  P_b_fac * InterWeights.interpolator(eats_function(ModVar , UseOp , Rad , nu[0] , None , [last_index] , [] , kappas , RS , True) , eats_function(ModVar , UseOp , Rad , nu[0] , None , [last_index+1] , [] , kappas , RS , True) , region)
                elif fast_cooling:
                    out[i] =  P_b_fac * InterWeights.interpolator(eats_function(ModVar , UseOp , Rad , nu[0] , None , [] , [last_index] , kappas , RS , True) , eats_function(ModVar , UseOp , Rad , nu[0] , None , [] , [last_index+1] , kappas , RS , True) , region)
                else:
                    raise NameError('variable regions is not properly used!')
            elif region == 'front':
                P_b_fac = 1e23 * thickness * InterWeights.interpolator(Dyn.R[first_index]**2 / Dyn.Gamma[first_index]**3 / (1-Dyn.beta[first_index]*np.cos(Phi[0]))**3  ,  Dyn.R[first_index+1]**2 / Dyn.Gamma[first_index+1]**3 / (1-Dyn.beta[first_index+1]*np.cos(Phi[-1]))**3  , region)
                if slow_cooling:
                    out[i] =  P_b_fac * InterWeights.interpolator(eats_function(ModVar , UseOp , Rad , nu[-1] , None , [first_index] , [] , kappas , RS , True) , eats_function(ModVar , UseOp , Rad , nu[-1] , None , [first_index+1] , [] , kappas , RS , True) , region)
                elif fast_cooling:
                    out[i] =  P_b_fac * InterWeights.interpolator(eats_function(ModVar , UseOp , Rad , nu[-1] , None , [] , [first_index] , kappas , RS , True) , eats_function(ModVar , UseOp , Rad , nu[-1] , None , [] , [first_index+1] , kappas , RS , True) , region)
                else:
                    raise NameError('variable regions is not properly used!')
        return out
                                                        
def self_absorption(Dyn, ModVar , absCon ,  Rad , NatCon , InterWeights , nuPrim , indeces , RS):
    import numpy as np
    from EATS_func import alphanu_func
    ### Class absCon is a rename of selfAbs, to avoid confusion. This is because we send the FS and RS classes separately

    ### We have to divide it into FS and RS

    alpha0Ffactor = 11.7 * absCon.phipF * absCon.XpF**(-3) * NatCon.qe / NatCon.mp


    if RS:
        where_slow_cooling = np.where(Rad.numRS[indeces] <= Rad.nucRS[indeces])
        slow_cooling_front = InterWeights.numRS_front <= InterWeights.nucRS_front
        slow_cooling_edge = InterWeights.numRS_edge <= InterWeights.nucRS_edge
        where_fast_cooling = np.where(Rad.nucRS[indeces] < Rad.numRS[indeces])
        alpha0Sfactor = 7.8 * absCon.phipS * absCon.XpS**(-(4+ModVar.pRS)/2.) * (ModVar.pRS+2)*(ModVar.pRS-1)* NatCon.qe / NatCon.mp / (ModVar.pRS+2/3.)

        rhoPrim = 4 * Dyn.rho4[indeces] * Dyn.Gamma[indeces]
        rhoPrim_front = 4 * InterWeights.rho4_front * InterWeights.Gamma_front
        rhoPrim_edge = 4 * InterWeights.rho4_edge * InterWeights.Gamma_edge

        alpha0F = alpha0Ffactor * rhoPrim / Dyn.BRS[indeces] * Dyn.gammacRS[indeces] ** (-5)
        alpha0F_front = alpha0Ffactor * rhoPrim_front / InterWeights.BRS_front * InterWeights.gammacRS_front ** (-5)
        alpha0F_edge = alpha0Ffactor * rhoPrim_edge / InterWeights.BRS_edge * InterWeights.gammacRS_edge ** (-5)

        alpha0S = alpha0Sfactor * rhoPrim * Dyn.gamma_min_RS[indeces] ** (-5) / Dyn.BRS[indeces] 
        alpha0S_front = alpha0Sfactor * rhoPrim_front * InterWeights.gamma_min_RS_front ** (-5) / InterWeights.BRS_front 
        alpha0S_edge = alpha0Sfactor * rhoPrim_edge * InterWeights.gamma_min_RS_edge ** (-5) / InterWeights.BRS_edge 

        alphanu = np.zeros(len(indeces)+2)
        alphanu[1:-1][where_fast_cooling] = alpha0F[where_fast_cooling] * alphanu_func(nuPrim[where_fast_cooling] , Rad.numRS[where_fast_cooling] , Rad.nucRS[where_fast_cooling], True , ModVar.pRS)
        alphanu[1:-1][where_slow_cooling] = alpha0S[where_slow_cooling] * alphanu_func(nuPrim[where_slow_cooling] , Rad.numRS[where_slow_cooling] , Rad.nucRS[where_slow_cooling], False , ModVar.pRS)
        ### front of EATS
        alphanu[-1] = alpha0F_front * alphanu_func(InterWeights.nuPrim_front , InterWeights.numRS_front , InterWeights.nucRS_front , slow_cooling_front==False , ModVar.pRS)
        alphanu[0] = alpha0F_edge * alphanu_func(InterWeights.nuPrim_edge , InterWeights.numRS_edge , InterWeights.nucRS_edge , slow_cooling_edge==False , ModVar.pRS)
        thickness = np.zeros(len(indeces+2))
        thickness[1:-1] = Dyn.thickness_RS[indeces]
        thickness[0] = InterWeights.thickness_RS_edge
        thickness[-1] = InterWeights.thickness_RS_front


    else:
        where_fast_cooling = np.where(Rad.nuc[indeces] < Rad.num[indeces])
        where_slow_cooling = np.where(Rad.num[indeces] <= Rad.nuc[indeces])
        slow_cooling_edge = InterWeights.num_edge <= InterWeights.nuc_edge
        slow_cooling_front = InterWeights.num_front <= InterWeights.nuc_front

        alpha0Ffactor = 11.7 * absCon.phipF * absCon.XpF**(-3) * NatCon.qe / NatCon.mp
        alpha0Sfactor = 7.8 * absCon.phipS * absCon.XpS**(-(4+ModVar.p)/2.) * (ModVar.p+2)*(ModVar.p-1)* NatCon.qe / NatCon.mp / (ModVar.p+2/3.)

        rhoPrim = 4 * Dyn.rho[indeces] * Dyn.Gamma[indeces]
        rhoPrim_edge = 4 * InterWeights.rho_edge * InterWeights.Gamma_edge
        rhoPrim_front = 4 * InterWeights.rho_front * InterWeights.Gamma_front

        alpha0F = alpha0Ffactor * rhoPrim / Dyn.B[indeces] * Dyn.gammac[indeces] ** (-5)
        alpha0F_front = alpha0Ffactor * rhoPrim_front / InterWeights.B_front * InterWeights.gammac_front ** (-5)
        alpha0F_edge = alpha0Ffactor * rhoPrim_edge / InterWeights.B_edge * InterWeights.gammac_edge ** (-5)

        alpha0S = alpha0Sfactor * rhoPrim * Dyn.gamma_min[indeces] ** (-5) / Dyn.B[indeces] 
        alpha0S_edge = alpha0Sfactor * rhoPrim_edge * InterWeights.gamma_min_edge ** (-5) / InterWeights.B_edge 
        alpha0S_front = alpha0Sfactor * rhoPrim_front * InterWeights.gamma_min_front ** (-5) / InterWeights.B_front 

        alphanu = np.zeros(len(indeces)+2)
        alphanu[1:-1][where_fast_cooling] = alpha0F[where_fast_cooling] * alphanu_func(nuPrim[where_fast_cooling] , Rad.num[where_fast_cooling] , Rad.nuc[where_fast_cooling], True , ModVar.p)
        alphanu[1:-1][where_slow_cooling] = alpha0S[where_slow_cooling] * alphanu_func(nuPrim[where_slow_cooling] , Rad.num[where_slow_cooling] , Rad.nuc[where_slow_cooling], False , ModVar.p)
        ### Front of EATS

        alphanu[0] = alpha0F_edge * alphanu_func(InterWeights.nuPrim_edge , InterWeights.num_edge , InterWeights.nuc_edge, slow_cooling_edge==False , ModVar.p)
        alphanu[-1] = alpha0F_front * alphanu_func(InterWeights.nuPrim_front , InterWeights.num_front , InterWeights.nuc_front, slow_cooling_front==False , ModVar.p)
        thickness = np.zeros(len(indeces+2))
        thickness[1:-1] = Dyn.thickness_FS[indeces]
        thickness[0] = InterWeights.thicknessFS_edge
        thickness[-1] = InterWeights.thickness_FS_front


    return alphanu * thickness / 2


class weights:
    def __init__(self,Dyn,UseOp,Rad,ModVar,NatCon,tobs,tobs_behind,tobs_before,first_index,last_index,onePzFreq):
         ### Weights to interpolate the center point of the EATSurface
        import numpy as np

        self.frontWeight1 = Dyn.tobs[first_index+1] - tobs
        self.frontWeight2 = tobs - Dyn.tobs[first_index]
        self.frontWeight = Dyn.tobs[first_index+1] - Dyn.tobs[first_index]

        ### Weights to interpolate the edge point of the EATSurface
        self.edgeWeight1 = tobs_before - tobs
        if self.edgeWeight1 < 0:
            print 'edgeWeight < 0'
            print self.edgeWeight1
            print tobs_behind
            print tobs_before
            print tburst
            print tobs_before
        self.edgeWeight2 = tobs - tobs_behind
        self.edgeWeight = tobs_before - tobs_behind

        ### Interpolating dynamics values at the edge

        self.R_edge = self.interpolator(Dyn.R[last_index] , Dyn.R[last_index+1] , 'edge')
        self.tburst_edge = self.interpolator(Dyn.tburst[last_index] , Dyn.tburst[last_index+1] , 'edge')
        self.tobs_edge = self.interpolator(Dyn.tobs[last_index] , Dyn.tobs[last_index+1] , 'edge')
        self.Gamma_edge = self.interpolator(Dyn.Gamma[last_index] , Dyn.Gamma[last_index+1] , 'edge')
        self.beta_edge = self.interpolator(Dyn.beta[last_index] , Dyn.beta[last_index+1] , 'edge')
        self.theta_edge = self.interpolator(Dyn.theta[last_index] , Dyn.theta[last_index+1] , 'edge')
        self.rho_edge = self.interpolator(Dyn.rho[first_index] , Dyn.rho[first_index+1] , 'edge')
        self.B_edge = self.interpolator(Dyn.B[last_index] , Dyn.B[last_index+1] , 'edge')
        self.num_edge = self.interpolator(Rad.num[last_index] , Rad.num[last_index+1] , 'edge')
        self.nuc_edge = self.interpolator(Rad.nuc[last_index] , Rad.nuc[last_index+1] , 'edge')
        self.gammac_edge = self.interpolator(Dyn.gammac[last_index] , Dyn.gammac[last_index+1] , 'edge')
        self.gamma_min_edge = self.interpolator(Dyn.gamma_min[last_index] , Dyn.gamma_min[last_index+1] , 'edge')
        self.thickness_FS_edge = self.interpolator(Dyn.thickness_FS[last_index] , Dyn.thickness_FS[last_index+1] , 'edge')
        ### Integration angle Phi
        cosPhi_edge = NatCon.c/self.R_edge*(self.tburst_edge - tobs / (1+ModVar.z))
        if (cosPhi_edge >= 1) or (cosPhi_edge <= -1):
            self.Phi_edge = np.pi/2
        else:
            self.Phi_edge = np.arccos(cosPhi_edge)
            
        self.nuPrim_edge = onePzFreq * self.Gamma_edge * (1-self.beta_edge * np.cos(self.Phi_edge))
        
        if UseOp.reverseShock:
            self.BRS_edge = self.interpolator(Dyn.BRS[last_index] , Dyn.BRS[last_index+1] , 'edge')
            self.numRS_edge = self.interpolator(Rad.numRS[last_index] , Rad.numRS[last_index+1] , 'edge')
            self.nucRS_edge = self.interpolator(Rad.nucRS[last_index] , Rad.nucRS[last_index+1] , 'edge')
            self.rho4_edge = self.interpolator(Dyn.rho4[first_index] , Dyn.rho4[first_index+1] , 'edge')            
            self.gammac_RS_edge = self.interpolator(Dyn.gammac_RS[last_index] , Dyn.gammac_RS[last_index+1] , 'edge')
            self.gamma_min_RS_edge = self.interpolator(Dyn.gamma_min_RS[last_index] , Dyn.gamma_min_RS[last_index+1] , 'edge')
            self.thickness_RS_edge = self.interpolator(Dyn.thickness_RS[last_index] , Dyn.thickness_RS[last_index+1] , 'edge')
        ### Interpolating dynamics values of the LoS part of the EATS
        
        self.R_front = self.interpolator(Dyn.R[first_index] , Dyn.R[first_index+1] , 'front')
        self.tburst_front = self.interpolator(Dyn.tburst[first_index] , Dyn.tburst[first_index+1] , 'front')
        self.tobs_front = self.interpolator(Dyn.tobs[first_index] , Dyn.tobs[first_index+1] , 'front')
        self.Gamma_front = self.interpolator(Dyn.Gamma[first_index] , Dyn.Gamma[first_index+1] , 'front')
        self.beta_front = self.interpolator(Dyn.beta[first_index] , Dyn.beta[first_index+1] , 'front')
        self.theta_front = self.interpolator(Dyn.theta[first_index] , Dyn.theta[first_index+1] , 'front')
        self.thickness_FS_front = self.interpolator(Dyn.thickness_FS[last_index] , Dyn.thickness_FS[last_index+1] , 'front')
        self.B_front = self.interpolator(Dyn.B[first_index] , Dyn.B[first_index+1] , 'front')
        self.rho_front = self.interpolator(Dyn.rho[first_index] , Dyn.rho[first_index+1] , 'front')
        self.num_front = self.interpolator(Rad.num[last_index] , Rad.num[last_index+1] , 'front')
        self.nuc_front = self.interpolator(Rad.nuc[last_index] , Rad.nuc[last_index+1] , 'front')
        self.gammac_front = self.interpolator(Dyn.gammac[last_index] , Dyn.gammac[last_index+1] , 'front')
        self.gamma_min_front = self.interpolator(Dyn.gamma_min[last_index] , Dyn.gamma_min[last_index+1] , 'front')
        self.nuPrim_front = onePzFreq * self.Gamma_front * (1-self.beta_front)
        if UseOp.reverseShock:
            self.rho4_front = self.interpolator(Dyn.rho4[first_index] , Dyn.rho4[first_index+1] , 'front')            
            self.BRS_front = self.interpolator(Dyn.BRS[last_index] , Dyn.BRS[last_index+1] , 'front')
            self.numRS_front = self.interpolator(Rad.numRS[last_index] , Rad.numRS[last_index+1] , 'front')
            self.nucRS_front = self.interpolator(Rad.nucRS[last_index] , Rad.nucRS[last_index+1] , 'front')
            self.gammac_RS_front = self.interpolator(Dyn.gammac_RS[last_index] , Dyn.gammac_RS[last_index+1] , 'front')
            self.gamma_min_RS_front = self.interpolator(Dyn.gamma_min_RS[last_index] , Dyn.gamma_min_RS[last_index+1] , 'front')
            self.thickness_RS_front = self.interpolator(Dyn.thickness_RS[last_index] , Dyn.thickness_RS[last_index+1] , 'front')
    def interpolator(self,lower,upper,region):
        if region == 'front':
            return (lower*self.frontWeight1 + upper*self.frontWeight2) / self.frontWeight
        elif region == 'edge':
            return (lower*self.edgeWeight1 + upper*self.edgeWeight2) / self.edgeWeight
        else:
            print region
            raise NameError('Something wrong! This else should not be entered')

class flux_allocation:
    def __init__(self,UseOp,nflux,Plot_Exceptions,npoints,freq):
        import numpy as np
        for i in range(nflux):
            if not Plot_Exceptions.RS_only:
                exec('self.FFS_%d = np.zeros(npoints)'%i)
            if not Plot_Exceptions.FS_only:
                exec('self.FRS_%d = np.zeros(npoints)'%i)
            exec('self.name_%d = \'%s\''%(i,freq[i]))


class exception_class:
    def __init__(self):
        self.FS_only = False
        self.RS_only = False
        self.thermal_only = False

class plot_details:
    def __init__(self):
        self.colourCycle = ['b','g','r','c','m','y','k']     #Cycle of colours in plotting. Matplotlib standard cycle
        self.scalePlotTime = {'d': 86400. , 'h' : 3600. , 'm' : 60. , 's' : 1. }
        self.scaleFluxAxis = {'mJy' : 1.e3 , 'Jy' : 1. }
