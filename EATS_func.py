class kappa_constants():
    def __init__(self,p):
        kappa1 = 2.37 - 0.3*p
        kappa2 = 14.7 - 8.68*p + 1.4*p**2
        kappa3 = 6.94 - 3.844*p + 0.62*p**2
        kappa4 = 3.5 - 0.2*p

        self.kappa13 = -kappa1/3
        self.kappa12 = kappa1/2
        self.kappa11 = -1/kappa1
        self.kappa2p = kappa2*(p-1)/2
        self.kappa12inv = -1/kappa2
        self.kappa33 = -kappa3/3
        self.kappa3p = kappa3*(p-1)/2
        self.kappa13inv = -1/kappa3
        self.kappa42 = kappa4/2
        self.kappa14 = -1/kappa4        

def eats_function(ModVar , UseOp , Rad , nu , elements ,  where_slow_cooling , where_fast_cooling , kappas , RS=False , single_value=False):
    import numpy as np

    if single_value: ### This option is used if we are interpolating and want a single value only
        ### If we are interpolating, we send either fast_cooling or slow_cooling as a single element array, and the other as an empty array (length 0)
        if len(where_slow_cooling) == 1: ### Slow cooling


            if RS:
                return Rad.PmaxS_RS[where_slow_cooling] * ((nu/Rad.numRS[where_slow_cooling])**(kappas.kappa33) + (nu/Rad.numRS[where_slow_cooling])**(kappas.kappa3p))**(kappas.kappa13inv) * (1+(nu/Rad.nucRS[where_slow_cooling])**(kappas.kappa42))**(kappas.kappa14)
            else:
                return Rad.PmaxS[where_slow_cooling] * ((nu/Rad.num[where_slow_cooling])**(kappas.kappa33) + (nu/Rad.num[where_slow_cooling])**(kappas.kappa3p))**(kappas.kappa13inv) * (1+(nu/Rad.nuc[where_slow_cooling])**(kappas.kappa42))**(kappas.kappa14)
        elif len(where_fast_cooling) == 1: ### Fast cooling
            if RS: 
                return Rad.PmaxF_RS[where_fast_cooling] * ((nu/Rad.nucRS[where_fast_cooling])**(kappas.kappa13) + (nu/Rad.nucRS[where_fast_cooling])**(kappas.kappa12)) ** (kappas.kappa11) * (1+(nu/Rad.numRS[where_fast_cooling])**(kappas.kappa2p))**(kappas.kappa12inv)
            else: ###Fast cooling, FS
                return Rad.PmaxF[where_fast_cooling] * ((nu/Rad.nuc[where_fast_cooling])**(kappas.kappa13) + (nu/Rad.nuc[where_fast_cooling])**(kappas.kappa12)) ** (kappas.kappa11) * (1+(nu/Rad.num[where_fast_cooling])**(kappas.kappa2p))**(kappas.kappa12inv)
        else:
            raise NameError('This else should not be entered!')
        

    slow_cooling = elements[where_slow_cooling]
    fast_cooling = elements[where_fast_cooling]
    P_out = np.zeros(len(slow_cooling) + len(fast_cooling)) ###Through this formulation we may use this function when interpolating as well
    

    if RS:
        try:
            P_out[where_fast_cooling] = Rad.PmaxF_RS[fast_cooling] * ((nu[where_fast_cooling]/Rad.nucRS[fast_cooling])**(kappas.kappa13) + (nu[where_fast_cooling]/Rad.nucRS[fast_cooling])**(kappas.kappa12)) ** (kappas.kappa11) * (1+(nu[where_fast_cooling]/Rad.numRS[fast_cooling])**(kappas.kappa2p))**(kappas.kappa12inv)


            P_out[where_slow_cooling] = Rad.PmaxS_RS[slow_cooling] * ((nu[where_slow_cooling]/Rad.numRS[slow_cooling])**(kappas.kappa33) + (nu[where_slow_cooling]/Rad.numRS[slow_cooling])**(kappas.kappa3p))**(kappas.kappa13inv) * (1+(nu[where_slow_cooling]/Rad.nucRS[slow_cooling])**(kappas.kappa42))**(kappas.kappa14)
        except:
            print where_fast_cooling
            print len(where_fast_cooling)
            print len(where_slow_cooling)
            print np.sum(where_fast_cooling)
            print len(P_out)
            print len(nu)
            raw_input("holding")
            raise NameError("slutar")

    else:

        P_out[where_fast_cooling] = Rad.PmaxF[fast_cooling] * ((nu[where_fast_cooling]/Rad.nuc[fast_cooling])**(kappas.kappa13) + (nu[where_fast_cooling]/Rad.nuc[fast_cooling])**(kappas.kappa12)) ** (kappas.kappa11) * (1+(nu[where_fast_cooling]/Rad.num[fast_cooling])**(kappas.kappa2p))**(kappas.kappa12inv)

        P_out[where_slow_cooling] =  Rad.PmaxS[slow_cooling] * ((nu[where_slow_cooling]/Rad.num[slow_cooling])**(kappas.kappa33) + (nu[where_slow_cooling]/Rad.num[slow_cooling])**(kappas.kappa3p))**(kappas.kappa13inv) * (1+(nu[where_slow_cooling]/Rad.nuc[slow_cooling])**(kappas.kappa42))**(kappas.kappa14)
                
    return P_out


class BL_constants:
    def __init__(self,p):
        self.phipF = 1.89 - 0.935*p + 0.17*p**2
        self.phipS = 0.54 + 0.08*p
        self.XpF = 0.455 + 0.08*p
        self.XpS = 0.06 + 0.28*p
        



def alphanu_func(nu , num , nuc , fast_cooling , p):
    ##########################################################################
    ### Returns broken power-law approximation of alpha_nu, not normalized ###
    ##########################################################################
    import numpy as np
    alpha_out = np.zeros(len(nu))  ###Allocating space
    if fast_cooling:
        if len(nu) == 1:
            if nu<=nuc:
                alpha_out = (nu/nuc)**(-5/3.)
            elif (nuc < nu) and (nu <= num):
                alpha_out = (nu/nuc)**(-3)
            elif num<nu:
                alpha_out = (num/nuc)**(-3) * (nu/num)**(-(p+5)/2)

        else:
            nu_less_nuc = np.where(nu<=nuc)
            nuc_less_nu_less_num = np.where((nuc<nu)&(nu<=num))
            num_less_nu = np.where(num<nu)
        
            alpha_out[nu_less_nuc] = (nu[nu_less_nuc]/nuc[nu_less_nuc])**(-5/3.)
            alpha_out[nuc_less_nu_less_num] = (nu[nuc_less_nu_less_num]/nuc[nuc_less_nu_less_num])**(-3)
            alpha_out[num_less_nu] = (num[num_less_nu]/nuc[num_less_nu])**(-3) * (nu[num_less_nu]/num[num_less_nu])**(-(p+5)/2)

    else:
        if len(nu) == 1:
            if nu<=num:
                alpha_out = (nu/num)**(-5/3.)
            elif (num<nu) and (nu<=nuc):
                alpha_out = (nu/num)**(-(p+4)/2)
            elif nuc<nu:
                alpha_out = (nuc/num)**(-(p+4)/2) * (nu/nuc)**(-(p+5)/2)

        else:
            nu_less_num = np.where(nu<=num)
            num_less_nu_less_nuc = np.where((num<nu)&(nu<=nuc))
            nuc_less_nu = np.where(nuc<nu)

            alpha_out[nu_less_num] = (nu[nu_less_num]/num[nu_less_num])**(-5/3.)
            alpha_out[num_less_nu_less_nuc] = (nu[num_less_nu_less_nuc]/num[num_less_nu_less_nuc])**(-(p+4)/2)
            alpha_out[nuc_less_nu] = (nuc[nuc_less_nu]/num[nuc_less_nu])**(-(p+4)/2) * (nu[nuc_less_nu]/nuc[nuc_less_nu])**(-(p+5)/2)
    
    return alpha_out
