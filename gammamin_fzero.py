def minimize_gammamin(gamma_c , gamma_max , p , epsilone , Gamma_minus_one , mue):
    import numpy as np
    from scipy import optimize

    ### Temporary
    from matplotlib import pyplot as plt
    ### End temporary

    gmin_func = lambda gmin: gmin_fzero(gmin , gamma_c , gamma_max , p , epsilone , Gamma_minus_one , mue)

    if np.sign(gmin_func(1e-5)) != np.sign(gmin_func(gamma_max*1.00001)):  ### Employ bisect method
        try:
            res = optimize.bisect(gmin_func , 1e-5 , gamma_max*1.00001)
        except:
            raise NameError('bisect')
        
        if res <1:
            return 1.
        else:
        
            return res
    else:  ### gamma_min > gamma_max
        #raise NameError('gamma_min > gamma_max')
        return gamma_max
        """
        #print 'does not have different signs'
        try:
            res,_,ler,msg = optimize.fsolve(gmin_func , 1.   ,   full_output=True)
            if ler != 1:
                print 'No solution! Message:\n'
                print msg
                print 'gamma_min =',res[0]
                print 'gamma_max =',gamma_max
                print 'gamma_c   =',gamma_c
                print gamma_max
                print gamma_max.shape
                print gmin_func(1e-3)
                print gmin_func(gamma_max)
                print gmin_func(gamma_max*1.00001)
                print gmin_func(609723672.392)
                gplot = np.logspace(-3,13,10000)
                output = np.zeros(len(gplot))
                for i,gp in enumerate(gplot):
                    output[i] = gmin_func(gp)
                plt.plot(gplot,output)
                plt.show()
            if res[0] > gamma_max:
                return gamma_max
            else:
                print res
                print 'gamma_min =',res[0]
                print 'gamma_max =',gamma_max
                print 'gamma_c   =',gamma_c
                print gamma_max
                print gamma_max.shape
                print gmin_func(1e-3)
                print gmin_func(gamma_max)
                print gmin_func(gamma_max*1.00001)
                print gmin_func(609723672.392)
                gplot = np.logspace(-3,13,10000)
                output = np.zeros(len(gplot))
                for i,gp in enumerate(gplot):
                    output[i] = gmin_func(gp)
                plt.plot(gplot,output)
                plt.show()
            
            return res[0]
        except:
            print 'failed'
            return (p-2)/(p-1)*(epsilone/mue*Gamma_minus_one+1)
        """

def gmin_fzero(gamma_min , gamma_c , gamma_max , p , epsilone , Gamma_minus_one , mue):
    import numpy as np
    """
    if special:
        numerator = gamma_min**(1-p)*(np.log(gamma_min/gamma_c) + gamma_min**-1+gamma_c**-1) + (gamma_max**(1-p) - gamma_min**(1-p))/(1-p) + (gamma_max**(-p)-gamma_min**(-p))/p
        denominator = gamma_min**(1-p)*gamma_c**-1 - gamma_min**-p - (gamma_max**(-p) - gamma_min**(-p))/p
    else:
    """
    #gamma_max_m1 = gamma_max - 1
    #gamma_c_m1 = gamma_c - 1
    #gamma_min_m1 = gamma_min - 1
    if gamma_c < gamma_min:  ### Fast cooling
        
        #if special:
            #print 'fast cooling'
        ### See note from 2016-06-16
        numerator = gamma_min**(1-p)*(np.log(gamma_min/gamma_c) + gamma_min**-1+gamma_c**-1) + (gamma_max**(1-p) - gamma_min**(1-p))/(1-p) + (gamma_max**(-p)-gamma_min**(-p))/p
        denominator = gamma_min**(1-p)*gamma_c**-1 - gamma_min**-p - (gamma_max**(-p) - gamma_min**(-p))/p
        
        #numerator = (gamma_max**(1-p) - gamma_min**(1-p))/(1-p) + (gamma_max**(-p)-gamma_min**(-p))/p
        #denominator = (gamma_min**(-p) - gamma_max**(-p))/p
        """
        ### Using dN/dgamma = (gamma-1)**(-p-1) instead of gamma**(-p-1)
        numerator = (gamma_max_m1**(1-p) - gamma_min_m1**(1-p)) * p
        denominator = (gamma_min_m1**-p-gamma_max_m1**-p) * (1-p)
        """
    else:
        if gamma_c > gamma_max:

            gamma_c = gamma_max

        numerator = (gamma_c**(2-p) - gamma_min**(2-p)) / (2-p) - (gamma_c**(1-p)-gamma_min**(1-p))/(1-p) + (gamma_c*gamma_max**(1-p) - gamma_c**(2-p))/(1-p) + (gamma_c*gamma_max**-p - gamma_c**(1-p))/p

        denominator = (gamma_c**(1-p) - gamma_min**(1-p))/(1-p) - (gamma_c*gamma_max**-p - gamma_c**(1-p))/p
        """
        ### Using dN/dgamma = (gamma-1)**(-p) instead of gamma**(-p)
        numerator = (gamma_c_m1**(2-p)-gamma_min_m1**(2-p))/(2-p) + (gamma_max_m1**(1-p)-gamma_c_m1**(1-p))/(1-p)
        if np.isnan(numerator) or np.isinf:
            print 'gamma_c_m1 =',gamma_c_m1
            print 'gamma_min_m1 =',gamma_min_m1
            raw_input('gamma_max_m1 = %s'%gamma_max_m1)
        denominator = (gamma_c_m1**(1-p)-gamma_min_m1**(1-p))/(1-p) - (gamma_max_m1**(-p)-gamma_c_m1**(-p))/p
        """

    #if denominator == 0:
        #return -float('inf')
    lhs = epsilone * Gamma_minus_one/mue
    return lhs - numerator / denominator
