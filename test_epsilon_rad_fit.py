def test_epsrad(epsilon_rad,B,Gamma,delta_m,p,gamma_max,epsilone,t_comoving_step):
    import numpy as np
    me =  9.1094e-28 
    c = 2.9979e10
    mp = 1.6726e-24
    mue = me/mp
    sigmaT = 6.6524e-25
    #gammae = (1-epsilon_rad) * (Gamma - 1) * epsilone/mue + 1                   
    gamma_min = (p-2)/(p-1)*(Gamma - 1) * epsilone/mue * (1-epsilon_rad)

    ### Using the scale time of radiative losses
    #tstep = 6*np.pi*c*me / sigmaT / epsilone / B**2 * (3-p)/(p-2) * gamma_min**(2-p) / (gamma_max**(3-p) - gamma_min**(3-p))
    tstep = t_comoving_step
    gamma_c_mean = 6*np.pi*me*c / sigmaT / B**2 / tstep

    dE2sh_el = (Gamma-1) * delta_m * c**2 * epsilone
    
    N0 = delta_m / me * (p-2) * gamma_min**(p-2) * (Gamma - 1) * epsilone

    #dErad2 = N0*me*c**2*((gamma_c_mean**(2-p)-gamma_max**(2-p))/(3-p)/(p-2) - (gamma_min**(3-p)*(gamma_c_mean**-1 - gamma_max**-1))/(3-p) + (gamma_c_mean**(2-p) - gamma_max**(2-p))/(2-p)**2   -  gamma_max**(2-p)*np.log(gamma_c_mean/gamma_max) / (2-p) )

    dErad2 = N0*me*c**2*((gamma_c_mean**(2-p))/(3-p)/(p-2) - (gamma_min**(3-p)*(gamma_c_mean**-1))/(3-p) + (gamma_c_mean**(2-p))/(2-p)**2   )

    return dErad2 / dE2sh_el - epsilon_rad
