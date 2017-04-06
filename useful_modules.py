class cgs_constants:
    def __init__(self):
        self.c =  2.9979e10
        self.mp = 1.6726e-24
        self.me = 9.1094e-28
        self.mpe = self.mp+self.me
        self.mue = self.me/self.mp
        self.hcgs = 6.6260755e-27   #Planck constant in cgs
        self.kB = 1.380658e-16
        self.sigmaT = 6.6524e-25
        self.qe = 4.803204e-10
        self.sigma_B = 5.6704e-5  ### Stephan Boltzann constant in cgs



def factor_of_ten(in_value):
    ### returns the factor of ten of in_value. For example, 542,3 returns 2
    import numpy as np
    if in_value == 0:
        return 0
    else:
        return int(np.floor(np.log10(in_value)))
    

def round_off(input_value,round_factor):
    import numpy as np
    ### Rounds input_value off to round_factor significant figures

    ### Treats all values as positive
    invalue_sign = np.sign(input_value)
    input_value = np.abs(input_value)

    try: ### Input is array
        out_array = np.copy(input_value)

        for i_round in range(len(input_value)): 
            if input_value[i_round] == 0:
                continue
        
            fot = int(np.log10(input_value[i_round]))
            fot -= (fot < 0)
            
            out_array[i_round] = round(input_value[i_round] / 10**fot,round_factor) * 10**fot
        return out_array
    except:
        fot = int(np.log10(input_value))
        fot -= (fot < 0)

        return round(input_value / 10**fot,round_factor) * 10**fot  *  invalue_sign

def binner(in_array,n_points,option):
    ### Bins in_array to n_points equally spaced bins. option='linear' returns linear sum, option='square' returns a geometric sum
    import numpy as np
    n_bins = len(in_array) / n_points  #Number of bins
    output = np.zeros(n_bins)
    for i in range(n_bins):
        if option == 'linear':
            output[i] = np.sum(in_array[i*n_points:(i+1)*n_points]) / n_points
        elif option == 'square':
            output[i] = np.sqrt(np.sum(in_array[i*n_points:(i+1)*n_points]**2)) / n_points
        else:
            raise NameError('Invalid bin option %s'%option)
    remaining_points = (len(in_array) % n_points)
    if remaining_points != 0: #Number of points is not devidable by number of bins. We have to construct the last bin separatly
        if option == 'linear':
            output = np.append(output,np.sum(in_array[n_bins*n_points:])/remaining_points)
        elif option == 'square':
            output = np.append(output,np.sqrt(np.sum(in_array[n_bins*n_points:]))/remaining_points)
                               
    return output


### Function that reduces number of ticks in plot. Matplotlib has a tendency to make plots look cramped
def reduce_ticks(in_plot,first_plot=False,last_plot=True,lower_xlim=float('-inf'),upper_xlim=float('inf'),include_y=False): 
    import numpy as np
    ### Boolean include_y is true if a parameter is on the y axis, false if normalized probability is

    if include_y:
        iterator = ['x','y']
    else:
        iterator = ['x']
        yticks = np.linspace(0,1.,3)

    for i_ticks in iterator:
        exec('in_ticks = in_plot.get_%sticks()'%i_ticks)


        in_ticks = in_ticks[np.where(in_ticks >= 0.)]

        ### If this is not the last plot in a panel (not last_plot), we don't want to few ticks because the last will be removed

        if len(in_ticks) % 3 == 0 or len(in_ticks) == 7: ### 3 ticks
            if not last_plot:
                if len(in_ticks) == 6:
                    out_ticks = [in_ticks[0] , (in_ticks[1]+in_ticks[2])/2 , in_ticks[3]]

                elif len(in_ticks) == 7:
                    out_ticks = in_ticks[:-1:2]
                    
                elif len(in_ticks) == 9:
                    out_ticks = in_ticks[::3]
                elif len(in_ticks) == 12:
                    out_ticks = [in_ticks[0],in_ticks[3],in_ticks[6],in_ticks[9]]
            else:
                if len(in_ticks) == 3: ### Odd number of ticks, picking centre value for centre tick, rounded to three digits
                    out_ticks = in_ticks
                elif len(in_ticks) == 9:
                    out_ticks = in_ticks[::4]
                else: ### Even number of ticks, only 2 ticks
                    out_ticks = [in_ticks[0],round_off(in_ticks[(len(in_ticks)-1)/2],2),in_ticks[-1]]

        elif len(in_ticks) % 4 == 0: ### 4 ticks. If there are four ticks going in, we don't alter it
            if not last_plot:
                if len(in_ticks) == 4:
                    out_ticks = in_ticks[:-1]
                elif len(in_ticks) == 8:
                    out_ticks = in_ticks[:-2:2]
                elif len(in_ticks) == 12:
                    out_ticks = in_ticks[::3]
            else:
                if len(in_ticks) == 4: 
                    out_ticks = np.copy(in_ticks)
                elif len(in_ticks) == 8: ### Three ticks, where the centre is the mean of the two nearest ticks, rounded to three digits
                    out_ticks = [in_ticks[0],round_off((in_ticks[3]+in_ticks[4])/2,3),in_ticks[-1]]
            
        elif len(in_ticks) % 5 == 0:
            if not last_plot:
                if len(in_ticks) == 5:
                    out_ticks = in_ticks[:-1]
                elif len(in_ticks) == 10:
                    out_ticks = in_ticks[:-1:3]
            else:
                if len(in_ticks) % 2 == 0: ### Even number of ticks: 4 ticks
                    out_ticks = [in_ticks[0],in_ticks[(len(in_ticks)-1)/3],in_ticks[(len(in_ticks)-1)*2/3],in_ticks[-1]]
                else: ### Odd number of ticks: 3 ticks
                    out_ticks = [in_ticks[0],round_off(in_ticks[(len(in_ticks)-1)/2],2),in_ticks[-1]]
            

        else: out_ticks = in_ticks
        
        exec('%sticks = np.copy(out_ticks)'%i_ticks)
    ### Rounding off; above routine may create floats infinitesimally close to incoming value


        

    in_plot.set_xticks(xticks)
    in_plot.set_yticks(yticks)        


