def plotMarginal():
    import numpy as np
    import pymultinest
    from matplotlib import pyplot as plt
    outputFiles_rawbase = raw_input('Type path to chains folder: ')
    prefix = '%s/chains/1-'%outputFiles_rawbase
    
    open_params = open('%s/chains/.tmp_myCode/parameters.txt'%outputFiles_rawbase)
    paramNames = open_params.read().split()
    open_params.close()
    
    n_params = len(paramNames)

    a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
    s = a.get_stats()

    p = pymultinest.PlotMarginalModes(a)

    for i in range(1,n_params):
#        plt.subplot(n_params, n_params, n_params * i + i + 1)
        #p.plot_marginal(i, ls='-', color='blue', linewidth=3)
        #p.plot_marginal(i, with_ellipses = True, dim2 = True, with_points = False, grid_points=50)
        
        
        for j in range(i):

            plt.subplot(n_params-1, n_params-1, (n_params-1) * j + i)
            p.plot_conditional(i, j, with_ellipses = False, with_points = False, grid_points=30)
            
        plt.subplot(n_params-1,n_params-1,(n_params-1)*(i-1) + i)
        plt.xlabel(paramNames[i])
        plt.ylabel(paramNames[i-1])

    plt.show()

def plot_contours(path_to_chains):
    import numpy as np
    import os
    from matplotlib import pyplot as plt
    
    ### Getting number of modes
    open_stats = open('%s/chains/1-stats.dat'%path_to_chains,'r')
    read_stats = open_stats.read().split('Total Modes Found:')
    open_stats.close()

    n_modes = int(read_stats[1].split()[0])

    ### Reading posterior files
    open_posterior = open('%s/chains/1-post_separate.dat'%path_to_chains,'r')
    read_posterior = open_posterior.read().split('\n')
    open_posterior.close()

    n_params = len(read_posterior[2].split())-2
    raw_input(n_params)

    while read_posterior[0] == '': read_posterior = read_posterior[1:]

    print 'Posterior has %d mode%s'%(n_modes,'s'*(n_modes>1))
    modes_split = np.where(read_posterior=='')

    ### probability has dimensions number of modes x number of data points
    probability = np.zeros([n_modes,len(read_posterior)])
    ### parameters has dimenstions number of modes x number of parameters x number of data points
    parameters = np.zeros([n_modes,n_params,len(read_posterior)])
    first_line = np.zeros(n_modes,dtype=int)
    
    n_points = 0

    ### Interpreting posterior files
    for i_modes in range(n_modes):
        
        for read_line in range(len(read_posterior)-first_line[i_modes]):
            this_line_in = read_posterior[read_line + first_line[i_modes]]
            if this_line_in == '': 
                if read_line+1 > n_points: n_points = np.copy(read_line) + 1
                print '%d points in mode %d'%(read_line+1,i_modes+1)
                break  #New mode
            this_line = map(float,this_line_in.split())
            probability[i_modes,read_line] = this_line[0]
            try: parameters[i_modes,:,read_line] = this_line[2:]
            except:
                raise NameError('Numbers of parameters should be %d, not %d'%(len(this_line[2:]),n_params))


        if i_modes != (n_modes - 1): 
            first_line[i_modes+1] = np.copy(read_line)
            while read_posterior[first_line[i_modes+1]] == '': first_line[i_modes+1] += 1

    ### Sanity check
    if np.abs((np.sum(probability) - n_modes))  > 1e-4:
        print '\n\nWarning! Sum of probability in posterior points is n_modes - p = %f\n\n'%(n_modes-np.sum(probability))


    return probability, parameters
