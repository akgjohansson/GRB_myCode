def set_subplots(params_list,plot_scale,x_labels,whereParam,paramTotalLength):
    import numpy as np
    from matplotlib import gridspec
    from matplotlib import pyplot as plt
    from matplotlib import ticker
    import matplotlib
    from matplotlib.ticker import ScalarFormatter 
    #Fontsizes
    xaxis_fontsize = 60
    yaxis_fontsize = 50

    n_params = len(params_list)
    #Plot labels       
    epsilon_label = r'$\epsilon_{\rm rad}$'
    epsilone_label = r'$\epsilon_{\rm e}$'
    epsilone_FS_label = r'$\epsilon_{\rm e,FS}$'
    epsilone_RS_label = r'$\epsilon_{\rm e,RS}$'
    E0_label = r'$E_0[{\rm erg}]$'
    n_label = r'$n_{\rm CM}[{\rm cm}^{-3}]$'
    A0_label = r'$A_0 [{\rm cm}^{-3+s}]$'
    R_ISM_label = r'$R_{\rm ISM}$'
    Gamma0_label = r'$\Gamma_0$'
    epsilonB_label = r'$\epsilon_{\rm B}$'
    epsilonB_FS_label = r'$\epsilon_{\rm B,FS}$'
    epsilonB_RS_label = r'$\epsilon_{\rm B,RS}$'
    p_label = r'$p$'
    p_FS_label = r'$p_{\rm FS}$'
    p_RS_label = r'$p_{\rm RS}$'
    theta0_label = r'$\theta_0$'
    alpha_label = r'$\alpha$'
    tprompt_label = r'$\Delta t_{\rm prompt}[s]$'
    T_coc_label = r'$T_{\rm coc}$ [K]'
    R_coc_label = r'$R_{\rm coc}$ [cm]'
    theta0_coc_label = r'$\theta_{\rm 0,coc}$ [$^{\circ}$]'
    Gamma0_coc_label = r'$\Gamma_{\rm 0,coc}$'
    N_coc_label = r'$N_{\rm coc}$'

    
    if n_params == 2:
        plt.figure(figsize=(18,12))
        G = gridspec.GridSpec(1,2)
        plot_grid = ['[0,0]','[0,1]']
        first = np.array([0])
        alone = []
    elif n_params == 3:
        plt.figure(figsize=(6*n_params,3*n_params))
        G = gridspec.GridSpec(1,3)
        plot_grid = ['[0,0]','[0,1]','[0,2]']
        first = np.array([0])
        alone = []
    elif n_params == 4:
        plt.figure(figsize=(6*n_params,6*n_params))
        G = gridspec.GridSpec(2,2)
        plot_grid = ['[0,0]','[0,1]','[1,0]','[1,1]']
        first = np.array([0,2])
        alone = []
    elif n_params == 5:
        plt.figure(figsize=(4*n_params,3*n_params))
        G = gridspec.GridSpec(2,6)
        plot_grid = ['[0,0:2]','[0,2:4]','[0,4:]','[1,:3]','[1,3:]']
        first = np.array([0,3])
        alone = []
    elif n_params == 6:
        plt.figure(figsize=(6*n_params,6*n_params))
        G = gridspec.GridSpec(2,3)
        plot_grid = ['[0,0]','[0,1]','[0,2]','[1,0]','[1,1]','[1,3]']
        first = np.array([0,3])
        alone = []
    elif n_params == 7:
        plt.figure(figsize=(2*n_params,3*n_params))
        G = gridspec.GridSpec(3,6)
        plot_grid = ['[0,:2]','[0,2:4]','[0,4:]','[1,:3]','[1,3:]','[2,:3]','[2,3:]']
        first = np.array([0,3,5])
        alone = []
    elif n_params == 8: ### Two columns, four rows
        plt.figure(figsize=(2*n_params,3*n_params))
        G = gridspec.GridSpec(3,6)
        plot_grid = ['[0,:2]','[0,2:4]','[0,4:]','[1,:2]','[1,2:4]','[1,4:]','[2,:3]','[2,3:]']
        first = np.array([0,3,6])
        alone = []
    elif n_params == 9: ### Three columns, three rows
        plt.figure(figsize=(2*n_params,3*n_params))
        G = gridspec.GridSpec(3,3)
        plot_grid = ['[0,0]','[0,1]','[0,2]','[1,0]','[1,1]','[1,2]','[2,0]','[2,1]','[2,2]']
        first = np.array([0,3,6])
        alone = []
    elif n_params == 10: ### Three
        plt.figure(figsize=(2*n_params, 3*n_params))
        G = gridspec.GridSpec(3,12)
        plot_grid = ['[0,:4]','[0,4:8]','[0,8:]' , '[1,:4]','[1,4:8]','[1,8:]' , '[2,:3]','[2,3:6]','[2,6:9]','[2,9:]']
        first = np.array([0,3,6])
        alone = []
        
    elif n_params == 11: ### Three rows
        plt.figure(figsize=(2*n_params, 3*n_params))
        G = gridspec.GridSpec(3,12)
        plot_grid = ['[0,:4]','[0,4:8]','[0,8:]' , '[1,:3]','[1,3:6]','[1,6:9]','[1,9:]' , '[2,:3]','[2,3:6]','[2,6:9]','[2,9:]']
        first = np.array([0,3,7])
        alone = []
    elif n_params == 12: ## Three columns, four rows
        plt.figure(figsize=(2*n_params , 3*n_params))
        G = gridspec.GridSpec(4,3)
        plot_grid = ['[0,0]','[0,1]','[0,2]','[1,0]','[1,1]','[1,2]','[2,0]','[2,1]','[2,2]','[3,0]','[3,1]','[3,2]']
        first = np.array([0,3,6,9])
        alone = []
    elif n_params == 13: ### Four rows
        plt.figure(figsize=(2*n_params , 3*n_params))
        G = gridspec.GridSpec(4,12)
        plot_grid = ['[0,:3]','[0,3:6]','[0,6:9]','[0,9:]','[1,:4]','[1,4:8]','[1,8:]','[2,:4]','[2,4:8]','[2,8:]','[3,:4]','[3,4:8]','[3,8:]']
        first = np.array([0,4,7,10])
        alone = []
    elif n_params == 14: ### Four rows
        plt.figure(figsize=(2*n_params , 3*n_params))
        G = gridspec.GridSpec(4,12)
        plot_grid = ['[0,:3]','[0,3:6]','[0,6:9]','[0,9:]','[1,:3]','[1,3:6]','[1,6:9]','[1,9:]','[2,:4]','[2,4:8]','[2,8:]','[3,:4]','[3,4:8]','[3,8:]']
        first = np.array([0,4,8,11])
        alone = []


    G.update(wspace=0)
    out_array = np.array([None]*paramTotalLength)
    for i_plot in range(n_params):
        
#        if whereParam[i_plot] == -1: continue ### Parameter user has choosen to omit
        if plot_scale[i_plot] == 'log': 
            scale_text = ',xscale=u\'log\''
            logscale_this = True
        else: 
            scale_text = ''
            logscale_this = False
        
        if np.sum(first == i_plot):
            share_text = ''
            visible_text = ''
            share_here_text = ',sharey=%s_subplot'%params_list[i_plot]
        else:
            share_text = share_here_text
            visible_text = '%s_subplot.yaxis.set_visible(False)'%params_list[i_plot]

        
        
        exec('%s_subplot = plt.subplot(G%s%s,xlabel=r\'%s\'%s,ylabel=\'Norm. Prob.\')'%(params_list[i_plot] , plot_grid[i_plot] , scale_text , x_labels[i_plot],share_text))
        
        #matplotlib.axes.Axes.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #formatter = matplotlib.ticker.ScalarFormatter()
        #formatter.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        #formatter.set_scientific(True)
        #formatter.set_powerlimits((-1,1))

        #if not logscale_this:
        #    exec('%s_subplot.ticklabel_format(style=\'sci\', axis=\'x\', scilimits=(0,2))'%(params_list[i_plot]))
        
        exec('out_array[whereParam[i_plot]] = %s_subplot'%params_list[i_plot])

        ### Making y-axes of all but the leftmost plots invisible
        if visible_text != '': exec(visible_text)

        exec('%s_subplot.xaxis.label.set_fontsize(xaxis_fontsize)'%params_list[i_plot])
        exec('%s_subplot.yaxis.label.set_fontsize(yaxis_fontsize)'%params_list[i_plot])

        ### Setting limits for scientific factoring of large or small ticks. Default is >1e-3 or <1e4

#        exec('formatter = %s_subplot.xaxis.ScalarFormatter()'%params_list[i_plot])
#        formatter.set_scientific(True)
#        formatter.set_powerlimits((-1,2))

        #
#        exec('%s_subplot.ticker.ScalarFormatter.set_powerlimits((-1,2))'%params_list[i_plot])
        

    return out_array , first

