def stats_subplot(this_plot,x_grid_in,y_grid_in,plotLineIn,plot_multinest_gaussian,plot_dim,param_dim,first_plot,last_plot,numberOfPlots,param_counter,mstats_text,which_plot,gauss_linestyle=None,midLinePos=None,midLineColor=None,plot_lims=[float('-inf'),float('inf')]):
    import numpy as np
    from useful_modules import factor_of_ten
    import matplotlib.ticker as ticker
    
    line_width = 5.0
    gauss_line_width = 5.0

    ### Normalizing all peaks to 1.0
    y_grid_in = y_grid_in / np.max(y_grid_in)

    ### Removing low-probability tails
    for i_tail in range(len(y_grid_in)-1,0,-1):
        if y_grid_in[i_tail] > 0.01:
            break
    for j_tail in range(len(y_grid_in)):
        if y_grid_in[j_tail] > 0.01:
            break

    y_grid_in = y_grid_in[j_tail:i_tail+1]
    x_grid_in = x_grid_in[j_tail:i_tail+1]

    yupper = 1.
    x_grid_selection = np.where((x_grid_in>plot_lims[0]) & (x_grid_in<plot_lims[1]))

    if param_dim == 'log':

        this_plot.plot(10**x_grid_in[x_grid_selection],y_grid_in[x_grid_selection],plotLineIn,linewidth=line_width)
        if last_plot: ### Plotting input line. last_plot is only True if this is the last plot, and user wants to plot an input parameter marker

            this_plot.plot(10**np.array([midLinePos,midLinePos]),[0,yupper*1.1],color=midLineColor,linewidth=line_width)
            this_plot.set_ylim([0,yupper])
            
        min_xvalue_out = min(10**x_grid_in[x_grid_selection])
        max_xvalue_out = max(10**x_grid_in[x_grid_selection])

    elif param_dim == 'lin':
        if not first_plot:
            x_lims = np.zeros(2)
            x_lims[:] = this_plot.get_xlim()
        this_plot.plot(x_grid_in[x_grid_selection],y_grid_in[x_grid_selection],plotLineIn,linewidth=line_width)
        if not first_plot:
            if x_grid_in[x_grid_selection][0] < x_lims[0]:
                x_lims[0] = np.copy(x_grid_in[x_grid_selection][0])
            if x_grid_in[x_grid_selection][1] > x_lims[1]:
                x_lims[1] = np.copy(x_grid_in[x_grid_selection][1])
            this_plot.set_xlim(x_lims)
        if last_plot:

            this_plot.plot([midLinePos,midLinePos],[0,yupper*1.1],color=midLineColor,linewidth=line_width)
            this_plot.set_ylim([0,yupper])
            deltaMin = (x_grid_in[x_grid_selection][-1] - x_grid_in[x_grid_selection][0])
            facTen = factor_of_ten(x_grid_in[0])
#            this_plot.set_xlim([x_grid_in[x_grid_selection][0]-deltaMin,x_grid_in[x_grid_selection][-1]+deltaMin])
        min_xvalue_out = min(x_grid_in[x_grid_selection])
        max_xvalue_out = max(x_grid_in[x_grid_selection])

    elif param_dim == 'deg':
        if not first_plot:
            x_lims = np.zeros(2)
            x_lims[:] = this_plot.get_xlim()

        this_plot.plot(x_grid_in[x_grid_selection]*180/np.pi,y_grid_in[x_grid_selection],plotLineIn,linewidth=line_width)
        #this_plot.plot(x_grid_in[x_grid_selection],y_grid_in[x_grid_selection],plotLineIn,linewidth=line_width)
        if not first_plot:
            if min(x_grid_in[x_grid_selection])*180/np.pi < x_lims[0]:
                x_lims[0] = np.copy(x_grid_in[x_grid_selection][np.argmin(x_grid_in[x_grid_selection])]*180/np.pi)
            if max(x_grid_in[x_grid_selection])*180/np.pi > x_lims[1]:
                x_lims[1] = np.copy(x_grid_in[x_grid_selection][np.argmax(x_grid_in[x_grid_selection])]*180/np.pi)

            this_plot.set_xlim(x_lims)

        if last_plot:

            xlower,xupper = this_plot.get_xlim()
            this_plot.plot([180/np.pi*midLinePos,180/np.pi*midLinePos],[0,yupper*1.1],color=midLineColor,linewidth=line_width)
            this_plot.set_xlim([0,xupper])
            this_plot.set_ylim([0,yupper])
        deltaMin = (x_grid_in[x_grid_selection][-1] - x_grid_in[x_grid_selection][0]) * 180 / np.pi
        facTen = factor_of_ten(x_grid_in[x_grid_selection][0])
        min_xvalue_out = min(x_grid_in[x_grid_selection] * 180 / np.pi)
        max_xvalue_out = max(x_grid_in[x_grid_selection] * 180 / np.pi)

    if plot_multinest_gaussian: plot_mstats(x_grid_in[x_grid_selection],y_grid_in[x_grid_selection],n_params,this_plot,gauss_linestyle,gauss_line_width,param_dim,mstats_text,param_counter)    
    for xtick in this_plot.xaxis.get_major_ticks(): xtick.label.set_fontsize(40)
    for ytick in this_plot.yaxis.get_major_ticks(): ytick.label.set_fontsize(40)


    return min_xvalue_out,max_xvalue_out
