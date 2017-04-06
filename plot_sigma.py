#This function reads in the posterior points, finds the 68% range and feeds it to the lightcurve generator, in order to plot the 1-sigma field on the lightcurve plot
def binner(x_array,y_array,bins):
    import numpy as np
    x_output , y_output = np.zeros(bins) , np.zeros(bins)
    bin_space = np.linspace(x_array[0],x_array[-1],bins+1)

    for i in range(bins):
        lower_lim = np.argmin(np.abs(x_array - bin_space[i]))
        upper_lim = np.argmin(np.abs(x_array - bin_space[i+1]))
        if lower_lim == upper_lim:
            x_output[i] = x_array[lower_lim]
            y_output[i] = 0
            continue
        else:
            x_output[i] = np.sum(x_array[lower_lim:upper_lim]) / (upper_lim - lower_lim)   #summing all elements  lower_lim <= i < upper_lim and normalizing
            y_output[i] = np.sum(y_array[lower_lim:upper_lim])    #Summing all elements, without normalizing
    
    return x_output , y_output , bin_space


def get_one_sigma(n_params,sigma_range,i_mode):

    from fitterFunction import modelFunc
    from options import userOptions
    import numpy as np
    from matplotlib import pyplot as plt
    from numpy.core.defchararray import split as splt

#Read in chains/1-post_separate.dat

    file_name = 'chains/1-post_separate.dat'

#Finding number of modes

    plot_area = True

    try:
        read_text = open(file_name,'r')
        dat_text = read_text.read().split('\n')
        read_text.close()
        raw_text = np.array(dat_text)
    except:
        print 'No such file %s. Now exiting'%file_name
        raise SystemExit(0)

    empties = np.where(raw_text=='')[0]
    
    break_mode = np.array([],dtype=int)
    print empties
    modes = 0

    ### Finding the breaks between modes
    for i in empties:
        try:
            if raw_text[i+1] == '': 
                modes += 1
                break_mode = np.append(break_mode,i+2)
        except:continue
    break_mode = np.append(break_mode,empties[-1]+1)
    print break_mode
    print '%s has %d nodes'%(file_name,modes)


    if i_mode == 0: ### i_mode == 0 means we want the top 68% from the entire selection
        in_data = np.zeros([break_mode[-1],n_params+2])
        mode_limits = [break_mode[0],break_mode[-1]]
    else:
        in_data = np.zeros([break_mode[i_mode]-break_mode[i_mode-1]-1,n_params+2])
        mode_limits = [break_mode[i_mode-1],break_mode[i_mode]]

    i_counter = 0
    for i_select in range(mode_limits[0],mode_limits[1]): ### Even if the last entry in mode_limits is not two lines from the last entry of the previous mode, the code allows us to run into a couple of empty lines
        if raw_text[i_select] == '': 
            continue
        else:
#            print len(map(float,raw_text[i_select].split()))
#            print len(in_data[i_counter])
#            print len(in_data)
            in_data[i_counter] = map(float,raw_text[i_select].split())
            
#            in_data[i_counter,0] = float(raw_text[i_select].split()[0])
#            in_data[i_counter,0] = float(raw_text[i_select].split()[1])
            i_counter += 1
#            print i_counter
#        raw_input(in_data)
        

    length = len(in_data)

    tot_prob = 0.
    
    which_points = np.array([],dtype=int)
    order = np.argsort(-in_data[:,0]) #Sorted on probability, highest first
    sorted_prob = in_data[order,0] 

    print '---\nMode no. %d\n---'%i_mode

    if i_mode == 0: print 'chains/1-post_separate.dat has %d points'%length
    else: print 'The %d%s mode of chains/1-post_separate.dat has %d points'%(i_mode,'st'*(i_mode==1)+'nd'*(i_mode==2)+'rd'*(i_mode==3)+'th'*(i_mode>3) , length)

    for i in range(length): ### Going through the points from the start, until the desired percentage is reached
        tot_prob += sorted_prob[i]
#        print tot_prob
        which_points = np.append(which_points,order[i])
        if tot_prob >= (sigma_range): break
    print 'Taking %d of these to make %.2f%%'%(len(which_points) , sigma_range*100)
    ### Finding difference in total chi2
    
        
    highest_chi2 , lowest_chi2 = np.max(in_data[which_points,1]) , np.min(in_data[which_points,1])
    print 'Highest chi2 in distribution: %s\nLowest chi2 in distribution: %s\nchi2 difference: %s\n'%(highest_chi2,lowest_chi2,highest_chi2-lowest_chi2)
    
    return in_data[which_points]




