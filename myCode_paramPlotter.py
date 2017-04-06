import numpy as np
from matplotlib import pyplot as plt
import os
import glob
from EATS_func import BL_constants

tobs_scaled = False
this_path_input = raw_input('Is the Parameters directory stored in this directory? ([y]/n): ')
use_this_path = (this_path_input == 'y') or (this_path_input == 'Y') or (this_path_input == '')
if use_this_path:
    this_path = os.path.abspath('Parameters/')

else:
    this_path = raw_input('Type path to Parameters directory: ')

### Loading input parameters
"""
load_parameters_txt_input = raw_input('Load input parameters from parameters.txt? (otherwise from options.py) ([y]/n): ')
load_parameters_txt = not ((load_parameters_txt_input == 'n') or (load_parameters_txt_input == 'N'))
"""
load_parameters_txt = True
if load_parameters_txt:
    txt_open = open('%s/../parameters.txt'%this_path)
    txt_read = txt_open.read()
    txt_open.close()

    inlines_temp = np.array(txt_read.split('\n'))

    ### Removing last line
    if inlines_temp[-1] == '':
        inlines_temp = inlines_temp[:-1]

    ### Splitting lines
    inlines = np.zeros([len(inlines_temp),2] , dtype='S256')
    for i_inl in range(len(inlines_temp)):
        inlines[i_inl] = inlines_temp[i_inl].split('=')
        

### Calculating parameters
reverseShock = os.path.isfile('%s/BRS.txt'%this_path)

parameter_list_input = glob.glob('%s/*.txt'%this_path)
parameter_list = np.zeros(len(parameter_list_input),dtype='S256')

normalize_input = raw_input('Do you want to normalize the plot? (y/[n]): ')
normalize = (normalize_input == 'Y') or (normalize_input == 'y')

for i_plist in range(len(parameter_list)):
    parameter_list[i_plist] = parameter_list_input[i_plist].split('/')[-1].split('.')[0]
    exec('loaded_parameter = np.loadtxt(\'%s\')'%parameter_list_input[i_plist])

    if len(np.where(loaded_parameter>0)[0]) == 0: ### Negative parameter
        loaded_parameter *= -1

    if normalize:
        loaded_parameter /= np.max(loaded_parameter)

    exec('%s = np.copy(loaded_parameter)'%(parameter_list[i_plist]))
    
### Constants
p_FS_line = np.where(inlines == 'p_FS')[0]
p_FS = float(inlines[p_FS_line,1][0])
if reverseShock:
    p_RS_line = np.where(inlines == 'p_RS')[0]
    p_RS = float(inlines[p_RS_line,1][0])

mp = 1.6726e-24
c = 2.9979e10
me = 9.1094e-28
sigmaT = 6.6524e-25
qe = 4.803204e-10

phipF , phipS , XpF , XpS = BL_constants(p_FS)

if reverseShock:
    phipFRS , phipSRS , XpFRS , XpSRS = BL_constants(p_RS)
    variables = np.array([
        phipS*11.17*(p_FS-1)*qe**3*(4*Gamma*rho/mp)*B/(3*p_FS-1)/me/(c**2)  ,
        phipF*2.234*qe**3*(4*Gamma*rho/mp)*B/me/c**2 ,
        phipS*11.17*(p_RS-1)*qe**3*(4*Gamma*rho4/mp)*BRS/(3*p_RS-1)/me/(c**2) ,
        phipF*2.234*qe**3*(4*Gamma*rho4/mp)*BRS/me/c**2 ,
        (6*np.pi*qe/sigmaT/B)**.5 , 
        (6*np.pi*qe/sigmaT/BRS)**.5 
    ])

else:
    variables = np.array([])
variable_names = np.array(['PmaxS_FS','PmaxF_FS','PmaxS_RS','PmaxF_RS','gamma_max_FS','gamma_max_RS'])
print reverseShock
print variables
print len(variables)

for i_set_var in range(len(variable_names)):
    exec('%s = variables[i_set_var]'%variable_names[i_set_var])

parameter_list = np.append(parameter_list , variable_names)

while True:

    ### Choosing what parameter or variable to plot
    print 'Choose parameter or variable to plot\n\nParameters:\n'
    for i_plotchoise in range(len(parameter_list)):
        print '(%d): %s'%(i_plotchoise+1 , parameter_list[i_plotchoise])
    


    choise = raw_input('Your choise (x,y): ')
    try:
        x_choise = int(choise.split(',')[0])-1
        y_choise = int(choise.split(',')[1])-1
    except:
        raw_input('Bad input. Press enter and try again: ')


    if parameter_list[x_choise] == 'tobs' and not tobs_scaled:

        tobs /= 86400
        tobs_scaled = True
    print 'plt.plot(%s,%s)'%(parameter_list[x_choise] , parameter_list[y_choise])
    exec('plt.plot(%s,%s)'%(parameter_list[x_choise] , parameter_list[y_choise]))

    plot_additional = raw_input('Plot additional curves? (y/[n]): ')
    if (plot_additional == 'y') or (plot_additional == 'Y'):
        continue

    if parameter_list[x_choise] == 'tobs':
        plt.xlabel('tobs [days]')
    else:
        plt.xlabel(parameter_list[x_choise])
    plt.ylabel(parameter_list[y_choise])
    plt.loglog()
    plt.show()
    

    break
