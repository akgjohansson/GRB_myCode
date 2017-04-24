 #This program runs the main fitter

import pymultinest
import numpy as np
import time 
import random
import sys
import os
import glob
import resource
from cosmocalc import cosmocalc
from fitterFunction import modelFunc
from options import model_variables
from options import userOptions
tidPre = time.gmtime()
from decimal import Decimal
from radiation_modules import plot_details
from radiation_modules import exception_class
global numRuns,surfRingsOut,tdata,FdataInput,errorbarInput,numberOfEmpties

np.set_printoptions(threshold=np.nan)


def myPrior(cube,ndims,nparam):

    cubeIndex = np.array(np.where(UseOp.parametrar))
    for i in range(ndims): cube[i] = UseOp.paramLimits[cubeIndex[0,i],0] +  cube[i] * (UseOp.paramLimits[cubeIndex[0,i],1]-UseOp.paramLimits[cubeIndex[0,i],0])
    

def logLikelihood(cube,ndims,nparam,cm_FdataInput=None,cm_tdata=None,cm_errorbarInput=None,cm_numberOfEmpties=None,cm_numberOfPoints=None):
    global lowestChi2, numberOfPoints , numRuns , surfRingsOut , tdata , FdataInput , errorbarInput , numberOfEmpties

    if (UseOp.runOption=='LC') and UseOp.createMock: FdataInput, tdata, errorbarInput, numberOfEmpties,numberOfEmpties = cm_FdataInput, cm_tdata, cm_errorbarInput, cm_numberOfEmpties, cm_numberOfEmpties

    if UseOp.runOption == 'fit': 
        for i in range(np.sum(UseOp.parametrar)): 
            ModVar.new_value(np.array(np.where(UseOp.parametrar))[0] , cube) #Assigning input value from prior
    else:
        startTime = time.time()
        print "Number of parameters: %d"%np.sum(UseOp.parametrar)
            
    
    if (UseOp.runOption == 'LC'): 
        
        lightcurve,startJetBreak,endJetBreak,tempGrid = modelFunc(R,ModVar,UseOp,PlotDetails,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,numberOfPoints,np.sum(UseOp.parametrar),Plot_Exceptions,plot_SED)
        endTime = time.time()
        print "Timeuse: %f s"%(endTime-startTime)
        return lightcurve , tempGrid


    ### Sampler
    elif UseOp.runOption == 'fit': 
        
        ### Setting number of rings in EATS integrator. The total number of rings in integrator is dependent on theta0
        surfRingsOut = 50

        chi2 , _,_ = modelFunc(R,ModVar,UseOp,PlotDetails,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,numberOfPoints,np.sum(UseOp.parametrar),Plot_Exceptions)

        if UseOp.printProcess & UseOp.allowPrint:

            try:
                if chi2 < lowestChi2:

                    print numRuns
                    
                    if True: #Set this to True to allow the programme to print the process in a file while running fitting routine
                        try: 
                            interLog = open('interactive_log.txt','r')
                            inText = interLog.read()
                            interLog.close()
                            lowestReadChi2 = (float((inText.split('\n')[1]).split(':')[1]))
                            if lowestReadChi2 < chi2 : #If the file interactive_log.txt has a lower chi2 than the routine, then don't rewrite file. This may occur when running on many cores
                                lowestChi2 = lowestReadChi2
                                numRuns += 1
                                print numRuns
                                return - chi2 / 2
                        except:
                            os.system('touch interactive_log.txt')
                            inText = ''

                    utText =  "New lowest chi2:\n"
                    paramOut = "%s=%s"%(paramNames[0],ModVar.echo_value(0))
                    for j in range(1,len(ModVar.const_names)):
                        paramOut = "%s\n%s=%s"%(paramOut,paramNames[j],ModVar.echo_value(j))
                        ### why are the parameters printed in the wrong order?
                        print j
                        print paramNames[j]
                        print ModVar.echo_value(j)
                        print '-'*20
                        print ModVar.const_names[j]
                        
                    paramOut = "%s\nsurfaceRings = %d"%(paramOut,surfRingsOut)
                    utParams = open('parameters.txt','w')
                    utParams.write(paramOut)
                    utParams.close()
                    
                    printTime = time.gmtime()

                    printParams = "At %d:%d:%d: \nNew lowest chi2: %s\n"%(printTime[3],printTime[4],printTime[5],chi2)
                    for i in range(ndims): 
                        #Adding exceptions for print-out that prints a variable on a different form than the fitter is using
                        printParams = '%s\n%s = %s %s'%(printParams,paramNamesShort[whereParam[i]], 10**ModVar.echo_value(whereParam[i])*(preferredScale[whereParam[i]]=='log') + ModVar.echo_value(whereParam[i])*((preferredScale[whereParam[i]]=='deg')*180/np.pi + (preferredScale[whereParam[i]]=='lin')) , 'degrees' * (preferredScale[whereParam[i]]=='deg'))
                    printParams = '%s\nnumRuns = %d'%(printParams,numRuns)
                    printParams = '%s\nsurfaceRings = %d'%(printParams,surfRingsOut)
                    printParams = '%s\nlog10(Likelihood) = %s\n\n'%(printParams,-chi2/2)

                    lowestChi2 = chi2

                    if True: #Set this to True to allow the programme to print the process in a file while running fitting routine
                        interLog = open('interactive_log.txt','w')
                        interLog.write("%s%s"%(printParams,inText))
                        interLog.close()
                        doOnce = False
            except: lowestChi2 = chi2
        numRuns += 1
        return - chi2 / 2



############################################
#              Loading data....            #
############################################



def loadData():
    global tobsmax, numberOfPoints
    from useful_modules import binner
    #Loading input data

    if UseOp.useData:     #If UseOp.useData == True, data is collected from the mock-data directory, otherwise from the real-data dir.
        os.chdir(inputData)
        if UseOp.allowPrint: print "Getting data files from %s"%inputData
    else: 
        os.chdir(inputData)
        if UseOp.allowPrint: print "Getting data files from %s"%inputData

    dataFilesTemp = glob.glob("*.%s"%UseOp.fileExtention)
    outTemp = np.zeros(len(dataFilesTemp))
    tdataTemp,FdataTemp,errorbarTemp,errorIs = np.array([]),np.array([]),np.array([]),np.array([],dtype=bool)
    maxIndex = 0 
    for giveName in range(len(dataFilesTemp)): outTemp[giveName] = float(dataFilesTemp[giveName].split(UseOp.breakName)[0])

    ### Chosing what lightcurves to load, if user specified the input option 'band=?' ###
    try:
        if bandLimit[1] == 0:   #No upper limit
            bandLimit[1] = 1e30 #A rediculously high limit to facilitate a 'no upper limit'
        inputNumber = len(dataFilesTemp)
        keep = map(bool,np.ones(inputNumber))  #boolean array 'keep' returns True if file should be loaded
        
        for rInd in range(inputNumber):
            if (outTemp[rInd] >= bandLimit[0]) and (outTemp[rInd] <= bandLimit[1]): keep[rInd] = True
            else: keep[rInd] = False

        out = np.zeros(np.sum(keep))
        dataFiles = ['']*np.sum(keep)
        for rKeep in range(inputNumber):
            if keep[rKeep]:
                out[np.sum(keep[:rKeep+1]) - 1] = outTemp[rKeep]
                dataFiles[np.sum(keep[:rKeep+1]) - 1] = dataFilesTemp[rKeep]

    except:
        out = outTemp
        dataFiles = dataFilesTemp

    
    if UseOp.allowPrint: print "Getting data from %d lightcurves"%len(dataFiles)

    dataIndeces = np.array([0]*(len(dataFiles)+1)) #dataIndeces: The index after the last entry of the input data (se concatenation below)


    numberOfPoints = 0
    for i in range(len(dataFiles)):
        
        #Determining unit in file name
        if UseOp.inputFreqUnit == 'GHz': out[i] *= 1e9
        elif UseOp.inputFreqUnit == 'MHz': out[i] *= 1e6
        elif UseOp.inputFreqUnit == 'cm': out[i] = c / out[i]
    
        try:
            input_data = np.loadtxt(dataFiles[i])
            tTemp = input_data[:,0]
            FRead = input_data[:,1]
            sigmaTemp = input_data[:,2]
            errorIsTemp = np.ones(len(tTemp),dtype=bool)
            if UseOp.inputTimeFormat == 'day': tTemp *= 86400
            elif UseOp.inputTimeFormat == 'hour': tTemp *= 3600
            elif UseOp.inputTimeFormat == 'min': tTemp *= 60
            elif UseOp.inputTimeFormat == 'auto':
                if np.sum(tTemp)/len(tTemp) < 100:  ### Time is in days
                    tTemp *= 86400
                ### If input time is in seconds, nothing is done to it


        except:
            a = open(dataFiles[i])
            inText = a.read().split('\n')
            while inText[-1] == '':
                inText = inText[:-1]
            a.close()
        
            if (len(inText[0].split())!=1): #if the first entry in the .dat file is not the number of data points
                inputLength = 0
                for k in range(len(inText)):
                    if (inText[inputLength] != ''):
                        try:
                            map(float,inText[inputLength].split())   #Is the row a row of data points?
                            inputLength += 1
                        except: break
                isFirstTheNumber=0
            else: #if the first entry in the .dat file IS the number of data points
                inputLength = int(inText[0])
                if inputLength != len(inText[1:]): raw_input("Input data is of other length than is stated in top of file!")  #Sanity check
                isFirstTheNumber = 1
        

        
            tTemp,FRead,sigmaTemp,errorIsTemp = np.zeros(inputLength),np.zeros(inputLength),np.zeros(inputLength),np.zeros(inputLength,dtype=bool) #errorIsTemp is an array to say if the errorbar contains true errorbars (True) or upper limits (False)
        #Converting string line to float array
            for j in range(isFirstTheNumber,len(inText)):
                inRow = map(float,inText[j].split())
                tTemp[j-1],FRead[j-1] = np.array(inRow[:2])
            #Converting input time
                if UseOp.inputTimeFormat == 'day': tTemp[j-1] *= 86400
                elif UseOp.inputTimeFormat == 'hour': tTemp[j-1] *= 3600
                elif UseOp.inputTimeFormat == 'min': tTemp[j-1] *= 60
                elif UseOp.inputTimeFormat == 'auto':
                    if tTemp[j-1] < 100: ### Input time format is in days
                        tTemp[j-1] *= 86400
                try: 
                    sigmaTemp[j-isFirstTheNumber] = np.array(inRow[2]) 
                    errorIsTemp[j-isFirstTheNumber] = True
                except: 
                    sigmaTemp[j-isFirstTheNumber] = np.array(inRow[1])   #If the error bar is missing, sigma = flux
                    errorIsTemp[j-isFirstTheNumber] = False
        
            if inputLength != len(tTemp): raw_input("Input data is of other length than is stated in top of file!")  #Sanity check
        
        ### Sorting and binning input data
        indexOrder = np.argsort(tTemp)
        if bin_data_bool:
            tTemp = binner(tTemp[indexOrder],bin_data,'linear')
            FRead = binner(FRead[indexOrder],bin_data,'linear')
            sigmaTemp = binner(sigmaTemp[indexOrder],bin_data,'square')
            errorIsTemp = binner(errorIsTemp,bin_data,'linear')
            indexOrder = np.arange(len(tTemp))
        
        for indSort in range(len(tTemp)):
            tdataTemp = np.append(tdataTemp,tTemp[indexOrder[indSort]])
            FdataTemp = np.append(FdataTemp,FRead[indexOrder[indSort]])
            errorbarTemp = np.append(errorbarTemp,sigmaTemp[indexOrder[indSort]])
            errorIs = np.append(errorIs,errorIsTemp[indexOrder[indSort]])

        inputLength = len(tTemp)
        numberOfPoints += inputLength
        dataIndeces[i+1] = int(dataIndeces[i] + inputLength)
        if (dataIndeces[i+1]-dataIndeces[i]) > maxIndex: maxIndex = dataIndeces[i+1]-dataIndeces[i]
        


    firstColumn,FdataOut,errorbarOut,errorIsOut = np.zeros([len(dataFiles),maxIndex]),np.zeros([len(dataFiles),maxIndex]),np.zeros([len(dataFiles),maxIndex]),np.zeros([len(dataFiles),maxIndex],dtype=bool)

    os.chdir(currentPath)

    for i in range(len(dataFiles)):
        firstColumn[i] = np.concatenate([tdataTemp[dataIndeces[i]:dataIndeces[i+1]],[-1]*(maxIndex-(dataIndeces[i+1]-dataIndeces[i]))])
        FdataOut[i] = np.concatenate([FdataTemp[dataIndeces[i]:dataIndeces[i+1]],[-1]*(maxIndex-(dataIndeces[i+1]-dataIndeces[i]))])
        errorbarOut[i] = np.concatenate([errorbarTemp[dataIndeces[i]:dataIndeces[i+1]],[-1]*(maxIndex-(dataIndeces[i+1]-dataIndeces[i]))])
        errorIsOut[i] = np.concatenate([errorIs[dataIndeces[i]:dataIndeces[i+1]],[-1]*(maxIndex-(dataIndeces[i+1]-dataIndeces[i]))])

        minFdataThis , maxFdataThis = FdataTemp[dataIndeces[i] + np.argmin(FdataTemp[dataIndeces[i]:dataIndeces[i+1]])] , FdataTemp[dataIndeces[i] + np.argmax(FdataTemp[dataIndeces[i]:dataIndeces[i+1]])]
        try: 
            if minFdataThis < minFdata: minFdata = minFdataThis
        except: minFdata = minFdataThis
        try: 
            if maxFdataThis > maxFdata: maxFdata = maxFdataThis
        except: maxFdata = maxFdataThis
        minTdataThis , maxTdataThis = tdataTemp[dataIndeces[i] + np.argmin(tdataTemp[dataIndeces[i]:dataIndeces[i+1]])] , tdataTemp[dataIndeces[i] + np.argmax(tdataTemp[dataIndeces[i]:dataIndeces[i+1]])]
        try: 
            if minTdataThis < minTdata: minTdata = minTdataThis
        except: minTdata = minTdataThis
        try: 
            if maxTdataThis > maxTdata: maxTdata = maxTdataThis
        except: maxTdata = maxTdataThis


    

    return firstColumn,FdataOut,errorbarOut,errorIsOut,out


#######################################
#   Analysing and plotting results    #
#######################################

def plot_mstats(x_grid,y_grid_in,n_params,this_plot,gauss_linetype,gauss_line_width,plot_dim,mstats_text,param_number):
    from matplotlib import mlab

    if plot_dim == 'log' or plot_dim == 'lin': x_grid_in = np.copy(x_grid)
    elif plot_dim == 'deg': x_grid_in = x_grid * 180 / np.pi

    ### mu and sigma arrays. One per parameter
    mstats_mu = np.zeros(n_params)
    mstats_sigma = np.zeros(n_params)


    this_line = map(float,mstats_text[param_number+1].split())


    if param_number + 1 != this_line[0]:# Sanity check. If this fails the layout of the 1-stats.dat file is probably changed                                        
        print 'Something occur while reading chains/1-stats.dat. Now exiting.'
        raise SystemExit(0)

    mstats_mu = this_line[1]
    mstats_sigma = this_line[2]

    if len(x_grid_in) == 0: #No points in this distribution.
        print 'Parameter %s without distribution! Skipping this one.'%paramNames[whereParam[param_number]]
    else:
        grid_area = np.sum(np.concatenate([[(x_grid_in[1]-x_grid_in[0])/2],(x_grid_in[2:]-x_grid_in[:-2])/2,[(x_grid_in[-1]-x_grid_in[-2])/2]]) * y_grid_in)
        ### Plotting gaussians
        xlims = this_plot.get_xlim()
        if plot_dim == 'log': 
            x_plotarray = np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),100)
            plot_pdf = mlab.normpdf(np.linspace(np.log10(min(x_plotarray)),np.log10(max(x_plotarray)),100),mstats_mu,mstats_sigma)
        elif plot_dim == 'lin':
            x_plotarray = np.linspace(xlims[0],xlims[1],100)
            plot_pdf = mlab.normpdf(x_plotarray,mstats_mu,mstats_sigma)
        elif plot_dim == 'deg': 
            x_plotarray = np.linspace(xlims[0],xlims[1],100)
            plot_pdf = mlab.normpdf(x_plotarray,mstats_mu*180/np.pi,mstats_sigma*180/np.pi)
        this_plot.plot(x_plotarray , plot_pdf * grid_area,gauss_linetype,linewidth=gauss_line_width)


##################################
### Plotting comparing results ###
##################################

def plot_comres(this_plot,fit_pars,color,scale):
    n_comres = np.shape(fit_pars)[0]
    x_lims = np.zeros(2)
    x_lims[:] = this_plot.get_xlim()
    
    ### fit_pars has first dimension the number of the model , and second dimension has length 3: mid-value, lower-limit and upper-limit

    for i_comres in range(n_comres):


        if fit_pars[i_comres,0] != 0.:
            if scale == 'log':
                lower_lim = 10**(np.log10(fit_pars[i_comres,1]) - np.log10(fit_pars[i_comres,2]/fit_pars[i_comres,1])/5)
                upper_lim = 10**(np.log10(fit_pars[i_comres,1]) + np.log10(fit_pars[i_comres,2]/fit_pars[i_comres,1])/5)
            else: ### lin or deg
                lower_lim = fit_pars[i_comres,1] - (fit_pars[i_comres,2] - fit_pars[i_comres,1]) / 5
                upper_lim = fit_pars[i_comres,1] + (fit_pars[i_comres,2] - fit_pars[i_comres,1]) / 5
            if lower_lim < x_lims[0]: 
                x_lims[0] = np.copy(lower_lim)
            if upper_lim > x_lims[1]: 
                x_lims[1] = np.copy(upper_lim)



            this_plot.plot([fit_pars[i_comres,0],fit_pars[i_comres,0]] , [0.,1.1],color=color[i_comres],linewidth=5.0) #Central value
            this_plot.plot([fit_pars[i_comres,1],fit_pars[i_comres,1]] , [0., 1.1],linestyle='--',color=color[i_comres],linewidth=5.0) #Lower 1-sigma
            this_plot.plot([fit_pars[i_comres,2],fit_pars[i_comres,2]] , [0., 1.1],linestyle='--',color=color[i_comres],linewidth=5.0) #Upper 1-sigma

            this_plot.set_ylim([0,1])
    return x_lims





def plotter():
    from set_subplots import set_subplots
    from matplotlib import pyplot as plt
    from stats_subplot import stats_subplot
    from useful_modules import reduce_ticks
    import matplotlib as mpl
    import kde_weights as kde
    #Analyse output

    
    #Making sure directory ../Figures exist
    pathName = '%s/Figures/%s'%(os.path.abspath('..'),UseOp.GRBlabel)
    if os.path.isdir('%s/Figures/'%(os.path.abspath('..'))) == False: os.system('mkdir ../Figures/')
    if os.path.isdir('%s/Figures/%s/'%(os.path.abspath('..'),UseOp.GRBlabel)) == False: os.system('mkdir ../Figures/%s/'%UseOp.GRBlabel)    
    os.system('mkdir %s'%(pathName))
    mpl.rc('xtick', labelsize=40) 
    mpl.rc('ytick', labelsize=40) 

    #Number of plots
    numberOfPlots = input('Number of plots: ')

    #Create average plots
    sum_plots_input = raw_input('Do you want to create average plots? (y/[n]): ')
    if sum_plots_input == 'y' or sum_plots_input == 'Y': sum_plots = True
    else: sum_plots = False

    compare_results_input = raw_input('Do you want to compare with and plot other fit parameters? (y/[n]): ')
    compare_results = (compare_results_input=='Y') or (compare_results_input=='y')

    #Print input parameter line?
    if not compare_results:
        printInputLineInput = raw_input('Plot an input value line (use for fitting mock observations) (y/[n]):')
        printInputLine = (printInputLineInput == 'y') or (printInputLineInput == 'Y')
    else:
        n_comres = input('Number of models to input and compare with: ')
        printInputLine = False
    if printInputLine: midLineColor = raw_input('Color of the input value line: ')
    else: midLineColor = None
    pathToPlot = ['']*numberOfPlots        
    plotLine = np.zeros(numberOfPlots,dtype='S10')
    currentPlotPath = os.path.abspath('.')
    #n_average_plots = np.ones(numberOfPlots,dtype=int)



    ### Looping over plot paths
    ### In this loop we also choose if a gaussian fit should be plotted
    ### and collects the list of what parameters were used to create the posterior,
    ### stored in chains/.tmp_myCode/parameters.txt
    plot_multinest_gaussian = np.zeros(numberOfPlots,dtype=bool)
    n_params = np.zeros(numberOfPlots,dtype=int)
    for plotI in range(numberOfPlots):
        if sum_plots:
            #n_average_plots[plotI] = input('Number of plots to average over in plot %d: '%(plotI+1))
            pathToPlot[plotI] = raw_input('Type general path (using asterix) to the %d%s location of all chains folders: '%(plotI+1,'st'*(plotI==0)+'nd'*(plotI==1)+'rd'*(plotI==2)+'th'*(plotI>2)))
            plot_multinest_gaussian[i] = False
            gauss_color = ''
            gauss_linetype = ''
        else:
            new_path_input = 2
            while new_path_input:
                pathToPlotInput = raw_input('Path where to find chains/ no. %d: '%(plotI+1))
                if pathToPlotInput[0] == '~':
                    pathToPlot[plotI] = '%s%s'%(homePath , pathToPlotInput[1:])
                else:
                    pathToPlot[plotI] = os.path.abspath(pathToPlotInput)
                if not os.path.isdir(pathToPlot[plotI]):
                    print 'No such path %.'
                    if new_path_input == 1:
                        print 'Now exiting.'
                        raise SystemExit(0)
                    else:
                        print 'Please try once more'
                    new_path_input -= 1
                else:
                    new_path_input = False

            
            #plot_multinest_gaussian_input = raw_input('Plot MultiNest output parameters as gaussian? ([y]/n): ')
            plot_multinest_gaussian_input = 'n'
            plot_multinest_gaussian[plotI] = plot_multinest_gaussian_input == 'y' or plot_multinest_gaussian_input == 'Y' or plot_multinest_gaussian_input == ''
    
            try:gauss_linetype
            except: 
                gauss_linetype = np.zeros(numberOfPlots,dtype='S20')
                gauss_color = np.zeros(numberOfPlots,dtype='S20')
            if plot_multinest_gaussian[plotI]: 
                gauss_linetype[plotI] = np.copy(raw_input('Gaussian linetype: ([-],--,-.,:): '))
                if gauss_linetype[plotI] == '': gauss_linetype = '-'
                gauss_color[plotI] = np.copy(raw_input('Gaussian color: (default=g) '))
                if gauss_color[plotI] == '': gauss_color[plotI] = 'g'
            
            ### Loading chains/.tmp_myCode/parameters.txt
            param_list_path = '%s/chains/.tmp_myCode/parameters.txt'%pathToPlot[plotI]

            if not os.path.isfile(param_list_path):
                print 'No such file %s. Posterior files are probably produced with a deprecated code. Please run the new version of the code on the same chains files to make them up to date. Now exiting'%param_list_path
                raise SystemExit(0)

            ### list_of_plot_params contains a list of all parameters to be plotted, counting all choosen paths to chains files
            open_param_list = open(param_list_path,'r')
            exec('list_of_plot_params_%d = open_param_list.read().split()'%plotI)
            open_param_list.close()
            exec('n_params[plotI] = len(list_of_plot_params_%d)'%plotI)
            if plotI == 0: ### First iteration, no previously stored list of parameters to append to
                list_of_plot_params = np.copy(list_of_plot_params_0)

            else:
                ### Adding new params to list
                exec('list_of_params_tmp = np.copy(list_of_plot_params_%d)'%plotI)


                for a_param_list in range(len(list_of_params_tmp)):
                    if np.sum(list_of_params_tmp[a_param_list] == list_of_plot_params) == 0: ### If this is true, the parameter does not already exist in the array list_of_plot_params
                        list_of_plot_params = np.append(list_of_plot_params,list_of_params_tmp[a_param_list])

    
        ### Line color of marginal distribution plot
        plotLine[plotI] = raw_input('Line color. Default=\'b\': ')
        if plotLine[plotI] == '': plotLine[plotI] = r'b'



    ### Getting total number of paramters for all collected posteriors

    number_of_plot_params = len(list_of_plot_params)



    ### Giving the user the posibility to omit parameters when plotting
    plot_all_parameters = raw_input('Plot all %d parameters? ([y]/n): '%number_of_plot_params)
    if plot_all_parameters == 'n' or plot_all_parameters == 'N':
        ### Printing available parameters
        for print_params in range(number_of_plot_params):
            print '(%d): %s'%(print_params+1 , list_of_plot_params[print_params])
        ### Letting the user choose what parameters to omit
        while True:
            try: 
                omit_params = np.array(map(int,raw_input('Which parameter/s do you want to omit? Write corresponding number above, separate by space for many: ').split()))
                break
            except: 
                do_cancel = raw_input('Bad input. Enter \'c\' to cancel, press enter to try again: ')
                if do_cancel == 'c': raise SystemExit(0)
        ### Omitting chosen parameters
        omit_these_parameters = np.zeros(len(omit_params),dtype='S256')
        
        print 'Omitting %s'%(' and '.join(np.array(list_of_plot_params[omit_params-1])))
        for i_omit in range(len(omit_params)):
            omit_these_parameters[i_omit] = list_of_plot_params[omit_params[i_omit] - 1] ### The -1 is introduced because the user chooses parameter from no. 1 and up
#            list_of_plot_params = np.append(list_of_plot_params[:omit_params[i_omit] - 1 - i_omit],list_of_plot_params[omit_params[i_omit] - i_omit:]) ### The '- i_omit' is introduced because the length of list_of_plot_params is decreased for each iteration
            list_of_plot_params[omit_params[i_omit]-1] = None

    
    ### Locating index in parametrar array corresponding to list_of_plot_params
        ### whereParam_short contains reduced whereParam, ordered as programme default, not after largest y-value
    whereParam_short = np.zeros(len(list_of_plot_params),dtype=int)
    for i_whereParam in range(number_of_plot_params):
        if list_of_plot_params[i_whereParam] == str(None):
            whereParam_short[i_whereParam] = -1
            continue
        for j_whereParam in range(len(UseOp.parametrar)):
            if list_of_plot_params[i_whereParam] == paramNamesShort[j_whereParam]: 
                whereParam_short[i_whereParam] = np.copy(j_whereParam)
                break
            
            ### Sanity check - if loop goes through entire loop without finding a match, something is wrong. Possibly has the chains files been run by a deprecated code
            if j_whereParam == (len(UseOp.parametrar)-1):
                print 'Could not match all parameters in chains/.tmp_myCode/parameters.txt. Possibly has the chains files been produced by a deprecated code. Run the new code on the existing chains files and try again. Now exiting.'
                raise SystemExit(0)


    if compare_results:   #Input fit parameters for comparison
        ### Allocating

        compare_input = np.zeros([n_comres,3,number_of_plot_params])
        comres_infile = np.zeros(n_comres,dtype='S256')
        comres_color = np.zeros(n_comres,dtype='S20')


        for i_comps in range(n_comres):
            comres_infile[i_comps] = raw_input('Name of file with fit parameters%s. Leave empty for manual input: '%((' for model %d'%(plotI+1))*(numberOfPlots>1)))
            comres_color[i_comps] = raw_input('Color: ')


            comres_out_text = ''

            if not comres_infile[i_comps] == '':
                if True:
                    open_comres = open(comres_infile[i_comps])
                    read_comres = open_comres.read().split('\n')
                    open_comres.close()


                    for i_read_comres in range(len(read_comres)):
                        if read_comres[i_read_comres] == '':
                            continue
                        this_comres = read_comres[i_read_comres].split()[0].split(':')[0]
                        try:
                            this_param = np.where(this_comres == list_of_plot_params)[0][0]
                        except:
                            print 'Parameter %s is not represented in %s.'%(this_comres , comres_infile[i_comps])
                        compare_input[i_comps,:,this_param] = np.array(read_comres[i_read_comres].split()[1:],dtype=float)
                else:
                    print 'Something wrong. Now exiting'
                    raise SystemExit(0)
            else:
                for i_comres in range(len(list_of_plot_params)):

                    ### index i points to the entry in the parametrar array corresponding to the parameter we are into in this stage of the for-loop
                    try:
                        this_param = np.where(list_of_plot_params[i_comres] == list_of_plot_params)[0][0]
                    except: ### this param is set to None and is omitted
                        continue

                    print '\nEnter fit parameters of model %d. Leave empty to skip the parameter. %s'%(i_comps+1,'Enter 0 as 1-sigma to not plot the uncertainty range')
                    try:
                        compare_input[i_comps,0,this_param] = input('Central value for %s: '%(list_of_plot_params[i_comres]))
                    except: 
                        compare_input[i_comps,0,this_param] = 0.
                        compare_input[i_comps,1,this_param] = 0.
                        compare_input[i_comps,2,this_param] = 0.
                    if compare_input[i_comps,0,this_param] != 0.:
                        compare_input[i_comps,1,this_param] = input('Lower 1-sigma for %s: '%(list_of_plot_params[i_comres]))
                        if compare_input[i_comps,1,this_param] != 0:  #Will not plot uncertainty range if <-- == 0
                            compare_input[i_comps,2,this_param] = input('Upper 1-sigma for %s: '%(list_of_plot_params[i_comres]))
                        else:
                            compare_input[i_comps,1,this_param] = compare_input[i_comps,0,this_param]
                            compare_input[i_comps,2,this_param] = compare_input[i_comps,0,this_param]
                        if compare_input[i_comps,1,this_param] < 0: #Transforms relative standard dev. to absolute
                            compare_input[i_comps,1,this_param] += compare_input[i_comps,0,this_param]
                            compare_input[i_comps,2,this_param] += compare_input[i_comps,0,this_param]
                    
                    comres_out_text += '%s: %f  %f  %f\n'%(list_of_plot_params[i_comres] , compare_input[i_comps,0,this_param] , compare_input[i_comps,1,this_param] , compare_input[i_comps,2,this_param])
                save_comres_filename = raw_input('File name to save fit parameters of model %d: '%(i_comps+1))
                open_save_comres = open(save_comres_filename,'w')
                open_save_comres.write(comres_out_text)
                open_save_comres.close()

            
    
            
            

    create_temp_catalogue = False


    ### Lower and upper limits for subplots
    min_xvalue = np.zeros(len(UseOp.parametrar))
    max_xvalue = np.zeros(len(UseOp.parametrar))

    y_max = np.zeros(number_of_plot_params)
    plot_list_range_raw = range(number_of_plot_params) ### The indeces to all parameters, when used in for-looping

    iteration_length = np.zeros(numberOfPlots,dtype=int)

    x_grid_mean = np.zeros(number_of_plot_params) ### Stores mean values of x_grid, used to correct order-of-ten of x_grid

    ### Ask user if the isotropic or collimated energy should be plotted
    collimate_energy_input = raw_input('Plot the collimated energy? ([y]/n): ')
    collimate_energy = (collimate_energy_input != 'n') and (collimate_energy_input != 'N')

    ### Looping over plot paths
    for j in range(numberOfPlots):
        if sum_plots: sum_average_paths = glob.glob(pathToPlot[j])
        ### Creating average plots
        if sum_plots and (len(sum_average_paths) > 1): 
            
            x_grid_lims = np.zeros([n_params[j],2])
            for sum_average in range(len(sum_average_paths)):
                os.chdir(currentPlotPath)
                os.chdir(sum_average_paths[sum_average])
                print 'Moving to %s'%sum_average_paths[sum_average]
        
                a = pymultinest.Analyzer(n_params = n_params)
                s = a.get_stats()
                p = pymultinest.PlotMarginalModes(a)
                i_counter = 0
                if os.path.isdir('%s/temp_catalogue'%currentPlotPath) and j == 0 and sum_average == 0: 
                    os.system('rm %s/temp_catalogue/*'%currentPlotPath)
                elif not os.path.isdir('%s/temp_catalogue'%currentPlotPath): 
                    create_temp_catalogue = True
                    os.system('mkdir %s/temp_catalogue'%currentPlotPath)
                for i in np.where(UseOp.parametrar)[0]:

                    x_grid,y_grid = p.plot_marginal(i_counter, with_ellipses = True, with_points = False, grid_points=50)
                    np.savetxt('%s/temp_catalogue/grid_%d%d.txt'%(currentPlotPath,i,sum_average),[x_grid,y_grid])

                    #Finding extremes for the x-axis
                    if sum_average == 0:
                        x_grid_lims[i_counter] = np.array([min(x_grid),max(x_grid)])
                    else:
                        if x_grid_lims[i_counter,0] > min(x_grid): x_grid_lims[i_counter,0] = np.copy(min(x_grid))
                        if x_grid_lims[i_counter,1] < max(x_grid): x_grid_lims[i_counter,1] = np.copy(max(x_grid))
                    i_counter += 1
            os.chdir(currentPlotPath)
            print 'Moving to %s'%currentPlotPath
            mstats_text = None

        else:
            
                                                            

            os.chdir(pathToPlot[j])
        
            a = pymultinest.Analyzer(n_params = n_params[j])
            s = a.get_stats()

            #p = pymultinest.PlotMarginalModes(a)
    
            #Reading mutlinest output from chains/1-stats.dat
            if plot_multinest_gaussian[j]:
                open_mstats = open('chains/1-stats.dat')
                read_mstats = open_mstats.read().split('Dim No.')
                open_mstats.close()
                mstats_text = read_mstats[1].split('\n')
            else:
                mstats_text = 'None'
        
        limCounter = 0
    
### At this stage we have read chains/.tmp_myCode/parameters.txt list, and found a list of all the total paramters.
        
        if True:
            ### Getting list of parameters for the current chains we are analysing
            ### Constructing plot_list_range for the specific set of posterior under the j:th index
            exec('i_plot_list = np.copy(list_of_plot_params_%d)'%j) 
            exec('plot_list_range_%d = np.array([],dtype=int)'%j)

            ### Picking out indeces from plot_list_range_raw that is found in this set of posterior files
            
            for i_include_these in plot_list_range_raw:
                if np.sum(list_of_plot_params[i_include_these] == i_plot_list) == 1:
                    iteration_length[j] += 1

            ### Getting data
            in_data = a.get_data()
            post_in = in_data[:,0]
            x_grid_total = in_data[:,2:]

            param_mean = np.zeros(n_params[j])
            onesigma = np.zeros([n_params[j],2])
            twosigma = np.zeros([n_params[j],2])
            threesigma = np.zeros([n_params[j],2])
            post_prob = np.sum(post_in)

            ### Correcting for collimated energy
            if collimate_energy:
                try:
                    E0_param = np.where(i_plot_list == 'E0')[0]
                    theta0_param = np.where(i_plot_list == 'theta0')[0]
                    lowest_theta = np.argmin(x_grid_total[:,theta0_param])
                    x_grid_total[:,E0_param] = np.log10(10**x_grid_total[:,E0_param] * (1-np.cos(x_grid_total[:,theta0_param])))
                    latexParamNamesLin[np.where(paramNamesShort=='E0')[0]] = r'$E_0$'

                except:
                    raise NameError('Could not locate E0 and theta0 params in posterior distribution.')
                    

            ### Looping over parameters in the j:th path to localize parameter indeces in array parametrar
            for plot_list_i in range(len(i_plot_list)):#range(iteration_length[j]):
                param_mean[plot_list_i] = s['marginals'][plot_list_i]['median']
                onesigma[plot_list_i] = s['marginals'][plot_list_i]['1sigma']
                twosigma[plot_list_i] = s['marginals'][plot_list_i]['2sigma']
                threesigma[plot_list_i] = s['marginals'][plot_list_i]['3sigma']

                ### Correcting 1, 2, and 3 sigma
                if collimate_energy and i_plot_list[plot_list_i] == 'E0':
                    theta0_1sigma = s['marginals'][theta0_param]['1sigma']
                    theta0_2sigma = s['marginals'][theta0_param]['2sigma']
                    theta0_3sigma = s['marginals'][theta0_param]['3sigma']
                    
                    onesigma[E0_param] = np.log10(10**onesigma[plot_list_i] * (1-np.cos(theta0_1sigma)))
                    twosigma[E0_param] = np.log10(10**twosigma[plot_list_i] * (1-np.cos(theta0_2sigma)))
                    threesigma[E0_param] = np.log10(10**threesigma[plot_list_i] * (1-np.cos(theta0_3sigma)))


                x_grid = np.linspace(threesigma[plot_list_i,0] , threesigma[plot_list_i,1] , n_plot_marginal)#                x_grid = np.linspace(min(x_grid_total[:,plot_list_i]) , max(x_grid_total[:,plot_list_i]) , n_plot_marginal)


                
                ### index i points to the entry in the parametrar array corresponding to the parameter we are into in this stage of the for-loop
                try:
                    this_param = np.where(i_plot_list[plot_list_i] == list_of_plot_params)[0][0]

                    exec('plot_list_range_%d = np.append(plot_list_range_%d , this_param)'%(j,j))
                except:
                    print 'Skipping %s in %s/chains/'%(i_plot_list[plot_list_i] , pathToPlot[j])
                    continue



#                exec('plot_list_range_%d = np.append(plot_list_range_%d , this_param)'%(plot_list_i,plot_list_i))
                ### CONTINUE HERE!!! We must make sure that this_param points to the right index!!!! ###


            
                ### Creating average plots
                if sum_plots and (len(sum_average_paths) > 1): #Creating average x-and-y grid
                    n_ave = 100
                    x_grid = np.linspace(x_grid_lims[limCounter,0],x_grid_lims[limCounter,1],n_ave)
                    y_grid = np.zeros(n_ave)
                    y_counter = 0
                    for average_y in range(len(sum_average_paths)):
                        read_x_grid, read_y_grid = np.loadtxt('%s/temp_catalogue/grid_%d%d.txt'%(currentPlotPath,i,average_y))
                        for y_index in range(n_ave):
                            #Interpolating
                            nearest_under = np.argmin(np.abs(read_x_grid - x_grid[y_index])) #Finding the index in the read-in x_grid closest to the active x_grid element
                            nearest_under -= (read_x_grid[nearest_under] > x_grid[y_index]) #Making sure we have the index just under the active x_grid index
                            if nearest_under < 0: continue
                            elif nearest_under > (len(read_y_grid)-2): break
                            y_grid[y_index] += (read_y_grid[nearest_under] * (read_x_grid[nearest_under+1] - x_grid[y_index]) + read_y_grid[nearest_under+1] * (x_grid[y_index] - read_x_grid[nearest_under])) / (read_x_grid[nearest_under+1] - read_x_grid[nearest_under])  #Interpolation to find y_grid value on the new grid
                            
                        y_counter += 1
                    y_grid /= y_counter   #Averaging

                ### Reading chains files directly
                else:
                    ### Using kde_weights to estimate multiple gaussian fit to posterior

                    x_grid_total[:,plot_list_i]

                    cv = kde.Covariator(x_grid_total[:,plot_list_i] , post_in)
                    ic,norm = cv()
                    g_kde = kde.gaussian_kde(x_grid_total[:,plot_list_i] , post_in , ic , norm)
                    y_grid = g_kde(x_grid)*(x_grid[1] - x_grid[0])
                    

                exec('x_grid_%d_%d,y_grid_%d_%d = np.copy(x_grid) , np.copy(y_grid)'%(j,this_param,j,this_param))
                
                ### Storing mean values for order-of-ten correction of x_grid
                if preferredScale[whereParam_short[this_param]] == 'log':
                    x_grid_mean_temp = (10**x_grid[0] + 10**x_grid[-1]) / 2
                elif preferredScale[whereParam_short[this_param]] == 'deg':
                    x_grid_mean_temp = (x_grid[0] + x_grid[-1]) * 90 / np.pi
                else:
                    x_grid_mean_temp = (x_grid[0] + x_grid[-1]) / 2
                    
                ### Decides whether to store x_grid_mean_temp. Conditions to store is if it no previous value is stored (it is 0.), or if previously stored mean values are larger than the current ones
                if x_grid_mean[this_param] == 0. or x_grid_mean[this_param] > x_grid_mean_temp:
                    x_grid_mean[this_param] = np.copy(x_grid_mean_temp)



    xlabel_name = np.zeros(number_of_plot_params,dtype='S256')
    for j_store in range(numberOfPlots):
        exec('plot_list_range = np.copy(plot_list_range_%d)'%j_store)

        for i_store in plot_list_range:

            if list_of_plot_params[i_store] != str(None): ### If this is false, the parameter should be omitted
                ### Storing maximum y-values
                if y_max[i_store] == 0.: 
                    exec('y_max[%d] = max(y_grid_%d_%d)'%(i_store,j_store,i_store))
                    
                else:
                    exec('y_max_past = np.copy(y_max[%d])'%(i_store))
                    exec('y_max_here = max(y_grid_%d_%d)'%(j_store,i_store))


                    if y_max_here > y_max_past:
                        y_max[i_store] = np.copy(y_max_here)

                ### Setting xlabels
                if j_store == (numberOfPlots - 1): ### The last fileset
                    #preferredPlotScale[whereParam_short[i_store]] == 'lin': 
                    xlabel_name[i_store] = latexParamNamesLin[whereParam_short[i_store]]
                    #else: xlabel_name[i_store] = latexParamNames[whereParam_short[i_store]]

    ### Sorting plots after maximum y-value, highest first to present sharp distributions first. Plots with similar y-values should be on the same row in the plot panels

    plot_max_order = np.argsort(-y_max)

    if np.count_nonzero(y_max) > 0: plot_max_order = plot_max_order[: np.count_nonzero(y_max)]  ### numpy counts None as a zero when running count_nonzero. This line removes entries with None, i.e. entries user has choosen to omit

    ### plot_max_order lists the order of array 

    ### Setting subplot-grid
    
    list_subplots , first_plot = set_subplots(list_of_plot_params[plot_max_order],np.array(UseOp.preferredPlotScale)[whereParam_short[plot_max_order]] , xlabel_name[plot_max_order] , whereParam_short[plot_max_order] , len(UseOp.parametrar))

    ### Setting array with indeces pointing to last plot in each panel row
    last_plot = np.zeros(len(first_plot),dtype=int)
    if len(first_plot) > 1:
        last_plot[0] = first_plot[1]-1
        for lplot in range(1,len(first_plot)-1):
            last_plot[lplot] = first_plot[lplot+1] - 1
                    
    last_plot[-1] = number_of_plot_params - 1

    
        ### Plotting the posterior

    ### Plotting distributions

    for j_plot_all in range(numberOfPlots):
        for which_i_pa , i_plot_all in enumerate(plot_max_order):
            if whereParam_short[i_plot_all] == -1: 
                raw_input('We shouldn\'t be here!')
                continue
            
            exec('plot_list = list_of_plot_params_%d'%j_plot_all)
            try: 
                this_param_plot_list = np.where(np.array(paramNamesShort[whereParam_short[i_plot_all]]) == plot_list)[0]
                this_param = np.where(np.array(paramNamesShort[whereParam_short[i_plot_all]]) == list_of_plot_params)[0]
                if len(this_param_plot_list) != 1: 
                    continue
                ### this param states the index of the parameter in the particular chains files we are in, counted from zero. It is used when reading chains/1-stats.dat file
                this_param = this_param[0]
            except: 
                print 'stepping by'
                continue

            exec('x_grid = x_grid_%d_%d'%(j_plot_all , this_param))
            exec('y_grid = y_grid_%d_%d'%(j_plot_all , this_param))

            ### Correcting order-of-ten
            ### Reducing orders of magnitude of x_grid                                                                                                                          
            if not printInputLine:
                inputLine_position=None
            if UseOp.preferredPlotScale[whereParam_short[this_param]] != 'log':
                param_dim = preferredScale[whereParam_short[this_param]]
                
                if np.log10(x_grid_mean[this_param]) < -1 or np.log10(x_grid_mean[this_param]) > 1:
                    correct_mag_factor = np.floor(np.log10(x_grid_mean[this_param]))
                    if param_dim == 'log':
                        x_grid -= correct_mag_factor
                        paramLimits_corr = UseOp.paramLimits[whereParam_short[i_plot_all]] - correct_mag_factor
                        if printInputLine:
                            inputLine_position = ModVar.echo_value(whereParam_short[i_plot_all]) - correct_mag_factor

                    else:
                        x_grid *= 10**-correct_mag_factor
                        paramLimits_corr = UseOp.paramLimits[whereParam_short[i_plot_all]] * 10** (-correct_mag_factor)
                        if printInputLine:
                            inputLine_position = ModVar.echo_value(whereParam_short[i_plot_all]) * 10**(-correct_mag_factor)


                    if compare_results and j_plot_all == 0:
                        compare_input[:,:,this_param] *= 10**-correct_mag_factor


                    ### Correcting label, if this is the last set of plots
                    if j_plot_all == numberOfPlots - 1:
                        correct_xlabel = list_subplots[whereParam_short[i_plot_all]].get_xlabel()
                        correct_xlabel = correct_xlabel[:-1] + '\\times 10^{%d}'%(-correct_mag_factor) + '$'
                        list_subplots[whereParam_short[i_plot_all]].set_xlabel(correct_xlabel)
                else:
                    paramLimits_corr = UseOp.paramLimits[whereParam_short[i_plot_all]]
                    if printInputLine:
                        inputLine_position = ModVar.echo_value(whereParam_short[i_plot_all])

                    
            else:
                paramLimits_corr = UseOp.paramLimits[whereParam_short[i_plot_all]]
                if printInputLine:
                    inputLine_position = ModVar.echo_value(whereParam_short[i_plot_all])
                        
            first_plot_set = (j_plot_all==0)
            if not first_plot_set: #first_plot_set can be true if this parameter has not been plotted yet
                first_plot_set = True
                for i_test_first in range(j_plot_all):
                    exec('temp_param_list = list_of_plot_params_%d'%i_test_first)
                    if np.sum(np.where(plot_list[plot_max_order[i_plot_all]] == temp_param_list)) != 0: ### Testing if this parameter has been plotted yet. If not, first_plot_set is True
                        first_plot_set = False
                        break

            
            if plot_multinest_gaussian[j_plot_all]: min_xvalue_in,max_xvalue_in = stats_subplot(list_subplots[whereParam_short[i_plot_all]],x_grid,y_grid,plotLine[j_plot_all],plot_multinest_gaussian[j_plot_all],preferredPlotScale[whereParam_short[i_plot_all]],preferredScale[whereParam_short[i_plot_all]],first_plot_set,(j_plot_all==(numberOfPlots-1)) and printInputLine,numberOfPlots,this_param_plot_list,mstats_text,'%s%s'%(gauss_color[j_plot_all],gauss_linetype[j_plot_all]),inputLine_position,midLineColor,paramLimits_corr)


            else: min_xvalue_in,max_xvalue_in = stats_subplot(list_subplots[whereParam_short[i_plot_all]],x_grid,y_grid,plotLine[j_plot_all],plot_multinest_gaussian[j_plot_all],preferredPlotScale[whereParam_short[i_plot_all]],preferredScale[whereParam_short[i_plot_all]],first_plot_set,(j_plot_all==(numberOfPlots-1)) and printInputLine,numberOfPlots,this_param_plot_list,mstats_text,j_plot_all,None,inputLine_position,midLineColor,paramLimits_corr)

            if j_plot_all == 0: ### First plot, extreme values are assigned as plot limit material
                min_xvalue[i_plot_all] = np.copy(min_xvalue_in)
                max_xvalue[i_plot_all] = np.copy(max_xvalue_in)
            else:
                if min_xvalue[i_plot_all] > min_xvalue_in: min_xvalue[i_plot_all] = np.copy(min_xvalue_in)
                if max_xvalue[i_plot_all] < max_xvalue_in: max_xvalue[i_plot_all] = np.copy(max_xvalue_in)
            


            ### Plotting input parameters, used to compare results
            if compare_results:   #Plot input fit parameters
                
                xlow_in , xhigh_in = plot_comres(list_subplots[whereParam_short[i_plot_all]],compare_input[:,:,i_plot_all],comres_color,UseOp.preferredPlotScale[whereParam_short[i_plot_all]])
                if xlow_in < min_xvalue[i_plot_all]:
                    min_xvalue[i_plot_all] = np.copy(xlow_in)
                if xhigh_in > max_xvalue[i_plot_all]:
                    max_xvalue[i_plot_all] = np.copy(xhigh_in)
                
    
                limCounter += 1



            if j_plot_all == (numberOfPlots - 1): ### Last plot
                #############
                ### Ticks ###
                #############

                if UseOp.preferredPlotScale[whereParam_short[i_plot_all]] == 'log':
                    if min_xvalue[i_plot_all] <= 0: 
                        lower_xlim = 0.
                    else:
                        lower_xlim = 10**(np.floor(np.log10(min_xvalue[i_plot_all])))

                    if max_xvalue[i_plot_all] <= 0: 
                        upper_xlim = 0.
                    else:
                        upper_xlim = 10**(np.floor(np.log10(max_xvalue[i_plot_all])) + 1)


                ### making sure the x-limit isn't outside prior range
                if preferredScale[whereParam_short[i_plot_all]] == 'log':

                    if lower_xlim*10 < 10**UseOp.paramLimits[whereParam_short[i_plot_all]][0]:
                        lower_xlim = 10**UseOp.paramLimits[whereParam_short[i_plot_all]][0]
                    if upper_xlim > 10**UseOp.paramLimits[whereParam_short[i_plot_all]][1]:
                        upper_xlim = 10**UseOp.paramLimits[whereParam_short[i_plot_all]][1]


                if UseOp.preferredPlotScale[whereParam_short[i_plot_all]] == 'log':
                    min_log_value = 10**np.floor(np.log10(min_xvalue_in))
                    max_log_value = 10**(np.floor(np.log10(max_xvalue_in))+1)
                    if min_log_value < lower_xlim * 10:
                        lower_xlim = min_log_value/10
                    if max_log_value > upper_xlim:
                        upper_xlim = max_log_value

                    list_subplots[whereParam_short[i_plot_all]].set_xlim([lower_xlim*10,upper_xlim])
                    list_subplots[whereParam_short[i_plot_all]].axes.tick_params(axis='x',pad=10)

                    ### If too many ticks, reducing by half
                    in_log_ticks = list_subplots[whereParam_short[i_plot_all]].get_xticks()
                    new_log_ticks = False
                    while len(in_log_ticks) > 4:
                        in_log_ticks = in_log_ticks[::2]
                        new_log_ticks = True

                    if not np.sum(last_plot==which_i_pa)>0: ### If not last plot in panel
                        in_log_range = list_subplots[whereParam_short[i_plot_all]].get_xlim()

                        while in_log_ticks[-1] >= in_log_range[-1]:
                            in_log_ticks = in_log_ticks[:-1]
                            new_log_ticks = True
                    if new_log_ticks:
                        list_subplots[whereParam_short[i_plot_all]].set_xticks(in_log_ticks)

                    list_subplots[whereParam_short[i_plot_all]].set_xlim([lower_xlim*10,upper_xlim])
                else:
                    lower_xlim = min_xvalue[i_plot_all]
                    upper_xlim = max_xvalue[i_plot_all]


                    reduce_ticks(list_subplots[whereParam_short[i_plot_all]],np.sum(first_plot==which_i_pa)>0,np.sum(last_plot==which_i_pa)>0,lower_xlim,upper_xlim)


                

                ##############
                ### Limits ###
                ##############

                    
                


    


                
    #Deleting temp_catalogue
    if create_temp_catalogue: os.system('rm %s/temp_catalogue -r'%currentPlotPath)
    os.chdir(currentPlotPath)

    
    plt.tight_layout()
    if os.path.isfile('%s/probability_%s.%s'%(pathName,UseOp.GRBlabel,UseOp.figTypes)):
        del_fig_file = raw_input('Overwrite probability_%s.%s? ([y]/n): '%(UseOp.GRBlabel,UseOp.figTypes))
        if (del_fig_file == 'n') or (del_fig_file == 'N'): 
            output_filename = raw_input('Type file name without suffix: ')
        else:
            output_filename = 'probability_%s'%(UseOp.GRBlabel)
    else:
        output_filename = 'probability_%s'%(UseOp.GRBlabel)

    plt.savefig('%s/%s.%s'%(pathName,output_filename,UseOp.figTypes))
    plt.savefig('%s/%s.pdf'%(pathName,output_filename))
    plt.close()
    #os.system('pdftk %s/*.pdf cat output %s/merged.pdf'%(pathName,pathName))
    if UseOp.allowPrint: print "Probablility plots saved in %s/%s.%s"%(pathName,output_filename,UseOp.figTypes)



### This function gathers the information from the chains/1-stats.dat file and prints an analysis of them compared to the input value
def gather_multinest_stats():
    stats_path_input = raw_input('Path where to find the chains folder containing the 1-stats.dat file. Type general path using asterix, separated by space to analyse many or a file containing a list of all paths separated by newlines: ')
    write_tex_files_input = raw_input('Write tex-files? ([y]/n): ')
    write_tex_files = (write_tex_files_input == '') or (write_tex_files_input == 'y') or (write_tex_files_input == 'Y')
    if not write_tex_files:
        read_output_matrix_input = raw_input('Analyse written files? ([y]/n): ')
        read_output_matrix = (read_output_matrix_input == '') or (read_output_matrix_input == 'Y') or (read_output_matrix_input == 'y')
    
    if len(stats_path_input.split()) > 1: ### Many inputs
        stats_path = np.zeros(0,dtype='S256')
        for i in range(len(stats_path_input.split())):
            stats_path = np.append(stats_path,glob.glob('%s/chains/1-stats.dat'%stats_path_input.split()[i]))
    else:
       try: 
           open_stats_path = open(stats_path_input,'r')
           stats_path = open_stats_path.read().split('\n')
           open_stats_path.close()
           while stats_path[-1] == '': stats_path = stats_path[:-1]
       except:
           stats_path = glob.glob('%s/chains/1-stats.dat'%stats_path_input)
    stats_path_to_chains = np.zeros(np.shape(stats_path),dtype='S256')

    true_input = ModVar.echo_value(whereParam)

    lower_limit , upper_limit = np.zeros([len(stats_path),n_params],dtype=bool) , np.zeros([len(stats_path),n_params],dtype=bool)

    dist_type = np.zeros((len(stats_path),n_params),dtype='S40')
    type_of_distribution = np.zeros((len(stats_path),n_params),dtype=int)
    
    stats_gaussian = np.zeros([len(stats_path),n_params,2])
    stats_max_like = np.zeros([len(stats_path),n_params])
#    stats_MAP = np.zeros([len(stats_path),n_params,2])
        
    sigma_percentage = np.zeros([len(stats_path),n_params])
    sigma_percentage_text = np.zeros([len(stats_path),n_params],dtype='S20')
    true_input_array = np.ones([len(stats_path),n_params]) * true_input
    mu_diff = np.zeros([len(stats_path),n_params])
    mu_diff_sign = np.zeros([len(stats_path),n_params],dtype=bool)

    bestfit_diff = np.zeros([len(stats_path),n_params])
    bestfit_diff_sign = np.zeros([len(stats_path),n_params],dtype=bool)
    same_side_sign = np.zeros([len(stats_path),n_params],dtype=bool)

    diff_sign_text = np.zeros(np.shape(same_side_sign),dtype='S20')
    inside_sigma = np.zeros([len(stats_path),n_params],dtype=bool)
    inside_sigma_text = np.zeros(np.shape(inside_sigma),dtype='S20')

#    write_cover_file_input = raw_input('Create cover file? ([y]/n): ')
#    write_cover_file = write_cover_file_input == '' or write_cover_file_input == 'y' or write_cover_file_input == 'Y'
#    if write_cover_file: cover_text = ''

    #if read_output_matrix:
        ### Creating a 3-D tensor with dimensions Number of files x Number of entries x Number of parameters
    #    input_data = np.zeros([1,6,n_params])

    for i_path in range(len(stats_path)):
        if write_tex_files:
            outtext = '\\begin{tabular}{l|c|c|c|c|c|c|c|c}\nParameter & Dist. type & T.I. loc. & Max. Like. & $\mu$ & T.I. & No. S.Ds. & Bias & Size of $\sigma$\\\\\n'
        stats_path_to_chains[i_path] = os.path.abspath('%s'%('/'.join(stats_path[i_path].split('/')[:-2])))

        if write_tex_files: os.chdir(stats_path_to_chains[i_path])

        path_name_split = os.path.abspath(stats_path_to_chains[i_path]).split('/')
        if len(path_name_split[-1]) > 4: path_name_split[-1] = path_name_split[-1][:4]


        if write_tex_files:
            output_filename = raw_input('Output file name without suffix [%s_%s_%s.tex]: '%(path_name_split[-3],path_name_split[-2],path_name_split[-1]))
        else:
            output_filename = '%s_%s_%s'%(path_name_split[-3],path_name_split[-2],path_name_split[-1])
        if output_filename == '': output_filename = '%s_%s_%s'%(path_name_split[-3],path_name_split[-2],path_name_split[-1])
        if os.path.isfile('%s.tex'%output_filename) and write_tex_files:
            ask_del_input = raw_input('Delete existing file %s.tex? ([y]/n): '%output_filename)
            ask_del = ask_del_input == 'y' or ask_del_input == 'Y' or ask_del_input == ''
            if not ask_del:
                print 'Will not overwrite file %s.tex. Now exiting.'%output_filename
                raise SystemExit(0)
            else: os.system('rm  %s.tex'%output_filename)


        ####################################
        ### Loading written matrix files ###
        ####################################

        if read_output_matrix:
            input_matrix_files = glob.glob('%s/%s_mode*.dat'%(stats_path_to_chains[i_path].split('chains/1-stats.dat')[0] , output_filename))

            for open_file_modes in range(len(input_matrix_files)):  ### Looping over all modes
                try:
                    input_data_temp = np.copy(input_data)
                    input_data = np.zeros([np.shape(input_data)[0]+1,np.shape(input_data)[1],np.shape(input_data)[2]])
                    input_data[-1] = np.loadtxt(input_matrix_files[open_file_modes])
                    input_data_path = np.append(input_data_path,'%s'%input_matrix_files[open_file_modes])
                    input_data_model_list = np.append(input_data_model_list,np.array(path_name_split[-1]))
                except: 
                    input_data = np.zeros([1,6,n_params])
                    input_data_path = np.zeros(1,dtype='S256')
                    input_data_path[0] = '%s'%input_matrix_files[open_file_modes]
                    input_data_model_list = np.array([path_name_split[-1]])
                    try:input_data[0] = np.loadtxt(input_matrix_files[open_file_modes])
                    except:
                        print 'Wrong number of parameters specified in options.py. Now exiting.'
                        raise SystemExit(0)
            
            

            continue

        ######################################
        ### Loading posterior distribution ###
        ######################################

        
        try:
            a = pymultinest.Analyzer(n_params = n_params)
            s = a.get_stats()
            p = pymultinest.PlotMarginalModes(a)
        except:
            print 'Could not read files in %s. Now exiting'%stats_path_to_chains[i_path]
            raise SystemExit(0)

        for read_x_grid in range(n_params):

            x_grid,y_grid = p.plot_marginal(read_x_grid, with_ellipses = True, with_points = False, grid_points=50)
            x_grid_selection = np.where((x_grid>UseOp.paramLimits[whereParam[read_x_grid],0]) & (x_grid<UseOp.paramLimits[whereParam[read_x_grid],1])) #Reducing grid to only contain points within user specified prior range
            x_grid,y_grid = x_grid[x_grid_selection] , y_grid[x_grid_selection]


            ### Detecting lower and upper limits ;  type_of_distribution: see note from 12/8 -15
            
            if y_grid[0] > (max(y_grid)*0.5): 
                upper_limit[i_path,read_x_grid] = True
                dist_type[i_path,read_x_grid] = 'UL'
                type_of_distribution[i_path,read_x_grid] = 2
            elif y_grid[0] > (max(y_grid)*0.2):
                upper_limit[i_path,read_x_grid] = True
                dist_type[i_path,read_x_grid] = 'PUL'
                type_of_distribution[i_path,read_x_grid] = 3
            if y_grid[-1] > (max(y_grid)*0.5): 
                lower_limit[i_path,read_x_grid] = True
                dist_type[i_path,read_x_grid] = 'LL'
                type_of_distribution[i_path,read_x_grid] = 4
            elif y_grid[-1] > (max(y_grid)*0.2): 
                lower_limit[i_path,read_x_grid] = True
                dist_type[i_path,read_x_grid] = 'PLL'
                type_of_distribution[i_path,read_x_grid] = 5
            if upper_limit[i_path,read_x_grid] == False and lower_limit[i_path,read_x_grid] == False:                 
                dist_type[i_path,read_x_grid] = '{\color{green} CD }'
                type_of_distribution[i_path,read_x_grid] = 0
            if upper_limit[i_path,read_x_grid] and lower_limit[i_path,read_x_grid]: 
                dist_type[i_path,read_x_grid] = '{\\color{red} UC}'
                type_of_distribution[i_path,read_x_grid] = 1
                





        ### Loading multinest outputs from chains/1-stats.dat
        os.chdir(currentPath)
        open_stats = open(os.path.abspath(stats_path[i_path]))
        read_stats = open_stats.read()
        open_stats.close()

        number_of_modes = int(read_stats.split('Total Modes Found:')[1].split()[0])
        print '\n%s has %d mode%s'%(os.path.abspath(stats_path[i_path]),number_of_modes,'s'*(number_of_modes>1))

        

        ### Reading in the stats ###
        local_log_evidence = np.zeros([number_of_modes,n_params])
        global_log_evidence = np.zeros(n_params)

        global_log_evidence[0] = float(read_stats.split('+/-')[0].split()[-1])
        global_log_evidence[1] = float(read_stats.split('+/-')[1].split()[0])

        #if number_of_modes == 1:
        for i_modes in range(number_of_modes):
            if number_of_modes > 1 and write_tex_files:
                outtext += 'mode %d &&&&&\\\\\n'%(i_modes+1)
            stats_split = read_stats.split('Dim No.')
            local_log_evidence[i_modes,0] = float(read_stats.split('Local Log-Evidence')[i_modes+1].split()[0])
            local_log_evidence[i_modes,1] = float(read_stats.split('Local Log-Evidence')[i_modes+1].split()[2])
            
            

                                                                                                           
                                                
#            print stats_split
            for i_stats in range(n_params):
                
                
                stats_gaussian[i_path,i_stats] = map(float,stats_split[3*i_modes+1].split('\n')[i_stats+1].split())[1:3]
                
#                print map(float,stats_split[1].split('\n')[i_stats+1].split())[1:3]
                
                stats_max_like[i_path,i_stats] = map(float,stats_split[3*i_modes+2].split('\n')[i_stats+1].split())[1]
#                stats_MAP[i_path,i_stats] = map(float,stats_split[number_of_modes*i_modes+3].split('\n')[i_stats+1].split())
#            if np.sum(stats_max_like==stats_MAP) == np.size(stats_max_like):
#                print ''#\nMaximum Likelihood Parameters equal the MAP Parameters'
#            else:
#                print ''#\n\n\nMaximum Likelihood Paramters does not equal the MAP Parameters! What to do now?\n\n\n'
                
                ### Percentage of the sigma/mu 
                ### ### For log-scale, print how many decades
                ### ### For lin-scale, print the percentage
        
                if preferredScale[whereParam[i_stats]] == 'lin' or preferredScale[whereParam[i_stats]] == 'deg': #Parameter is fitted on its linear
                    if preferredScale[whereParam[i_stats]] == 'deg': value_factor = 180/np.pi  #Factor to write out degrees in output table
                    else: value_factor = 1.
                    maxlike_exp = int(np.log10(stats_max_like[i_path,i_stats]*value_factor)) * ((np.log10(stats_max_like[i_path,i_stats]*value_factor) > 3) or (np.log10(stats_max_like[i_path,i_stats]*value_factor) < -3))
                    maxlike_base = 10**(np.log10(stats_max_like[i_path,i_stats]*value_factor) - maxlike_exp)  #Gives power-of-ten base if too large or too small to write out
                    
                    mu_exp = int(np.log10(stats_gaussian[i_path,i_stats,0]*value_factor)) * ((np.log10(stats_gaussian[i_path,i_stats,0]*value_factor) > 3) or (np.log10(stats_gaussian[i_path,i_stats,0]*value_factor) < -3))
                    mu_base = 10**(np.log10(stats_gaussian[i_path,i_stats,0]*value_factor)-mu_exp)  #Gives power-of-ten base if too large or too small to write out

                    trueinput_exp = int(np.log10(true_input[i_stats]*value_factor)) * ((np.log10(true_input[i_stats]*value_factor) > 3) or (np.log10(true_input[i_stats]*value_factor) < -3))
                    trueinput_base = 10**(np.log10(true_input[i_stats]*value_factor) - trueinput_exp)


                    if stats_gaussian[i_path,i_stats,1] == 0.: sigma_percentage_text[i_path,i_stats] = 'inf'
                    else:
                        sigma_percentage[i_path,i_stats] = stats_gaussian[i_path,i_stats,1] / stats_gaussian[i_path,i_stats,0] * 100
                        sigma_percentage_text[i_path,i_stats] = '%s \\%% '%(round_off(sigma_percentage[i_path,i_stats]*100,2))
#                elif preferredScale[whereParam[i_stats]] == 'deg':
#                    if stats_gaussian[i_path,i_stats,1] == 0.: sigma_percentage_text[i_path,i_stats] = 'inf'
#                    else:
#                        sigma_percentage_text[i_path,i_stats] = '%s $^{\circ}$'%(round_off( stats_gaussian[i_path,i_stats,1]*180/np.pi ,2))
                elif preferredScale[whereParam[i_stats]] == 'log':
                    maxlike_exp = int(stats_max_like[i_path,i_stats]) * ((stats_max_like[i_path,i_stats] > 3) or (stats_max_like[i_path,i_stats] < -3))
                    maxlike_base = 10**(stats_max_like[i_path,i_stats] - maxlike_exp)  #Gives power-of-ten base if too large or too small to write out
                    
                    mu_exp = int(stats_gaussian[i_path,i_stats,0]) * ((stats_gaussian[i_path,i_stats,0] > 3) or (stats_gaussian[i_path,i_stats,0] < -3))
                    mu_base = 10**(stats_gaussian[i_path,i_stats,0]-mu_exp)  #Gives power-of-ten base if too large or too small to write out

                    trueinput_exp = int(true_input[i_stats]) * ((true_input[i_stats] > 3) or (true_input[i_stats] < -3))
                    trueinput_base = 10**(true_input[i_stats] - trueinput_exp)

                    if stats_gaussian[i_path,i_stats,1] == 0.: sigma_percentage_text[i_path,i_stats] = 'inf'
                    else:
                        sigma_percentage[i_path,i_stats] = stats_gaussian[i_path,i_stats,1]   #How many decades does the standard deviation cover?
                        sigma_percentage_text[i_path,i_stats] = '%s dec'%(((-1)**(sigma_percentage[i_path,i_stats]<0))*round_off(np.abs(sigma_percentage[i_path,i_stats]),2))
                else:
                    print 'Scale is not defined for parameter %s. Now exiting.'%(paramNames[whereParam[i_stats]])
                    raise SystemExit(0)




#        else:
#            raw_input('\n\n\nMore than one mode nb nb nb!!!! Now write the code to analyse this file!\n\n\n')
                             
    ####################################################
    ### Comparing multinest output to the true input ###
    ####################################################

    

    ### Difference between true input and mean value. Negative - over estimated; positive - under estimated
                mu_diff[i_path,i_stats] = true_input[i_stats] - stats_gaussian[i_path,i_stats,0]
                mu_diff_sign[i_path,i_stats] = mu_diff[i_path,i_stats] > 0
    ### Difference between true input and the best-fit value. Negative - over estimated; positive - under estimated
                bestfit_diff[i_path,i_stats] = true_input[i_stats] - stats_max_like[i_path,i_stats]
                bestfit_diff_sign[i_path,i_stats] = bestfit_diff[i_path,i_stats] > 0

                same_side_sign[i_path,i_stats] = mu_diff_sign[i_path,i_stats] == bestfit_diff_sign[i_path,i_stats] #Are mu and the best-fit on the same side of the true input?

    

#                where_diff_sign[i_path,i_stats] = np.where(same_side_sign[i_path,i_stats])
#                where_diff_sign_false = np.where(same_side_sign==False)

                if same_side_sign[i_path,i_stats]:
#                for i_bestfit in range(len(where_diff_sign[0])):
                    diff_sign_text[i_path,i_stats] = mu_diff_sign[i_path,i_stats] * '{\\color{red}UE}' + (mu_diff_sign[i_path,i_stats]==False) * '{\\color{blue}OE}'
                else:
# i_bestfit_false in range(len(where_diff_sign_false[0])):

                    diff_sign_text[i_path,i_stats] = '$\mu$%sb%s'%(mu_diff_sign[i_path,i_stats] * 'U' + (mu_diff_sign[i_path,i_stats]==False) * 'O'  ,  bestfit_diff_sign[i_path,i_stats] * 'U' + (bestfit_diff_sign[i_path,i_stats]==False) * 'O')
    
    ### Is the true input inside the standard deviation?
                inside_sigma[i_path,i_stats] = np.array((true_input[i_stats] > (stats_gaussian[i_path,i_stats,0] - stats_gaussian[i_path,i_stats,1])) & (true_input[i_stats] < (stats_gaussian[i_path,i_stats,0] + stats_gaussian[i_path,i_stats,1])))


                if inside_sigma[i_path,i_stats]: inside_sigma_text[i_path,i_stats] = '{\color{green} IS}'
                else: inside_sigma_text[i_path,i_stats] = '{\color{red} OS}'

 
   


    
    #######################
    ### Printing output ###
    #######################


                no_s_d_ =  round_off(np.abs(mu_diff[i_path,i_stats]),2)
                
                maxlike_text = '$%s%s%s$'%(round_off(maxlike_base,2) , (' \\times 10^{%d}'%maxlike_exp)*(maxlike_exp != 0) , ' ^{\circ}'*(preferredScale[whereParam[i_stats]]=='deg')) #Only writes out the power of ten if it is needed
                
                mu_text = '$%s%s%s$'%(round_off(mu_base,2) , (' \\times 10^{%d}'%mu_exp)*(mu_exp!=0) , ' ^{\circ}'*(preferredScale[whereParam[i_stats]]=='deg')) #Only writes out the power of ten if it is needed
                
                trueinput_text = '$%s%s%s$'%(round_off(trueinput_base,2) , (' \\times 10^{%d}'%trueinput_exp)*(trueinput_exp!=0) , ' ^{\circ}'*(preferredScale[whereParam[i_stats]]=='deg')) #Only writes out the power of ten if it is needed


                if write_tex_files:
                    outtext += '%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\\n'%(latexParamNames[whereParam[i_stats]],dist_type[i_path,i_stats], inside_sigma_text[i_path,i_stats] ,maxlike_text , mu_text , trueinput_text  , no_s_d_ , diff_sign_text[i_path,i_stats],sigma_percentage_text[i_path,i_stats])
             ############################
             ### Writing matrix files ###
             ############################
                
                #######
                ### matrix layout: mu, sigma, maximum likelihood, type of distribution (see note from 12/8 -15
                #######

                ### Saves file in the same directory as where the chains folder is found
                
            out_array = np.zeros([6,n_params])
            out_array[0] = stats_gaussian[i_path,:,0]
            out_array[1] = stats_gaussian[i_path,:,1]
            out_array[2] = stats_max_like[i_path]
            out_array[3] = type_of_distribution[i_path]
            out_array[4] = local_log_evidence[i_modes]
            out_array[5] = global_log_evidence
            np.savetxt('%s/%s_mode_%d.dat'%(stats_path_to_chains[i_path].split('chains/1-stats.dat')[0], output_filename , i_modes), out_array )
        if write_tex_files:
            outtext += '\\end{tabular}'
        
            write_output = open('%s.tex'%(output_filename),'w')
            write_output.write(outtext)
            write_output.close()



    #######################
    ### Analysing input ###
    #######################

    if read_output_matrix:
        print np.shape(input_data)
        print input_data_path
        print len(input_data_model_list)
        ### Degrees of Freedom for each type of coverage
        dict = {'full_log':162 , 'full_inconsistent':162 , 'optical':152 , 'xrays':12 , 'radio':138 , 'prebreak':152 , 'postbreak':152}
        
        true_input_list = [[51.,53.],[-2.,0.],[5.,15.],[1.,7.]]

        ########################
        ### Analysing biases ###
        ########################

        bias_after_model = np.zeros([7,32])

        for i_model in range(2):
            for j_model in range(2):
                for k_model in range(2):
                    for l_model in range(2):
                        for s_model in [0,2]:
                            ### Order in input_data: mu, sigma, maximum likelihood, type of distribution, local log-evidence (central value and 1 sigma), global log-evidence (central value and 1 sigma)
                            
                            blabla
                            ###################
####################3

#CONTINUE HERE!!!!
        
        


###############################################
#         Beginning of programme              #
###############################################  

#Importing user options and constants from options.py

UseOp = userOptions()

### Issuing warnings to avoid confusion in turned-off switches
if np.sum((UseOp.radiativeLosses==False)  +  UseOp.fixedRSFSratio + UseOp.fixed_epsilon) != 0 or (np.sum(UseOp.parametrar)>= 12 and not UseOp.reverseShock):
    warning_out = '!!! Warning! Following switches are on/off: !!!'
    print '-'*len(warning_out)
    print warning_out
    print '-'*len(warning_out)
    print '\n\n\n'

    if not UseOp.radiativeLosses:
        print '- radiative losses are turned off'
    if UseOp.fixedRSFSratio:
        print '- RS microphysical parameters are fixed to FS values!'
    if UseOp.fixed_epsilon:
        print '- Dynamics module is assuming a fixed radiative loss fraction'
    if (np.sum(UseOp.parametrar)>= 12 and not UseOp.reverseShock):
        print '- there are %d free parameters, but RS is turned off!'%np.sum(UseOp.parametrar)
    print '\n\n\n'

ModVar = model_variables()
PlotDetails = plot_details()



#constants = np.array(options.inputConstants()) 

### Correcting boolean vector parametrar to automatically prior range for the number density n, depending if CM or WM model
   ### Variable s is constants[8], density prior range is parametrar[6]. Now adding another boolean into parametrar[6]
if UseOp.parametrar[6]:
    paramet_tail = np.copy(UseOp.parametrar[7:])
    if ModVar.s == 0: ### CM
        UseOp.parametrar = np.concatenate([UseOp.parametrar[:6],[True,False],UseOp.parametrar[7:]])
        UseOp.parametrar[9] = False
    else:
        UseOp.parametrar = np.concatenate([UseOp.parametrar[:6],[False,True],UseOp.parametrar[7:]])
        


if UseOp.allowPrint: print "Program started %dh %dm %ds"%(tidPre[3:6])

n_params = int(np.sum(UseOp.parametrar))
currentPath = os.path.abspath('.')
abovePath = os.path.abspath('../')
homePath = os.path.expanduser('~')
whereParam = np.where(UseOp.parametrar)[0]
n_plot_marginal = 50

surfRingsOut = 300
numRuns = 0
loadInputConstants = True #When this is true, the options.py file is loaded for constants. If not, the parameters.txt file is loaded
printPlot = False
plot_area = False
plot_SED = False
paramNames = np.array(['log10(epsilon_rad)','log10(epsilon_rad_RS)','log10(epsilone_FS)','log10(epsilonp_FS)','log10(epsilone_RS)','log10(epsilonp_RS)','log10(E0)','log10(n_CM)','log10(A0)','s','R_ISM','log10(Gamma0)','log10(epsilonB_FS)','log10(epsilonB_RS)','p_FS','t0','theta0','alpha','t_prompt','p_RS','z'])
latexParamNames = [r'$\log_{10}(\epsilon_{\rm rad})$',r'$\log_{10}(\epsilon_{\rm e})$',r'$\log_{10}(\epsilon_{\rm p})$',r'$\log_{10}(\epsilon_{\rm e,RS})$',r'$\log_{10}(\epsilon_{\rm p,RS})$',r'$\log_{10}(E_0)$',r'$\log_{10}(n_{\rm CM})$',r'$\log_{10}(A_0) [{\rm cm}^{-3+s}]$',r'$s$',r'$log_{10}(R_{\rm ISM})$',r'$\log_{10}(\Gamma_0)$',r'$\log_{10}(\epsilon_{\rm B})$',r'$\log_{10}(\epsilon_{\rm B,RS})$',r'$p_{\rm FS}$',r'$t_0$',r'$\theta_0$',r'$\alpha$',r'$\Delta t_{\rm of}$',r'$p_{\rm RS}$',r'$z$']
latexParamNamesLin = [r'$\epsilon_{\rm rad}$',r'$\epsilon_{\rm e}$',r'$\epsilon_{\rm p}$',r'$\epsilon_{\rm e,RS}$',r'$\epsilon_{\rm p,RS}$',r'$E_{\rm 0,iso}$',r'$n_{\rm CM}$',r'$A_0$',r'$s$',r'$R_{\rm ISM}$',r'$\Gamma_0$',r'$\epsilon_{\rm B}$',r'$\epsilon_{\rm B,RS}$',r'$p_{\rm FS}$',r'$t_0$',r'$\theta_0 (^{\circ})$',r'$\alpha (^{\circ})$',r'$\Delta t_{\rm of}$',r'$p_{\rm RS}$',r'$z$']
paramNamesShort = np.array(['epsilon','epsilon_e','epsilon_p','epsilon_e_RS','epsilon_p_RS','E0','n_CM','A_0','s','R_ISM','Gamma0','epsilon_B','epsilon_B_RS','p','t0','theta0','alpha','t_prompt','p_RS','z'])
#                            e_rad e_e2  e_p2  e_e3  e_p3  E0    n_CM   A0    s    R_ISM Gam0  e_B2  e_B3   p2       theta0 alpha tprompt p3  
preferredScale = np.array(['log','log','log','log','log','log','log','log','lin','log','log','log','log','lin','None','deg','deg','log','lin','lin'])
constantsLength = len(paramNames) #The number of elements in contants list



### Handling input arguments ###


nArgs = len(sys.argv)
Plot_Exceptions = exception_class()

bin_data_bool = False
if nArgs > 1: 
    inArg = sys.argv
    
    #Help input argument
    if (inArg[1] == "h") | (inArg[1] == "help"):
        print "This is the help section, called by commandline arguments h or help.\n\n--------------------------------\nAvailable commandline arguments: \n\nburst= Type the date of the GRB\noption= \n    LC: lightcurve production\n    fit: run the fitting routine\n    marginal: plot 2D marginal plots\n    plot-marginal: advanced marginal plotter. Lets the user plot many posteriors on top of each other\n    print-stats: Produce a LaTeX table with stats from MultiNest output\n    read-stats: Produce a LaTeX table with the properties of the probability distributions. LaTeX output requires package \'color\'\n    the default option is set in options.py, named runOption=''\n\nload=\n    output: load constants from the parameters.txt file, printed by the fitting routine\n    chains: Load best-fit parameters from MultiNest output\n    Default is to load the parameters from the options.py file\n\nexcept=\n    FS: Plot the forward-shock only\n    RS: Plot the reverse-shock component only\n    thermal: Plot the thermal component only\n    These above may be combined by calling except= many times. Option expept is best suited when printing lightcurve output files (print-write). Save each component in a specific file, and load all files when running print-read. See manual for mor details, by entering commandline argument man. This option is still under construction\n\nbin=\n    Enter number of data points per bin. Default is no binning.\n\nband=\n    xrays/x-rays: Load X-ray lightcurves >0.1 keV (2.4e16 Hz)\n    UV: Load UV lightcurves 3.2 eV - 0.1 keV (7.7e14 - 2.4e16 Hz)\n    UV-optical/optical-UV: Load optical and UV lightcurves 1.6eV - 3.2eV\n    optical: Load optical lightcurves 0.4eV - 4eV\n    radio-mm: Load radio-mm lightcurves < 1.6eV\n    radio: Loading radio lightcurves (<600 GHz)\n    Manual input in GHz, range limits separated with a \'-\'. For example band=1e5-1e8 for the range 1e5 to 1e8 GHz\n\n\nProgramme is called by typing \'(i)python mainFit.py ->commandline arguments<-\'\n\nFor other questions, mail Andreas Johansson on johansson.mcquack@gmail.com\nNow exiting.  "
        raise SystemExit(0)
    elif inArg[1] == 'man':
        print 'There is no manual pages writted yet. Please contact Andreas Johansson via johansson.mcquack@gmail.com\nNow exiting.'
        raise SystemExit(0)
    #Other input arguments
    for i in range(1,nArgs):
        inputArg = inArg[i].split('=')
        #Setting burst name

        if inputArg[0] == 'except':
            if inputArg[1] == 'FS': Plot_Exceptions.FS_only = True
            elif inputArg[1] == 'RS': Plot_Exceptions.RS_only = True
            elif inputArg[1] == 'thermal': Plot_Exceptions.thermal_only = True
            else: print 'No except I can understand!\n\nPossible expects:\n\'FS\' - plot forward shock only\n\'RS\' - plot reverse shock only'

        elif inputArg[0] == 'bin':
            try: 
                bin_data = int(inputArg[1])
                bin_data_bool = True  #Boolean variable to determine whether to bin input data or not.
            except: 
                print 'Invalid number of points per bin. Type bin=number_of_points_per_bin. Now exiting'
                raise SystemExit(0)
        elif inputArg[0] == 'burst':
            UseOp.GRBlabel = inputArg[1]
        #Setting runOption
        elif inputArg[0] == 'option':
            if inputArg[1] == 'print-write':
                UseOp.runOption = 'LC'
                printPlot = True
                UseOp.plotComponents = False
            elif inputArg[1] == 'print-read':
                UseOp.runOption = 'LC'
                UseOp.plotComponents = True
                print_read_input = True
                printPlot = True
            elif inputArg[1] == 'print-read-ni': ### No input when print-reading
                UseOp.runOption = 'LC'
                UseOp.plotComponents = True
                print_read_input = False
                printPlot = True
            elif inputArg[1] == 'area':
                UseOp.runOption = 'LC'
                plot_area = True
            elif inputArg[1] == 'SED': #Plotting an SED
                UseOp.runOption = 'LC'
                plot_SED = True
            else:
                UseOp.runOption = inputArg[1]
 
        elif inputArg[0] == 'load':
            print "Loading parameters.txt"
            if inputArg[1] == 'output': #Loading the output file parameters.txt created in the lowest-chi2 routine
                try:
                    os.path.isfile('parameters.txt')
                    ladda = open('parameters.txt')
                    loadInText = ladda.read().split('\n')
                    ladda.close()
#                    constants = np.zeros(constantsLength)
                    for j in range(constantsLength):
                        ModVar.new_value(j , float(loadInText[j].split('=')[1]))
                    surfRingsOut = int(loadInText[constantsLength].split('=')[1])
                    loadInputConstants = False
                except: print "Could not find file parameters.txt. Loading the options.py file instead"
            elif inputArg[1] == 'chains':
                ### Reading best-fit values from chains/1-stats.dat
                print 'Loading best-fit parameters from MultiNest output'
                try:
                    open_stats = open('chains/1-stats.dat')
                    read_stats = open_stats.read().split('Dim No.')[2].split('\n')
                    open_stats.close()
                    for i_chains in range(1,n_params+1):
                        ModVar.new_value(whereParam[i_chains-1] , float(read_stats[i_chains].split()[1]))
                except:
                    print 'Failed when reading chains/1-stats.dat. Either file does not exist, or the layout of the MultiNest output has been changed. Now exiting.'
                    raise SystemExit(0)
            else: print "WARNING! No argument %s"%inputArg[1]
        elif inputArg[0] == 'band':
            bandLimit = np.zeros(2)   #bandLimit[0] lower limit, bandLimit[1] upper limit
            if (inputArg[1] == 'x-ray') or (inputArg[1] == 'X-ray') or (inputArg[1] == 'x-rays') or (inputArg[1] == 'X-rays') or (inputArg[1] == 'xrays'):
                print "Loading X-ray lightcurves >0.1 keV (2.4e16 Hz)"
                bandLimit[0] = 2.4e16
                bandLimit[1] = 0
                bandLabel = "X-rays >0.1 keV (>2.4e7 GHz)"
            elif (inputArg[1] == 'UV') or (inputArg[1] == 'uv'):
                print "Loading UV lightcurves 3.2 eV - 0.1 keV (7.7e14 - 2.4e16 Hz)"
                bandLimit[0] = 7.7e14
                bandLimit[1] = 2.4e16
                bandLabel = "UV 3.2-100 eV (7.7e5 - 2.4e7 GHz)"
            elif (inputArg[1] == 'optical-UV') or (inputArg[1] == 'UV-optical'):
                print "Loading optical and UV lightcurves 1.6eV - 3.2eV"
                bandLimit[0] = 2.e14
                bandLimit[1] = 2.4e16
                bandLabel = "UV and optical 1.6-100 eV (2e5 - 2.4e7 GHz)"
            elif (inputArg[1] == 'optical'):
                print "Loading optical lightcurves 0.4eV - 4eV"
                bandLimit[0] = 1.e14
                bandLimit[1] = 1.e15
                bandLabel = "optical 0.4-4 eV (1e5/3cm - 1e6/300nm GHz)"
            elif (inputArg[1] == 'radio-mm'):
                print "Loading radio-mm lightcurves < 1.6eV"
                bandLimit[0] = 0
                bandLimit[1] = 2.e14
                bandLabel = "radio - submm - IR < 1.6 eV (<2.4e5 GHz)"
            elif (inputArg[1] == 'radio'):
                print "Loading radio lightcurves"
                bandLimit[0] = 0
                bandLimit[1] = 6e11 #lambda = 2mm
                bandlabel = 'radio - submm (<2mm ; <600 GH)'
            elif (len(inputArg[1].split('-'))==2):
                bandLimit[0] = float(inputArg[1].split('-')[0])
                bandLimit[1] = float(inputArg[1].split('-')[1])
                printBand1, printBand1e, printBand2 , printBand2e = bandLimit[0]/1e9 , 0 , bandLimit[1]/1e9 , 0
                while printBand1 >= 10:
                    printBand1 /= 10
                    printBand1e += 1
                while printBand2 >= 10:
                    printBand2 /= 10
                    printBand2e += 1
                print "Loading %se%d - %se%d GHz"%(printBand1,printBand1e,printBand2,printBand2e)
                bandLabel = "%se%d - %se%d GHz"%(printBand1,printBand1e,printBand2,printBand2e)
            else:
                try:
                    
                    bandLimit[0] = float(inputArg[1])
                    bandLimit[1] = float(inputArg[1])
                    print "Loading %s"%bandLimit[0]
                except:
                    print "Bad band input. Loading all lightcurves..."

    print_read = printPlot and UseOp.plotComponents
    print_write = printPlot and (not UseOp.plotComponents)

        

else:
    if UseOp.allowPrint: print "No input arguments"

### Save parameters, active if save_params==True and UseOp.runOption=LC
if UseOp.save_params and (UseOp.runOption=='LC'): 
    UseOp.save_params = True
    if not os.path.isdir('Parameters'): os.system('mkdir Parameters')
else: UseOp.save_params = False


if UseOp.printStats: ### To avoid flooding run log
    if UseOp.runOption == 'fit': UseOp.printStats = False
if (UseOp.createMock) & (UseOp.runOption == 'LC'): 
    UseOp.useData = False    #Overrun
    UseOp.plotMock = True
if (UseOp.createMock == False) & (UseOp.useData == False): UseOp.plotMock = False   #Overrun
if UseOp.runOption != 'fit': UseOp.plotOutput == False
if (UseOp.createMock) & (UseOp.runOption != 'LC'): UseOp.createMock = False



#Creating pathways
inputData = '%s%s%sfit/'%(UseOp.inputData , '/'*(UseOp.inputData.split('/')[-1]!='') , UseOp.GRBlabel)
if not os.path.isdir(inputData):
    #if (UseOp.runOption == 'LC') & UseOp.createMock: 
    #    os.system('mkdir %s'%inputData)
    if UseOp.runOption == 'fit': 
        print 'No input data found. Input data should be saved in %s. Now exiting'%inputData
        raise SystemExit(0)
elif UseOp.createMock and (UseOp.runOption=='LC'):
    deleteDir = raw_input('Save files in existing directory %s? ([y]/n)'%inputData)
    if deleteDir == '': deleteDir = 'y'
    if (deleteDir != 'y') and (deleteDir != 'Y'): 
        print "Relabel the burst in file option.py or by adding flag burst=\nNow exiting"
        raise SystemExit(0)
        

color_68='#616161'
color_95='#cfcfcf'


#mockPath = '%sMock/'%inputData
#if (UseOp.runOption == 'LC') & UseOp.createMock:
#    if os.path.isdir(inputData) == False: os.system('mkdir %s'%inputData)


if os.path.isfile('interactive_log.txt') & (UseOp.runOption=='fit'): os.system('rm interactive_log.txt')
while UseOp.allowDataRemoval & UseOp.allowPrint:
    removePrevious = raw_input('OK removing previous data? [default=no, y=yes]: ') 
    if removePrevious == 'y':
        os.system('rm -r chains/') 
        useArchive = False #Whether files should be saved to archive after run or not. False: write to archive
        break
    elif removePrevious == '':
        useArchive = True
        break
    else: print "Bad input"

if os.path.isdir('chains/') == False: os.system('mkdir chains')
os.system('chmod a+rw chains/')


R = np.logspace(UseOp.gridStart,UseOp.gridEnd,UseOp.gridStep)


parameterIndeces,parIndCount = np.array([0]*np.sum(UseOp.parametrar)),0
colourCycle = ['b','g','r','c','m','y','k']     #Cycle of colours in plotting. Matplotlib standard cycle
scalePlotTime = {'d': 86400. , 'h' : 3600. , 'm' : 60. , 's' : 1. }
scaleFluxAxis = {'mJy' : 1.e3 , 'Jy' : 1. }

for parI in range(len(UseOp.parametrar)):
    if UseOp.parametrar[parI] == True: 
        parameterIndeces[parIndCount] = parI
        parIndCount += 1

if UseOp.runOption == 'plot':
    plotter()
    raise SystemExit(0)
elif UseOp.runOption == 'read-stats':
    gather_multinest_stats()
    raise SystemExit(0)
elif UseOp.runOption == 'plot-marginal':
    from plotMarginal import plot_contours
    from plotMarginal import plotMarginal
    from matplotlib import pyplot as plt
    from matplotlib.colors import LogNorm
    from useful_modules import reduce_ticks
    import scipy.ndimage as ndimage
    import kde_weights as kde
    import matplotlib as mpl

    choise = input('Plot marginal plots in a gathered sheet (1) or separately with contours (2)?: ')
    mpl.rc('xtick', labelsize=20) 
    mpl.rc('ytick', labelsize=20)    
    if choise == 1:
        plotMarginal()
    elif choise == 2: 
        ### Method:
        ### - Read in all path ways
        ### - Read in all parameters in all path ways
        ### - Create a dictionary containing, for each path way, the name of each parameter and corresponding index in the respective path way
        ### - Iterate over all combinations within the path ways N.B. Not combinations between path ways


        number_of_plots = input('Number of plots: ')
        path_to_chains = np.zeros(number_of_plots,dtype='S256')
        contour_colors = np.zeros(number_of_plots,dtype='S5')
        cmaps = np.zeros(number_of_plots,dtype=tuple)
        for set_path in range(number_of_plots):
            path_to_chains[set_path] = raw_input('Path to %d%s chains folder: '%(set_path+1,'st'*(set_path==0)+'nd'*(set_path==1)+'rd'*(set_path==2)+'th'*(set_path>2)))
            contour_colors[set_path] = raw_input('Line colour ([b],{r,g}): ')
            if contour_colors[set_path] == '': contour_colors[set_path] = 'b'
            if contour_colors[set_path] == 'b':
                cmaps[set_path] = plt.cm.Blues
            elif contour_colors[set_path] == 'r':
                cmaps[set_path] = plt.cm.Reds
            elif contour_colors[set_path] == 'g':
                cmaps[set_path] = plt.cm.Greens
            else:
                input_cmaps=raw_input('I don\'t have a preferred color map in mind. Please enter one (c.f. matplotlib colormaps): ')
                try:
                    exec('cmaps[set_path] = plt.cm.%s'%input_cmaps)
                except:
                    print 'Colormap matplotlib.pyplot.cm.%s does not exist. Now exiting.'%input_cmaps
                    raise SystemExit(0)

        if not os.path.isdir('Figures'): 
            os.system('mkdir Figures')
            print 'Created directory Figures'
        marginal_prefix = raw_input('Choose name prefix for these files: ')
        if not os.path.isdir('Figures/%s'%marginal_prefix): 
            os.system('mkdir Figures/%s'%marginal_prefix)
            print 'Created directory Figures/%s'%marginal_prefix
        
        ### Preparing to read data

        n_bins = 120
        
        number_of_params = np.zeros(number_of_plots , dtype=int)
        path_comb = np.array([])
        ### Reading path ways
        this_path = os.path.abspath('.')
        for read_path in range(number_of_plots):
            os.chdir(path_to_chains[read_path])
            abs_path = os.path.abspath('.')
            ### Reading parameters
            read_parameters = open('chains/.tmp_myCode/parameters.txt')
            these_params = read_parameters.read().split()
            read_parameters.close()

            number_of_params[read_path] = len(these_params)
                

            ### Analysing output using PyMultiNest built-in fuctions
            exec('analyze_%d = pymultinest.Analyzer(n_params = number_of_params[read_path], outputfiles_basename = \'%s/chains/1-\')'%(read_path,abs_path))

            exec('stats_in = analyze_%d.get_stats()'%read_path)
            
            ### Getting 3-sigma points from posterior
            exec('local_extremes_%d = np.zeros([number_of_params[read_path],2])'%read_path)
            for i_read in range(number_of_params[read_path]):
                which_param = np.where(these_params[i_read] == paramNamesShort)[0]
                exec('local_extremes_%d[i_read] = stats_in[\'marginals\'][i_read][\'3sigma\']'%read_path)
                
                ### Expanding a factor of 2 outside the local extremes
                if preferredScale[which_param] == 'log':
                    exec('local_extremes_%d[i_read,0] -= np.log(2)'%read_path)
                    exec('local_extremes_%d[i_read,1] += np.log(2)'%read_path)
                else:
                    exec('local_extremes_%d[i_read,0] /= 2'%read_path)
                    exec('local_extremes_%d[i_read,1] *= 2'%read_path)


            if False:

                exec('s_%d = analyze_%d.get_stats()'%(read_path,read_path))
                
                exec('p_%d = pymultinest.PlotMarginalModes(analyze_%d)'%(read_path,read_path))



            #try:
            if True:
                ### Creating dictionaries binding indeces and parameter names
                dict_text = ''
                for print_dict_text in range(len(these_params)):
                    if dict_text != '':
                        dict_text += ' , '
                    dict_text += '\'%s\' : %d'%(these_params[print_dict_text] , print_dict_text)
                exec('path_%d_dict = {%s}'%(read_path , dict_text))

                ### Saving combinations
                number_of_combos = len(these_params)*(len(these_params)-1)/2



                ### Collecting combinations
                for i_hshake in range(len(these_params)):
                    for j_hshake in range(len(these_params)):
                        if i_hshake >= j_hshake:
                            continue
                            
                        ### Checking if this combinations already is saved. If not, save it
                        if not np.sum(path_comb[np.where(path_comb == these_params[i_hshake])[0]] == these_params[j_hshake]):
                            path_comb = np.append(path_comb , [these_params[i_hshake],these_params[j_hshake]])
                            if len(path_comb) > 1:
                                path_comb = np.array(np.split(path_comb,len(path_comb)/2))
                            
                                
                print 'Path %d contains following %d parameters: \n\n%s\n'%(read_path+1 ,number_of_params[read_path] , ', '.join(these_params))
            #except:
            else:
                print 'Could not find file %s/chains/.tmp_myCode/parameters.txt. Now exiting'%path_to_chains[read_path]
                raise SystemExit(0)

        os.chdir(this_path)
        use_all_combos_input = raw_input('Data contains in total %d parameter combinations. Do you want to plot all combinations? ([y]/n): '%len(path_comb))
        use_all_combos = (use_all_combos_input == '') or (use_all_combos_input == 'Y') or (use_all_combos_input == 'y')

        ### Selecting combinations for subplots
        if not use_all_combos:
            print 'Select combinations by stating corresponding number. Separate with space for many.\n\n'
            for i_print_combos in range(len(path_comb)):
                print '(%d): %s - %s'%(i_print_combos+1 , path_comb[i_print_combos,0] , path_comb[i_print_combos,1])
            good_input,input_counter = False,0
            while not good_input:
                try:
                    newcomb_selection = np.array(raw_input('Combinations: ').split(),dtype=int)
                    good_input = True
                except:
                    input_counter += 1
                    if input_counter > 4:
                        print 'Now exiting'
                        raise SystemExit(0)
                    print 'Bad input. Try again'
                    
            path_comb = path_comb[newcomb_selection-1]

            ### Asking user to plot subplots
            if len(path_comb) == 1:
                plot_subplots_input = 'n'
            else:
                plot_subplots_input = raw_input('Plot as subplots? ([y]/n): ')
            plot_subplots = not ((plot_subplots_input == 'n') or (plot_subplots_input == 'N'))
            if plot_subplots:
                good_input,input_counter = False,0
                ### Asking user to specify subplot grid
                while not good_input:
                    if input_counter > 4:
                        print 'Now exiting'
                        raise SystemExit(0)

                    try:
                        dimensions = np.array(raw_input('Dimensions (rows x cols): ').split('x'),dtype=int)
                        if len(dimensions) != 2:
                            raise NameError()
                        good_input = True
                    except:
                        input_counter +=1
                        print 'Bad input. Type input values separated by \'x\''
                

                plot_ctrl_text = 'axes'
                
            else:
                plot_ctrl_text = 'gca()'
        else:
            plot_subplots = False
            plot_ctrl_text = 'gca()'
                
        density_temp = np.zeros([n_bins,n_bins])  ### Used when finding bin limits        
        n_params_c = np.zeros(number_of_plots)
        ### Fetching parameter and probability data
        for i_fetch in range(number_of_plots):
            exec('in_data = analyze_%d.get_data()'%i_fetch)
            

            exec('probability_%d = in_data[:,0]'%i_fetch)
            exec('parameters_%d = in_data[:,2:]'%i_fetch)
            n_params_c[i_fetch] = len(in_data[0])-2
            exec('density_%d = np.zeros([n_params_c[i_fetch],n_params_c[i_fetch],n_bins,n_bins])'%i_fetch)
            
        ### Plotting combinations
        if plot_subplots:
        ### Checking if one parameter occurs in all selected combinations. If so, plot these on x-axis
        

            if np.sum(path_comb == path_comb[0,0])==len(path_comb):
                share_x_axis = path_comb[0,0]
            elif np.sum(path_comb == path_comb[0,1])==len(path_comb):
                share_x_axis = path_comb[0,1]
            else: 
                share_x_axis = False
            if share_x_axis:
                for rect_comb in range(len(path_comb)):
                    if path_comb[rect_comb,1] == share_x_axis:
                        path_comb[rect_comb] = path_comb[rect_comb][::-1]
                    
            ### Setting subplots

            subplot_array = ''
            next_comb = 0
            sharex_owner = np.zeros(dimensions[1],dtype='S256')
            for subplot_combs_i in range(dimensions[0]):
                if dimensions[0] > 1:
                    subplot_array += '('
                for subplot_combs_j in range(dimensions[1]):
                    try:
                        subplot_array += 'fig_%s_%s%s'%(path_comb[next_comb,0],path_comb[next_comb,1], ' , '*(subplot_combs_j != (dimensions[1]-1)))
                        ### Assigning which plots to own the share pyplot.sharex object
                    
                        if (len(path_comb)-next_comb) <= dimensions[1]:
                            sharex_owner[len(path_comb)-next_comb-dimensions[1]] = 'fig_%s_%s'%(path_comb[next_comb,0],path_comb[next_comb,1])
                    except:
                        subplot_array += 'empty%s'%(' , '*(subplot_combs_j != (dimensions[1]-1)))
                    

                    next_comb += 1
                    
                if dimensions[0] > 1:
                    subplot_array += ')%s'%(','*(subplot_combs_i!=(dimensions[0]-1)))

            exec('fig,(%s) = plt.subplots(dimensions[0],dimensions[1],figsize=(6*dimensions[0],6*dimensions[1]))'%subplot_array)
            
            if share_x_axis:
                fig.subplots_adjust(wspace=0.)
        else:
            share_x_axis = False
            
                


        for plot_combs in range(len(path_comb)):
            
            if not plot_subplots:
                exec('fig_%s_%s = plt.figure(figsize=(6,6))'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))            

            for plot_path in range(number_of_plots):


                #try:
                ### Trying to fetch parameters from the path plot_patch. If it fails, the combination plot_combs isn't there

                exec('iparamx = path_%d_dict[path_comb[plot_combs,0]]'%plot_path)
                exec('iparamy = path_%d_dict[path_comb[plot_combs,1]]'%plot_path)

                ### Indeces pointing to the code arrays
                codeparam_x = np.where(path_comb[plot_combs,0]==paramNamesShort)[0]
                codeparam_y = np.where(path_comb[plot_combs,1]==paramNamesShort)[0]


                    

                ### Setting limits separating the bins



                exec('bin_limits_x = np.linspace(local_extremes_%d[iparamx,0],local_extremes_%d[iparamx,1],n_bins+1)'%(plot_path,plot_path))
                exec('bin_limits_y = np.linspace(local_extremes_%d[iparamy,0],local_extremes_%d[iparamy,1],n_bins+1)'%(plot_path,plot_path))


                ### Binning. Assumes linear bin spacing
#                for x_bins in range(n_bins):
#                    for y_bins in range(n_bins):
#                        exec('density_%d[iparamx,iparamy,y_bins,x_bins] = np.sum(probability_%d[np.where((parameters_%d[:,iparamx]>bin_limits_x[x_bins]) & (parameters_%d[:,iparamx]<=bin_limits_x[x_bins+1]) & (parameters_%d[:,iparamy]>bin_limits_y[y_bins]) & (parameters_%d[:,iparamy]<=bin_limits_y[y_bins+1]))])'%(plot_path , plot_path , plot_path , plot_path , plot_path , plot_path))



                ### Kernel density estimation, from kde_weights.py
                exec('cv = kde.Covariator(np.array([parameters_%d[:,iparamx] , parameters_%d[:,iparamy]]) , probability_%d)'%(plot_path , plot_path , plot_path))
                ic,norm = cv()
                exec('gkde = kde.gaussian_kde(np.array([parameters_%d[:,iparamx] , parameters_%d[:,iparamy]]) , probability_%d , ic , norm)'%(plot_path , plot_path , plot_path))



                plot_y , plot_x = np.meshgrid(np.linspace((bin_limits_y[0]+bin_limits_y[1])/2,(bin_limits_y[-1]+bin_limits_y[-2])/2,n_bins) , np.linspace((bin_limits_x[0]+bin_limits_x[1])/2,(bin_limits_x[-1]+bin_limits_x[-2])/2,n_bins) )

                ### Getting prob density for each point in the plot_x-plot_y grid from the gaussian kde grid
                exec('density_%d[iparamx,iparamy] = np.reshape(gkde(np.vstack([plot_x.flatten() , plot_y.flatten()])),plot_x.shape)'%plot_path)
                
                ### Locating 1- and 2- sigma levels

                exec('density_sorted = np.sort(density_%d[iparamx,iparamy].flatten())'%plot_path)
                prob_cum_sum = np.cumsum(density_sorted[::-1]*(plot_x[1,0]-plot_x[0,0])*(plot_y[0,1]-plot_y[0,0]))[::-1]
                i1sigma = density_sorted[np.argmin(np.abs(prob_cum_sum-0.68))]
                i2sigma = density_sorted[np.argmin(np.abs(prob_cum_sum-0.95))]
                if i2sigma == 0.:
                    sigma_border = 1
                    levels = '[i1sigma]'
                else:
                    sigma_border = 2
                    levels = '[i1sigma,i2sigma]'
                print i2sigma
                

                    

                ### Plotting
                        
                prior_x_limits = np.array([UseOp.paramLimits[codeparam_x,0][0],UseOp.paramLimits[codeparam_x,1][0]])
                prior_y_limits = np.array([UseOp.paramLimits[codeparam_y,0][0],UseOp.paramLimits[codeparam_y,1][0]])
                if preferredScale[np.where(paramNamesShort == path_comb[plot_combs,0])] == 'log':
                    plot_x = 10**plot_x

                    if plot_path == 0:
                        prior_x_limits = 10**prior_x_limits
                elif preferredScale[np.where(paramNamesShort == path_comb[plot_combs,0])] == 'deg':
                    plot_x *= 180 / np.pi
                if preferredScale[np.where(paramNamesShort == path_comb[plot_combs,1])] == 'log':
                    plot_y = 10**plot_y

                    if plot_path == 0:
                        prior_y_limits = 10**prior_y_limits
                elif preferredScale[np.where(paramNamesShort == path_comb[plot_combs,1])] == 'deg':
                    plot_y *= 180 / np.pi

                ### Setting limits 5% outside the 2-sigma lines
                   ### Approach - find lowest and highest value of each row to find limits of the y-axis, and of each column to find the limits of the x-axis
                ### x-axis

                
                exec('iteration_l = np.shape(density_%s[iparamx,iparamy])[0]'%plot_path)
                if plot_path == 0: ### Else we have an array from previous plotted densities in this plot/subplot
                    x_limits = np.zeros(2)
                    y_limits = np.zeros(2)

                ### x-limits
                if sigma_border == 1:
                    lim_border = i1sigma
                else:
                    lim_border = i2sigma
                    
                for i_xlim in range(iteration_l):
                    exec('this_array = density_%s[iparamx,iparamy,:,i_xlim]'%plot_path)
                    from_below = 0
                    from_above = iteration_l-1
                    while this_array[from_below] < lim_border and from_below < iteration_l-1:
                        from_below += 1
                    if plot_x[from_below,0] < x_limits[0] or i_xlim == 0:
                        x_limits[0] = plot_x[from_below,0]
                    while this_array[from_above] < lim_border and from_above > 0:
                        from_above -= 1
                    if plot_x[from_above,0] > x_limits[1] or i_xlim == 0:
                        x_limits[1] = plot_x[from_above,0]

                exec('print local_extremes_%d[iparamx]'%plot_path)

                ### y-limits
                for i_ylim in range(iteration_l):
                    exec('this_array = density_%s[iparamx,iparamy,i_ylim,:]'%plot_path)
                    from_below = 0
                    from_above = iteration_l-1
                    while this_array[from_below] < lim_border and from_below < iteration_l-1:
                        from_below += 1
                    if plot_y[0,from_below] < y_limits[0] or i_ylim == 0:
                        y_limits[0] = plot_y[0,from_below]
                    while this_array[from_above] < lim_border and from_above > 0:
                        from_above -= 1
                    if plot_y[0,from_above] > y_limits[1] or i_ylim == 0:
                        y_limits[1] = plot_y[0,from_above]
                


                ### Plotting contours and smudge


                #exec('plt.pcolormesh(plot_x,plot_y,np.ma.masked_less(density_%d[iparamx,iparamy],i2sigma) , alpha=0.5, shading=\'gouraud\', cmap=cmaps[plot_path])'%plot_path)
                
                if plot_subplots:
                    exec('cont = fig_%s_%s.contour(plot_x,plot_y,density_%d[iparamx,iparamy],levels=%s,colors=contour_colors[plot_path] , linestyles=[\'-\',\'--\'])'%(path_comb[plot_combs,0] , path_comb[plot_combs,1],plot_path,levels))
                    exec('fig_%s_%s.pcolormesh(plot_x,plot_y,density_%d[iparamx,iparamy] , alpha=0.5, shading=\'gouraud\', cmap=cmaps[plot_path])'%(path_comb[plot_combs,0],path_comb[plot_combs,1],plot_path))
                else:
                    exec('cont = plt.contour(plot_x,plot_y,density_%d[iparamx,iparamy],levels=[i1sigma,i2sigma],colors=contour_colors[plot_path] , linestyles=[\'-\',\'--\'])'%(plot_path))
                    exec('plt.pcolormesh(plot_x,plot_y,density_%d[iparamx,iparamy] , alpha=0.5, shading=\'gouraud\', cmap=cmaps[plot_path])'%plot_path)
                plt.clabel(cont,inline=True,fmt={i1sigma:r'$1\sigma$' , i2sigma:r'$2\sigma$'})
                    
                

                        
                if np.array(UseOp.preferredPlotScale)[np.where(paramNamesShort == path_comb[plot_combs,0])] == 'log':
                    exec('fig_%s_%s.%s.set_xscale("log")'%(path_comb[plot_combs,0] , path_comb[plot_combs,1],plot_ctrl_text))

                if np.array(UseOp.preferredPlotScale)[np.where(paramNamesShort == path_comb[plot_combs,1])] == 'log':
                    exec('fig_%s_%s.%s.set_yscale("log")'%(path_comb[plot_combs,0] , path_comb[plot_combs,1],plot_ctrl_text))
                print x_limits
                print y_limits
                if plot_path == (number_of_plots - 1):
                    
                    if UseOp.preferredPlotScale[iparamx] == 'log':
                        lim_distance = x_limits[1] / x_limits[0]
                        x_limits[0] /= 10**(lim_distance*5e-2)
                        x_limits[1] *= 10**(lim_distance*5e-2)
                    else:
                        lim_distance = x_limits[1]-x_limits[0]
                        x_limits[0] -= lim_distance*5e-2
                        x_limits[1] += lim_distance*5e-2


                    if UseOp.preferredPlotScale[iparamy] == 'log':
                        lim_distance = y_limits[1] / y_limits[0]
                        y_limits[0] /= 10**(lim_distance*5e-2)
                        y_limits[1] *= 10**(lim_distance*5e-2)
                    else:
                        lim_distance = y_limits[1]-y_limits[0]
                        y_limits[0] -= lim_distance*5e-2
                        y_limits[1] += lim_distance*5e-2
                    
                    print x_limits
                    print y_limits

                    exec('fig_%s_%s.%s.set_xlim(x_limits)'%(path_comb[plot_combs,0] , path_comb[plot_combs,1],plot_ctrl_text))
                    exec('fig_%s_%s.%s.set_ylim(y_limits)'%(path_comb[plot_combs,0] , path_comb[plot_combs,1],plot_ctrl_text))

                #else:
                    #exec('p_%d.plot_conditional(iparamx, iparamy, with_ellipses = False, with_points = False, grid_points=n_bins)'%plot_path)
                #except:
                
            ### Correcting plot limits to be constrained within prior limits

            """
            x_lims = np.zeros(2)
            y_lims = np.zeros(2)
            
            exec('x_lims[:] = fig_%s_%s.%s.get_xlim()'%(path_comb[plot_combs,0],path_comb[plot_combs,1],plot_ctrl_text))
            exec('y_lims[:] = fig_%s_%s.%s.get_ylim()'%(path_comb[plot_combs,0],path_comb[plot_combs,1],plot_ctrl_text))
            
            
            if x_limits[0] < x_lims[0]:
                x_limits[0] = x_lims[0]
            if x_limits[1] > x_lims[1]:
                x_limits[1] = x_lims[1]
            if y_limits[0] < y_lims[0]:
                y_limits[0] = y_lims[0]
            if y_limits[1] > y_lims[1]:
                y_limits[1] = y_lims[1]


            ### Making sure plotlims are within parameter limits stated in options.py
            
            if x_limits[0] < prior_x_limits[0]:
                x_limits[0] = prior_x_limits[0]
            if x_limits[1] > prior_x_limits[1]:
                x_limits[1] = prior_x_limits[1]
            if y_limits[0] < prior_y_limits[0]:
                y_limits[0] = prior_y_limits[0]
            if y_limits[1] > prior_y_limits[1]:
                y_limits[1] = prior_y_limits[1]

            """



            if plot_subplots:
                if not (share_x_axis and ((len(path_comb) - plot_combs)>dimensions[1])): ### If all subplots have the same x-label, only print it on the bottom two
                    print 'fig_%s_%s.axes.set_xlabel(latexParamNamesLin[np.where(paramNamesShort == path_comb[plot_combs,0])[0]])'%(path_comb[plot_combs,0],path_comb[plot_combs,1])
                    exec('fig_%s_%s.axes.set_xlabel(latexParamNamesLin[np.where(paramNamesShort == path_comb[plot_combs,0])[0]])'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                    
                    exec('fig_%s_%s.xaxis.label.set_fontsize(50)'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                    exec('fig_%s_%s.yaxis.label.set_fontsize(40)'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                    if  UseOp.preferredPlotScale[codeparam_x] != 'log':
                        exec('fig_%s_%s.ticklabel_format(style=\'sci\', axis=\'x\', scilimits=(0,2))'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                        
                    if  UseOp.preferredPlotScale[codeparam_y] != 'log':
                        exec('fig_%s_%s.ticklabel_format(style=\'sci\', axis=\'y\', scilimits=(0,2))'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                    exec('xlims = fig_%s_%s.axes.get_xlim()'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                    exec('ylims = fig_%s_%s.axes.get_ylim()'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                    exec('reduce_ticks(fig_%s_%s,include_y=True)'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                    exec('fig_%s_%s.axes.set_xlim(xlims)'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                    exec('fig_%s_%s.axes.set_ylim(ylims)'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
                elif share_x_axis:
                    print len(path_comb)
                    print plot_combs
                    print dimensions[1]
                    exec('fig_%s_%s.get_shared_x_axes().join(fig_%s_%s,%s)'%(path_comb[plot_combs,0],path_comb[plot_combs,1],path_comb[plot_combs,0],path_comb[plot_combs,1],sharex_owner[(plot_combs+dimensions[0]*dimensions[1]-len(path_comb))%dimensions[1]]))
                    print 'fig_%s_%s.get_shared_x_axes().join(fig_%s_%s,%s)'%(path_comb[plot_combs,0],path_comb[plot_combs,1],path_comb[plot_combs,0],path_comb[plot_combs,1],sharex_owner[(plot_combs+dimensions[0]*dimensions[1]-len(path_comb))%dimensions[1]])
            
                exec('fig_%s_%s.axes.set_ylabel(latexParamNamesLin[np.where(paramNamesShort == path_comb[plot_combs,1])[0]])'%(path_comb[plot_combs,0],path_comb[plot_combs,1]))
            else:
                plt.xlabel(latexParamNamesLin[np.where(paramNamesShort == path_comb[plot_combs,0])[0]],fontsize=40)
            
                plt.ylabel(latexParamNamesLin[np.where(paramNamesShort == path_comb[plot_combs,1])[0]],fontsize=40)
                ### Setting limits
#                plt.xlim(x_limits)
#                plt.ylim(y_limits)
            if (plot_subplots and (plot_combs==(len(path_comb)-1))) or not plot_subplots:
                if plot_subplots:
                    filename = 'subplots'
                else:
                    filename = '%s_%s'%(path_comb[plot_combs,0],path_comb[plot_combs,1])
                plt.tight_layout()
                exec('plt.savefig(\'Figures/%s/%s.%s\')'%(marginal_prefix,filename,UseOp.figTypes))
                print 'Saved figure Figures/%s/%s.%s'%(marginal_prefix,filename,UseOp.figTypes)
                exec('plt.savefig(\'Figures/%s/%s.%s\')'%(marginal_prefix,filename,'pdf'))
                print 'Saved figure Figures/%s/%s.%s'%(marginal_prefix,filename,'pdf')
                plt.close()

                    ### Recording maximum densities
#                    if np.max(density[i_modes,x,y]) > z_max_contour[x,y]: z_max_contour[x,y] = np.max(density[i_modes,x,y])  
                       
                        ### Setting limits
#                        if plot_path == (number_of_plots-1):
#                            if UseOp.preferredPlotScale[whereParam[x]] == 'log':
#                                plt.xlim([10**(np.floor(x_lolim[x])) , 10**(np.floor(x_hilim[x])+1)])
#                            else: plt.xlim([min(plot_x[0]) , max(plot_x[0])])
#                            if UseOp.preferredPlotScale[whereParam[y]] == 'log':
#                                plt.ylim([10**(np.floor(y_lolim[y])) , 10**(np.floor(y_hilim[y])+1)])
#                            else: plt.ylim([min(plot_y[:,0]) , max(plot_y[:,0])])

#                        plt.xlabel(latexParamNamesLin[whereParam[x]])
#                        plt.ylabel(latexParamNamesLin[whereParam[y]])
#                        cons_plot_num += 1
                        
#                        np.savetxt('Figures/%s_%s_xaxis_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , plot_x)
#                        np.savetxt('Figures/%s_%s_yaxis_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , plot_y)
#                        np.savetxt('Figures/%s_%s_density_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , density[i_modes,x,y])
                        
#                        if plot_path == (number_of_plots-1):
#                            if n_modes > 1: mode_text = 'mode_%d/'%(i_modes+1)
#                            else: mode_text = ''
                            
#                            exec('plt.savefig(\'Figures/%s/%s%s_%s.%s\')'%(marginal_prefix,mode_text,paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]],UseOp.figTypes))
#                            print 'Saved figure Figures/%s/%s%s_%s.%s'%(marginal_prefix,mode_text,paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]],UseOp.figTypes)



                    




#            probabiliy , parameters = plot_contours(os.path.abspath(path_to_chains[plot_path]))

#            probability , parameters = plot_contours(n_params,os.path.abspath(path_to_chains[plot_path]))

#            n_modes = np.shape(probability)[0]
#            density = np.zeros([n_modes,n_params,n_params,n_bins,n_bins])
        
#            for i_modes in range(n_modes):
#                if n_modes > 1: os.system('mkdir Figures/%s/mode_%d'%(marginal_prefix,i_modes+1))

#                for rescale_par in range(n_params):
#                    if preferredScale[whereParam[rescale_par]] == 'deg':
#                        parameters[i_modes,rescale_par] *= 180 / np.pi
#                sum_prob = 0.
#                prob_temp = np.copy(probability[i_modes])


                ### The arrays parameters and probability may contain zeros (if there are multiple modes), and low-probability outliers. We pick out the points within two sigma, and use them as the bin limits
#                while sum_prob < 0.98:
#                    prob_max = np.argmax(prob_temp)
#                    sum_prob += prob_temp[prob_max]
#                    prob_temp[prob_max] -= 1
#                    prob_ind = np.where(prob_temp < 0)

                ### z_max_contour contains the maximum value of each parameter pair, first dimension x, second y
#                z_max_contour = np.zeros([n_params,n_params]) 

#                for x in range(n_params):
#                    for y in range(n_params):
#                        if x >= y: continue #We do not want to plot the parameter on itself                                                                                              
#                        if plot_path == 0:
#                            plot_num[x,y] = cons_plot_num
#                        exec('fig_%d_%d = plt.figure(%d)'%(x,y,plot_num[x,y]))
#                        if preferredScale[whereParam[x]] == 'log':
#                            bin_limits_x = np.logspace(min(parameters[i_modes,x][prob_ind]),max(parameters[i_modes,x][prob_ind]),n_bins+1)

                        ### Setting limits separating the bins
#                        bin_limits_x = np.linspace(min(parameters[i_modes,x][prob_ind]),max(parameters[i_modes,x][prob_ind]),n_bins+1)
#                        bin_limits_y = np.linspace(min(parameters[i_modes,y][prob_ind]),max(parameters[i_modes,y][prob_ind]),n_bins+1)
                        
                        
                        ### Binning

#                        for x_bins in range(n_bins):
#                            for y_bins in range(n_bins):
#                                density[i_modes,x,y,x_bins,y_bins] = np.sum(probability[i_modes,np.where((parameters[i_modes,x]>bin_limits_x[x_bins]) & (parameters[i_modes,x]<bin_limits_x[x_bins+1]) & ((parameters[i_modes,y]>bin_limits_y[y_bins])) & ((parameters[i_modes,y]<bin_limits_y[y_bins+1])))])
#                        plot_x , plot_y = np.meshgrid(np.linspace(bin_limits_x[0],bin_limits_x[-1],n_bins) , np.linspace(bin_limits_y[0],bin_limits_y[-1],n_bins))

                        ### Plotting

#                        exec('fig_%d_%d = plt.figure()'%(x,y))
                        
                        ### I have done the following test. I plotted all parameters linearly without rescaling from logarithm to it's power-of-ten. The density distribution in the contour curves then looked identical to when the logarithmic scales are rescaled to it's power-of-ten and then plotted on it's logarithm.
                        
##                        if preferredScale[whereParam[x]] == 'log':
#                            plot_x = 10**plot_x
#                            plt.xscale('log')
#                        if preferredScale[whereParam[y]] == 'log':
#                            plot_y = 10**plot_y
#                            plt.yscale('log')
#                        if paramNames[whereParam[x]] == 'Gamma0 (cocoon)':
#                            plot_x = np.sqrt(plot_x**2 - 1)
#                        elif paramNames[whereParam[y]] == 'Gamma0 (cocoon)':
#                            plot_y = np.sqrt(plot_y**2 - 1)
#                        exec('fig_%d_%d = plt.contour(plot_x,plot_y, density[i_modes,x,y],n_contours,colors=contour_colors[plot_path])'%(x,y))
#                        if UseOp.preferredPlotScale[whereParam[x]] == 'log':
#                            exec('fig_%d_%d.ax.set_xscale("log")'%(x,y))
#                        if UseOp.preferredPlotScale[whereParam[y]] == 'log':
#                            exec('fig_%d_%d.ax.set_yscale("log")'%(x,y))

                        ### Correcting limits
#                        if plot_path == 0:
#                            if preferredScale[whereParam[x]] == 'log':
#                                x_lolim[x] = np.log10(min(plot_x[0]))
#                                x_hilim[x] = np.log10(max(plot_x[0]))
#                            if preferredScale[whereParam[y]] == 'log':
#                                y_lolim[y] = np.log10(min(plot_y[:,0]))#
#                                y_hilim[y] = np.log10(max(plot_y[:,0])#)


#                        else:
#                            if preferredScale[whereParam[x]] == 'log':
#                                if x_lolim[x] > np.log10(min(plot_x[0])): x_lolim[x] = np.log10(min(plot_x[0]))
#                                if x_hilim[x] < np.log10(max(plot_x[0])): x_hilim[x] = np.log10(max(plot_x[0]))
#                            if preferredScale[whereParam[y]] == 'log':
#                                if y_lolim[y] > np.log10(min(plot_y[:,0])): y_lolim[y] = np.log10(min(plot_y[:,0]))
#                                if y_hilim[y] < np.log10(max(plot_y[:,0])): y_hilim[y] = np.log10(max(plot_y[:,0]))

                        ### Recording maximum densities
#                        if np.max(density[i_modes,x,y]) > z_max_contour[x,y]: z_max_contour[x,y] = np.max(density[i_modes,x,y])  
                       
                        ### Setting limits
#                        if plot_path == (number_of_plots-1):
#                            if UseOp.preferredPlotScale[whereParam[x]] == 'log':
#                                plt.xlim([10**(np.floor(x_lolim[x])) , 10**(np.floor(x_hilim[x])+1)])
#                            else: plt.xlim([min(plot_x[0]) , max(plot_x[0])])
#                            if UseOp.preferredPlotScale[whereParam[y]] == 'log':
#                                plt.ylim([10**(np.floor(y_lolim[y])) , 10**(np.floor(y_hilim[y])+1)])
#                            else: plt.ylim([min(plot_y[:,0]) , max(plot_y[:,0])])

#                        plt.xlabel(latexParamNamesLin[whereParam[x]])
#                        plt.ylabel(latexParamNamesLin[whereParam[y]])
#                        cons_plot_num += 1
                        
#                        np.savetxt('Figures/%s_%s_xaxis_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , plot_x)
#                        np.savetxt('Figures/%s_%s_yaxis_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , plot_y)
#                        np.savetxt('Figures/%s_%s_density_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , density[i_modes,x,y])
                        
#                        if plot_path == (number_of_plots-1):
#                            if n_modes > 1: mode_text = 'mode_%d/'%(i_modes+1)
#                            else: mode_text = ''
                            
#                            exec('plt.savefig(\'Figures/%s/%s%s_%s.%s\')'%(marginal_prefix,mode_text,paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]],UseOp.figTypes))
#                            print 'Saved figure Figures/%s/%s%s_%s.%s'%(marginal_prefix,mode_text,paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]],UseOp.figTypes)
        
    else:
        print 'Invalid option %d. Now exiting.'%choise
    raise SystemExit(0)


### Load data to plot uncertainty range
def load_sigma(sigma_prob_range,sigma_prob_input,load_area_name):
    for sigma_prob in sigma_prob_range:
        load_name_lowest = '%s/lowest_%d.txt'%(load_area_name,int(100*sigma_prob))
        load_name_highest = '%s/highest_%d.txt'%(load_area_name,int(100*sigma_prob))
        load_name_tempgrid = '%s/temp_grid_%d.txt'%(load_area_name,int(100*sigma_prob))
        if sigma_prob == 0.68:
            lowest_LC_68 = np.loadtxt(load_name_lowest)
            highest_LC_68 = np.loadtxt(load_name_highest)
            sigmaTempGrid_68 = np.loadtxt(load_name_tempgrid)
        elif sigma_prob == 0.95: 
            lowest_LC_95 = np.loadtxt(load_name_lowest)
            highest_LC_95 = np.loadtxt(load_name_highest)
            sigmaTempGrid_95 = np.loadtxt(load_name_tempgrid)
    if sigma_prob_input == 3: return lowest_LC_68, highest_LC_68, sigmaTempGrid_68, lowest_LC_95, highest_LC_95, sigmaTempGrid_95 #Plot both 68% and 95 %
    elif sigma_prob_input == 1: return lowest_LC_68, highest_LC_68, sigmaTempGrid_68,[],[],[]  #Plot only 68%
    else: return [],[],[], lowest_LC_95, highest_LC_95, sigmaTempGrid_95 #Plot only 95%


#Loading data. If we are creating mock observations, data will not be loaded
if ((UseOp.runOption != 'LC') | (UseOp.createMock == False)) and (UseOp.runOption != 'plot'):
    if (UseOp.createMock == False) & ((UseOp.runOption == 'LC') & (UseOp.useData == False)):  #If we only want to plot a certain model
        if UseOp.mockDistribution == 'log': tdata = 10**np.array([np.linspace(np.log10(UseOp.mockIntervalTime[0]),np.log10(UseOp.mockIntervalTime[1]),UseOp.numberOfMock)]*len(UseOp.freqGrid)) 
        elif UseOp.mockDistribution == 'lin': tdata = np.array([np.linspace(UseOp.mockIntervalTime[0],UseOp.mockIntervalTime[1],UseOp.numberOfMock)]*len(UseOp.freqGrid)) 
        else: 
            try:tdata = np.array[UseOp.mockDistribution]
            except: 
                if UseOp.allowPrint: print "Bad mockDistribution input. Exiting."
                raise SystemExit(0)
        freqInput =  np.array(UseOp.freqGrid)
    elif plot_SED:
        try: SED_time = raw_input('Enter times, separated by space (days): ')
        except:
            print 'Bad input. Now exiting.'
            raise SystemExit(0)
        freqGridInput = SED_time.split()
        numberOfPoints = len(UseOp.freqGridInput)
        UseOp.freqGrid = np.logspace(4,20,100)
        tdata = np.zeros([100,numberOfPoints])
        print freqGridInput
        for iFreq in range(numberOfPoints):
            try:tdata[:,iFreq] = float(freqGridInput[iFreq])*86400.#np.array([[SED_time]]*len(UseOp.freqGrid))
            except: exec('tdata[:,iFreq] = %s*86400.'%freqGridInput[iFreq])
        freqInput = np.copy(UseOp.freqGrid)
#        FdataInput = 
            

    else:tdata,FdataInput,errorbarInput,errorIsArr,freqInput = loadData()
    
    if UseOp.createMock: freq = freqMock
        
            
    else: freq = freqInput
    print "Number of data points = %d"%numberOfPoints
    

#When wished output is in the temporal grid, the frequency might be a scalar, and the same for the time when the frequency is resolved. To keep the code general, we now make sure that these quantities alwas are arrays
if UseOp.mockDim == 'E': 
    if UseOp.createMock: tdata = np.array(UseOp.tgrid)
    try: iterationLength = len(tdata)    #Makes sure that input tdata is an array
    except: 
        tdata = np.array([tdata])
        iterationLength = 1
    numberOfEmpties = np.zeros(iterationLength)
    if UseOp.createMock == False:
        for i in range(iterationLength):
            numberOfEmpties[i] = np.count_nonzero(freq[i]+1)
elif UseOp.mockDim == 'T':
    if UseOp.createMock: freq = np.array(UseOp.freqGrid)
    try: iterationLength = len(freq)
    except: 
        freq = np.array([freq])  #Same as above
        iterationLength = 1
    numberOfEmpties = np.zeros(iterationLength, dtype=int)
    if UseOp.createMock == False:
        for i in range(iterationLength):
            numberOfEmpties[i] = np.count_nonzero(tdata[i]+1.)



### Plotting the 1-sigma area around the lightcurves

def plot_area_func():
    from plot_sigma import get_one_sigma
    from matplotlib import pyplot as plt

    if print_read_input:
        load_area_name = raw_input('Type GRB label to store uncertainty area files. If files exist for this name, those will be loaded. Leave blank to use the input burst name: ') #This input is used when plotting the areas
    else:
        load_area_name = ''

    ### Checking number of modes in distribution
    stats_text_open = open('chains/1-stats.dat')
    stats_text = stats_text_open.read()
    stats_text_open.close()
    
    n_modes = int(stats_text.split('Total Modes Found:')[1].split()[0])
    print n_modes
    

    ### Picking what uncertainty range to plot

    if load_area_name == '': load_area_name = UseOp.GRBlabel
    try:
        if print_read_input:
            sigma_prob_input = input('Plot 1-sigma ([1]) or 1-sigma and 2-sigma(2)?: ')
        else:
            sigma_prob_input = 1
    except: #Input must be an integer
        sigma_prob_input = 1
    if sigma_prob_input == 1: sigma_prob_range = [0.68]
    elif sigma_prob_input == 2: sigma_prob_range = [0.68,0.95]
    else: #If bad input
        print 'Bad input %s! Now exiting'%sigma_prob_input
        raise SystemExit(0)
    if os.path.isdir(load_area_name):
        lowest_LC_68, highest_LC_68, sigmaTempGrid_68, lowest_LC_95, highest_LC_95, sigmaTempGrid_95 = load_sigma(sigma_prob_range,sigma_prob_input,load_area_name)
    else:
        for sigma_prob in sigma_prob_range: ### Picking points within choosen sigma range
            for i_modes in range(1,n_modes+1):
                if i_modes == 1: all_sigma_points = get_one_sigma(n_params,sigma_prob,i_modes)
                else: all_sigma_points = np.append(all_sigma_points,get_one_sigma(n_params,sigma_prob,i_modes),0)
                
            ### Analysing total chi-square range
            chi2_max = max(all_sigma_points[:,1])
            chi2_min = min(all_sigma_points[:,1])
            print 'Total chi-square difference: %s'%(chi2_max-chi2_min)
            

        ### Picking a selection of the points at random
            ### It may happen that the fitter hasn't saved enough posterior points, e.g. if the run isn't finished. Then this is stated and plotting cancelled
            cancel_run = False
            if len(all_sigma_points) < 100:
                cancel_run_input = raw_input('Posterior only has %d points within %d sigma. Continue anyway? (y/[n]): '%(len(all_sigma_points),sigma_prob_input))
                cancel_run = (cancel_run_input == 'n') or (cancel_run_input == 'N') or (cancel_run_input == '')
            if cancel_run:
                print 'Now exiting'
                raise SystemExit(0)
            if print_read_input:
                plot_all_points = raw_input('Sample contains %d points. Plot all of them? (y/[n]): '%len(all_sigma_points))
            else:
                plot_all_points = 'y'
            plot_all_points = (plot_all_points == 'y' or plot_all_points == 'Y')
            if plot_all_points:
                div_factor = 1
            else:
                if len(all_sigma_points) < 400: 
                    div_factor = 1
                elif len(all_sigma_points) < 800: 
                    div_factor = 2
                elif len(all_sigma_points) < 4000: 
                    div_factor = 5
                else: 
                    div_factor = 10
            if div_factor == 1: 
                sigma_points = np.copy(all_sigma_points)
            else:
                sigma_points = all_sigma_points[random.sample(range(len(all_sigma_points)),len(all_sigma_points)/div_factor)]
            print 'Getting %d points from data set, using %d points for plotting uncertainty range'%(len(all_sigma_points),len(sigma_points))
            last_fraq_completion = 0.
            for sigma_LC in range(len(sigma_points)):

        ### Setting constants from the one-sigma points
                
                ModVar.new_value(whereParam , sigma_points[sigma_LC,2:])
                    
                runOption_temp = np.copy(UseOp.runOption)
                UseOp.runOption = 'one-sigma'
                lightcurve_sigma,hej,hej,sigmaTempGrid = modelFunc(R,ModVar,UseOp,PlotDetails,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,UseOp.numberOfPoints,np.sum(UseOp.parametrar),Plot_Exceptions,Plot_SED)
                UseOp.runOption = np.copy(runOption_temp)

                if sigma_LC == 0:
                    lowest_LC = np.copy(lightcurve_sigma)
                    highest_LC = np.copy(lightcurve_sigma)
                else:                

                ### Assigning values if extremes are reached
                    for extremes in range(len(lightcurve_sigma)):

                        where_lower = np.where(lightcurve_sigma[extremes] < lowest_LC[extremes])[0]
                        where_higher = np.where(lightcurve_sigma[extremes] > highest_LC[extremes])[0]
                        lowest_LC[extremes,where_lower] = lightcurve_sigma[extremes,where_lower]
                        highest_LC[extremes,where_higher] = lightcurve_sigma[extremes,where_higher]

                ### Printing completion percentage

                fraq_completion = sigma_LC / float(len(sigma_points))
                if fraq_completion > last_fraq_completion:
                    print '%d%% complete'%(int(100*last_fraq_completion))
                    last_fraq_completion += 0.1

            ### Storing data
            if not os.path.isdir(load_area_name): os.system('mkdir %s/'%load_area_name)
            lowest_filename = '%s/lowest_%d.txt'%(load_area_name,int(100*sigma_prob))
            highest_filename = '%s/highest_%d.txt'%(load_area_name,int(100*sigma_prob))
            tempgrid_filename = '%s/temp_grid_%d.txt'%(load_area_name,int(100*sigma_prob))
            np.savetxt(lowest_filename,lowest_LC)
            np.savetxt(highest_filename,highest_LC)
            np.savetxt(tempgrid_filename,sigmaTempGrid)
            print 'Files %s, %s and %s are written'%(lowest_filename,highest_filename,tempgrid_filename)
                
            if sigma_prob == 0.68:
                lowest_LC_68 = np.copy(lowest_LC)
                highest_LC_68 = np.copy(highest_LC)
                sigmaTempGrid_68 = np.copy(sigmaTempGrid)
            elif sigma_prob == 0.95:
                lowest_LC_95 = np.copy(lowest_LC)
                highest_LC_95 = np.copy(highest_LC)
                sigmaTempGrid_95 = np.copy(sigmaTempGrid)
                    
        np.savetxt('%s/freq.txt'%load_area_name,freq)

    return sigma_prob_range,sigma_prob_input,load_area_name



def lightcurve_production(freq,UseOp,numberOfEmpties,FdataInput,tdata,errorbarInput):
    from matplotlib import pyplot as plt
    from matplotlib import gridspec
    import matplotlib as mpl
    from plot_sigma import get_one_sigma
    from plotMarginal import plotMarginal
#    global FdataInput,errorbarInput,numberOfPoints,tdata
    ccgs = 2.99792458e10   #Speed of light in CGS units
    if loadInputConstants: surfRingsOut = 200   #To assure we have enough with rings in EATS integrator
    if UseOp.createMock:
        FdataInput,errorbarInput,numberOfPoints = [],[],0
    #Creating temporal grid...
        if UseOp.mockDistribution == 'log': #Logarithmic distance between mock data points
            if UseOp.mockDim == 'T':     tdata = np.array([10**np.linspace(np.log10(UseOp.mockIntervalTime[0]),np.log10(UseOp.mockIntervalTime[1]),UseOp.numberOfMock)]*len(freq))
            elif UseOp.mockDim == 'E':   freq = np.array([10**np.linspace(np.log10(UseOp.mockIntervalFreq[0]),np.log10(UseOp.mockIntervalFreq[1]),UseOp.numberOfMock)]*len(tdata))
        elif UseOp.mockDistribution == 'lin':
            if UseOp.mockDim == 'T':     tdata = np.array([np.linspace(UseOp.mockIntervalTime[0],UseOp.mockIntervalTime[1],UseOp.numberOfMock)]*len(freq))
            elif UseOp.mockDim == 'E':   freq = np.array([np.linspace(UseOp.mockIntervalFreq[0],UseOp.mockIntervalFreq[1],UseOp.numberOfMock)]*len(tdata))
        else: 
            if UseOp.mockDim == 'T':     tdata = UseOp.mockDistribution
            elif UseOp.mockDim == 'E':   freq = UseOp.mockDistribution
            if nl > 1: UseOp.numberOfMock = len(UseOp.mockDistribution[0])
            else: UseOp.numberOfMock = len(UseOp.mockDistribution)
        numberOfEmpties = np.array([UseOp.numberOfMock]*iterationLength)

        
    nl = len(freq)

    
    


    ### Plotting the lightcurves

    if not print_read:    # If this is false, the user want's to read the printed data and plot it
        if UseOp.createMock:lightcurve , tempGrid = logLikelihood([],[],[],FdataInput,tdata,errorbarInput,numberOfEmpties)
        else: lightcurve , tempGrid = logLikelihood([],[],[])

    
    ### Create and save mock observations

    if UseOp.createMock:
        if UseOp.mockDim == 'T':    mockGrid,mockLC = tdata,lightcurve   #In janskys
        elif UseOp.mockDim == 'E':  mockGrid,mockLC = freq,lightcurve
        FdataInput = lightcurve
        mockError = mockLC * UseOp.gaussianError      #For now
        #Gaussian offset
        if UseOp.gaussianOffset: 
            store_seed = np.zeros([len(mockLC),UseOp.numberOfMock])
            for i in range(len(mockLC)):
                
                for j in range(UseOp.numberOfMock):
                    #Set seed
                    time_now = time.time()
                    microseconds = int((time_now - int(time_now)) * 1e6)
                    store_seed[i,j] = microseconds
                    random.seed(microseconds)
                    #Create data points
                    if UseOp.offsetType == 'lin':
                        mockLC[i,j] = random.gauss(mockLC[i,j],mockLC[i,j]*UseOp.gaussianError)
                        while mockLC[i,j] <= 0: mockLC[i,j] = random.gauss(mockLC[i,j],mockLC[i,j]*UseOp.gaussianError) #Preventing negative flux
                    elif UseOp.offsetType == 'log':                        
                        mockLC[i,j] = np.random.lognormal(np.log(mockLC[i,j]),np.e**mockLC[i,j]*UseOp.gaussianError)
                        while mockLC[i,j] <= 0: mockLC[i,j] = random.gauss(np.log10(mockLC[i,j]),np.log10(np.abs(mockLC[i,j]))*UseOp.gaussianError) #Preventing negative flux
                    else:
                        print 'Bad offset type \'%s\'. Now exiting'
                        raise SystemError(0)
            
        print os.path.isdir(inputData)
        print inputData
        if os.path.isdir(inputData):
            loopAgain = True
            while loopAgain:
                optIn = raw_input('Erase files in %s before saving files there? ([y]/n): '%inputData)
                if optIn == '': optIn = 'y'
                if (optIn == 'y') | (optIn == 'Y'):
                    os.system('rm -r %s*'%inputData)
                    loopAgain = False
                elif (optIn == 'n') | (optIn == 'N'): loopAgain = False
        for iNu in range(iterationLength):
            outText = "%d\n"%(UseOp.numberOfMock-(len(mockLC[iNu]) - np.count_nonzero(mockLC[iNu])))
            for oi in range(np.count_nonzero(mockLC[iNu]+1)): 
                #if mockLC[iNu,oi] != 0.:outText += "%s %s %s\n"%(mockGrid[iNu,oi],mockLC[iNu,oi],mockLC[iNu,oi]*0.1)#mockError[iNu,oi]) #Producing file with a fix errorbar (not necessarily the correct errorbar)
                if mockLC[iNu,oi] != 0.:outText += "%s %s %s\n"%(mockGrid[iNu,oi],mockLC[iNu,oi],mockError[iNu,oi])
            if UseOp.mockDim == 'T':
                if (os.path.isdir('%s'%inputData)==False): os.system('mkdir %s'%inputData)
                if os.path.isdir('%s'%inputData) == False: raw_input('Couldn\'nt create directory %s%stimeResolved! Create it manually before continuing!'%(dataFiles,'/'*(dataFiles.split('/')[-1] != '')  ) )
                os.system('touch %s%slc.%s'%(inputData,freq[iNu],UseOp.fileExtention))
                a = open('%s%slc.%s'%(inputData,freq[iNu],UseOp.fileExtention),'w')
            elif UseOp.mockDim == 'E': 
                if os.path.isdir('%s%sfrequencyResolved/'%(inputData, '/'*(inputData.split('/')[-1] != '') )) == False: os.system('mkdir %s%sfrequencyResolved'%(inputData,'/'*(inputData.split('/')[-1] != '') ))
                a = open('%s%sfrequencyResolved/%slc.%s'%(inputData,'/'*(inputData.split('/')[-1] != ''),tdata[iNu],UseOp.fileExtention),'w')
            a.write(outText)
            a.close()
        np.savetxt('%s/seed.txt'%inputData,store_seed)  #Storing seed used to produce mock data
        print "Mock data produced. Now exiting"
        raise SystemExit(0)
       

    #####################
    ### Plotting data ###
    #####################

    if UseOp.numberOfMock == 1: plotDims = raw_input('Do you want to plot the lightcurve with frequency on the x-axis? (y/n): ')
    else: plotDims = 'n'
            
     ### Writing print-outs. If user input option=print-write  -  print plot to a txt file

    if print_read:  #Loading print-outs
        if print_read_input:
            fileNameIn = raw_input('Type file name or list of file names, separated with a space [print_plot.txt]: ')
        else:
            fileNameIn = ''
        if fileNameIn == '': fileNameIn = 'print_plot.txt'
        fileName = fileNameIn.split()
    elif print_write: #Writing print-outs
        fileNameIn = raw_input('Enter name of file to save plot data in [default=print_plot.txt]: ')
        if fileNameIn == '': fileName = 'print_plot.txt'
        else: fileName = fileNameIn

        if os.path.isfile(fileName):
            print "File %s exists. Adding plot data after existing data."%fileName
            readPlotPrint = open(fileName)
            readPrint = readPlotPrint.read()
            readPlotPrint.close()
        else: 
            readPrint = ''
            os.system('touch %s'%fileName)


    ### Plotting print-outs ###

    if not print_write: #If the option is not to write a print-out file
        if print_read: #Read the print-out
            
            bandNames = ['U-band','B-band','V-band','G-band','R-band','I-band','Z-band','Y-band','J-band','H-band','K-band','L-band','M-band','N-band','Q-band']
            bandWaveLengths = [365.,445.,551.,605.,658.,806.,900.,1020.,1220.,1630.,2190.,3450.,4750.,10500.,2100.]
            if print_read_input:
                use_standard_colors_input = raw_input('Use standard color scheme? ([y]/n): ')
            else:
                use_standard_colors_input = 'y'
            if (use_standard_colors_input=='y') or (use_standard_colors_input=='Y') or (use_standard_colors_input==''): use_standard_colors = True
            else: use_standard_colors = False
            bandStandardColors = ['k','c','m','g','b','y','r','c','k','y','r','g','m','b','c','k']
            standardFreqs = [1.9e9,4.8e9,8.4e9,22.5e9,43e9,100e9,300e9,1.4e14,2.5e14,3.7e14,4.6e14,4.8e14,5.4e14,6.7e14,8.2e14,2.4e18]
            legend = ['']*iterationLength

            lineFreq = np.array([])  # Each frequency in print-out file gets an entry in this array
            lineColorMem = np.array([])
            scaleLC = np.array([],dtype=float)
            lineTypeMem = np.array([],dtype=float)
            #legend = np.array([])
            bandName = np.array([],dtype=str)
            freqMem = np.array([[]])
            tempMem = np.array([[],[]])
            fluxMem = np.array([])
            tempMemIndex = np.array([0],dtype=int)

            if print_read_input:
                plot_area_input = raw_input('Plot uncertainty areas? ([y]/n): ')
            else:
                plot_area_input = 'y'
            if (plot_area_input == 'y') or (plot_area_input == 'Y') or (plot_area_input == ''):
                sigma_prob_range,sigma_prob_input,load_area_name = plot_area_func()
                
                lowest_LC_68, highest_LC_68, sigmaTempGrid_68, lowest_LC_95, highest_LC_95, sigmaTempGrid_95 = load_sigma(sigma_prob_range,sigma_prob_input,load_area_name)
                if os.path.isfile('%s/freq.txt'%load_area_name): sigma_freq = np.loadtxt('%s/freq.txt'%load_area_name)
                else: sigma_freq = np.loadtxt('%s/freq.txt'%load_area_name)

                plot_area = True
            else:
                plot_area = False

            scaleLightcurves = raw_input("Scale the lightcurves? ([y]/n): ")
            if scaleLightcurves == '': scaleLightcurves = 'y'
            if (scaleLightcurves == 'y') or (scaleLightcurves == 'Y'):
                ### Loading data file with scale factors for lightcurve scaling
                loadScaleData = raw_input('Name of list file with scale input. Leave blank for manual input: ')
                if loadScaleData != '':
                    try:
                        openScaleData = open(loadScaleData)
                        scaleDataInput = openScaleData.read().split('\n')
                        openScaleData.close()
                    
                        ### Reading in scale data to variable scaleData (nx2 matrix)
                        scaleData = np.zeros([len(scaleDataInput)-1,2],dtype=float)
                        for iScale in range(len(scaleData)):
                            scaleData[iScale] = map(float,scaleDataInput[iScale].split())
                    except:
                        print 'No such file. Now exiting.'
                        raise SystemExit(0)
                else: #If input is left empty, code will assume manual input and giving the opportunity to save input to a text file for later use
                    saveScaleData = ''
            for iFile in range(len(fileName)):
                readPrintOpen = open(fileName[iFile])
                readPrint = readPrintOpen.read().split('\n')
                readPrintOpen.close()
                        
                readLength = len(readPrint) / 3
                readLineLength = len(readPrint[1].split())
                tempGrid,lightcurve,readFreq,linecolor = np.zeros([readLength,readLineLength]),np.zeros([readLength,readLineLength]),np.zeros(readLength),np.zeros(readLength,dtype=str)
                
                lineTypeInput = raw_input('Choose line type for file %s ([-],--,:,-.): '%fileName[iFile])
                if lineTypeInput == '': lineType = '-'
                else: lineType = np.copy(lineTypeInput)

                while not ((lineType == '-') | (lineType == '-.') | (lineType == ':') | (lineType == '--')):
                    lineType = raw_input("Bad line type \'%s\'! Try again: "%lineType)

                for iLine in range(0,len(readPrint),3):  #Looping over input frequencies
                    if readPrint[iLine] == '': break
                    readFreq[iLine/3] = float(readPrint[iLine].split(':')[1])
                    readoutFreq , readFreqExp = readFreq[iLine/3] , 0

                    while readoutFreq > 10: #Reducing frequency to a reader-friendly form
                        readoutFreq /= 10
                        readFreqExp += 1

                    tempGrid[iLine/3] = map(float,readPrint[iLine+1].split())
                    lightcurve[iLine/3] = map(float,readPrint[iLine+2].split())
                    lineTypeMem = np.append(lineTypeMem,lineType)

                    
                    tempMem = np.append(tempMem,tempGrid[iLine/3])
                    fluxMem = np.append([fluxMem],[[lightcurve[iLine/3]]])
                    freqMem = np.append([freqMem],[readFreq[iLine/3]])
                    
                    tempMemIndex = np.append(tempMemIndex,len(tempMem))

                    if not np.count_nonzero(lineFreq == readFreq[iLine/3]): #If no choise has been done for this frequency
                        lineFreq = np.append(lineFreq,readFreq[iLine/3])
                        wave_length = ccgs / lineFreq[-1] * 1e7   #Wave length in nanometers
                        if wave_length < 10: bandNameInput = 'X-Rays'
                        elif (wave_length > 30000): 
                            if readFreqExp <= 11:  #Writing out numbers smaller than 1000 without power-of-ten
                                bandNameInputTemp = readoutFreq * 10 ** (readFreqExp-9)
                                if bandNameInputTemp == int(bandNameInputTemp):
                                    bandNameInputOut = int(readoutFreq * 10 ** (readFreqExp-9))
                                else:
                                    bandNameInputOut = readoutFreq * 10 ** (readFreqExp-9)
                                bandNameInput = '%s GHz'%bandNameInputOut
                            else:
                                bandNameInput = r'$%s \times 10^{%d}$ GHz'%(readoutFreq,readFreqExp-9)
                            
                        else:
                            bandNameInput = np.copy(bandNames[np.argmin(np.abs(bandWaveLengths - wave_length))])

                        if use_standard_colors:
                            linecolor = bandStandardColors[np.argmin(np.abs(standardFreqs-lineFreq[-1]))]
                        else:
                            linecolor = raw_input('Colour for %s (matplotlib standard, html hex or gray scale): '%(bandNameInput))
                        lineColorMem = np.append(lineColorMem,linecolor)


                        ### Scaling lightcurves

                        if (scaleLightcurves == 'y') or (scaleLightcurves == 'Y'):
                            if loadScaleData == '':
                                try: scaleLCinput = input("Scale factor for %s: "%bandNameInput)
                                except: scaleLCinput = input("Bad input! Scale factor for %s: "%bandNameInput)
                                saveScaleData = '%s%s %f\n'%(saveScaleData,lineFreq[-1],scaleLCinput)
                            else:
                                ### Finding the line in scaleData variable corresponging to the desired frequency
                                scaleDataIndex = np.argmin(np.abs(scaleData[:,0] - lineFreq[-1]))
                                if scaleData[scaleDataIndex,0] != lineFreq[-1]:
                                    new_scaleLCinput = raw_input('Data with scale factors has no input corresponding to the frequency %s. Type new scale factor  or leave blank if you want to use scale factor for the closest frequency (%s Hz): '%(lineFreq[-1],scaleData[scaleDataIndex,0]))
                                    if new_scaleLCinput != '':
                                        scaleLCinput = float(new_scaleLCinput)
                                    else:                                        
                                        scaleLCinput = scaleData[scaleDataIndex,1]
                                else:                                        
                                    scaleLCinput = scaleData[scaleDataIndex,1]

                            if scaleLCinput < 1: 
                                if (int(1./scaleLCinput) == (Decimal(str(1./scaleLCinput)).quantize(Decimal('0.01')))): #Removing numeric deviations
                                    bandName = np.append(bandName,'%s / %d'%(bandNameInput,int(1./scaleLCinput)))
                                elif ((int(1./scaleLCinput)+1) == (Decimal(str(1./scaleLCinput)).quantize(Decimal('0.01')))):
                                    bandName = np.append(bandName,'%s / %d'%(bandNameInput,int(1./scaleLCinput)+1))
                                else: bandName = np.append(bandName,'%s / %s'%(bandNameInput,1./scaleLCinput))
                            else: 
                                if float(scaleLCinput) == int(scaleLCinput):bandName = np.append(bandName,'%s x %d'%(bandNameInput,int(scaleLCinput)))
                                else: np.append(bandName,'%s x %d'%(bandNameInput,scaleLCinput))
                            scaleLC = np.append(scaleLC,np.copy(scaleLCinput))
                                
                                
                        else: 
                            scaleLC = np.append(scaleLC,1.)
                            bandName = np.append(bandName,bandNameInput)
                        legend[iLine/3] = bandName[-1]
                        
                    else: 
                        lineIndex = np.where(lineFreq==readFreq[iLine/3])
                        linecolor = lineColorMem[lineIndex][0]
                        lineColorMem = np.append(lineColorMem,linecolor)
                        scaleLC = np.append(scaleLC,scaleLC[lineIndex])


            ### Saving file with scale factors
            if (scaleLightcurves == 'y') or (scaleLightcurves == 'Y'):
                if (loadScaleData == ''):
                    saveScaleDataName = raw_input('Name of file to save scale factors in. Leave blank to not save: ')
                    if saveScaleDataName != '':
                        writeScaleData = open(saveScaleDataName,'w')
                        writeScaleData.write(saveScaleData)
                        writeScaleData.close()
                
            
            ### Printing legend with corresponding colors
            for iFreqLegend in range(len(lineFreq)):
                plt.plot(0,0,lineColorMem[iFreqLegend])
            plt.legend(legend,loc=input('Position of legend (1-4): '))


            for iPlot in range(len(freqMem)):
                plt.plot(tempMem[tempMemIndex[iPlot]:tempMemIndex[iPlot+1]]/scalePlotTime[UseOp.daysOrSec],fluxMem[tempMemIndex[iPlot]:tempMemIndex[iPlot+1]] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot],lineTypeMem[iPlot],color=lineColorMem[iPlot])

                ### Plotting uncertainty areas
                if plot_area:
                    try:
                        len(sigma_freq)
                        sigma_freq_is_scalar = False
                    except:
                        sigma_freq_is_scalar = True
                    if sigma_freq_is_scalar:
                        #which_sigma = np.where(sigma_freq == freqMem[iPlot])[0]
                        if (sigma_prob_input == 1) or (sigma_prob_input == 2):
                            plt.fill_between(sigmaTempGrid_68/scalePlotTime[UseOp.daysOrSec] , highest_LC_68 * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot] , lowest_LC_68 * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot],color=color_68)
                        if (sigma_prob_input == 2):
                            plt.fill_between(sigmaTempGrid_95/scalePlotTime[UseOp.daysOrSec] , highest_LC_95 * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot] , highest_LC_68 * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot],color=color_95)
                            plt.fill_between(sigmaTempGrid_95/scalePlotTime[UseOp.daysOrSec] , lowest_LC_95 * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot] , lowest_LC_68 * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot],color=color_95)
                    else:
                        which_sigma = np.where(sigma_freq == freqMem[iPlot])[0]
                        if (sigma_prob_input == 1) or (sigma_prob_input == 2):
                            plt.fill_between(sigmaTempGrid_68[which_sigma][0]/scalePlotTime[UseOp.daysOrSec] , highest_LC_68[which_sigma][0] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot] , lowest_LC_68[which_sigma][0] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot],color=color_68)
                        if (sigma_prob_input == 2):
                            plt.fill_between(sigmaTempGrid_95[which_sigma][0]/scalePlotTime[UseOp.daysOrSec] , highest_LC_95[which_sigma][0] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot] , highest_LC_68[which_sigma][0] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot],color=color_95)
                            plt.fill_between(sigmaTempGrid_95[which_sigma][0]/scalePlotTime[UseOp.daysOrSec] , lowest_LC_95[which_sigma][0] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot] , lowest_LC_68[which_sigma][0] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[iPlot],color=color_95)
            plot_area = False
        
    if printPlot:# ### Ordering input data. Useful when writing print-out and plotting data from print-out ###
#        printFreqTemp = np.copy(freq)
#        printFreqOrd = np.zeros(iterationLength,dtype=int)
        
#        for orderFreq in range(iterationLength):
#            printFreqOrd[orderFreq] = np.argmin(printFreqTemp)
#            printFreqTemp[printFreqOrd[orderFreq]] += 1e20
        printFreqOrd = np.argsort(freq)

        printFreqOut = freq[printFreqOrd]
        tempGridOut = tempGrid[printFreqOrd]
        lightcurveOut = lightcurve[printFreqOrd]
        
        FdataPrint = np.copy(FdataInput[printFreqOrd])
        tdataPrint = np.copy(tdata[printFreqOrd])
        errorbarPrint = errorbarInput[printFreqOrd]

            
    for i in range(iterationLength):  #Looping over all frequencies
        ### Index i gets the items from the load data, and from the lightcurve output. Index printFreq[i] finds corresponding index for the saved deviation and scaling data

        ### Writing print-out

        if print_write: 
            printFreq, printFreqIte = printFreqOut[i], 0
            while printFreq >= 10:
                printFreq /= 10
                printFreqIte += 1

            
            readPrint += 'Band:%se%d\n%s\n'%(printFreq,printFreqIte,'\n'.join([' '.join(map(str,tempGridOut[i])),' '.join(map(str,lightcurveOut[i]))]))
            

        elif not print_write:   #If not read print-out
            if (plotDims == 'y') | (plotDims == 'Y'):
                plt.plot(UseOp.freqGrid,FdataInput[:,0] * scaleFluxAxis[UseOp.fluxAxis],'%so'%colourCycle[i%len(colourCycle)])
            
            if (not print_read) and  (not UseOp.createMock) and (not plot_SED):  ### Plotting LC
                plt.plot(tempGrid[i]/scalePlotTime[UseOp.daysOrSec],lightcurve[i] * scaleFluxAxis[UseOp.fluxAxis],colourCycle[i%len(colourCycle)])

            if not plot_SED:
                if print_read: n_o_e_iterator = np.copy(numberOfEmpties[printFreqOrd[i]])
                else: n_o_e_iterator = np.copy(numberOfEmpties[i])
                for plotLC in range(n_o_e_iterator):
                    #Finding limits for plot
                
                    minFdataThis , maxFdataThis = FdataInput[i,np.argmin(FdataInput[i,:numberOfEmpties[i]])] , FdataInput[i,np.argmax(FdataInput[i,:numberOfEmpties[i]])]
                    try: 
                        if minFdataThis < minFdata: minFdata = minFdataThis
                    except: minFdata = minFdataThis
                    try: 
                        if maxFdataThis > maxFdata: maxFdata = maxFdataThis
                    except: maxFdata = maxFdataThis
                    minTdataThis , maxTdataThis = tdata[i,np.argmin(tdata[i,:numberOfEmpties[i]])] , tdata[i,np.argmax(tdata[i,:numberOfEmpties[i]])]
                    try: 
                        if minTdataThis < minTdata: minTdata = minTdataThis
                    except: minTdata = minTdataThis
                    try: 
                        if maxTdataThis > maxTdata: maxTdata = maxTdataThis
                    except: maxTdata = maxTdataThis

                    if print_read: is_error = errorIsArr[printFreqOrd[i],plotLC]
                    else: is_error = errorIsArr[i,plotLC]
                    if is_error: #Plotting dot and errorbar. Else: plotting upper limit

                        if print_read:
#                            yerror = FdataPrint[i,plotLC]
#                            dz = yerror / FdataInput[printFreqOrd[i],plotLC] / np.log(10)  ### Transform from linear to logarithmic standard deviation. No need to include scaling, since the scale enters both in numerator and denominator
#                            yerr_upper = 10**(np.log10(FdataInput[printFreqOrd[i],plotLC]) + dz)
#                            yerr_lower = 10**(np.log10(FdataInput[printFreqOrd[i],plotLC]) - dz)

                            plt.plot(tdata[printFreqOrd[i],plotLC]/scalePlotTime[UseOp.daysOrSec],FdataInput[printFreqOrd[i],plotLC] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[i],'%so'%lineColorMem[i])
                            plt.errorbar(tdata[printFreqOrd[i],plotLC]/scalePlotTime[UseOp.daysOrSec],FdataInput[printFreqOrd[i],plotLC] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[i],yerr=errorbarInput[printFreqOrd[i],plotLC]*scaleLC[i]*scaleFluxAxis[UseOp.fluxAxis],ecolor='%s'%lineColorMem[i])
#                        plt.plot(tdataPrint[i,plotLC]/scalePlotTime[daysOrSec],FdataPrint[i,plotLC] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[i],'%so'%lineColorMem[i])
#                        plt.errorbar(tdataPrint[i,plotLC]/scalePlotTime[daysOrSec],FdataPrint[i,plotLC] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[i],yerr=errorbarPrint[i,plotLC]*scaleLC[i],ecolor='%s'%lineColorMem[i])#,elinewidth=1)
                        else:
                            yerror = errorbarInput[i,plotLC]
                            dz = yerror / FdataInput[i,plotLC] / np.log(10)  ### Transform from linear to logarithmic standard deviation. No need to include scaling, since the scale enters both in numerator and denominator
                            yerr_upper = 10**(np.log10(FdataInput[i,plotLC]) + dz)
                            yerr_lower = 10**(np.log10(FdataInput[i,plotLC]) - dz)
                        
                            plt.plot(tdata[i,plotLC]/scalePlotTime[UseOp.daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[UseOp.fluxAxis],'%so'%colourCycle[i%len(colourCycle)])
                            plt.errorbar(tdata[i,plotLC]/scalePlotTime[UseOp.daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[UseOp.fluxAxis],yerr=errorbarInput[i,plotLC]*scaleFluxAxis[UseOp.fluxAxis],ecolor='%s'%colourCycle[i%len(colourCycle)])#,elinewidth=1)
                    else:
                        if print_read:
                            plt.errorbar(tdata[printFreqOrd[i],plotLC]/scalePlotTime[UseOp.daysOrSec],FdataInput[printFreqOrd[i],plotLC] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[i]*0.75,yerr=0.25*FdataInput[printFreqOrd[i],plotLC] * scaleFluxAxis[UseOp.fluxAxis] * scaleLC[i] , lolims=True,ecolor='%s'%lineColorMem[i])  ### Plotting upper limits in the print-read option. Factors 0.75 and 0.25 are used to make the arrow stem from the FdataInput[printFreqOrd[i],plotLC] value
                        else:
                            plt.errorbar(tdata[i,plotLC]/scalePlotTime[UseOp.daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[UseOp.fluxAxis]*0.75,yerr=FdataInput[i,plotLC] * scaleFluxAxis[UseOp.fluxAxis]*0.25,lolims=True,ecolor='%s'%colourCycle[i%len(colourCycle)]) 

                        
                #else:plt.errorbar(tdata[i,plotLC]/scalePlotTime[daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[UseOp.fluxAxis]*0.75,yerr=FdataInput[i,plotLC] * scaleFluxAxis[UseOp.fluxAxis]*0.25,lolims=True,ecolor='%s'%colourCycle[i%len(colourCycle)]) 

    if print_write:
        writePlotPrint = open(fileName,'w')
        writePlotPrint.write(readPrint)
        writePlotPrint.close()
    elif not plot_SED:
        plt.xlim([0.5*minTdata/scalePlotTime[UseOp.daysOrSec],2*maxTdata/scalePlotTime[UseOp.daysOrSec]])
        plt.ylim([0.5*minFdata * scaleFluxAxis[UseOp.fluxAxis],2*maxFdata * scaleFluxAxis[UseOp.fluxAxis]])
                
    #    elif mockDim == 'E': #option 'E' is obsolete
    #        plt.plot(freq[i,:numberOfEmpties[i]]/scalePlotTime,lightcurve[i,:numberOfEmpties[i]] * scaleFluxAxis,colourCycle[i%len(colourCycle)])
    #        if UseOp.useData: plt.plot(freq[i,:numberOfEmpties[i]]/scalePlotTime,FdataInput[i,:numberOfEmpties[i]] * scaleFluxAxis,'%so'%colourCycle[i%len(colourCycle)])
    #    plt.loglog()
   # 

        if (plotDims == 'y') | (plotDims == 'Y'):  plt.xlabel('Frequency (Hz)') 
        else: 
            plt.xlabel('Observing time (%s)'%('days'*(UseOp.daysOrSec=='d') + 'hours'*(UseOp.daysOrSec=='h') + 'minutes'*(UseOp.daysOrSec=='m') + 'seconds'*(UseOp.daysOrSec=='s')))
        plt.title(raw_input('Plot title: '))
        
#        if raw_input('Plot legend? (y/n)') == 'y':
#            legend = ['']*iterationLength
#            for iLegend in range(iterationLength):legend[iLegend] = 
#            plt.legend(legend)
        plt.ylabel('Flux (%s)'%UseOp.fluxAxis)
        plt.loglog()
        plt.show()
    else: ### Plotting SED
        print len(lightcurve)
        for iPlotSED in range(np.shape(lightcurve)[1]):#numberOfPoints):
            print np.shape(lightcurve)
            plt.plot(UseOp.freqGrid,lightcurve[:,iPlotSED],colourCycle[iPlotSED%len(colourCycle)])
        plt.loglog()
        plt.show()



    raise SystemExit(0)





#Produce surface plot of log-likelihood
if UseOp.runOption == 'surf': 
    irange,jrange = np.arange(0,1,.05),np.arange(0,1,.05)
    xGrid,yGrid = np.zeros(len(irange)),np.zeros(len(jrange))
    loglike = np.zeros([len(irange),len(jrange)])

    #Loading mock observations
    try:loadMock = np.loadtxt('mockData.txt')
    except:
        if UseOp.allowPrint: print "No previously saved mock observations (./mockData.txt) found!"
        raise SystemExit(0)
    tdata,Fdata = [],[]
    for imock in range(0,len(loadMock),2): 
        tdata.append(loadMock[imock])
        Fdata.append(loadMock[imock+1])
    errorbar = Fdata / 10. #Errorbar estimation for mock data


    for iSurfGrid in range(len(irange)):
        for jSurfGrid in range(len(jrange)):
            cube = [irange[iSurfGrid],jrange[jSurfGrid]]
            myPrior(cube,2,2)
            if iSurfGrid == 0: xGrid[jSurfGrid] = cube[1]
            if jSurfGrid == 0: yGrid[iSurfGrid] = cube[0]
            nl = len(freq)
            loglike[iSurfGrid][jSurfGrid] = logLikelihood(cube,n_params,n_params)
    xAxis,yAxis = np.meshgrid(xGrid,yGrid)
    tredfig = plt.figure()
    ax = Axes3D(tredfig)
    ax.plot_surface(xAxis,yAxis,np.log10(-loglike), rstride=1, cstride=1, cmap=plt.cm.hot)
    #ax.plot_surface(xAxis,yAxis,-loglike*2, rstride=1, cstride=1, cmap=plt.cm.hot)
    nameCount,j = 0,0
    titel = []
    for i in UseOp.parametrar:
        if i: 
            titel.append(paramNames[j])
            nameCount += 1
        j += 1
    plt.xlabel(titel[1])
    plt.ylabel(titel[0])
    #plt.zlabel('Probability')
    plt.savefig('3Dplot.jpg')
    plt.show(ax)


    
    raise SystemExit(0)


elif UseOp.runOption == 'marginal':
    plotMarginal(n_params)


############################################
#           Fitting routine                #
############################################

elif UseOp.runOption == 'fit':
    UseOp.createMock = False
    #Run MultiNest
    nl = len(freq)
    
    ### Writing a list of parameters, used for plotting marginal probabilities
    try:
        if not os.path.isdir('chains/.tmp_myCode/'): 
            try:
                os.mkdir('chains/.tmp_myCode/')
            except:
                os.mkdir('chains/')
                os.mkdir('chains/.tmp_myCode/')
        if not os.path.isfile('chains/.tmp_myCode/parameters.txt'):
            os.system('touch chains/.tmp_myCode/parameters.txt')

            write_parameter_list = open('chains/.tmp_myCode/parameters.txt','w')
            write_parameter_list.write(' '.join(np.array(paramNamesShort)[whereParam]))
            write_parameter_list.close()
    except: 
        None

    pymultinest.run(logLikelihood, myPrior, n_params, importance_nested_sampling = False, resume = True, verbose = True, sampling_efficiency = 'model', n_live_points = UseOp.livePoints,evidence_tolerance=0.5)

    
    tidPost = np.array(time.gmtime())
    tidDiff = np.array(tidPost)-np.array(tidPre)
    tidMod = [24,60,60]
    if (tidDiff[1] != 0) & UseOp.allowPrint: print "Crossing of month during simulation occured!"
    for tidCorr in range(5,2,-1): 
        if tidDiff[tidCorr] < 0: 
            tidDiff[tidCorr] += tidMod[tidCorr-3]
            tidDiff[tidCorr-1] -= 1
    if UseOp.printCrucial: print "Fitting done! \n\n\nThat was exhausting, thank heaven it's done! Time elapsed: %dh %dm %ds. \n\nAwaiting orders!"%(tidDiff[3],tidDiff[4],tidDiff[5])


### Printing best fits ###

elif UseOp.runOption=='print-stats':
     ### PRINTING STATS FROM chains/1-stats.dat ###
    path_to_chains = raw_input('Path to where to find chains directory. Use asterix analyse many paths: ')
    file_paths = glob.glob('%s/chains/1-stats.dat'%path_to_chains)
    column_title_input = np.zeros(len(file_paths),dtype='S256')
    outtext_first = '\\begin{tabular}{l'

    for i_column_titles in range(len(file_paths)):
        column_title_input[i_column_titles] = raw_input('Title for column of path %s: '%os.path.abspath(file_paths[i_column_titles]))
    chi2 = np.array([])
    chi2_sigma = np.array([])
    chi2_text = ''
    for i_file in range(len(file_paths)):
        readStats = open(file_paths[i_file],'r')
        statsIn = readStats.read().split('Dim No.')
        number_of_modes = int(statsIn[0].split('Total Modes Found:')[1].split()[0])
        print 'Posterior has %d mode%s'%(number_of_modes,'s'*(number_of_modes>1))
    #    gaussian_raw = statsIn[1].split('\n')
        number_of_parameters = np.sum(UseOp.parametrar)
        
        gaussian = np.zeros([number_of_modes,number_of_parameters,2])
        best_fit = np.zeros([number_of_modes,number_of_parameters])
        stats_diff = np.zeros([number_of_modes,number_of_parameters,2])
        #if number_of_modes > 1: outtext = '\\begin{tabular}{l|l}\nMode 1 & \\\\\nParameter & Value\\\\\n'
        #else: outtext = '\\begin{tabular}{l|l}\nParameter & %s\\\\\n'%(' & '.join(column_title))
        round_factor = 2  #How many decimals
        
        for i_modes in range(number_of_modes):
            if i_file == 0 and i_modes == 0: column_title = '%s'%column_title_input[i_file]
            else: column_title += ' & %s'%column_title_input[i_file]

            if number_of_modes > 1:  column_title += ' (Mode %d)'%(i_modes+1)

            outtext_first += '|l'
#            if i_modes > 0: outtext += '%s Mode %d'%(column_titles[i_file],i_modes+1)
            if number_of_modes > 1:
                chi2 = np.append(chi2,-2*float(statsIn[0+i_modes*3].split('Strictly Local Log-Evidence')[1].split('+/-')[0].split()[-1]))
                chi2_sigma = np.append(chi2_sigma,2*float(statsIn[0+i_modes*3].split('Strictly Local Log-Evidence')[1].split('+/-')[1].split()[0]))
            else:
                chi2 = np.append(chi2,-2*float(statsIn[0].split('+/-')[0].split()[-1]))
                chi2_sigma = np.append(chi2_sigma,2*float(statsIn[0].split('+/-')[1].split()[0]))
            chi2_text += ' & $%s\\pm%s$'%(round_off(chi2[-1] / numberOfPoints,round_factor) , round_off(chi2_sigma[-1] / numberOfPoints,round_factor))
            for i_stats in range(number_of_parameters):
                gaussian[i_modes,i_stats] = np.array(map(float,statsIn[1+3*i_modes].split('\n')[i_stats+1].split()[1:]))
                best_fit[i_modes,i_stats] = float(statsIn[2+3*i_modes].split('\n')[i_stats+1].split()[1])
                if preferredScale[whereParam[i_stats]] == 'log':
                    stats_diff[i_modes,i_stats,0] = 10**(gaussian[i_modes,i_stats,0] + gaussian[i_modes,i_stats,1]) - 10**best_fit[i_modes,i_stats]
                    stats_diff[i_modes,i_stats,1] = 10**best_fit[i_modes,i_stats] - 10**(gaussian[i_modes,i_stats,0] - gaussian[i_modes,i_stats,1])
                    gaussian[i_modes,i_stats] = 10**gaussian[i_modes,i_stats]
                    best_fit[i_modes,i_stats] = 10**best_fit[i_modes,i_stats]
                else:
                    stats_diff[i_modes,i_stats,0] = gaussian[i_modes,i_stats,0] + gaussian[i_modes,i_stats,1] - best_fit[i_modes,i_stats]
                    stats_diff[i_modes,i_stats,1] = best_fit[i_modes,i_stats] - gaussian[i_modes,i_stats,0] + gaussian[i_modes,i_stats,1]
                if preferredScale[whereParam[i_stats]] == 'deg':
                    stats_diff[i_modes,i_stats,0] = 180 / np.pi * stats_diff[i_modes,i_stats,0]
                    stats_diff[i_modes,i_stats,1] = 180 / np.pi * stats_diff[i_modes,i_stats,1]
                    gaussian[i_modes,i_stats] = 180 / np.pi * gaussian[i_modes,i_stats]
                    best_fit[i_modes,i_stats] = 180 / np.pi * best_fit[i_modes,i_stats]

                if (stats_diff[i_modes,i_stats,0] < 0) or (stats_diff[i_modes,i_stats,1] < 0):
                    print 'best-fit value of parameter %s is outside of the standard deviation!'%paramNames[whereParam[i_stats]]
                    stats_diff[i_modes,i_stats] = np.abs(stats_diff[i_modes,i_stats])

                stats_fot = int(np.log10(best_fit[i_modes,i_stats])) #Factor of ten of best fit
                stats_fot -= (stats_fot < 0)  #Making sure we get the nearest factor of ten below the value
                if i_file == 0 and i_modes==0:
                    exec('outtext_line%d = \'%%s \'%%latexParamNamesLin[whereParam[i_stats]]'%(i_stats+i_modes*n_params))
#                if i_file == 0 and i_modes > 0:
#                    exec('outtext_line%d = \'%%s \'%%latexParamNamesLin[whereParam[i_stats]]'%(i_stats+i_modes*n_params))
                if (stats_fot >2) or (stats_fot < -2):
                    exec('outtext_line%d += \' & $%s^{+%s}_{-%s}\\\\times 10^{%d}$\''%(i_stats,round_off(best_fit[i_modes,i_stats]*10**(-stats_fot),round_factor) , round_off(stats_diff[i_modes,i_stats,0]*10**(-stats_fot),round_factor) , round_off(stats_diff[i_modes,i_stats,1]*10**(-stats_fot),round_factor) , stats_fot))
                else:
                    exec('outtext_line%d += \' & $%s^{+%s}_{-%s}$\''%(i_stats , round_off(best_fit[i_modes,i_stats],round_factor)  , round_off(stats_diff[i_modes,i_stats,0],round_factor) , round_off(stats_diff[i_modes,i_stats,1],round_factor)))
    outtext_lines = ''
    for join_lines in range(n_params):
        exec('outtext_lines += outtext_line%s'%join_lines)
        outtext_lines += '\\\\\n'
    outtext = '%s}\\\\\nParameters & %s\\\\\n%s$\\chi^2_{\\rm red}$ %s\\n\\end{tabular}'%(outtext_first,column_title,outtext_lines,chi2_text)
    output_file_name = raw_input('Output file name without suffix [output.tex]: ')
    if output_file_name == '': output_file_name = 'output'
    if os.path.isfile('%s.tex'%output_file_name):
        overwrite_input = raw_input('Overwrite file %s.tex? ([y]/n): '%output_file_name)
        if not ((overwrite_input == '') or (overwrite_input == 'y') or (overwrite_input == 'Y')):
            print 'Will not overwrite file %s.tex. Now exiting'%output_file_name
            raise SystemExit(0)
    write_outtext = open('%s.tex'%output_file_name,'w')
    write_outtext.write(outtext)
    write_outtext.close()
    

            
            


#############################################
#        Lightcurve production              #
#############################################


elif UseOp.runOption == 'LC': #Producing light curves of choosen frequencies and constants from constants.txt
    if UseOp.createMock and (UseOp.runOption=='LC'): 
        FdataInput , tdata, errorbarInput, numberOfPoints  = [],[],[],0
    elif plot_SED:
        FdataInput, errorbarInput = [],[]
    lightcurve_production(freq,UseOp,numberOfEmpties,FdataInput,tdata,errorbarInput)



 
#Prints passed time
