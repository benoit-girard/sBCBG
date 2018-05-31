#!/usr/bin/python
# -*- coding: utf-8 -*-

##
## testPlausibility.py
##
## This script tests the 14 electrophysiological constraints from LG14
## Works both in single-channel and multi-channels cases

from iniBG import *
from modelParams import *
import nest.raster_plot
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import math

restFR = {} # this will be populated with firing rates of all nuclei, at rest
oscilPow = {} # Oscillations power and frequency at rest
oscilFreq = {}

#------------------------------------------
# Checks whether the BG model respects the electrophysiological constaints (firing rate at rest).
# If testing for a given antagonist injection experiment, specifiy the injection site in antagInjectionSite, and the type of antagonists used in antag.
# Returns [score obtained, maximal score]
# params possible keys:
# - nb{MSN,FSI,STN,GPi,GPe,CSN,PTN,CMPf} : number of simulated neurons for each population
# - Ie{GPe,GPi} : constant input current to GPe and GPi
# - G{MSN,FSI,STN,GPi,GPe} : gain to be applied on LG14 input synaptic weights for each population
#------------------------------------------
def checkAvgFR(showRasters=False,params={},antagInjectionSite='none',antag='',logFileName=''):
  nest.ResetNetwork()
  initNeurons()

  showPotential = False # Switch to True to graph neurons' membrane potentials - does not handle well restarted simulations

  dataPath='log/'
  nest.SetKernelStatus({"overwrite_files":True}) # when we redo the simulation, we erase the previous traces

  nstrand.set_seed(params['nestSeed'], params['pythonSeed']) # sets the seed for the simulation

  simulationOffset = nest.GetKernelStatus('time')
  print('Simulation Offset: '+str(simulationOffset))
  offsetDuration = 1000.
  simDuration = params['tSimu'] # ms

  # single or multi-channel?
  if params['nbCh'] == 1:
    connect_detector = lambda N: nest.Connect(Pop[N], spkDetect[N])
    disconnect_detector = lambda N: nest.Disconnect(Pop[N], spkDetect[N])
    connect_multimeter = lambda N: nest.Connect(multimeters[N], [Pop[N][0]])
  else:
    connect_detector= lambda N: [nest.Connect(Pop[N][i], spkDetect[N]) for i in range(len(Pop[N]))]
    disconnect_detector= lambda N: [nest.Disconnect(Pop[N][i], spkDetect[N]) for i in range(len(Pop[N]))]
    connect_multimeter = lambda N: nest.Connect(multimeters[N], [Pop[N][0][0]])

  #-------------------------
  # measures
  #-------------------------
  spkDetect={} # spike detectors used to record the experiment
  multimeters={} # multimeters used to record one neuron in each population
  expeRate={}

  antagStr = ''
  if antagInjectionSite != 'none':
    antagStr = antagInjectionSite+'_'+antag+'_'

  for N in NUCLEI:
    # 1000ms offset period for network stabilization
    spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": antagStr+N, "to_file": storeGDF, 'start':offsetDuration+simulationOffset,'stop':offsetDuration+simDuration+simulationOffset})
    connect_detector(N)
    if showPotential:
      # multimeter records only the last 200ms in one neuron in each population
      multimeters[N] = nest.Create('multimeter', params = {"withgid": True, 'withtime': True, 'interval': 0.1, 'record_from': ['V_m'], "label": antagStr+N, "to_file": False, 'start':offsetDuration+simulationOffset+simDuration-200.,'stop':offsetDuration+simDuration+simulationOffset})
      connect_multimeter(N)

  #-------------------------
  # Simulation
  #-------------------------
  nest.Simulate(simDuration+offsetDuration)

  score = 0

  text=[]
  frstr = "#" + str(params['LG14modelID'])+ " , " + antagInjectionSite + ', '
  s = '----- RESULTS -----'
  print s
  text.append(s+'\n')
  if antagInjectionSite == 'none':
    validationStr = "\n#" + str(params['LG14modelID']) + " , "
    frstr += "none , "
    for N in NUCLEI:
      strTestPassed = 'NO!'
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
      if expeRate[N] <= FRRNormal[N][1] and expeRate[N] >= FRRNormal[N][0]:
        # if the measured rate is within acceptable values
        strTestPassed = 'OK'
        score += 1
        validationStr += N + "=OK , "
      else:
      # out of the ranges
        if expeRate[N] > FRRNormal[N][1] :
          difference = expeRate[N] - FRRNormal[N][1]
          validationStr += N + "=+%.2f , " % difference
        else:
          difference = expeRate[N] - FRRNormal[N][0]
          validationStr += N + "=%.2f , " % difference

      frstr += '%f , ' %(expeRate[N])
      s = '* '+N+' - Rate: '+str(expeRate[N])+' Hz -> '+strTestPassed+' ('+str(FRRNormal[N][0])+' , '+str(FRRNormal[N][1])+')'
      print s
      text.append(s+'\n')
      restFR[N] = str(expeRate[N])

      oscilPow[N] = -1.
      oscilFreq[N] = -1.
      try:
        spikes_N = nest.GetStatus(spkDetect[N], keys="events")[0]['times'] # get the timing of all spikes
        data = np.bincount([int(i-offsetDuration-simulationOffset) for i in spikes_N], minlength=int(simDuration)) # discretize them in bins of 1ms
        ps = np.abs(np.fft.fft(data))**2
        time_step = 1 / 1000. # 1000 ms
        freqs = np.fft.fftfreq(data.size, time_step)
        idx = np.argsort(freqs)
        posi_spectrum = np.where((freqs[idx]>0) & (freqs[idx]<200)) # restrict the analysis to freqs < 200 Hz
        oscilPow[N] = np.max(ps[idx][posi_spectrum])
        oscilFreq[N] = freqs[idx][posi_spectrum][np.where(oscilPow[N] == ps[idx][posi_spectrum])[0][0]]
        #pl.plot(freqs[idx][posi_spectrum], ps[idx][posi_spectrum]) # simple plot
        #pl.show()
      except:
        print("Power spectrum computation failed - skipping")
  else:
    validationStr = ""
    frstr += str(antag) + " , "
    for N in NUCLEI:
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
      if N == antagInjectionSite:
        strTestPassed = 'NO!'
        if expeRate[N] <= FRRAnt[N][antag][1] and expeRate[N] >= FRRAnt[N][antag][0]:
          # if the measured rate is within acceptable values
          strTestPassed = 'OK'
          score += 1
          validationStr += N + "_" + antag + "=OK , "
        else:
        # out of the ranges
          if expeRate[N] > FRRNormal[N][1] :
            difference = expeRate[N] - FRRNormal[N][1]
            validationStr += N + "_" + antag + "=+%.2f , " % difference
          else:
            difference = expeRate[N] - FRRNormal[N][0]
            validationStr += N + "_" + antag + "=%.2f , " % difference

        s = '* '+N+' with '+antag+' antagonist(s): '+str(expeRate[N])+' Hz -> '+strTestPassed+' ('+str(FRRAnt[N][antag][0])+' , '+str(FRRAnt[N][antag][1])+')'
        print s
        text.append(s+'\n')
      else:
        s = '* '+N+' - Rate: '+str(expeRate[N])+' Hz'
        print s
        text.append(s+'\n')
      frstr += '%f , ' %(expeRate[N])

  s = '-------------------'
  print s
  text.append(s+'\n')

  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','a')
  firingRatesFile.writelines(frstr)
  firingRatesFile.close()

  #print "************************************** file writing",text
  #res = open(dataPath+'OutSummary_'+logFileName+'.txt','a')
  res = open(dataPath+'OutSummary.txt','a')
  res.writelines(text)
  res.close()

  validationFile = open("validationArray.csv",'a')
  validationFile.write(validationStr)
  validationFile.close()
  #-------------------------
  # Displays
  #-------------------------
  if showRasters and interactive:
    displayStr = ' ('+antagStr[:-1]+')' if (antagInjectionSite != 'none') else ''
    for N in NUCLEI:
      # histograms crash in the multi-channels case
      nest.raster_plot.from_device(spkDetect[N], hist=(params['nbCh'] == 1), title=N+displayStr)

    if showPotential:
      pl.figure()
      nsub = 231
      for N in NUCLEI:
        pl.subplot(nsub)
        nest.voltage_trace.from_device(multimeters[N],title=N+displayStr+' #0')
        disconnect_detector(N)
        pl.axhline(y=BGparams[N]['V_th'], color='r', linestyle='-')
        nsub += 1
    pl.show()

  return score, 5 if antagInjectionSite == 'none' else 1

# -----------------------------------------------------------------------------
# This function verify if their is pauses in the GPe and if the caracteristiques
# of theses pauses are relevant with the data of the elias paper 2007 
# It is run after CheckAVGFR because it uses the gdf files of the simulation.
# -----------------------------------------------------------------------------

#---------------------------- begining getSpikes ------------------------------
# return an ordered dictionnary of the spikes occurences by neuron and in the time
def getSpikes(Directory, Nuclei):
    spikesDict = {}
    spikesList = []
    gdfList = os.listdir(Directory + '/NoeArchGdf')
    
    for f in gdfList:
        if f.find(Nuclei) != -1 and f[-4:] == ".gdf" :
            spikeData = open(Directory +'/NoeArchGdf/' + f)
            for line in spikeData: # take the spike and put it in neuronRecording
                spk = line.split('\t')
                spk.pop()
                spikesList.append(float(spk[1]))
                if spk[0] in spikesDict:
                    spikesDict[spk[0]].append(float(spk[1]))
                else:
                    spikesDict[spk[0]] = [float(spk[1])]
        
    for neuron in spikesDict:
        spikesDict[neuron] = sorted(spikesDict[neuron])
    
    return spikesDict, spikesList
#---------------------------- end getSpikes -----------------------------------
    
#--------------------------- begining getISIs ---------------------------------
# return ISIs ordered by neuron in a dictionnary
def getISIs(spikesDict):
    ISIsDict = {}
    for neuron in spikesDict:
        ISIsDict[neuron] = []
        for i in range(len(spikesDict[neuron]) - 1):
            ISIsDict[neuron].append(round(spikesDict[neuron][i+1] - spikesDict[neuron][i], 1))
    ISIsList = []
    for neuron in ISIsDict:
        for isi in ISIsDict[neuron]:
            ISIsList.append(isi)       
    return ISIsDict, ISIsList
#----------------------------- end getISIs ------------------------------------ 
    
#--------------------------- begining rasterPlot ------------------------------
# plot rasters figures in the directory /raster
def rasterPlot(spikesDict, Nuclei, Directory):
    rasterList = []
    
    if not os.path.exists(Directory + '/rasterPlot'):
        os.makedirs(Directory + '/rasterPlot')

    for neuron in spikesDict:
        rasterList.append(spikesDict[neuron])  
    plt.figure(figsize=(40,15))
    plt.eventplot(rasterList, linelengths = 0.8, linewidths = 0.6)
    plt.title('Spike raster plot ' + Nuclei)
    plt.grid()
    plt.savefig(Directory + '/rasterPlot/' + 'RasterPlot_' + Nuclei + '.png')
#----------------------------- end rasterPlot ---------------------------------
    
#--------------------------- begining BarPlot ---------------------------------
# plot the nuclei histogram of ISIs
def activityHistPlot(spikesList, Nuclei, Directory):
    
    if not os.path.exists(Directory + '/activityHistPlot'):
        os.makedirs(Directory + '/activityHistPlot')

    plt.figure(figsize=(40,5))
    plt.hist(spikesList, bins=200, normed=0.5)
    plt.title('Histogram of the activity' + Nuclei)
    plt.grid()
    plt.savefig(Directory + '/activityHistPlot/'+ 'activityHistPlot_' + Nuclei + '.png')
#----------------------------- end BarPlot ------------------------------------
    
#--------------------------- begining BarPlot ---------------------------------
# plot the nuclei histogram of ISIs
def HistPlot(ISIsList, Nuclei, Directory):
    
    if not os.path.exists(Directory + '/histPlot'):
        os.makedirs(Directory + '/histPlot')
        
    plt.figure()
    plt.hist(ISIsList, bins=20, normed=0.5)
    plt.title('Histogram ' + Nuclei)
    plt.grid()
    plt.savefig(Directory + '/histPlot/'+ 'HistPlot_' + Nuclei + '.png')
#----------------------------- end BarPlot ------------------------------------
    
#--------------------------- begining poisson ---------------------------------
# compute the poissonian probability that n or less spike occure during T ms
def poisson(n, r, T): # Tsum of 2 isi or 3 ? n = 2
    P = 0
    for i in range(n):
        P += math.pow(r*T, i)/ math.factorial(i)

    return P*math.exp(-r*T)
#----------------------------- end poisson ------------------------------------

#----------------------- begining Pause Analysis ------------------------------
def PauseAnalysis(ISIsDict,ISIsList): # Tsum of 2 isi or 3 ? n = 2
    simuSpecs = {'meanISI': np.mean(ISIsList),}
    
    r = 1/float(simuSpecs['meanISI'])
    pausesDict = {}
    pausesList = []
    coreIList = []
    
    isiThreshold = 0

    if max(ISIsList) >= 250:
        isiThreshold = 250
    elif max(ISIsList) >= 200:
        isiThreshold = 200
    elif max(ISIsList) >= 150:
        isiThreshold = 150
    elif max(ISIsList) >= 100:
        isiThreshold = 100
    elif max(ISIsList) >= 80:
        isiThreshold = 80
    elif max(ISIsList) >= 60:
        isiThreshold = 60
    elif max(ISIsList) >= 40:
        isiThreshold = 40
    else:
        isiThreshold = 20
          
    for neuron in ISIsDict:
        skip = False
        for i in range(1,len(ISIsDict[neuron])-1):
            if ISIsDict[neuron][i] >= isiThreshold and not skip :
                coreI = ISIsDict[neuron][i]
                pause = coreI
                s = -math.log10(poisson(1, r, coreI))
                s2 = -math.log10(poisson(2, r, coreI+ISIsDict[neuron][i-1]))
                s3 = -math.log10(poisson(2, r, coreI+ISIsDict[neuron][i+1]))
                if s2 > s and s2 >= s3:
                    s = s2
                    pause += ISIsDict[neuron][i-1]
                elif s3 > s:
                    s = s3
                    pause += ISIsDict[neuron][i+1]
                    skip = True
        
                if neuron in pausesDict:
                    pausesDict[neuron].append(pause)
                    pausesList.append(pause)
                    coreIList.append(coreI)
                else:
                    pausesDict[neuron] = [pause]
                    pausesList.append(pause)
                    coreIList.append(coreI)
            else:
                skip = False
        
        pausersFRRList = []
        correctedFRRList = []
        for neuron in pausesDict:
            pausersFRRList.append((len(ISIsDict[neuron])+1)*1000/float(params['tSimu']))
            pausesLength = sum(pausesDict[neuron])
            correctedFRRList.append((len(ISIsDict[neuron])-len(pausesDict[neuron])+1)*1000/(float(params['tSimu'])-pausesLength))
            
            
    
    simuSpecs['isiThreshold'] = isiThreshold
    simuSpecs['percentagePausers'] = len(pausesDict)/float(len(ISIsDict))*100
    simuSpecs['nbPausersNeurons'] = len(pausesDict)
    simuSpecs['meanPausesDuration'] = round(np.mean(pausesList),2)
    simuSpecs['meanCoreI'] = round(np.mean(coreIList),2)
    simuSpecs['nbPausesPerMin'] = round(len(pausesList)/float(len(pausesDict)*params['tSimu'])*60000,2)
    simuSpecs['nbPauses'] = len(pausesList)
    simuSpecs['meanISI'] = round(np.mean(ISIsList),2)
    simuSpecs['pausersFRR'] = round(np.mean(pausersFRRList),2)
    simuSpecs['minPausersFRR'] = round(min(pausersFRRList),2)
    simuSpecs['correctedPausersFRR'] = round(np.mean(correctedFRRList),2)

    return simuSpecs
#-------------------------- end Pause Analysis --------------------------------
    
#------------------------- begining gdf exploitation --------------------------
# call the function and plot results
def gdfExploitation(Directory):
    for N in NUCLEI:
        a = getSpikes(Directory, N)
        spikesDict = a[0]
        activityHistPlot(a[1], N, Directory)
        rasterPlot(spikesDict, N, Directory)
        
        if N == 'Arky' or N == 'Prot' or N == 'GPe':
            
            simuSpecs = PauseAnalysis(getISIs(spikesDict)[0], getISIs(spikesDict)[1])
            
            text = "\n################# Pause Results " + N + " #################"
            text += "\n ISI threshold       = " + str(simuSpecs['isiThreshold']) + " ms    | 250 ms"
            text += "\n Mean coreI duration = " + str(simuSpecs['meanCoreI']) + " ms | [200 - 600]"
            text += "\n Mean pause duration = " + str(simuSpecs['meanPausesDuration']) + " ms | 620 ms"
            text += "\n Mean ISI            = " + str(simuSpecs['meanISI']) + "  ms | 15 ms"
            text += "\n total Pauses Nb     = " + str(simuSpecs['nbPauses']) 
            text += "\n pause/min/neuron    = " + str(simuSpecs['nbPausesPerMin']) + "    | [13 - 24]"
            text += "\n Pauser neurons Nb   = " + str(simuSpecs['nbPausersNeurons'])
            text += "\n % Pauser neurons    = " + str(simuSpecs['percentagePausers'])  + "     | [60 - 100]\n"
            text += "\n pausersFRR          = " + str(simuSpecs['pausersFRR']) + "    | [37 - 54]"
            text += "\n corr pausers FRR    = " + str(simuSpecs['correctedPausersFRR']) + "     | [44 - 62]"
            text += "\n Min pausers FRR     = " + str(simuSpecs['minPausersFRR']) + "      | [30 - 54] \n"
            text += "\n#####################################################\n"
            
            res = open(Directory+'/log/OutSummary.txt','a')
            res.writelines(text)
            res.close()
            print text

#---------------------------- end gdf exploitation ----------------------------

pausesDATA = {'percentagePausers':  [40. ,    100.,   75.],        # percentage of pauser neurons in GPe [low value, High value, perfect value]
              'shortPercentageISI': [0 ,     0.70,    0.2],         # percentage of Interspike intervals inferior to 2 ms
              'meanPausesDuration': [450. ,  730.,   620.],     # change to [0.45, 0.73, 0.62] are the  extreme recorded values if it is too selective
              'nbPausesPerMin':     [8. ,     23.,    13.],            # change to [8, 23, 13] are the  extreme recorded values if it is too selective
              'meanIPI':            [2.63 ,  8.74,   6.19],     # InterPauses Inteval | [2.63, 8.74, 6.19]are the  extreme recorded values if it is too selective
              'pausersFRR':         [37.48 , 71.25, 54.37],  # change to [21.47, 76.04, 54.13] which are the  extreme recorded values if it is too selective
              'correctedPausersFRR':[44.04 , 81.00, 62.52],  # change to [22.60, 86.63, 62.52] which are the  extreme recorded values if it is too selective
              'nonPausersFRR':      [37.10 , 75.75, 56.43],} # change to [31.37, 91.70, 56.43] which are the  extreme recorded values if it is too selective

#------------------------------------------------------------------------------

def main():
  if len(sys.argv) >= 2:
    print "Command Line Parameters"
    paramKeys = ['LG14modelID',
                 'nbMSN',
                 'nbFSI',
                 'nbSTN',
                 'nbGPe',
                 'nbGPi',
                 'nbCSN',
                 'nbPTN',
                 'nbCMPf',
                 'GMSN',
                 'GFSI',
                 'GSTN',
                 'GGPe',
                 'GGPi', 
                 'IeGPe',
                 'IeGPi',
                 'inDegCSNMSN',
                 'inDegPTNMSN',
                 'inDegCMPfMSN',
                 'inDegFSIMSN',
                 'inDegMSNMSN', 
                 'inDegCSNFSI',
                 'inDegPTNFSI',
                 'inDegSTNFSI',
                 'inDegGPeFSI',
                 'inDegCMPfFSI',
                 'inDegFSIFSI',
                 'inDegPTNSTN',
                 'inDegCMPfSTN',
                 'inDegGPeSTN',
                 'inDegCMPfGPe',
                 'inDegSTNGPe',
                 'inDegMSNGPe',
                 'inDegGPeGPe',
                 'inDegMSNGPi',
                 'inDegSTNGPi',
                 'inDegGPeGPi',
                 'inDegCMPfGPi',
                 ]
    if len(sys.argv) == len(paramKeys)+1:
      print "Using command line parameters"
      print sys.argv
      i = 0
      for k in paramKeys:
        i+=1
        params[k] = float(sys.argv[i])
    else :
      print "Incorrect number of parameters:",len(sys.argv),"-",len(paramKeys),"expected"

  nest.set_verbosity("M_WARNING")
  
  instantiate_BG(params, antagInjectionSite='none', antag='')
  score = np.zeros((2))
  #mapTopology2D(show=True)
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=True)
  
  Directory = os.getcwd()
  os.system('mkdir NoeArchGdf')  # save the .gdf files before antagonist desaster 
  os.system('cp log/MSN* log/STN* log/Arky* log/Prot* log/GPi* log/FSI* NoeArchGdf/ ')
  os.system('rm log/MSN* log/STN* log/Arky* log/Prot* log/FSI* log/GPi*')
  gdfExploitation(Directory)

  # don't bother with deactivation tests if activities at rest are not within plausible bounds
  if score[0] < score[1]:
    print("Activities at rest do not match")
  
  #-------------------------
  print "******************"
  print "* Score:",score[0],'/',score[1]
  print "******************"

  #-------------------------
  # log the results in a file
  #-------------------------
  res = open('log/OutSummary.txt','a')
  for k,v in params.iteritems():
    res.writelines(k+' , '+str(v)+'\n')
  res.writelines("Score: "+str(score[0])+' , '+str(score[1]))
  res.close()

  res = open('score.txt','w')
  res.writelines(str(score[0])+'\n')
  res.close()

  # combined params+score output, makes it quicker to read the outcome of many experiments
  params['sim_score'] = score[0]
  params['max_score'] = score[1]
  with open('params_score.csv', 'wb') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in params.items():
       writer.writerow([key, value])
    for key, value in restFR.items():
       writer.writerow([key+'_Rate', value])
    for key, value in oscilPow.items():
       writer.writerow([key+'_Pow', value])
    for key, value in oscilFreq.items():
       writer.writerow([key+'_Freq', value])

#---------------------------
if __name__ == '__main__':
  main()
