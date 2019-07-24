#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
## testReactionToInput.py
##
## This script tests the evolution of the firing rates in the BG when the inputs (CSN, PTN or CM/Pf)
## gradually increase (number of recruited neurons & level of activity of these recruited neurons)

from iniBG import *

import matplotlib                 # utile ?
matplotlib.use('Agg')             # utile ?
import matplotlib.pyplot as plt   # utile ?
import pylab as pl                # utile ?
import nest.raster_plot as raster # utile ?
from iniBG import *
# from spikeProcessing import FanoFactor, OscIndex
# from filter import lowpass
import os
import numpy as np
from modelParams import *
restFR = {} # this will be populated with firing rates of all nuclei, at rest
#oscilPow = {} # Oscillations power and frequency at rest
#oscilFreq = {}

#------------------------------------------
#
# Returns dictionary of the Firing Rates observed in each nucleus, at each step
# Variables:
# - inputName: input that will be evaluated (CSN, PTN or CMPf)
# - nbInNeurons: list of nb of activated neurons to be tested (between 0 and 12000)
# - activityLevels: list of the levels of activation to be tested (between 0 and 1, in percentage)
#------------------------------------------
def ReactionToInput(showRasters=False, params={}, nbInNeurons=[4000], activityLevels=[0., 1.], logFileName=''):

  if 'inputPop' in params:
    inputName=params['inputPop']
  else:
    print('testReactionToInput.py: inputPop not specified in the params dictionary')
    exit()

  if inputName not in ['CSN','PTN','CMPf']:
    print('testReactionToInput.py: '+inputName+' not an input')
    exit()
  else:
    print ("\nTesting the reactions of model "+str(params['LG14modelID'])+" to varying inputs on the "+inputName+" input.")
    print ("=================================")
    print ("- numbers of activated neurons: "+str(nbInNeurons))
    print ("- levels of activation: "+str(activityLevels))
  #nest.ResetNetwork() # pas sur qu'on veuille faire ca
  #initNeurons() # pas sur qu'on veuille faire ca


  dataPath='log/'
  nest.SetKernelStatus({"overwrite_files":True}) # when we redo the simulation, we erase the previous traces

  nstrand.set_seed(params['nestSeed'], params['pythonSeed']) # sets the seed for the simulation

  simulationOffset = nest.GetKernelStatus('time')
  print('Simulation Offset: '+str(simulationOffset))
  if 'offsetDuration' not in params.keys():
    offsetDuration = 1000.
  else:
    offsetDuration = params['offsetDuration']

  if 'simDuration' not in params.keys():
    simDuration = 1000. # step duration period
  else:
    simDuration = params['simDuration']

  #-------------------------
  # measures
  #-------------------------
  expeRate={}
  spkDetect={} # spike detectors used to record the experiment

  # single or multi-channel?
  if params['nbCh'] != 1:
    print("ReactionToInput is to be operated with a 1 channel model only.")
    exit()

  connect_detector = lambda N: nest.Connect(Pop[N], spkDetect[N])

  # for debug purposes:
  #CSNspkDetect = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label":N, "to_file": storeGDF, 'start':offsetDuration+simulationOffset,'stop':offsetDuration+simDuration+simulationOffset})

  #-------------------------
  # prepare the firing rates of the inputs for all the steps of the experiment
  #-------------------------
  # ranges of activation for each input, first value: 0% activation, last value: 100% activation
  FR={'CSN':[2.,20.],'PTN':[15.,46.],'CMPf':[4.,34.]}

  g = FR[inputName][1] - FR[inputName][0]
  rate = g * np.array(activityLevels) + FR[inputName][0]*np.ones((len(activityLevels)))

  #-------------------------
  # and prepare the lists of neurons that will be affected by these activity changes
  #-------------------------
  src = Pop[inputName]
  if 'Fake' in globals():
    if inputName in Fake:
      src = Fake[inputName]
    if nbInNeurons==params['nb'+inputName]:
      ActPop=src
    else:
      ActPop = tuple(rnd.choice(a=np.array(src),size=max(nbInNeurons),replace=False))

  #-------------------------
  # write header in firingRate summary file
  #-------------------------
  firingRatesFile=open(dataPath+'firingRates.csv','a')
  frstr = 'nb'+inputName+ ' , act. level, FR MSN, FR FSI, FR STN, FR GPe, FR GPi\n' # not adapted to arky/proto at the moment
  firingRatesFile.writelines(frstr)

  #-------------------------
  # Simulation
  #-------------------------
  for nbN in nbInNeurons:
    for i in range(len(activityLevels)):
      print('* '+str(nbN)+' '+inputName+' neurons at '+str(rate[i])+' Hz')

      #print(ActPop)
      for N in NUCLEI:
        spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label":N+'-'+str(nbN)+'-'+str(activityLevels[i]), "to_file": storeGDF, 'start':offsetDuration+simulationOffset,'stop':offsetDuration+simDuration+simulationOffset})
        connect_detector(N)



      nest.SetStatus(ActPop,{'rate':FR[inputName][0]})
      if nbN>0:
        print(nbN,rate[i])
        nest.SetStatus(ActPop[:nbN],{'rate':rate[i]})

      nest.Simulate(simDuration+offsetDuration)

      frstr=str(nbN)+' , '+str(activityLevels[i])+' , '
      for N in NUCLEI:
        expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
        frstr += '%f , ' %(expeRate[N])


      simulationOffset = nest.GetKernelStatus('time')

      frstr+='\n'
      firingRatesFile.writelines(frstr)
  firingRatesFile.close()

  #print "************************************** file writing",text
  #res = open(dataPath+'OutSummary.txt','a')
  #res.writelines(text)
  #res.close()

#-----------------------------------------------------------------------
def main():
  if len(sys.argv) >= 2:
    print("Command Line Parameters")
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
      print("Using command line parameters")
      print(str(sys.argv))
      i = 0
      for k in paramKeys:
        i+=1
        params[k] = float(sys.argv[i])
    else :
      print("Incorrect number of parameters")

  nest.set_verbosity("M_WARNING")

  # the number of channels is expected to be 1, but not enforced here...
  instantiate_BG(params, antagInjectionSite='none', antag='')

  ReactionToInput(params=params, nbInNeurons=[250,500,1000,2000,4000], activityLevels=[0., 0.2, 0.4, 0.6, 0.8, 1.])

  #score = np.zeros((2))
  #mapTopology2D(show=True)
  #score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=True)

#---------------------------
if __name__ == '__main__':
  main()
