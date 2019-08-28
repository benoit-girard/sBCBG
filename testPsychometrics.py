#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
## testSelection.py
##
## This script evaluates the selection ability of a model with a systematic exploration of
## the evolution of the GPi output of a multiple channel model where 2 channels compete for selection.
## An adaptation of the tests initially proposed in Gurney et al., 2001b

from iniBG import *

import nest.raster_plot as raster # utile ?
from iniBG import *
# from spikeProcessing import FanoFactor, OscIndex
# from filter import lowpass
import numpy as np
from modelParams import *
restFR = {} # this will be populated with firing rates of all nuclei, at rest

#------------------------------------------
#
#
#------------------------------------------
def twoChannelPsychometricCompetition(params={}, nbInNeurons=1000, nbSteps=11, CMPfbackground=4.):

  if params['nbCh'] < 3:
    print('testSelection.py: operates on 3 channels minimum. You asked for ',params['nbCh'])
    exit()
  elif nbSteps < 2:
    print('testSelection.py: tests 3 input activity levels minimum [0, 1]. You asked for ',nbSteps)
    exit()
  else:
    print ("\nTesting the selection of model "+str(params['LG14modelID'])+", with "+str(nbSteps)+" exploration steps.")
    print ("=================================")
    print ("- numbers of activated input neurons: "+str(nbInNeurons))
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
  for n in NUCLEI:
    spkDetect[n] = [0,0,0]
    expeRate[n] = [0,0,0]

  print(expeRate,spkDetect)

  connect_detector = lambda N,C: nest.Connect(Pop[N][C], spkDetect[N][C])

  # for debug purposes:
  #CSNspkDetect = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label":N, "to_file": storeGDF, 'start':offsetDuration+simulationOffset,'stop':offsetDuration+simDuration+simulationOffset})

  #-------------------------
  # prepare the firing rates of the inputs for all the steps of the experiment
  #-------------------------
  # ranges of activation for each input, first value: 0% activation, last value: 100% activation
  FR={'CSN':[2.,20.],'PTN':[15.,46.]}
  g = {}
  g['CSN'] = FR['CSN'][1] - FR['CSN'][0]
  g['PTN'] = FR['PTN'][1] - FR['PTN'][0]
  activityLevels=[]
  for i in range(nbSteps-1):
    activityLevels.append(1./(nbSteps-1) * i)
  activityLevels.append(1.)
  #print('activityLevels[1]: ',activityLevels[1])
  rate={}
  rate['CSN'] = g['CSN'] * np.array(activityLevels) + FR['CSN'][0]*np.ones((len(activityLevels)))
  rate['PTN'] = g['PTN'] * np.array(activityLevels) + FR['PTN'][0]*np.ones((len(activityLevels)))
  #print('rate CSN',rate['CSN'])

  #-------------------------
  # and prepare the lists of neurons that will be affected by these activity changes
  #-------------------------
  ActPop= {'CSN':[(),()],'PTN':[(),()]}
  for inputName in ['CSN','PTN']:
    src = Pop[inputName]
    for i in [0,1]:
      if 'Fake' in globals():
        if inputName in Fake:
          src = Fake[inputName][i]
        if nbInNeurons==params['nb'+inputName]:
          ActPop[inputName][i] = src
        else:
          ActPop[inputName][i] = tuple(rnd.choice(a=np.array(src),size=nbInNeurons,replace=False))

  #print('actpop CSN',ActPop['CSN'][0])

  #-------------------------
  # write header in firingRate summary files (one per pop in the model)
  # made of 3 lines :
  # * the input activity level in [0,1]
  # * the correspondign firing rates from the CSNs
  # * the correspondign firing rates from the PTNs
  #-------------------------
  firingRatesFiles = {}
  frstr=''
  for i in range(nbSteps):
    frstr += str(activityLevels[i])+' , '
  frstr+='\n'
  for i in range(nbSteps):
    frstr += str(rate['CSN'][i])+' , '
  frstr+='\n'
  for i in range(nbSteps):
    frstr += str(rate['PTN'][i])+' , '
  frstr+='\n'

  for popName in NUCLEI:
    firingRatesFiles[popName]=[]
    for C in [0,1,2]:
      firingRatesFiles[popName].append(open(dataPath+'firingRates_'+popName+'_Ch'+str(C)+'_Selection.csv','a'))
      firingRatesFiles[popName][C].writelines(frstr)

  for C in range(params['nbCh']):
    if 'Fake' in globals():
      if 'CMPf' in Fake:
        nest.SetStatus(Fake['CMPf'][C],{'rate':CMPfbackground})
    else:
      nest.SetStatus(Pop['CMPf'][C],{'rate':CMPfbackground})

  #-------------------------
  # Simulation
  #-------------------------
  for i1 in range(nbSteps):
      print('STEP ',i1*nbSteps,'/',nbSteps)
      print('* channel 1 act. level '+str(activityLevels[i1])+' (CSN: '+str(rate['CSN'][i1])+' Hz - PTN: '+str(rate['PTN'][i1])+' Hz)')
      print('* channel 2 act. level '+str(activityLevels[nbSteps-i1-1])+' (CSN: '+str(rate['CSN'][nbSteps-i1-1])+' Hz - PTN: '+str(rate['PTN'][nbSteps-i1-1])+' Hz)')

      # create spike detectors for the trial, and connect them to the populations
      # necessary to do that at each trial ? Would be better to reset the existing ones...
      start=offsetDuration+simulationOffset
      stop=offsetDuration+simDuration+simulationOffset
      for N in NUCLEI:
        for C in [0,1,2]:
          spkDetect[N][C] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label":N+'-'+str(activityLevels[i1])+'-'+str(activityLevels[i2]), "to_file": storeGDF, 'start':start,'stop':stop})
          connect_detector(N,C)

      nest.SetStatus(ActPop['CSN'][0],{'rate':rate['CSN'][i1]})
      nest.SetStatus(ActPop['CSN'][1],{'rate':rate['CSN'][nbSteps-i1-1]})
      nest.SetStatus(ActPop['PTN'][0],{'rate':rate['PTN'][i1]})
      nest.SetStatus(ActPop['PTN'][1],{'rate':rate['PTN'][nbSteps-i1-1]})

      nest.Simulate(simDuration+offsetDuration)
      for N in NUCLEI:
        for C in [0,1,2]:
          expeRate[N][C] = nest.GetStatus(spkDetect[N][C], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
          frstr='%f , ' %(expeRate[N][C])
          if i1 == nbSteps-1:
            frstr+='\n'
          firingRatesFiles[N][C].write(frstr)

      simulationOffset = nest.GetKernelStatus('time')
# TODO
# 1) visualiser les gdf pour voir si tout s'est bien passé
# 2) revoir la création et l'écriture dans les firingrate.csv
  for N in NUCLEI:
    for C in [0,1,2]:
      firingRatesFiles[N][C].close()

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

  #ReactionToInput(params=params, nbInNeurons=[250,500,1000,2000,4000], activityLevels=[0., 0.2, 0.4, 0.6, 0.8, 1.])
  twoChannelPsychometricCompetition(params=params, nbInNeurons=500, nbSteps=5, CMPfbackground=4.)

  #score = np.zeros((2))
  #mapTopology2D(show=True)
  #score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=True)

#---------------------------
if __name__ == '__main__':
  main()
