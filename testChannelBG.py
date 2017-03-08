#!/apps/free/python/2.7.10/bin/python
# -*- coding: utf-8 -*-    
from LGneurons import *
from modelParams import *
import nest.raster_plot
#import time
import sys

# params possible keys:
# - nb{MSN,FSI,STN,GPi,GPe,CSN,PTN,CMPf} : number of simulated neurons for each population
# - Ie{GPe,GPi} : constant input current to GPe and GPi
# - G{MSN,FSI,STN,GPi,GPe} : gain to be applied on LG14 input synaptic weights for each population

#------------------------------------------
# Creates the populations of neurons necessary to simulate a BG circuit
#------------------------------------------
def createBG_MC():
  #==========================
  # Creation of neurons
  #-------------------------
  print '\nCreating neurons\n================'

  nbSim['MSN'] = params['nbMSN']
  createMC('MSN',params['nbCh'])

  nbSim['FSI'] = params['nbFSI']
  createMC('FSI',params['nbCh'])

  nbSim['STN'] = params['nbSTN']
  createMC('STN',params['nbCh'])

  nbSim['GPe'] = params['nbGPe']
  createMC('GPe',params['nbCh'])
  for i in range(len(Pop['GPe'])):
    nest.SetStatus(Pop['GPe'][i],{"I_e":params['IeGPe']})

  nbSim['GPi'] = params['nbGPi']
  createMC('GPi',params['nbCh'])
  for i in range(len(Pop['GPi'])):
    nest.SetStatus(Pop['GPi'][i],{"I_e":params['IeGPi']})

  parrot = True # switch to False at your risks & perils...                                                                                                                   
  nbSim['CSN'] = params['nbCSN']
  createMC('CSN',params['nbCh'], fake=True, parrot=parrot)

  nbSim['PTN'] = params['nbPTN']
  createMC('PTN',params['nbCh'], fake=True, parrot=parrot)

  nbSim['CMPf'] = params['nbCMPf']
  createMC('CMPf',params['nbCh'], fake=True, parrot=False)

  print "Number of simulated neurons:", nbSim

#------------------------------------------
# Connects the populations of a previously created multi-channel BG circuit 
#------------------------------------------
def connectBG_MC(antagInjectionSite,antag):
  G = {'MSN': params['GMSN'],
       'FSI': params['GFSI'],
       'STN': params['GSTN'],
       'GPe': params['GGPe'],
       'GPi': params['GGPi'],
      }

  print "Gains on LG14 syn. strength:", G

  #-------------------------
  # connection of populations
  #-------------------------
  print '\nConnecting neurons\n================'
  print "**",antag,"antagonist injection in",antagInjectionSite,"**"
  print '* MSN Inputs'
  connectMC('ex','CSN','MSN', 'focused', inDegree= min(params['inDegCSNMSN'],nbSim['CSN']),gain=G['MSN'])
  connectMC('ex','PTN','MSN', 'focused', inDegree= min(params['inDegPTNMSN'],nbSim['PTN']),gain=G['MSN'])
  connectMC('ex','CMPf','MSN','focused', inDegree= min(params['inDegCMPfMSN'],nbSim['CMPf']),gain=G['MSN'])
  connectMC('in','MSN','MSN', 'focused', inDegree= min(params['inDegMSNMSN'],nbSim['MSN']),gain=G['MSN'])
  connectMC('in','FSI','MSN', 'diffuse', inDegree= min(params['inDegFSIMSN'],nbSim['FSI']),gain=G['MSN']) # diffuse ? focused ?                                               
  # some parameterizations from LG14 have no STN->MSN or GPe->MSN synaptic contacts
  if alpha['STN->MSN'] != 0:
    print "alpha['STN->MSN']",alpha['STN->MSN']
    connectMC('ex','STN','MSN', 'diffuse', inDegree= min(params['inDegSTNMSN'],nbSim['STN']),gain=G['MSN'])
  if alpha['GPe->MSN'] != 0:
    print "alpha['GPe->MSN']",alpha['GPe->MSN']
    connectMC('in','GPe','MSN', 'diffuse', inDegree= min(params['inDegGPeMSN'],nbSim['GPe']),gain=G['MSN']) # diffuse ? focused ?                                             

  print '* FSI Inputs'
  connectMC('ex','CSN','FSI', 'focused', inDegree= min(params['inDegCSNFSI'],nbSim['CSN']),gain=G['FSI'])
  connectMC('ex','PTN','FSI', 'focused', inDegree= min(params['inDegPTNFSI'],nbSim['PTN']),gain=G['FSI'])
  if alpha['STN->FSI'] != 0:
    connectMC('ex','STN','FSI', 'diffuse', inDegree= min(params['inDegSTNFSI'],nbSim['STN']),gain=G['FSI'])
  connectMC('in','GPe','FSI', 'focused', inDegree= min(params['inDegGPeFSI'],nbSim['GPe']),gain=G['FSI'])
  connectMC('ex','CMPf','FSI','focused', inDegree= min(params['inDegCMPfFSI'],nbSim['CMPf']),gain=G['FSI'])
  connectMC('in','FSI','FSI', 'diffuse',inDegree= min(params['inDegFSIFSI'],nbSim['FSI']),gain=G['FSI'])

  print '* STN Inputs'
  connectMC('ex','PTN','STN', 'focused',inDegree= min(params['inDegPTNSTN'],nbSim['PTN']), gain=G['STN'])
  connectMC('ex','CMPf','STN','focused',inDegree= min(params['inDegCMPfSTN'],nbSim['CMPf']), gain=G['STN'])
  connectMC('in','GPe','STN', 'focused',inDegree= min(params['inDegGPeSTN'],nbSim['GPe']), gain=G['STN']) # or diffuse, to be in line with the 2008 model?                    

  print '* GPe Inputs'
  if antagInjectionSite == 'GPe':
    if   antag == 'AMPA':
      connectMC('NMDA','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connectMC('NMDA','STN','GPe','diffuse', inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
      connectMC('in','MSN','GPe','focused', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
      connectMC('in','GPe','GPe','diffuse', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe']) # diffuse or focused?                                           
    elif antag == 'NMDA':
      connectMC('AMPA','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connectMC('AMPA','STN','GPe','diffuse', inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
      connectMC('in','MSN','GPe','focused', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
      connectMC('in','GPe','GPe','diffuse', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe'])
    elif antag == 'AMPA+GABAA':
      connectMC('NMDA','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connectMC('NMDA','STN','GPe','diffuse',inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    elif antag == 'GABAA':
      connectMC('ex','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connectMC('ex','STN','GPe', 'diffuse',inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connectMC('ex','CMPf','GPe','focused',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
    connectMC('ex','STN','GPe', 'diffuse',inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    connectMC('in','MSN','GPe', 'focused',inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
    connectMC('in','GPe','GPe', 'diffuse',inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe'])

  print '* GPi Inputs'
  if antagInjectionSite =='GPi':
    if   antag == 'All':
      pass
    elif antag == 'NMDA':
      connectMC('in','MSN','GPi',   'focused',inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connectMC('AMPA','STN','GPi', 'diffuse', inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connectMC('in','GPe','GPi',   'diffuse',inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
      connectMC('AMPA','CMPf','GPi','focused',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    elif antag == 'NMDA+AMPA':
      connectMC('in','MSN','GPi', 'focused',inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connectMC('in','GPe','GPi', 'diffuse',inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
    elif antag == 'AMPA':
      connectMC('in','MSN','GPi',   'focused',inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connectMC('NMDA','STN','GPi', 'diffuse',inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connectMC('in','GPe','GPi',   'diffuse',inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
      connectMC('NMDA','CMPf','GPi','focused',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    elif antag == 'GABAA':
      connectMC('ex','STN','GPi', 'diffuse',inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connectMC('ex','CMPf','GPi','focused',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connectMC('in','MSN','GPi', 'focused',inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
    connectMC('ex','STN','GPi', 'diffuse',inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
    connectMC('in','GPe','GPi', 'diffuse',inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
    connectMC('ex','CMPf','GPi','focused',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])

#------------------------------------------
# Checks that the BG model parameterization defined by the "params" dictionary can respect the electrophysiological constaints (firing rate at rest).
# If testing for a given antagonist injection experiment, specifiy the injection site in antagInjectionSite, and the type of antagonists used in antag.
# Returns [score obtained, maximal score]
# params possible keys:
# - nb{MSN,FSI,STN,GPi,GPe,CSN,PTN,CMPf} : number of simulated neurons for each population
# - Ie{GPe,GPi} : constant input current to GPe and GPi
# - G{MSN,FSI,STN,GPi,GPe} : gain to be applied on LG14 input synaptic weights for each population
#------------------------------------------

def checkAvgFR(showRasters=False,params={},antagInjectionSite='none',antag='',logFileName=''):
  nest.ResetKernel()
  dataPath='log/'
  nest.SetKernelStatus({'local_num_threads': params['nbcpu'] if ('nbcpu' in params) else 2, "data_path": dataPath})
  initNeurons()

  offsetDuration = 1000.
  simDuration = 5000. # ms
  # nest.SetKernelStatus({"overwrite_files":True}) # Thanks to use of timestamps, file names should now 
                                                   # be different as long as they are not created during the same second

  print '/!\ Using the following LG14 parameterization',params['LG14modelID']
  loadLG14params(params['LG14modelID'])

  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','inDegCSNMSN','inDegPTNMSN','inDegCMPfMSN','inDegMSNMSN','inDegFSIMSN','inDegSTNMSN','inDegGPeMSN','inDegCSNFSI','inDegPTNFSI','inDegSTNFSI','inDegGPeFSI','inDegCMPfFSI','inDegFSIFSI','inDegPTNSTN','inDegCMPfSTN','inDegGPeSTN','inDegCMPfGPe','inDegSTNGPe','inDegMSNGPe','inDegGPeGPe','inDegMSNGPi','inDegSTNGPi','inDegGPeGPi','inDegCMPfGPi',]
  for np in necessaryParams:
    if np not in params:
      print "Missing parameter:",np 
      exit()

  #-------------------------
  # creation and connection of the neural populations
  #-------------------------
  createBG_MC()

  connectBG_MC(antagInjectionSite,antag)

  #-------------------------
  # measures
  #-------------------------
  spkDetect={} # spike detectors used to record the experiment
  expeRate={}

  antagStr = ''
  if antagInjectionSite != 'none':
    antagStr = antagInjectionSite+'_'+antag+'_'

  for N in NUCLEI:
    # 1000ms offset period for network stabilization
    spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": antagStr+N, "to_file": True, 'start':offsetDuration,'stop':offsetDuration+simDuration})
    for i in range(len(Pop[N])):
      nest.Connect(Pop[N][i], spkDetect[N])

  #-------------------------
  # Simulation
  #-------------------------
  nest.Simulate(simDuration+offsetDuration)

  score = 0

  text=[]
  frstr = antagInjectionSite + ', '
  s = '----- RESULTS -----'
  print s
  text.append(s+'\n')
  if antagInjectionSite == 'none':
    for N in NUCLEI:
      strTestPassed = 'NO!'
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
      if expeRate[N] <= FRRNormal[N][1] and expeRate[N] >= FRRNormal[N][0]:
        # if the measured rate is within acceptable values
        strTestPassed = 'OK'
        score += 1
      frstr += '%f , ' %(expeRate[N])
      s = '* '+N+' - Rate: '+str(expeRate[N])+' Hz -> '+strTestPassed+' ('+str(FRRNormal[N][0])+' , '+str(FRRNormal[N][1])+')'
      print s
      text.append(s+'\n')
  else:
    for N in NUCLEI:
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
      if N == antagInjectionSite:
        strTestPassed = 'NO!'
        if expeRate[N] <= FRRAnt[N][antag][1] and expeRate[N] >= FRRAnt[N][antag][0]:
          # if the measured rate is within acceptable values
          strTestPassed = 'OK'
          score += 1
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

  #-------------------------
  # Displays
  #-------------------------
  if showRasters and interactive:
    displayStr = ' ('+antagStr[:-1]+')' if (antagInjectionSite != 'none') else ''
    for N in NUCLEI:
      nest.raster_plot.from_device(spkDetect[N],hist=True,title=N+displayStr)

    nest.raster_plot.show()

  return score, 5 if antagInjectionSite == 'none' else 1

#-----------------------------------------------------------------------
# PActiveCNS/PTN : proportion of "active" neurons in the CSN/PTN populations (in [0.,1.])
#
#-----------------------------------------------------------------------
def checkGurneyTest(showRasters=False,params={},CSNFR=[2.,10.], PActiveCSN=1., PTNFR=[15.,35], PActivePTN=1., antagInjectionSite='none',antag=''):
  nest.ResetKernel()
  dataPath='log/'
  nest.SetKernelStatus({'local_num_threads': params['nbcpu'] if ('nbcpu' in params) else 2, "data_path": dataPath})
  initNeurons()

  offsetDuration = 200.
  simDuration = 800. # ms
  loadLG14params(params['LG14modelID'])

  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','inDegCSNMSN','inDegPTNMSN','inDegCMPfMSN','inDegMSNMSN','inDegFSIMSN','inDegSTNMSN','inDegGPeMSN','inDegCSNFSI','inDegPTNFSI','inDegSTNFSI','inDegGPeFSI','inDegCMPfFSI','inDegFSIFSI','inDegPTNSTN','inDegCMPfSTN','inDegGPeSTN','inDegCMPfGPe','inDegSTNGPe','inDegMSNGPe','inDegGPeGPe','inDegMSNGPi','inDegSTNGPi','inDegGPeGPi','inDegCMPfGPi',]
  for nepa in necessaryParams:
    if nepa not in params:
      print "Missing parameter:",nepa
      exit()

  nbRecord=2 # number of channels whose activity will be recorded
  if params['nbCh']<2:
    print 'need at least 2 channels to perform Gurney test'
    exit()
  elif params['nbCh']>2:
    nbRecord = 3 # if we have more than 2 channels, we will also record one of the neutral channels

  #-------------------------
  # creation and connection of the neural populations
  #-------------------------
  createBG_MC()
  connectBG_MC(antagInjectionSite,antag)

  #-------------------------
  # prepare the firing rates of the inputs for the 5 steps of the experiment
  #-------------------------  
  gCSN = CSNFR[1]-CSNFR[0]
  gPTN = PTNFR[1]-PTNFR[0]
  activityLevels = np.array([[0,0.4,0.4,0.6,0.4], [0.,0.,0.6,0.6,0.6]]) 

  CSNrate= gCSN * activityLevels + np.ones((5)) * CSNFR[0]
  PTNrate= gPTN * activityLevels + np.ones((5)) * PTNFR[0]

  #-------------------------
  # and prepare the lists of neurons that will be affected by these activity changes
  #-------------------------
  ActPop = {'CSN':[(),()],'PTN':[(),()]}
  if 'Fake' in globals():
    if 'CSN' in Fake:
      if PActiveCSN==1.:
       ActPop['CSN']=Fake['CSN']
      else:
        for i in range(2):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Fake['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    else:
      if PActiveCSN==1.:
       ActPop['CSN']=Pop['CSN']
      else:
        for i in range(2):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if 'PTN' in Fake:
      if PActivePTN==1.:
        ActPop['PTN']=Fake['PTN']
      else:
        for i in range(2):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Fake['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
    else:
      if PActivePTN==1.:
        ActPop['PTN']=Pop['PTN']
      else:
        for i in range(2):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
  else:
    if PActiveCSN==1.:
     ActPop['CSN']=Pop['CSN']
    else:
      for i in range(2):
        ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if PActivePTN==1.:
      ActPop['PTN']=Pop['PTN']
    else:
      for i in range(2):
        ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))

  #-------------------------
  # log-related variales
  #-------------------------
  score = 0
  expeRate={}
  for N in NUCLEI:
    expeRate[N]=-1. * np.ones((nbRecord,5))

  inspector = {}
  for N in NUCLEI:
    inspector[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": N, "to_file": False})
    for i in range(nbRecord):
      nest.Connect(Pop[N][i],inspector[N])

  #-------------------------
  # write header in firingRate summary file
  #-------------------------
  frstr = 'step , '
  for i in range(nbRecord):
    for N in NUCLEI:
      frstr += N+' ('+str(i)+') , '
  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','w')
  firingRatesFile.writelines(frstr)

  #----------------------------------
  # Loop over the 5 steps of the test
  #----------------------------------
  for timeStep in range(5):
    #-------------------------
    # measures                                                                                                                                                  
    #-------------------------
    spkDetect=[{},{},{}] # list of spike detector dictionaries used to record the experiment in the first 3 channels

    antagStr = ''
    if antagInjectionSite != 'none':
      antagStr = antagInjectionSite+'_'+antag+'_'

    for i in range(nbRecord):
      for N in NUCLEI:
        spkDetect[i][N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": str(timeStep)+'_'+antagStr+N, "to_file": True, 'start':offsetDuration + timeStep*(offsetDuration+simDuration),'stop':(timeStep+1)*(offsetDuration+simDuration)})
        nest.Connect(Pop[N][i], spkDetect[i][N])

    frstr = str(timeStep) + ', '

    #-------------------------
    # Simulation
    #-------------------------
    print '====== Step',timeStep,'======'
    print 'Channel 0:',CSNrate[0,timeStep],PTNrate[0,timeStep]
    print 'Channel 1:',CSNrate[1,timeStep],PTNrate[1,timeStep]

    '''
    if 'Fake' in globals():
      if 'CSN' in Fake:
        nest.SetStatus(Fake['CSN'][0],{'rate':CSNrate[0,timeStep]})
        nest.SetStatus(Fake['CSN'][1],{'rate':CSNrate[1,timeStep]})
      else:
        nest.SetStatus(Pop['CSN'][0],{'rate':CSNrate[0,timeStep]})
        nest.SetStatus(Pop['CSN'][1],{'rate':CSNrate[1,timeStep]})
      if 'PTN' in Fake:
        nest.SetStatus(Fake['PTN'][0],{'rate':PTNrate[0,timeStep]})
        nest.SetStatus(Fake['PTN'][1],{'rate':PTNrate[1,timeStep]})
      else:
        nest.SetStatus(Pop['PTN'][0],{'rate':PTNrate[0,timeStep]})
        nest.SetStatus(Pop['PTN'][1],{'rate':PTNrate[1,timeStep]})
    else:
      nest.SetStatus(Pop['CSN'][0],{'rate':CSNrate[0,timeStep]})
      nest.SetStatus(Pop['CSN'][1],{'rate':CSNrate[1,timeStep]})
      nest.SetStatus(Pop['PTN'][0],{'rate':PTNrate[0,timeStep]})
      nest.SetStatus(Pop['PTN'][1],{'rate':PTNrate[1,timeStep]})
    '''

    nest.SetStatus(ActPop['CSN'][0],{'rate':CSNrate[0,timeStep]})
    nest.SetStatus(ActPop['CSN'][1],{'rate':CSNrate[1,timeStep]})
    nest.SetStatus(ActPop['PTN'][0],{'rate':PTNrate[0,timeStep]})
    nest.SetStatus(ActPop['PTN'][1],{'rate':PTNrate[1,timeStep]})

    nest.Simulate(simDuration+offsetDuration)

    for i in range(nbRecord):
      print '------ Channel',i,'-------'
      for N in NUCLEI:
        #strTestPassed = 'NO!'
        expeRate[N][i,timeStep] = nest.GetStatus(spkDetect[i][N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
        print 't('+str(timeStep)+')',N,':',expeRate[N][i,timeStep],'Hz'
        frstr += '%f , ' %(expeRate[N][i,timeStep])

    strTestPassed = 'YES!'
    if timeStep == 0:
      for i in range(params['nbCh']):
        if expeRate['GPi'][0,timeStep]<FRRNormal['GPi'][0]:
          strTestPassed = 'NO!'
      meanRestGPi = expeRate['GPi'][:,timeStep].mean()
    elif timeStep == 1:
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9:
        strTestPassed = 'NO!'
    elif timeStep == 2:
      if expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][0,timeStep] < expeRate['GPi'][1,timeStep]:
        strTestPassed = 'NO!'
    elif timeStep == 3:
      if expeRate['GPi'][0,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 :
        strTestPassed = 'NO!'
    elif timeStep == 4:
      if expeRate['GPi'][1,timeStep] > meanRestGPi*0.9 or expeRate['GPi'][0,timeStep] < expeRate['GPi'][1,timeStep]:
        strTestPassed = 'NO!'

    if strTestPassed == 'YES!':
      score +=1

    print '------ Result ------'
    print expeRate['GPi'][0,timeStep],'Hz',expeRate['GPi'][1,timeStep],'Hz',strTestPassed

    # write measured firing rates in csv file
    frstr+='\n'
    firingRatesFile.writelines(frstr)

    #-------------------------
    # Displays
    #-------------------------
    '''
    if showRasters and interactive:
      for i in range(nbRecord):
        displayStr = ' Channel '+str(i)
        displayStr+=' ('+antagStr[:-1]+')' if (antagInjectionSite != 'none') else ''
        #for N in NUCLEI:
        for N in ['MSN','STN']:
          #nest.raster_plot.from_device(spkDetect[i][N],hist=True,title=N+displayStr)
          nest.raster_plot.from_device(spkDetect[i][N],hist=False,title=N+displayStr)

      nest.raster_plot.show()
    '''

  if showRasters and interactive:
    for N in NUCLEI:
      nest.raster_plot.from_device(inspector[N],hist=True,title=N)
    nest.raster_plot.show()

  firingRatesFile.close()

  return score,5

#-----------------------------------------------------------------------
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
  
  #execTime = time.localtime()
  #timeStr = str(execTime[0])+'_'+str(execTime[1])+'_'+str(execTime[2])+'_'+str(execTime[3])+':'+str(execTime[4])+':'+str(execTime[5])

  score = np.zeros((2))
  '''
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=True)

  for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
    score += checkAvgFR(params=params,antagInjectionSite='GPe',antag=a)

  for a in ['All','AMPA','NMDA+AMPA','NMDA','GABAA']:
    score += checkAvgFR(params=params,antagInjectionSite='GPi',antag=a)
  '''

  score += checkGurneyTest(showRasters=True,params=params,PActiveCSN=0.2,PActivePTN=0.2)
  #score += checkGurneyTest(showRasters=True,params=params,PActiveCSN=1.,PActivePTN=1.)

  #-------------------------
  print "******************"
  print "* Score:",score[0],'/',score[1]
  print "******************"

  #-------------------------
  # log the results in a file
  #-------------------------
  #res = open('log/OutSummary_'+timeStr+'.txt','a')
  res = open('log/OutSummary.txt','a')
  for k,v in params.iteritems():
    res.writelines(k+' , '+str(v)+'\n')
  res.writelines("Score: "+str(score[0])+' , '+str(score[1]))
  res.close()

  res = open('score.txt','w')
  res.writelines(str(score[0])+'\n')
  res.close()

#---------------------------
if __name__ == '__main__':
  main()
