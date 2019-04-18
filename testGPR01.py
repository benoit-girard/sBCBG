#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
## testGPR01.py
##
## This script implements the selection test from Gurney, Prescott and Redrave, 2001b

from iniBG import *


#-----------------------------------------------------------------------
# PActiveCNS/PTN : proportion of "active" neurons in the CSN/PTN populations (in [0.,1.])
#
#-----------------------------------------------------------------------
def checkGurneyTest(showRasters=False,params={},CSNFR=[2.,10.], PActiveCSN=1., PTNFR=[15.,35], PActivePTN=1., antagInjectionSite='none',antag='', CMPfFR=[4., 4.], PActiveCMPf=1., transientCMPf=0.1):

  nest.ResetNetwork()
  initNeurons()

  dataPath='log/'
  nest.SetKernelStatus({"overwrite_files":True}) # when we redo the simulation, we erase the previous traces

  offsetDuration = 200.
  simDuration = 800. # ms

  nbRecord=2 # number of channels whose activity will be recorded
  if params['nbCh']<2:
    print 'need at least 2 channels to perform Gurney test'
    exit()
  elif params['nbCh']>2:
    nbRecord = 3 # if we have more than 2 channels, we will also record one of the neutral channels

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
  def activate_pop(N, PActive, nbCh=2):
    src = Pop[n]
    if 'Fake' in globals():
      if N in Fake:
        src = Fake[N]
    if PActive==1.:
     ActPop[N]=src
    else:
      for i in range(nbCh):
        ActPop[N][i] = tuple(rnd.choice(a=np.array(src[i]),size=int(nbSim[N]*PActive),replace=False))
  activate_pop('CSN', PActiveCSN)
  activate_pop('PTN', PActivePTN)
  activate_pop('CMPf', PActiveCMPf, nbCh=params['nbCh'])

  #-------------------------
  # log-related variales
  #-------------------------
  score = 0
  expeRate={}
  for N in NUCLEI:
    expeRate[N]=-1. * np.ones((nbRecord,5))

  antagStr = ''
  if antagInjectionSite != 'none':
    antagStr = antagInjectionSite+'_'+antag+'_'

  inspector = {}
  for N in NUCLEI:
    inspector[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": antagStr+N, "to_file": True, 'start':0. ,'stop':offsetDuration+simDuration*6})
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

  #-------------------------
  # measures
  #-------------------------
  spkDetect=[{},{},{}] # list of spike detector dictionaries used to record the experiment in the first 3 channels


  for i in range(nbRecord):
    for N in NUCLEI:
      spkDetect[i][N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": antagStr+N, "to_file": False, 'start':offsetDuration ,'stop':offsetDuration+simDuration*6})
      nest.Connect(Pop[N][i], spkDetect[i][N])

  # stimulation without inputs, to make sure that the NN is at rest
  nest.Simulate(offsetDuration)

  #----------------------------------
  # Loop over the 5 steps of the test
  #----------------------------------
  for timeStep in range(5):

    frstr = str(timeStep) + ', '

    #-------------------------
    # Simulation
    #-------------------------
    print '====== Step',timeStep,'======'
    print 'Channel 0:',CSNrate[0,timeStep],PTNrate[0,timeStep]
    print 'Channel 1:',CSNrate[1,timeStep],PTNrate[1,timeStep]

    nest.SetStatus(ActPop['CSN'][0],{'rate':CSNrate[0,timeStep]})
    nest.SetStatus(ActPop['CSN'][1],{'rate':CSNrate[1,timeStep]})
    nest.SetStatus(ActPop['PTN'][0],{'rate':PTNrate[0,timeStep]})
    nest.SetStatus(ActPop['PTN'][1],{'rate':PTNrate[1,timeStep]})

    if PActiveCMPf == 0. or CMPfFR[0] == CMPfFR[1] or transientCMPf == 0.:
      # vanilla GPR01 scenario
      nest.Simulate(simDuration)
    else:
      # CMPf activity during selection
      print('CMPf activity increased to ' + str(CMPfFR[1]) + ' for ' + str((simDuration+offsetDuration)*transientCMPf) + ' ms\n')
      for Ch in range(params['nbCh']):
        nest.SetStatus(ActPop['CMPf'][Ch],{'rate': CMPfFR[1]})
      nest.Simulate(simDuration*transientCMPf)
      for Ch in range(params['nbCh']):
        nest.SetStatus(ActPop['CMPf'][Ch],{'rate': CMPfFR[0]})
      nest.Simulate(simDuration*(1.-transientCMPf))
        

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

  print('Summary of activity')
  import pprint
  pprint.pprint(expeRate['GPi'])

  if showRasters and interactive:
    for N in NUCLEI:
      nest.raster_plot.from_device(inspector[N],hist=True,title=N)
    nest.raster_plot.show()

  firingRatesFile.close()

  return score,5

#-----------------------------------------------------------------------
# PActiveCNS/PTN : proportion of "active" neurons in the CSN/PTN populations (in [0.,1.])
#
#-----------------------------------------------------------------------
def checkGeorgopoulosTest(showRasters=False,params={},CSNFR=[2.,10.], PActiveCSN=1., PTNFR=[15.,35], PActivePTN=1., antagInjectionSite='none',antag=''):

  offsetDuration = 500.
  simDuration = 1000. # ms
  loadLG14params(params['LG14modelID'])

  #-------------------------
  # prepare the firing rates of the inputs
  #-------------------------
  gCSN = CSNFR[1]-CSNFR[0]
  gPTN = PTNFR[1]-PTNFR[0]
  activityLevels = np.ones((params['nbCh']))
  for i in range(params['nbCh']):
    activityLevels[i] = 2 * np.pi / params['nbCh'] * i
  activityLevels = (np.cos(activityLevels)+1)/2.

  CSNrate= gCSN * activityLevels + np.ones((params['nbCh'])) * CSNFR[0]
  PTNrate= gPTN * activityLevels + np.ones((params['nbCh'])) * PTNFR[0]

  #-------------------------
  # and prepare the lists of neurons that will be affected by these activity changes
  #-------------------------
  ActPop = {'CSN':[() for i in range(params['nbCh'])],'PTN':[() for i in range(params['nbCh'])]}
  if 'Fake' in globals():
    if 'CSN' in Fake:
      if PActiveCSN==1.:
       ActPop['CSN']=Fake['CSN']
      else:
        for i in range(params['nbCh']):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Fake['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    else:
      if PActiveCSN==1.:
       ActPop['CSN']=Pop['CSN']
      else:
        for i in range(params['nbCh']):
          ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if 'PTN' in Fake:
      if PActivePTN==1.:
        ActPop['PTN']=Fake['PTN']
      else:
        for i in range(params['nbCh']):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Fake['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
    else:
      if PActivePTN==1.:
        ActPop['PTN']=Pop['PTN']
      else:
        for i in range(params['nbCh']):
          ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))
  else:
    if PActiveCSN==1.:
     ActPop['CSN']=Pop['CSN']
    else:
      for i in range(params['nbCh']):
        ActPop['CSN'][i] = tuple(rnd.choice(a=np.array(Pop['CSN'][i]),size=int(nbSim['CSN']*PActiveCSN),replace=False))
    if PActivePTN==1.:
      ActPop['PTN']=Pop['PTN']
    else:
      for i in range(params['nbCh']):
        ActPop['PTN'][i] = tuple(rnd.choice(a=np.array(Pop['PTN'][i]),size=int(nbSim['PTN']*PActivePTN),replace=False))

  #-------------------------
  # log-related variables
  #-------------------------
  GPiRestRate= -1*np.ones((params['nbCh']))
  expeRate={}
  for N in NUCLEI:
    expeRate[N]=-1. * np.ones((params['nbCh']))

  inspector = {}
  for N in NUCLEI:
    inspector[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": N, "to_file": False})
    for i in range(params['nbCh']):
      nest.Connect(Pop[N][i],inspector[N])

  #-------------------------
  # write header in firingRate summary file
  #-------------------------
  frstr = 'channel , '
  for N in NUCLEI:
    frstr += N+', '
  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','w')
  firingRatesFile.writelines(frstr)

  #-------------------------
  # measures
  #-------------------------
  spkDetect=[{} for i in range(params['nbCh'])] # list of spike detector dictionaries used to record the experiment in all the channels

  antagStr = ''
  if antagInjectionSite != 'none':
    antagStr = antagInjectionSite+'_'+antag+'_'

  for i in range(params['nbCh']):
    for N in NUCLEI:
      spkDetect[i][N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": antagStr+N, "to_file": True, 'start':2*offsetDuration+simDuration,'stop':2*(offsetDuration+simDuration)})
      nest.Connect(Pop[N][i], spkDetect[i][N])

  GPiRestSpkDetect = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": antagStr+'GPiRest', "to_file": True, 'start':offsetDuration,'stop':offsetDuration+simDuration})
  for i in range(params['nbCh']):
      nest.Connect(Pop['GPi'][i], GPiRestSpkDetect)

  #-------------------------
  # Simulation
  #-------------------------

  # step 1 : stimulation without inputs, to calibrate the GPi activity at rest :
  nest.Simulate(simDuration+offsetDuration)

  print '------ Rest Period ------'  
  frstr = 'rest, , , , ,' # only GPi is recorded at rest, and on all channels
  GPiRestRate = nest.GetStatus(GPiRestSpkDetect, 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
  print "GPi rate at rest:",GPiRestRate;"Hz"
  frstr += '%f \n' %GPiRestRate
  firingRatesFile.writelines(frstr)

  for i in range(params['nbCh']):
    nest.SetStatus(ActPop['CSN'][i],{'rate':CSNrate[i]})
    nest.SetStatus(ActPop['PTN'][i],{'rate':PTNrate[i]})

  nest.Simulate(simDuration+offsetDuration)

  for i in range(params['nbCh']):
    print '------ Channel',i,'------'
    frstr = str(i)+', '
    for N in NUCLEI:
      expeRate[N][i] = nest.GetStatus(spkDetect[i][N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
      print N,':',expeRate[N][i],'Hz'
      frstr += '%f , ' %(expeRate[N][i])
    frstr += '\n'

    firingRatesFile.writelines(frstr)

  firingRatesFile.close()

  #-------------------------
  # Displays
  #-------------------------
  if showRasters and interactive:
    for N in NUCLEI:
      nest.raster_plot.from_device(inspector[N],hist=True,title=N)
    nest.raster_plot.show()

    pylab.plot(expeRate['GPi'])
    pylab.show()

  contrast = 2. / CSNrate * GPiRestRate / expeRate['GPi']

  return contrast

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
  
  instantiate_BG(params, antagInjectionSite='none', antag='')
  score = np.zeros((2))

  #proportion = 0.1
  proportion = 0.5
  score += checkGurneyTest(showRasters=True,params=params,PActiveCSN=proportion,PActivePTN=proportion,CSNFR=[2.,4.], PTNFR=[15.,15], CMPfFR=[4., 4.], transientCMPf=1.)
  #score += checkGeorgopoulosTest(params=params,PActiveCSN=proportion,PActivePTN=proportion)
  # Test Liénard-style:
 # score += checkGeorgopoulosTest(params=params,CSNFR=[2.,4.], PTNFR=[15.,15])
  #score += checkGurneyTest(showRasters=True,params=params,PActiveCSN=1.,PActivePTN=1.)

  #-------------------------
  print "******************"
  #print "* Score:",score[0],'/',score[1]
  print score
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
