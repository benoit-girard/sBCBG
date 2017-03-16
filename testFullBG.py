#!/apps/free/python/2.7.10/bin/python
from LGneurons import *
from modelParams import *
import nest.raster_plot
#import time
import sys


#------------------------------------------
# Creates the populations of neurons necessary to simulate a BG circuit
#------------------------------------------
def createBG():
  #==========================
  # Creation of neurons
  #-------------------------
  print '\nCreating neurons\n================'

  nbSim['MSN'] = params['nbMSN']
  create('MSN')

  nbSim['FSI'] = params['nbFSI']
  create('FSI')

  nbSim['STN'] = params['nbSTN']
  create('STN')

  nbSim['GPe'] = params['nbGPe']
  create('GPe')
  nest.SetStatus(Pop['GPe'],{"I_e":params['IeGPe']})

  nbSim['GPi'] = params['nbGPi']
  create('GPi')
  nest.SetStatus(Pop['GPi'],{"I_e":params['IeGPi']})

  parrot = True # switch to False at your risks & perils...                                                                                                                                     
  nbSim['CSN'] = params['nbCSN']
  create('CSN', fake=True, parrot=parrot)

  nbSim['PTN'] = params['nbPTN']
  create('PTN', fake=True, parrot=parrot)

  nbSim['CMPf'] = params['nbCMPf']
  create('CMPf', fake=True, parrot=False)

  print "Number of simulated neurons:", nbSim

#------------------------------------------
# Connects the populations of a previously created multi-channel BG circuit
#------------------------------------------
def connectBG(antagInjectionSite,antag):
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
  connect('ex','CSN','MSN', inDegree= min(params['inDegCSNMSN'],nbSim['CSN']),gain=G['MSN'])
  connect('ex','PTN','MSN', inDegree= min(params['inDegPTNMSN'],nbSim['PTN']),gain=G['MSN'])
  connect('ex','CMPf','MSN', inDegree= min(params['inDegCMPfMSN'],nbSim['CMPf']),gain=G['MSN'])
  connect('in','MSN','MSN', inDegree= min(params['inDegMSNMSN'],nbSim['MSN']),gain=G['MSN'])
  connect('in','FSI','MSN', inDegree= min(params['inDegFSIMSN'],nbSim['FSI']),gain=G['MSN'])
  # some parameterizations from LG14 have no STN->MSN or GPe->MSN synaptic contacts
  if alpha['STN->MSN'] != 0:
    print "alpha['STN->MSN']",alpha['STN->MSN']
    connect('ex','STN','MSN', inDegree= min(params['inDegSTNMSN'],nbSim['STN']),gain=G['MSN'])
  if alpha['GPe->MSN'] != 0:
    print "alpha['GPe->MSN']",alpha['GPe->MSN']
    connect('in','GPe','MSN', inDegree= min(params['inDegGPeMSN'],nbSim['GPe']),gain=G['MSN']) # diffuse ? focused ?                                                               

  print '* FSI Inputs'
  connect('ex','CSN','FSI',  inDegree= min(params['inDegCSNFSI'],nbSim['CSN']),gain=G['FSI'])
  connect('ex','PTN','FSI',  inDegree= min(params['inDegPTNFSI'],nbSim['PTN']),gain=G['FSI'])
  if alpha['STN->FSI'] != 0:
    connect('ex','STN','FSI',inDegree= min(params['inDegSTNFSI'],nbSim['STN']),gain=G['FSI'])
  connect('in','GPe','FSI',  inDegree= min(params['inDegGPeFSI'],nbSim['GPe']),gain=G['FSI'])
  connect('ex','CMPf','FSI', inDegree= min(params['inDegCMPfFSI'],nbSim['CMPf']),gain=G['FSI'])
  connect('in','FSI','FSI',  inDegree= min(params['inDegFSIFSI'],nbSim['FSI']),gain=G['FSI'])

  print '* STN Inputs'
  connect('ex','PTN','STN', inDegree= min(params['inDegPTNSTN'],nbSim['PTN']), gain=G['STN'])
  connect('ex','CMPf','STN',inDegree= min(params['inDegCMPfSTN'],nbSim['CMPf']), gain=G['STN'])
  connect('in','GPe','STN', inDegree= min(params['inDegGPeSTN'],nbSim['GPe']), gain=G['STN']) # or diffuse, to be in line with the 2008 model?                                      

  print '* GPe Inputs'
  if antagInjectionSite == 'GPe':
    if   antag == 'AMPA':
      connect('NMDA','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connect('NMDA','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
      connect('in','MSN','GPe', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
      connect('in','GPe','GPe', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe']) # diffuse or focused?                                                             
    elif antag == 'NMDA':
      connect('AMPA','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connect('AMPA','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
      connect('in','MSN','GPe', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
      connect('in','GPe','GPe', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe'])
    elif antag == 'AMPA+GABAA':
      connect('NMDA','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connect('NMDA','STN','GPe',inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    elif antag == 'GABAA':
      connect('ex','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
      connect('ex','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connect('ex','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']), gain=G['GPe'])
    connect('ex','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']), gain=G['GPe'])
    connect('in','MSN','GPe', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']), gain=G['GPe'])
    connect('in','GPe','GPe', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']), gain=G['GPe'])

  print '* GPi Inputs'
  if antagInjectionSite =='GPi':
    if   antag == 'All':
      pass
    elif antag == 'NMDA':
      connect('in','MSN','GPi',   inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connect('AMPA','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connect('in','GPe','GPi',   inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
      connect('AMPA','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    elif antag == 'NMDA+AMPA':
      connect('in','MSN','GPi', inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connect('in','GPe','GPi', inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
    elif antag == 'AMPA':
      connect('in','MSN','GPi',   inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
      connect('NMDA','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connect('in','GPe','GPi',   inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
      connect('NMDA','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    elif antag == 'GABAA':
      connect('ex','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
      connect('ex','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connect('in','MSN','GPi', inDegree= min(params['inDegMSNGPi'],nbSim['MSN']),gain=G['GPi'])
    connect('ex','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']),gain=G['GPi'])
    connect('in','GPe','GPi', inDegree= min(params['inDegGPeGPi'],nbSim['GPe']),gain=G['GPi'])
    connect('ex','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']),gain=G['GPi'])

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
  #------------------------
  # creation and connection of the neural populations
  #------------------------

  createBG()

  connectBG(antagInjectionSite,antag)

  #-------------------------
  # measures
  #-------------------------
  spkDetect={} # spike detectors used to record the experiment
  expeRate={}

  #if logFileName == '':
  #  execTime = time.localtime()
  #  logFileName = str(execTime[0])+'_'+str(execTime[1])+'_'+str(execTime[2])+'_'+str(execTime[3])+':'+str(execTime[4])+':'+str(execTime[5])

  antagStr = ''
  if antagInjectionSite != 'none':
    antagStr = antagInjectionSite+'_'+antag+'_'

  for N in NUCLEI:
    # 1000ms offset period for network stabilization
    #spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": logFileName+'_'+antagStr+N, "to_file": True, 'start':offsetDuration,'stop':offsetDuration+simDuration})
    spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": antagStr+N, "to_file": True, 'start':offsetDuration,'stop':offsetDuration+simDuration})
    nest.Connect(Pop[N], spkDetect[N])

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
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
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
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
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
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='')#,showRasters=True)


  for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
    score += checkAvgFR(params=params,antagInjectionSite='GPe',antag=a)

  for a in ['All','AMPA','NMDA+AMPA','NMDA','GABAA']:
    score += checkAvgFR(params=params,antagInjectionSite='GPi',antag=a)

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
