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
  createMC('CMPf',params['nbCh'], fake=True, parrot=parrot)

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

  #IMPROVE
  params['nbCh']=8

  score = np.zeros((2))
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=True)

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
