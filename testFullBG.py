#!/apps/free/python/2.7.10/bin/python
from LGneurons import *
from modelParams import *
import nest.raster_plot
import nest.voltage_trace
import pylab as pl
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
  nest.SetStatus(Pop['MSN'],{"I_e":params['IeMSN']})

  nbSim['FSI'] = params['nbFSI']
  create('FSI')
  nest.SetStatus(Pop['FSI'],{"I_e":params['IeFSI']})

  nbSim['STN'] = params['nbSTN']
  create('STN')
  nest.SetStatus(Pop['STN'],{"I_e":params['IeSTN']})

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
  connect('ex','CSN','MSN', inDegree= params['inDegCSNMSN'], gain=G['MSN'])
  connect('ex','PTN','MSN', inDegree= params['inDegPTNMSN'], gain=G['MSN'])
  connect('ex','CMPf','MSN', inDegree=params['inDegCMPfMSN'],gain=G['MSN'])
  connect('in','MSN','MSN', inDegree= params['inDegMSNMSN'], gain=G['MSN'])
  connect('in','FSI','MSN', inDegree= params['inDegFSIMSN'], gain=G['MSN'])
  # some parameterizations from LG14 have no STN->MSN or GPe->MSN synaptic contacts
  if alpha['STN->MSN'] != 0:
    print "alpha['STN->MSN']",alpha['STN->MSN']
    connect('ex','STN','MSN', inDegree= params['inDegSTNMSN'],gain=G['MSN'])
  if alpha['GPe->MSN'] != 0:
    print "alpha['GPe->MSN']",alpha['GPe->MSN']
    connect('in','GPe','MSN', inDegree= params['inDegGPeMSN'],gain=G['MSN'])

  print '* FSI Inputs'
  connect('ex','CSN','FSI',  inDegree= params['inDegCSNFSI'], gain=G['FSI'])
  connect('ex','PTN','FSI',  inDegree= params['inDegPTNFSI'], gain=G['FSI'])
  if alpha['STN->FSI'] != 0:
    connect('ex','STN','FSI',inDegree= params['inDegSTNFSI'], gain=G['FSI'])
  connect('in','GPe','FSI',  inDegree= params['inDegGPeFSI'], gain=G['FSI'])
  connect('ex','CMPf','FSI', inDegree= params['inDegCMPfFSI'],gain=G['FSI'])
  connect('in','FSI','FSI',  inDegree= params['inDegFSIFSI'], gain=G['FSI'])

  print '* STN Inputs'
  connect('ex','PTN','STN', inDegree= params['inDegPTNSTN'], gain=G['STN'])
  connect('ex','CMPf','STN',inDegree= params['inDegCMPfSTN'],gain=G['STN'])
  connect('in','GPe','STN', inDegree= params['inDegGPeSTN'], gain=G['STN']) 

  print '* GPe Inputs'
  if antagInjectionSite == 'GPe':
    if   antag == 'AMPA':
      connect('NMDA','CMPf','GPe',inDegree=params['inDegCMPfGPe'],gain=G['GPe'])
      connect('NMDA','STN','GPe', inDegree=params['inDegSTNGPe'], gain=G['GPe'])
      connect('in','MSN','GPe', inDegree=  params['inDegMSNGPe'], gain=G['GPe'])
      connect('in','GPe','GPe', inDegree=  params['inDegGPeGPe'], gain=G['GPe']) 
    elif antag == 'NMDA':
      connect('AMPA','CMPf','GPe',inDegree= params['inDegCMPfGPe'],gain=G['GPe'])
      connect('AMPA','STN','GPe', inDegree= params['inDegSTNGPe'], gain=G['GPe'])
      connect('in','MSN','GPe', inDegree= params['inDegMSNGPe'],   gain=G['GPe'])
      connect('in','GPe','GPe', inDegree= params['inDegGPeGPe'],   gain=G['GPe'])
    elif antag == 'AMPA+GABAA':
      connect('NMDA','CMPf','GPe',inDegree= params['inDegCMPfGPe'],gain=G['GPe'])
      connect('NMDA','STN','GPe',inDegree= params['inDegSTNGPe'],  gain=G['GPe'])
    elif antag == 'GABAA':
      connect('ex','CMPf','GPe',inDegree= params['inDegCMPfGPe'], gain=G['GPe'])
      connect('ex','STN','GPe', inDegree= params['inDegSTNGPe'],  gain=G['GPe'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connect('ex','CMPf','GPe',inDegree= params['inDegCMPfGPe'],gain=G['GPe'])
    connect('ex','STN','GPe', inDegree= params['inDegSTNGPe'], gain=G['GPe'])
    connect('in','MSN','GPe', inDegree= params['inDegMSNGPe'], gain=G['GPe'])
    connect('in','GPe','GPe', inDegree= params['inDegGPeGPe'], gain=G['GPe'])

  print '* GPi Inputs'
  if antagInjectionSite =='GPi':
    if   antag == 'AMPA+NMDA+GABAA':
      pass
    elif antag == 'NMDA':
      connect('in','MSN','GPi',   inDegree= params['inDegMSNGPi'], gain=G['GPi'])
      connect('AMPA','STN','GPi', inDegree= params['inDegSTNGPi'], gain=G['GPi'])
      connect('in','GPe','GPi',   inDegree= params['inDegGPeGPi'], gain=G['GPi'])
      connect('AMPA','CMPf','GPi',inDegree= params['inDegCMPfGPi'],gain=G['GPi'])
    elif antag == 'NMDA+AMPA':
      connect('in','MSN','GPi', inDegree= params['inDegMSNGPi'], gain=G['GPi'])
      connect('in','GPe','GPi', inDegree= params['inDegGPeGPi'], gain=G['GPi'])
    elif antag == 'AMPA':
      connect('in','MSN','GPi',   inDegree= params['inDegMSNGPi'], gain=G['GPi'])
      connect('NMDA','STN','GPi', inDegree= params['inDegSTNGPi'], gain=G['GPi'])
      connect('in','GPe','GPi',   inDegree= params['inDegGPeGPi'], gain=G['GPi'])
      connect('NMDA','CMPf','GPi',inDegree= params['inDegCMPfGPi'],gain=G['GPi'])
    elif antag == 'GABAA':
      connect('ex','STN','GPi', inDegree= params['inDegSTNGPi'], gain=G['GPi'])
      connect('ex','CMPf','GPi',inDegree= params['inDegCMPfGPi'],gain=G['GPi'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connect('in','MSN','GPi', inDegree= params['inDegMSNGPi'], gain=G['GPi'])
    connect('ex','STN','GPi', inDegree= params['inDegSTNGPi'], gain=G['GPi'])
    connect('in','GPe','GPi', inDegree= params['inDegGPeGPi'], gain=G['GPi'])
    connect('ex','CMPf','GPi',inDegree= params['inDegCMPfGPi'],gain=G['GPi'])

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
  if 'nbcpu' in params:
    nest.SetKernelStatus({'local_num_threads': params['nbcpu']})
  nest.SetKernelStatus({"data_path": dataPath})
  initNeurons()

  offsetDuration = 1000.
  simDuration = 5000. # ms
  # nest.SetKernelStatus({"overwrite_files":True}) # Thanks to use of timestamps, file names should now 
                                                   # be different as long as they are not created during the same second

  print '/!\ Using the following LG14 parameterization',params['LG14modelID']
  loadLG14params(params['LG14modelID'])
  loadThetaFromCustomparams(params)

  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeMSN','IeFSI','IeSTN','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','inDegCSNMSN','inDegPTNMSN','inDegCMPfMSN','inDegMSNMSN','inDegFSIMSN','inDegSTNMSN','inDegGPeMSN','inDegCSNFSI','inDegPTNFSI','inDegSTNFSI','inDegGPeFSI','inDegCMPfFSI','inDegFSIFSI','inDegPTNSTN','inDegCMPfSTN','inDegGPeSTN','inDegCMPfGPe','inDegSTNGPe','inDegMSNGPe','inDegGPeGPe','inDegMSNGPi','inDegSTNGPi','inDegGPeGPi','inDegCMPfGPi',]
  for np in necessaryParams:
    if np not in params:
      raise KeyError('Missing parameter: '+np)

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
  frstr = "#" + str(params['LG14modelID'])+ " , " + antagInjectionSite + ', '
  s = '----- RESULTS -----'
  print s
  text.append(s+'\n')
  if antagInjectionSite == 'none':
    validationStr = "\n#" + str(params['LG14modelID']) + " , "
    frstr += "none , "
    for N in NUCLEI:
      strTestPassed = 'NO!'
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
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
  else:
    validationStr = ""
    frstr += str(antag) + " , "
    for N in NUCLEI:
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
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
  
  score = np.zeros((2))
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=True)


  for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
    score += checkAvgFR(params=params,antagInjectionSite='GPe',antag=a)

  for a in ['AMPA+NMDA+GABAA','AMPA','NMDA+AMPA','NMDA','GABAA']:
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
