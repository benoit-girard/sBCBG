from LGneurons import *
import nest.raster_plot
import time

# params possible keys:
# - nb{MSN,FSI,STN,GPi,GPe,CSN,PTN,CMPf} : number of simulated neurons for each population
# - Ie{GPe,GPi} : constant input current to GPe and GPi
# - G{MSN,FSI,STN,GPi,GPe} : gain to be applied on LG14 input synaptic weights for each population

def checkAvgFR(showRasters=True,params={},antagInjectionSite='none',antag=''):
  nest.ResetKernel()
  nest.SetKernelStatus({'local_num_threads':2, "data_path": "log/"})
  initNeurons()

  simDuration = 3000. # ms
  # nest.SetKernelStatus({"overwrite_files":True}) # Thanks to use of timestamps, file names should now 
                                                   # be different as long as they are not created during the same second

  #==========================
  # Creation of neurons
  #-------------------------
  print '\nCreating neurons\n================'

  nbSim['MSN'] = params['nbMSN'] if ('nbMSN' in params) else 2644.
  create('MSN')

  nbSim['FSI'] = params['nbFSI'] if ('nbFSI' in params) else 53.
  create('FSI')

  nbSim['STN'] = params['nbSTN'] if ('nbSTN' in params) else 8.
  create('STN')

  nbSim['GPe'] = params['nbGPe'] if ('nbGPe' in params) else 25.
  create('GPe')
  nest.SetStatus(Pop['GPe'],{"I_e":params['IeGPe'] if ('IeGPe' in params) else 13.})

  nbSim['GPi'] = params['nbGPi'] if ('nbGPi' in params) else 14.
  create('GPi')
  nest.SetStatus(Pop['GPi'],{"I_e":params['IeGPi'] if ('IeGPi' in params) else 9.})
  
  nbSim['CSN'] = params['nbCSN'] if ('nbCSN' in params) else 3000
  create('CSN', fake=True)
  
  nbSim['PTN'] = params['nbPTN'] if ('nbPTN' in params) else 100.
  create('PTN', fake=True)

  nbSim['CMPf'] = params['nbCMPf'] if ('nbCMPf' in params) else 9.
  create('CMPf', fake=True)

  print "Number of simulated neurons:", nbSim

  G = {'MSN': params['GMSN'] if ('GMSN' in params) else 3.9, # 3.9
       'FSI': params['GFSI'] if ('GFSI' in params) else 1.1,
       'STN': params['GSTN'] if ('GSTN' in params) else 1.35,
       'GPe': params['GGPe'] if ('GGPe' in params) else 1.3,
       'GPi': params['GGPi'] if ('GGPi' in params) else 1.3,
      }

  print "Gains on LG14 syn. strength:", G

  #-------------------------
  # connection of populations
  #-------------------------
  print '\nConnecting neurons\n================'
  print "**",antag,"antagonist injection in",antagInjectionSite,"**"
  print '* MSN Inputs'
  print 'DEBUG',nbSim['CSN'], min(params['inDegCSNMSN'],nbSim['CSN']) if ('inDegCSNMSN' in params) else min(100,nbSim['CSN']), G['MSN']
  connect('ex','CSN','MSN', inDegree= min(params['inDegCSNMSN'],nbSim['CSN']) if ('inDegCSNMSN' in params) else min(100,nbSim['CSN']),gain=G['MSN'])
  connect('ex','PTN','MSN', inDegree= min(params['inDegPTNMSN'],nbSim['PTN']) if ('inDegPTNMSN' in params) else 1,gain=G['MSN'])
  connect('ex','CMPf','MSN',inDegree= min(params['inDegCMPfMSN'],nbSim['CMPf']) if ('inDegCMPfMSN' in params) else 1,gain=G['MSN'])
  connect('in','MSN','MSN', inDegree= min(params['inDegMSNMSN'],nbSim['MSN']) if ('inDegMSNMSN' in params) else 50,gain=G['MSN'])
  connect('in','FSI','MSN', inDegree= min(params['inDegFSIMSN'],nbSim['FSI']) if ('inDegFSIMSN' in params) else 1,gain=G['MSN'])
  # some parameterizations from LG14 have no STN->MSN or GPe->MSN synaptic contacts
  if alpha['STN->MSN'] != 0:
    connect('ex','STN','MSN', inDegree= min(params['inDegSTNMSN'],nbSim['STN']) if ('inDegSTNMSN' in params) else 1,gain=G['MSN'])
  if alpha['GPe->MSN'] != 0:
    connect('in','GPe','MSN', inDegree= min(params['inDegGPeMSN'],nbSim['GPe']) if ('inDegGPeMSN' in params) else 1,gain=G['MSN'])

  print '* FSI Inputs'
  connect('ex','CSN','FSI', inDegree= min(params['inDegCSNFSI'],nbSim['CSN']) if ('inDegCSNFSI' in params) else 50,gain=G['FSI'])
  connect('ex','PTN','FSI', inDegree= min(params['inDegPTNFSI'],nbSim['PTN']) if ('inDegPTNFSI' in params) else 1,gain=G['FSI'])
  connect('ex','STN','FSI', inDegree= min(params['inDegSTNFSI'],nbSim['STN']) if ('inDegSTNFSI' in params) else 2,gain=G['FSI'])
  connect('in','GPe','FSI', inDegree= min(params['inDegGPeFSI'],nbSim['GPe']) if ('inDegGPeFSI' in params) else nbSim['GPe'],gain=G['FSI'])
  connect('ex','CMPf','FSI',inDegree= min(params['inDegCMPfFSI'],nbSim['CMPf']) if ('inDegCMPfFSI' in params) else nbSim['CMPf'],gain=G['FSI'])
  connect('in','FSI','FSI', inDegree= min(params['inDegFSIFSI'],nbSim['FSI']) if ('inDegFSIFSI' in params) else min(25,nbSim['FSI']),gain=G['FSI'])

  print '* STN Inputs'
  connect('ex','PTN','STN', inDegree= min(params['inDegPTNSTN'],nbSim['PTN']) if ('inDegPTNSTN' in params) else min(25,nbSim['PTN']), gain=G['STN'])
  connect('ex','CMPf','STN',inDegree= min(params['inDegCMPfSTN'],nbSim['CMPf']) if ('inDegCMPfSTN' in params) else min(84,nbSim['CMPf']), gain=G['STN'])
  connect('in','GPe','STN', inDegree= min(params['inDegGPeSTN'],nbSim['GPe']) if ('inDegGPeSTN' in params) else min(30,nbSim['GPe']), gain=G['STN'])

  print '* GPe Inputs'
  if antagInjectionSite == 'GPe':
    if   antag == 'AMPA':
      connect('NMDA','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']) if ('inDegCMPfGPe' in params) else min(32/2,nbSim['CMPf']), gain=G['GPe'])
      connect('NMDA','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']) if ('inDegSTNGPe' in params) else min(107/2,nbSim['STN']), gain=G['GPe'])
      connect('in','MSN','GPe', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']) if ('inDegMSNGPe' in params) else min(14723/2,nbSim['MSN']), gain=G['GPe'])
      connect('in','GPe','GPe', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']) if ('inDegGPeGPe' in params) else min(32,nbSim['GPe']), gain=G['GPe'])
    elif antag == 'NMDA':
      connect('AMPA','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']) if ('inDegCMPfGPe' in params) else min(32/2,nbSim['CMPf']), gain=G['GPe'])
      connect('AMPA','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']) if ('inDegSTNGPe' in params) else min(107/2,nbSim['STN']), gain=G['GPe'])
      connect('in','MSN','GPe', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']) if ('inDegMSNGPe' in params) else min(14723/2,nbSim['MSN']), gain=G['GPe'])
      connect('in','GPe','GPe', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']) if ('inDegGPeGPe' in params) else min(32,nbSim['GPe']), gain=G['GPe'])
    elif antag == 'AMPA+GABAA':
      connect('NMDA','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']) if ('inDegCMPfGPe' in params) else min(32/2,nbSim['CMPf']), gain=G['GPe'])
      connect('NMDA','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']) if ('inDegSTNGPe' in params) else min(107/2,nbSim['STN']), gain=G['GPe'])
    elif antag == 'GABAA':
      connect('ex','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']) if ('inDegCMPfGPe' in params) else min(32/2,nbSim['CMPf']), gain=G['GPe'])
      connect('ex','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']) if ('inDegSTNGPe' in params) else min(107/2,nbSim['STN']), gain=G['GPe'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else: 
    connect('ex','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']) if ('inDegCMPfGPe' in params) else min(32/2,nbSim['CMPf']), gain=G['GPe'])
    connect('ex','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']) if ('inDegSTNGPe' in params) else min(107/2,nbSim['STN']), gain=G['GPe'])
    connect('in','MSN','GPe', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']) if ('inDegMSNGPe' in params) else min(14723/2,nbSim['MSN']), gain=G['GPe'])
    connect('in','GPe','GPe', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']) if ('inDegGPeGPe' in params) else min(32,nbSim['GPe']), gain=G['GPe'])

  print '* GPi Inputs'
  if antagInjectionSite =='GPi':
    if   antag == 'All':
      pass
    elif antag == 'NMDA':
      connect('in','MSN','GPi', inDegree= min(params['inDegMSNGPi'],nbSim['MSN']) if ('inDegMSNGPi' in params) else min(14723/2,nbSim['MSN']),gain=G['GPi'])
      connect('AMPA','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']) if ('inDegSTNGPi' in params) else min(107/2,nbSim['STN']),gain=G['GPi'])
      connect('in','GPe','GPi', inDegree= min(params['inDegGPeGPi'],nbSim['GPe']) if ('inDegGPeGPi' in params) else nbSim['GPe'],gain=G['GPi'])
      connect('AMPA','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']) if ('inDegCMPfGPi' in params) else min(30,nbSim['CMPf']),gain=G['GPi'])
    elif antag == 'NMDA+AMPA':
      connect('in','MSN','GPi', inDegree= min(params['inDegMSNGPi'],nbSim['MSN']) if ('inDegMSNGPi' in params) else min(14723/2,nbSim['MSN']),gain=G['GPi'])
      connect('in','GPe','GPi', inDegree= min(params['inDegGPeGPi'],nbSim['GPe']) if ('inDegGPeGPi' in params) else nbSim['GPe'],gain=G['GPi'])
    elif antag == 'AMPA':
      connect('in','MSN','GPi', inDegree= min(params['inDegMSNGPi'],nbSim['MSN']) if ('inDegMSNGPi' in params) else min(14723/2,nbSim['MSN']),gain=G['GPi'])
      connect('NMDA','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']) if ('inDegSTNGPi' in params) else min(107/2,nbSim['STN']),gain=G['GPi'])
      connect('in','GPe','GPi', inDegree= min(params['inDegGPeGPi'],nbSim['GPe']) if ('inDegGPeGPi' in params) else nbSim['GPe'],gain=G['GPi'])
      connect('NMDA','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']) if ('inDegCMPfGPi' in params) else min(30,nbSim['CMPf']),gain=G['GPi'])
    elif antag == 'GABAA':
      connect('ex','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']) if ('inDegSTNGPi' in params) else min(107/2,nbSim['STN']),gain=G['GPi'])
      connect('ex','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']) if ('inDegCMPfGPi' in params) else min(30,nbSim['CMPf']),gain=G['GPi'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connect('in','MSN','GPi', inDegree= min(params['inDegMSNGPi'],nbSim['MSN']) if ('inDegMSNGPi' in params) else min(14723/2,nbSim['MSN']),gain=G['GPi'])
    connect('ex','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']) if ('inDegSTNGPi' in params) else min(107/2,nbSim['STN']),gain=G['GPi'])
    connect('in','GPe','GPi', inDegree= min(params['inDegGPeGPi'],nbSim['GPe']) if ('inDegGPeGPi' in params) else nbSim['GPe'],gain=G['GPi'])
    connect('ex','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']) if ('inDegCMPfGPi' in params) else min(30,nbSim['CMPf']),gain=G['GPi'])

  #-------------------------
  # measures
  #-------------------------
  spkDetect={} # spike detectors used to record the experiment
  expeRate={}

  execTime = time.localtime()
  timeStr = str(execTime[0])+'_'+str(execTime[1])+'_'+str(execTime[2])+'_'+str(execTime[3])+':'+str(execTime[4])+':'+str(execTime[5])
  antagStr = ''
  if antagInjectionSite != 'none':
    antagStr = antagInjectionSite+'_'+antag+'_'

  for N in NUCLEI:
    # 500ms offset period for network stabilization
    spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": timeStr+'_'+antagStr+N, "to_file": True, 'start':500.,'stop':500.+simDuration})
    nest.Connect(Pop[N], spkDetect[N])

  #-------------------------
  # Simulation
  #-------------------------
  nest.Simulate(simDuration+500.)

  score = 0

  if antagInjectionSite == 'none':
    for N in NUCLEI:
      strTestPassed = 'NO!'
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
      if expeRate[N] <= FRRNormal[N][1] and expeRate[N] >= FRRNormal[N][0]:
        # if the measured rate is within acceptable values
        strTestPassed = 'OK'
        score += 1
      print '*',N,'- Rate:',expeRate[N],'Hz -> '+strTestPassed
  else:
    for N in NUCLEI:
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
      if N == antagInjectionSite:
        strTestPassed = 'NO!'
        if expeRate[N] <= FRRAnt[N][antag][1] and expeRate[N] >= FRRAnt[N][antag][0]:
          # if the measured rate is within acceptable values
          strTestPassed = 'OK'
          score += 1
          print '*',N,'with',antag,'antagonist(s):', expeRate[N],'Hz ->',strTestPassed
      else:
        print '*',N,'- Rate:',expeRate[N],'Hz'

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
def main(showRasters=True,params={}):
  score = np.zeros((2))
  score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='none',antag='')

  for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
    score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPe',antag=a)

  for a in ['All','AMPA','AMPA+NMDA','NMDA','GABAA']:
    score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPi',antag=a)

  # GPe inactivations:
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPe',antag='AMPA')
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPe',antag='AMPA+GABAA')
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPe',antag='NMDA')
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPe',antag='GABAA')
  # GPi inactivations:
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPi',antag='All')
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPi',antag='NMDA')
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPi',antag='NMDA+AMPA')
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPi',antag='AMPA')
  #score += checkAvgFR(showRasters=True,params=params,antagInjectionSite='GPi',antag='GABAA')


  #-------------------------
  print "******************"
  print "* Score:",score[0],'/',score[1]
  print "******************"

#---------------------------
if __name__ == '__main__':
  params = {'GMSN':4.37,
            'GFSI':1.3,
            'GSTN':1.35,
            'GGPe':1.3,
            'GGPi':1.3,
            'inDegFSIMSN':30, # according to Humphries et al. 2010, 30-150 FSIs->MSN
            'inDegMSNMSN':70, # according to Koos et al. 2004, cited by Humphries et al., 2010, on avg 3 synpase per MSN-MSN connection
            'inDegFSIFSI':15, # according to Humphries et al., 2010, 13-63 FSIs->FSI
            'inDegGPeGPi':23
            } 
  main(showRasters=True,params=params)
