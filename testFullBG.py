from LGneurons import *
import nest.raster_plot
import time

# params possible keys:
# - nb{MSN,FSI,STN,GPi,GPe,CSN,PTN,CMPf} : number of simulated neurons for each population
# - Ie{GPe,GPi} : constant input current to GPe and GPi
# - G{MSN,FSI,STN,GPi,GPe} : gain to be applied on LG14 input synaptic weights for each population

def main(showRasters=True,params={}):

  simDuration = 3000. # ms
  # nest.SetKernelStatus({"overwrite_files":True}) # Thanks to use of timestamps, file names should now 
                                                   #be different as long as they are not created during the same second

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
  print '* MSN Inputs'
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
  connect('ex','CMPf','GPe',inDegree= min(params['inDegCMPfGPe'],nbSim['CMPf']) if ('inDegCMPfGPe' in params) else min(32/2,nbSim['CMPf']), gain=G['GPe'])
  connect('ex','STN','GPe', inDegree= min(params['inDegSTNGPe'],nbSim['STN']) if ('inDegSTNGPe' in params) else min(107/2,nbSim['STN']), gain=G['GPe'])
  connect('in','MSN','GPe', inDegree= min(params['inDegMSNGPe'],nbSim['MSN']) if ('inDegMSNGPe' in params) else min(14723/2,nbSim['MSN']), gain=G['GPe'])
  connect('in','GPe','GPe', inDegree= min(params['inDegGPeGPe'],nbSim['GPe']) if ('inDegGPeGPe' in params) else min(32,nbSim['GPe']), gain=G['GPe'])

  print '* GPi Inputs'
  connect('in','MSN','GPi', inDegree= min(params['inDegMSNGPi'],nbSim['MSN']) if ('inDegMSNGPi' in params) else min(14723/2,nbSim['MSN']),gain=G['GPi'])
  connect('ex','STN','GPi', inDegree= min(params['inDegSTNGPi'],nbSim['STN']) if ('inDegSTNGPi' in params) else min(107/2,nbSim['STN']),gain=G['GPi'])
  connect('in','GPe','GPi', inDegree= min(params['inDegGPeGPi'],nbSim['GPe']) if ('inDegGPeGPi' in params) else nbSim['GPe'],gain=G['GPi'])
  connect('ex','CMPf','GPi',inDegree= min(params['inDegCMPfGPi'],nbSim['CMPf']) if ('inDegCMPfGPi' in params) else min(30,nbSim['CMPf']),gain=G['GPi'])

  #-------------------------
  # measures
  #-------------------------

  spkDetect={} # spike detector used to record the whole experiment
  expeRate={}

  execTime = time.localtime()
  timeStr = str(execTime[0])+'_'+str(execTime[1])+'_'+str(execTime[2])+'_'+str(execTime[3])+':'+str(execTime[4])+':'+str(execTime[5])

  for N in NUCLEI:
    spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": timeStr+'_'+N, "to_file": True})
    nest.Connect(Pop[N], spkDetect[N])

  #-------------------------
  # Simulation
  #-------------------------
  tOffset = 0
  #nest.ResetKernel()
  nest.Simulate(simDuration)

  score = 0

  for N in NUCLEI:
    strTestPassed = 'NO!'
    # print '\n Spike Detector n_events',nest.GetStatus(spkDetect, 'n_events')[0]
    expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
    # print N,':',expeRate[N],' in ',FRRNormal[N][0], FRRNormal[N][1],'?'
    if expeRate[N] <= FRRNormal[N][1] and expeRate[N] >= FRRNormal[N][0]:
      # if the measured rate is within acceptable values
      strTestPassed = 'OK'
      score += 1
    print '\n*',N,'- Rate:',expeRate[N],'Hz -> '+strTestPassed

  #-------------------------
  # Antagonist injection simulation in GPe
  #-------------------------
  print "\n* GPe DEACTIVATIONS:"
  print "   - GABA Antagonist"

  tOffset += simDuration + 500.

  # spikedetector for deactivated period
  spkDetGPeNoGABA = nest.Create('spike_detector', params={"withgid": True, "withtime": True, "start": tOffset, "stop": tOffset+simDuration})
  nest.Connect(Pop['GPe'],spkDetGPeNoGABA)

  # GABA deactivation
  GPeInConnGPe = nest.GetConnections(target=Pop['GPe'], source=Pop['GPe'])
  GPeInConnMSN = nest.GetConnections(target=Pop['GPe'], source=Pop['MSN'])
  nest.SetStatus(GPeInConnGPe, {'weight':0.0}) 
  nest.SetStatus(GPeInConnMSN, {'weight':0.0}) 

  nest.Simulate(simDuration+500.)

  rate = nest.GetStatus(spkDetGPeNoGABA, 'n_events')[0] / float(nbSim['GPe']*simDuration) * 1000
  strTestPassed = 'NO!'
  if rate <= FRRGPe['GABAA'][1] and rate >= FRRGPe['GABAA'][0]:
    strTestPassed = 'OK'
    score += 1    
  print "  * No GABA_A:", rate,'Hz ->',strTestPassed

  print "******************"
  print "* Score:",score
  print "******************"

  #for N in NUCLEI:
  #  expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration) * 1000
  #  print '\n*',N,'- Rate:',expeRate[N],'Hz'


  #-------------------------
  # Displays
  #-------------------------
  '''
  dSD = nest.GetStatus(spkDetect,keys="events")[0]
  evs = dSD["senders"]
  ts = dSD["times"]
  
  pylab.figure(testedNucleus+' spikes')
  pylab.plot(ts, evs, ".")

  pylab.show()
  '''
  if showRasters and interactive:
    for N in NUCLEI:
      nest.raster_plot.from_device(spkDetect[N],hist=True,title=N)
    nest.raster_plot.from_device(spkDetGPeNoGABA,hist=True,title='GPe antag GABA')

    nest.raster_plot.show()

#---------------------------
if __name__ == '__main__':
  params = {'GMSN':4.37,
            'GFSI':1.3,
            'GSTN':1.35,
            'GGPi':1.3,
            'inDegFSIMSN':30, # according to Humphries et al. 2010, 30-150 FSIs->MSN
            'inDegMSNMSN':70, # according to Koos et al. 2004, cited by Humphries et al., 2010, on avg 3 synpase per MSN-MSN connection
            'inDegFSIFSI':15, # according to Humphries et al., 2010, 13-63 FSIs->FSI
            'inDegGPeGPi':23
            } 
  main(showRasters=True,params=params)
