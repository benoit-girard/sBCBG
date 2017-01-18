from LGneurons import *
import nest.raster_plot

testedNucleus = 'MSN'
simDuration = 5000. # ms
nest.SetKernelStatus({"overwrite_files":True})

if testedNucleus == 'STN':
#==========================
  print 'Creating neurons\n================'
  # creation of STN neurons
  #-------------------------
  nbSim['STN']=30.
  create('STN')

  # PTN
  #-------------------------
  nbSim['PTN'] = 150
  create('PTN', fake=True)
  #print '* PTN:',nbSim['PTN'],'Poisson generators with avg rate:',rate['PTN']
  #Pop['PTN']  = nest.Create('poisson_generator',int(nbSim['PTN']))
  #nest.SetStatus(Pop['PTN'],{'rate':rate['PTN']})

  # CMPf
  #-------------------------
  nbSim['CMPf']=33
  create('CMPf', fake=True)
  #print '* CMPf:',nbSim['CMPf'],'Poisson generators with avg rate:',rate['CMPf']
  #Pop['CMPf'] = nest.Create('poisson_generator',int(nbSim['CMPf']))
  #nest.SetStatus(Pop['CMPf'],{'rate': rate['CMPf']})

  # Fake GPe
  #-------------------------
  nbSim['GPe'] = 96
  create('GPe', fake=True)
  #print '* GPe:',nbSim['GPe'],'Poisson generators with avg rate:',rate['GPe']
  #Pop['GPe'] = nest.Create('poisson_generator',int(nbSim['GPe']))
  #nest.SetStatus(Pop['GPe'],{'rate':rate['GPe']})


  # connection of populations
  #-------------------------
  print 'Connecting neurons\n================'
  G = 1.35
  connect('ex','PTN','STN', inDegree=25, gain=G)
  connect('ex','CMPf','STN', inDegree=min(84,nbSim['CMPf']), gain=G)
  connect('in','GPe','STN', inDegree= 30, gain=G)

elif testedNucleus == 'GPe':
#==========================
  print 'Creating neurons\n================'
  # creation of GPe neurons
  #-------------------------
  nbSim['GPe']=30.
  create('GPe')
  nest.SetStatus(Pop['GPe'],{"I_e":13.}) # external current, necessary for tau_m=14 (not 20)

  # CMPf
  #-------------------------
  nbSim['CMPf']=int(nbSim['GPe']*0.34) # respects the approximate proportions
  create('CMPf', fake=True)

  # STN
  #-------------------------
  nbSim['STN']=int(nbSim['GPe']*0.3)
  create('STN', fake=True)

  # MSN
  #-------------------------
  nbSim['MSN']=200
  create('MSN', fake=True)

  print nbSim

  # connection of populations
  #---------------------------
  print 'Connecting neurons\n================'

  G = 1.3
  connect('ex','CMPf','GPe', inDegree=min(32/2,nbSim['CMPf']), gain=G)
  connect('ex','STN','GPe', inDegree=min(107/2,nbSim['STN']), gain=G)
  connect('in','MSN','GPe', inDegree=min(14723/2,nbSim['MSN']), gain=G)
  connect('in','GPe','GPe', inDegree=min(32,nbSim['GPe']), gain=G)

elif testedNucleus == 'MSN':
#==========================
  print 'Creating neurons\n================'

  # creation of MSN neurons
  #-------------------------
  nbSim['MSN']=50.
  create('MSN')

  # CSN
  #-------------------------
  nbSim['CSN'] = 200.
  create('CSN',fake=True)

  # PTN
  #-------------------------
  nbSim['PTN'] = 10.
  create('PTN',fake=True)

  # FSI
  #-------------------------
  nbSim['FSI'] = 1.
  create('FSI',fake=True)

  # CMPf
  #-------------------------
  nbSim['CMPf'] = 1.
  create('CMPf',fake=True)

  # connection of populations
  #---------------------------
  print 'Connecting neurons\n================'
  G = 3.9
  connect('ex','CSN','MSN',inDegree=100,gain=G)
  connect('ex','PTN','MSN',inDegree=1,gain=G)
  connect('ex','CMPf','MSN',inDegree=1,gain=G)
  connect('in','FSI','MSN',inDegree=1,gain=G)

elif testedNucleus == 'FSI':
#===========================
  print 'Creating neurons\n================'

  nbSim['FSI']=50.
  create('FSI')
  nbSim['STN']=7.
  create('STN',fake=True)
  nbSim['GPe']=25.
  create('GPe',fake=True)
  nbSim['CSN']=200
  create('CSN',fake=True)
  nbSim['PTN']=10.
  create('PTN',fake=True)
  nbSim['CMPf']=8.
  create('CMPf',fake=True)

  print 'Connecting neurons\n================'
  G = 1.1
  connect('ex','CSN','FSI',inDegree=50,gain=G)
  connect('ex','PTN','FSI',inDegree=1,gain=G)
  connect('ex','STN','FSI',inDegree=7,gain=G)
  connect('in','GPe','FSI',inDegree=25,gain=G)
  connect('ex','CMPf','FSI',inDegree=8,gain=G)

elif testedNucleus == 'GPi':
#===========================
  print 'Creating neurons\n================'

  nbSim['GPi']=50.
  create('GPi')
  # add a constant input current to the GPi so that event without excitatory inputs, we have activity
  nest.SetStatus(Pop['GPi'],{"I_e":9.})
  nbSim['MSN']=1000.
  create('MSN',fake=True)
  nbSim['STN']=30.
  create('STN',fake=True)
  nbSim['GPe']=88.
  create('GPe',fake=True)
  nbSim['CMPf']=30.
  create('CMPf',fake=True)

  print 'Connecting neurons\n================'
  G = 1.3
  connect('ex','STN','GPi',inDegree=30,gain=G)
  connect('ex','CMPf','GPi',inDegree=30,gain=G)
  connect('in','MSN','GPi',inDegree=50,gain=G)
  connect('in','GPe','GPi',inDegree=10,gain=G)


#-------------------------
# measures
#-------------------------

#mSTN = nest.Create("multimeter")
#nest.SetStatus(mSTN, {"withtime":True, "record_from":["V_m","currents"]})
#nest.Connect(mSTN, Pop['STN'])
  
spkDetect = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label":testedNucleus, "to_file": True})
nest.Connect(Pop[testedNucleus], spkDetect)

#-------------------------
# disconnection tests
#-------------------------

#GPeInConn = nest.GetConnections(target=Pop['GPe'])
#GPeInConnSTN = nest.GetConnections(target=Pop['GPe'],source=Pop['STN'])
#GPeInConnCMPf = nest.GetConnections(target=Pop['GPe'],source=Pop['CMPf'])
#print "GPe Exitatory inputs disconnected"
#print "NB CONNECT",len(GPeInConn),len(GPeInConnSTN)
#print GPeInConnSTN
#nest.SetStatus(GPeInConnSTN, {'weight':0.0}) # how to be specific to AMPA or NMDA ?
#nest.SetStatus(GPeInConnCMPf, {'weight':0.0}) # how to be specific to AMPA, NMDA or GABA inputs ?                                                                
#-------------------------
# Simulation
#-------------------------
#nest.ResetKernel()
#nest.SetKernelStatus({'local_num_threads': 2})
nest.Simulate(simDuration)
  
print '\n Spike Detector n_events',nest.GetStatus(spkDetect, 'n_events')[0]
expeRate = nest.GetStatus(spkDetect, 'n_events')[0] / float(nbSim[testedNucleus]*simDuration)
print '\n ',testedNucleus,'Rate:',expeRate*1000,'Hz'

# Displays
#-------------------------

#dmSTN = nest.GetStatus(mSTN)[0]
#VmSTN = dmSTN["events"]["V_m"]
#ImSTN = dmSTN["events"]["currents"]
#tSTN = dmSTN["events"]["times"]
'''  
dSD = nest.GetStatus(spkDetect,keys="events")[0]
evs = dSD["senders"]
ts = dSD["times"]
  
pylab.figure(testedNucleus+' spikes')
pylab.plot(ts, evs, ".")

pylab.show()
'''

nest.raster_plot.from_device(spkDetect,hist=True,title=testedNucleus)
nest.raster_plot.show()
