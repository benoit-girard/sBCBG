from LGneurons import *

simDuration = 5000. # ms

#==========================
print '\nCreating neurons\n================'

nbSim['MSN'] = 2644.
create('MSN')

nbSim['FSI'] = 53.
create('FSI')

nbSim['STN'] = 8.
create('STN')

nbSim['GPe'] = 25.
create('GPe')
nest.SetStatus(Pop['GPe'],{"I_e":13.})

nbSim['GPi'] = 14.
create('GPi')
nest.SetStatus(Pop['GPi'],{"I_e":9.})

nbSim['CSN'] = 3000
create('CSN', fake=True)

nbSim['PTN'] = 100.
create('PTN', fake=True)

nbSim['CMPf']=9.
create('CMPf', fake=True)

print "Number of simulated neurons:", nbSim

G = {'MSN': 3.9,
     'FSI': 1.1,
     'STN': 1.35,
     'GPe': 1.3,
     'GPi': 1.3,
    }
# connection of populations
#-------------------------
print '\nConnecting neurons\n================'
print '* MSN Inputs'
connect('ex','CSN','MSN', inDegree=100,gain=G['MSN'])
connect('ex','PTN','MSN', inDegree=1,gain=G['MSN'])
#connect('ex','CMPf','MSN',inDegree=nbSim['CMPf'],gain=G['MSN'])
connect('ex','CMPf','MSN',inDegree=1,gain=G['MSN'])
#connect('in','FSI','MSN', inDegree=nbSim['FSI'],gain=G['MSN'])
connect('in','FSI','MSN', inDegree=1,gain=G['MSN'])

print '* FSI Inputs'
connect('ex','CSN','FSI', inDegree=50,gain=G['FSI'])
connect('ex','PTN','FSI', inDegree=1,gain=G['FSI'])
#connect('ex','STN','FSI', inDegree=nbSim['STN'],gain=G['FSI'])
connect('ex','STN','FSI', inDegree=2,gain=G['FSI'])
connect('in','GPe','FSI', inDegree=nbSim['GPe'],gain=G['FSI'])
connect('ex','CMPf','FSI',inDegree=nbSim['CMPf'],gain=G['FSI'])

print '* STN Inputs'
connect('ex','PTN','STN', inDegree=min(25,nbSim['PTN']), gain=G['STN'])
connect('ex','CMPf','STN',inDegree=min(84,nbSim['CMPf']), gain=G['STN'])
connect('in','GPe','STN', inDegree=min(30,nbSim['GPe']), gain=G['STN'])

print '* GPe Inputs'
connect('ex','CMPf','GPe',inDegree=min(32/2,nbSim['CMPf']), gain=G['GPe'])
connect('ex','STN','GPe', inDegree=min(107/2,nbSim['STN']), gain=G['GPe'])
connect('in','MSN','GPe', inDegree=min(14723/2,nbSim['MSN']), gain=G['GPe'])
connect('in','GPe','GPe', inDegree=min(32,nbSim['GPe']), gain=G['GPe'])

print '* GPi Inputs'
connect('in','MSN','GPi', inDegree=min(14723/2,nbSim['MSN']),gain=G['GPi'])
connect('ex','STN','GPi', inDegree=min(107/2,nbSim['STN']),gain=G['GPi'])
connect('in','GPe','GPi', inDegree=nbSim['GPe'],gain=G['GPi'])
connect('ex','CMPf','GPi',inDegree=min(30,nbSim['CMPf']),gain=G['GPi'])

#-------------------------
# measures
#-------------------------

spkDetect={}
expeRate={}

for N in NUCLEI:
  spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
  nest.Connect(Pop[N], spkDetect[N])

# Simulation
#-------------------------
#nest.ResetKernel()
nest.Simulate(simDuration)

for N in NUCLEI:  
  # print '\n Spike Detector n_events',nest.GetStatus(spkDetect, 'n_events')[0]
  expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration)
  print '\n*',N,'- Rate:',expeRate[N]*1000,'Hz'

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
