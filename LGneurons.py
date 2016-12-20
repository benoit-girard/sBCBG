# -*- coding: utf-8 -*-
import pylab
import nest
import numpy
import csv
from math import sqrt, cosh, exp, pi


#-------------------------------------------------------------------------------
# Establishes a connexion between two populations, following the results of LG14
# type : a string 'ex' or 'in', defining whether it is excitatory or inhibitory
# nameTgt, nameSrc : strings naming the populations, as defined in NUCLEI list
# srcPop : population emiting the connexions (nest node ID)
# tgtPop : population receiving the connexions (nest node ID)
#-------------------------------------------------------------------------------
def connect(type,nameSrc,nameTgt):
  print "*connect() DEBUG: connecting ",nameSrc,"to",nameTgt,"with",type,"connection"
  # process receptor types
  if type == 'ex':
    lRecType = ['AMPA','NMDA']
  elif type == 'in':
    lRecType = ['GABA']
  else:
    print "Undefined connexion type:",type

  # compute connexion strength
  # nu is the average total synaptic inputs a neuron of tgt receives from different neurons of src
  if nameSrc=='CSN' or nameSrc=='PTN':
    nu = alpha[nameSrc+'->'+nameTgt]
    print '\t',nameSrc+' -> '+nameTgt+': unknown number of different input neurons'
    print '\t',str(nbSim[nameSrc]),"neurons from",nameSrc,"will provide inputs"
  else:
    nu = neuronCounts[nameSrc] / neuronCounts[nameTgt] * P[nameSrc+'->'+nameTgt] * alpha[nameSrc+'->'+nameTgt]
    print '\t',nameSrc+' -> '+nameTgt+': avg number of different input neurons: '+str(neuronCounts[nameSrc] / neuronCounts[nameTgt] * P[nameSrc+'->'+nameTgt])
    print '\tcompare with the effective nb of neurons simulated '+str(nbSim[nameSrc]/float(nbSim[nameTgt])*P[nameSrc+'->'+nameTgt])

  # attenuation due to the distance from the receptors to the soma of tgt:
  attenuation = cosh(LX[nameTgt]*(1-p[nameSrc+'->'+nameTgt])) / cosh(LX[nameTgt])

  # final weight computation, modulated by the type of receptor
  # construction of the NEST connection
  for r in lRecType:
    print '\t',nu,nbSim[nameSrc],attenuation,r,recType[r]-1,wPSP # verifier /nbsim(nameTGT)
    w = nu / float(nbSim[nameSrc]) * attenuation * wPSP[recType[r]-1]
    nest.Connect(Pop[nameSrc],Pop[nameTgt],syn_spec={'receptor_type':recType[r],'weight':w})



NbNeurons = 1
dt = 0.01
simDuration = 1000. # in ms

# imported from Chadoeuf "connexweights"
# All the parameters needed to replicate Lienard model
#
#-------------------------

# Load the file with the Lienard solutions:
# solutions = csv.DictReader(open("solutions_simple.csv"))

# fixed parameters
A_GABA=-0.25 # mV
A_AMPA= 1.
A_NMDA= 0.025
D_GABA=5./exp(1)   # ms ; /e because Dn is peak half-time in LG14, while it is supposed to be tau_peak in NEST
D_AMPA=5./exp(1)
D_NMDA=100./exp(1)
Ri=200.E-2   # Ohms.m
Rm=20000.E-4 # Ohms.m^2

NUCLEI=['MSN','FSI','STN','GPe','GPi']

# Number of neurons in the real macaque brain
neuronCounts={'MSN': 26448.0E3,
              'FSI':   532.0E3,
              'STN':    77.0E3,
              'GPe':   251.0E3,
              'GPi':   143.0E3,
              'CMPf':   86.0E3
             }

# Number of neurons that will be simulated
nbSim = {'MSN': 0.,
         'FSI': 0.,
         'STN': 1.,
         'GPe': 0.,
         'GPi': 0.,
         'CMPf':0.,
         'CSN': 0.,
         'PTN': 0.,
        }

# P(X->Y): probability that a given neuron from X projects to at least neuron of Y
P = {'MSN->GPe': 0.82,
     'MSN->GPi': 1.,
     'FSI->MSN': 1.,
     'FSI->FSI': 1.,
     'STN->GPe': 0.83,
     'STN->GPi': 0.72,
     'STN->MSN': 0.17,
     'STN->FSI': 0.17,
     'GPe->STN': 1.,
     'GPe->GPe': 0.84,
     'GPe->GPi': 0.84,
     'GPe->MSN': 0.16,
     'GPe->FSI': 0.16,
     'CSN->MSN': 1.,
     'CSN->FSI': 1.,
     'PTN->MSN': 1.,
     'PTN->FSI': 1.,
     'PTN->STN': 1.,
     'CMPf->STN': 1.,
     'CMPf->MSN': 1.,
     'CMPf->FSI': 1.,
     'CMPf->GPe': 1.,
     'CMPf->GPi': 1.
    }

# alpha X->Y: average number of synaptic contacts made by one neuron of X to one neuron of Y, when there is a connexion
# for the moment set from one specific parameterization, should be read from Jean's solution file 
alpha = {'MSN->GPe':   171,
         'MSN->GPi':   210,
         'FSI->MSN':  4362,
         'FSI->FSI':   116,
         'STN->GPe':   428,
         'STN->GPi':   233,
         'STN->MSN':     0,
         'STN->FSI':    91,
         'GPe->STN':    19,
         'GPe->GPe':    38,
         'GPe->GPi':    16,
         'GPe->MSN':     0,
         'GPe->FSI':   353,
         'CSN->MSN':   342, # here, represents directly \nu
         'CSN->FSI':   250, # here, represents directly \nu
         'PTN->MSN':   342, # here, represents directly \nu ; TO BE CHECKED
         'PTN->FSI':   250, # here, represents directly \nu ; TO BE CHECKED
         'PTN->STN':   259, # here, represents directly \nu
         'CMPf->MSN': 4965,
         'CMPf->FSI': 1053,
         'CMPf->STN':   76,
         'CMPf->GPe':   79,
         'CMPf->GPi':  131
        }

# p(X->Y): probability that a given neuron from X projects to at least neuron of Y
# Warning: p is not P!
p = {'MSN->GPe':  0.48,
     'MSN->GPi':  0.59,
     'FSI->MSN':  0.19,
     'FSI->FSI':  0.16,
     'STN->GPe':  0.30,
     'STN->GPi':  0.59,
     'STN->MSN':  0.16,
     'STN->FSI':  0.41,
     'GPe->STN':  0.58,
     'GPe->GPe':  0.01,
     'GPe->GPi':  0.13,
     'GPe->MSN':  0.06,
     'GPe->FSI':  0.58,
     'CSN->MSN':  0.95,
     'CSN->FSI':  0.82,
     'PTN->MSN':  0.98,
     'PTN->FSI':  0.70,
     'PTN->STN':  0.97,
     'CMPf->STN': 0.46,
     'CMPf->MSN': 0.27,
     'CMPf->FSI': 0.06,
     'CMPf->GPe': 0.0,
     'CMPf->GPi': 0.48
    }

# electrotonic constant L computation:
dx={'MSN':1.E-6,'FSI':1.5E-6,'STN':1.5E-6,'GPe':1.7E-6,'GPi':1.2E-6}
lx={'MSN':619E-6,'FSI':961E-6,'STN':750E-6,'GPe':865E-6,'GPi':1132E-6}
LX={}
for n in NUCLEI:
    LX[n]=lx[n]*sqrt((4*Ri)/(dx[n]*Rm))

# P_X_Y: probability that a neuron from X projects to Y
P_MSN_GPe = 0.82
P_STN_GPe = 0.83
P_STN_GPi = 0.72
P_STN_MSN = 0.17
P_STN_FSI = 0.17
P_GPe_GPe = 0.84
P_GPe_GPi = 0.84
P_GPe_MSN = 0.16
P_GPe_FSI = 0.16
P_MSN_GPi = 0.82

# setting the 3 input ports for AMPA, NMDA and GABA receptor types
#-------------------------

nbPorts = 3
recType = {'AMPA':1,'NMDA':2,'GABA':3}
tau_syn = [D_AMPA, D_NMDA, D_GABA]
wPSP = [A_AMPA, A_NMDA, A_GABA]  # PSP amplitude (mV) ; A in LG14 notation

# parameterization of each neuronal type
#-------------------------

CommonParams = {'t_ref':         2.0,
                'V_m':           0.0,
                'V_th':         10.0, # dummy value to avoid NEST complaining about identical V_th and V_reset values
                'E_L':           0.0,
                'V_reset':       0.0,
                'I_e':           0.0,
                'V_min':       -20.0, # as in HSG06
                'tau_syn':   tau_syn
               }
nest.SetDefaults("iaf_psc_alpha_multisynapse", CommonParams)

Strparams = {'tau_m':        13.0, # according to SBE12
             'V_th':         30.0, # value of the LG14 example model, table 9
             'C_m':          13.0  # so that R_m=1, C_m=tau_m
            }

FSIparams = {'tau_m':         3.1, # from http://www.neuroelectro.org/article/75165/
             'V_th':         16.0, # value of the LG14 example model, table 9
             'C_m':           3.1  # so that R_m=1, C_m=tau_m
            }

STNparams = {'tau_m':         6.0, # as in HSG06 (but they model rats...)
             'V_th':         26.0, # value of the LG14 example model, table 9
             'C_m':           6.0  # so that R_m=1, C_m=tau_m
            }

GPeparams = {'tau_m':        20.0,
             'V_th':         11.0, # value of the LG14 example model, table 9
             'C_m':          20.0  # so that R_m=1, C_m=tau_m
            }

GPiparams = {'tau_m':        20.0,
             'V_th':          6.0, # value of the LG14 example model, table 9
             'C_m':          20.0, # so that R_m=1, C_m=tau_m
            }


# dictionary of the parameterizations of each neuronal type
#-------------------------
 
BGparams = {'Str':Strparams,
            'FSI':FSIparams,
            'STN':STNparams,
            'GPe':GPeparams,
            'GPi':GPiparams}


# Pop is the dictionary that will contain the Nest IDs of all populations in the model
#-------------------------
Pop = {}

# creation of STN neurons
#-------------------------
Pop['STN'] = nest.Create("iaf_psc_alpha_multisynapse",NbNeurons,params=BGparams['STN'])


#-------------------------
# creation of external inputs (ctx, CMPf)
#-------------------------

# CSN
#-------------------------
#Pop['CSN']  = nest.Create('poisson_generator')
#nest.SetStatus(Pop['CSN'],{'rate': 2.0})


# PTN
#-------------------------
nbSim['PTN'] = 5
Pop['PTN']  = nest.Create('poisson_generator',nbSim['PTN'])
nest.SetStatus(Pop['PTN'],{'rate':15.})

connect('ex','PTN','STN')

#r = recType['AMPA']
#print 'w PTN->STN, attenuation:', attenuate('STN', 0.97)
#w = 259. / nbPTN * attenuate('STN',0.97) * wPSP[r-1]# a terme, remplacer 259. par la valeur lue dans le fichier
#print 'w PTN->STN, weight:', w
#nest.Connect(PTN,STN,syn_spec={'receptor_type':r,'weight':w})

#r = recType['NMDA']
#w = 259. / nbPTN * attenuate('STN',0.97) * wPSP[r-1]# a terme, remplacer 259. par la valeur lue dans le fichier
#nest.Connect(PTN,STN,syn_spec={'receptor_type':r,'weight':w})

# CMPf
#-------------------------
nbSim['CMPf']=1.
Pop['CMPf'] = nest.Create('poisson_generator')
nest.SetStatus(Pop['CMPf'],{'rate': 4.})

# one STN neuron should receive N_CMPf/N_STN=1.1 inputs from the CMPf
# thus we simulate 1 CMPf neuron with a weight of 1.1

connect('ex','CMPf','STN')

#r = recType['AMPA']
#print 'w CMPf->STN, attenuation:', attenuate('STN', 0.46)
#w = N_CMPf / N_STN * attenuate('STN',0.46) * wPSP[r-1]
#print 'w CMPf->STN, weight:', w
#nest.Connect(CMPf,STN,syn_spec={'receptor_type':r,'weight':w})

#r = recType['NMDA']
#w = N_CMPf / N_STN * attenuate('STN',0.46) * wPSP[r-1]
#nest.Connect(CMPf,STN,syn_spec={'receptor_type':r,'weight':w})

# Fake GPe
#-------------------------
nbSim['GPe'] = int(neuronCounts['GPe']/neuronCounts['STN'])
Pop['GPe'] = nest.Create('poisson_generator',nbSim['GPe'])
nest.SetStatus(Pop['GPe'],{'rate':62.6})

connect('in','GPe','STN')
#r = recType['GABA']
#w = N_GPe/N_STN/nbFakeGPe * 19. * attenuate('STN',0.58) * wPSP[r-1]
#print 'w fakeGPe->STN, nb synapses per input neuron:',N_GPe/N_STN/nbFakeGPe * 19.
#print 'w fakeGPe->STN, attenuation:',attenuate('STN',0.58)
#print 'w fakeGPe->STN, weight:', w
#nest.Connect(fakeGPe,STN,syn_spec={'receptor_type':r,'weight':w})

#-------------------------
# Stimulator neuron & connection
#-------------------------

#Stimulator  = nest.Create("iaf_neuron")
#nest.SetStatus(Stimulator, {'tau_m':20.0})
#nest.SetStatus(Stimulator, {"I_e": 376.0}) #376.0

# print 'Stimulator param ',nest.GetStatus(Stimulator,'C_m')

# print 'STN param ',nest.GetStatus(STN,'recordables')

#r = recType['AMPA']
#nest.Connect(Stimulator,STN,syn_spec={'receptor_type':r,'weight':wPSP[r-1]})

# measures
#-------------------------

mSTN = nest.Create("multimeter")
nest.SetStatus(mSTN, {"withtime":True, "record_from":["V_m","currents"]})
nest.Connect(mSTN, Pop['STN'])

spkDetect = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
nest.Connect(Pop['STN'], spkDetect)


#mStim = nest.Create("multimeter")
#nest.SetStatus(mStim, {"withtime":True, "record_from":["V_m"]})
#nest.Connect(mStim, Stimulator)

# Simulation
#-------------------------
nest.Simulate(simDuration)

# Displays
#-------------------------

showStimV = False
showSynCurr = False

#dmStim = nest.GetStatus(mStim)[0]
#VmStim = dmStim["events"]["V_m"]
#tStim = dmStim["events"]["times"]

dmSTN = nest.GetStatus(mSTN)[0]
VmSTN = dmSTN["events"]["V_m"]
ImSTN = dmSTN["events"]["currents"]
tSTN = dmSTN["events"]["times"]

dSD = nest.GetStatus(spkDetect,keys="events")[0]
evs = dSD["senders"]
ts = dSD["times"]

pylab.figure('STN spikes')
pylab.plot(ts, evs, ".")

pylab.figure('STN Voltage')
if (showStimV):
  pylab.plot(tStim, VmStim)
pylab.plot(tSTN, VmSTN)

if (showSynCurr):
  pylab.figure("STN input PSPs")
  pylab.plot(tSTN, ImSTN)

pylab.show()
