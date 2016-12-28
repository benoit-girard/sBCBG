# -*- coding: utf-8 -*-
import pylab
import nest
import numpy as np
import numpy.random as rnd
import csv
from math import sqrt, cosh, exp, pi

rnd.seed(17)

#-------------------------------------------------------------------------------
# Creates a population of neurons
# name: string naming the population, as defined in NUCLEI list
# fake: if fake is True, the neurons will be replaced by Poisson generators, firing 
#       at the rate indicated in the "rate" dictionary
#-------------------------------------------------------------------------------
def create(name,fake=False):
  if nbSim[name] == 0:
    print 'ERROR: create(): nbSim['+name+'] = 0'
    exit()
  if fake:
    if rate[name] == 0:
      print 'ERROR: create(): rate['+name+'] = 0 Hz'
    print '* '+name+'(fake):',nbSim[name],'Poisson generators with avg rate:',rate[name]
    Pop[name]  = nest.Create('poisson_generator',int(nbSim[name]))
    nest.SetStatus(Pop[name],{'rate':rate[name]})
  else:
    print '* '+name+':',nbSim[name],'neurons with parameters:',BGparams[name]
    Pop[name] = nest.Create("iaf_psc_alpha_multisynapse",int(nbSim[name]),params=BGparams[name])

#-------------------------------------------------------------------------------
# Establishes a connexion between two populations, following the results of LG14
# type : a string 'ex' or 'in', defining whether it is excitatory or inhibitory
# nameTgt, nameSrc : strings naming the populations, as defined in NUCLEI list
# inDegree : number of neurons from Src project to a single Tgt neuron
# gain : allows to amplify the weight normally deduced from LG14
#-------------------------------------------------------------------------------
def connect(type,nameSrc,nameTgt,inDegree,delay=1.,gain=1.):
  print "* connecting ",nameSrc,"->",nameTgt,"with",type,"connection and",inDegree,"inputs"
  # process receptor types
  if type == 'ex':
    lRecType = ['AMPA','NMDA']
  elif type == 'AMPA':
    lRecType = ['AMPA']
  elif type == 'NMDA':
    lRecType = ['NMDA']
  elif type == 'in':
    lRecType = ['GABA']
  else:
    print "Undefined connexion type:",type

  # compute connexion strength
  # nu is the average total synaptic inputs a neuron of tgt receives from different neurons of src
  if nameSrc=='CSN' or nameSrc=='PTN':
    nu = alpha[nameSrc+'->'+nameTgt]
    print '\tnu',nu
    print '\t',nameSrc+' -> '+nameTgt+': unknown number of different input neurons'
    print '\t',str(inDegree),"neurons from",nameSrc,"will provide inputs"
  else:
    nu = neuronCounts[nameSrc] / neuronCounts[nameTgt] * P[nameSrc+'->'+nameTgt] * alpha[nameSrc+'->'+nameTgt]
    print '\tMaximal number of distinct input neurons (nu):',nu
    print '\tMinimal number of distinct input neurons     :',str(neuronCounts[nameSrc] / neuronCounts[nameTgt] * P[nameSrc+'->'+nameTgt])
    print '\tCompare with the effective chosen inDegree   :'+str(inDegree)

  # attenuation due to the distance from the receptors to the soma of tgt:
  attenuation = cosh(LX[nameTgt]*(1-p[nameSrc+'->'+nameTgt])) / cosh(LX[nameTgt])

  # To ensure that for excitatory connections, Tgt neurons receive AMPA and NMDA projections from the same Src neurons, 
  # we have to handle the "indegree" connectivity ourselves:
  for nTgt in range(int(nbSim[nameTgt])):
    inputTable = rnd.choice(int(nbSim[nameSrc]),size=int(inDegree),replace=False)

    for input in inputTable:
      # final weight computation, modulated by the type of receptor and
      # construction of the NEST connection
      for r in lRecType:
        #print '\t',nu,inDegree,attenuation,r,recType[r]-1,wPSP # verifier /nbsim(nameTGT)
        w = nu / float(inDegree) * attenuation * wPSP[recType[r]-1] * gain
        #print '\t',Pop[nameSrc], Pop[nameTgt]
        #print '\t',Pop[nameSrc][input], Pop[nameTgt][nTgt],w
        nest.Connect(pre=(Pop[nameSrc][input],), post=(Pop[nameTgt][nTgt],),syn_spec={'receptor_type':recType[r],'weight':w,'delay':delay})

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Establishes a connexion between two populations, following the results of LG14
# type : a string 'ex' or 'in', defining whether it is excitatory or inhibitory
# nameTgt, nameSrc : strings naming the populations, as defined in NUCLEI list
# outDegree : number of Tgt neurons targeted by one single Src neuron 
# 
#-------------------------------------------------------------------------------
'''
def connectOut(type,nameSrc,nameTgt,outDegree,delay=1.):
  print "* connecting ",nameSrc,"->",nameTgt,"with",type,"connection and",inDegree,"inputs"
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
    print '\tnu',nu
    print '\t',nameSrc+' -> '+nameTgt+': unknown number of different input neurons'
    print '\t',str(inDegree),"neurons from",nameSrc,"will provide inputs"
  else:
    nu = neuronCounts[nameSrc] / neuronCounts[nameTgt] * P[nameSrc+'->'+nameTgt] * alpha[nameSrc+'->'+nameTgt]
    print '\tnu',nu
    print '\t',nameSrc+' -> '+nameTgt+': avg number of different input neurons: '+str(neuronCounts[nameSrc] / neuronCounts[nameTgt] * P[nameSrc+'->'+nameTgt])
    print '\tcompare with the effective chosen inDegree '+str(inDegree)

  # attenuation due to the distance from the receptors to the soma of tgt:
  attenuation = cosh(LX[nameTgt]*(1-p[nameSrc+'->'+nameTgt])) / cosh(LX[nameTgt])

  # To ensure that for excitatory connections, Tgt neurons receive AMPA and NMDA projections from the same Src neurons, 
  # we have to handle the "indegree" connectivity ourselves:
  for nTgt in range(int(nbSim[nameTgt])):
    inputTable = rnd.choice(int(nbSim[nameSrc]),size=inDegree,replace=False)

    for input in inputTable:
      # final weight computation, modulated by the type of receptor and
      # construction of the NEST connection
      for r in lRecType:
        #print '\t',nu,inDegree,attenuation,r,recType[r]-1,wPSP # verifier /nbsim(nameTGT)
        w = nu / float(inDegree) * attenuation * wPSP[recType[r]-1]
        #print '\t',Pop[nameSrc], Pop[nameTgt]
        #print '\t',Pop[nameSrc][input], Pop[nameTgt][nTgt],w
        nest.Connect(pre=(Pop[nameSrc][input],), post=(Pop[nameTgt][nTgt],),syn_spec={'receptor_type':recType[r],'weight':w,'delay':delay})
'''
#-------------------------------------------------------------------------------

dt = 0.01 # ms
simDuration = 10000. # in ms

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
         'STN': 0.,
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
         'PTN->MSN':     5, # here, represents directly \nu ; TO BE CHECKED
         'PTN->FSI':     5, # here, represents directly \nu ; TO BE CHECKED
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

MSNparams = {'tau_m':        13.0, # according to SBE12
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

GPeparams = {'tau_m':        14.0, # 20 -> 14 based on Johnson & McIntyre 2008, JNphy)
             'V_th':         11.0, # value of the LG14 example model, table 9
             'C_m':          14.0  # so that R_m=1, C_m=tau_m
            }

GPiparams = {'tau_m':        14.0, # 20 -> 14 based on Johnson & McIntyre 2008, JNphy)
             'V_th':          6.0, # value of the LG14 example model, table 9
             'C_m':          14.0  # so that R_m=1, C_m=tau_m
            }


# dictionary of the parameterizations of each neuronal type
#-------------------------
 
BGparams = {'MSN':MSNparams,
            'FSI':FSIparams,
            'STN':STNparams,
            'GPe':GPeparams,
            'GPi':GPiparams}

Pop = {}

# the dictionary used to store the desired discharge rates of the various Poisson generators that will be used as external inputs
rate = {'CSN':   2.  ,
        'PTN':  15.  ,
        'CMPf':  4.  ,
        'MSN':   0.25, # MSN and the following will be used when the corresponding nucleus is not explicitely simulated
        'FSI':  16.6 ,
        'STN':  14.3 ,
        'GPe':  62.6 ,
        'GPi':  64.2
        } 

#---------------------------
def main():

  # Pop is the dictionary that will contain the Nest IDs of all populations in the model
  #-------------------------
  print 'Creating neurons'

  # creation of STN neurons
  #-------------------------
  nbSim['STN']=10.
  print '* STN:',nbSim['STN'],'neurons with parameters:',BGparams['STN']

  Pop['STN'] = nest.Create("iaf_psc_alpha_multisynapse",int(nbSim['STN']),params=BGparams['STN'])

  #-------------------------
  # creation of external inputs (ctx, CMPf)
  #-------------------------
  rate = {} # the dictionary used to store the desired discharge rates of the various Poisson generators that will be used as external inputs

  # CSN
  #-------------------------
  #Pop['CSN']  = nest.Create('poisson_generator')
  #nest.SetStatus(Pop['CSN'],{'rate': 2.0})


  # PTN
  #-------------------------
  nbSim['PTN'] = 5*nbSim['STN']
  rate['PTN'] =  15.
  print '* PTN:',nbSim['PTN'],'Poisson generators with avg rate:',rate['PTN']
  Pop['PTN']  = nest.Create('poisson_generator',int(nbSim['PTN']))
  nest.SetStatus(Pop['PTN'],{'rate':rate['PTN']})

  connect('ex','PTN','STN', inDegree=5)

  # CMPf
  #-------------------------
  nbSim['CMPf']=nbSim['STN']
  rate['CMPf']=  4.
  print '* CMPf:',nbSim['CMPf'],'Poisson generators with avg rate:',rate['CMPf']
  Pop['CMPf'] = nest.Create('poisson_generator',int(nbSim['CMPf']))
  nest.SetStatus(Pop['CMPf'],{'rate': rate['CMPf']})

  connect('ex','CMPf','STN', inDegree=1)

  # Fake GPe
  #-------------------------
  nbSim['GPe'] = int(neuronCounts['GPe']/neuronCounts['STN']) * nbSim['STN']
  rate['GPe']= 62.6
  print '* GPe:',nbSim['GPe'],'Poisson generators with avg rate:',rate['GPe']
  Pop['GPe'] = nest.Create('poisson_generator',int(nbSim['GPe']))
  nest.SetStatus(Pop['GPe'],{'rate':rate['GPe']})

  connect('in','GPe','STN', inDegree= int(neuronCounts['GPe']/neuronCounts['STN'])) 

  #-------------------------
  # measures
  #-------------------------

  mSTN = nest.Create("multimeter")
  nest.SetStatus(mSTN, {"withtime":True, "record_from":["V_m","currents"]})
  nest.Connect(mSTN, Pop['STN'])
  
  spkDetect = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
  nest.Connect(Pop['STN'], spkDetect)

  # Simulation
  #-------------------------
  nest.Simulate(simDuration)


  # Experimental estimation of the firing rate:
  print '\n Spike Detector n_events',nest.GetStatus(spkDetect, 'n_events')[0]
  expeRate = nest.GetStatus(spkDetect, 'n_events')[0] / float(nbSim['STN']*simDuration)
  print '\n Rate:',expeRate*1000,'Hz'


  # Displays
  #-------------------------

  showSynCurr = False
  showSTNVolt = False

  #print nest.GetStatus(mSTN)
  #print "============================="
  #print nest.GetStatus(mSTN)[0]

  dmSTN = nest.GetStatus(mSTN)[0]
  VmSTN = dmSTN["events"]["V_m"]
  ImSTN = dmSTN["events"]["currents"]
  tSTN = dmSTN["events"]["times"]
  
  dSD = nest.GetStatus(spkDetect,keys="events")[0]
  evs = dSD["senders"]
  ts = dSD["times"]

  pylab.figure('STN spikes')
  pylab.plot(ts, evs, ".")

  if showSTNVolt:
    pylab.figure('STN Voltage')
    pylab.plot(tSTN, VmSTN)

  if (showSynCurr):
    pylab.figure("STN input PSPs")
    pylab.plot(tSTN, ImSTN)

  pylab.show()

#---------------------------
if __name__ == '__main__':
  main()
