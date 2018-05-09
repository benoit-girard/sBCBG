#!/usr/bin/python
# -*- coding: utf-8 -*-

interactive = False # avoid loading X dependent things
                   # set to False for simulations on Sango
storeGDF = True # unless overriden by run.py, keep spike rasters

import nstrand

import pandas as pd
import pylab

import nest
import numpy as np
import numpy.random as rnd
import csv
from math import sqrt, cosh, exp, pi

AMPASynapseCounter = 0 # counter variable for the fast connect

#-------------------------------------------------------------------------------
# Loads a given LG14 model parameterization
# ID must be in [0,14]
#-------------------------------------------------------------------------------
def loadLG14params(ID):
  # Load the file with the Lienard solutions:
  LG14SolutionsReader = csv.DictReader(open("solutions_simple_unique.csv"),delimiter=';')
  LG14Solutions = []
  for row in LG14SolutionsReader:
    LG14Solutions.append(row)

  print '### Parameterization #'+str(ID)+' from (Lienard & Girard, 2014) is used. ###'

  for k,v in alpha.iteritems():
    try:
      alpha[k] = round(float(LG14Solutions[ID]['ALPHA_'+k.replace('->','_')]),0)
    except:
      print('Could not find LG14 parameters for connection `'+k+'`, trying to run anyway.')

  for k,v in p.iteritems():
    try:
      p[k] = round(float(LG14Solutions[ID]['DIST_'+k.replace('->','_')]),2)
    except:
      print('Could not find LG14 parameters for connection `'+k+'`, trying to run anyway.')

  for k,v in BGparams.iteritems():
    try:
      BGparams[k]['V_th'] = round(float(LG14Solutions[ID]['THETA_'+k]),1)
    except:
      print('Could not find LG14 parameters for connection `'+k+'`, trying to run anyway.')


def loadThetaFromCustomparams(params):
  for k,v in BGparams.iteritems():
    try:
      newval = round(float(params['THETA_'+k]), 2)
      print("WARNING: overwriting LG14 value for theta in "+k+" from original value of "+str(BGparams[k]['V_th'])+" to new value: "+str(newval))
      BGparams[k]['V_th'] = newval # firing threshold
    except:
      print("INFO: keeping LG14 value for theta in "+k+" to its original value of "+str(BGparams[k]['V_th']))
      pass

#-------------------------------------------------------------------------------
# Changes the default of the iaf_psc_alpha_multisynapse neurons
# Very important because it defines the 3 types of receptors (AMPA, NMDA, GABA) that will be needed
# Has to be called after any KernelReset
#-------------------------------------------------------------------------------
def initNeurons():
  nest.SetDefaults("iaf_psc_alpha_multisynapse", CommonParams)

#-------------------------------------------------------------------------------
# Creates a population of neurons
# name: string naming the population, as defined in NUCLEI list
# fake: if fake is True, the neurons will be replaced by Poisson generators, firing
#       at the rate indicated in the "rate" dictionary
# parrot: do we use parrot neurons or not? If not, there will be no correlations in the inputs, and a waste of computation power...
#-------------------------------------------------------------------------------
def create(name,fake=False,parrot=True):
  if nbSim[name] == 0:
    print 'ERROR: create(): nbSim['+name+'] = 0'
    exit()
  if fake:
    if rate[name] == 0:
      print 'ERROR: create(): rate['+name+'] = 0 Hz'
    print '* '+name+'(fake):',nbSim[name],'Poisson generators with avg rate:',rate[name]
    if not parrot:
      print "/!\ /!\ /!\ /!\ \nWARNING: parrot neurons not used, no correlations in inputs\n"
      Pop[name]  = nest.Create('poisson_generator',int(nbSim[name]))
      nest.SetStatus(Pop[name],{'rate':rate[name]})
    else:
      Fake[name]  = nest.Create('poisson_generator',int(nbSim[name]))
      nest.SetStatus(Fake[name],{'rate':rate[name]})
      Pop[name]  = nest.Create('parrot_neuron',int(nbSim[name]))
      nest.Connect(pre=Fake[name],post=Pop[name],conn_spec={'rule':'one_to_one'})

  else:
    print '* '+name+':',nbSim[name],'neurons with parameters:',BGparams[name]
    Pop[name] = nest.Create("iaf_psc_alpha_multisynapse",int(nbSim[name]),params=BGparams[name])

#-------------------------------------------------------------------------------
# Creates a popolation of neurons subdivided in Multiple Channels
#
# name: string naming the population, as defined in NUCLEI list
# nbCh: integer stating the number of channels to be created
# fake: if fake is True, the neurons will be replaced by Poisson generators, firing
#       at the rate indicated in the "rate" dictionary
# parrot: do we use parrot neurons or not? If not, there will be no correlations in the inputs, and a waste of computation power...
#-------------------------------------------------------------------------------
def createMC(name,nbCh,fake=False,parrot=True):
  if nbSim[name] == 0:
    print 'ERROR: create(): nbSim['+name+'] = 0'
    exit()

  Pop[name]=[]

  if fake:
    Fake[name]=[]
    if rate[name] == 0:
      print 'ERROR: create(): rate['+name+'] = 0 Hz'
    print '* '+name+'(fake):',nbSim[name]*nbCh,'Poisson generators (divided in',nbCh,'channels) with avg rate:',rate[name]
    if not parrot:
      print "/!\ /!\ /!\ /!\ \nWARNING: parrot neurons not used, no correlations in inputs\n"
      for i in range(nbCh):
        Pop[name].append(nest.Create('poisson_generator',int(nbSim[name])))
        nest.SetStatus(Pop[name][i],{'rate':rate[name]})
    else:
      for i in range(nbCh):
        Fake[name].append(nest.Create('poisson_generator',int(nbSim[name])))
        nest.SetStatus(Fake[name][i],{'rate':rate[name]})
        Pop[name].append(nest.Create('parrot_neuron',int(nbSim[name])))
        nest.Connect(pre=Fake[name][i],post=Pop[name][i],conn_spec={'rule':'one_to_one'})

  else:
    print '* '+name+':',nbSim[name]*nbCh,'neurons (divided in',nbCh,'channels) with parameters:',BGparams[name]
    for i in range(nbCh):
      Pop[name].append(nest.Create("iaf_psc_alpha_multisynapse",int(nbSim[name]),params=BGparams[name]))


#------------------------------------------------------------------------------
# Routine to perform the fast connection using nest built-in `connect` function
# - `source` & `dest` are lists defining Nest IDs of source & target population
# - `synapse_label` is used to tag connections and be able to find them quickly
#   with function `mass_mirror`, that adds NMDA on top of AMPA connections
# - `inDegree`, `receptor_type`, `weight`, `delay` are Nest connection params
#------------------------------------------------------------------------------
def mass_connect(source, dest, synapse_label, inDegree, receptor_type, weight, delay, stochastic_delays=None, verbose=False):
  def printv(text):
    if verbose:
      print(text)

  # potential initialization of stochastic delays
  if 'stochastic_delays' != None and delay > 0:
    printv('Using stochastic delays in mass-connect')
    low = delay * 0.5
    high = delay * 1.5
    sigma = delay * stochastic_delays
    delay =  {'distribution': 'normal_clipped', 'low': low, 'high': high, 'mu': delay, 'sigma': sigma}

  # The first `fixed_indegree` connection ensures that all neurons in `dest`
  # are targeted by the same number of axons (an integer number)
  integer_inDegree = np.floor(inDegree)
  if integer_inDegree>0:
    printv('Adding '+str(int(integer_inDegree*len(dest)))+' connections with rule `fixed_indegree`\n')
    nest.Connect(source,
                 dest,
                 {'rule': 'fixed_indegree', 'indegree': int(integer_inDegree)},
                 {'model': 'static_synapse_lbl', 'synapse_label': synapse_label, 'receptor_type': receptor_type, 'weight': weight, 'delay':delay})

  # The second `fixed_total_number` connection distributes remaining axonal
  # contacts at random (i.e. the remaining fractional part after the first step)
  float_inDegree = inDegree - integer_inDegree
  remaining_connections = np.round(float_inDegree * len(dest))
  if remaining_connections > 0:
    printv('Adding '+str(remaining_connections)+' remaining connections with rule `fixed_total_number`\n')
    nest.Connect(source,
                 dest,
                 {'rule': 'fixed_total_number', 'N': int(remaining_connections)},
                 {'model': 'static_synapse_lbl', 'synapse_label': synapse_label, 'receptor_type': receptor_type, 'weight': weight, 'delay':delay})

#------------------------------------------------------------------------------
# Routine to duplicate a connection made with a specific receptor, with another
# receptor (typically to add NMDA connections to existing AMPA connections)
# - `source` & `synapse_label` should uniquely define the connections of
#   interest - typically, they are the same as in the call to `mass_connect`
# - `receptor_type`, `weight`, `delay` are Nest connection params
#------------------------------------------------------------------------------
def mass_mirror(source, synapse_label, receptor_type, weight, delay, stochastic_delays, verbose=False):
  def printv(text):
    if verbose:
      print(text)

  # find all AMPA connections for the given projection type
  printv('looking for AMPA connections to mirror with NMDA...\n')
  ampa_conns = nest.GetConnections(source=source, synapse_label=synapse_label)
  # in rare cases, there may be no connections, guard against that
  if ampa_conns:
    # extract just source and target GID lists, all other information is irrelevant here
    printv('found '+str(len(ampa_conns))+' AMPA connections\n')
    if 'stochastic_delays' != None and delay > 0:
      printv('Using stochastic delays in mass-miror')
      delay = np.array(nest.GetStatus(ampa_conns, keys=['delay'])).flatten()
    src, tgt, _, _, _ = zip(*ampa_conns)
    nest.Connect(src, tgt, 'one_to_one',
                 {'model': 'static_synapse_lbl',
                  'synapse_label': synapse_label, # tag with the same number (doesn't matter)
                  'receptor_type': receptor_type, 'weight': weight, 'delay':delay})

#-------------------------------------------------------------------------------
# Establishes a connexion between two populations, following the results of LG14
# type : a string 'ex' or 'in', defining whether it is excitatory or inhibitory
# nameTgt, nameSrc : strings naming the populations, as defined in NUCLEI list
# redundancy : value that characterizes the number of repeated axonal contacts from one neuron of Src to one neuron of Tgt (see RedundancyType for interpretation of this value)
# RedundancyType : string
#   if 'inDegreeAbs': `redundancy` is the number of neurons from Src that project to a single Tgt neuron
#   if 'outDegreeAbs': `redundancy` is number of axonal contacts between each neuron from Src onto a single Tgt neuron
#   if 'outDegreeCons': `redundancy` is a scaled proportion of axonal contacts between each neuron from Src onto a single Tgt neuron given arithmetical constraints, ranging from 0 (minimal number of contacts to achieve required axonal bouton counts) to 1 (maximal number of contacts with respect to population numbers)
# LCGDelays: shall we use the delays obtained by (Liénard, Cos, Girard, in prep) or not (default = True)
# gain : allows to amplify the weight normally deduced from LG14
#-------------------------------------------------------------------------------
def connect(type, nameSrc, nameTgt, redundancy, RedundancyType, LCGDelays=True, gain=1., stochastic_delays=None, verbose=False, projType=''):

  def printv(text):
    if verbose:
      print(text)

  printv("* connecting "+nameSrc+" -> "+nameTgt+" with "+type+" connection")

  if RedundancyType == 'inDegreeAbs':
    # inDegree is already provided in the right form
    inDegree = float(redundancy)
  elif RedundancyType == 'outDegreeAbs':
    #### fractional outDegree is expressed as a fraction of max axo-dendritic contacts
    inDegree = get_frac(1./redundancy, nameSrc, nameTgt, neuronCounts[nameSrc], neuronCounts[nameTgt], verbose=verbose)
  elif RedundancyType == 'outDegreeCons':
    #### fractional outDegree is expressed as a ratio of min/max axo-dendritic contacts
    inDegree = get_frac(redundancy, nameSrc, nameTgt, neuronCounts[nameSrc], neuronCounts[nameTgt], useMin=True, verbose=verbose)
  else:
    raise KeyError('`RedundancyType` should be one of `inDegreeAbs`, `outDegreeAbs`, or `outDegreeCons`.')

  # check if in degree acceptable (not larger than number of neurons in the source nucleus)
  if inDegree  > nbSim[nameSrc]:
    printv("/!\ WARNING: required 'in degree' ("+str(inDegree)+") larger than number of neurons in the source population ("+str(nbSim[nameSrc])+"), thus reduced to the latter value")
    inDegree = nbSim[nameSrc]

  if inDegree == 0.:
    printv("/!\ WARNING: non-existent connection strength, will skip")
    return

  global AMPASynapseCounter

  # process receptor types
  if type == 'ex':
    lRecType = ['AMPA','NMDA']
    AMPASynapseCounter = AMPASynapseCounter + 1
    lbl = AMPASynapseCounter # needs to add NMDA later
  elif type == 'AMPA':
    lRecType = ['AMPA']
    lbl = 0
  elif type == 'NMDA':
    lRecType = ['NMDA']
    lbl = 0
  elif type == 'in':
    lRecType = ['GABA']
    lbl = 0
  else:
    raise KeyError('Undefined connexion type: '+type)

  W = computeW(lRecType, nameSrc, nameTgt, inDegree, gain, verbose=False)

  printv("  W="+str(W)+" and inDegree="+str(inDegree))

  #if nameSrc+'->'+nameTgt in ConnectMap:
  #  loadConnectMap = True
  #else:
  #  loadConnectMap = False
  #  ConnectMap[nameSrc+'->'+nameTgt] = []

  # determine which transmission delay to use:
  if LCGDelays:
    delay= tau[nameSrc+'->'+nameTgt]
  else:
    delay= 1.

  mass_connect(Pop[nameSrc], Pop[nameTgt], lbl, inDegree, recType[lRecType[0]], W[lRecType[0]], delay, stochastic_delays = stochastic_delays)
  if type == 'ex':
    # mirror the AMPA connection with similarly connected NMDA connections
    mass_mirror(Pop[nameSrc], lbl, recType['NMDA'], W['NMDA'], delay, stochastic_delays = stochastic_delays)

  return W


#-------------------------------------------------------------------------------
# Establishes a connexion between two populations, following the results of LG14, in a MultiChannel context
# type : a string 'ex' or 'in', defining whether it is excitatory or inhibitory
# nameTgt, nameSrc : strings naming the populations, as defined in NUCLEI list
# projType : type of projections. For the moment: 'focused' (only channel-to-channel connection) and
#            'diffuse' (all-to-one with uniform distribution)
# redundancy, RedundancyType : contrains the inDegree - see function `connect` for details
# LCGDelays : shall we use the delays obtained by (Liénard, Cos, Girard, in prep) or not (default = True)
# gain : allows to amplify the weight normally deduced from LG14
# source_channels : By default with `source_channels=None`, the connection is implemented using all source channels
#                   Specify a custom list of channels to implement connections only from these channels
#                   For example, calling successively `connectMC(...,projType='focused',source_channels=[0])` and then `connectMC(...,projType='diffuse',source_channels=[1])` would implement first a focused projection using only source channel 0 and then a diffuse connection using only source channel 1:
#                   Src channels:   (0) (1)
#                                    | / |
#                   Tgt channels:   (0) (1)
#-------------------------------------------------------------------------------
def connectMC(type, nameSrc, nameTgt, projType, redundancy, RedundancyType, LCGDelays=True, gain=1., source_channels = None, stochastic_delays=None, verbose=False):

  def printv(text):
    if verbose:
      print(text)

  printv("* connecting "+nameSrc+" -> "+nameTgt+" with "+projType+" "+type+" connection")

  if source_channels == None:
    # if not specified, assume that the connection originates from all channels
    source_channels = range(len(Pop[nameSrc]))

  if RedundancyType == 'inDegreeAbs':
    # inDegree is already provided in the right form
    inDegree = float(redundancy)
  elif RedundancyType == 'outDegreeAbs':
    #### fractional outDegree is expressed as a fraction of max axo-dendritic contacts
    inDegree = get_frac(1./redundancy, nameSrc, nameTgt, neuronCounts[nameSrc], neuronCounts[nameTgt], verbose=verbose)
  elif RedundancyType == 'outDegreeCons':
    #### fractional outDegree is expressed as a ratio of min/max axo-dendritic contacts
    inDegree = get_frac(redundancy, nameSrc, nameTgt, neuronCounts[nameSrc], neuronCounts[nameTgt], useMin=True, verbose=verbose)
  else:
    raise KeyError('`RedundancyType` should be one of `inDegreeAbs`, `outDegreeAbs`, or `outDegreeCons`.')

  # check if in degree acceptable (not larger than number of neurons in the source nucleus)
  if projType == 'focused' and inDegree > nbSim[nameSrc]:
    printv("/!\ WARNING: required 'in degree' ("+str(inDegree)+") larger than number of neurons in individual source channels ("+str(nbSim[nameSrc])+"), thus reduced to the latter value")
    inDegree = nbSim[nameSrc]
  if projType == 'diffuse' and inDegree  > nbSim[nameSrc]*len(source_channels):
    printv("/!\ WARNING: required 'in degree' ("+str(inDegree)+") larger than number of neurons in the overall source population ("+str(nbSim[nameSrc]*len(source_channels))+"), thus reduced to the latter value")
    inDegree = nbSim[nameSrc]*len(source_channels)

  if inDegree == 0.:
    printv("/!\ WARNING: non-existent connection strength, will skip")
    return

  global AMPASynapseCounter

  inDegree = inDegree * (float(len(source_channels)) / float(len(Pop[nameSrc])))

  # prepare receptor type lists:
  if type == 'ex':
    lRecType = ['AMPA','NMDA']
    AMPASynapseCounter = AMPASynapseCounter + 1
    lbl = AMPASynapseCounter # needs to add NMDA later
  elif type == 'AMPA':
    lRecType = ['AMPA']
    lbl = 0
  elif type == 'NMDA':
    lRecType = ['NMDA']
    lbl = 0
  elif type == 'in':
    lRecType = ['GABA']
    lbl = 0
  else:
    raise KeyError('Undefined connexion type: '+type)

  # compute the global weight of the connection, for each receptor type:
  W = computeW(lRecType, nameSrc, nameTgt, inDegree, gain, verbose=False)

  printv("  W="+str(W)+" and inDegree="+str(inDegree))

  ## check whether a connection map has already been drawn or not:
  #if nameSrc+'->'+nameTgt in ConnectMap:
  #  #print "Using existing connection map"
  #  loadConnectMap = True
  #else:
  #  #print "Will create a connection map"
  #  loadConnectMap = False
  #  ConnectMap[nameSrc+'->'+nameTgt] = [[] for i in range(len(Pop[nameTgt]))]

  # determine which transmission delay to use:
  if LCGDelays:
    delay = tau[nameSrc+'->'+nameTgt]
  else:
    delay = 1.

  if projType == 'focused': # if projections focused, input come only from the same channel as tgtChannel
     for src_channel in source_channels: # for each relevant channel of the Source nucleus
       mass_connect(Pop[nameSrc][src_channel], Pop[nameTgt][src_channel-source_channels[0]], lbl, inDegree, recType[lRecType[0]], W[lRecType[0]], delay, stochastic_delays = stochastic_delays)
  elif projType == 'diffuse': # if projections diffused, input connections are shared among each possible input channel equally
    for src_channel in source_channels: # for each relevant channel of the Source nucleus
      for tgt_channel in range(len(Pop[nameTgt])): # for each channel of the Target nucleus
        mass_connect(Pop[nameSrc][src_channel], Pop[nameTgt][tgt_channel], lbl, inDegree/len(Pop[nameTgt]), recType[lRecType[0]], W[lRecType[0]], delay, stochastic_delays = stochastic_delays)

  if type == 'ex':
    # mirror the AMPA connection with similarly connected NMDA connections
    for src_channel in source_channels: # for each relevant channel of the Source nucleus
      mass_mirror(Pop[nameSrc][src_channel], lbl, recType['NMDA'], W['NMDA'], delay, stochastic_delays = stochastic_delays)

  return W

#-------------------------------------------------------------------------------
# returns the minimal & maximal numbers of distinct input neurons for one connection
#-------------------------------------------------------------------------------
def get_input_range(nameSrc, nameTgt, cntSrc, cntTgt, verbose=False):
  if nameSrc=='CSN' or nameSrc=='PTN':
    nu = alpha[nameSrc+'->'+nameTgt]
    nu0 = 0
    if verbose:
      print('\tMaximal number of distinct input neurons (nu): '+str(nu))
      print('\tMinimal number of distinct input neurons     : unknown (set to 0)')
  else:
    nu = cntSrc / float(cntTgt) * P[nameSrc+'->'+nameTgt] * alpha[nameSrc+'->'+nameTgt]
    nu0 = cntSrc / float(cntTgt) * P[nameSrc+'->'+nameTgt]
    if verbose:
      print('\tMaximal number of distinct input neurons (nu): '+str(nu))
      print('\tMinimal number of distinct input neurons     : '+str(nu0))
  return [nu0, nu]

#-------------------------------------------------------------------------------
# computes the inDegree as a fraction of maximal possible inDegree
# FractionalOutDegree: outDegree, expressed as a fraction
#-------------------------------------------------------------------------------
def get_frac(FractionalOutDegree, nameSrc, nameTgt, cntSrc, cntTgt, useMin=False, verbose=False):
  if useMin == False:
    # 'FractionalOutDegree' is taken to be relative to the maximal number of axo-dendritic contacts
    inDegree = get_input_range(nameSrc, nameTgt, cntSrc, cntTgt, verbose=verbose)[1] * FractionalOutDegree
  else:
    # 'FractionalOutDegree' is taken to be relative to the maximal number of axo-dendritic contacts and their minimal number
    r = get_input_range(nameSrc, nameTgt, cntSrc, cntTgt, verbose=verbose)
    inDegree = (r[1] - r[0]) * FractionalOutDegree + r[0]
  if verbose:
    print('\tConverting the fractional outDegree of '+nameSrc+' -> '+nameTgt+' from '+str(FractionalOutDegree)+' to inDegree neuron count: '+str(round(inDegree, 2))+' (relative to minimal value possible? '+str(useMin)+')')
  return inDegree

#-------------------------------------------------------------------------------
# computes the weight of a connection, based on LG14 parameters
#-------------------------------------------------------------------------------
def computeW(listRecType, nameSrc, nameTgt, inDegree, gain=1.,verbose=False):
  nu = get_input_range(nameSrc, nameTgt, neuronCounts[nameSrc], neuronCounts[nameTgt], verbose=verbose)[1]
  if verbose:
    print '\tCompare with the effective chosen inDegree   :',str(inDegree)

  # attenuation due to the distance from the receptors to the soma of tgt:
  attenuation = cosh(LX[nameTgt]*(1-p[nameSrc+'->'+nameTgt])) / cosh(LX[nameTgt])

  w={}
  for r in listRecType:
    w[r] = nu / float(inDegree) * attenuation * wPSP[recType[r]-1] * gain

  return w

#-------------------------------------------------------------------------------

dt = 0.01 # ms
simDuration = 10000. # in ms

# Acceptable firing rate ranges (FRR) in normal and deactivation experiments
# extracted from LG14 Table 5

FRRNormal = {'MSN': [0,1],
             'FSI': [7.8,14.0], # the refined constraint of 10.9 +/- 3.1 Hz was extracted from the following papers: Adler et al., 2016; Yamada et al., 2016 (summarizing date from three different experiments); and Marche and Apicella, 2017
             'STN': [15.2,22.8],
             'GPe': [55.7,74.5],
             'GPi': [59.1,79.5],
             }
FRRGPi = {'AMPA+NMDA+GABAA':[53.4,96.8],
          'NMDA':[27.2451,78.6255],
          'NMDA+AMPA':[6.811275,52.364583],
          'AMPA':[5.7327,66.0645],
          'GABAA':[44.1477,245.8935],
          }
FRRGPe = {'AMPA':[4.2889,58.7805],
          'AMPA+GABAA':[10.0017148,137.076126],
          'NMDA':[29.5767,61.1645],
          'GABAA':[74.8051,221.4885],
          }
FRRAnt = {'GPe':FRRGPe,'GPi':FRRGPi}

# imported from Chadoeuf "connexweights"
# All the parameters needed to replicate Lienard model
#
#-------------------------


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
# one hemisphere only, based on Hardman et al. 2002 paper, except for striatum & CM/Pf
neuronCounts={'MSN': 26448.0E3,
              'FSI':   532.0E3,
              'STN':    77.0E3,
              'GPe':   251.0E3,
              'GPi':   143.0E3,
              'CMPf':   86.0E3,
              'CSN': None, 'PTN': None # prevents key error
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
P = {'MSN->GPe': 1.,
     'MSN->GPi': 0.82,
     'MSN->MSN': 1.,
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
         'MSN->MSN':   210,
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
         'PTN->MSN':     5, # here, represents directly \nu
         'PTN->FSI':     5, # here, represents directly \nu
         'PTN->STN':   259, # here, represents directly \nu
         'CMPf->MSN': 4965,
         'CMPf->FSI': 1053,
         'CMPf->STN':   76,
         'CMPf->GPe':   79,
         'CMPf->GPi':  131
        }

# p(X->Y): relative distance on the dendrite from the soma, where neurons rom X projects to neurons of Y
# Warning: p is not P!
p = {'MSN->GPe':  0.48,
     'MSN->GPi':  0.59,
     'MSN->MSN':  0.77,
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

# tau: communication delays
tau = {'MSN->GPe':    7.,
       'MSN->GPi':   11.,
       'MSN->MSN':    1.,
       'FSI->MSN':    1.,
       'FSI->FSI':    1.,
       'STN->GPe':    3.,
       'STN->GPi':    3.,
       'STN->MSN':    3.,
       'STN->FSI':    3.,
       'GPe->STN':   10.,
       'GPe->GPe':    1.,
       'GPe->GPi':    3.,
       'GPe->MSN':    3.,
       'GPe->FSI':    3.,
       'CSN->MSN':    7.,
       'CSN->FSI':    7.,
       'PTN->MSN':    3.,
       'PTN->FSI':    3.,
       'PTN->STN':    3.,
       'CMPf->MSN':   7.,
       'CMPf->FSI':   7.,
       'CMPf->STN':   7.,#4
       'CMPf->GPe':   7.,#5
       'CMPf->GPi':   7.,#6
       }


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
initNeurons() # sets the default params of iaf_alpha_psc_mutisynapse neurons to CommonParams

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
Fake= {} # Fake contains the Poisson Generators, that will feed the parrot_neurons, stored in Pop
ConnectMap = {} # when connections are drawn, in "create()", they are stored here so as to be re-usable

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

  if interactive:
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
