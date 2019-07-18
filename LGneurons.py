#!/usr/bin/python
# -*- coding: utf-8 -*-

interactive = False # avoid loading X dependent things
                   # set to False for simulations on Sango
storeGDF = True # unless overriden by run.py, keep spike rasters

import nstrand

import pandas as pd
from modelParams import *
import nest
import numpy as np
import numpy.random as rnd
import csv
from math import sqrt, cosh, exp, pi

AMPASynapseCounter = 0 # counter variable for the fast connect

import nest.topology as nesttopo

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

  print('### Parameterization #'+str(ID)+' from (Lienard & Girard, 2014) is used. ###')

  for k,v in alpha.items():
    try:
      #print k,v,round(float(LG14Solutions[ID]['ALPHA_'+k.replace('->','_')]),0)
      alpha[k] = round(float(LG14Solutions[ID]['ALPHA_'+k.replace('->','_')]),0)
    except:
      print('Could not find LG14 alpha parameters for connection `'+k+'`, trying to run anyway.')

  for k,v in p.items():
    try:
      #print 'dist:',k,v,round(float(LG14Solutions[ID]['DIST_'+k.replace('->','_')]),2)
      p[k] = round(float(LG14Solutions[ID]['DIST_'+k.replace('->','_')]),2)
    except:
      print('Could not find LG14 distance (p) parameters for connection `'+k+'`, trying to run anyway.')

  for k,v in BGparams.items():
    try:
      BGparams[k]['V_th'] = round(float(LG14Solutions[ID]['THETA_'+k]),1)
    except:
      print('Could not find LG14 theta parameters for connection `'+k+'`, trying to run anyway.')

#-------------------------------------------------------------------------------
# Overrides the neuron threshold parameters from LG14 with those defined in parameters
# Added by Jean during his stay at OIST, should not be used anymore.
#-------------------------------------------------------------------------------
#def loadThetaFromCustomparams(params):
#  for k,v in BGparams.items():
#    try:
#      newval = round(float(params['THETA_'+k]), 2)
#      print("WARNING: overwriting LG14 value for theta in "+k+" from original value of "+str(BGparams[k]['V_th'])+" to new value: "+str(newval))
#      BGparams[k]['V_th'] = newval # firing threshold
#    except:
#      print("INFO: keeping LG14 value for theta in "+k+" to its original value of "+str(BGparams[k]['V_th']))
#      pass


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
    print('ERROR: create(): nbSim['+name+'] = 0')
    exit()
  if fake:
    if rate[name] == 0:
      print ('ERROR: create(): rate['+name+'] = 0 Hz')
    print('* '+name+'(fake):',nbSim[name],'Poisson generators with avg rate:',rate[name])
    if not parrot:
      print("/!\ /!\ /!\ /!\ \nWARNING: parrot neurons not used, no correlations in inputs\n")
      Pop[name]  = nest.Create('poisson_generator',int(nbSim[name]))
      nest.SetStatus(Pop[name],{'rate':rate[name]})
    else:
      Fake[name]  = nest.Create('poisson_generator',int(nbSim[name]))
      nest.SetStatus(Fake[name],{'rate':rate[name]})
      Pop[name]  = nest.Create('parrot_neuron',int(nbSim[name]))
      nest.Connect(pre=Fake[name],post=Pop[name],conn_spec={'rule':'one_to_one'})

  else:
    print('* '+name+':',nbSim[name],'neurons with parameters:',BGparams[name])
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
    print('ERROR: create(): nbSim['+name+'] = 0')
    exit()

  Pop[name]=[]

  if fake:
    Fake[name]=[]
    if rate[name] == 0:
      print('ERROR: create(): rate['+name+'] = 0 Hz')
    print('* '+name+'(fake):',nbSim[name]*nbCh,'Poisson generators (divided in',nbCh,'channels) with avg rate:',rate[name])
    if not parrot:
      print ("/!\ /!\ /!\ /!\ \nWARNING: parrot neurons not used, no correlations in inputs\n")
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
    print ('* '+name+':',nbSim[name]*nbCh,'neurons (divided in',nbCh,'channels) with parameters:',BGparams[name])
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

  sigmaDependentInterval = True # Hugo's method

  # potential initialization of stochastic delays
  if stochastic_delays != None and delay > 0 and stochastic_delays > 0.:
    printv('Using stochastic delays in mass-connect')
    sigma = delay * stochastic_delays
    if sigmaDependentInterval:
      n = 2 # number of standard deviation to include in the distribution
      if stochastic_delays >= 1./n:
        print('Error : stochastic_delays >= 1/n and the distribution of delays therefore includes 0 which is not possible -> Jean\'s method is used')
        sigmaDependentInterval = False
      else:
        low = delay - n*sigma
        high = delay + n*sigma


    if not sigmaDependentInterval:
      low = .5*delay
      high = 1.5*delay

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
    if stochastic_delays != None and delay > 0:
      printv('Using stochastic delays in mass-mirror')
      delay = np.array(nest.GetStatus(ampa_conns, keys=['delay'])).flatten()
    src, tgt, _, _, _ = zip(*ampa_conns)
    nest.Connect(src, tgt, 'one_to_one',
                 {'model': 'static_synapse_lbl',
                  'synapse_label': synapse_label, # tag with the same number (doesn't matter)
                  'receptor_type': receptor_type, 'weight': weight, 'delay':delay})

#-------------------------------------------------------------------------------
# Creates a topological population of neurons subdivided in Multiple Channels
#
# name: string naming the population, as defined in NUCLEI list
# nbCh: integer stating the number of channels to be created
# c: distance to the center (small distance means more channels in competition)
# r: radius of each channel (leading to larger overlap and thus broader competition)
# fake: if fake is True, the neurons will be replaced by Poisson generators, firing
#       at the rate indicated in the "rate" dictionary
# parrot: do we use parrot neurons or not? If not, there will be no correlations in the inputs, and a waste of computation power...
#-------------------------------------------------------------------------------
def createTopoMC(name, nbCh, layout, c=0.3, r=0.25, fake=False, parrot=True):
  if nbSim[name] == 0:
    print('ERROR: create(): nbSim['+name+'] = 0')
    exit()

  Pop[name]=[]

  # helper function that gives the channel center
  def circular_center(nbCh, c, Ch=None):
    # equi-distant points on a circle
    if Ch == None:
      indices = np.arange(0, nbCh, dtype=float) + 0.5
    else:
      indices = np.array(Ch) + 0.5
    angles = (1. - indices/nbCh) * 2. * np.pi
    x, y = np.cos(angles)*c, np.sin(angles)*c
    ## re-project in [0,1]x[0,1]
    #x = (x + 1.) / 2.
    #y = (y + 1.) / 2.
    return {'x': x, 'y': y}

  def circular_positions(nbCh, c, r, sim_pts, Ch=None):
    if Ch == None:
      Ch = range(nbCh)
    center_xy = circular_center(nbCh, c, Ch=Ch)
    xSim = []
    ySim = []
    for i in range(len(Ch)):
      angleSim = np.random.uniform(0., 2.*np.pi, int(sim_pts))
      rSim = np.random.uniform(0., r, int(sim_pts))
      xSim = xSim + (np.cos(angleSim)*rSim + center_xy['x'][i]).tolist()
      ySim = ySim + (np.sin(angleSim)*rSim + center_xy['y'][i]).tolist()
    return (xSim, ySim)

  def grid_positions(nbCh, sim_pts):
    n = int(sim_pts*nbCh)
    n_squared = np.ceil(np.sqrt(n))
    coord = [[x/n_squared*2.-1., y/n_squared*2.-1.] for x in np.arange(0,n_squared, dtype=float) for y in np.arange(0,n_squared, dtype=float)]
    # too many points due to square root rounding? remove at random
    if len(coord) > n:
      coord = np.array(coord)[np.sort(np.random.choice(range(len(coord)), size=n, replace=False))].tolist()
    return ([coord[i][0] for i in range(len(coord))], [coord[i][1] for i in range(len(coord))])

  # compute the neuron coordinates
  if layout == 'circular':
    positions = circular_positions(nbCh, c, r, nbSim[name])
    edge_wrap = False
  elif layout == 'grid':
    positions = grid_positions(nbCh, nbSim[name])
    edge_wrap = True
  else:
    raise KeyError('`layout` must be `circular` or `grid`.')

  if fake:
    Fake[name]=[]
    if rate[name] == 0:
      print('ERROR: create(): rate['+name+'] = 0 Hz')
    print ('* '+name+'(fake):',nbSim[name]*nbCh,'Poisson generators (divided in',nbCh,'channels) with avg rate:',rate[name])
    if not parrot:
      print("/!\ /!\ /!\ /!\ \nWARNING: parrot neurons not used, no correlations in inputs\n")
      Topo[name] = nesttopo.CreateLayer({'positions': [[positions[0][i], positions[1][i]] for i in range(len(positions[0]))], 'elements': 'poisson_generator', 'extent': [2., 2.], 'center':[0., 0.], 'edge_wrap': edge_wrap})
      all_nodes = nest.GetNodes(Topo[name])
      for i in range(nbCh):
        Pop[name].append([all_nodes[0][j] for j in np.arange(int(i*nbSim[name]), int((i+1)*nbSim[name]))])
        nest.SetStatus(Pop[name][i],{'rate':rate[name]})
    else:
      Topo[name] = nesttopo.CreateLayer({'positions': [[positions[0][i], positions[1][i]] for i in range(len(positions[0]))], 'elements': 'parrot_neuron', 'extent': [2., 2.], 'center':[0., 0.], 'edge_wrap': edge_wrap})
      all_nodes = nest.GetNodes(Topo[name])
      for i in range(nbCh):
        Fake[name].append(nest.Create('poisson_generator',int(nbSim[name])))
        nest.SetStatus(Fake[name][i],{'rate':rate[name]})
      for i in range(nbCh):
        Pop[name].append([all_nodes[0][j] for j in np.arange(int(i*nbSim[name]), int((i+1)*nbSim[name]))])
        nest.Connect(pre=Fake[name][i],post=Pop[name][i],conn_spec={'rule':'one_to_one'})

  else:
    print('* '+name+':',nbSim[name]*nbCh,'neurons (divided in',nbCh,'channels) with parameters:',BGparams[name])
    nest.SetDefaults('iaf_psc_alpha_multisynapse', BGparams[name])
    Topo[name] = nesttopo.CreateLayer({'positions': [[positions[0][i], positions[1][i]] for i in range(len(positions[0]))], 'elements': 'iaf_psc_alpha_multisynapse', 'extent': [2., 2.], 'center':[0., 0.], 'edge_wrap': edge_wrap})
    all_nodes = nest.GetNodes(Topo[name])
    for i in range(nbCh):
      Pop[name].append([all_nodes[0][j] for j in np.arange(int(i*nbSim[name]), int((i+1)*nbSim[name]))])

  # writes the layout to a file
  dataPath='log/'
  topoPositions = 'ID, Ch, X, Y, Z\n'
  topoFile=open(dataPath+'topo_'+name+'.csv','w',1)
  ni = 0
  for i in range(nbCh):
    for j in range(int(nbSim[name])):
      topoPositions += str(Pop[name][i][j])+', '+str(i)+', '+str(positions[0][ni])+', '+str(positions[1][ni])+', '+str(0)+'\n'
      ni += 1
  topoFile.write(topoPositions)
  topoFile.close()

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

  sigmaDependentInterval = True # using Hugo's method

  # potential initialization of stochastic delays
  if stochastic_delays != None and delay > 0 and stochastic_delays > 0.:
    printv('Using stochastic delays in mass-connect')
    sigma = delay * stochastic_delays
    if sigmaDependentInterval:
      n = 2 # number of standard deviation to include in the distribution
      if stochastic_delays >= 1./n:
        print ('Error : stochastic_delays >= 1/n and the distribution of delays therefore includes 0 which is not possible -> Jean\'s method is used')
        sigmaDependentInterval = False
      else:
        low = delay - n*sigma
        high = delay + n*sigma
    else:
      low = .5*delay
      high = 1.5*delay
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
# Routine to perform the fast connection using nest built-in `connect` function
# And in the topological case
# - `sourceName` & `destName` are names of two different layers
# - `synapse_label` is used to tag connections and be able to find them quickly
#   with function `mass_mirror`, that adds NMDA on top of AMPA connections
# - `inDegree`, `receptor_type`, `weight`, `delay` are Nest connection params
# - `spread` is a parameter that affects the diffusion level of the connection
#------------------------------------------------------------------------------
def mass_connect_topo(sourceName, destName, synapse_label, inDegree, receptor_type, weight, delay, spread, stochastic_delays=None, verbose=False):
  def printv(text):
    if verbose:
      print(text)

  # potential initialization of stochastic delays
  if stochastic_delays != None and delay > 0:
    printv('Using stochastic delays in mass-connect')
    low = delay * 0.5
    high = delay * 1.5
    sigma = delay * stochastic_delays
    delay =  {'distribution': 'normal_clipped', 'low': low, 'high': high, 'mu': delay, 'sigma': sigma}

  ## creation of the synapse model
  #nest.CopyModel('static_synapse_lbl', 'mass_connected_'+sourceName+'_'+destName, {'synapse_label': synapse_label, 'receptor_type': receptor_type, 'weight': weight, 'delay':delay})

  ## creation of the topological connection dict
  #conndict = {'connection_type': 'convergent',
  #            'mask': {'circular': {'radius': spread}},
  #            'synapse_model': 'mass_connected_'+sourceName+'_'+destName
  #           }

  #nest.CopyModel('static_synapse_lbl', 'mass_connected_'+sourceName+'_'+destName, {'synapse_label': synapse_label, 'receptor_type': receptor_type, 'weight': weight, 'delay':delay})

  nest.SetDefaults('static_synapse_lbl', {'synapse_label': synapse_label, 'receptor_type': receptor_type})

  # creation of the topological connection dict
  conndict = {'connection_type': 'convergent',
              'mask': {'circular': {'radius': spread}},
              'synapse_model': 'static_synapse_lbl', 'weights': weight, 'delays':delay,
              'allow_oversized_mask': True, 'allow_multapses': True}

  # The first call ensures that all neurons in `destName`
  # have at least `int(inDegree)` incoming connections
  integer_inDegree = np.floor(inDegree)
  if integer_inDegree>0:
    printv('Adding '+str(int(integer_inDegree*len(Pop[destName])))+' connections with rule `fixed_indegree`\n')
    integer_conndict = conndict.copy()
    integer_conndict.update({'number_of_connections': int(integer_inDegree)})
    nesttopo.ConnectLayers(Topo[sourceName], Topo[destName], integer_conndict)

  # The second call distributes the approximate number of remaining axonal
  # contacts at random (i.e. the remaining fractional part after the first step)
  # Why "approximate"? Because with pynest layers, there are only two ways to specify
  # the number of axons in a connection:
  #    1) with an integer, specified with respect to each source (alt. target) neurons
  #    2) as a probability
  # Here, we have a fractional part - not an integer number - so that leaves us option 2.
  # However, because the new axonal contacts are drawn at random, we will not have the
  # exact number of connections
  float_inDegree = inDegree - integer_inDegree
  remaining_connections = np.round(float_inDegree * len(Pop[destName]))
  if remaining_connections > 0:
    printv('Adding '+str(remaining_connections)+' remaining connections with rule `fixed_total_number`\n')
    float_conndict = conndict.copy()
    float_conndict.update({'kernel': 1. / (nbSim[sourceName] * float(remaining_connections))})
    nesttopo.ConnectLayers(Topo[sourceName], Topo[destName], float_conndict)

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
    if stochastic_delays != None and delay > 0:
      printv('Using stochastic delays in mass-mirror')
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
#   if 'outDegreeAbs': `redundancy` is number of axonal contacts from an individual Src neuron onto a single Tgt neuron
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
# LCGDelays: shall we use the delays obtained by (Liénard, Cos, Girard, in prep) or not (default = True)
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
# Establishes a topological connection between two populations
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
def connectTopoMC(type, nameSrc, nameTgt, projType, redundancy, RedundancyType, LCGDelays=True, gain=1., source_channels = None, stochastic_delays=None, spreads=None, verbose=False):

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
    mass_connect_topo(nameSrc, nameTgt, lbl, inDegree, recType[lRecType[0]], W[lRecType[0]], delay, spread=spreads[0], stochastic_delays = stochastic_delays) # 0.5 spread for now
  elif projType == 'diffuse': # if projections diffused, input connections are shared among each possible input channel equally
    mass_connect_topo(nameSrc, nameTgt, lbl, inDegree, recType[lRecType[0]], W[lRecType[0]], delay, spread=spreads[1], stochastic_delays = stochastic_delays) # for now, arbitrary high spread

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
    print('\tCompare with the effective chosen inDegree   : '+str(inDegree))

  # attenuation due to the distance from the receptors to the soma of tgt:
  attenuation = cosh(LX[nameTgt]*(1-p[nameSrc+'->'+nameTgt])) / cosh(LX[nameTgt])

  w={}
  for r in listRecType:
    w[r] = nu / float(inDegree) * attenuation * wPSP[recType[r]-1] * gain

  return w

#-------------------------------------------------------------------------------

#rnd.seed(17)
#nest.SetKernelStatus({'local_num_threads':2, "data_path": "log/", "overwrite_files":True})
#nest.SetKernelStatus({'local_num_threads':2, "data_path": "log/"})

dt = 0.01 # ms
simDuration = 10000. # in ms

# Acceptable firing rate ranges (FRR) in normal and deactivation experiments
# extracted from LG14 Table 5

FRRNormal = {'MSN': [0,1],
             'FSI': [7.8,14.0], # the refined constraint of 10.9 +/- 3.1 Hz was extracted from the following papers: Adler et al., 2016; Yamada et al., 2016 (summarizing date from three different experiments); and Marche and Apicella, 2017 # old values: [0,20]
             'STN': [15.2,22.8],
             'GPe': [55.7,74.5],
             'Arky': [55.7,74.5],
             'Prot': [55.7,74.5],
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

FRRAnt = {'Arky':FRRGPe,'Prot':FRRGPe,'GPe':FRRGPe,'GPi':FRRGPi}

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

if params['splitGPe']:
  NUCLEI=['MSN','FSI','STN','Arky','Prot','GPi']
else:
  NUCLEI=['MSN','FSI','STN','GPe','GPi']

# Number of neurons in the real macaque brain
# one hemisphere only, based on Hardman et al. 2002 paper, except for striatum & CM/Pf
neuronCounts={'MSN': 26448.0E3,
              'FSI':   532.0E3,
              'STN':    77.0E3,
              'GPe':   251.0E3,
              'Arky':  251.0E3,
              'Prot':  251.0E3,
              'GPi':   143.0E3,
              'CMPf':   86.0E3,
              'CSN': None, 'PTN': None # prevents key error
             }

# Number of neurons that will be simulated
nbSim = {'MSN': 0.,
         'FSI': 0.,
         'STN': 0.,
         'GPe': 0.,
         'Arky': 0.,
         'Prot': 0.,
         'GPi': 0.,
         'CMPf':0.,
         'CSN': 0.,
         'PTN': 0.,}

# P(X->Y): probability that a given neuron from X projects to at least neuron of Y
P = {'MSN->GPe': 1.,
     'MSN->Arky': 1.,
     'MSN->Prot': 1.,
     'MSN->GPi': 0.82,
     'MSN->MSN': 1.,

     'FSI->MSN': 1.,
     'FSI->FSI': 1.,

     'STN->GPe':  0.83,
     'STN->Arky': 0.83,
     'STN->Prot': 0.83,
     'STN->GPi':  0.72,
     'STN->MSN':  0.17,
     'STN->FSI':  0.17,

     'GPe->STN': 1.,
     'GPe->GPe': 0.84,
     'GPe->GPi': 0.84,
     'GPe->MSN': 0.16,
     'GPe->FSI': 0.16,

     'Arky->Arky': 0.84,
     'Arky->Prot': 0.84,
     'Arky->MSN': 0.16,
     'Arky->FSI': 0.16,

     'Prot->STN': 1.,
     'Prot->Arky': 0.84,
     'Prot->Prot': 0.84,
     'Prot->GPi': 0.84,

     'CSN->MSN': 1.,
     'CSN->FSI': 1.,

     'PTN->MSN': 1.,
     'PTN->FSI': 1.,
     'PTN->STN': 1.,

     'CMPf->STN': 1.,
     'CMPf->MSN': 1.,
     'CMPf->FSI': 1.,
     'CMPf->GPe': 1.,
     'CMPf->Arky': 1.,
     'CMPf->Prot': 1.,
     'CMPf->GPi': 1.,}

# alpha X->Y: average number of synaptic contacts made by one neuron of X to one neuron of Y, when there is a connexion
# for the moment set from one specific parameterization, should be read from Jean's solution file
alpha = {'MSN->GPe':   171,
         'MSN->Arky':   171,
         'MSN->Prot':   171,
         'MSN->GPi':   210,
         'MSN->MSN':   210,

         'FSI->MSN':  4362,
         'FSI->FSI':   116,

         'STN->GPe':   428,
         'STN->Arky':   428,
         'STN->Prot':   428,
         'STN->GPi':   233,
         'STN->MSN':     0,
         'STN->FSI':    91,

         'GPe->STN':    19,
         'GPe->GPe':    38,
         'GPe->GPi':    16,
         'GPe->MSN':     0,
         'GPe->FSI':   353,

         'Arky->Arky':    38,
         'Arky->Prot':    38,
         'Arky->MSN':     0,
         'Arky->FSI':   353,

         'Prot->STN':    19,
         'Prot->Arky':    38,
         'Prot->Prot':    38,
         'Prot->GPi':    16,

         'CSN->MSN':   342, # here, represents directly \nu
         'CSN->FSI':   250, # here, represents directly \nu

         'PTN->MSN':     5, # here, represents directly \nu
         'PTN->FSI':     5, # here, represents directly \nu
         'PTN->STN':   259, # here, represents directly \nu

         'CMPf->MSN': 4965,
         'CMPf->FSI': 1053,
         'CMPf->STN':   76,
         'CMPf->GPe':   79,
         'CMPf->Arky':   79,
         'CMPf->Prot':   79,
         'CMPf->GPi':  131,}

# p(X->Y): relative distance on the dendrite from the soma, where neurons rom X projects to neurons of Y
# Warning: p is not P!
p = {'MSN->GPe':  0.48,
     'MSN->Arky':  0.48,
     'MSN->Prot':  0.48,
     'MSN->GPi':  0.59,
     'MSN->MSN':  0.77,

     'FSI->MSN':  0.19,
     'FSI->FSI':  0.16,

     'STN->GPe':  0.30,
     'STN->Prot':  0.30,
     'STN->Arky':  0.30,
     'STN->GPi':  0.59,
     'STN->MSN':  0.16,
     'STN->FSI':  0.41,

     'GPe->STN':  0.58,
     'GPe->GPe':  0.01,
     'GPe->GPi':  0.13,
     'GPe->MSN':  0.06,
     'GPe->FSI':  0.58,

     'Arky->Arky':  0.01,
     'Arky->Prot':  0.01,
     'Arky->MSN':  0.06,
     'Arky->FSI':  0.58,

     'Prot->STN':  0.58,
     'Prot->Arky':  0.01,
     'Prot->Prot':  0.01,
     'Prot->GPi':  0.13,

     'CSN->MSN':  0.95,
     'CSN->FSI':  0.82,

     'PTN->MSN':  0.98,
     'PTN->FSI':  0.70,
     'PTN->STN':  0.97,

     'CMPf->STN': 0.46,
     'CMPf->MSN': 0.27,
     'CMPf->FSI': 0.06,
     'CMPf->GPe': 0.00,
     'CMPf->Arky': 0.00,
     'CMPf->Prot': 0.00,
     'CMPf->GPi': 0.48,}

# electrotonic constant L computation:
dx={'MSN':1.E-6,'FSI':1.5E-6,'STN':1.5E-6,'GPe':1.7E-6,'Arky':1.7E-6,'Prot':1.7E-6,'GPi':1.2E-6}
lx={'MSN':619E-6,'FSI':961E-6,'STN':750E-6,'GPe':865E-6,'Arky':865E-6,'Prot':865E-6,'GPi':1132E-6}
LX={}
for n in NUCLEI:
    LX[n]=lx[n]*sqrt((4*Ri)/(dx[n]*Rm))

# tau: communication delays
tau = {'MSN->GPe':    7.,
       'MSN->Arky':    7.,
       'MSN->Prot':    7.,
       'MSN->GPi':   11.,
       'MSN->MSN':    1.,

       'FSI->MSN':    1.,
       'FSI->FSI':    1.,

       'STN->GPe':    3.,
       'STN->Arky':    3.,
       'STN->Prot':    3.,
       'STN->GPi':    3.,
       'STN->MSN':    3.,
       'STN->FSI':    3.,

       'GPe->STN':   10.,
       'GPe->GPe':    1.,
       'GPe->GPi':    3.,
       'GPe->MSN':    3.,
       'GPe->FSI':    3.,

       'Arky->Arky':    1.,
       'Arky->Prot':    1.,
       'Arky->MSN':    3.,
       'Arky->FSI':    3.,

       'Prot->STN':   10.,
       'Prot->Arky':    1.,
       'Prot->Prot':    1.,
       'Prot->GPi':    3.,

       'CSN->MSN':    7.,
       'CSN->FSI':    7.,

       'PTN->MSN':    3.,
       'PTN->FSI':    3.,
       'PTN->STN':    3.,

       'CMPf->MSN':   7.,
       'CMPf->FSI':   7.,
       'CMPf->STN':   7.,#4
       'CMPf->GPe':   7.,#5
       'CMPf->Arky':   7.,#5
       'CMPf->Prot':   7.,#5
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
                'tau_syn':   tau_syn,}

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

Arkyparams = {'tau_m':        14.0, # 20 -> 14 based on Johnson & McIntyre 2008, JNphy)
             'V_th':         11.0, # value of the LG14 example model, table 9
             'C_m':          14.0  # so that R_m=1, C_m=tau_m
            }

Protparams = {'tau_m':        14.0, # 20 -> 14 based on Johnson & McIntyre 2008, JNphy)
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
            'Arky':Arkyparams,
            'Prot':Protparams,
            'GPi':GPiparams}

Pop = {}
Fake= {} # Fake contains the Poisson Generators, that will feed the parrot_neurons, stored in Pop
ConnectMap = {} # when connections are drawn, in "create()", they are stored here so as to be re-usable
Topo = {} # stores the neuron coordinates (topological map)

# the dictionary used to store the desired discharge rates of the various Poisson generators that will be used as external inputs
rate = {'CSN':   2.  ,
        'PTN':  15.  ,
        'CMPf':  4.  ,
        'MSN':   0.25, # MSN and the following will be used when the corresponding nucleus is not explicitely simulated
        'FSI':  16.6 ,
        'STN':  14.3 ,
        'GPe':  62.6 ,
        'Arky':  62.6 ,
        'Prot':  62.6 ,
        'GPi':  64.2 ,
        }

#---------------------------
def main():

  # Pop is the dictionary that will contain the Nest IDs of all populations in the model
  #-------------------------
  print('Creating neurons')

  # creation of STN neurons
  #-------------------------
  nbSim['STN']=10.
  print('* STN: '+str(nbSim['STN'])+' neurons with parameters: '+str(BGparams['STN']))

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
  print('* PTN: '+str(nbSim['PTN'])+' Poisson generators with avg rate: '+str(rate['PTN']))
  Pop['PTN']  = nest.Create('poisson_generator',int(nbSim['PTN']))
  nest.SetStatus(Pop['PTN'],{'rate':rate['PTN']})

  connect('ex','PTN','STN', inDegree=5)

  # CMPf
  #-------------------------
  nbSim['CMPf']=nbSim['STN']
  rate['CMPf']=  4.
  print('* CMPf: '+str(nbSim['CMPf'])+' Poisson generators with avg rate: '+str(rate['CMPf']))
  Pop['CMPf'] = nest.Create('poisson_generator',int(nbSim['CMPf']))
  nest.SetStatus(Pop['CMPf'],{'rate': rate['CMPf']})

  connect('ex','CMPf','STN', inDegree=1)

  # Fake GPe
  #-------------------------
  nbSim['GPe'] = int(neuronCounts['GPe']/neuronCounts['STN']) * nbSim['STN']
  rate['GPe']= 62.6
  print('* GPe: '+str(nbSim['GPe'])+' Poisson generators with avg rate: '+str(rate['GPe']))
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
  print('\n Spike Detector n_events '+str(nest.GetStatus(spkDetect, 'n_events')[0]))
  expeRate = nest.GetStatus(spkDetect, 'n_events')[0] / float(nbSim['STN']*simDuration)
  print('\n Rate:'+str(expeRate*1000)+' Hz')


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
