#!/usr/bin/python
# -*- coding: utf-8 -*-

import nstrand

from LGneurons import *
from modelParams import *
import nest.raster_plot
import nest.voltage_trace
import pylab as pl
import sys

import csv


#------------------------------------------
# Creates the populations of neurons necessary to simulate a BG circuit
#------------------------------------------
def createBG():
  #==========================
  # Creation of neurons
  #-------------------------
  print '\nCreating neurons\n================'

  # single or multi-channel?
  if params['nbCh'] == 1:
    create_pop = lambda *args, **kwargs: create(*args, **kwargs)
    update_Ie = lambda p: nest.SetStatus(Pop[p],{"I_e":params['Ie'+p]})
  else:
    create_pop = lambda *args, **kwargs: createMC(nbCh=params['nbCh'], *args, **kwargs)
    update_Ie = lambda p: [nest.SetStatus(Pop[p][i],{"I_e":params['Ie'+p]}) for i in range(len(Pop[p]))]

  nbSim['MSN'] = params['nbMSN']
  create_pop('MSN')
  update_Ie('MSN')

  nbSim['FSI'] = params['nbFSI']
  create_pop('FSI')
  update_Ie('FSI')

  nbSim['STN'] = params['nbSTN']
  create_pop('STN')
  update_Ie('STN')

  nbSim['GPe'] = params['nbGPe']
  create_pop('GPe')
  update_Ie('GPe')

  nbSim['GPi'] = params['nbGPi']
  create_pop('GPi')
  update_Ie('GPi')

  parrot = True # switch to False at your risks & perils...
  nbSim['CSN'] = params['nbCSN']
  create_pop('CSN', fake=True, parrot=parrot)

  nbSim['PTN'] = params['nbPTN']
  create_pop('PTN', fake=True, parrot=parrot)

  nbSim['CMPf'] = params['nbCMPf']
  create_pop('CMPf', fake=True, parrot=params['parrotCMPf']) # was: False

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

  # single or multi-channel?
  if params['nbCh'] == 1:
    connect_pop = lambda *args, **kwargs: connect(*args, **kwargs)
  else:
    connect_pop = lambda *args, **kwargs: connectMC(*args, **kwargs)

  #-------------------------
  # connection of populations
  #-------------------------
  print '\nConnecting neurons\n================'
  print "**",antag,"antagonist injection in",antagInjectionSite,"**"
  print '* MSN Inputs'
  connect_pop('ex','CSN','MSN', projType=params['cTypeCSNMSN'], inDegree= params['inDegCSNMSN'], gain=G['MSN'])
  connect_pop('ex','PTN','MSN', projType=params['cTypePTNMSN'], inDegree= params['inDegPTNMSN'], gain=G['MSN'])
  connect_pop('ex','CMPf','MSN',projType=params['cTypeCMPfMSN'],inDegree= params['inDegCMPfMSN'],gain=G['MSN'])
  connect_pop('in','MSN','MSN', projType=params['cTypeMSNMSN'], inDegree= params['inDegMSNMSN'], gain=G['MSN'])
  connect_pop('in','FSI','MSN', projType=params['cTypeFSIMSN'], inDegree= params['inDegFSIMSN'], gain=G['MSN'])
  # some parameterizations from LG14 have no STN->MSN or GPe->MSN synaptic contacts
  if alpha['STN->MSN'] != 0:
    print "alpha['STN->MSN']",alpha['STN->MSN']
    connect_pop('ex','STN','MSN', projType=params['cTypeSTNMSN'], inDegree= params['inDegSTNMSN'], gain=G['MSN'])
  if alpha['GPe->MSN'] != 0:
    print "alpha['GPe->MSN']",alpha['GPe->MSN']
    connect_pop('in','GPe','MSN', projType=params['cTypeGPeMSN'], inDegree= params['inDegGPeMSN'], gain=G['MSN'])

  print '* FSI Inputs'
  connect_pop('ex','CSN','FSI', projType=params['cTypeCSNFSI'], inDegree= params['inDegCSNFSI'], gain=G['FSI'])
  connect_pop('ex','PTN','FSI', projType=params['cTypePTNFSI'], inDegree= params['inDegPTNFSI'], gain=G['FSI'])
  if alpha['STN->FSI'] != 0:
    connect_pop('ex','STN','FSI', projType=params['cTypeSTNFSI'],inDegree= params['inDegSTNFSI'],gain=G['FSI'])
  connect_pop('in','GPe','FSI', projType=params['cTypeGPeFSI'], inDegree= params['inDegGPeFSI'], gain=G['FSI'])
  connect_pop('ex','CMPf','FSI',projType=params['cTypeCMPfFSI'],inDegree= params['inDegCMPfFSI'],gain=G['FSI'])
  connect_pop('in','FSI','FSI', projType=params['cTypeFSIFSI'], inDegree= params['inDegFSIFSI'], gain=G['FSI'])

  print '* STN Inputs'
  connect_pop('ex','PTN','STN', projType=params['cTypePTNSTN'], inDegree= params['inDegPTNSTN'],  gain=G['STN'])
  connect_pop('ex','CMPf','STN',projType=params['cTypeCMPfSTN'],inDegree= params['inDegCMPfSTN'], gain=G['STN'])
  connect_pop('in','GPe','STN', projType=params['cTypeGPeSTN'], inDegree= params['inDegGPeSTN'],  gain=G['STN'])

  print '* GPe Inputs'
  if antagInjectionSite == 'GPe':
    if   antag == 'AMPA':
      connect_pop('NMDA','CMPf','GPe',projType=params['cTypeCMPfGPe'],inDegree= params['inDegCMPfGPe'],gain=G['GPe'])
      connect_pop('NMDA','STN','GPe', projType=params['cTypeSTNGPe'], inDegree= params['inDegSTNGPe'], gain=G['GPe'])
      connect_pop('in','MSN','GPe',   projType=params['cTypeMSNGPe'], inDegree= params['inDegMSNGPe'], gain=G['GPe'])
      connect_pop('in','GPe','GPe',   projType=params['cTypeGPeGPe'], inDegree= params['inDegGPeGPe'], gain=G['GPe'])
    elif antag == 'NMDA':
      connect_pop('AMPA','CMPf','GPe',projType=params['cTypeCMPfGPe'],inDegree= params['inDegCMPfGPe'],gain=G['GPe'])
      connect_pop('AMPA','STN','GPe', projType=params['cTypeSTNGPe'], inDegree= params['inDegSTNGPe'], gain=G['GPe'])
      connect_pop('in','MSN','GPe',   projType=params['cTypeMSNGPe'], inDegree= params['inDegMSNGPe'], gain=G['GPe'])
      connect_pop('in','GPe','GPe',   projType=params['cTypeGPeGPe'], inDegree= params['inDegGPeGPe'], gain=G['GPe'])
    elif antag == 'AMPA+GABAA':
      connect_pop('NMDA','CMPf','GPe',projType=params['cTypeCMPfGPe'],inDegree= params['inDegCMPfGPe'],gain=G['GPe'])
      connect_pop('NMDA','STN','GPe', projType=params['cTypeSTNGPe'], inDegree= params['inDegSTNGPe'], gain=G['GPe'])
    elif antag == 'GABAA':
      connect_pop('ex','CMPf','GPe',projType=params['cTypeCMPfGPe'],inDegree= params['inDegCMPfGPe'],gain=G['GPe'])
      connect_pop('ex','STN','GPe', projType=params['cTypeSTNGPe'], inDegree= params['inDegSTNGPe'], gain=G['GPe'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connect_pop('ex','CMPf','GPe',projType=params['cTypeCMPfGPe'],inDegree= params['inDegCMPfGPe'],gain=G['GPe'])
    connect_pop('ex','STN','GPe', projType=params['cTypeSTNGPe'], inDegree= params['inDegSTNGPe'], gain=G['GPe'])
    connect_pop('in','MSN','GPe', projType=params['cTypeMSNGPe'], inDegree= params['inDegMSNGPe'], gain=G['GPe'])
    connect_pop('in','GPe','GPe', projType=params['cTypeGPeGPe'], inDegree= params['inDegGPeGPe'], gain=G['GPe'])

  print '* GPi Inputs'
  if antagInjectionSite =='GPi':
    if   antag == 'AMPA+NMDA+GABAA':
      pass
    elif antag == 'NMDA':
      connect_pop('in','MSN','GPi',   projType=params['cTypeMSNGPi'], inDegree= params['inDegMSNGPi'], gain=G['GPi'])
      connect_pop('AMPA','STN','GPi', projType=params['cTypeSTNGPi'], inDegree= params['inDegSTNGPi'], gain=G['GPi'])
      connect_pop('in','GPe','GPi',   projType=params['cTypeGPeGPi'], inDegree= params['inDegGPeGPi'], gain=G['GPi'])
      connect_pop('AMPA','CMPf','GPi',projType=params['cTypeCMPfGPi'],inDegree= params['inDegCMPfGPi'],gain=G['GPi'])
    elif antag == 'NMDA+AMPA':
      connect_pop('in','MSN','GPi', projType=params['cTypeMSNGPi'],inDegree= params['inDegMSNGPi'], gain=G['GPi'])
      connect_pop('in','GPe','GPi', projType=params['cTypeGPeGPi'],inDegree= params['inDegGPeGPi'], gain=G['GPi'])
    elif antag == 'AMPA':
      connect_pop('in','MSN','GPi',   projType=params['cTypeMSNGPi'], inDegree= params['inDegMSNGPi'], gain=G['GPi'])
      connect_pop('NMDA','STN','GPi', projType=params['cTypeSTNGPi'], inDegree= params['inDegSTNGPi'], gain=G['GPi'])
      connect_pop('in','GPe','GPi',   projType=params['cTypeGPeGPi'], inDegree= params['inDegGPeGPi'], gain=G['GPi'])
      connect_pop('NMDA','CMPf','GPi',projType=params['cTypeCMPfGPi'],inDegree= params['inDegCMPfGPi'],gain=G['GPi'])
    elif antag == 'GABAA':
      connect_pop('ex','STN','GPi', projType=params['cTypeSTNGPi'], inDegree= params['inDegSTNGPi'], gain=G['GPi'])
      connect_pop('ex','CMPf','GPi',projType=params['cTypeCMPfGPi'],inDegree= params['inDegCMPfGPi'],gain=G['GPi'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connect_pop('in','MSN','GPi', projType=params['cTypeMSNGPi'], inDegree= params['inDegMSNGPi'], gain=G['GPi'])
    connect_pop('ex','STN','GPi', projType=params['cTypeSTNGPi'], inDegree= params['inDegSTNGPi'], gain=G['GPi'])
    connect_pop('in','GPe','GPi', projType=params['cTypeGPeGPi'], inDegree= params['inDegGPeGPi'], gain=G['GPi'])
    connect_pop('ex','CMPf','GPi',projType=params['cTypeCMPfGPi'],inDegree= params['inDegCMPfGPi'],gain=G['GPi'])


#------------------------------------------
# Re-weight a specific connection, characterized by a source, a target, and a receptor
# Returns the previous value of that connection (useful for 'reactivating' after a deactivation experiment)
#------------------------------------------
def alter_connection(src, tgt, tgt_receptor, altered_weight):
  if params['nbCh'] != 1:
    raise NotImplementedError('Altering connection is implemented only in the one-channel case')
  recTypeEquiv = {'AMPA':1,'NMDA':2,'GABA':3, 'GABAA':3} # adds 'GABAA'
  # check that we have this connection in the current network
  conns_in = nest.GetConnections(source=Pop[src], target=Pop[tgt])
  if len(conns_in):
    receptors = nest.GetStatus(conns_in, keys='receptor')
    previous_weights = nest.GetStatus(conns_in, keys='weight')
    rec_nb = recTypeEquiv[tgt_receptor]
    if isinstance(altered_weight, int):
      altered_weights = [altered_weight] * len(receptors)
    elif len(altered_weight) == len(receptors):
      altered_weights = altered_weight # already an array
    else:
      raise LookupError('Wrong size for the `altered_weights` variable (should be scalar or a list with as many items as there are synapses in that connection - including non-targeted receptors)')
    new_weights = [{'weight': float(previous_weights[i])} if receptors[i] != rec_nb else {'weight': float(altered_weights[i])} for i in range(len(receptors))] # replace the weights for the targeted receptor
    nest.SetStatus(conns_in, new_weights)
    return previous_weights
  return None

#------------------------------------------
# gets the nuclei involved in deactivation experiments in GPe/GPi
#------------------------------------------
def get_afferents(a):
  GABA_afferents = ['MSN', 'GPe'] # afferents with gabaergic connections
  GLUT_afferents = ['STN', 'CMPf'] # afferents with glutamatergic connections
  if a == 'GABAA':
    afferents = GABA_afferents
  elif a == 'AMPA+GABAA':
    afferents = GABA_afferents + GLUT_afferents
  elif a == 'AMPA+NMDA+GABAA':
    afferents = GABA_afferents + GLUT_afferents
  else:
    afferents = GLUT_afferents
  return afferents

#------------------------------------------
# deactivate connections based on antagonist experiment
#------------------------------------------
def deactivate(site, a):
  ww = {}
  for src in get_afferents(a):
    ww[src] = None
    for rec in a.split('+'):
      w = alter_connection(src, site, rec, 0)
      if ww[src] == None:
        ww[src] = w # keep the original weights only once
  return ww

#------------------------------------------
# reactivate connections based on antagonist experiment
#------------------------------------------
def reactivate(site, a, ww):
  for src in get_afferents(a):
    for rec in a.split('+'):
      alter_connection(src, site, rec, ww[src])

#------------------------------------------
# Instantiate the BG network according to the `params` dictionnary
# For now, this instantiation respects the hardcoded antagonist injection sites
# In the future, these will be handled by changing the network weights
#------------------------------------------
def instantiate_BG(params={}, antagInjectionSite='none', antag=''):
  nest.ResetKernel()
  dataPath='log/'
  if 'nbcpu' in params:
    nest.SetKernelStatus({'local_num_threads': params['nbcpu']})

  nstrand.set_seed(params['nestSeed'], params['pythonSeed']) # sets the seed for the BG construction

  nest.SetKernelStatus({"data_path": dataPath})
  initNeurons()

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

