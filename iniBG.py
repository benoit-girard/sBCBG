#!/usr/bin/python
# -*- coding: utf-8 -*-

import nstrand

from LGneurons import *
from modelParams import *
import os
import nest.raster_plot
import nest.voltage_trace
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
    def create_pop(*args, **kwargs):
      if 'nbCh' in kwargs.keys():
        # remove the extra arg
        kwargs.pop("nbCh", None)
      create(*args, **kwargs)
    update_Ie = lambda p: nest.SetStatus(Pop[p],{"I_e":params['Ie'+p]})
  else:
    def create_pop(*args, **kwargs):
      if 'nbCh' not in kwargs.keys():
        # enforce the default
        kwargs['nbCh'] = params['nbCh']
      createMC(*args, **kwargs)
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

  if params['splitGPe']:
    nbSim['Arky'] = params['nbArky']
    create_pop('Arky')
    update_Ie('Arky')
  
    nbSim['Prot'] = params['nbProt']
    create_pop('Prot')
    update_Ie('Prot')
  else:
    nbSim['GPe'] = params['nbGPe']
    create_pop('GPe')
    update_Ie('GPe')
      
      
  nbSim['GPi'] = params['nbGPi']
  create_pop('GPi')
  update_Ie('GPi')

  parrot = True # switch to False at your risks & perils...
  nbSim['CSN'] = params['nbCSN']
  if 'nbCues' in params.keys():
    # cue channels are present
    CSNchannels = params['nbCh']+params['nbCues']
  else:
    CSNchannels = params['nbCh']
  create_pop('CSN', nbCh=CSNchannels, fake=True, parrot=parrot)

  nbSim['PTN'] = params['nbPTN']
  create_pop('PTN', fake=True, parrot=parrot)

  nbSim['CMPf'] = params['nbCMPf']
  create_pop('CMPf', fake=True, parrot=params['parrotCMPf']) # was: False

  print "Number of simulated neurons:", nbSim

#------------------------------------------
# Connects the populations of a previously created multi-channel BG circuit
#------------------------------------------
def connectBG(antagInjectionSite,antag):

  # single or multi-channel?
  if params['nbCh'] == 1:
    connect_pop = lambda *args, **kwargs: connect(*args, RedundancyType=params['RedundancyType'], stochastic_delays=params['stochastic_delays'], **kwargs)
  else:
    def connect_pop(*args, **kwargs):
      if 'source_channels' not in kwargs.keys():
        # enforce the default
        kwargs['source_channels'] = range(params['nbCh'])
      return connectMC(*args, RedundancyType=params['RedundancyType'], stochastic_delays=params['stochastic_delays'], **kwargs)

  #-------------------------
  # connection of populations
  #-------------------------
  print '\nConnecting neurons\n================'
  print "**",antag,"antagonist injection in",antagInjectionSite,"**"
  print '* MSN Inputs'
  if 'nbCues' not in params.keys():
    # usual case: CSN have as the same number of channels than the BG nuclei
    CSN_MSN = connect_pop('ex','CSN','MSN', projType=params['cTypeCSNMSN'], redundancy=params['redundancyCSNMSN'], gain=params['GCSNMSN'])
  else:
    # special case: extra 'cue' channels that target MSN
    CSN_MSN = connect_pop('ex','CSN','MSN', projType=params['cTypeCSNMSN'], redundancy=params['redundancyCSNMSN'], gain=Params['GCSNMSN']/2., source_channels=range(params['nbCh']))
    connect_pop('ex','CSN','MSN', projType='diffuse', redundancy=params['redundancyCSNMSN'], gain=params['GCSNMSN']/2., source_channels=range(params['nbCh'], params['nbCh']+params['nbCues']))
  PTN_MSN = connect_pop('ex','PTN','MSN', projType=params['cTypePTNMSN'], redundancy= params['redundancyPTNMSN'], gain=params['GPTNMSN'])
  CMPf_MSN = connect_pop('ex','CMPf','MSN',projType=params['cTypeCMPfMSN'],redundancy= params['redundancyCMPfMSN'],gain=params['GCMPfMSN'])
  connect_pop('in','MSN','MSN', projType=params['cTypeMSNMSN'], redundancy= params['redundancyMSNMSN'], gain=params['GMSNMSN'])
  connect_pop('in','FSI','MSN', projType=params['cTypeFSIMSN'], redundancy= params['redundancyFSIMSN'], gain=params['GFSIMSN'])
  # some parameterizations from LG14 have no STN->MSN or GPe->MSN synaptic contacts
  if alpha['STN->MSN'] != 0:
    print "alpha['STN->MSN']",alpha['STN->MSN']
    connect_pop('ex','STN','MSN', projType=params['cTypeSTNMSN'], redundancy= params['redundancySTNMSN'], gain=params['GSTNMSN'])
  if alpha['GPe->MSN'] != 0:
    if params['splitGPe']:
      print "alpha['Arky->MSN']",alpha['Arky->MSN']
      connect_pop('in','Arky','MSN', projType=params['cTypeArkyMSN'], redundancy= params['redundancyArkyMSN'], gain=params['GArkyMSN'])
    else:
      print "alpha['GPe->MSN']",alpha['GPe->MSN']
      connect_pop('in','GPe','MSN', projType=params['cTypeGPeMSN'], redundancy= params['redundancyGPeMSN'], gain=params['GGPeMSN'])

  print '* FSI Inputs'
  connect_pop('ex','CSN','FSI', projType=params['cTypeCSNFSI'], redundancy= params['redundancyCSNFSI'], gain=params['GCSNFSI'])
  connect_pop('ex','PTN','FSI', projType=params['cTypePTNFSI'], redundancy= params['redundancyPTNFSI'], gain=params['GPTNFSI'])
  if alpha['STN->FSI'] != 0:
    connect_pop('ex','STN','FSI', projType=params['cTypeSTNFSI'],redundancy= params['redundancySTNFSI'],gain=params['GSTNFSI'])
  if params['splitGPe']:
    connect_pop('in','Arky','FSI', projType=params['cTypeArkyFSI'], redundancy= params['redundancyArkyFSI'], gain=params['GArkyFSI'])
  else:
    connect_pop('in','GPe','FSI', projType=params['cTypeGPeFSI'], redundancy= params['redundancyGPeFSI'], gain=params['GGPeFSI'])
    
  connect_pop('ex','CMPf','FSI',projType=params['cTypeCMPfFSI'],redundancy= params['redundancyCMPfFSI'],gain=params['GCMPfFSI'])
  connect_pop('in','FSI','FSI', projType=params['cTypeFSIFSI'], redundancy= params['redundancyFSIFSI'], gain=params['GFSIFSI'])

  print '* STN Inputs'
  connect_pop('ex','PTN','STN', projType=params['cTypePTNSTN'], redundancy= params['redundancyPTNSTN'],  gain=params['GPTNSTN'])
  connect_pop('ex','CMPf','STN',projType=params['cTypeCMPfSTN'],redundancy= params['redundancyCMPfSTN'], gain=params['GCMPfSTN'])
  if params['splitGPe']:
    connect_pop('in','Prot','STN', projType=params['cTypeProtSTN'], redundancy= params['redundancyProtSTN'],  gain=params['GProtSTN'])
  else:
    connect_pop('in','GPe','STN', projType=params['cTypeGPeSTN'], redundancy= params['redundancyGPeSTN'],  gain=params['GGPeSTN'])
  
  if params['splitGPe']:
      print '* Arky Inputs'
      if 'fakeArkyRecurrent' not in params.keys():
        # usual case: Arky's recurrent collaterals are handled normally
        Arky_recurrent_source = 'Arky'
      else:
        # here collaterals are simulated with Poisson train spikes firing at the frequency given by params['fakeArkyRecurrent']
        rate['Fake_Arky'] = float(params['fakeArkyRecurrent'])
        for nucleus_dict in [nbSim, neuronCounts]:
          nucleus_dict['Fake_Arky'] = nucleus_dict['Arky']
        for connection_dict in [P, alpha, p, tau]:
          connection_dict['Fake_Arky->Arky'] = connection_dict['Arky->Arky']
        if params['nbCh'] == 1:
          create('Fake_Arky', fake=True, parrot=True)
        else:
          createMC('Fake_Arky', params['nbCh'], fake=True, parrot=True)
        Arky_recurrent_source = 'Fake_Arky'
      if antagInjectionSite == 'GPe':
        if   antag == 'AMPA':
          connect_pop('NMDA','CMPf','Arky',projType=params['cTypeCMPfArky'],redundancy= params['redundancyCMPfArky'],gain=params['GCMPfArky'])
          connect_pop('NMDA','STN','Arky', projType=params['cTypeSTNArky'], redundancy= params['redundancySTNArky'], gain=params['GSTNArky'])
          connect_pop('in','MSN','Arky',   projType=params['cTypeMSNArky'], redundancy= params['redundancyMSNArky'], gain=params['GMSNArky'])
          connect_pop('in','Prot','Arky', projType=params['cTypeProtArky'], redundancy= params['redundancyProtArky'], gain=params['GProtArky'])
          connect_pop('in', Arky_recurrent_source, 'Arky', projType=params['cTypeArkyArky'], redundancy= params['redundancyArkyArky'], gain=params['GArkyArky'])
        elif antag == 'NMDA':
          connect_pop('AMPA','CMPf','Arky',projType=params['cTypeCMPfArky'],redundancy= params['redundancyCMPfArky'],gain=params['GCMPfArky'])
          connect_pop('AMPA','STN','Arky', projType=params['cTypeSTNArky'], redundancy= params['redundancySTNArky'], gain=params['GSTNArky'])
          connect_pop('in','MSN','Arky',   projType=params['cTypeMSNArky'], redundancy= params['redundancyMSNArky'], gain=params['GMSNArky'])
          connect_pop('in','Prot','Arky', projType=params['cTypeProtArky'], redundancy= params['redundancyProtArky'], gain=params['GProtArky'])
          connect_pop('in', Arky_recurrent_source, 'Arky', projType=params['cTypeArkyArky'], redundancy= params['redundancyArkyArky'], gain=params['GArkyArky'])
        elif antag == 'AMPA+GABAA':
          connect_pop('NMDA','CMPf','Arky',projType=params['cTypeCMPfArky'],redundancy= params['redundancyCMPfArky'],gain=params['GCMPfArky'])
          connect_pop('NMDA','STN','Arky', projType=params['cTypeSTNArky'], redundancy= params['redundancySTNArky'], gain=params['GSTNArky'])
        elif antag == 'GABAA':
          connect_pop('ex','CMPf','Arky',projType=params['cTypeCMPfArky'],redundancy= params['redundancyCMPfArky'],gain=params['GCMPfArky'])
          connect_pop('ex','STN','Arky', projType=params['cTypeSTNArky'], redundancy= params['redundancySTNArky'], gain=params['GSTNArky'])
        else:
          print antagInjectionSite,": unknown antagonist experiment:",antag
      else:
        connect_pop('ex','CMPf','Arky',projType=params['cTypeCMPfArky'],redundancy= params['redundancyCMPfArky'],gain=params['GCMPfArky'])
        connect_pop('ex','STN','Arky', projType=params['cTypeSTNArky'], redundancy= params['redundancySTNArky'], gain=params['GSTNArky'])
        connect_pop('in','MSN','Arky', projType=params['cTypeMSNArky'], redundancy= params['redundancyMSNArky'], gain=params['GMSNArky'])
        connect_pop('in','Prot','Arky', projType=params['cTypeProtArky'], redundancy= params['redundancyProtArky'], gain=params['GProtArky'])
        connect_pop('in', Arky_recurrent_source, 'Arky', projType=params['cTypeArkyArky'], redundancy= params['redundancyArkyArky'], gain=params['GArkyArky'])
        
      print '* Prot Inputs'
      if 'fakeProtRecurrent' not in params.keys():
        # usual case: Prot's recurrent collaterals are handled normally
        Prot_recurrent_source = 'Prot'
      else:
        # here collaterals are simulated with Poisson train spikes firing at the frequency given by params['fakeProtRecurrent']
        rate['Fake_Prot'] = float(params['fakeProtRecurrent'])
        for nucleus_dict in [nbSim, neuronCounts]:
          nucleus_dict['Fake_Prot'] = nucleus_dict['Prot']
        for connection_dict in [P, alpha, p, tau]:
          connection_dict['Fake_Prot->Prot'] = connection_dict['Prot->Prot']
        if params['nbCh'] == 1:
          create('Fake_Prot', fake=True, parrot=True)
        else:
          createMC('Fake_Prot', params['nbCh'], fake=True, parrot=True)
        Prot_recurrent_source = 'Fake_Prot'
      if antagInjectionSite == 'GPe':
        if   antag == 'AMPA':
          connect_pop('NMDA','CMPf','Prot',projType=params['cTypeCMPfProt'],redundancy= params['redundancyCMPfProt'],gain=params['GCMPfProt'])
          connect_pop('NMDA','STN','Prot', projType=params['cTypeSTNProt'], redundancy= params['redundancySTNProt'], gain=params['GSTNProt'])
          connect_pop('in','MSN','Prot',   projType=params['cTypeMSNProt'], redundancy= params['redundancyMSNProt'], gain=params['GMSNProt'])
          connect_pop('in','Arky','Prot',   projType=params['cTypeArkyProt'], redundancy= params['redundancyArkyProt'], gain=params['GArkyProt'])
          connect_pop('in', Prot_recurrent_source, 'Prot', projType=params['cTypeProtProt'], redundancy= params['redundancyProtProt'], gain=params['GProtProt'])
        elif antag == 'NMDA':
          connect_pop('AMPA','CMPf','Prot',projType=params['cTypeCMPfProt'],redundancy= params['redundancyCMPfProt'],gain=params['GCMPfProt'])
          connect_pop('AMPA','STN','Prot', projType=params['cTypeSTNProt'], redundancy= params['redundancySTNProt'], gain=params['GSTNProt'])
          connect_pop('in','MSN','Prot',   projType=params['cTypeMSNProt'], redundancy= params['redundancyMSNProt'], gain=params['GMSNProt'])
          connect_pop('in','Arky','Prot',   projType=params['cTypeArkyProt'], redundancy= params['redundancyArkyProt'], gain=params['GArkyProt'])
          connect_pop('in', Prot_recurrent_source, 'Prot', projType=params['cTypeProtProt'], redundancy= params['redundancyProtProt'], gain=params['GProtProt'])
        elif antag == 'AMPA+GABAA':
          connect_pop('NMDA','CMPf','Prot',projType=params['cTypeCMPfProt'],redundancy= params['redundancyCMPfProt'],gain=params['GCMPfProt'])
          connect_pop('NMDA','STN','Prot', projType=params['cTypeSTNProt'], redundancy= params['redundancySTNProt'], gain=params['GSTNProt'])
        elif antag == 'GABAA':
          connect_pop('ex','CMPf','Prot',projType=params['cTypeCMPfProt'],redundancy= params['redundancyCMPfProt'],gain=params['GCMPfProt'])
          connect_pop('ex','STN','Prot', projType=params['cTypeSTNProt'], redundancy= params['redundancySTNProt'], gain=params['GSTNProt'])
        else:
          print antagInjectionSite,": unknown antagonist experiment:",antag
      else:
        connect_pop('ex','CMPf','Prot',projType=params['cTypeCMPfProt'],redundancy= params['redundancyCMPfProt'],gain=params['GCMPfProt'])
        connect_pop('ex','STN','Prot', projType=params['cTypeSTNProt'], redundancy= params['redundancySTNProt'], gain=params['GSTNProt'])
        connect_pop('in','MSN','Prot', projType=params['cTypeMSNProt'], redundancy= params['redundancyMSNProt'], gain=params['GMSNProt'])
        connect_pop('in','Arky','Prot',   projType=params['cTypeArkyProt'], redundancy= params['redundancyArkyProt'], gain=params['GArkyProt'])
        connect_pop('in', Prot_recurrent_source, 'Prot', projType=params['cTypeProtProt'], redundancy= params['redundancyProtProt'], gain=params['GProtProt'])
        
  else:
      print '* GPe Inputs'
      if 'fakeGPeRecurrent' not in params.keys():
        # usual case: GPe's recurrent collaterals are handled normally
        GPe_recurrent_source = 'GPe'
      else:
        # here collaterals are simulated with Poisson train spikes firing at the frequency given by params['fakeGPeRecurrent']
        rate['Fake_GPe'] = float(params['fakeGPeRecurrent'])
        for nucleus_dict in [nbSim, neuronCounts]:
          nucleus_dict['Fake_GPe'] = nucleus_dict['GPe']
        for connection_dict in [P, alpha, p, tau]:
          connection_dict['Fake_GPe->GPe'] = connection_dict['GPe->GPe']
        if params['nbCh'] == 1:
          create('Fake_GPe', fake=True, parrot=True)
        else:
          createMC('Fake_GPe', params['nbCh'], fake=True, parrot=True)
        GPe_recurrent_source = 'Fake_GPe'
      if antagInjectionSite == 'GPe':
        if   antag == 'AMPA':
          connect_pop('NMDA','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=params['GCMPfGPe'])
          connect_pop('NMDA','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=params['GSTNGPe'])
          connect_pop('in','MSN','GPe',   projType=params['cTypeMSNGPe'], redundancy= params['redundancyMSNGPe'], gain=params['GMSNGPe'])
          connect_pop('in', GPe_recurrent_source, 'GPe', projType=params['cTypeGPeGPe'], redundancy= params['redundancyGPeGPe'], gain=params['GGPeGPe'])
        elif antag == 'NMDA':
          connect_pop('AMPA','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=params['GCMPfGPe'])
          connect_pop('AMPA','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=params['GSTNGPe'])
          connect_pop('in','MSN','GPe',   projType=params['cTypeMSNGPe'], redundancy= params['redundancyMSNGPe'], gain=params['GMSNGPe'])
          connect_pop('in', GPe_recurrent_source, 'GPe', projType=params['cTypeGPeGPe'], redundancy= params['redundancyGPeGPe'], gain=params['GGPeGPe'])
        elif antag == 'AMPA+GABAA':
          connect_pop('NMDA','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=params['GCMPfGPe'])
          connect_pop('NMDA','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=params['GSTNGPe'])
        elif antag == 'GABAA':
          connect_pop('ex','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=params['GCMPfGPe'])
          connect_pop('ex','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=params['GSTNGPe'])
        else:
          print antagInjectionSite,": unknown antagonist experiment:",antag
      else:
        connect_pop('ex','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=params['GCMPfGPe'])
        connect_pop('ex','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=params['GSTNGPe'])
        connect_pop('in','MSN','GPe', projType=params['cTypeMSNGPe'], redundancy= params['redundancyMSNGPe'], gain=params['GMSNGPe'])
        connect_pop('in', GPe_recurrent_source, 'GPe', projType=params['cTypeGPeGPe'], redundancy= params['redundancyGPeGPe'], gain=params['GGPeGPe'])

  print '* GPi Inputs'
  if antagInjectionSite =='GPi':
    if   antag == 'AMPA+NMDA+GABAA':
      pass
    elif antag == 'NMDA':
      connect_pop('in','MSN','GPi',   projType=params['cTypeMSNGPi'], redundancy= params['redundancyMSNGPi'], gain=params['GMSNGPi'])
      connect_pop('AMPA','STN','GPi', projType=params['cTypeSTNGPi'], redundancy= params['redundancySTNGPi'], gain=params['GSTNGPi'])
      if params['splitGPe']:
        connect_pop('in','Prot','GPi',   projType=params['cTypeProtGPi'], redundancy= params['redundancyProtGPi'], gain=params['GProtGPi'])
      else:
        connect_pop('in','GPe','GPi',   projType=params['cTypeGPeGPi'], redundancy= params['redundancyGPeGPi'], gain=params['GGPeGPi'])
      connect_pop('AMPA','CMPf','GPi',projType=params['cTypeCMPfGPi'],redundancy= params['redundancyCMPfGPi'],gain=params['GCMPfGPi'])
    elif antag == 'NMDA+AMPA':
      connect_pop('in','MSN','GPi', projType=params['cTypeMSNGPi'],redundancy= params['redundancyMSNGPi'], gain=params['GMSNGPi'])
      if params['splitGPe']:
        connect_pop('in','Prot','GPi',   projType=params['cTypeProtGPi'], redundancy= params['redundancyProtGPi'], gain=params['GProtGPi'])
      else:
        connect_pop('in','GPe','GPi',   projType=params['cTypeGPeGPi'], redundancy= params['redundancyGPeGPi'], gain=params['GGPeGPi'])
    elif antag == 'AMPA':
      connect_pop('in','MSN','GPi',   projType=params['cTypeMSNGPi'], redundancy= params['redundancyMSNGPi'], gain=params['GMSNGPi'])
      connect_pop('NMDA','STN','GPi', projType=params['cTypeSTNGPi'], redundancy= params['redundancySTNGPi'], gain=params['GSTNGPi'])
      if params['splitGPe']:
        connect_pop('in','Prot','GPi',   projType=params['cTypeProtGPi'], redundancy= params['redundancyProtGPi'], gain=params['GProtGPi'])
      else:
        connect_pop('in','GPe','GPi',   projType=params['cTypeGPeGPi'], redundancy= params['redundancyGPeGPi'], gain=params['GGPeGPi'])
      connect_pop('NMDA','CMPf','GPi',projType=params['cTypeCMPfGPi'],redundancy= params['redundancyCMPfGPi'],gain=params['GCMPfGPi'])
    elif antag == 'GABAA':
      connect_pop('ex','STN','GPi', projType=params['cTypeSTNGPi'], redundancy= params['redundancySTNGPi'], gain=params['GSTNGPi'])
      connect_pop('ex','CMPf','GPi',projType=params['cTypeCMPfGPi'],redundancy= params['redundancyCMPfGPi'],gain=params['GCMPfGPi'])
    else:
      print antagInjectionSite,": unknown antagonist experiment:",antag
  else:
    connect_pop('in','MSN','GPi', projType=params['cTypeMSNGPi'], redundancy= params['redundancyMSNGPi'], gain=params['GMSNGPi'])
    connect_pop('ex','STN','GPi', projType=params['cTypeSTNGPi'], redundancy= params['redundancySTNGPi'], gain=params['GSTNGPi'])
    if params['splitGPe']:
      connect_pop('in','Prot','GPi',   projType=params['cTypeProtGPi'], redundancy= params['redundancyProtGPi'], gain=params['GProtGPi'])
    else:
      connect_pop('in','GPe','GPi',   projType=params['cTypeGPeGPi'], redundancy= params['redundancyGPeGPi'], gain=params['GGPeGPi'])
    connect_pop('ex','CMPf','GPi',projType=params['cTypeCMPfGPi'],redundancy= params['redundancyCMPfGPi'],gain=params['GCMPfGPi'])

  base_weights = {'CSN_MSN': CSN_MSN, 'PTN_MSN': PTN_MSN, 'CMPf_MSN': CMPf_MSN}

  return base_weights

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
  if params['splitGPe']:
    GABA_afferents = ['MSN', 'Arky', 'Prot'] # afferents with gabaergic connections
  else:
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
  #nest.SetKernelStatus({"resolution": 0.005}) # simulates with a higher precision
  initNeurons()

  print '/!\ Using the following LG14 parameterization',params['LG14modelID']
  loadLG14params(params['LG14modelID'])
  loadThetaFromCustomparams(params)

  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  if params['splitGPe']:
    necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbArky','nbProt','nbGPi','nbCSN','nbPTN','nbCMPf',
                     'IeMSN','IeFSI','IeSTN','IeGPe','IeArky','IeProt','IeGPi',
                     'GCSNMSN','GPTNMSN','GCMPfMSN','GMSNMSN','GFSIMSN','GSTNMSN','GGPeMSN','GArkyMSN',
                     'GCSNFSI','GPTNFSI','GSTNFSI','GGPeFSI','GArkyFSI','GCMPfFSI','GFSIFSI',
                     'GPTNSTN','GCMPfSTN','GGPeSTN','GProtSTN',
                     'GCMPfGPe','GSTNGPe','GMSNGPe','GGPeGPe',
                     'GCMPfArky','GSTNArky','GMSNArky','GArkyArky','GProtArky',
                     'GCMPfProt','GSTNProt','GMSNProt','GProtProt','GArkyProt',
                     'GMSNGPi','GSTNGPi','GGPeGPi','GProtGPi','GCMPfGPi',
                     'redundancyCSNMSN','redundancyPTNMSN','redundancyCMPfMSN','redundancyMSNMSN','redundancyFSIMSN','redundancySTNMSN','redundancyGPeMSN','redundancyArkyMSN',
                     'redundancyCSNFSI','redundancyPTNFSI','redundancySTNFSI','redundancyGPeFSI','redundancyArkyFSI','redundancyCMPfFSI','redundancyFSIFSI',
                     'redundancyPTNSTN','redundancyCMPfSTN','redundancyGPeSTN','redundancyProtSTN',
                     'redundancyCMPfGPe','redundancySTNGPe','redundancyMSNGPe','redundancyGPeGPe',
                     'redundancyCMPfArky','redundancySTNArky','redundancyMSNArky','redundancyArkyArky','redundancyProtArky',
                     'redundancyCMPfProt','redundancySTNProt','redundancyMSNProt','redundancyProtProt','redundancyArkyProt',
                     'redundancyMSNGPi','redundancySTNGPi','redundancyGPeGPi','redundancyProtGPi','redundancyCMPfGPi',]
  else:
    necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf',
                     'IeMSN','IeFSI','IeSTN','IeGPe','IeGPi',
                     'GCSNMSN','GPTNMSN','GCMPfMSN','GMSNMSN','GFSIMSN','GSTNMSN','GGPeMSN',
                     'GCSNFSI','GPTNFSI','GSTNFSI','GGPeFSI','GCMPfFSI','GFSIFSI',
                     'GPTNSTN','GCMPfSTN','GGPeSTN',
                     'GCMPfGPe','GSTNGPe','GMSNGPe','GGPeGPe',
                     'GMSNGPi','GSTNGPi','GGPeGPi','GCMPfGPi',
                     'redundancyCSNMSN','redundancyPTNMSN','redundancyCMPfMSN','redundancyMSNMSN','redundancyFSIMSN','redundancySTNMSN','redundancyGPeMSN',
                     'redundancyCSNFSI','redundancyPTNFSI','redundancySTNFSI','redundancyGPeFSI','redundancyCMPfFSI','redundancyFSIFSI',
                     'redundancyPTNSTN','redundancyCMPfSTN','redundancyGPeSTN',
                     'redundancyCMPfGPe','redundancySTNGPe','redundancyMSNGPe','redundancyGPeGPe',
                     'redundancyMSNGPi','redundancySTNGPi','redundancyGPeGPi','redundancyCMPfGPi',]
  for np in necessaryParams:
    if np not in params:
      raise KeyError('Missing parameter: '+np)

  #------------------------
  # creation and connection of the neural populations
  #------------------------

  createBG()
  return connectBG(antagInjectionSite,antag)




