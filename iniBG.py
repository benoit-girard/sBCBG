#!/usr/bin/env python
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
  print('\nCreating neurons\n================')

  if params['nbCh'] == 1:
    # single-channel nuclei
    def create_pop(*args, **kwargs):
      if 'nbCh' in kwargs.keys():
        # remove the extra arg
        kwargs.pop("nbCh", None)
      create(*args, **kwargs)
    update_Ie = lambda p: nest.SetStatus(Pop[p],{"I_e":params['Ie'+p]})
  elif 'topo' not in params.keys() or params['topo'] == False:
    # multi-channel nuclei
    def create_pop(*args, **kwargs):
      if 'nbCh' not in kwargs.keys():
        # enforce the default
        kwargs['nbCh'] = params['nbCh']
      createMC(*args, **kwargs)
    update_Ie = lambda p: [nest.SetStatus(Pop[p][i],{"I_e":params['Ie'+p]}) for i in range(len(Pop[p]))]
  else:
    # topological nuclei
    def create_pop(*args, **kwargs):
      if 'nbCh' not in kwargs.keys():
        # enforce the default
        kwargs['nbCh'] = params['nbCh']
      createTopoMC(*args, layout=params['topo'], c=params['channel_center'], r=params['channel_radius'], **kwargs)
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

  print("Number of simulated neurons:"+str(nbSim))

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

  print("Gains on LG14 syn. strength:"+str(G))

  if params['nbCh'] == 1:
    # single-channel nuclei
    connect_pop = lambda *args, **kwargs: connect(*args, RedundancyType=params['RedundancyType'], stochastic_delays=params['stochastic_delays'], **kwargs)
  elif 'topo' not in params.keys() or params['topo'] == False:
    # multi-channel nuclei
    def connect_pop(*args, **kwargs):
      if 'source_channels' not in kwargs.keys():
        # enforce the default
        kwargs['source_channels'] = range(params['nbCh'])
      return connectMC(*args, RedundancyType=params['RedundancyType'], stochastic_delays=params['stochastic_delays'], **kwargs)
  else:
    # topological nuclei
    def connect_pop(*args, **kwargs):
      if 'source_channels' not in kwargs.keys():
        # enforce the default
        kwargs['source_channels'] = range(params['nbCh'])
      topo_spreads = [params['spread_focused'], params['spread_diffuse']]
      return connectTopoMC(*args, RedundancyType=params['RedundancyType'], stochastic_delays=params['stochastic_delays'], spreads=topo_spreads, **kwargs)

  #-------------------------
  # connection of populations
  #-------------------------
  print('\nConnecting neurons\n================')
  print("** "+antag+" antagonist injection in "+antagInjectionSite+" **")
  print('* MSN Inputs')
  if 'nbCues' not in params.keys():
    # usual case: CSN have as the same number of channels than the BG nuclei
    CSN_MSN = connect_pop('ex','CSN','MSN', projType=params['cTypeCSNMSN'], redundancy=params['redundancyCSNMSN'], gain=G['MSN'])
  else:
    # special case: extra 'cue' channels that target MSN
    CSN_MSN = connect_pop('ex','CSN','MSN', projType=params['cTypeCSNMSN'], redundancy=params['redundancyCSNMSN'], gain=G['MSN']/2., source_channels=range(params['nbCh']))
    connect_pop('ex','CSN','MSN', projType='diffuse', redundancy=params['redundancyCSNMSN'], gain=G['MSN']/2., source_channels=range(params['nbCh'], params['nbCh']+params['nbCues']))
  PTN_MSN = connect_pop('ex','PTN','MSN', projType=params['cTypePTNMSN'], redundancy= params['redundancyPTNMSN'], gain=G['MSN'])
  CMPf_MSN = connect_pop('ex','CMPf','MSN',projType=params['cTypeCMPfMSN'],redundancy= params['redundancyCMPfMSN'],gain=G['MSN'])
  connect_pop('in','MSN','MSN', projType=params['cTypeMSNMSN'], redundancy= params['redundancyMSNMSN'], gain=G['MSN'])
  connect_pop('in','FSI','MSN', projType=params['cTypeFSIMSN'], redundancy= params['redundancyFSIMSN'], gain=G['MSN'])
  # some parameterizations from LG14 have no STN->MSN or GPe->MSN synaptic contacts
  if alpha['STN->MSN'] != 0:
    print("Note: alpha['STN->MSN'] = "+str(alpha['STN->MSN']))
    connect_pop('ex','STN','MSN', projType=params['cTypeSTNMSN'], redundancy= params['redundancySTNMSN'], gain=G['MSN'])
  if alpha['GPe->MSN'] != 0:
    if params['splitGPe']:
      print("Note: alpha['Arky->MSN']"+str(alpha['Arky->MSN']))
      connect_pop('in','Arky','MSN', projType=params['cTypeArkyMSN'], redundancy=params['redundancyArkyMSN'], gain=params['GArkyMSN'])
    else:
      print("Note: alpha['GPe->MSN'] = "+str(alpha['GPe->MSN']))
      connect_pop('in','GPe','MSN', projType=params['cTypeGPeMSN'], redundancy= params['redundancyGPeMSN'], gain=G['MSN'])

  print('* FSI Inputs')
  connect_pop('ex','CSN','FSI', projType=params['cTypeCSNFSI'], redundancy= params['redundancyCSNFSI'], gain=G['FSI'])
  connect_pop('ex','PTN','FSI', projType=params['cTypePTNFSI'], redundancy= params['redundancyPTNFSI'], gain=G['FSI'])
  if alpha['STN->FSI'] != 0:
      connect_pop('ex','STN','FSI', projType=params['cTypeSTNFSI'],redundancy= params['redundancySTNFSI'],gain=G['FSI'])
  if params['splitGPe']:
      connect_pop('in','Arky','FSI', projType=params['cTypeArkyFSI'], redundancy= params['redundancyArkyFSI'], gain=params['GArkyFSI'])
  else:
      connect_pop('in','GPe','FSI', projType=params['cTypeGPeFSI'], redundancy= params['redundancyGPeFSI'], gain=G['FSI'])
  connect_pop('ex','CMPf','FSI',projType=params['cTypeCMPfFSI'],redundancy= params['redundancyCMPfFSI'],gain=G['FSI'])
  connect_pop('in','FSI','FSI', projType=params['cTypeFSIFSI'], redundancy= params['redundancyFSIFSI'], gain=G['FSI'])

  print('* STN Inputs')
  connect_pop('ex','PTN','STN', projType=params['cTypePTNSTN'], redundancy= params['redundancyPTNSTN'],  gain=G['STN'])
  connect_pop('ex','CMPf','STN',projType=params['cTypeCMPfSTN'],redundancy= params['redundancyCMPfSTN'], gain=G['STN'])
  if params['splitGPe']:
      connect_pop('in','Prot','STN', projType=params['cTypeProtSTN'], redundancy= params['redundancyProtSTN'],  gain=params['GProtSTN'])
  else:
      connect_pop('in','GPe','STN', projType=params['cTypeGPeSTN'], redundancy= params['redundancyGPeSTN'],  gain=G['STN'])

  print('* GPe Inputs')
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
         if antag == 'AMPA':
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
         if antag == 'AMPA':
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
             connect_pop('NMDA','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=G['GPe'])
             connect_pop('NMDA','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=G['GPe'])
             connect_pop('in','MSN','GPe',   projType=params['cTypeMSNGPe'], redundancy= params['redundancyMSNGPe'], gain=G['GPe'])
             connect_pop('in', GPe_recurrent_source, 'GPe', projType=params['cTypeGPeGPe'], redundancy= params['redundancyGPeGPe'], gain=G['GPe'])
         elif antag == 'NMDA':
             connect_pop('AMPA','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=G['GPe'])
             connect_pop('AMPA','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=G['GPe'])
             connect_pop('in','MSN','GPe',   projType=params['cTypeMSNGPe'], redundancy= params['redundancyMSNGPe'], gain=G['GPe'])
             connect_pop('in', GPe_recurrent_source, 'GPe', projType=params['cTypeGPeGPe'], redundancy= params['redundancyGPeGPe'], gain=G['GPe'])
         elif antag == 'AMPA+GABAA':
             connect_pop('NMDA','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=G['GPe'])
             connect_pop('NMDA','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=G['GPe'])
         elif antag == 'GABAA':
             connect_pop('ex','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=G['GPe'])
             connect_pop('ex','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=G['GPe'])
         else:
             print(antagInjectionSite+": unknown antagonist experiment: "+antag)
      else:
         connect_pop('ex','CMPf','GPe',projType=params['cTypeCMPfGPe'],redundancy= params['redundancyCMPfGPe'],gain=G['GPe'])
         connect_pop('ex','STN','GPe', projType=params['cTypeSTNGPe'], redundancy= params['redundancySTNGPe'], gain=G['GPe'])
         connect_pop('in','MSN','GPe', projType=params['cTypeMSNGPe'], redundancy= params['redundancyMSNGPe'], gain=G['GPe'])
         connect_pop('in', GPe_recurrent_source, 'GPe', projType=params['cTypeGPeGPe'], redundancy= params['redundancyGPeGPe'], gain=G['GPe'])

  print('* GPi Inputs')
  if antagInjectionSite =='GPi':
    if   antag == 'AMPA+NMDA+GABAA':
      pass
    elif antag == 'NMDA':
      connect_pop('in','MSN','GPi',   projType=params['cTypeMSNGPi'], redundancy= params['redundancyMSNGPi'], gain=G['GPi'])
      connect_pop('AMPA','STN','GPi', projType=params['cTypeSTNGPi'], redundancy= params['redundancySTNGPi'], gain=G['GPi'])
      if params['splitGPe']:
          connect_pop('in','Prot','GPi',   projType=params['cTypeProtGPi'], redundancy= params['redundancyProtGPi'], gain=params['GProtGPi'])
      else:
          connect_pop('in','GPe','GPi',   projType=params['cTypeGPeGPi'], redundancy= params['redundancyGPeGPi'], gain=G['GPi'])
      connect_pop('AMPA','CMPf','GPi',projType=params['cTypeCMPfGPi'],redundancy= params['redundancyCMPfGPi'],gain=G['GPi'])
    elif antag == 'NMDA+AMPA':
      connect_pop('in','MSN','GPi', projType=params['cTypeMSNGPi'],redundancy= params['redundancyMSNGPi'], gain=G['GPi'])
      if params['splitGPe']:
          connect_pop('in','Prot','GPi',   projType=params['cTypeProtGPi'], redundancy= params['redundancyProtGPi'], gain=params['GProtGPi'])
      else:
          connect_pop('in','GPe','GPi', projType=params['cTypeGPeGPi'],redundancy= params['redundancyGPeGPi'], gain=G['GPi'])
    elif antag == 'AMPA':
      connect_pop('in','MSN','GPi',   projType=params['cTypeMSNGPi'], redundancy= params['redundancyMSNGPi'], gain=G['GPi'])
      connect_pop('NMDA','STN','GPi', projType=params['cTypeSTNGPi'], redundancy= params['redundancySTNGPi'], gain=G['GPi'])
      if params['splitGPe']:
          connect_pop('in','Prot','GPi',   projType=params['cTypeProtGPi'], redundancy= params['redundancyProtGPi'], gain=params['GProtGPi'])
      else:
          connect_pop('in','GPe','GPi',   projType=params['cTypeGPeGPi'], redundancy= params['redundancyGPeGPi'], gain=G['GPi'])
      connect_pop('NMDA','CMPf','GPi',projType=params['cTypeCMPfGPi'],redundancy= params['redundancyCMPfGPi'],gain=G['GPi'])
    elif antag == 'GABAA':
      connect_pop('ex','STN','GPi', projType=params['cTypeSTNGPi'], redundancy= params['redundancySTNGPi'], gain=G['GPi'])
      connect_pop('ex','CMPf','GPi',projType=params['cTypeCMPfGPi'],redundancy= params['redundancyCMPfGPi'],gain=G['GPi'])
    else:
      print(antagInjectionSite+": unknown antagonist experiment: "+antag)
  else:
     connect_pop('in','MSN','GPi', projType=params['cTypeMSNGPi'], redundancy= params['redundancyMSNGPi'], gain=G['GPi'])
     connect_pop('ex','STN','GPi', projType=params['cTypeSTNGPi'], redundancy= params['redundancySTNGPi'], gain=G['GPi'])
     if params['splitGPe']:
         connect_pop('in','Prot','GPi', projType=params['cTypeProtGPi'], redundancy= params['redundancyProtGPi'], gain=params['GProtGPi'])
     else:
         connect_pop('in','GPe','GPi', projType=params['cTypeGPeGPi'], redundancy= params['redundancyGPeGPi'], gain=G['GPi'])
     connect_pop('ex','CMPf','GPi',projType=params['cTypeCMPfGPi'],redundancy= params['redundancyCMPfGPi'],gain=G['GPi'])

  base_weights = {'CSN_MSN': CSN_MSN, 'PTN_MSN': PTN_MSN, 'CMPf_MSN': CMPf_MSN}

  return base_weights


#------------------------------------------
# Creates a "fake" copy of a nucleus with the same number of neurons and which fires with Poisson spike trains of a specified frequency
# /original_nuc/ is the name of the original nucleus to copy, /new_nuc/ is the desired name of the Poisson nucleus to create
# /poisson_rate/ specifies the desired firing frequency
#------------------------------------------
def Poisson_copy(original_nuc, new_nuc, poisson_rate):
  rate[new_nuc] = poisson_rate
  for nucleus_dict in [nbSim, neuronCounts]:
    nucleus_dict[new_nuc] = nucleus_dict[original_nuc]
  if params['nbCh'] == 1:
    create(new_nuc, fake=True, parrot=True)
  else:
    createMC(new_nuc, params['nbCh'], fake=True, parrot=True)


#------------------------------------------
# Replaces the specified connection in the normal BG circuitry by a "bypass" from another nucleus
# The arguments specify respectively the previous source nucleus (/original_nuc_from/ to be disconnected), the new source nucleus (/new_nuc_from/ to be connected instead), with respect to the target nucleus /nuc_to/
# Note: this is really useful in conjonction with a fake nucleus created through a call to copy_nucleus_Poisson
#------------------------------------------
def bypass_connection(original_nuc_from, new_nuc_from, nuc_to):
  # retrieve the existing connections
  existing_conns = nest.GetConnections(source=np.array(Pop[original_nuc_from]).flatten().tolist(), target=np.array(Pop[nuc_to]).flatten().tolist())
  if len(existing_conns) == 0:
    print('skipping non-existent connection '+original_nuc_from+'->'+nuc_to)
  else:
    print("Replacing connection "+original_nuc_from+"->"+nuc_to+" by "+new_nuc_from+"->"+nuc_to+"...")
    # retrieving the receptor, weight and delay parameters of the current connection
    weights_delays_recs = np.array(nest.GetStatus(existing_conns, keys=['weight', 'delay', 'receptor'])).transpose()
    old_src, tgt, _, _, _ = zip(*existing_conns)
    # get input neurons from the fake nucleus that match the current source nucleus neurons
    new_src = np.array(sorted(old_src, key=lambda x: x)) - min(old_src) + Pop[new_nuc_from][0][0] # remark: this relies on sequential ordering of pynest created neurons
    # destroy the current connection by setting its weight to 0
    nest.SetStatus(existing_conns, [{'weight': 0.} for i in range(len(existing_conns))])
    # connect the fake nucleus instead, looping over the required receptors
    for r in np.unique(weights_delays_recs[2]):
      # get all connection indices using this receptor
      same_r = np.where(weights_delays_recs[2] == r)[0]
      # creates the new connection for these indices
      nest.Connect(pre=new_src[same_r].tolist(),
                   post=np.array(tgt)[same_r].tolist(),
                   conn_spec='one_to_one',
                   syn_spec={'model': 'static_synapse',
                             #'model': 'static_synapse_lbl', 'synapse_label': 0,
                             'receptor_type': int(r),
                             'weight': weights_delays_recs[0][same_r].tolist(),
                             'delay': weights_delays_recs[1][same_r].tolist()})

#------------------------------------------
# Re-weight a specific connection, characterized by a source, a target, and a receptor
# Returns the previous value of that connection (useful for 'reactivating' after a deactivation experiment)
#    _
#   / \
#  / | \    Careful: this routine seems to hang forever with nest-5g (development version built on July 20, 2018)
# /  o  \
# -------
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

  print('/!\ Using the following LG14 parameterization'+str(params['LG14modelID']))
  loadLG14params(params['LG14modelID'])
  #loadThetaFromCustomparams(params)

  # We check that all the necessary parameters have been defined. They should be in the modelParams.py file.
  # If one of them misses, we exit the program.
  if params['splitGPe']:
    necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbArky','nbProt','nbGPi','nbCSN','nbPTN','nbCMPf','IeMSN','IeFSI','IeSTN','IeArky','IeProt','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','GArkyMSN','GArkyFSI','GProtSTN','GCMPfArky','GMSNArky','GCMPfProt','GMSNProt','GSTNProt','GSTNArky','GProtGPi','GArkyArky','GProtArky','GArkyProt','GProtProt','redundancyCSNMSN','redundancyPTNMSN','redundancyCMPfMSN','redundancyMSNMSN','redundancyFSIMSN','redundancySTNMSN','redundancyArkyMSN','redundancyCSNFSI','redundancyPTNFSI','redundancySTNFSI','redundancyArkyFSI','redundancyCMPfFSI','redundancyFSIFSI','redundancyPTNSTN','redundancyCMPfSTN','redundancyProtSTN','redundancyCMPfArky','redundancySTNArky','redundancyMSNArky','redundancyArkyArky','redundancyProtArky','redundancyCMPfProt','redundancySTNProt','redundancyMSNProt','redundancyArkyProt','redundancyProtProt','redundancyMSNGPi','redundancySTNGPi','redundancyProtGPi','redundancyCMPfGPi',]
  else:
    necessaryParams=['nbCh','nbMSN','nbFSI','nbSTN','nbGPe','nbGPi','nbCSN','nbPTN','nbCMPf','IeMSN','IeFSI','IeSTN','IeGPe','IeGPi','GMSN','GFSI','GSTN','GGPe','GGPi','redundancyCSNMSN','redundancyPTNMSN','redundancyCMPfMSN','redundancyMSNMSN','redundancyFSIMSN','redundancySTNMSN','redundancyGPeMSN','redundancyCSNFSI','redundancyPTNFSI','redundancySTNFSI','redundancyGPeFSI','redundancyCMPfFSI','redundancyFSIFSI','redundancyPTNSTN','redundancyCMPfSTN','redundancyGPeSTN','redundancyCMPfGPe','redundancySTNGPe','redundancyMSNGPe','redundancyGPeGPe','redundancyMSNGPi','redundancySTNGPi','redundancyGPeGPi','redundancyCMPfGPi',]
  for np in necessaryParams:
    if np not in params:
      raise KeyError('Missing parameter: '+np)

  #------------------------
  # creation and connection of the neural populations
  #------------------------

  createBG()
  return connectBG(antagInjectionSite,antag)
