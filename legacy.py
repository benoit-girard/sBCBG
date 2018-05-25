#!/usr/bin/python
# -*- coding: utf-8 -*-

###################
# ABOUT THIS FILE #
###################

# legacy.py was created as an effort to keep the codebase clean while not throwing away code that might be useful again.
# The functions contained here should not be used directly, but should rather be incorporated into the main codebase


import nest
import numpy as np
import numpy.random as rnd

#-------------------------------------------------------------------------------
# Simple queue mechanism
# This function makes connections by chunk, until the chunk size is too small
#-------------------------------------------------------------------------------
def empty_queue_legacy(l, W, delay, lRecType, n, do_empty = False):
  for connect_chunk in [[l[0][i:i + n], l[1][i:i + n]] for i in xrange(0, len(l[0]), n)]:
    if len(connect_chunk[0]) < n and not do_empty:
      # return the incomplete chunk, it will processed next time
      return connect_chunk
    else:
      # call the connect function
      for r in lRecType:
        w = W[r]
        nest.Connect(pre=tuple(connect_chunk[1]), post=tuple(connect_chunk[0]), conn_spec='one_to_one', syn_spec={'receptor_type':recType[r],'weight':w,'delay':delay})
  return np.array([[],[]], dtype=int)


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
def connect_legacy(type,nameSrc,nameTgt,redundancy,RedundancyType,LCGDelays=True,gain=1., verbose=True, projType=''):

  def printv(text):
    if verbose:
      print(text)

  if RedundancyType == 'inDegreeAbs':
    printv("not implemented")
    return
  elif RedundancyType == 'outDegreeAbs':
    inDegree = get_frac(1./redundancy, nameSrc, nameTgt, verbose=False)
  elif RedundancyType == 'outDegreeCons':
    #### fractional inDegree is expressed as a fraction of max number of neurons
    ###printv('Using fractional inDegree (value supplied: '+str(inDegree)+')')
    ###inDegree = get_frac(inDegree, nameSrc, nameTgt, verbose=False)
    printv("not implemented")
    return

  # check if in degree acceptable (not larger than number of neurons in the source nucleus)
  if inDegree  > nbSim[nameSrc]:
    printv("/!\ WARNING: required 'in degree' ("+str(inDegree)+") larger than number of neurons in the source population ("+str(nbSim[nameSrc])+"), thus reduced to the latter value")
    inDegree = nbSim[nameSrc]

  if inDegree == 0.:
    printv("/!\ WARNING: non-existent connection strength, will skip")
    return

  printv("* connecting "+nameSrc+" -> "+nameTgt+" with "+type+" connection and "+str(inDegree)+ " inputs")

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
    raise KeyError('Undefined connexion type: '+type)
  
  W = computeW(lRecType,nameSrc,nameTgt,inDegree,gain,verbose=verbose)

  if nameSrc+'->'+nameTgt in ConnectMap:
    loadConnectMap = True
  else:
    loadConnectMap = False
    ConnectMap[nameSrc+'->'+nameTgt] = []

  # determine which transmission delay to use:
  if LCGDelays:
    delay= tau[nameSrc+'->'+nameTgt]
  else:
    delay= 1.

  pop_array = np.array(Pop[nameSrc]) # convert to numpy array to allow indexing
  connect_queue = np.array([[],[]], dtype=int) # the connection queue
  max_chunk_size = 10000 # what is the max number of connections to add simultaneously?

  # To ensure that for excitatory connections, Tgt neurons receive AMPA and NMDA projections from the same Src neurons, 
  # we have to handle the "indegree" connectivity ourselves:
  for nTgt in range(int(nbSim[nameTgt])):
    #node_info   = nest.GetStatus([Pop[nameTgt][nTgt]])[0] # use the process-specific random number generator
    #nrnd = nstrand.pyRngs[node_info['vp']]                # ^
    nrnd = nstrand.pyMasterRng # use the master python seed
    if loadConnectMap:
      # use previously created connectivity map
      inputTable = ConnectMap[nameSrc+'->'+nameTgt][nTgt]
      inDeg = len(inputTable)
    else:
      # if no connectivity map exists between the two populations, let's create one
      r = inDegree - int(inDegree)
      inDeg = int(inDegree) if nrnd.rand() > r else int(inDegree)+1
      inputTable = pop_array[nrnd.choice(int(nbSim[nameSrc]), size=inDeg, replace=False)]
      ConnectMap[nameSrc+'->'+nameTgt] += [tuple(inputTable)]

    connect_queue = np.concatenate((connect_queue, np.array([[Pop[nameTgt][nTgt]]*inDeg, inputTable])), axis=1)

    if connect_queue.shape[1] > max_chunk_size:
      # space for at least one chunk? => empty the queue
      connect_queue = empty_queue_legacy(connect_queue, W, delay, lRecType, max_chunk_size, do_empty=False)

    # old way of connecting neurons
    #for r in lRecType:
    #  w = W[r]
    #  nest.Connect(pre=ConnectMap[nameSrc+'->'+nameTgt][nTgt], post=(Pop[nameTgt][nTgt],), syn_spec={'receptor_type':recType[r],'weight':w,'delay':delay})

  empty_queue_legacy(connect_queue, W, delay, lRecType, max_chunk_size, do_empty=True)


#-------------------------------------------------------------------------------
# Establishes a connexion between two populations, following the results of LG14, in a MultiChannel context
# type : a string 'ex' or 'in', defining whether it is excitatory or inhibitory
# nameTgt, nameSrc : strings naming the populations, as defined in NUCLEI list
# projType : type of projections. For the moment: 'focused' (only channel-to-channel connection) and 
#            'diffuse' (all-to-one with uniform distribution)
# inDegree : number of neurons from Src project to a single Tgt neuron
# LCGDelays: shall we use the delays obtained by (Liénard, Cos, Girard, in prep) or not (default = True)
# gain : allows to amplify the weight normally deduced from LG14
#-------------------------------------------------------------------------------
def connectMC_legacy(type,nameSrc,nameTgt,projType,redundancy,RedundancyType,LCGDelays=True,gain=1., verbose=True):

  def printv(text):
    if verbose:
      print(text)


  if RedundancyType == 'inDegreeAbs':
    printv("not implemented")
    return
  elif RedundancyType == 'outDegreeAbs':
    inDegree = get_frac(1./redundancy, nameSrc, nameTgt, verbose=False)
  elif RedundancyType == 'outDegreeCons':
    #### fractional inDegree is expressed as a fraction of max number of neurons
    ###printv('Using fractional inDegree (value supplied: '+str(inDegree)+')')
    ###inDegree = get_frac(inDegree, nameSrc, nameTgt, verbose=False)
    printv("not implemented")
    return

  printv("* connecting "+nameSrc+" -> "+nameTgt+" with "+projType+" "+type+" connection and "+str(inDegree)+" inputs")

  # check if in degree acceptable (not larger than number of neurons in the source nucleus)
  if projType == 'focused' and inDegree > nbSim[nameSrc]:
    printv("/!\ WARNING: required 'in degree' ("+str(inDegree)+") larger than number of neurons in the source population ("+str(nbSim[nameSrc])+"), thus reduced to the latter value")
    inDegree = nbSim[nameSrc]
  if projType == 'diffuse' and inDegree  > nbSim[nameSrc]*len(Pop[nameSrc]):
    printv("/!\ WARNING: required 'in degree' ("+str(inDegree)+") larger than number of neurons in the source population ("+str(nbSim[nameSrc]*len(Pop[nameSrc]))+"), thus reduced to the latter value")
    inDegree = nbSim[nameSrc]*len(Pop[nameSrc])

  if inDegree == 0.:
    printv("/!\ WARNING: non-existent connection strength, will skip")
    return

  # prepare receptor type lists:
  if type == 'ex':
    lRecType = ['AMPA','NMDA']
  elif type == 'AMPA':
    lRecType = ['AMPA']
  elif type == 'NMDA':
    lRecType = ['NMDA']
  elif type == 'in':
    lRecType = ['GABA']
  else:
    raise KeyError('Undefined connexion type: '+type)

  # compute the global weight of the connection, for each receptor type:
  W = computeW(lRecType,nameSrc,nameTgt,inDegree,gain,verbose=verbose)

  # check whether a connection map has already been drawn or not:
  if nameSrc+'->'+nameTgt in ConnectMap:
    #print "Using existing connection map"
    loadConnectMap = True
  else:
    #print "Will create a connection map"
    loadConnectMap = False
    ConnectMap[nameSrc+'->'+nameTgt] = [[] for i in range(len(Pop[nameTgt]))]

  # determine which transmission delay to use:
  if LCGDelays:
    delay = tau[nameSrc+'->'+nameTgt]
  else:
    delay = 1.

  # To ensure that for excitatory connections, Tgt neurons receive AMPA and NMDA projections from the same Src neurons,
  # we have to handle the "indegree" connectivity ourselves:
  for tgtChannel in range(len(Pop[nameTgt])): # for each channel of the Target nucleus
    for nTgt in range(int(nbSim[nameTgt])): # for each neuron in this channel 
      if not loadConnectMap:
      # if no connectivity map exists between the two populations, let's create one
        if projType =='focused': # if projections focused, input come only from the same channel as tgtChannel
          r = inDegree - int(inDegree)
          inputTable = rnd.choice(int(nbSim[nameSrc]),size=int(inDegree) if rnd.rand() > r else int(inDegree)+1,replace=False)
          inputPop = []
          for i in inputTable:
            inputPop.append(Pop[nameSrc][tgtChannel][i])
          inputPop = tuple(inputPop)

          ConnectMap[nameSrc+'->'+nameTgt][tgtChannel].append(inputPop)
        elif projType=='diffuse': # if projections diffused, input connections are shared among each possible input channel equally
          n = int(inDegree)/int(len(Pop[nameSrc]))
          r = float(inDegree)/float(len(Pop[nameSrc])) - n
          inputPop = []
          #print nameSrc,'->',nameTgt,'#input connections:',n,'(',r,')'
          for srcChannel in range(len(Pop[nameSrc])):
            if rnd.rand() < r:
              nbInPerChannel = n + 1
            else:
              nbInPerChannel = n
            #print '   ',nbInPerChannel
            inputTable = rnd.choice(int(nbSim[nameSrc]),size=nbInPerChannel,replace=False)
            for i in inputTable:
              inputPop.append(Pop[nameSrc][srcChannel][i])

          inputPop = tuple(inputPop)
          ConnectMap[nameSrc+'->'+nameTgt][tgtChannel].append(inputPop)
        else:
          print "Unknown multiple channel connection method",projType
      else:
      #otherwise, use the existing one
        #print nameSrc,"->",nameTgt,"using previously defined connection map"
        inputPop = ConnectMap[nameSrc+'->'+nameTgt][tgtChannel][nTgt]

      for r in lRecType:
        w = W[r]

        nest.Connect(pre=inputPop, post=(Pop[nameTgt][tgtChannel][nTgt],),syn_spec={'receptor_type':recType[r],'weight':w,'delay':delay})
