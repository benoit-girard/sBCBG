# Base parameters and their meaning
# these defaults can be overrided via a custom python parameter file or via the commandline (in this order: commandline arguments take precedence over customParams.py, which take precendence over the defaults defined here)

####################################################################################
# This file should not be modified! Use commandline arguments or custom param file #
####################################################################################

params = {
'splitGPe':                  False, # if True, the GPe will be splitted in prototypical and arkypallidal neurons
'nbCh':                          1, # number of concurrent channels to simulate
'LG14modelID':                   9, # LG 2014 parameterization used (default: model #9)
'nbMSN':                     2644., # population size (default: 1/1000 of the BG)
'nbFSI':                       53., # ^
'nbSTN':                        8., # ^
'nbGPe':                       25., # ^
'nbArky':                       5., # part of the GPe which projects to the striatum
'nbProt':                      20., # part of the GPe which projects to the STN and GPi/SNr
'nbGPi':                       14., # ^
'nbCSN':                     3000., # ^
'nbPTN':                      100., # ^
'nbCMPf':                       9., # ^
'GMSN':                         1., # gain on all synapses (default: no gain)
'GFSI':                         1., # ^
'GSTN':                         1., # ^
'GGPe':                         1., # ^
'GGPi':                         1., # ^

'GArkyMSN':                     1., # the following gains can be used to test various models of GPe splitting
'GArkyFSI':                     1., # ^
'GProtSTN':                     1., # ^
'GCMPfArky':                    1., # ^
'GMSNArky':                     1., # ^
'GCMPfProt':                    1., # ^
'GMSNProt':                     1., # ^
'GSTNProt':                     1., # ^
'GSTNArky':                     1., # ^
'GProtGPi':                     1., # ^
'GArkyArky':                  1/5., # ^
'GProtArky':                  4/5., # ^
'GArkyProt':                  1/5., # ^
'GProtProt':                  4/5., # ^

'IeMSN':                        0., # tonic input currents (default: no input current)
'IeFSI':                        0., # ^
'IeSTN':                        0., # ^
'IeGPe':                        0., # ^
'IeArky':                       0., # ^
'IeProt':                       0., # ^
'IeGPi':                        0., # ^
# There are 3 different format for setting the inDegreeXY (=number of different incoming neurons from X that target one neuron in Y)
# - 'inDegreeAbs': specify the absolute number of different input neurons from X that target one neuron of Y --- be careful when using this setting, as the population size varies widly between the striatum and the downstream nuclei
# - 'outDegreeAbs': specify the absolute number of contacts between an axon from X to each dendritic tree in Y
# - 'outDegreeCons': specify the constrained number of contacts between an axon from X to each dendritic tree in Y as a fraction between 0 (minimal number of contacts to achieve required axonal bouton counts) and 1 (maximal number of contacts with respect to population numbers)
'RedundancyType':   'outDegreeAbs', # by default all axons are hypothesized to target each dendritic tree at 3 different locations
'redundancyCSNMSN':              3, # ^
'redundancyPTNMSN':              3, # ^
'redundancyCMPfMSN':             3, # ^
'redundancyMSNMSN':              3, # ^
'redundancyFSIMSN':              3, # ^
'redundancySTNMSN':              3, # ^
'redundancyGPeMSN':              3, # ^
'redundancyArkyMSN':             3, # ^
'redundancyCSNFSI':              3, # ^
'redundancyPTNFSI':              3, # ^
'redundancySTNFSI':              3, # ^
'redundancyGPeFSI':              3, # ^
'redundancyArkyFSI':             3, # ^
'redundancyCMPfFSI':             3, # ^
'redundancyFSIFSI':              3, # ^
'redundancyPTNSTN':              3, # ^
'redundancyCMPfSTN':             3, # ^
'redundancyGPeSTN':              3, # ^
'redundancyProtSTN':             3, # ^
'redundancyCMPfGPe':             3, # ^
'redundancySTNGPe':              3, # ^
'redundancyMSNGPe':              3, # ^
'redundancyGPeGPe':              3, # ^

'redundancyCMPfArky':            3, # ^
'redundancySTNArky':             3, # ^
'redundancyMSNArky':             3, # ^
'redundancyArkyArky':            3, # ^
'redundancyProtArky':            3, # ^

'redundancyCMPfProt':            3, # ^
'redundancySTNProt':             3, # ^
'redundancyMSNProt':             3, # ^
'redundancyArkyProt':            3, # ^
'redundancyProtProt':            3, # ^

'redundancyMSNGPi':              3, # ^
'redundancySTNGPi':              3, # ^
'redundancyGPeGPi':              3, # ^
'redundancyProtGPi':             3, # ^
'redundancyCMPfGPi':             3, # ^
'cTypeCSNMSN':           'focused', # defining connection types for channel-based models (focused or diffuse) based on LG14 - refer to this paper for justification
'cTypePTNMSN':           'focused', # ^
'cTypeCMPfMSN':          'diffuse', # ^
'cTypeMSNMSN':           'diffuse', # ^
'cTypeFSIMSN':           'diffuse', # ^
'cTypeSTNMSN':           'diffuse', # ^
'cTypeGPeMSN':           'diffuse', # ^
'cTypeCSNFSI':           'focused', # ^
'cTypePTNFSI':           'focused', # ^
'cTypeCMPfFSI':          'diffuse', # ^
'cTypeFSIFSI':           'diffuse', # ^
'cTypeSTNFSI':           'diffuse', # ^
'cTypeGPeFSI':           'diffuse', # ^
'cTypePTNSTN':           'focused', # ^
'cTypeCMPfSTN':          'diffuse', # ^
'cTypeGPeSTN':           'focused', # ^
'cTypeCMPfGPe':          'diffuse', # ^
'cTypeMSNGPe':           'focused', # ^
'cTypeSTNGPe':           'diffuse', # ^
'cTypeGPeGPe':           'diffuse', # ^
'cTypeCMPfGPi':          'diffuse', # ^
'cTypeMSNGPi':           'focused', # ^
'cTypeSTNGPi':           'diffuse', # ^
'cTypeGPeGPi':           'diffuse', # LG14: no data available to decide; setting to diffuse improve selection properties
'cTypeArkyMSN':          'diffuse', # ^
'cTypeArkyFSI':          'diffuse', # ^
'cTypeProtSTN':          'focused', # ^
'cTypeCMPfArky':         'diffuse', # ^
'cTypeMSNArky':          'focused', # ^
'cTypeSTNArky':          'diffuse', # ^
'cTypeArkyArky':         'diffuse', # ^
'cTypeProtArky':         'diffuse', # ^
'cTypeCMPfProt':         'diffuse', # ^
'cTypeMSNProt':          'focused', # ^
'cTypeSTNProt':          'diffuse', # ^
'cTypeArkyProt':         'diffuse', # ^
'cTypeProtProt':         'diffuse', # ^
'cTypeProtGPi':          'diffuse', # ^
'topo':                      False, # Enable the topological version of the BG? not by default
'spread_focused':              0.1, # connection spread of focused projections (used if topo==True)
'spread_diffuse':               4., # connection spread of diffuse projections (used if topo==True)
'channel_center':              0.3, # distance of channels to the center (used if topo==True)
'channel_radius':             0.25, # radius of channels (used if topo==True)
'parrotCMPf' :                True, # Should the CMPf be simulated using parrot neurons?
'stochastic_delays':          None, # If specified, gives the relative sd of a clipped Gaussian distribution for the delays
# For convenience, a few simulator variables are also set here
'whichTest':    'testPlausibility', # task to be run (default: test the plausibility through deactivation simulations)
'offsetDuration':             1000, # non-recorded stabilization period, in ms
'simDuration':                1000, # simulation duration, in ms
'nestSeed':                     20, # nest seed (affects input poisson spike trains)
'pythonSeed':                   10, # python seed (affects connection map)
'nbcpu':                         1, # number of CPUs to be used by nest
'durationH':                  '08', # max duration of a simulation, used by Sango & ISIR clusters
'nbnodes':                     '1', # number of nodes, used by K computer
}
