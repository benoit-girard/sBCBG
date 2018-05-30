# Base parameters and their meaning
# these defaults can be overrided via a custom python parameter file or via the commandline (in this order: commandline arguments take precedence over customParams.py, which take precendence over the defaults defined here)

####################################################################################
# This file should not be modified! Use commandline arguments or custom param file #
####################################################################################

params = {
'splitGPe':                  False,
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

'GCSNMSN':           1., # defining connection types for channel-based models (focused or diffuse) based on LG14 - refer to this paper for justification
'GPTNMSN':           1., # ^
'GCMPfMSN':          1., # ^
'GMSNMSN':           1., # ^
'GFSIMSN':           1., # ^
'GSTNMSN':           1., # ^
'GGPeMSN':           1., # ^
'GArkyMSN':          1., # ^

'GCSNFSI':           1., # ^
'GPTNFSI':           1., # ^
'GCMPfFSI':          1., # ^
'GFSIFSI':           1., # ^
'GSTNFSI':           1., # ^
'GGPeFSI':           1., # ^
'GArkyFSI':          1., # ^

'GPTNSTN':           1., # ^
'GCMPfSTN':          1., # ^
'GGPeSTN':           1., # ^
'GProtSTN':          1., # ^

'GCMPfGPe':          1., # ^
'GMSNGPe':           1., # ^
'GSTNGPe':           1., # ^
'GGPeGPe':           1., # ^

'GCMPfArky':         1., # ^
'GMSNArky':          1., # ^
'GSTNArky':          1., # ^
'GArkyArky':         1/5., # ^
'GProtArky':         4/5., # ^

'GCMPfProt':         1., # ^
'GMSNProt':          1., # ^
'GSTNProt':          1., # ^
'GArkyProt':         1/5., # ^
'GProtProt':         4/5., # ^

'GCMPfGPi':          1., # ^
'GMSNGPi':           1., # ^
'GSTNGPi':           1., # ^
'GGPeGPi':           1., # LG14: no data available to decide; setting to diffuse improve selection properties
'GProtGPi':          1., #

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
'redundancyCMPfFSI':             3, # ^
'redundancyFSIFSI':              3, # ^
'redundancySTNFSI':              3, # ^
'redundancyGPeFSI':              3, # ^
'redundancyArkyFSI':             3, # ^

'redundancyPTNSTN':              3, # ^
'redundancyCMPfSTN':             3, # ^
'redundancyGPeSTN':              3, # ^
'redundancyProtSTN':              3, # ^

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

'redundancyCMPfGPi':             3, # ^
'redundancyMSNGPi':              3, # ^
'redundancySTNGPi':              3, # ^
'redundancyGPeGPi':              3, # ^
'redundancyProtGPi':             3, # ^

'cTypeCSNMSN':           'focused', # defining connection types for channel-based models (focused or diffuse) based on LG14 - refer to this paper for justification
'cTypePTNMSN':           'focused', # ^
'cTypeCMPfMSN':          'diffuse', # ^
'cTypeMSNMSN':           'diffuse', # ^
'cTypeFSIMSN':           'diffuse', # ^
'cTypeSTNMSN':           'diffuse', # ^
'cTypeGPeMSN':           'diffuse', # ^
'cTypeArkyMSN':          'diffuse', # ^

'cTypeCSNFSI':           'focused', # ^
'cTypePTNFSI':           'focused', # ^
'cTypeCMPfFSI':          'diffuse', # ^
'cTypeFSIFSI':           'diffuse', # ^
'cTypeSTNFSI':           'diffuse', # ^
'cTypeGPeFSI':           'diffuse', # ^
'cTypeArkyFSI':          'diffuse', # ^

'cTypePTNSTN':           'focused', # ^
'cTypeCMPfSTN':          'diffuse', # ^
'cTypeGPeSTN':           'focused', # ^
'cTypeProtSTN':          'focused', # ^

'cTypeCMPfGPe':          'diffuse', # ^
'cTypeMSNGPe':           'focused', # ^
'cTypeSTNGPe':           'diffuse', # ^
'cTypeGPeGPe':           'diffuse', # ^

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

'cTypeCMPfGPi':          'diffuse', # ^
'cTypeMSNGPi':           'focused', # ^
'cTypeSTNGPi':           'diffuse', # ^
'cTypeGPeGPi':           'diffuse', # LG14: no data available to decide; setting to diffuse improve selection properties
'cTypeProtGPi':          'diffuse', #

'parrotCMPf' :                True, # Should the CMPf be simulated using parrot neurons?
'stochastic_delays':          None, # If specified, gives the relative sd of a clipped Gaussian distribution for the delays
# For convenience, a few simulator variables are also set here
'whichTest':          'testFullBG', # task to be run (default: test the plausibility through deactivation simulations)
'nestSeed':                     20, # nest seed (affects input poisson spike trains)
'pythonSeed':                   10, # python seed (affects connection map)
'nbcpu':                         1, # number of CPUs to be used by nest
'durationH':                  '08', # max duration of a simulation, used by Sango cluster
'nbnodes':                     '1', # number of nodes, used by K computer
'tsimu':                     5000., # time duration of one simulation
}

