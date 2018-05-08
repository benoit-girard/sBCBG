# test for RedundancyType=outDegreeAbs
# Running this file should give the exact same firing rates as params/params9_inDegreeAbs.py

the_scale = 4.

params = {'LG14modelID':9 ,
          'IeMSN':  24.5,
          'IeFSI':  8.  ,
          'IeSTN':  9.5 ,
          'IeGPe':  12. ,
          'IeGPi':  11. ,
          'nbMSN':  2644.*the_scale, # original number of neurons, possibly scaled
          'nbFSI':  53.*the_scale  , # ^
          'nbSTN':  8.*the_scale   , # ^
          'nbGPe':  25.*the_scale  , # ^
          'nbGPi':  14.*the_scale  , # ^
          'nbCSN':  3000.*the_scale, # large pool of CSN neurons (split per channel)
          'nbPTN':  3000.*the_scale , # large pool of PTN neurons (possibly split per channel)
          'nbCMPf': 3000.*the_scale, # large pool of thalamic neurons (not split per channel)
'RedundancyType':   'outDegreeAbs', # by default all axons are hypothesized to target each dendritic tree at 3 different locations
'redundancyCSNMSN':              3, # ^
'redundancyPTNMSN':              3, # ^
'redundancyCMPfMSN':             3, # ^
'redundancyMSNMSN':              3, # ^
'redundancyFSIMSN':              3, # ^
'redundancySTNMSN':              3, # ^
'redundancyGPeMSN':              3, # ^
'redundancyCSNFSI':              3, # ^
'redundancyPTNFSI':              3, # ^
'redundancySTNFSI':              3, # ^
'redundancyGPeFSI':              3, # ^
'redundancyCMPfFSI':             3, # ^
'redundancyFSIFSI':              3, # ^
'redundancyPTNSTN':              3, # ^
'redundancyCMPfSTN':             3, # ^
'redundancyGPeSTN':              3, # ^
'redundancyCMPfGPe':             3, # ^
'redundancySTNGPe':              3, # ^
'redundancyMSNGPe':              3, # ^
'redundancyGPeGPe':              3, # ^
'redundancyMSNGPi':              3, # ^
'redundancySTNGPi':              3, # ^
'redundancyGPeGPi':              3, # ^
'redundancyCMPfGPi':             3, # ^
          }
