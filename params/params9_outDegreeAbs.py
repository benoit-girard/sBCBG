# test for RedundancyType=outDegreeAbs
# Running this file should give the exact same firing rates as params/params9_inDegreeAbs.py

the_scale = 4.

params = {'tSimu':  5000.,
          'LG14modelID':9 ,
          'IeMSN':  24.5,
          'IeFSI':  8.  ,
          'IeSTN':  9.5 ,
          'IeGPe':  12. ,
          'IeArky': 12. ,
          'IeProt': 12. ,
          'IeGPi':  11. ,
          
          'nbMSN':  2644.*the_scale, # original number of neurons, possibly scaled
          'nbFSI':  53.*the_scale  , # ^
          'nbSTN':  8.*the_scale   , # ^
          'nbGPe':  25.*the_scale  , # ^
          'nbArky': 5.*the_scale   , # ^
          'nbProt': 20.*the_scale  , # ^
          'nbGPi':  14.*the_scale  , # ^
          'nbCSN':  3000.*the_scale, # large pool of CSN neurons (split per channel)
          'nbPTN':  3000.*the_scale, # large pool of PTN neurons (possibly split per channel)
          'nbCMPf': 3000.*the_scale, # large pool of thalamic neurons (not split per channel)
          
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
'redundancyProtProt':            3, # ^
'redundancyArkyProt':            3, # ^

'redundancyMSNGPi':              3, # ^
'redundancySTNGPi':              3, # ^
'redundancyGPeGPi':              3, # ^
'redundancyProtGPi':             3, # ^
'redundancyCMPfGPi':             3,}
