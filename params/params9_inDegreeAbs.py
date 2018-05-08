# test for RedundancyType=inDegreeAbs
# Running this file should give the exact same firing rates as params/params9_outDegreeAbs.py

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
'RedundancyType':   'inDegreeAbs', # 
'redundancyCSNMSN':  114,
'redundancyPTNMSN':  1.66666666667,
'redundancyCMPfMSN': 5.38150332728,
'redundancyMSNMSN':  70.0,
'redundancyFSIMSN':  29.2471264368,
'redundancySTNMSN':  0, # connection not defined
'redundancyGPeMSN':  0, # connection not defined
'redundancyCSNFSI':  83.3333333333,
'redundancyPTNFSI':  1.66666666667,
'redundancySTNFSI':  0.746359649123,
'redundancyGPeFSI':  8.88250626566,
'redundancyCMPfFSI': 56.7406015038,
'redundancyFSIFSI':  38.6666666667,
'redundancyPTNSTN':  86.3333333333,
'redundancyCMPfSTN': 28.2943722944,
'redundancyGPeSTN':  20.645021645,
'redundancyCMPfGPe': 9.02257636122,
'redundancySTNGPe':  36.326002656,  # should be reduced to 32 when using scale=4
'redundancyMSNGPe':  6006.11952191,
'redundancyGPeGPe':  10.36,
'redundancyMSNGPi':  10616.1902098, # should be reduced to 10576.0 when using scale=4
'redundancySTNGPi':  30.1107692308,
'redundancyGPeGPi':  7.8634965035,
'redundancyCMPfGPi': 26.2610722611,
          }
