# test for RedundancyType=outDegreeCons
# Running this file should give slightly different results from params/params9_outDegreeAbs.py or params/params9_inDegreeAbs.py
# Reason: taking 1/3 ratio among the possible min/max range of axo-dendritic contacts is not exactly the same as taking an outDegree of 3 (=each target dendritic tree is contacted in three places by the same source axon)
# Yet these parameters should fulfill all 14 plausibility objectives

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
'RedundancyType':   'outDegreeCons', # 
'redundancyCSNMSN':     0.333333333,
'redundancyPTNMSN':     0.333333333,
'redundancyCMPfMSN':    0.333333333,
'redundancyMSNMSN':     0.333333333,
'redundancyFSIMSN':     0.333333333,
'redundancySTNMSN':     0.333333333,
'redundancyGPeMSN':     0.333333333,
'redundancyCSNFSI':     0.333333333,
'redundancyPTNFSI':     0.333333333,
'redundancySTNFSI':     0.333333333,
'redundancyGPeFSI':     0.333333333,
'redundancyCMPfFSI':    0.333333333,
'redundancyFSIFSI':     0.333333333,
'redundancyPTNSTN':     0.333333333,
'redundancyCMPfSTN':    0.333333333,
'redundancyGPeSTN':     0.333333333,
'redundancyCMPfGPe':    0.333333333,
'redundancySTNGPe':     0.333333333,
'redundancyMSNGPe':     0.333333333,
'redundancyGPeGPe':     0.333333333,
'redundancyMSNGPi':     0.333333333,
'redundancySTNGPi':     0.333333333,
'redundancyGPeGPi':     0.333333333,
'redundancyCMPfGPi':    0.333333333,
          }
