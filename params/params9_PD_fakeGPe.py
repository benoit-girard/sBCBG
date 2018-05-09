# Median parameters giving a plausibility score of 14/14 with all inDegree values set to 1/3

the_scale = 4. # should be 4

params = {'LG14modelID':9 ,
          'IeMSN':  24.5,
          'IeFSI':  8.  ,
          'IeSTN':  9.5 ,
          'IeGPe':  12. ,
          'IeGPi':  11. ,
          'GSTN':   1.5, # increased coupling gain
          'GGPe':   1.5, # ^
          'fakeGPeRecurrent': 65, # "fake" GPe recurrent collaterals coming from an artificial nucleus firing Poisson spike trains at 65 Hz
          'nbMSN':  2644.*the_scale, # original number of neurons, possibly scaled
          'nbFSI':  53.*the_scale  , # ^
          'nbSTN':  8.*the_scale   , # ^
          'nbGPe':  25.*the_scale  , # ^
          'nbGPi':  14.*the_scale  , # ^
          'nbCSN':  3000.*the_scale, # large pool of CSN neurons (split per channel)
          'nbPTN':  3000.*the_scale , # large pool of PTN neurons (possibly split per channel)
          'nbCMPf': 3000.*the_scale, # large pool of thalamic neurons (not split per channel)
}
