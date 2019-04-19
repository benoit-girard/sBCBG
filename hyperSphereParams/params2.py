# Maximal hypersphere of parameters giving a plausibility score of 14/14

the_scale = 4.

params = {
  # original LG14 parameterization
  'LG14modelID':2 ,
  # hypersphere parameterization
  'IeMSN': 24.5 ,
  'IeFSI': 5.5 ,
  'IeSTN': 7.5 ,
  'IeGPe': 12. ,
  'IeGPi': 9.5 ,
  # original number of neurons, possibly scaled
  'nbMSN':  the_scale * 2644. ,
  'nbFSI':  the_scale * 53. ,
  'nbSTN':  the_scale * 8. ,
  'nbGPe':  the_scale * 25. ,
  'nbGPi':  the_scale * 14. ,
  'nbCSN':  the_scale * 3000. ,
  'nbPTN':  the_scale * 3000. ,
  'nbCMPf': the_scale * 3000. ,
  # duration of a simulation (in ms)
  'simDuration': 1000.
}
