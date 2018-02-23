# Hand-tuned parameters with custom inDegree values for original parameterization #9
# Degrees of freedom explored: theta for MSN, FSI, STN and GPe + tonic input for GPe and GPi

params = {'LG14modelID':9,
          'THETA_MSN':  5.0, # overrides the parameter theta from mean-field models
          'THETA_FSI': 12.5, # ^
          'THETA_STN': 19.5, # ^
          'THETA_GPe':  10., # ^
          'IeGPe':   13.,
          'IeGPi':   11.,
          'nbMSN':            2644., # original number of neurons, possibly scaled
          'nbFSI':              53., # ^
          'nbSTN':               8., # ^
          'nbGPe':              25., # ^
          'nbGPi':              14., # ^
          'nbCSN':            3000., # large pool of CSN neurons (split per channel)
          'nbPTN':             100., # large pool of PTN neurons (possibly split per channel)
          'nbCMPf':           3000., # large pool of thalamic neurons (not split per channel)
          'inDegCMPfMSN':  9., # inDegree from CMPf to * are 9
          'inDegCMPfFSI':  9., # ^
          'inDegCMPfSTN':  9., # ^
          'inDegCMPfGPe':  9., # ^
          'inDegCMPfGPi':  9., # ^
          'inDegPTNMSN':   3., # inDegree from PTN to striatum are 3 (not a sensitive parameter)
          'inDegPTNFSI':   3., # ^
          'inDegCSNMSN':       100., # otherwise all other inDegrees are as in custom_noparrot_params9.py
          'inDegFSIMSN':        30.,
          'inDegMSNMSN':        70.,
          'inDegSTNMSN':         0.,
          'inDegGPeMSN':         0.,
          'inDegCSNFSI':        50.,
          'inDegSTNFSI':         2.,
          'inDegGPeFSI':        25.,
          'inDegFSIFSI':        15.,
          'inDegPTNSTN':        25.,
          'inDegGPeSTN':        25.,
          'inDegSTNGPe':         8.,
          'inDegMSNGPe':      2644.,
          'inDegGPeGPe':        25.,
          'inDegMSNGPi':      2644.,
          'inDegSTNGPi':         8.,
          'inDegGPeGPi':        23.,
          'parrotCMPf': True, # use parrot neurons for CMPf as well
          }
