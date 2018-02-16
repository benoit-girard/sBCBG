#!/apps/free/python/2.7.10/bin/python 

# defines the value of the parameters that will be used by testFullbG.py

params = {'LG14modelID':9,
          'GMSN':              4.37, # Same as in baseParams
          'GFSI':               1.3, # ^
          'GSTN':              1.38, # ^
          'GGPe':               1.3, # ^
          'GGPi':                1., # ^
          'IeGPe':              13., # ^
          'IeGPi':              11., # ^
          'nbMSN':            2644., # original number of neurons, possibly scaled
          'nbFSI':              53., # ^
          'nbSTN':               8., # ^
          'nbGPe':              25., # ^
          'nbGPi':              14., # ^
          'nbCSN':            3000., # large pool of CSN neurons (split per channel)
          'nbPTN':             100., # large pool of thalamic neurons (possibly split per channel)
          'nbCMPf':           3000., # large pool of thalamic neurons (not split per channel)
          'inDegCMPfMSN':  9., # inDegree from CMPf to * are 9
          'inDegCMPfFSI':  9., # ^
          'inDegCMPfSTN':  9.,  # ^
          'inDegCMPfGPe':  9.,  # ^
          'inDegCMPfGPi':  9.,  # ^
          'inDegPTNMSN':   3., # inDegree from PTN to striatum are 3 (not a sensitive parameter)
          'inDegPTNFSI':   3., # ^
          'parrotCMPf': True, # use parrot neurons for CMPf as well
          }

