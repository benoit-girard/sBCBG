# Hand-tuned parameters with custom inDegree values for original parameterization #9
# CM/Pf are NOT modeled with parrot neurons
# Degrees of freedom explored: gains for all nuclei + tonic input for GPe and GPi


params   =       {'LG14modelID':          9,
                  'nbMSN':            2644.,
                  'nbFSI':              53.,
                  'nbSTN':               8.,
                  'nbGPe':              25.,
                  'nbGPi':              14.,
                  'nbCSN':            3000.,
                  'nbPTN':             100.,
                  'nbCMPf':              9.,
                  'GMSN':              4.37,
                  'GFSI':               1.3,
                  'GSTN':              1.38,
                  'GGPe':               1.3,
                  'GGPi':                1.,
                  'IeGPe':              13.,
                  'IeGPi':              11.,
                  'inDegCSNMSN':       100.,
                  'inDegPTNMSN':         1.,
                  'inDegCMPfMSN':        1.,
                  'inDegFSIMSN':        30., # 30 : according to Humphries et al. 2010, 30-150 FSIs->MSN
                  'inDegMSNMSN':        70., # 70 = 210/3 : according to Koos et al. 2004, cited by Humphries et al., 2010, on avg 3 synpase per MSN-MSN connection
                  'inDegSTNMSN':         0.,
                  'inDegGPeMSN':         0.,
                  'inDegCSNFSI':        50.,
                  'inDegPTNFSI':         1.,
                  'inDegSTNFSI':         2.,
                  'inDegGPeFSI':        25.,
                  'inDegCMPfFSI':        9.,
                  'inDegFSIFSI':        15., # 15 : according to Humphries et al., 2010, 13-63 FSIs->FSI
                  'inDegPTNSTN':        25.,
                  'inDegCMPfSTN':        9.,
                  'inDegGPeSTN':        25.,
                  'inDegCMPfGPe':        9.,
                  'inDegSTNGPe':         8.,
                  'inDegMSNGPe':      2644.,
                  'inDegGPeGPe':        25.,
                  'inDegMSNGPi':      2644.,
                  'inDegSTNGPi':         8.,
                  'inDegGPeGPi':        23.,
                  'inDegCMPfGPi':        9.,
                  'cTypeCMPfMSN': 'focused', # LG14: diffuse
                  'cTypeMSNMSN':  'focused', # LG14: diffuse
                  'cTypeGPeFSI':  'focused', # LG14: diffuse
                  'cTypeCMPfFSI': 'focused', # LG14: diffuse
                  'cTypeCMPfSTN': 'focused', # LG14: diffuse
                  'cTypeGPeSTN':  'diffuse', # LG14: focused
                  'cTypeCMPfGPe': 'focused', # LG14: diffuse
                  'cTypeCMPfGPi': 'focused', # LG14: diffuse
                  'parrotCMPf' :      False,
                  }

