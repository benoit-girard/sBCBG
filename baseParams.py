# Simulation defaults
# these defaults can be overrided via a custom python parameter file or via the commandline (in this order: commandline arguments take precedence over customParams.py, which take precendence over the defaults defined here)
params   =       {'durationH':         '04', # used by Sango
                  'durationMin':       '00', # used by Sango
                  'nbnodes':            '1', # used by K
                  'nestSeed':           20, # nest seed (affects input poisson spike trains)
                  'pythonSeed':         10, # python seed (affects connection map)
                  'nbcpu':                1,
                  'whichTest': 'testFullBG',
                  'nbCh':                 1,
                  'LG14modelID':          9,
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
                  'IeMSN':              0.,
                  'IeFSI':              0.,
                  'IeSTN':              0.,
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
                  'cTypeCSNMSN':  'focused', # defining connection types for channel-based models (focused or diffuse)
                  'cTypePTNMSN':  'focused',
                  'cTypeCMPfMSN': 'focused', # LG14: diffuse
                  'cTypeFSIMSN':  'diffuse',
                  'cTypeMSNMSN':  'focused', # LG14: diffuse
                  'cTypeSTNMSN':  'diffuse',
                  'cTypeGPeMSN':  'diffuse',
                  'cTypeCSNFSI':  'focused',
                  'cTypePTNFSI':  'focused',
                  'cTypeSTNFSI':  'diffuse',
                  'cTypeGPeFSI':  'focused', # LG14: diffuse
                  'cTypeCMPfFSI': 'focused', # LG14: diffuse
                  'cTypeFSIFSI':  'diffuse',
                  'cTypePTNSTN':  'focused',
                  'cTypeCMPfSTN': 'focused', # LG14: diffuse
                  'cTypeGPeSTN':  'diffuse', # LG14: focused
                  'cTypeCMPfGPe': 'focused', # LG14: diffuse
                  'cTypeSTNGPe':  'diffuse',
                  'cTypeMSNGPe':  'focused',
                  'cTypeGPeGPe':  'diffuse',
                  'cTypeMSNGPi':  'focused',
                  'cTypeSTNGPi':  'diffuse',
                  'cTypeGPeGPi':  'diffuse', # LG14: no data available to decide; setting to diffuse improve selection properties
                  'cTypeCMPfGPi': 'focused', # LG14: diffuse
                  'parrotCMPf' :      False,
                  }

