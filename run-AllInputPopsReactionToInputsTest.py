#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time

for modelID in [10,11,12,13,14]:
#for modelID in [3]:
  for p in ['CSN','PTN','CMPf']:
    commande = "python run.py --platform Local --LG14modelID "+str(modelID)+" --whichTest testReactionToInput --inputPop "+p+" --tag "+p+" --nbcpu -1 --nbCh 1 --gdf --custom hyperSphereParams/params"+str(modelID)+".py"
    print("### Testing Inputs from "+p+" ### \n"+commande)
    print("Starting from: "+time.asctime( time.localtime() ))
    os.system(commande)