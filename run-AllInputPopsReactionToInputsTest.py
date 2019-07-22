#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time

modelID=0

for p in ['CSN','PTN','CMPf']:
    commande = "python run.py --platform Local --LG14modelID "+str(modelID)+" --whichTest testReactionToInput --inputPop "+p+" --nbcpu -1 --nbCh 1 --gdf --custom hyperSphereParams/params"+str(modelID)+".py"
    print("### Testing Inputs from "+p+" ### \n"+commande)
    print("Starting from: "+time.asctime( time.localtime() ))
    os.system(commande)
