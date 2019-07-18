#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time

for i in range(6,15):
    commande = "python run.py --platform Local --LG14modelID "+str(i)+" --whichTest testPlausibility --nbcpu -1 --nbCh 1 --gdf --custom hyperSphereParams/params"+str(i)+".py"
    print("### Running model "+str(i)+" ### \n"+commande)
    print("Starting from: "+time.asctime( time.localtime() ))
    os.system(commande)
