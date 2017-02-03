#!/usr/bin/python
# -*- coding: utf-8 -*-

#import shlex
# import subprocess
import os
import time

execTime = time.localtime()
timeString = str(execTime[0])+'_'+str(execTime[1])+'_'+str(execTime[2])+'_'+str(execTime[3])+':'+str(execTime[4])

print 'Time:', timeString

# job for one parameterization test:
#----------------------------------------- 
def launchOneParameterizedRun(i):
  print "Ready to launch job",i

  IDstring = timeString+'_%05d' %(i)

  print 'Create subdirectory:',IDstring
  os.system('mkdir '+IDstring)
  os.system('cp LGneurons.py '+IDstring+'/')
  os.system('cp testFullBG.py '+IDstring+'/')
  os.system('cp solutions_simple_unique.csv '+IDstring+'/')
  os.system('cp __init__.py '+IDstring+'/')
  os.chdir(IDstring)
  os.system('mkdir log')

  # creation of the modelParams.py file that will correspond to the run at hand
  mltstr = '''#!/apps/free/python/2.7.10/bin/python

# defines the value of the parameters that will be used by testFullbG.py
# generated by sangoScript.py

interactive = True

params = {'LG14modelID': %2d,
          'nbMSN': 2644.,
          'nbFSI':   53.,
          'nbSTN':    8.,
          'nbGPe':   25.,
          'nbGPi':   14.,
          'nbCSN': %5.f,
          'nbPTN': %4.f,
          'nbCMPf':   9.,
          'GMSN':     %4.2f,
          'GFSI':     %4.2f,
          'GSTN':     %4.2f,
          'GGPe':     %4.2f,
          'GGPi':     %4.2f, 
          'IeGPe':    %3.1f,
          'IeGPi':    %3.1f,
          'inDegCSNMSN': 100.,
          'inDegPTNMSN':   1.,
          'inDegCMPfMSN':  1.,
          'inDegFSIMSN':  30., # according to Humphries et al. 2010, 30-150 FSIs->MSN
          'inDegMSNMSN':  70., # according to Koos et al. 2004, cited by Humphries et al., 2010, on avg 3 synpase per MSN-MSN connection
          'inDegCSNFSI':  50.,
          'inDegPTNFSI':   1.,
          'inDegSTNFSI':   2.,
          'inDegGPeFSI':  25.,
          'inDegCMPfFSI':  9.,
          'inDegFSIFSI':  15., # according to Humphries et al., 2010, 13-63 FSIs->FSI
          'inDegPTNSTN':  25.,
          'inDegCMPfSTN':  9.,
          'inDegGPeSTN':  25.,
          'inDegCMPfGPe':  9.,
          'inDegSTNGPe':   8.,
          'inDegMSNGPe':2644.,
          'inDegGPeGPe':  25.,
          'inDegMSNGPi':2644.,
          'inDegSTNGPi':   8.,
          'inDegGPeGPi':  23.,
          'inDegCMPfGPi':  9.,
          }
''' %(lg14modelid,nbcsn,nbptn,gmsn,gfsi,gstn,ggpe,ggpi,iegpe,iegpi)

  print 'Write modelParams.py'
  paramsFile = open('modelParams.py','w')
  paramsFile.writelines(mltstr)
  paramsFile.close()

  # execute the script file
  command = 'python testFullBG.py'
  os.system(command)

  os.chdir('..')

#===============================

'''
for iegpi in [10.,11.,12.]:
  launchOneParameterizedRun(i)
  i+=1
'''

i = 0                                                                                                                                                                           
# which LG14 parameterization to use?
lg14modelid = 2

# with which additional parameters?
nbcsn = 3000.
nbptn = 100.

gmsn=4.
gfsi=1.
gstn=1.4
ggpe=1.
ggpi=1.
iegpe=9.
iegpi=8.

'''
gmsn=4.
gfsi=1.
gstn=1.4
ggpe=1.
ggpi=1.
iegpe=11.
iegpi=10.
'''

#nbcsn=12000.
#nbptn=400.
launchOneParameterizedRun(i)
i+=1
