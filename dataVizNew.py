#!/apps/free/python/2.7.10/bin/python
#import sys
import nest.raster_plot
import argparse
import os

# This script loads the data in the gdf files produced by nest,
# and then uses the raster_plot module to plot the neural activity.
# first argument : name of the files without nest ID and cpu number
# second argument: nest ID
# third argument : nb of CPUs used, i.e., number of gdf files

dataFileName=[]
nbCPUs = 2

parser = argparse.ArgumentParser(description="Raster visualization.", formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=27))
parser._action_groups.pop()

RequiredNamed = parser.add_argument_group('mandatory arguments')
RequiredNamed.add_argument('--fileName', type=str, help='Name of the files without nestID and CPU number', default=None)
RequiredNamed.add_argument('--NestID', type=int, help='ID given by Nest.', default=None)

Optional = parser.add_argument_group('optional arguments')
Optional.add_argument('--path', type=str, help='Data path if not current directory', default='')
Optional.add_argument('--nbCPU', type=int, help='Number of CPUs', default=2)
Optional.add_argument('--nbChannels', type=int, help='Number of Channels (if more than one, give the smallestNestID for the NestID argument).', default=1)

cmd_args = parser.parse_args()

fstr    =cmd_args.fileName
NestID  =cmd_args.NestID

dataPath=cmd_args.path
nbCPUs  =cmd_args.nbCPU
nbCh    =cmd_args.nbChannels

for i in range(nbCPUs):
  for j in range(nbCh):
    fullName = dataPath+fstr+str(NestID+j)+'-'+str(i)+'.gdf'
    if os.stat(fullName).st_size >0:
      dataFileName.append(fullName)
    else:
      print(fullName,'excluded because empty')
print(dataFileName)

#-----------

nest.raster_plot.from_file(dataFileName,hist=True,title=fstr)

nest.raster_plot.show()
