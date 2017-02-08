#!/apps/free/python/2.7.10/bin/python
import sys
import nest.raster_plot

# This script loads the data in the gdf files produced by nest, 
# and then uses the raster_plot module to plot the neural activity.
# first argument : ID of the files, first four digits chosen by nest when saving the gdf files
#                  example: file MSN-89621-00.gdf => ID = 8962
# second argument: nb of CPUs used, i.e., number of gdf files


NUCLEI={'MSN':1,'FSI':2,'STN':3,'GPe':4,'GPi':5}

deactivationList = []
for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
  deactivationList.append('GPe_'+a+'_')
for a in ['All','AMPA','NMDA+AMPA','NMDA','GABAA']:
  deactivationList.append('GPi_'+a+'_')

dataFileNames={'MSN':{},'FSI':{},'STN':{},'GPe':{},'GPi':{}}

nbCPUs = 10
if len(sys.argv) >= 2:
  fstr = '-'+sys.argv[1]
  if len(sys.argv) >= 3:
    nbCPUs = int(sys.argv[2])
else:
  print "please provide file ID"
  exit()

for N,I in NUCLEI.iteritems():
  dataFileNames[N]['none']=[]
  for i in range(nbCPUs):
    dataFileNames[N]['none'].append(N+fstr+str(I)+'-0'+str(i)+'.gdf')
  print dataFileNames[N]['none']

for d in deactivationList:
  for N,I in NUCLEI.iteritems():
    dataFileNames[N][d]=[]
    for i in range(nbCPUs):
      dataFileNames[N][d].append(d+N+fstr+str(I)+'-0'+str(i)+'.gdf')

#-----------

for N,I in NUCLEI.iteritems():
  nest.raster_plot.from_file(dataFileNames[N]['none'],hist=True,title=N)

nest.raster_plot.show()
