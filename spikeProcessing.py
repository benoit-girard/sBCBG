#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy
import pylab
import nest.raster_plot as raster

filePath = '2017_2_13_14:42_00000/log/'
NUCLEI = ['STN','GPe']
fileList = {}
fileList['STN'] = ['STN-33652-0.gdf','STN-33652-1.gdf']
fileList['GPe'] = ['GPe-33653-0.gdf','GPe-33653-1.gdf']

# read files & combine data :
#----------------------------
data = {}
ts = {}
gids = {}
for N in NUCLEI:
  data[N] = None
  for f in fileList[N]:
    newdata = numpy.loadtxt(filePath+f)
    if data[N] is None:
      data[N] = newdata
    else:
      data[N] = numpy.concatenate((data[N],newdata))

  ts[N] = data[N][:,1] # complete time series of spiking events
  gids[N] = data[N][:,0] # corresponding list of firing neuron ID

# prepare signal : histogram of spiking events
#---------------------------------------------
signal = {}
t_bins = numpy.arange(numpy.amin(ts[NUCLEI[0]]),numpy.amax(ts[NUCLEI[0]]),1.0)
for N in NUCLEI:
  signal[N],bins = raster._histogram(ts[N], bins=t_bins)

# show histogram
#---------------
#for N in NUCLEI:
#  nbNeurons = len(numpy.unique(gids[N]))
#  heights = 1000 * signal[N] / (1.0 * nbNeurons)
#  pylab.bar(t_bins, heights, width=1.0, color="black")
#  pylab.show()

# compute FFT
#------------
T = 0.001 #1 ms, sampling period
for N in NUCLEI:
  Nb = signal[N].size

  #ps = numpy.abs(numpy.fft.fft(signal[N]))**2
  #freqs = numpy.fft.fftfreq(Nb,T)
  #idx = numpy.argsort(freqs)
  #idx = numpy.fft.fftshift(freqs)

  yf = numpy.fft.fft(signal[N])
  xf = numpy.linspace(0.,1./(2.*T),Nb/2)

  # show FFT
  #---------
  #pylab.plot(freqs[idx], ps[idx])
  #pylab.show()

  pylab.plot(xf[1:], 2.0/Nb * numpy.abs(yf[1:Nb//2])**2)
  #pylab.plot(xf[1:], numpy.abs(yf[1:Nb//2])**2)
  pylab.show()



