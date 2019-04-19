#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import csv
import matplotlib.pyplot as plt
import nest.raster_plot as raster
from scipy import signal as sig
import os
import glob
import re
import pylab as pl


def process(filePath, nbCpu):

  print '\n\n ===================\n  Processing spikes \n ===================\n\n'

  NUCLEI = ['MSN','FSI','STN','GPe','GPi']
  fileList = {}
  rootFID = int(float(re.findall(r'\d+', glob.glob('log/*'+NUCLEI[0]+'*.gdf')[0])[0]))
  FID = 0
  for N in NUCLEI:
    fileList[N] = []
    cpu = 0
    while cpu<nbCpu:
      fileList[N].append(N+'-'+str(rootFID+FID)+'-'+str(cpu)+'.gdf')
      cpu += 1
    FID += 1

  showFFT = False

  # read files & combine data :
  #----------------------------
  data = {}
  ts = {}
  gids = {}

  for N in NUCLEI:
    cpt=0
    data[N] = None
    for f in fileList[N]:
      cpt+=1
      newdata = np.loadtxt('log/'+f, ndmin=2)
      if newdata.any():
        if data[N] is None:
          data[N] = newdata
        else:
          data[N] = np.concatenate((data[N],newdata))
    ts[N] = []
    gids[N] = []
    try:
      ts[N] = [x - 1000 for x in data[N][:,1]] # complete time series of spiking events (-1000 ms to remove initial delay)
      gids[N] = data[N][:,0] # corresponding list of firing neuron ID
    except:
      TypeError('No spikes')

  print 'Continuing'
  # prepare signal : histogram of spiking events
  #---------------------------------------------
  signal = {}
  binPeriod = 1.0 # in ms
  t_bins = np.arange(np.amin(ts[NUCLEI[0]]),np.amax(ts[NUCLEI[0]]),binPeriod)
  for N in NUCLEI:
    signal[N],bins = raster._histogram(ts[N], bins=t_bins)

  '''
  # show/save histogram
  #---------------
  plt.rcParams["figure.figsize"] = (16,6)
  ax = {}
  i = 0
  for N in NUCLEI:
    i += 1
    nbNeurons = len(np.unique(gids[N]))
    #heights = 1000 * signal[N] / nbNeurons # Normalize
    heights = signal[N]
    ax[N] = plt.subplot(3,2,i)

    ax[N].bar(t_bins, heights, width=2.0, color="black")
    plt.subplots_adjust(left=None, bottom=None, right=None, top=2.0, wspace=None, hspace=0.8)
    plt.ylabel('Spikes nb')
    plt.xlabel('Time [ms]')
    plt.title(N)
  if not os.path.exists("plots/"):
      os.makedirs("plots/")
  plt.savefig('plots/ActHisto.pdf', bbox_inches='tight')

  plt.clf()'''

  if not os.path.exists("plots/"):
    os.makedirs("plots/")

  # compute FFT
  #------------
  T = binPeriod/1000 #0.1 ms, sampling period
  OI = {}
  FF = {}
  PS = {}
  for N in NUCLEI:
    Nb = signal[N].size # number of sample points
    yf = np.fft.fft(signal[N])
    xf = np.linspace(0.,1./(2.*T),Nb//2)

    # show FFT
    #---------
    PS[N] = 2.0/Nb * np.abs(yf[1:Nb//2])**2
    if showFFT and N in ['STN', 'GPe']:
      pl.plot(xf[1:], PS[N], linewidth=0.5, color='red') # simple plot
      pl.show()
      '''plt.plot(xf[1:], PS[N])
      plt.xlabel('Freq. [Hz]')
      plt.show()
      plt.savefig("plots/"+N+'_PwSpec.pdf', bbox_inches='tight')
      plt.close()'''



  # Oscillation index computation :
  #--------------------------------
  #binPeriod = 5.
  #t_bins = np.arange(np.amin(ts[NUCLEI[0]]),np.amax(ts[NUCLEI[0]]),binPeriod)

  #Frequencies of interest
  a = 15
  b = 30

  fieldnames = []
  for N in NUCLEI:
    fieldnames.append(N+'_FF')
    fieldnames.append(N+'_OI'+str(a)+"-"+str(b))
    signal[N],bins = raster._histogram(ts[N], bins=t_bins)
    FF[N] = FanoFactor(signal[N])
    OI[N] = OscIndex(PS[N], a, b, binPeriod)
    print N, FF[N], OI[N]

  if not os.path.exists("report/"):
    os.makedirs("report/")
  if os.path.isfile('report/oscillations.csv'):
    os.remove('report/oscillations.csv')
  with open('report/oscillations.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile, delimiter=';',
                        quotechar="'", quoting=csv.QUOTE_MINIMAL)
    report = [[],[]]
    for N in NUCLEI:
      report[0] += [N+'_FF', N+'_OI('+str(a)+"-"+str(b)+')']
      report[1] += [FF[N], OI[N]]
    writer.writerow(report[0])
    writer.writerow(report[1])
    csvfile.close()


#------------------------------------------
# Computes Fano Factor from Kumar 2011
# FF[pop] = Var[pop] / E[pop]
# returns 0 if E[pop] == 0
#
# t_bins = np.arange(np.amin(ts[NUCLEI[0]]),np.amax(ts[NUCLEI[0]]),5.0)
#------------------------------------------
def FanoFactor(raster):
  mean = np.mean(raster)
  var = np.var(raster)
  if mean != 0:
    return var/mean
  else:
   return 0

#------------------------------------------
# Computes Oscillation index from Kumar 2011
# Integral from a to b (spectrum) / Integral from 0 to (sampling freq. / 2) (spectrum)
# Note : the given spectrum must already be truncated (i.e. with x < sampling freq. / 2)
# returns 0 if denominator == 0
#------------------------------------------
def OscIndex(PS, freqs, a=15, b=30):
  tot = PS.sum()

  if tot != 0:
    idx = np.argsort(freqs)
    posi_spectrum = np.where((freqs[idx]>a) & (freqs[idx]<b)) # restrict the analysis to freqs [a-b] Hz
    return PS[posi_spectrum].sum()/tot
  else:
   return 0

def main():
  process(os.path.split(os.getcwd())[-1])
  pass

#---------------------------
if __name__ == '__main__':
  main()
