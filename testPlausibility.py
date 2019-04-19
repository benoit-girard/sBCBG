#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
## testPlausibility.py
##
## This script tests the 14 electrophysiological constraints from LG14
## Works both in single-channel and multi-channels cases

from iniBG import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
import nest.raster_plot as raster
from iniBG import *
from spikeProcessing import FanoFactor, OscIndex
from filter import lowpass
import os
import numpy as np
from modelParams import *
restFR = {} # this will be populated with firing rates of all nuclei, at rest
oscilPow = {} # Oscillations power and frequency at rest
oscilFreq = {}

#------------------------------------------
# Checks whether the BG model respects the electrophysiological constaints (firing rate at rest).
# If testing for a given antagonist injection experiment, specifiy the injection site in antagInjectionSite, and the type of antagonists used in antag.
# Returns [score obtained, maximal score]
# params possible keys:
# - nb{MSN,FSI,STN,GPi,GPe,CSN,PTN,CMPf} : number of simulated neurons for each population
# - Ie{GPe,GPi} : constant input current to GPe and GPi
# - G{MSN,FSI,STN,GPi,GPe} : gain to be applied on LG14 input synaptic weights for each population
#------------------------------------------
def checkAvgFR(showRasters=False,params={},antagInjectionSite='none',antag='',logFileName='', computeLFP=True, computePS=True):
  nest.ResetNetwork()
  initNeurons()

  showPotential = False # Switch to True to graph neurons' membrane potentials - does not handle well restarted simulations

  dataPath='log/'
  nest.SetKernelStatus({"overwrite_files":True}) # when we redo the simulation, we erase the previous traces

  nstrand.set_seed(params['nestSeed'], params['pythonSeed']) # sets the seed for the simulation

  simulationOffset = nest.GetKernelStatus('time')
  print('Simulation Offset: '+str(simulationOffset))
  if 'offsetDuration' not in params.keys():
    offsetDuration = 1000.
  else:
    offsetDuration = params['offsetDuration']

  if 'simDuration' not in params.keys():
    simDuration = 1000. # step duration period
    #simDuration = 50000. # ms VERY LONG
    #simDuration = 10000. # ms LONG
    #simDuration = 3000. # ms NORMAL
    #simDuration = 1000. * 100. # ms DISTANCE BASED
  else:
    simDuration = params['simDuration']

  # single or multi-channel?
  if params['nbCh'] == 1:
    connect_detector = lambda N: nest.Connect(Pop[N], spkDetect[N])
    disconnect_detector = lambda N: nest.Disconnect(Pop[N], spkDetect[N])
    connect_multimeter = lambda N: nest.Connect(multimeters[N], [Pop[N][0]])
  else:
    connect_detector= lambda N: [nest.Connect(Pop[N][i], spkDetect[N]) for i in range(len(Pop[N]))]
    disconnect_detector= lambda N: [nest.Disconnect(Pop[N][i], spkDetect[N]) for i in range(len(Pop[N]))]
    connect_multimeter = lambda N: nest.Connect(multimeters[N], [Pop[N][0][0]])

  #-------------------------
  # measures
  #-------------------------
  spkDetect={} # spike detectors used to record the experiment
  multimeters={} # multimeters used to record one neuron in each population
  expeRate={}

  if computeLFP:
    samplingRate = 2000.
    voltmeters = {} # voltmeters used to compute LFP
    maxRecord = 100 # 100 is greater than len(STN) and len(GPe) -> recording all neurons of both pops
    time_step = 1 / samplingRate
  focusNuclei = ['GPe','STN', 'GPi']

  antagStr = ''
  if antagInjectionSite != 'none':
    antagStr = antagInjectionSite+'_'+antag+'_'

  for N in NUCLEI:
    # 1000ms offset period for network stabilization
    spkDetect[N] = nest.Create("spike_detector", params={"withgid": True, "withtime": True, "label": antagStr+N, "to_file": storeGDF, 'start':offsetDuration+simulationOffset,'stop':offsetDuration+simDuration+simulationOffset})
    connect_detector(N)
    #spkDetect[N] = nest.Create("spike_detector", len([Pop[N][x] for x in range(len(Pop[N])) if x<maxRecord]), params={"withgid": True, "withtime": True, "label": antagStr+N, "to_file": storeGDF, 'start':offsetDuration+simulationOffset,'stop':offsetDuration+simDuration+simulationOffset})
    #nest.Connect([Pop[N][0]], spkDetect[N])
    if showPotential:
      # multimeter records only the last 200ms in one neuron in each population
      multimeters[N] = nest.Create('multimeter', params = {"withgid": True, 'withtime': True, 'interval': 0.1, 'record_from': ['V_m'], "label": antagStr+N, "to_file": False, 'start':offsetDuration+simulationOffset+simDuration-200.,'stop':offsetDuration+simDuration+simulationOffset})
      connect_multimeter(N)
    if computeLFP and N in focusNuclei:
      #voltmeters used to compute LFP
      voltmeters[N] = nest.Create("voltmeter", len([Pop[N][x] for x in range(len(Pop[N])) if x<maxRecord]), params = {'to_accumulator':True, 'start':offsetDuration+simulationOffset,'stop':offsetDuration+simDuration+simulationOffset})
      nest.SetStatus(voltmeters[N], {"withtime":True, 'interval': time_step*1000})
      nest.Connect(voltmeters[N], [Pop[N][x] for x in range(len(Pop[N])) if x<maxRecord])


  #-------------------------
  # Simulation
  #-------------------------
  nest.Simulate(simDuration+offsetDuration)

  score = 0

  text=[]
  frstr = "#" + str(params['LG14modelID'])+ " , " + antagInjectionSite + ', '
  s = '----- RESULTS -----'

  print(s)
  text.append(s+'\n')
  if antagInjectionSite == 'none':
    validationStr = "\n#" + str(params['LG14modelID']) + " , "
    frstr += "none , "

    if not os.path.exists("plots/"):
      os.makedirs("plots/")
    if not os.path.exists("data/"):
      os.makedirs("data/")
    if computePS:
      #Frequencies of interest (Hz)
      a = 15
      b = 30
      PSmetrics = {}

    for N in NUCLEI:
      strTestPassed = 'NO!'
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
      if expeRate[N] <= FRRNormal[N][1] and expeRate[N] >= FRRNormal[N][0]:
        # if the measured rate is within acceptable values
        strTestPassed = 'OK'
        score += 1
        validationStr += N + "=OK , "
      else:
      # out of the ranges
        if expeRate[N] > FRRNormal[N][1] :
          difference = expeRate[N] - FRRNormal[N][1]
          validationStr += N + "=+%.2f , " % difference
        else:
          difference = expeRate[N] - FRRNormal[N][0]
          validationStr += N + "=%.2f , " % difference

      frstr += '%f , ' %(expeRate[N])
      s = '* '+N+' - Rate: '+str(expeRate[N])+' Hz -> '+strTestPassed+' ('+str(FRRNormal[N][0])+' , '+str(FRRNormal[N][1])+')'
      print(s)
      text.append(s+'\n')
      restFR[N] = str(expeRate[N])

      oscilPow[N] = -1.
      oscilFreq[N] = -1.
      try:
        spikes_N = nest.GetStatus(spkDetect[N], keys="events")[0]['times'] # get the timing of all spikes
        data = np.bincount([int(i-offsetDuration-simulationOffset) for i in spikes_N], minlength=int(simDuration)) # discretize them in bins of 1ms
        ps = np.abs(np.fft.fft(data))**2
        time_step = 1 / 1000. # 1000 ms
        freqs = np.fft.fftfreq(data.size, time_step)
        idx = np.argsort(freqs)
        posi_spectrum = np.where((freqs[idx]>0) & (freqs[idx]<200)) # restrict the analysis to freqs < 200 Hz
        oscilPow[N] = np.max(ps[idx][posi_spectrum])
        oscilFreq[N] = freqs[idx][posi_spectrum][np.where(oscilPow[N] == ps[idx][posi_spectrum])[0][0]]
        #pl.plot(freqs[idx][posi_spectrum], ps[idx][posi_spectrum]) # simple plot
        #pl.show()
      except:
        print("Power spectrum computation failed - skipping")

      if N in focusNuclei:
        #--------------------
        #  Compute pseudo-LFP
        #--------------------
        if computeLFP:
          nbNeurons = len([Pop[N][x] for x in range(len(Pop[N])) if x<maxRecord])
          Vm = nest.GetStatus(voltmeters[N])[0]['events']['V_m']/nbNeurons

          # save LFP plots
          if True:
            plt.plot(range(len(Vm[500:1000])), Vm[500:1000],color='black', linewidth=.5) # plot from 500 to 1000ms
            plt.title(N+' LFP ('+str(nbNeurons)+' neurons)')
            plt.savefig('plots/'+N+' LFP.pdf')
            plt.close()

          filter = False
          if filter:
            freqCut = 250. # Cutoff frequency (Hz)
            freqTrans = 100. # Transition band size (Hz)
            hamming = 0
            fVm = lowpass(Vm, freqCut / samplingRate, freqTrans / samplingRate)
            fVm = fVm[(len(fVm)-len(Vm)) // 2 + hamming : -(len(fVm)-len(Vm)) // 2 - hamming]





        #--------------------
        #  Compute power spectrum
        #--------------------
        if computePS:
          PSmetrics[N] = []
          #--------------------
          #  Compute spikes histogram for PS
          #--------------------
          data = nest.GetStatus(spkDetect[N], keys="events")[0]['times']
          data = [round(x-offsetDuration-simulationOffset,4) for x in data]
          bins = np.arange(0, simDuration, 1000*time_step)
          data,bins = raster._histogram(data, bins=bins)
          PSmetrics[N] += [{'name' : 'spikes', 'data' : data}]

          if computeLFP:
            #--------------------
            #  Prepare LFP data for PS
            #--------------------
            if filter:
              PSmetrics[N] += [{'name' : 'filteredLFP', 'data' : fVm}]
            PSmetrics[N] += [{'name' : 'LFP', 'data' : Vm}]

          hanning = False # compute the Hanning window of the time series
          for metric in PSmetrics[N]:
            data = metric['data']
            if metric['name'] == 'spikes':
              #--------------------
              #  Compute Fano Factor from spikes
              #--------------------
              metric['FF'] = FanoFactor(data)


            #--------------------
            #  Compute PS
            #--------------------
            if hanning:
              data = data * np.hanning(len(data))
            ps = np.abs(np.fft.fft(data))**2
            freqs = np.fft.fftfreq(data.size, time_step)
            idx = np.argsort(freqs)
            posi_spectrum = np.where((freqs[idx]>1) & (freqs[idx]<150)) # restrict the analysis to specified freqs

            plot = True
            if plot:
              plt.plot(freqs[idx][posi_spectrum], ps[idx][posi_spectrum], color='black', linewidth=.5)
              plt.xlabel('Freq. [Hz]')
              plt.gca().get_yaxis().set_visible(False)
              plt.title(N+' Power Spectrum.pdf')
              plt.savefig("plots/"+N+'_PwSpec_'+metric['name']+'.pdf', bbox_inches='tight',)
              plt.close()

            saveData = True
            if saveData:
              with open("data/"+N+'_PwSpec_'+metric['name']+'.py', 'w') as file:
                file.write('ps = '+str({'power': np.ndarray.tolist(ps[idx][posi_spectrum]), 'freqs': np.ndarray.tolist(freqs[idx][posi_spectrum])}))
                file.close

            #--------------------
            #  Compute Oscillation Index and absolute power levels
            #--------------------
            metric['OI'] = OscIndex(ps[idx][posi_spectrum], freqs[idx][posi_spectrum], a, b)
            metric['power'] = np.prod(ps[idx][np.where((freqs[idx]>a) & (freqs[idx]<b))])




    if computePS:
      os.makedirs("report/")
      with open('report/OI.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=';',
                            quotechar="'", quoting=csv.QUOTE_MINIMAL)
        report = [[],[]]
        for N in focusNuclei:
          tmp = []
          for metric in PSmetrics[N]:
            report[0] += [N+'_OI_'+metric['name']]
            report[1] += [metric['OI']]
        writer.writerow(report[0])
        writer.writerow(report[1])
        csvfile.close()

      with open('report/FF.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=';',
                            quotechar="'", quoting=csv.QUOTE_MINIMAL)
        report = [[],[]]
        for N in focusNuclei:
          for metric in PSmetrics[N]:
            if 'FF' in metric:
              report[0] += [N+'_FF'+metric['name']]
              report[1] += [metric['FF']]
        writer.writerow(report[0])
        writer.writerow(report[1])
        csvfile.close()

      with open('report/power.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=';',
                            quotechar="'", quoting=csv.QUOTE_MINIMAL)
        report = [[],[]]
        for N in focusNuclei:
          tmp = []
          for metric in PSmetrics[N]:
            report[0] += [N+'_power_'+metric['name']]
            report[1] += [metric['power']]
        writer.writerow(report[0])
        writer.writerow(report[1])
        csvfile.close()


  else:
    validationStr = ""
    frstr += str(antag) + " , "
    for N in NUCLEI:
      expeRate[N] = nest.GetStatus(spkDetect[N], 'n_events')[0] / float(nbSim[N]*simDuration*params['nbCh']) * 1000
      if N == antagInjectionSite:
        strTestPassed = 'NO!'
        if expeRate[N] <= FRRAnt[N][antag][1] and expeRate[N] >= FRRAnt[N][antag][0]:
          # if the measured rate is within acceptable values
          strTestPassed = 'OK'
          score += 1
          validationStr += N + "_" + antag + "=OK , "
        else:
        # out of the ranges
          if expeRate[N] > FRRNormal[N][1] :
            difference = expeRate[N] - FRRNormal[N][1]
            validationStr += N + "_" + antag + "=+%.2f , " % difference
          else:
            difference = expeRate[N] - FRRNormal[N][0]
            validationStr += N + "_" + antag + "=%.2f , " % difference

        s = '* '+N+' with '+antag+' antagonist(s): '+str(expeRate[N])+' Hz -> '+strTestPassed+' ('+str(FRRAnt[N][antag][0])+' , '+str(FRRAnt[N][antag][1])+')'
        print(s)
        text.append(s+'\n')
      else:
        s = '* '+N+' - Rate: '+str(expeRate[N])+' Hz'
        print(s)
        text.append(s+'\n')
      frstr += '%f , ' %(expeRate[N])

  s = '-------------------'
  print(s)
  text.append(s+'\n')

  frstr+='\n'
  firingRatesFile=open(dataPath+'firingRates.csv','a')
  firingRatesFile.writelines(frstr)
  firingRatesFile.close()

  #print "************************************** file writing",text
  #res = open(dataPath+'OutSummary_'+logFileName+'.txt','a')
  res = open(dataPath+'OutSummary.txt','a')
  res.writelines(text)
  res.close()

  validationFile = open("validationArray.csv",'a')
  validationFile.write(validationStr)
  validationFile.close()



  #-------------------------
  # Displays
  #-------------------------
  if showRasters:
    for N in NUCLEI:
        dSD = nest.GetStatus(spkDetect[N],keys="events")[0]
        evs = dSD["senders"]
        ts = dSD["times"]
        times = [[ts[t] for t in range(len(ts)) if evs[t]==neuron] for neuron in Pop[N]]
        colors = [['black' for t in range(len(ts)) if evs[t]==neuron] for neuron in Pop[N]]

        try:
          plt.plot(ts, evs, "|", color='black', markersize=2)
          plt.title(N+ ' raster plot')
          plt.savefig("plots/"+N+'_rastePlot.pdf')
          plt.close()
        except AttributeError:
          print 'A neuron of the '+N+' didn\'t spike which caused an error in the raster plot --> skipping'


  if showRasters and interactive:
    displayStr = ' ('+antagStr[:-1]+')' if (antagInjectionSite != 'none') else ''
    for N in NUCLEI:
      # histograms crash in the multi-channels case
      nest.raster_plot.from_device(spkDetect[N], hist=(params['nbCh'] == 1), title=N+displayStr)

    if showPotential:
      pl.figure()
      nsub = 231
      for N in NUCLEI:
        pl.subplot(nsub)
        nest.voltage_trace.from_device(multimeters[N],title=N+displayStr+' #0')
        disconnect_detector(N)
        pl.axhline(y=BGparams[N]['V_th'], color='r', linestyle='-')
        nsub += 1
    pl.show()

  return score, 5 if antagInjectionSite == 'none' else 1


#-----------------------------------------------------------------------
def main():
  if len(sys.argv) >= 2:
    print("Command Line Parameters")
    paramKeys = ['LG14modelID',
                 'nbMSN',
                 'nbFSI',
                 'nbSTN',
                 'nbGPe',
                 'nbGPi',
                 'nbCSN',
                 'nbPTN',
                 'nbCMPf',
                 'GMSN',
                 'GFSI',
                 'GSTN',
                 'GGPe',
                 'GGPi',
                 'IeGPe',
                 'IeGPi',
                 'inDegCSNMSN',
                 'inDegPTNMSN',
                 'inDegCMPfMSN',
                 'inDegFSIMSN',
                 'inDegMSNMSN',
                 'inDegCSNFSI',
                 'inDegPTNFSI',
                 'inDegSTNFSI',
                 'inDegGPeFSI',
                 'inDegCMPfFSI',
                 'inDegFSIFSI',
                 'inDegPTNSTN',
                 'inDegCMPfSTN',
                 'inDegGPeSTN',
                 'inDegCMPfGPe',
                 'inDegSTNGPe',
                 'inDegMSNGPe',
                 'inDegGPeGPe',
                 'inDegMSNGPi',
                 'inDegSTNGPi',
                 'inDegGPeGPi',
                 'inDegCMPfGPi',
                 ]
    if len(sys.argv) == len(paramKeys)+1:
      print("Using command line parameters")
      print(str(sys.argv))
      i = 0
      for k in paramKeys:
        i+=1
        params[k] = float(sys.argv[i])
    else :
      print("Incorrect number of parameters")

  nest.set_verbosity("M_WARNING")

  instantiate_BG(params, antagInjectionSite='none', antag='')
  score = np.zeros((2))
  #mapTopology2D(show=True)
  score += checkAvgFR(params=params,antagInjectionSite='none',antag='',showRasters=True)

  # don't bother with deactivation tests if activities at rest are not within plausible bounds
  if score[0] < score[1]:
    print("Activities at rest do not match: skipping deactivation tests")
  else:
      #if params['nbCh'] == 1:
      #  # The following implements the deactivation tests without re-wiring the BG (faster but implemented only in single-channel case)
      #  for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
      #    ww = deactivate('GPe', a)
      #    score += checkAvgFR(params=params,antagInjectionSite='GPe',antag=a)
      #    reactivate('GPe', a, ww)

      #  for a in ['AMPA+NMDA+GABAA','AMPA','NMDA+AMPA','NMDA','GABAA']:
      #    ww = deactivate('GPi', a)
      #    score += checkAvgFR(params=params,antagInjectionSite='GPi',antag=a)
      #    reactivate('GPi', a, ww)
      #else:
      # The following implements the deactivation tests with re-creation of the entire BG every time (slower but also implemented for multi-channels)
      for a in ['AMPA','AMPA+GABAA','NMDA','GABAA']:
        instantiate_BG(params, antagInjectionSite='GPe', antag=a)
        score += checkAvgFR(params=params,antagInjectionSite='GPe',antag=a)

      for a in ['AMPA+NMDA+GABAA','AMPA','NMDA+AMPA','NMDA','GABAA']:
        instantiate_BG(params, antagInjectionSite='GPi', antag=a)
        score += checkAvgFR(params=params,antagInjectionSite='GPi',antag=a)

  #-------------------------
  print "******************"
  print "* Score:",score[0],'/',score[1]
  print "******************"

  #-------------------------
  # log the results in a file
  #-------------------------
  res = open('log/OutSummary.txt','a')
  for k,v in params.iteritems():
    res.writelines(k+' , '+str(v)+'\n')
  res.writelines("Score: "+str(score[0])+' , '+str(score[1]))
  res.close()

  res = open('score.txt','w')
  res.writelines(str(score[0])+'\n')
  res.close()

  # combined params+score output, makes it quicker to read the outcome of many experiments
  params['sim_score'] = score[0]
  params['max_score'] = score[1]
  with open('params_score.csv', 'wb') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in params.items():
       writer.writerow([key, value])
    for key, value in restFR.items():
       writer.writerow([key+'_Rate', value])
    for key, value in oscilPow.items():
       writer.writerow([key+'_Pow', value])
    for key, value in oscilFreq.items():
       writer.writerow([key+'_Freq', value])

#---------------------------
if __name__ == '__main__':
  main()
