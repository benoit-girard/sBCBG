'''
The purpose of this file is to process the results of Parkinson Disease from a sequence of simulations (Hugo's work)
'''

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.ndimage.filters as fil
import pickle
from textwrap import wrap
import sys


def extract(directory, LFP = True):
  os.chdir(directory)
  
  sim = []
  parameters = {}
  METRICS = ['FF', 'OI']
  #METRICS = ['FF', 'OI', 'power']
  cpt=0
  for filename in os.listdir(os.getcwd()): # foreach xp folder
    cpt+=1
    print cpt
    if(os.path.isdir(filename) and filename not in ['results', 'restrictedResults', 'plots']): # not the results folder
      try:
        exec(open(filename+'/modelParams.py').read()) 
        parameters.update(params)
        sim.append({'SimID': filename, 'Gain' : parameters["GSTN"], 'LG14modelID' : parameters["LG14modelID"],
                    'THETA_STN' : parameters['THETA_STN'], 'nestSeed' : parameters['nestSeed'], 'pythonSeed' : parameters['pythonSeed'], 'Delays_sd' : parameters['stochastic_delays']})
          
      except:
        ImportError('Parameters not found.')
      
      # read simulation results (oscillations)
      for metric in METRICS:
        with open(filename+'/report/'+metric+'.csv', 'rb') as csvfile:
          spamreader = csv.reader(csvfile, delimiter=';', quotechar="'")
          results = []
          for row in spamreader:
            results.append(row)
          for i in range(len(results[0])):
            sim[-1][results[0][i]] = results[1][i]

      with open(filename+'/params_score.csv', mode='r') as ratesFile:
        fastMethod = True
        reader = csv.reader(ratesFile)
        for rows in reader:
          k = rows[0]
          if '_Rate' in k:
            sim[-1][k] = rows[1]
  return sim


def analysis(results, directory, outPath = 'results', plotstyle=True):
  if plotstyle:
    cmap='coolwarm'
    interpolation='bilinear'
  else:
    cmap='Greys'
    interpolation='None'
  NUCLEI = ['GPe', 'STN']

  ##----------------------
  # sequence parameters
  ##----------------------
  params = []
  params.append({'name' : 'Nucleus',        'min' : 0.,     'max' : 2.,     'delta' : 1.})

  params.append({'name' : 'G_alpha',       'min' : 1.,     'max' : 1.6,   'delta' : 0.1})
  params.append({'name' : 'Theta_offset',    'min' : 0.,     'max' : 12.,   'delta' : 2.})
  params.append({'name' : 'Model',          'min' : 0,     'max' : 15,    'delta' : 1})
  params.append({'name' : 'Common_seed',      'min' : 1,     'max' : 11,     'delta' : 1})
  params.append({'name' : 'Nest_seed',      'min' : 1,     'max' : 2,     'delta' : 1})
  params.append({'name' : 'Py_seed',        'min' : 1,     'max' : 2,     'delta' : 1})
  params.append({'name' : 'Delays_sd',        'min' : 0.,     'max' : 0.15,     'delta' : 0.05})



  # Find dimension of result matrix
  dim = []
  for param in params:
    dim += (int(round((param['max']-param['min'])/param['delta'],2)),)
  print dim

  METRICS = []
  METRICS.append({'name' : 'FFspikes', 'mat' : np.zeros(tuple(dim)), 'NUCLEI':NUCLEI[:]})
  #METRICS.append({'name' : 'OI_filteredLFP', 'mat' : np.zeros(tuple(dim)), 'NUCLEI':NUCLEI[:]})
  METRICS.append({'name' : 'OI_LFP', 'mat' : np.zeros(tuple(dim)), 'NUCLEI':NUCLEI[:]})
  METRICS.append({'name' : 'OI_spikes', 'mat' : np.zeros(tuple(dim)), 'NUCLEI':NUCLEI[:]})
  #METRICS.append({'name' : 'power_spikes', 'mat' : np.zeros(tuple(dim)), 'NUCLEI':NUCLEI[:]})
  #METRICS.append({'name' : 'power_LFP', 'mat' : np.zeros(tuple(dim)), 'NUCLEI':NUCLEI[:]})
  addNUCLEI = ['GPi']
  dim[0]+=len(addNUCLEI)
  METRICS.append({'name' : 'Rate', 'mat' : np.zeros(tuple(dim)), 'NUCLEI':NUCLEI[:]+addNUCLEI})
  
  max = [{M['name']: {key : [0. for x in NUCLEI] for key in ['max', 'THETA_STN', 'Gain', 'xp']} for M in METRICS} for i in range(15)]

  with open('../solutions_simple_unique.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar="'")
    solutions = []
    for row in spamreader:
      solutions.append(row)

  computePS=[]

  for dict in results:
    open(dict['SimID']+'/data/__init__.py','w').write('#!/usr/bin/env python\n# -*- coding: utf-8 -*-\n')
    open(dict['SimID']+'/__init__.py','w').write('#!/usr/bin/env python\n# -*- coding: utf-8 -*-\n')
    source_thr = float(solutions[dict['LG14modelID']+1][50]) # with 50 the index of THETA_STN
    

    if dict['nestSeed'] == 0 and dict['LG14modelID'] != 99:
      if round(dict['Gain'], 3) == round(params[1]['max']-params[1]['delta'], 3) and round(dict['THETA_STN']-source_thr, 3) == round(params[2]['min'], 3):
        computePS.append({'condition':'PD','SimID':dict['SimID'], 'commonSeed':dict['nestSeed'], 'Model':dict['LG14modelID']})
      elif round(dict['Gain'], 3) == round(params[1]['min'], 3) and round(dict['THETA_STN']-source_thr, 3) == round(params[2]['min'], 3):
        computePS.append({'condition':'normal','SimID':dict['SimID'], 'commonSeed':dict['nestSeed'], 'Model':dict['LG14modelID']})
    
    for M in METRICS:
      for N in range(len(M["NUCLEI"])):     
      
        #if dict['pythonSeed'] < params[4]['max']:
          if dict['Gain'] == 1. and int(round(((round(dict['THETA_STN']-source_thr,2))/params[2]['delta']),1)) == 0 and dict['LG14modelID'] == 9 and dict['pythonSeed']==1:
            print dict['SimID'], dict['Delays_sd']
          x = float(dict[M["NUCLEI"][N]+'_'+M['name']])
          '''if x > max[int(dict['LG14modelID'])][M['name']]['max'][N]:
            max[int(dict['LG14modelID'])][M['name']]['max'][N] = x
            max[int(dict['LG14modelID'])][M['name']]['THETA_STN'][N] = dict['THETA_STN']
            max[int(dict['LG14modelID'])][M['name']]['Gain'][N] = dict['Gain']'''
          
          M['mat'][N,int(round((dict['Gain']-params[1]['min'])/params[1]['delta'],2)),int(round(((round(dict['THETA_STN']-source_thr,2))/params[2]['delta']),1)),int(dict['LG14modelID'])-params[3]['min'],int(round((dict['nestSeed']-params[4]['min'])/params[4]['delta'],2)),0,0, int(round((dict['Delays_sd']-params[7]['min'])/params[7]['delta'],2))] = x

  recomputePS = False
  if recomputePS:
    # recompute PS with the right y axis range  
    if not os.path.exists('results/PS_oneseed'):
      os.makedirs('results/PS_oneseed')
    
    for M in METRICS:
      if M['name'] in ['OI_LFP','OI_spikes']:
        if not os.path.exists('results/PS_oneseed/'+M['name']):
          os.makedirs('results/PS_oneseed/'+M['name'])
        for N in NUCLEI:
          if not os.path.exists('results/PS_oneseed/'+M['name']+'/'+N):
            os.makedirs('results/PS_oneseed/'+M['name']+'/'+N)

    
    modules = {N:{} for N in NUCLEI}
    for M in METRICS:
      if M['name'] in ['OI_LFP','OI_spikes']:
        max = {N:-1 for N in NUCLEI}
        for sim in computePS:
          sys.path.insert(0, sim['SimID']+'/data')
          for N in NUCLEI:
            try:
              reload(modules[N][M['name'].replace('OI_','')])
            except KeyError:
              modules[N][M['name'].replace('OI_','')] = __import__(N+'_PwSpec_'+M['name'].replace('OI_',''))
            tmp = np.max(modules[N][M['name'].replace('OI_','')].ps['power'])
            if tmp > max[N]:
              max[N] = tmp
          sys.path.pop(0)
        
        for N in NUCLEI:
          M['maxPS_'+N] = max[N]
   
    for sim in computePS:
      sys.path.insert(0, sim['SimID']+'/data')
      for M in METRICS:
        if M['name'] in ['OI_LFP','OI_spikes']:
          for N in NUCLEI:
            reload(modules[N][M['name'].replace('OI_','')])
            plt.plot(modules[N][M['name'].replace('OI_','')].ps['freqs'], modules[N][M['name'].replace('OI_','')].ps['power'], linewidth=.5)
            plt.xlabel('Freq. [Hz]')
            plt.ylim(0,M['maxPS_'+N])
            #plt.gca().get_yaxis().set_visible(False)
            plt.title('Model'+str(sim['Model'])+N+'_'+sim['condition']+'_'+M['name']+'_'+'_seed'+str(sim['commonSeed']))       

            plt.savefig('results/PS_oneseed/'+M['name']+'/'+N+'/model'+str(sim['Model'])+'_seed'+str(sim['commonSeed'])+'_'+sim['condition']+'.png', bbox_inches='tight',)          
            plt.close() 
      sys.path.pop(0)

    print 'all PS computed'

  from pylab import rcParams
  rcParams['figure.figsize'] = 8, 9



  if not os.path.exists(outPath):
    os.makedirs(outPath)

  matrices = []
  
  for M in METRICS:
    
    if not os.path.exists("results/"+M['name']):
      os.makedirs(outPath+'/'+M['name'])
    
    for sd in range(int(round((params[7]['max']-params[7]['min'])/params[7]['delta'],2))):
      for N in range(len(M['NUCLEI'])):  
        if not os.path.exists(outPath+'/'+M['name']+'/'+M['NUCLEI'][N]):
          os.makedirs(outPath+'/'+M['name']+'/'+M['NUCLEI'][N])  
        for model in range(int(params[3]['min']), int(params[3]['max'])):
          #each model
          mat = np.mean(M['mat'][N,:,:,model-int(params[3]['min']), :, 0,0,sd], axis=2)
          normalizedMat = normalize(mat)
          matrices.append({'mat' : normalizedMat, 'model' : model, 'nucleus' : N, 'metric' : M['name']})
          plt.matshow(normalizedMat, interpolation=interpolation, cmap=cmap, origin='lower')

          for coordY in np.linspace(-1.29,6.21,num=10,endpoint=False)[1:]:
            plt.text(7.25,coordY,'|',fontsize=11)

          plt.text(6.3,6.55,'Norm.',fontsize=11)
          plt.text(7.1,6.55,'Raw',fontsize=11)
          plt.text(7.1,6.21,str(round(np.max(mat),2)),fontsize=11)
          plt.text(7.1,-1.29,str(round(np.min(mat),2)),fontsize=11)
          plt.xticks(range(int(round((params[2]['max']-params[2]['min'])/params[2]['delta'],2))),np.linspace(params[2]['min'], params[2]['max']-params[2]['delta'], int(round((params[2]['max']-params[2]['min'])/params[2]['delta'],2))))
          plt.yticks(range(int(round((params[1]['max']-params[1]['min'])/params[1]['delta'],2))),np.linspace(params[1]['min'], params[1]['max']-params[1]['delta'], int(round((params[1]['max']-params[1]['min'])/params[1]['delta'],2))))
          plt.xlabel('STN THETA offset')
          plt.ylabel('Gain offset')
          plt.title("\n".join(wrap(M['name']+' computed over 0-150Hz freqs in model '+str(model)+'\'s '+M['NUCLEI'][N]+' (delays sd='+str(sd)+') averaged across ('+str(int(round((params[4]['max']-params[4]['min'])/params[4]['delta'],2)))+' seeds)', 60)), y=1.15)
          plt.gca().xaxis.tick_bottom()
          plt.colorbar()
          
          for (i, j), z in np.ndenumerate(mat):
            plt.text(j, i, round(z,2), ha='center', va='center')

          plt.savefig(outPath+"/"+M['name']+'/'+M['NUCLEI'][N]+'/model_'+str(model)+'sd_'+str(sd)+'.png')
          plt.close()

        #model averaging
        mat = np.mean(np.mean(M['mat'][N,:,:,:,:, 0,0,sd], axis=3), axis=2)
        normalizedMat = normalize(mat)
        plt.matshow(normalizedMat, interpolation=interpolation, cmap=cmap, origin='lower')
        plt.text(6.3,6.55,'Norm.',fontsize=11)
        
        for coordY in np.linspace(-1.29,6.21,num=10,endpoint=False)[1:]:
          plt.text(7.25,coordY,'|',fontsize=11)
        
        plt.text(7.1,6.55,'Raw',fontsize=11)
        plt.text(7.1,6.21,str(round(np.max(mat),2)),fontsize=11)
        plt.text(7.1,-1.29,str(round(np.min(mat),2)),fontsize=11)
        plt.xticks(range(int(round((params[2]['max']-params[2]['min'])/params[2]['delta'],2))),np.linspace(params[2]['min'], params[2]['max']-params[2]['delta'], int(round((params[2]['max']-params[2]['min'])/params[2]['delta'],2))))
        plt.yticks(range(int(round((params[1]['max']-params[1]['min'])/params[1]['delta'],2))),np.linspace(params[1]['min'], params[1]['max']-params[1]['delta'], int(round((params[1]['max']-params[1]['min'])/params[1]['delta'],2))))
        plt.xlabel('STN THETA offset')
        plt.ylabel('Gain offset')
        plt.title("\n".join(wrap(M['name']+' computed over 0-150Hz freqs in '+M['NUCLEI'][N]+' (delays sd='+str(sd)+') averaged across '+str(int(round((params[3]['max']-params[3]['min'])/params[3]['delta'],2)))+'models ('+str(int(round((params[4]['max']-params[4]['min'])/params[4]['delta'],2)))+' seeds)', 60)), y=1.15)
        plt.gca().xaxis.tick_bottom()
        plt.colorbar()

        for (i, j), z in np.ndenumerate(mat):
          plt.text(j, i, round(z,2), ha='center', va='center')
      
        plt.savefig(outPath+"/"+M['name']+'/'+M['NUCLEI'][N]+'sd_'+str(sd)+'.png')
        plt.close() 
'''
  file = open(directory+'.txt','w')
  pickle.dump(METRICS, file)
  file.close()

  if not os.path.exists('results/ranked'):
    os.makedirs('results/ranked')
  matrices = sortByDistance(matrices)
  modelsDistances = [{'name':'Model '+str(model), 'tot':0, 'cpt':0} for model in range(15)]
  metricsDistances = [{'name':M['name'], 'tot':0, 'cpt':0} for M in METRICS]
  modSTN = {metric['name'] : [0 for i in range(15)] for metric in METRICS}
  modGPe = {metric['name'] : [0 for i in range(15)] for metric in METRICS}
  IDmodel = {metric['name'] : [i for i in range(15)] for metric in METRICS}
  for i in range(len(matrices)):
    mat = matrices[i]
    plt.matshow(mat['mat'], interpolation='nearest', cmap='Greys', origin='lower', extent=[params[2]['min'], params[2]['max'], params[1]['min']*10, params[1]['max']*10])
    plt.xlabel('STN THETA offset')
    plt.ylabel('Gain offset')
    plt.gca().xaxis.tick_bottom()
    plt.colorbar()
  
    plt.savefig('results/ranked/'+str(i+1)+'.png')
    plt.close() 

    filter(lambda x: x['name'] == NUCLEI[mat['nucleus']], nucleiDistances)[0]['tot'] += mat['distance']
    filter(lambda x: x['name'] == NUCLEI[mat['nucleus']], nucleiDistances)[0]['cpt'] += 1.
   
    if  NUCLEI[mat['nucleus']] == 'GPe':
      modGPe[mat['metric']][mat['model']] = mat['distance']
    elif  NUCLEI[mat['nucleus']] == 'STN':
      modSTN[mat['metric']][mat['model']] = mat['distance']
    
    
    filter(lambda x: x['name'] == 'Model '+str(mat['model']), modelsDistances)[0]['tot'] += mat['distance']
    filter(lambda x: x['name'] == 'Model '+str(mat['model']), modelsDistances)[0]['cpt'] += 1.
    
    filter(lambda x: x['name'] == mat['metric'], metricsDistances)[0]['tot'] += mat['distance']
    filter(lambda x: x['name'] == mat['metric'], metricsDistances)[0]['cpt'] += 1.

  colors = ['green', 'blue', 'red']
  for M in METRICS:
    plt.scatter(modSTN[M['name']], modGPe[M['name']], color = colors[0], label = M['name'], alpha=0.5, edgecolors='none')
    del colors[0]
    for i, txt in enumerate(IDmodel[M['name']]):
      plt.annotate(txt, (modSTN[M['name']][i],modGPe[M['name']][i]))
  plt.legend()
  plt.xlabel('Distance in STN')
  plt.ylabel('Distance in GPe')
  #plt.show()
  plt.close()



  x_labels = [model['name'] for model in modelsDistances]
  x = range(len(x_labels))
  plt.bar(x, [model['tot']/model['cpt'] for model in modelsDistances])
  plt.xticks(x, x_labels)
  plt.ylabel('Distance from reference matrix')
  #plt.show()
  plt.close()



  x_labels = [N['name'] for N in nucleiDistances]
  x = range(len(x_labels))
  plt.bar(x, [N['tot']/N['cpt'] for N in nucleiDistances])
  plt.xticks(x, x_labels)
  plt.ylabel('Distance from reference matrix')
  plt.show()
  plt.close()
  root = os.getcwd()
  parameters = {}



  x_labels = [M['name'] for M in metricsDistances]
  x = range(len(x_labels))
  plt.bar(x, [M['tot']/M['cpt'] for M in metricsDistances])
  plt.xticks(x, x_labels)
  plt.ylabel('Distance from reference matrix')
  #plt.show()
  plt.close()


'''

'''

  root = os.getcwd()
  parameters = {}
  if not os.path.exists("results/maxOI_RasterPlot"):
    os.makedirs("results/maxOI_RasterPlot")
  for N in range(len(NUCLEI)):
    if not os.path.exists("results/maxOI_RasterPlot/"+NUCLEI[N]):
      os.makedirs("results/maxOI_RasterPlot/"+NUCLEI[N])

  for file in os.listdir(os.getcwd()):
    if os.path.isdir(file) and file != 'results':
      os.chdir(file)
      
      try:
        exec(open('modelParams.py').read()) 
        parameters.update(params)
      except:
        raise ImportError('Couldn\'t find the file')
      for model in range(15):
        if int(parameters["LG14modelID"]) == model:

          for N in range(len(NUCLEI)):
            if round(parameters["GSTN"],2) == round(max[model]['Gain'][N],2) and round(parameters['THETA_STN'],2) == round(max[model]['THETA_STN'][N],2):
              max[model]['xp'][N] = file
              rasterPlot(getSpikes(os.getcwd()+'/log/',NUCLEI[N]), NUCLEI[N], os.getcwd()+'/log/', DirOut='../results/maxOI_RasterPlot/'+NUCLEI[N]+'/model_'+str(model)+'_Gain_'+str(round(parameters["GSTN"],2))+'_THETA_'+str(round(max[model]['THETA_STN'][N],2))+'.png')
      os.chdir(root)

'''

def normalize(data):
  return (data-np.min(data))/(np.max(data)-np.min(data))

##########################################################################################
##    Sort the matrices according to their distance with the BG & JL 2017 pattern       ##
##    Plot plt.matshow(reference, cmap='Greys', origin='lower') for an illustration     ##
##########################################################################################
def sortByDistance(matrices):
  reference = np.rot90([[1 if x+y > 5 else 0. for y in range(6)] for x in range(6)], k=3)
  reference = fil.gaussian_filter(reference, sigma=1.2) #apply gaussian filter
  reference = normalize(reference)

  [dic.update({'distance':np.sum(np.square(reference-dic['mat']))}) for dic in matrices]
  return sorted(matrices, key=lambda x: x['distance'])

def main():
  path = 'newvar'
  outPath = 'results'
  analysis(extract(path), path, outPath)

if __name__ == '__main__':
  main()

