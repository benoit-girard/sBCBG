import os
import csv
import numpy as np
import re
import shutil

# the base parameters from which the sequence will be generated
source = {'THETA_STN':1., 'cTypeFSIFSI': 'diffuse', 'cTypeGPeGPe': 'diffuse', 'cTypeGPeSTN': 'focused', 'inDegGPeFSI': 0.33333, 'nbCSN': 12000.0, 'nbnodes': '1', 'cTypeCMPfSTN': 'diffuse', 'inDegMSNGPi': 0.33333, 'inDegMSNGPe': 0.33333, 'cTypeCMPfGPe': 'diffuse', 'nbGPe': 100.0, 'inDegGPeMSN': 0.33333, 'nbGPi': 56.0, 'cTypeCMPfGPi': 'diffuse', 'cTypePTNMSN': 'focused', 'inDegCMPfGPe': 0.33333, 'inDegSTNGPi': 0.33333, 'inDegSTNGPe': 0.33333, 'inDegCMPfGPi': 0.33333, 'inDegFSIMSN': 0.33333, 'nestSeed': 1, 'LG14modelID': 0, 'inDegSTNFSI': 0.33333, 'nbCMPf': 12000.0, 'inDegCSNMSN': 0.33333, 'nbPTN': 400.0, 'inDegMSNMSN': 0.33333, 'cTypeMSNMSN': 'diffuse', 'IeFSI': [5.00, 5.00, 5.50, 7.0, 7.5, 4.00, 5.75, 6.00, 11.25, 8.00, 7.0, 7.5, 9.75, 5.50, 6.50], 'parrotCMPf': True, 'cTypePTNFSI': 'focused', 'cTypeSTNFSI': 'diffuse', 'inDegFSIFSI': 0.33333, 'IeMSN': [23.75, 23.75, 24.50, 24.0, 23.5, 23.75, 24.25, 24.25, 24.50, 24.75, 24.0, 24.5, 24.50, 24.50, 24.50], 'cTypeGPeGPi': 'diffuse', 'nbMSN': 10576.0, 'nbcpu': 8, 'inDegCMPfMSN': 0.33333, 'inDegGPeGPe': 0.33333, 'cTypeGPeFSI': 'diffuse', 'GGPe': 1.0, 'GMSN': 1.0, 'inDegGPeGPi': 0.33333, 'nbCh': 1, 'cTypeCSNMSN': 'focused', 'GGPi': 1.0, 'cTypeMSNGPi': 'focused', 'IeSTN': [8.75, 7.50, 7.75, 9.0, 7.5, 9.00, 9.00, 9.00, 9.50, 9.25, 9.0, 9.5, 9.00, 9.25, 8.75], 'durationH': '08', 'cTypeMSNGPe': 'focused', 'cTypeSTNGPe': 'diffuse', 'cTypeCMPfMSN': 'diffuse', 'inDegPTNFSI': 0.33333, 'GFSI': 1.0, 'nbSTN': 32.0, 'nbFSI': 212.0, 'cTypeGPeMSN': 'diffuse', 'email': '', 'inDegPTNMSN': 0.33333, 'inDegCMPfSTN': 0.33333, 'inDegPTNSTN': 0.33333, 'cTypeCSNFSI': 'focused', 'cTypeCMPfFSI': 'diffuse', 'cTypeFSIMSN': 'diffuse', 'cTypePTNSTN': 'focused', 'IeGPe': [12.50, 9.00, 12.00, 13.0, 12.5, 10.00, 10.50, 12.00, 11.50, 12.00, 13.0, 12.0, 12.00, 12.00, 10.00], 'GSTN': 1.0, 'IeGPi': [10.00, 9.00, 9.50, 9.5, 9.0, 9.50, 10.50, 11.00, 11.00, 11.00, 10.5, 11.0, 11.00, 11.00, 12.00], 'inDegCMPfFSI': 0.33333, 'pythonSeed': 1, 'cTypeSTNGPi': 'diffuse', 'cTypeSTNMSN': 'diffuse', 'whichTest': 'testFullBG', 'inDegSTNMSN': 0.33333, 'inDegCSNFSI': 0.33333, 'inDegGPeSTN': 0.33333}



#-------------------------------------------------------------------------------
# Generates a sequence of parameter files
#-------------------------------------------------------------------------------
def generate(path):

  ##------------------------
  # Create destination folder
  ##------------------------ 
  if os.path.exists(path):
    shutil.rmtree(path)
  os.makedirs(path)



  ##----------------------
  # sequence parameters
  ##----------------------
  params = []
  params.append({'name' : 'G_offset',       'min' : 1.,     'max' : 1.6,   'delta' : 0.1}) # max value is excluded
  params.append({'name' : 'Theta_offset',    'min' : 0.,     'max' : 12.,   'delta' : 2.}) # max value is excluded
  params.append({'name' : 'Model',          'min' : 0,      'max' : 15,     'delta' : 1})   # max value is excluded
  params.append({'name' : 'Common_seed',    'min' : 1,      'max' : 2,     'delta' : 1})   # max value is excluded
  params.append({'name' : 'Nest_seed',      'min' : 1,     'max' : 2,     'delta' : 1})   # max value is excluded
  params.append({'name' : 'Py_seed',        'min' : 1,     'max' : 2,     'delta' : 1})   # max value is excluded 
  
  # parameter to introduce vatiation in transmission delays 
  stochastic_delays = [.1]

  print '===================\n   Parameters :\n==================='
  
  # Find number of values of each parameter
  dim = ()
  for param in params:
    nb_values = (int(round((param['max']-param['min'])/param['delta'],2)),)
    print round((param['max']-param['min']),1)/param['delta']
    dim += nb_values
    print param['name']+':\n', '['+str(param['min'])+'-'+str(param['max'])+'[', '/ step:', param['delta'], '/ '+re.sub("[^0-9]", "", str(nb_values))+' value(s)\n'
  
  print 'Generating', np.prod(dim), 'parameter files', dim

  
  ##----------------------
  # generate files
  ##---------------------- 

  # find values in Lienard solutions file for some parameters (e.g. THETA)
  with open('../solutions_simple_unique.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar="'")
    rows = []
    for row in spamreader:
      rows.append(row)

  NUCLEI = ['MSN', 'FSI', 'STN', 'GPe', 'GPi']
  id_config = 0
  base = 'params = '
  
  # save Ie parameter for all nuclei and all models
  Ie_sauv = {'Ie'+N : source['Ie'+N] for N in NUCLEI} 
  
  #iterate stochastic_delays
  for stoch in stochastic_delays:
    source['stochastic_delays'] = stoch
    id_model = params[2]['min']
     
    #iterate models
    while id_model < params[2]['max']:
      source_thr = float(rows[int(id_model+1)][50]) # with 50 the index of THETA_STN    
      source['LG14modelID'] = id_model
      for N in NUCLEI:
        source['Ie'+N] = Ie_sauv['Ie'+N][int(id_model)]

      G_offset = params[0]['min']   
      
      #iterate gain
      while G_offset < params[0]['max']:
        source["GGPe"] = G_offset
        source["GSTN"] = G_offset
        Theta_offset = params[1]['min']
        
        #iterate theta
        while Theta_offset < params[1]['max']:
          source["THETA_STN"] = round(source_thr + Theta_offset,2)

          Common_seed = params[3]['min']
          
          #iterate common seed
          while Common_seed < params[3]['max']:
            source['pythonSeed'] = Common_seed
            source['nestSeed'] = Common_seed

            #write file
            new_param = open(path+'/config_'+str(id_config)+'.py', 'w+')
            new_param.write(base + str(source))

            #increment config ID
            id_config += 1


            Common_seed += params[3]['delta']
          
          Theta_offset = round(Theta_offset + params[1]['delta'],2)
        
        G_offset = round(G_offset + params[0]['delta'],2)
      
      id_model += params[2]['delta']

  #print summary
  print 'Generated', id_config, 'parameter files'


def main():
  generate('sequence')

if __name__ == '__main__':
  main()

