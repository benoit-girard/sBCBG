'''
The purpose of this file is to plot the power spectrums of various conditions
'''


import os
import shutil
from shutil import copyfile
import numpy as np
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
	focus = ['15-30', '0-250']
	shutil.rmtree("results", ignore_errors=True)
	linestyles = [[['blue', '--'], ['blue', '-']], [['red', ':'], ['red', '--']]]
	modelsParams = {n : [] for n in range(15)}
	for file in os.listdir('.'):
		if os.path.isdir(file):
			try:
				exec(open(file+'/modelParams.py').read()) 
				params['xpID'] = file
				modelsParams[params['LG14modelID']].append(params.copy())	 
			except:
				ImportError('Parameters not found.')
	NUCLEI = ['GPe', 'STN']
	METRICS = ['filteredLFP', 'LFP', 'spikes']
	os.makedirs("results")
	for foc in focus:
		os.makedirs("results/"+foc)
	for model in modelsParams:
		os.makedirs("results/"+foc+"/model_"+str(model))

		tmpGmin = []
		tmpGmax = []
		tmp = [[0,0],[0,0]]
		maxTH = -1
		maxTHind = []
		maxG = -1
		maxGind = []
		for i in range(len(modelsParams[model])):
			tmpDict = modelsParams[model][i]
			if tmpDict['GGPe'] >= maxG:
				if tmpDict['GGPe'] > maxG:
					maxGind = []
				maxG = tmpDict['GGPe']
				maxGind.append(i)
			if tmpDict['THETA_STN'] >= maxTH:
				if tmpDict['THETA_STN'] > maxTH:
					maxTHind = []
				maxTH = tmpDict['THETA_STN']
				maxTHind.append(i)
		indices = [0,1,2,3]
		tmp[0][0] = modelsParams[model][[x for x in indices if not (x in maxTHind or x in maxGind)][0]]
		tmp[0][1] = modelsParams[model][[x for x in indices if x in maxGind and not x in maxTHind][0]]
		tmp[1][0] = modelsParams[model][[x for x in indices if x in maxTHind and not x in maxGind][0]]
		tmp[1][1] = modelsParams[model][[x for x in maxTHind if x in maxGind][0]]
		modelsParams[model] = tmp
		
		for N in NUCLEI:
			for foc in focus:
				os.makedirs("results/"+foc+"/model_"+str(model)+'/'+N)
				for M in METRICS:
					lines = []
					for x in range(len(tmp)):
						for y in range(len(tmp[x])):
							exec(open(tmp[x][y]['xpID']+"/data/"+N+'_PwSpec_'+M+'.py').read())
							

							ps['freqs'] = np.asarray(ps['freqs'])
							ps['power'] = np.asarray(ps['power'])
							freqs = re.findall(r'\d+', foc)
							posi_spectrum = np.where((ps['freqs']>int(freqs[0])) & (ps['freqs']<int(freqs[1])))
							#fPS = lowpass(ps['power'], 0.025, 0.01)
							#fPS = fPS[(len(fPS)-len(ps['power'])) // 2 : -(len(fPS)-len(ps['power'])) // 2]
							fPS = ps['power']
							lines.append(plt.plot(ps['freqs'][posi_spectrum], fPS[posi_spectrum], linewidth=1.4, color=linestyles[x][y][0], linestyle = linestyles[x][y][1]))
							
							'''
							fPS = lowpass(ps['power'], 0.02, 0.05)
							fPS = fPS[(len(fPS)-len(ps['power'])) // 2 : -(len(fPS)-len(ps['power'])) // 2]
							lines.append(plt.plot(ps['freqs'], fPS, linewidth=1.4, color=linestyles[x][y][0], linestyle = linestyles[x][y][1]))
							#lines.append(plt.plot(ps['freqs'], ps['power'], linewidth=.4, color=linestyles[x][y][0], linestyle = linestyles[x][y][1]))
							'''
					plt.xlabel('Freq. [Hz]')
					#plt.legend(handles=lines, ['Gain = 1, THETA = 0', 'Gain = 1.6, THETA = 0', 'Gain = 1, THETA = 10', 'Gain = 1.6, THETA = 10'])
					plt.savefig("results/"+foc+"/model_"+str(model)+'/'+N+'/'+M+'.png')
					plt.close()
					

def lowpass(series, fc, b):

	#fc = 0.25  # Cutoff frequency as a fraction of the sampling rate (in (0, 0.5)).
	#b = 0.2  # Transition band, as a fraction of the sampling rate (in (0, 0.5)).
	N = int(np.ceil((4 / b)))
	if not N % 2: N += 1  # Make sure that N is odd.
	n = np.arange(N)
	 
	# Compute sinc filter.
	h = np.sinc(2 * fc * (n - (N - 1) / 2.))
	 
	# Compute Blackman window.
	w = 0.42 - 0.5 * np.cos(2 * np.pi * n / (N - 1)) + 0.08 * np.cos(4 * np.pi * n / (N - 1))
	 
	# Multiply sinc filter with window.
	h = h * w
	 
	# Normalize to get unity gain.
	h = h / np.sum(h)

	return np.convolve(series, h)

if __name__ == '__main__':
	main()