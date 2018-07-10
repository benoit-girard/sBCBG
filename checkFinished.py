'''
The purpose of this file is to count the number of finished and failed simulations of the running simulation sequence
'''

import os
import subprocess
from shutil import copyfile

def main():

	#number of simulations
	nbSim = 5400

	nbFinished = 0
	nbBug = 0
	ID = []
	for i in range(nbSim):
		root = '2018_07_03_13_53_27_572185_xp'
		os.chdir(root+('%06d' % i))
		if os.path.isdir('report'):
			nbFinished += 1
		if os.path.isfile('err'):
			if not open('err', 'r').read() == '':
				nbBug += 1
		os.chdir('..')
	print nbFinished, 'simulations finished over', nbSim
	print nbBug, 'simulations bugged over', nbSim

if __name__ == '__main__':
    main()