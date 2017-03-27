import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import spline

r = np.arange(0, 1.125, 0.125)
theta = 2 * np.pi *r

#print theta

#maxCSN = 10.
maxCSN = 4.
restCtx = 2. * np.ones((9))
activityLevels = np.ones((9))
for i in range(9):
  activityLevels[i] = 2 * np.pi / 8 * i
activityLevels = 2. + (maxCSN-2.) * ( (np.cos(activityLevels)+1)/2.)

#print activityLevels

#restGPi = 72. * np.ones((9))
#GPi = np.array([3.500000 , 14.142857 , 55.071429 , 80.285714 , 84.857143 , 84.000000 , 48.714286 , 12.500000, 3.5 ])
restGPi = 71.8 * np.ones((9))
GPi = np.array([0.0, 1.07, 32.4, 81.2, 89.7, 79.1, 32.4, 0.0, 0.0])

# attemps to smooth the graph:
thetanew = np.linspace(0,2*np.pi,100)
ALSmooth = spline(theta, activityLevels, thetanew)
GPISmooth = spline(theta, GPi, thetanew)
rGPISmooth = spline(theta, restGPi, thetanew)
rCTXSmooth = spline(theta, restCtx, thetanew)

# --- Plotting ---
ax1 = plt.subplot(121,projection = 'polar')
#ax1.plot(theta,activityLevels)
ax1.plot(thetanew,ALSmooth)
#ax1.plot(theta,restCtx)
ax1.plot(thetanew,rCTXSmooth)
ax1.set_rticks([2, 4, 6, 8, 10])
ax1.set_title("CSN input")

ax2 = plt.subplot(122,projection = 'polar')
#ax2.plot(theta,GPi)
ax2.plot(thetanew,GPISmooth)
#ax2.plot(theta,restGPi)
ax2.plot(thetanew,rGPISmooth)
ax2.set_rticks([20, 40, 60, 80, 100])
ax2.set_title("GPi output")

plt.show()
