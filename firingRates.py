from pylab import *
from numpy import random as rnd

xmin = 0.0
xmax = 0.100

ymin = 0.0
ymax = 300.

tau_m = 6E-3
tau_ref = 2E-3
th = 26E-3
sigma = 7E-3

S_max = 300.

def LG14(x):
  return S_max / (1 + exp((th-x)/3.8E-3) )

def iaf_IF(x,refr,membr,thres):
  if x > thres:
    return 1 / (refr + membr * log(x / (x-thres)))
  else:
    return 0

th_distrib = rnd.normal(th,sigma,500)

def popAvg(x,refr,membr):
  y =0
  for i in range(len(th_distrib)):
    y += iaf_IF(x,refr,membr,th_distrib[i])
  return y / float(len(th_distrib))


t = arange(xmin, xmax, (xmax-xmin)/300.)

ax = subplot(111)

# en gros, t_ref controle S_max
# et t_m, la pente au point d'inflexion

data=[]
for i in t:
  data.append(popAvg(i,3.33E-3,.1E-3))

ax.plot(t, data,'g-',label='avg iaf response')

data=[]
for i in t:
  data.append(popAvg(i,2E-3,1.E-3))

ax.plot(t, data,'k-')

data=[]
for i in t:
  data.append(popAvg(i,tau_ref,tau_m))

ax.plot(t, data,'b-')

#ax.plot(t, data,'k-',label='avg response '+str(len(th_distrib))+' iaf (t_r='+str(3.33)+' ; t_m='+str(0.3))

#ax.plot(t, STNavg(t,1E-3,tau_m,th),'g-')

#ax.plot(t, STNavg(t,0,tau_m,th),'g-')

#ax.plot(t, STNavg(t,tau_ref,3E-3,th),'-')

data=[]
for i in t:
  data.append(iaf_IF(i,tau_ref,tau_m,th))
ax.plot(t, data,'r-')

data=[]
for i in t:
  data.append(iaf_IF(i,tau_ref,tau_m,th+sigma))
ax.plot(t, data,'r:')

data=[]
for i in t:
  data.append(iaf_IF(i,tau_ref,tau_m,th-sigma))
ax.plot(t, data,'r:')

ax.plot(t, LG14(t), 'b-')

ax.grid(True)

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

#ax.fill_between(t, 0, function(t))

#plt.savefig('STN_IF.png')
show()
