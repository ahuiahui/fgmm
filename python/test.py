# 
# test script for python interface to fgmm 
#
# should plot noisy sine data (as black dots) 
# and try sampling regression (red line ) 
# and sampling from the sine data ( blue dots ) 
# 
import numpy as np
import fgmm
from matplotlib import pyplot

nstates = 8

t = np.linspace(0,6,1000)
#t = np.hstack([t[:300], t[600:] ])
#print t.shape
y = np.sin(t) + np.random.random(1000)*.5

pyplot.plot(t,y,'k.')

dat = np.vstack([t,y]).T

dat = dat.astype(np.dtype('f4')) # this is the messy stuff .. 

g = fgmm.GMM(nstates,2)

g.kmeans(dat)

g.Em(dat,1e-3,fgmm.COVARIANCE_FULL)

g.Dump()

a=1

g.InitRegression(a)

inp = np.linspace(0,6,100) 
r = []
gmr = []
print "jahldfuablfjhbl"

for t in inp :                                                
    r.append(g.DoSamplingRegression(np.array([t])))
    gmr.append(g.DoRegression(np.array([t])))

s=np.vstack([g.Draw() for _ in range(1000)])
pyplot.plot(inp,r,'r-',linewidth=3)
pyplot.plot(inp,gmr,'g-',linewidth=5)
pyplot.plot(s[:,0],s[:,1],'b.')

_t = np.linspace(0,2*np.pi,50)
circle = np.vstack( [np.cos(_t), np.sin(_t)])
for i in range(nstates ) :
    sig = g.GetCovariance(i)
    print sig
    el = np.dot(np.linalg.cholesky(sig),circle).T
    print el.shape
    print g.GetMean(i)
    el += g.GetMean(i)
    pyplot.plot(el[:,0],el[:,1],'k-',linewidth=3)

pyplot.show()


