
import numpy as np
import fgmm
from matplotlib import pyplot

t = np.linspace(0,6,1000)
#t = np.arange(2000)
y = np.sin(t) + np.random.random(1000)*.5

pyplot.plot(t,y,'k.')

dat = np.vstack([t,y]).T

dat = dat.astype(np.dtype('f4')) # this is the messy stuff .. 

g = fgmm.GMM(20,2)

g.init(dat)

g.Em(dat)

g.Dump()

a=1

g.InitRegression(a)

inp = np.linspace(0,6,100) 
r = []
print "jahldfuablfjhbl"

for t in inp :                                                
    r.append(g.DoSamplingRegression(np.array([t],dtype=np.dtype('f4'))))

s=np.vstack([g.Draw() for _ in range(1000)])
pyplot.plot(inp,r,'r-',linewidth=3)
pyplot.plot(s[:,0],s[:,1],'b.')

print g.GetMean(1)
print g.GetCovariance(1)
pyplot.show()


