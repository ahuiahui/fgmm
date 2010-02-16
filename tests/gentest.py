import numpy


STATES = 3

#means = numpy.random.rand(STATES,3)
means = numpy.array([[1.,0.,0.],
                  [0.,1.,0.],
                  [0.,0.,1.]])
print means
f = open("test.txt",'w')

for k in range(10000):
    cm = numpy.random.randint(STATES)
    x_ = numpy.random.randn(3)*0.2 + means[cm,:]
    f.write("%f  %f  %f \n"%(x_[0],x_[1],x_[2]))
f.close()
