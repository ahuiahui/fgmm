import numpy

def generate_data(states=3,dim=3,npoints=10000) :
    
    means = numpy.empty((states,dim))
    hic = 1.

    for i in range(states) :
        mu = numpy.zeros(dim)
        mu[(i+1)%dim] = hic
        means[i,:] = mu
        if (i+1)%dim == 0 :
            hic += 1

    print means
    f = open("test.txt",'w')
    for k in range(npoints):
        cm = numpy.random.randint(states)
        x_ = numpy.random.randn(dim)*0.2 + means[cm,:]
        for v in x_ :
            f.write("%f  "%v)
        f.write("\n")
    f.close()

generate_data(3,2)
