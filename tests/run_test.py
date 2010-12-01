import numpy
import subprocess

def generate_data(states=3,dim=3,npoints=100) :    

    means = numpy.empty((states,dim))
    hic = 1.

    for i in range(states) :
        mu = numpy.zeros(dim)
        mu[(i+1)%dim] = hic
        means[i,:] = mu
        if (i+1)%dim == 0 :
            hic += 1

    f = open("test.txt",'w')
    for k in range(npoints):
        cm = numpy.random.randint(states)
        x_ = numpy.random.randn(dim)*0.2 + means[cm,:]
        for v in x_ :
            f.write("%f  "%v)
        f.write("\n")
    f.close()

def run_test(path) :
    print path
    p = subprocess.Popen(path, stdout = subprocess.PIPE)
    out,err = p.communicate()
    if(p.returncode != 0) :
#        print out
        exit(1)
    print "pass"

generate_data(6,4,10000)


#subprocess.call(["./test_smat"])
run_test(["./test_smat"])
run_test(["./test_gaussian"])
run_test(["./test_em","test.txt","6"])
run_test(["./test_kmeans","test.txt","6"])
run_test(["./test_regression"])
run_test(["./test_pdf"])
