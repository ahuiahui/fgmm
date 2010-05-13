import subprocess
from matplotlib import pyplot

f = open("timings.txt",'w')
x,fg,MG = [],[],[]
for i in range(12) :
    npoints = 10**(i*.5+2)
    p = subprocess.Popen(["./bench","4","10",str(npoints)],stdout = subprocess.PIPE)
    out,_ = p.communicate()
    lines = out.splitlines()
    tfgmm = float(lines[1].split(' ')[0])
    tGMR = float(lines[3].split(' ')[0])

    x.append(npoints) 
    fg.append(tfgmm)
    MG.append(tGMR)

    f.write(" ".join(map(str,(npoints, tfgmm,tGMR))))
    f.write("\n")

f.close()

pyplot.semilogx(x,fg,'b',label="fgmm")
pyplot.semilogx(x,MG,'b',label="GMR c++ lib")
pyplot.title("time per EM iteration (4 states 10 D)")
pyplot.xlabel("number of data point")
pyplot.ylabel("[ms]")
pyplot.show()
