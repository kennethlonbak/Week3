import pylab as py

data = py.loadtxt("TIMING.dat")

fig, ax = py.subplots(1,1)
ax.plot(data[:,1]/(1e6),data[:,0])
#ax.set_xlim(0,0.0003)

py.show()