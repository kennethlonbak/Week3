import pylab as py
data = {}
n_run = 2
n_rank = 2
for run in range(1,n_run):
    data[run] = {}
    for rank in range(0,n_rank,2):
        print("Loading rank:%d for run:%d"%(rank,run))
        data[run][rank] = py.loadtxt("TIMING_rank%d_run%02d.dat"%(rank,run))
nr = 1
#mean = py.zeros(shape=py.shape(data[0][0][nr:,0]))
n=0
fig, ax = py.subplots(1,1)
for run in range(1,n_run):
    for rank in range(0,n_rank,2):
        n += 1
        print("Plotting rank:%d for run:%d" % (rank, run))
        t_rmean = py.convolve(data[run][rank][:,0],py.ones([nr,1]))
        b_rmean = py.convolve(data[run][rank][:,1]/1e6,py.ones([nr,1]))
        ax.plot(b_rmean,t_rmean,label="rank%d run%d"%(rank,run))
        #mean += t_rmean

#ax.plot(b_rmean,mean/n,label='Mean')
#ax.set_xlim(0,0.0003)
ax.legend()
py.show()