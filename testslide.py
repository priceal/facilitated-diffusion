"""
D is diffusion coefficient in 1D
k = off rate. 
These are used in scaling output of slide()

"""

D = 100.0
k = 0.05

numsamples = 20000
nsteps = 100

sam = []
for i in range(numsamples):
    delx, delt = 0.,0.
    for n in range(nsteps):
        dx, dt = fd.slide()
        delx += dx*sqrt(2*D/k)
        delt += dt/k
    sam.append((delx,delt))
sam = array(sam)

z = sam[:,0] # scale distance 
t = sam[:,1]          # scale time

figure(1)
scatter(sam[:,0],sam[:,1])

figure(2)

histtimes = histogram(t,bins=200)
#histnorm = histsam[0]/histsam[0].max()
thist = 0.5*(histtimes[1][:-1]+histtimes[1][1:])
plot(thist,histtimes[0],'o')

figure(3)
histz = histogram(z,bins=200)
#histnorm = histsam[0]/histsam[0].max()
zhist = 0.5*(histz[1][:-1]+histz[1][1:])
plot(zhist,histz[0],'o')

print(0.5*(z*z/t).mean(), 1.0/t.mean())