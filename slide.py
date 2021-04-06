"""
D is diffusion coefficient in 1D
k = off rate, or the probability per unit time of dissociation
These are used in scaling output of slide()

"""
D = 1.0
k = 1.0
numsamples = 2000

#################################################################
# import modules if not imported already
try:
    fd
    plt
    np
except NameError:
    import facilitatedDiffusion as fd
    import pylab as plt
    import numpy as np
    
# draw samples from sliding distribution and store in samples 
samples = []
for i in range(numsamples):
    dx, dt = fd.slide()
    samples.append((dx,dt))
samples = np.array(samples)

# scale displacements and residence times and store in z, t
z = samples[:,0] * np.sqrt(2*D/k)
t = samples[:,1] / k      

plt.figure(1)
plt.hist(t,bins=80)

plt.figure(2)
plt.hist(z,bins=80)

print("calculated D = {}, k = {}".format(0.5*(z*z/t).mean(), 1.0/t.mean()))