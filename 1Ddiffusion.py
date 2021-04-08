"""
D is diffusion coefficient in 1D
k = off rate, or the probability per unit time of dissociation
These are used in scaling output of slide()

"""

# 
D1 = 1.0
k = 10.0
D3 = 100.0

#
a = 1.0
R = 100.0
ro = 1.1

#
numsamples = 20000
numsteps = 500

#################################################################
runfile('init.py', current_namespace=True)

# draw samples from sliding distribution and store in samples 
t, z = [], []
for i in range(numsamples):

    delz, delt = 0.0, 0.0
    for n in range(numsteps):    
        
        # first, calculate the sliding displacement
        dz_slide, dt_slide = fd.slide()
    
        # now calculate the hopping displacement        
        x, y, dz_hop, dt_hop = fd.hop( ro/a, R/a )  # call it once
        while (x*x+y*y) > 1: # repeat if needed until captured
            x, y, dz_hop, dt_hop = fd.hop( ro/a, R/a )  # needs scaled inputs
        
        # scale the displacements and add to total displacement
        delt += dt_slide/k + dt_hop*a*a/D3/2.0    
        delz += dz_slide*math.sqrt(2.0*D1/k) + dz_hop*a

    # store total displacements in lists
    t.append(delt)
    z.append(delz)
        
z, t = np.array(z), np.array(t)

plt.figure(1)
plt.hist(t,bins=80)
plt.title('total time spent in 1D diffusion')

plt.figure(2)
plt.hist(z,bins=80)
plt.title('total displacement from starting position')
print("calculated effective coefficient for {} hopping/sliding cycles:".format(numsteps))
print("D_eff = {:3.5f}".format(0.5*(z*z/t).mean()))


