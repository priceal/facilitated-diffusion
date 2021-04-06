# simulate dissociation for a protein that binds at zo
# and has an offrate = k
# units of time = L^2/(2D)


numsamples = 1000000
zo = 0.5
k = 10.0


residencetime = []
position = []
for i in range(numsamples):
    timesample = wp.randexp()/k
    residencetime.append(timesample)
    positionsample = wp.randpdf1d(zo,timesample)
    position.append(positionsample)
residencetime = array(residencetime)
position = array(position)

#hist(residencetime,bins=20)
#hist(position,bins=100)
histpos = histogram(position,bins=100)
normpos = histpos[0]/histpos[0].sum()
zhist = 0.5*(histpos[1][:-1]+histpos[1][1:])
plot(zhist,normpos,'-')