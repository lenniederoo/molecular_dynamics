from molecularfunctionsOOP import particleClass
import molecularPhysicalQuantities as PQ
import numpy as np
#import matplotlib.pyplot as plt
#import time
from JosPlotPy import AnimatedScatter

print "HOI :D"
#set global constants
Np=108
deltat=.004
mass = 1
dens = 0.85
temp = 0.9
amountoftimesteps=1000

particles = particleClass(Np, dens, temp,mass)
particles.changeForces()
energies= np.zeros((amountoftimesteps+1,1),dtype = float)
momenta=np.zeros((amountoftimesteps+1,3),dtype = float)
Ekin=np.zeros((amountoftimesteps+1,1),dtype = float)
temperature=np.zeros((amountoftimesteps+1,1),dtype = float)
pot=np.zeros((amountoftimesteps+1,1),dtype = float)

energies[0] = particles.checkEnergy ()
#momenta[0] =particles.checkMomenta()
pot[0] =particles.checkPotential()
Ekin[0]=particles.checkKinEnergy()
temperature[0]=PQ.calc_Temp(particles)

for i in xrange(amountoftimesteps):
  #starttime= time.time()

  particles.update(deltat)
  if i%20 ==0:
    PQ.Temp_correction(particles,temp)

  energies[i+1]=particles.checkEnergy ()
  #momenta[i+1]=particles.checkMomenta()
  Ekin[i+1]=particles.checkKinEnergy()
  pot[i+1] =particles.checkPotential()
  temperature[i+1]=PQ.calc_Temp(particles)
  
  #timeremaining= (time.time()-starttime)*(amountoftimesteps-i)
  #print "time left: ", timeremaining
#
#t= np.arange(amountoftimesteps+1)
#targetarray=np.ones((amountoftimesteps+1,1),dtype = float)*temp
#figure()
#title('Simulation of %s particles'%(Np))
#subplot(141)
#title('Kinetic energy')
#plot(t,Ekin)
#subplot(142)
#plot(t, energies)
#title('total energy')
#subplot(143)
#plot(t, pot)
#title('potential energy')
#subplot(144)
##title('momenta')
##plot(t,momenta)
#plot(t, temperature)
#plot(t,targetarray)
#title('temperature')
#print ":)"
#show()
print ":)"
Animation=AnimatedScatter(particles,deltat)

