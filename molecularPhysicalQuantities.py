import numpy as np
import math
import matplotlib.pyplot as plt

      
def calc_av_Mom_sq(momenta,Np):
  av_Mom=np.sum(momenta**2)/Np
  return av_Mom

def calc_Temp(particles):
  Temp=1/3.0*calc_av_Mom_sq(particles.momenta,particles.Np)
  return Temp

def Temp_correction(particles,target):
  temp=calc_Temp(particles)
  labda=math.sqrt(target/temp)
  particles.momenta=labda*particles.momenta
  
def calc_Press(particles,deltat):
  indices1=[((particles.positions[:,0] + particles.momenta[:,0]*(deltat/particles.mass))>particles.L)]
  indices2=[((particles.positions[:,1] + particles.momenta[:,1]*(deltat/particles.mass))>particles.L)]
  indices3=[((particles.positions[:,2] + particles.momenta[:,2]*(deltat/particles.mass))>particles.L)]
  #print indices
  press1=np.sum(particles.momenta[indices1,0])/(particles.L**2*deltat)
  press2=np.sum(particles.momenta[indices2,1])/(particles.L**2*deltat)
  press3=np.sum(particles.momenta[indices3,2])/(particles.L**2*deltat)
  press=(press1+press2+press3)/3
  return press
  
  
class PlotPQs:
  def __init__(self, particles,amountoftimesteps,deltat):
    self.target=particles.temp
    self.n_t=amountoftimesteps
    self.deltat=deltat
    self.energies= np.zeros((self.n_t+1,1),dtype = float)
    self.mom=np.zeros((self.n_t+1,3),dtype = float)
    self.Ekin=np.zeros((self.n_t+1,1),dtype = float)
    self.temperature=np.zeros((self.n_t+1,1),dtype = float)
    self.pot=np.zeros((self.n_t+1,1),dtype = float)
    self.energies[0] = particles.checkEnergy()
    self.mom[0] =particles.checkMomenta()
    self.pot[0] =particles.checkPotential()
    self.Ekin[0]=particles.checkKinEnergy()
    self.temperature[0]=calc_Temp(particles)
    self.pressure=np.zeros((self.n_t+1,1),dtype=float)
    self.pressure[0]=particles.temp/particles.L**3
  
  def PlotThings(self,particles):
    for i in xrange(0,self.n_t):  
      particles.update(self.deltat)
      self.energies[i+1]=particles.checkEnergy ()
      self.mom[i+1]=particles.checkMomenta()
      self.Ekin[i+1]=particles.checkKinEnergy()
      self.pot[i+1] =particles.checkPotential()
      self.temperature[i+1]=calc_Temp(particles)  
      self.pressure[i+1]=calc_Press(particles,self.deltat)
    t= np.arange(self.n_t+1)
    targetarray=np.ones((self.n_t+1,1),dtype = float)*self.target
    plt.figure()
    plt.subplot(131)
    plt.title('Kinetic energy')
    plt.plot(t,self.Ekin)
    plt.subplot(132)
    plt.plot(t, self.energies)
    plt.title('total energy')
    plt.subplot(133)
    plt.plot(t, self.pot)
    plt.title('potential energy')
    plt.show()
    plt.figure()
    plt.subplot(121)
    plt.title('momenta')
    plt.plot(t,self.mom)
    plt.subplot(122)
    plt.plot(t, self.temperature)
    plt.plot(t,targetarray)
    plt.title('temperature')
    plt.show()
    plt.figure()
    plt.plot(t,self.pressure)
    plt.title('pressure')
    plt.show()