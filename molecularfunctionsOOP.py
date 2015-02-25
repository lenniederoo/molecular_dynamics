import numpy as np
import f90force
import f90pot
import molecularPhysicalQuantities as PQ

print "hello world!"

class particleClass:
  def __init__(self, Np, density, temp,mass):
    self.Np=Np
    self.temp=temp
    self.mass=mass
    self.forces=np.zeros((self.Np,3),dtype=float)
    self.positions, self.L = self.fill_init_pos(density)
    self.momenta = self.fill_init_mom(self.temp)
    self.changeForces()


  def show(self):
    print "Amount of particles: ",self.Np," box length: ",self.L
    print "Positions: ", self.positions
    print "Momenta: ", self.momenta
    print "Forces: ",self.forces
  
  def fill_init_pos(self,density):
    m = int((self.Np/4)**(1.0/3)+0.01) #amount of unit cells per direction  
    positions = np.zeros((self.Np,3),dtype=float)
    positions[0]=[0.0,0.0,0.0]
    positions[1]=[0.5,0.5,0.0]
    positions[2]=[0.5,0.0,0.5]
    positions[3]=[0.0,0.5,0.5]
    counter=0
    for i in xrange(0,m):
      for j in xrange(0,m):
        for k in xrange(0,m):
          positions[counter:counter+4]=positions[0:4]+[i,j,k]
          counter+=4 
    L = (self.Np/density)**(1.0/3)
    positions = positions*L/m
    return positions, L    

  def fill_init_mom(self,temp):
    mean,standdev=0,np.sqrt(temp)
    momenta = np.zeros((self.Np,3),dtype = float)
    momenta[:,0]=np.random.normal(mean,standdev,self.Np)
    momenta[:,1]=np.random.normal(mean,standdev,self.Np)
    momenta[:,2]=np.random.normal(mean,standdev,self.Np)
    return momenta

  def changeForces(self):
    # print f90force.calc_force.__doc__
    self.forces=f90force.calc_force(self.positions,self.forces,self.L, [self.Np])

  def changemom(self,deltat):
    self.momenta += self.forces*deltat
    
  def changepos(self,deltat):
    self.positions = (self.positions + self.momenta*(deltat/self.mass)) % self.L

  def update(self,deltat):
    self.changemom(deltat)
    self.changepos(deltat)
    self.changeForces()
    PQ.Temp_correction(self,self.temp)
    
  def checkMomenta(self):
    momtot=np.sum(self.momenta,axis=0)
    return momtot
  
  def checkPotential(self):
    totpot=f90pot.calc_pot(self.positions,self.L, [self.Np])
    return totpot
  
  def checkKinEnergy(self):
    Ekin=np.sum(self.momenta**2)*0.5/self.mass
    return Ekin

  def checkEnergy(self):
    Etot=self.checkKinEnergy()+self.checkPotential()
    return Etot    
