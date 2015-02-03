import numpy as np
import math

T = 1.0
dens = 0.9
sigma=1
epsilon=1

print "hello world!"

def fill_init_pos(n, positions):
  m = int((n/4)**(1.0/3)+0.01) #amount of unit cells per direction
  print "m: ",m
  a=(4.0/n)**(1./3)
  #do things
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
  return positions

def fill_init_mom(n,momenta):
  mean,standdev=0,5
  momenta[:,0]=np.random.normal(mean,standdev,n)
  momenta[:,1]=np.random.normal(mean,standdev,n)
  momenta[:,2]=np.random.normal(mean,standdev,n)
  return momenta

class particleClass:
  def __init__(self, n):
    self.positions=np.zeros((n,3),dtype=float)
    self.momenta=np.zeros((n,3),dtype=float)
    self.forces=np.zeros((n,3),dtype=float)
    self.positions = fill_init_pos(n, self.positions)
    self.momenta = fill_init_mom(n, self.momenta)
    self.forces=sum_forces(self.positions,self.forces,n)
  def show(self):
    print "Positions: ", self.positions
    print "Momenta: ", self.momenta
    print "Forces: ",self.forces

def calc_force(particle1,particle2):
  deltax=particle2[0]-particle1[0]
  deltay=particle2[1]-particle1[1]
  deltaz=particle2[2]-particle1[2]
  r=math.sqrt(deltax**2+deltay**2+deltaz**2)
  F=4*epsilon*((12*sigma**12)/r**13 - (6*sigma**6)/r**7)
  Fx=F*deltax/r
  Fy=F*deltay/r
  Fz=F*deltaz/r
  return [Fx,Fy,Fz]

def sum_forces(positions,forces,M):
	for i in xrange(0,M):
		for j in xrange(0,M):
			if i!=j:
				forces[i]+=calc_force(positions[i],positions[j]) 
 	return forces
def plotthings(positions)
  
particles = particleClass(864)
particles.show()
print calc_force(particles.positions[0],particles.positions[1])

  
