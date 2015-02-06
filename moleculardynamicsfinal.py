import numpy as np
import matplotlib.pyplot as plt
import math 
import mpl_toolkits.mplot3d.axes3d as p3
# m = int((Np/4)**(1.0/3)+0.01) #amount of unit cells per direction
# density = 0.85

print "hello world!"

def fill_init_pos(N, density):
#  print "m: ",m
  m = int((N/4)**(1.0/3)+0.01) 
#amount of unit cells per direction  #do things
  positions = np.zeros((N,3),dtype=float)
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
  L = (N/density)**(1.0/3)
  positions = positions*L/m
  return positions, L

def fill_init_mom(N,temp):
  mean,standdev=0,np.sqrt(temp)
  momenta = np.zeros((N,3),dtype = float)
  momenta[:,0]=np.random.normal(mean,standdev,N)
  momenta[:,1]=np.random.normal(mean,standdev,N)
  momenta[:,2]=np.random.normal(mean,standdev,N)
  return momenta

class particleClass:
  def __init__(self, Np, density, temp):
 #   self.momenta=np.zeros((Np,3),dtype=float)
    self.forces=np.zeros((Np,3),dtype=float)
    self.positions, self.L = fill_init_pos(Np, density )
    self.momenta = fill_init_mom(Np, temp)
    self.forces=sum_forces(self.positions,self.forces,Np, self.L)
  def show(self):
    print "Positions: ", self.positions
    print "Momenta: ", self.momenta
    print "Forces: ",self.forces

def calc_force(particle1,particle2, L):
  deltax= particle2[0]-particle1[0]
  deltax=deltax-round(abs(deltax/L))*L
  deltay=particle2[1]-particle1[1]
  deltay=deltay-round(abs(deltay/L))*L
  deltaz=particle2[2]-particle1[2]
  deltaz=deltaz-round(abs(deltaz/L))*L
  r2=(deltax**2+deltay**2+deltaz**2)
  F=24*(2/r2**7 - 1/r2**4)
  if (F>10.0):
    print  r2, L
    #raw_input(' ')
  Fx=F*deltax
  Fy=F*deltay 
  Fz=F*deltaz
  return [Fx,Fy,Fz]

def sum_forces(positions,forces,Np, L):
  forces = np.zeros((Np,3), dtype=float)
  for i in xrange(0,Np):
    for j in xrange(i+1,Np):
      forcebetweenparticles=calc_force(positions[i],positions[j], L) 
      forces[i]-=forcebetweenparticles
      forces[j]+=forcebetweenparticles
  return forces

def plotthings(positions,Np):
    fig= plt.figure()
    ax = fig.gca(projection='3d')
    #ax=Axes3D(fig)
    #plotaxis(5)
    for i in range (0,Np):
        xs,ys,zs = positions[i]
        ax.scatter(xs,ys,zs)
    plt.show()
    
    
def changemom(deltat,momenta,forces):
    momenta += forces*deltat
    return momenta

def changepos(deltat,momenta,positions,mass, L):
#    for i in range(0,m):
    positions = (positions + momenta*(deltat/mass)) % L

    return positions

def checkPotential(particle1,particle2,L):
  deltar= particle1-particle2
  deltar[0]=deltar[0]-round(deltar[0]/L)*L
  deltar[1]=deltar[1]-round(deltar[1]/L)*L
  deltar[2]=deltar[2]-round(deltar[2]/L)*L
  r2=np.sum(deltar**2)
  V=4*((1/r2)**6 - (1/r2)**3)
  #print "V: ",V
  return V

def checkpotentialSUM(positions,Np,L):
  potentialSUM=0.0
  for i in xrange(0,Np):
    for j in xrange(0,Np):
      if j!=i:
        potentialSUM+=0.5*checkPotential(positions[i],positions[j], L)
  return potentialSUM

def checkEnergy(Np,momenta,positions,mass,L):
  Etot=np.sum(momenta**2)*0.5
  Etot+=checkpotentialSUM(positions,Np,L)
#  print 'epot', checkpotentialSUM(positions,Np,L)
#  for i in xrange(0,Np):
#    Etot+=sum(0.5*momenta[i]**2/mass)
#    Etot=Etot+math.sqrt(sum(forces[i]**2))
  return Etot

def checkMomenta(momenta):
  momtot=np.sum(momenta)
  return momtot
"""
particles = particleClass(Np)
particles.show()
plotthings(particles.positions,Np)

print calc_force(particles.positions[0],particles.positions[1])
"""
