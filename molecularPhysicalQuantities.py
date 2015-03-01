import numpy as np
import math
import matplotlib.pyplot as plt
import f90press

      
def calc_av_Mom_sq(momenta,Np):
  av_Mom=np.sum(momenta**2)/Np
  return av_Mom

def calc_Temp(particles):
  Temp=1/3.0*calc_av_Mom_sq(particles.momenta,particles.Np)/particles.mass
  return Temp

def Temp_correction(particles,target):
  temp=calc_Temp(particles)
  labda=math.sqrt(target/temp)
  particles.momenta=labda*particles.momenta
  temp_new=calc_Temp(particles)
  #print 'target',target
  #print temp, 'temp oud'
  #print temp_new, 'temp nieuw'
  
def calc_Press(particles,Ekin):
    press=f90press.calc_press(particles.positions,particles.L, [particles.Np])
    press=particles.dens/(3*particles.Np)*(2*Ekin+press)
    return press
    
def calc_Corr_velocity(particles,tau,n_t,deltat):
  corr=0
  for j in xrange(0,n_t):
    for i in xrange(0,particles.Np):
      corr+=np.dot(particles.momenta[i],particles.momenta[i-tau])*deltat
    particles.update(deltat)
  corr=corr/(particles.Np*n_t*deltat)
  return corr
  
def plotcorr(particles,inittau,endtau,amountoftau,n_t,deltat):
  tau=np.linspace(inittau,endtau,amountoftau)
  corr=np.zeros((amountoftau,1),dtype = float)
  for i in xrange(0,amountoftau):
    corr[i]=calc_Corr_velocity(particles,tau[i],n_t,deltat)
  plt.figure()
  plt.plot(tau*deltat,corr,marker='o',linestyle='--')
  plt.title('Correlation')
  
def calc_SpecificHeat(particles,Ekin):
  #print error_calc(Ekin,20)
  Cv=((3*particles.Np)/2)/(1-(3*particles.Np*np.var(Ekin))/(2*Ekin**2))
  #print Cv
  return Cv
  
  
def error_calc(quantity,number_of_blocks):
    error=np.ones(np.shape(quantity),dtype=float)
    sigma=np.zeros(np.shape(quantity),dtype=float)
    #number_of_blocks=20
    block=np.size(quantity,0)/number_of_blocks
    #print 'block',block
    for i in xrange(0,number_of_blocks):
        #sigma[i*block:(i+1)*block]=np.mean(quantity[i*block:(i+1)*block]**2,axis=0)-np.mean(quantity[i*block:(i+1)*block],axis=0)**2
        sigma[i*block:(i+1)*block]=np.var(quantity[i*block:(i+1)*block],axis=0)    
        #print np.mean(quantity[i*block:(i+1)*block]**2),'mean of squared'
        #print np.mean(quantity[i*block:(i+1)*block])**2, 'squared of mean'
    #sigma=np.mean(quantity**2,axis=0)-np.mean(quantity,axis=0)**2    
    error=error*sigma/np.sqrt(block)
    #print 'error',error
    #print 'sigma',sigma
    #print quantity, 'quantity'
    return error
  
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
    self.krachten=np.zeros((self.n_t+1,3),dtype = float)
    self.krachten[0]=particles.checkforces()
  
  def PlotThings(self,particles,deltat):
    press=np.zeros((self.n_t+1,1),dtype = float)
    pressav=50
    for i in xrange(0,self.n_t):  
      particles.update(self.deltat)
      self.energies[i+1]=particles.checkEnergy ()
      self.mom[i+1]=particles.checkMomenta()
      self.Ekin[i+1]=particles.checkKinEnergy()
      self.pot[i+1] =particles.checkPotential()
      self.temperature[i+1]=calc_Temp(particles)  
      self.pressure[i+1]=calc_Press(particles,self.Ekin[i+1])
      self.krachten[i+1]=particles.checkforces()
      if i>pressav:
        press[i]=np.sum(self.pressure[i-pressav:i])/pressav
    press[0:pressav]=self.pressure[0:pressav]
    t= np.arange(self.n_t+1)
    targetarray=np.ones((self.n_t+1,1),dtype = float)*self.target
    #presserror=error_calc(press)
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
    #plt.plot(t, self.temperature)
    plt.errorbar(t,np.reshape(self.temperature,(self.temperature.shape[0], )),marker='o',linestyle='',yerr=np.reshape(error_calc(self.temperature,20),(self.temperature.shape[0], )))
    plt.plot(t,targetarray)
    plt.title('temperature')
    plt.show()
    plt.figure()
    #plt.plot(t,press,marker='o',linestyle='--')
    #print error_calc(press),'errorarray'
    #print np.reshape(error_calc(press),(press.shape[0], )),' npreshape' 
    plt.errorbar(t,np.reshape(self.pressure,(self.pressure.shape[0], )),marker='o',linestyle='',yerr=np.reshape(error_calc(self.pressure,20),(press.shape[0], )))
    plt.title('pressure')
    plt.show()
    plt.figure()
    plt.plot(t[0:self.n_t],calc_SpecificHeat(particles,self.Ekin[0:self.n_t]),marker='o',linestyle='--')
    plt.title('Specific Heat by constant Volume')
    plt.show()
    plt.figure()
    plt.plot(t,self.krachten)
    plt.title('forces')
    plt.show()
    return press