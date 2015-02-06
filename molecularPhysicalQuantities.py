import numpy as np
import math
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

