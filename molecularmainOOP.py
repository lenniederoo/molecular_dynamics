from molecularfunctionsOOP import particleClass
from molecularPhysicalQuantities import PlotPQs
import numpy as np
from JosPlotPy import AnimatedScatter
import molecularPhysicalQuantities as PQ

print "HOI :D here goes the main program"
matplotlib.pyplot.close("all") #closing all the figures

#set global constants
Np=108
deltat=.004
mass = 1
dens = 0.85
temp = 0.9
amountoftimesteps=1000

particles = particleClass(Np, dens, temp,mass)
plots=PlotPQs(particles,amountoftimesteps,deltat)
plots.PlotThings(particles)

print ":)"
print PQ.calc_Corr_velocity(particles,5,amountoftimesteps,deltat)
PQ.plotcorr(particles,0,100,20,amountoftimesteps,deltat)
Animation=AnimatedScatter(particles,deltat)
