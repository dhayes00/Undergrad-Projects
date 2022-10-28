# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:26:19 2022

@author: David Hayes (Python 3.8)

Description: This programme models Multi-Tone Continuous Wave (MTCW) Ranging. 
This was used in collaboration with another python programme which generated
and measured a physical ultrasonic MTCW, where the phase and amplitude of the
received wave were extracted and a distance measurement (with a sub-millimeter
accuracy) was obtained.
 
"""


import numpy as np
import math
from collections.abc import Iterable
import matplotlib.pyplot as plt
import time
startTime = time.time()
class Component:
    
    
    def __init__(self,*,f=40000, A=1.0, phi=0.0):
        self.f = f
        self.A = A
        self.phi = phi
        
    def __str__(self):
        return "{} amplitude, {} phase @ {} Hz".format(self.A, self.phi, self.f)
    
class Enviroment:
    
    def __init__(self,*,c=344):
        self.c = c
        
    def propagate(self, components):
        if isinstance(components, Iterable): #If there is more than 1 instance
            return [self.propagate1(component) for component in components]
        else:
            return self.propagate1(components)
    
    def propagate1(self, component):
        return component
        
        
class SimpleEnviroment(Enviroment):
    
    def __init__(self,*,c=344, d=0.0):
        super().__init__(c=c)
        self.d = d
        
    def propagate1(self, component):
        if self.d>0:
            deltaphi = 2*np.pi*(component.f/self.c)*self.d
            #print(deltaphi)
            return Component(f=component.f,
                             A = component.A/self.d,
                             phi = component.phi + deltaphi)
        else:
            return component
    
class NoisyEnviroment(SimpleEnviroment):
    
    def __init__(self,*,c=344, d=0.0, Anoise = 0.0, phinoise = 0.0):
        super().__init__(c=c, d=d)
        self.Anoise = Anoise
        self.phinoise = phinoise
        
    def propagate1(self, component):
        comp = super().propagate1(component)
        comp.A *= np.random.normal(1.0,self.Anoise)
        comp.phi += np.random.normal(0.0,self.phinoise)
        return comp

def print_components(components):
    if isinstance(components, Iterable):
        for i, component in enumerate(components):
            print(i,':',component)
    else:
        print(components)
        
        
class Measurement:
    
    def __init__(self,*,N,fs,t0=0.0, snoise = 0.0):
        self.fs = fs
        self.t0 = t0
        self.N = N
        self.snoise = snoise
        
    def measure(self,components):
        if isinstance(components, Iterable):
            total = self.measure1(components[0])
            for component in components[1:]:
                another = self.measure1(component)
                total.samples += another.samples
            return total
        else:
            return self.measure1(components)
        
    def measure1(self,comp):
        t = np.arange(self.N)/self.fs + self.t0 # Times signal is sampled
        samples =( comp.A*np.cos(-2*np.pi*comp.f*t + comp.phi)
                + np.random.normal(0,self.snoise,self.N))
        return SampledSignal(samples,fs=self.fs,t0=self.t0)
        
class SampledSignal:
    
    def __init__(self,samples,*,fs,t0 = 0.0):
        self.samples = samples
        self.fs = fs
        self.t0 = t0
        
    def times(self):
        t = np.arange(len(self.samples))/self.fs + self.t0 # Times signal is sampled
        return t
    
    def __str__(self):
        return "Signal with {} samples @ {} sps ; t0 = {} s".format(len(self.samples)
                , self.fs, self.t0)
        
    
def analyse_signal(f, sampledsignal):
    
    N = len(sampledsignal.samples)
    t = np.arange(N)/sampledsignal.fs + sampledsignal.t0 # Times signal is sampled
    
    swave = np.sin(-2*np.pi*f*t)
    cwave = np.cos(-2*np.pi*f*t)
    
    ss=np.dot(swave,swave)
    sc=np.dot(swave,cwave)
    cs=np.dot(cwave,swave)#NB: same as sc
    sd=np.dot(swave,sampledsignal.samples)
    cc=np.dot(cwave,cwave)
    cd=np.dot(cwave,sampledsignal.samples)
    
    A=np.mat(((ss,sc) , (cs,cc)))
    a=np.asarray(A.I*np.mat(((sd,) , (cd,))))

    samp = a[0,0]
    camp = a[1,0]
    phase = math.atan2(-samp, camp)
    amp=np.sqrt(samp*samp+camp*camp)
    
    return amp, phase

#================================= def main ==================================
from scipy.stats import linregress

freqs = np.linspace(38e3,41e3,15) # Frequencies of multi-tone signal
# Initialising amplitude, phase & Unwraped phase
Amp, phase = np.empty(len(freqs)), np.empty(len(freqs)) 
phaseU = np.empty(len(freqs))  

x = 0.5 #Distance = 0.5 m
  
M = Measurement(N = 750_000, fs = 1_000_000, snoise = 0.0) # Measuring
    
#env = SimpleEnviroment(d=x) # Propagating through an enviroment with no noise
env = NoisyEnviroment(d=x,Anoise = 0.1, phinoise = 0.05) # With noise

for fs in range(len(freqs)): # Creating a loop to analyse diff freq each time            
    A = Component(f = freqs[fs], A = 3.5)
    
    sampled = M.measure(env.propagate(A)) 
    Amp[fs], phase[fs] = analyse_signal(freqs[fs],sampled)
   

phaseU[:] = np.unwrap(phase) # Unwrapping phase
   
S, y0, r, _, sigmaS = linregress(freqs[:],phaseU[:]) # Model fitting 

distance = (S*(344/(2*np.pi))) # Calculating distance 

print("\nDistance =", distance*(10**2),'cm\n')


# Displaying Multi-tone signal
t = np.linspace(0,120*2.5e-5,1500)

signal = 0
for l in range(len(freqs)):
    
    signal = signal + Amp[l]*np.sin(2*np.pi*freqs[l]*t)

plt.figure(figsize = (8,6), dpi=1000, facecolor = 'black')

plt.plot(t, signal,'orange')
ax = plt.axes()
ax.set_facecolor('Black')
ax.spines['bottom'].set_color('white')
ax.spines['left'].set_color('white')

ax.xaxis.label.set_color('white')
ax.yaxis.label.set_color('white')
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
plt.xlabel('time [s]'), plt.ylabel('Amplitude [V]')

plt.grid(True)

executeTime = (time.time() - startTime)
print('Execution time :',str(executeTime),'\n')