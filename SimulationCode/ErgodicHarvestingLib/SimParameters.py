# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 09:32:15 2017

@author: Chen Chen
"""
import numpy as np

class SimParameters:
    res = 101      # Spatial Resolution of Ergodic Simulation, [# of points]
    tRes = 101     # Time Resolution of Ergodic Simulation, [# of samples per belief update]
    dt = 0.02      # Time step
    def __init__(self, dt):
        self.dt = dt
        
class ErgodicParameters(SimParameters):
    nIter = 14     # Number of iterations in ergodic optimization loop
    nFourier = 15  # Number of Fourier Coefficients used in Fisher Information EID
    # Cost Weights
    wBarrCost = 100  # Weight of the Barr Cost
    wErgCost  = 5    # Weight of the Ergodic Cost
    wControl  = 0.05 # Weight of the Control Cost
    wInitCtrl = 0    # Weight of the Initial Condition of the Control
    # ODE Integrater
    #   Using explicit runge-kutta method choose either dopri5 (5th order)
    #   or dop853 (8th order)
    odeIntegrator = 'dopri5'
    tRes = 101 # Time resolution for ODE
    # Armijo Line Search Parameters
    alpha = 0.4
    beta  = 0.025
    # Time horizon and delta time used in ergodic optimizations, including ndsolver
    timeHorizon = 1   

    def __init__(self):
        # Initialize Parameters
        self.dt = super().dt
        self.res = super().res
        # Initialize Time Horizon
        self.time = np.linspace(0, self.timeHorizon, self.tRes)
        
class EIDParameters(SimParameters):
    maxIter = 40  # Number of iterations
    SNR = 30      # SNR of the sensor measurement
    Sigma = 0.05  # Sigma of Sensor measurement model
    mass = 1      # Mass of the simulated sensor
    objAmp  = 0.10   # Amplitude of Object Oscillation
    objCenter = 0.5  # Center of Object Oscillation
    sInitPos = 0.4   # Initial position of the sensor
    maxT = 80        # Length of the simulation in seconds
    pLHistDepth = 1  # depth of the queue used to update likelihood function
    procNoiseSigma = 0.01 # Sigma of the process noise
    # Indexes of simulation iterations where the sensor will be blinded with randomized noise
    blindIdx = False
    # DeltaT for trajectory simulation
    dt = 0.02
    
    def __init__(self):
        self.res = super().res
        self.tRes = super().tRes
        self.maxIter = np.int(np.ceil(self.maxT / (self.dt * self.tRes)))
        
    def UpdateDeltaT(self, dt):
        self.dt = dt
        self.maxIter = np.int(np.ceil(self.maxT / (self.dt * self.tRes)))