# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 17:56:17 2017

@author: Chen Chen
"""
import numpy as np
from numpy.random import random, randn
from scipy.stats import norm
from scipy.signal import convolve#, firwin, filtfilt
from scipy.interpolate import interp1d
from scipy.io import loadmat
from ErgodicHarvestingLib.EntropyEID import EntropyEID
#import matplotlib.pyplot as plt

class EID(object):
    def __init__(self, eidParam):
        self.dt = eidParam.dt
        self.it = 0
        self.SNR = eidParam.SNR
        self.eidParam = eidParam
        
        self.pLHistDepth = eidParam.pLHistDepth
        self.pLPosQueue = np.array([])
        self.pLVQueue = np.array([])
        self.pLSigmaAmpBayesian = eidParam.pLSigmaAmpBayesian
        self.pLSigmaAmpEID = eidParam.pLSigmaAmpEID
        
        self.eidType = eidParam.eidType
        
        # State Transition Probability
        if eidParam.res % 2:
            self.pStateTransition = norm.pdf(np.linspace(0.0,1.0,eidParam.res), 0.5, eidParam.procNoiseSigma)
            self.pStateTransition = self.pStateTransition.reshape(eidParam.res,1)
        else:
            # Make sure it's symmetric
            self.pStateTransition = norm.pdf(np.linspace(0.0,1.0,eidParam.res+1), 0.5, eidParam.procNoiseSigma)
            self.pStateTransition = self.pStateTransition.reshape(eidParam.res+1,1)
        # Normalize
        self.pStateTransition /= np.sum(self.pStateTransition)

        # Internal Objects
        self.sensor = self.Sensor(snr=self.eidParam.SNR, sigma=self.eidParam.Sigma, res=self.eidParam.res, 
                                  maxT=self.eidParam.maxT, objAmp=self.eidParam.objAmp, objCenter=self.eidParam.objCenter, 
                                  rawTraj=self.eidParam.rawTraj)


    def UpdateTraj(self, traj, pLast):
        self.max_iter = len(traj)
        
        # Internal Objects
        self.sensor = self.Sensor(snr=self.eidParam.SNR, sigma=self.eidParam.Sigma, res=self.eidParam.res, 
                                  maxT=self.eidParam.maxT, objAmp=self.eidParam.objAmp, objCenter=self.eidParam.objCenter, 
                                  rawTraj=self.eidParam.rawTraj)
        
        # Sensor Position, shape = (max_iter, 1)
        self.sPos = traj
        
        # Sensor Belief, shape = (SpatialResol, max_iter)
        self.pB = np.zeros([self.sensor.SpatialResol, self.max_iter])

        # Prior Belief
        self.pLast = pLast

        # EID, shape = (SpatialResol, max_iter)
        self.phi = np.zeros([self.sensor.SpatialResol, self.max_iter])


    def Simulate(self, extIter):
        # Determine sensor blind behavior
        if type(self.eidParam.blindIdx) != bool and \
            type(self.eidParam.blindIdx) != int and extIter in self.eidParam.blindIdx:
            self.sensor.blind = self.eidParam.blind
            self.pLast = np.ones([self.sensor.SpatialResol,1])
            print("Enter blind sensor mode at iter #{0} ({1:.2f})".format(extIter, extIter * self.max_iter * self.dt), end='')
        else:
            self.sensor.blind = 0
        
        objTraj = []
        for it, sPos in enumerate(self.sPos):
            # Update Time Frame
            time = it * self.dt + extIter * self.max_iter * self.dt
    
            # STEP #1 Take New Measurement
            objTraj.append( self.sensor.Measure(sPos, time) )
            
            # STEP #2 Update Belief
            self.BayesUpdate(self.sensor, sPos, it)
            
            # STEP #3 Calculate EID
            if self.eidType == "FI":
                self.CalcFIEID(self.sensor, self.pB[:, it], it)
            elif self.eidType == "Entropy":
                self.CalcEntropyEID(self.pB[:, it], it)
            else:
                # By default, use Fisher Information (EEDI)
                self.CalcFIEID(self.sensor, self.pB[:, it], it)
        
        # Return most recent EID, sensor trajectory
        return self.phi, self.pB, objTraj
    
    def CalcFIEID(self, sensor, pB, it):
        phi = convolve(sensor.FIPow2, pB, mode='same')
        # Normalize
        phi = phi / phi.sum()
        
        print("Var(EID) = {0}".format(np.var(phi)))
        
        self.phi[:, it] = phi
        
    def CalcEntropyEID(self, pB, it=None):
        sigmaL = self.sensor.noiseSigma * self.pLSigmaAmpEID
        eid = EntropyEID(prior=pB, wsRes=self.sensor.SpatialResol, sigmaM=self.sensor.Sigma, sigmaL=sigmaL)
        if it != None:
            self.phi[:, it] = eid
            return self.phi[:, it]
        else:
            return eid.reshape(len(eid),1)

    def PlanNextMove(self, eid, sPos, wsLattice, stepSize=0.01, planDepth=1):
        # Loop through possible destinations (left, right, still)
        # and choose the one maximizes the EID
        # (minimizes entropy)
        moveEID = []
        # Interpret EID through workspace
        eidInterp = interp1d(wsLattice.flatten(), eid.flatten(), kind='cubic')
        
        # Compute planning window
        moves = np.arange(sPos-planDepth*stepSize, sPos+planDepth*stepSize, stepSize)
        moves = np.delete(moves, np.argwhere(moves < 0))
        moves = np.delete(moves, np.argwhere(moves > 1))
        
        for _, move in enumerate(moves):          
            if move < 0.05:
                # set EID to 0 if the candidate move is out of the workspace
                moveEID.append(0)
            elif move > 0.95:
                # set EID to 0 if the candidate move is out of the workspace
                moveEID.append(0)
            else:
                moveEID.append(eidInterp(move))
        # Global InfoMax position within the planning window
        posBest = moves[np.argmax(moveEID)]
        # Next move given stepsize
        pos = sPos + np.sign(posBest-sPos) * stepSize
        if pos < 0.05:
            pos = 0.05
        elif pos > 0.95:
            pos = 0.95
        return np.array([pos])

    def BayesUpdate(self, sensor, sPos, it):
        # Calculate Current Belief
        self.Likelihood(sensor, sPos)

        # State Transition Probability
        self.pStateTransition = self.pStateTransition.reshape(self.pLast.shape)
        self.pLast = convolve(self.pLast, self.pStateTransition, mode='same')
        self.pLast /= np.sum(self.pLast)
        
        # Update Current Belief
        if self.sensor.blind == 'zero':
            pB = np.ones([self.sensor.SpatialResol,1])
        else:
            pB = self.pL * self.pLast
        
        # Normalize
        pB /= np.sum(pB)
        
        # Update Belief
        self.pB[:, it] = pB.flatten()
        
        # Update Prior Belief
        self.pLast = self.pB[:, it].reshape(self.pLast.shape)
        
        pass

    def Likelihood(self, sensor, sPos):
        sigma = sensor.noiseSigma * self.pLSigmaAmpBayesian
        
        # Update history queue
        if len(self.pLPosQueue) < self.pLHistDepth:
            self.pLPosQueue = np.append(self.pLPosQueue, sPos)
            self.pLVQueue = np.append(self.pLVQueue, sensor.V)
        else:
            self.pLPosQueue[:-1] = self.pLPosQueue[1:]
            self.pLPosQueue[-1] = sPos
            self.pLVQueue[:-1] = self.pLVQueue[1:]
            self.pLVQueue[-1] = sensor.V

        pL = np.ones_like(sensor.RefMM)
        
        for s, v in zip(self.pLPosQueue, self.pLVQueue):  
            likelihood = (1.0 / (np.sqrt(2*np.pi) * sigma)) * \
                 np.exp( -(sensor.MeasureModel(s)-v)**2 / (2.0 * sigma**2) )
            likelihood /= sum(likelihood) 
            
            pL *= likelihood
            pL /= np.sum(pL)

        # Normalize pL
        #pL = pL / np.sum(pL)
        self.pL = pL.reshape(self.pLast.shape)
        return pL

    # FIR Low Pass Filter
    def LPF(self, x, sensor):
        # Flatten input and record shape
        #xShape = x.shape
        x = x.flatten()

        # Phase-invariant filtering using filtfilt
        #return filtfilt(b=self.firH, a=1.0, x, padlen=sensor.SpatialResol-1).reshape(xShape)
        #return filtfilt(b=self.firH, a=1.0, x=x, padlen=len(x)-1).reshape(xShape)

    # Performance Evaluation
    class Perf(object):
        def __init__(self, max_iter):
            self.MaxBelief = np.zeros([max_iter, 1])
            self.Variance = np.zeros([max_iter, 1])
            self.MeanObjPos = np.zeros([max_iter, 1])

    class Sensor(object):
        def __init__(self, snr, sigma, res, maxT, objAmp, objCenter, rawTraj):
            self.SNR = snr
            self.SpatialResol = res # Samples per measure
            self.Scale = 1.0
            self.Sigma = sigma
            self.sampGrid = np.linspace(0.0,1.0,self.SpatialResol)
            
            self.maxT = maxT
            self.objAmp  = objAmp
            self.objCenter = objCenter

            self.RefMM = self.MeasureModel(0.5)
            self.FI = self.FI1D(self.RefMM)
            self.FIPow2 = self.FI**2

            # is sensor blinded? used for simulating global search behavior
            self.blind = 0

            # Calculate Measurement Noise Amplitude
            mmPower = (self.RefMM**2).sum() / self.RefMM.size
            self.noisePower = 10.0 * np.log10(mmPower) - self.SNR
            self.noisePower = 10.0 ** (self.noisePower / 10.0)
            self.noiseSigma = np.sqrt( self.noisePower )
            #print("Noise power = {0}".format(self.noisePower))
            self.NoiseAmp = self.RefMM.max() / (10.0**(self.SNR/20.0))

            # Current Sensor Measurement
            self.V = 0.0
            self.Vp = interp1d(self.sampGrid, self.RefMM)
        
            
            # Load Object Trajectory Data
            if type(rawTraj) != bool:
                self.rawTraj = rawTraj
                # Create interpolate object
                #.rawTraj = interp1d(np.linspace(0.0,self.maxT,num=rawTraj.size),rawTraj.flatten(),kind='cubic')
                #print('Pre-defined object trajectory loaded!')
                self.ReScaleTraj()
            else:
                self.rawTraj = False
            
        def ReScaleTraj(self):
            rawTraj = self.rawTraj
            
            # Normalize and remove offset
            rawTraj = rawTraj / rawTraj.max()
            rawTraj = rawTraj - rawTraj.mean()
            
            # Scale
            rawTraj = rawTraj * self.objAmp / np.abs(rawTraj).max()
            
            # Recenter
            rawTraj = rawTraj + self.objCenter
            
            # Create interpolate object
            self.rawTraj = interp1d(np.linspace(0.0,self.maxT,num=rawTraj.size),rawTraj.flatten(),kind='cubic')
            
        def MeasureModel(self, miu):
            self.mm = norm.pdf(self.sampGrid, miu * self.Scale, self.Sigma)
            self.mm /= max(self.mm)
            return self.mm

        def FI1D(self, x):
            fi = np.diff(np.insert(x,0,0))
            fi = np.abs(fi)
            fi = fi / fi.sum()
            return fi

        def ObjPos(self, t):
            if type(self.rawTraj) == bool: 
                self.objFreq = 0.05
                return self.objAmp * np.sin(2.0 * np.pi * self.objFreq * t) + self.objCenter
            else:
                return self.rawTraj(t)

        def Measure(self, sPos, t):
            noise = self.NoiseAmp * randn(1)
            #noise = self.NoiseAmp * normal(0.0, self.Sigma, self.SpatialResol)
            #V = self.RefMM + noise

            # Interpolate over sample grid
            Vp = self.Vp
            
            offset = self.ObjPos(t) - sPos + 0.5
            if offset > 1:
                offset = 1.0
            elif offset < 0:
                offset = 0.0     
                
            if self.blind == 'zero':
                # Sensor blinded, returing either all zero or random readings
                self.V = 0.0
            elif self.blind == 'noise':
                # Sensor blinded, returing random sampled measurement
                self.V = Vp(random(1))
            else:
                # normal sensor reading
                self.V = Vp(offset) + noise
                
            return self.ObjPos(t)