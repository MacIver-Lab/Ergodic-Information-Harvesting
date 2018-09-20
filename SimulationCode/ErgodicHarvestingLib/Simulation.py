# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:53:25 2017

@author: Chen Chen @ MacIver Lab
"""

# Import basic packages
import numpy as np
from scipy.interpolate import interp1d
# Import Ergodic packages
from ErgodicHarvestingLib.EID import EID
from ErgodicHarvestingLib.Ergodicity import Ergodicity
from ErgodicHarvestingLib.EID_Opt import ergoptimize
# Time
from time import strftime
# Save to .mat file
from ErgodicHarvestingLib.save2mat import save2mat
from time import sleep
# Entropy
from scipy.stats import entropy


def EIDSim(ergParam, eidParam, showPlot=True, showLivePlot=False, plt_UseBlit=True, showMsg=True):    
    # Reseed the numpy random module for deterministic behavior cross machines
    np.random.seed(eidParam.randSeed)

    # Initialize
    eid = EID(eidParam)
    erg = Ergodicity(ergParam)
    
    # Initial control input is zero 
    uinit = 0

    # Initialize
    if eidParam.simType == 'EH':
        eidList = np.ones([ergParam.res, eidParam.maxIter])
        eidList[:, 1] /= eidList[:, 1].sum()
        phi = eidList[:, 1].reshape(len(eidList[:, 1]), 1)
        sTrajList = np.array([eidParam.sInitPos])
        oTrajList = np.empty(0)
        phiList = np.ones([ergParam.res, ergParam.res, eidParam.maxIter])
        pBList = np.ones([ergParam.res, ergParam.res, eidParam.maxIter])
        ergList = np.empty(eidParam.maxIter)
        enpList = np.empty(eidParam.maxIter)
        pLast = np.ones([ergParam.res, 1])
    elif eidParam.simType == 'IF':
        eidParam.maxIter = np.int(np.round(eidParam.maxT / ergParam.dt))
        eidList = np.ones([eidParam.res, eidParam.maxIter])
        sTrajList = np.array([eidParam.sInitPos])
        oTrajList = np.empty(0)
        phiList = np.ones([eidParam.res, 1, eidParam.maxIter])
        pBList = np.ones([eidParam.res, 1, eidParam.maxIter])
        ergList = np.empty(0)
        enpList = np.empty(eidParam.maxIter)
        pLast = np.ones([eidParam.res, 1])
        # Workspace lattice
        wsLattice = np.linspace(0, 1, eidParam.res)
        pLast /= pLast.sum()
        phiList[:, :, 0] = eid.CalcEntropyEID(pLast)
    
    # Initialize Plot
    if showPlot:
        #plt.ion() # Turn on interactive mode
        plt.close('all')
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        ax1.set_xlim([0, 1])
        ax1.set_ylim([0, 1])
        ax1.set_xlabel('Time')
        ax1.set_title('Simulated Trajectory')
        ax2.set_xlim([0, 1])
        ax2.set_ylim([0, 1])
        ax2.set_xlabel('Location')
        ax2.set_ylabel('Normalized EER')
        ax2.set_title('Current EER')
        trajGrid = np.linspace(0, 1, eidParam.res)
        h_sTraj, = ax1.plot(trajGrid, trajGrid)
        h_oTraj, = ax1.plot(trajGrid, trajGrid)
        ax1.legend(['Sensor Trajectory', 'Object Trajectory'])
        h_oldBelief, = ax2.plot(trajGrid, trajGrid)
        h_newBelief, = ax2.plot(trajGrid, trajGrid)
        ax2.legend(['Prior EER', 'Current EER'])
        fig.canvas.draw()
        #plt.pause(0.00001)
        if plt_UseBlit:
            ax1_bg = fig.canvas.copy_from_bbox(ax1.bbox)
            ax2_bg = fig.canvas.copy_from_bbox(ax2.bbox)
    
    # Main Loop
    for it in range(eidParam.maxIter-1):
        if eidParam.simType == 'EH':
            # Use ergodic optimization to derive a ergodic trajectory
            [rawTraj, uinit] = ergoptimize(
                    phi[:, -1], sTrajList[-1], control_init=uinit, 
                    ergParam=ergParam, showStats=False, showMsg=showMsg)
            # [Optional] Filter the wiggle to evaluate causation
            traj = interp1d(ergParam.time, rawTraj)(ergParam.eidTime)
            # Update measurement and EID with current trajectory
            eid.UpdateTraj(traj, pLast)
            phi, pB, objTraj = eid.Simulate(it)
            # Update belief and EID map
            eidList[:, it] = phi[:, -1]
            phiList[:, :, it] = phi
            pBList[:, :, it] = pB
            pLast = pB[:, -1]
            # Evaluate ergodicity
            ergList[it] = erg.computeErgMeasure(traj, pLast)
        elif eidParam.simType == 'IF':
            # Plan for next move
            traj = eid.PlanNextMove(phiList[:, :, it], sTrajList[it], wsLattice, 
                                    stepSize=eidParam.stepSize,
                                    planDepth=eidParam.planDepth)
            # Update measurement and EID with current move
            eid.UpdateTraj(traj, pLast)
            phi, pB, objTraj = eid.Simulate(it)
            # Update belief and EID map
            eidList[:, it+1] = phi[:, -1]
            phiList[:, :, it+1] = phi
            pBList[:, :, it+1] = pB
            pLast = pB[:, -1]
            # Evaluate Entropy
            enpList[it] = entropy(pLast)
            
        # Update trajectory list
        sTrajList = np.append(sTrajList, traj)
        oTrajList = np.append(oTrajList, objTraj)

        if showLivePlot and showLivePlot:
            # Update plot
            h_sTraj.set_ydata(traj)
            h_oTraj.set_ydata(objTraj)
            h_oldBelief.set_ydata(eidList[:, it]/max(eidList[:, it]))
            h_newBelief.set_ydata(phi[:, -1]/max(phi[:, -1]))
            if plt_UseBlit:
                # Restore background
                fig.canvas.restore_region(ax1_bg)
                fig.canvas.restore_region(ax2_bg)
                # Redraw just the points
                ax1.draw_artist(h_sTraj)
                ax1.draw_artist(h_oTraj)
                ax2.draw_artist(h_oldBelief)
                ax2.draw_artist(h_newBelief)
                # Fill in the axes rectangle
                fig.canvas.blit(ax1.bbox)
                fig.canvas.blit(ax2.bbox)
                #plt.pause(0.00001)
            else:
                # Redraw everything
                fig.canvas.draw()
            display.clear_output(wait=True)
            display.display(fig)
            sleep(0.0001)
        if showMsg:
            print('*** Iteration #{0} ***'.format(it))
    if showPlot:
        # Plot the entire trajectory
        plt.clf()
        plt.plot(sTrajList)
        plt.plot(oTrajList)
        plt.draw()
    
    # Save data in MATLAB mat format
    matDict = {}
    matDict['oTrajList'] = oTrajList
    matDict['sTrajList'] = sTrajList
    matDict['gAtten'] = 0
    matDict['timeHorizon'] = ergParam.timeHorizon
    matDict['SNR'] = eidParam.SNR
    matDict['dt'] = eidParam.dt
    matDict['maxTime'] = eidParam.maxT
    matDict['wControl'] = ergParam.wControl
    matDict['sInitPos'] = eidParam.sInitPos
    matDict['objCenter'] = eidParam.objCenter
    matDict['simType'] = eidParam.simType
    matDict['eidType'] = eidParam.eidType
    matDict['randSeed'] = eidParam.randSeed
    matDict['pLSigmaAmp'] =  eidParam.pLSigmaAmp
    matDict['pLSigmaAmpEID'] = eidParam.pLSigmaAmpEID
    matDict['pLSigmaAmpBayesian'] = eidParam.pLSigmaAmpBayesian
    matDict['procNoiseSigma'] = eidParam.procNoiseSigma
    matDict['tRes'] = ergParam.tRes
    matDict['objAmp'] = eidParam.objAmp
    matDict['Sigma'] = eidParam.Sigma
    matDict['eidList'] = eidList
    matDict['pB'] = pBList
    matDict['phi'] = phiList
    matDict['blindIdx'] = eidParam.blindIdx
    matDict['blindType'] = eidParam.blind
    matDict['timeDone'] = strftime("%b-%d-%H-%M-%S")
    matDict['enpList'] = enpList
    if eidParam.simType == 'EH':
        matDict['ergList'] = ergList
    elif eidParam.simType == 'IF':
        matDict['StepSize'] = eidParam.stepSize
        matDict['planDepth'] = eidParam.planDepth
    
    save2mat(eidParam.filename,matDict, path=eidParam.saveDir, appendIdx=False)