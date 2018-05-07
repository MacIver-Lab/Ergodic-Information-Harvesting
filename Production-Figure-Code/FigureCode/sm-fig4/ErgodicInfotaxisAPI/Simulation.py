# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:53:25 2017

@author: MacIver Lab
"""

# Import basic packages
import numpy as np
import matplotlib.pyplot as plt
#from numpy.random import rand
# Import Ergodic packages
from ErgodicInfotaxisAPI.EID import EID
from ErgodicInfotaxisAPI.Ergodicity import Ergodicity
from ErgodicInfotaxisAPI.EID_Opt import ergoptimize
# Time
from time import strftime
# Save to .mat file
from ErgodicInfotaxisAPI.save2mat import save2mat
from IPython import display
from time import sleep



def EIDSim(ergParam, eidParam, showPlot=True, showLivePlot=False, plt_UseBlit=True, showMsg=True):
    # Initialize
    eid = EID(eidParam)
    erg = Ergodicity(ergParam)

    eidList = np.ones([ergParam.res, eidParam.maxIter])
    sTrajList = np.array([eidParam.sInitPos])
#    sTrajListRaw = np.array([eidParam.sInitPos])
    oTrajList = np.empty(0)
    uinit = 0
    phiList = np.ones([ergParam.res, ergParam.res, eidParam.maxIter])
    pBList = np.ones([ergParam.res, ergParam.res, eidParam.maxIter])
    ergList = np.empty(eidParam.maxIter)
    pLast = np.ones([ergParam.res, 1])

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
        ax2.set_ylabel('Normalized EID')
        ax2.set_title('Current EID')

        trajGrid = np.linspace(0, 1, ergParam.res)
        h_sTraj, = ax1.plot(trajGrid, trajGrid)
        h_oTraj, = ax1.plot(trajGrid, trajGrid)
        ax1.legend(['Sensor Trajectory', 'Object Trajectory'])

        h_oldBelief, = ax2.plot(trajGrid, trajGrid)
        h_newBelief, = ax2.plot(trajGrid, trajGrid)
        ax2.legend(['Prior EID', 'Current EID'])

        fig.canvas.draw()
        #plt.pause(0.00001)

        if plt_UseBlit:
            ax1_bg = fig.canvas.copy_from_bbox(ax1.bbox)
            ax2_bg = fig.canvas.copy_from_bbox(ax2.bbox)


    # Main Loop
    for it in range(eidParam.maxIter-1):
        # Use ergodic optimization to derive a ergodic trajectory
        [rawTraj, uinit] = ergoptimize(
                eidList[:, it], sTrajList[-1], control_init=uinit,
                ergParam=ergParam, showStats=False, showMsg=showMsg)

        # [Optional] Filter the wiggle to evaluate causation
        if ergParam.filterWiggle and it >= ergParam.filterWiggleStartIdx:
            traj = ergParam.FilterWiggle(rawTraj)
        else:
            traj = rawTraj

        # Update measurement and EID with current trajectory
        eid.UpdateTraj(traj, pLast)
        phi, pB, objTraj = eid.Simulate(it)

        # Update belief and EID map
        eidList[:, it+1] = phi[:, -1]
        phiList[:, :, it+1] = phi
        pBList[:, :, it+1] = pB
        pLast = pB[:, -1]

        # Evaluate ergodicity
        ergList[it] = erg.computeErgMeasure(traj, phi[:, -1])

        # Update trajectory list
        sTrajList = np.append(sTrajList, traj)
#        sTrajListRaw = np.append(sTrajListRaw, rawTraj)
        oTrajList = np.append(oTrajList, objTraj)

        if showLivePlot and showLivePlot:
            # Update plot
            h_sTraj.set_ydata(traj)
            h_oTraj.set_ydata(objTraj)
            h_oldBelief.set_ydata(eidList[:, it]/max(eidList[:, it]))
            h_newBelief.set_ydata(eidList[:, it+1]/max(eidList[:, it+1]))

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

    # Save data to MATLAB
    matDict = {}
    matDict['oTrajList'] = oTrajList
    matDict['sTrajList'] = sTrajList
    matDict['gAtten'] = 0
    matDict['timeHorizon'] = ergParam.timeHorizon
    matDict['SNR'] = eidParam.SNR
    matDict['dt'] = eidParam.dt
    matDict['maxTime'] = eidParam.maxT
    matDict['wControl'] = ergParam.wControl
    matDict['objAmp'] = eidParam.objAmp
    matDict['Sigma'] = eidParam.Sigma
    matDict['eidList'] = eidList
    matDict['pB'] = pBList
    matDict['phi'] = phiList
    matDict['ergList'] = ergList
    matDict['blindIdx'] = eidParam.blindIdx
    matDict['blindType'] = eidParam.blind
    matDict['timeDone'] = strftime("%b-%d-%H-%M-%S")

    #filename = 'SNR-' + str(eidParam.SNR) + \
    #           '-Sigma-' + str(eidParam.Sigma) + \
    #           '-oAmp-' + str(eidParam.objAmp)+ \
    #           '-wCtrl-' + str(ergParam.wControl) + \
    #           '-dt-' + str(ergParam.dt) + \
    #           '-blind-' + eidParam.blind
           #'-timeH-' + str(ergParam.timeHorizon)

    save2mat(eidParam.filename,matDict, path=eidParam.saveDir)
