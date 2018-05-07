# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:24:51 2017

@author: Chen Chen
"""
from ErgodicInfotaxisAPI.ergodic import ErgodicOpt
import numpy as np

def ergoptimize(pdf, state_init,
                control_init,
                ergParam,
                showStats=False,
                showMsg=True):

    state_init = np.array([state_init])
    
    solver = ErgodicOpt(1, 1, ergParam, control_init)
    solver.set_pdf(pdf)
    
    
    U0 = np.zeros([ergParam.tRes, 1])
    U0[0] = control_init
    
    X0 = solver.simulate(state_init,U0)
    solver.update_traj(X0,U0)

    #set up some containers
    costs = solver.evalcost()
    trajlist = [X0,U0]
    
    # Armijo line search parameters
    alpha = ergParam.alpha
    beta  = ergParam.beta

    if showMsg:
        print("Ergodic optimization loop - k = ", end='')
        
    for k in range(ergParam.nIter):
        if showStats:
            print("{0}".format(k))
        else:
            if showMsg:
                print("{0} ".format(k), end='')
            
        descdir = solver.descentdirection()
        newdcost = solver.dcost(descdir)
        #print("Cost =", newdcost)
        gamma = 1.0
        newXU = [trajlist[0] + gamma*descdir[0], trajlist[1] + gamma*descdir[1]]
        newtraj = solver.project(state_init,newXU)
        solver.update_traj(newtraj[0],newtraj[1])
        newcost = solver.evalcost()
        
        while newcost > (costs + alpha*gamma*newdcost) and gamma>.00000001:
            gamma = beta*gamma
            #print(gamma)
            stepDirection = [trajlist[0]+gamma*descdir[0],trajlist[1]+gamma*descdir[1]]
            newtraj = solver.project(state_init,stepDirection)
            solver.update_traj(newtraj[0],newtraj[1])
            newcost = solver.evalcost()
            
        trajlist = [solver.X_current,solver.U_current]
        
        if showStats:
            print("J = {0}, diff = {1}".format(newcost, costs - newcost))
            
        if costs - newcost < 0.00001:
            break
        if k > 9:
            break
        costs = newcost
        
    return [solver.X_current.flatten(), solver.U_current.flatten()[-1]]

''' For testing
xgrid = np.linspace(0,1,100)
#pdf = norm.pdf(xgrid, 0.8, 0.1) + norm.pdf(xgrid, 0.4, 0.1)
pdf = norm.pdf(xgrid, 0.5, 0.1)
pdf /= pdf.sum()
[traj,uinit] = ergoptimize(pdf, 0.1, simRes=len(xgrid))
plt.plot(traj)
plt.plot(pdf/pdf.max())
plt.pause(1)
'''