# -*- coding: utf-8 -*-
import numpy as np

from ErgodicHarvestingLib.ergodic import ErgodicOpt


def ergoptimize(pdf, state_init, control_init, ergParam, showStats=False, showMsg=True):

    state_init = np.array([state_init])
    solver = ErgodicOpt(1, 1, ergParam, control_init)
    solver.set_pdf(pdf)
    U0 = np.zeros([ergParam.tRes, 1])
    U0[0] = control_init
    X0 = solver.simulate(state_init, U0)
    solver.update_traj(X0, U0)
    # set up some containers
    costs = solver.evalcost()
    trajlist = [X0, U0]
    # Armijo line search parameters
    alpha = ergParam.alpha
    beta = ergParam.beta

    for k in range(ergParam.nIter):
        descdir = solver.descentdirection()
        newdcost = solver.dcost(descdir)
        gamma = 1.0
        newXU = [trajlist[0] + gamma * descdir[0], trajlist[1] + gamma * descdir[1]]
        newtraj = solver.project(state_init, newXU)
        solver.update_traj(newtraj[0], newtraj[1])
        newcost = solver.evalcost()

        while newcost > (costs + alpha * gamma * newdcost) and gamma > 1e-8:
            gamma = beta * gamma
            stepDirection = [
                trajlist[0] + gamma * descdir[0],
                trajlist[1] + gamma * descdir[1],
            ]
            newtraj = solver.project(state_init, stepDirection)
            solver.update_traj(newtraj[0], newtraj[1])
            newcost = solver.evalcost()

        trajlist = [solver.X_current, solver.U_current]
        if showStats:
            print("J = {0}, diff = {1}".format(newcost, costs - newcost))
        if (costs - newcost) < 1e-5 or k > 9:
            break
        costs = newcost

    if showStats or showMsg:
        print(f"Ergodic optimization loop - k = {k}")
    return [solver.X_current.flatten(), solver.U_current.flatten()[-1]]
