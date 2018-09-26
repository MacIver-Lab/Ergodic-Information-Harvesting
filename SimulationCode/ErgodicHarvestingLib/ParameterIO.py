import json
import numpy as np
from ErgodicHarvestingLib.SimParameters import ErgodicParameters, EIDParameters

def loadParams(dat):
    with open(dat, 'r') as f:
        pDict = json.load(f)

    # Parameter Objects
    ergParam = ErgodicParameters()
    ergParam.alpha = pDict['ergParam']['alpha']
    ergParam.beta = pDict['ergParam']['beta']
    ergParam.timeHorizon = pDict['ergParam']['timeHorizon']
    ergParam.time = np.array(pDict['ergParam']['time'])
    ergParam.wControl = pDict['ergParam']['wControl']

    eidParam = EIDParameters()
    eidParam.res = pDict['eidParam']['res']
    eidParam.eidType = pDict['eidParam']['eidType']
    eidParam.simType = pDict['eidParam']['simType']
    eidParam.pLHistDepth = pDict['eidParam']['pLHistDepth']
    eidParam.pLSigmaAmpBayesian = pDict['eidParam']['pLSigmaAmpBayesian']
    eidParam.pLSigmaAmpEID = pDict['eidParam']['pLSigmaAmpEID']
    eidParam.pLSigmaAmp = pDict['eidParam']['pLSigmaAmp']
    if type(pDict['eidParam']['rawTraj']) == bool or type(pDict['eidParam']['rawTraj']) == str:
        eidParam.rawTraj = pDict['eidParam']['rawTraj']
    else:
        eidParam.rawTraj = np.array(pDict['eidParam']['rawTraj'])
    eidParam.maxT = pDict['eidParam']['maxT']
    eidParam.objCenter = pDict['eidParam']['objCenter']
    eidParam.sInitPos = pDict['eidParam']['sInitPos']
    eidParam.blindIdx = pDict['eidParam']['blindIdx']
    eidParam.blind = pDict['eidParam']['blind']
    eidParam.saveDir = pDict['eidParam']['saveDir']
    if eidParam.simType == 'IF':
        eidParam.stepSize = pDict['eidParam']['stepSize']
        eidParam.planDepth = pDict['eidParam']['planDepth']
    # Conditions to simulate
    SNR = pDict['SNR']
    wControl = pDict['wControl']
    Sigma = pDict['Sigma']
    objAmp = pDict['objAmp']
    figName = pDict['figName']
    filename = pDict['filename']
    procNoiseSigma = pDict['procNoiseSigma']
    pLSigmaAmp = pDict['pLSigmaAmp']
    dt = pDict['dt']
    randSeed = pDict['randSeed']
    tRes = pDict['tRes']
    
    # Load parameters into global variable
    outDict = {}
    outDict['ergParam'] = ergParam
    outDict['eidParam'] = eidParam
    outDict['SNR'] = SNR
    outDict['wControl'] = wControl
    outDict['Sigma'] = Sigma
    outDict['objAmp'] = objAmp
    outDict['figName'] = figName
    outDict['filename'] = filename
    outDict['procNoiseSigma'] = procNoiseSigma
    outDict['pLSigmaAmp'] = pLSigmaAmp
    outDict['dt'] = dt
    outDict['randSeed'] = randSeed
    outDict['tRes'] = tRes
    return outDict
    
    
# Support encoder class to serialize numpy objects
class ErgodicParamEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)