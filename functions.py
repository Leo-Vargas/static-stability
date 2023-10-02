import numpy as np
from scipy import linalg
import math
from classes import Measurement, StateData, busTypes


# ----------------------- Data Processing Related Functions -----------------------------
def hSetup(gridTopology: dict, angle, voltage):
    hValues = np.ones(((len(gridTopology)*2)), dtype=float)
    hDataDict = {}
    
    for i in range((len(gridTopology)*2)):
        if i < len(gridTopology):
            hValues[i]=angle[i]*math.pi/180
            hDataDict.update({'a'+str(gridTopology[i]): StateData(0, i, gridTopology[i])})
        else:
            hValues[i]=voltage[i-(len(gridTopology))]
            hDataDict.update({'v'+str(gridTopology[i-(len(gridTopology))]): StateData(1, i, gridTopology[i-(len(gridTopology))])})

    hDataDictUpdate(hValues, hDataDict)

    return [hValues, hDataDict]

def hDataDictUpdate(hvalues: np.ndarray, hDataDict: dict):
   
    for i, (key, hData) in enumerate(hDataDict.items()):
        hData.addValue(hvalues[i])
        


def is_symmetric(a, rtol=1e-05, atol=1e-08):
    if np.allclose(a, a.T, rtol=rtol, atol=atol) == True:
        print('ok')
    else:
        print('not ok')

def checkEquality(a, b, rtol=1e-05, atol=1e-08):
    if np.allclose(a, b, rtol=rtol, atol=atol) == True:
        print('ok')
    else:
        print('not ok')

def checkSensitivity(reducedJacobian: np.ndarray):
    sensitivity = np.diag(linalg.inv(reducedJacobian))
    
    for i in range(len(sensitivity)):
        print(sensitivity[i])

# ---------------------- YBUS --------------------
def makeYbus(linedata: np.ndarray):
    nl = np.array(linedata[:, 0], dtype=int)
    nr = np.array(linedata[:, 1], dtype=int)
    R = np.array(linedata[:, 2], dtype=complex)
    X = np.array(linedata[:, 3], dtype=complex)
    Bc = np.array(1j*linedata[:, 4], dtype=complex)
    a = np.array(linedata[:, 5])

    nbr = len(nl)
    nbus = np.max([np.max(nl), np.max(nr)])
    Z = R + X*1j

    y = np.divide(np.ones(nbr), Z)

    Ybus=np.zeros((int(nbus), int(nbus)), dtype=complex)

    for i in range(nbr):
        if a[i] <= 0:
            a[i] = 1
        
        nl[i]-=1
        nr[i]-=1

        Ybus[nl[i], nr[i]] = Ybus[nl[i], nr[i]] - y[i]/a[i]
        Ybus[nr[i], nl[i]] = Ybus[nl[i], nr[i]]

    for n in range(nbus):
        for k in range(nbr):
            if nl[k] == n:
                Ybus[n, n] = Ybus[n, n]+y[k]/(a[k]**2) + Bc[k]
            elif nr[k] == n:
                Ybus[n, n] = Ybus[n, n] + y[k] + Bc[k]

    return Ybus

def makeBusConfiguration(measurementDict: dict):
    busConfiguration = {
        'slack': busTypes(),
        'PV': busTypes(),
        'PQ': busTypes() 
    }
    
    for measurement in measurementDict.values():
        if measurement.category == 0:
            if measurement.busType == 0:
                busConfiguration['slack'].addBus(measurement.bus)
            elif measurement.busType == 1:
                busConfiguration['PV'].addBus(measurement.bus)
            elif measurement.busType == 2:
                busConfiguration['PQ'].addBus(measurement.bus)
        else:
            return busConfiguration
        


# ------------------------ PARTICIPATION ---------------
def calcBusParticipation(rightEvectors, leftEvectors):
    busParticipation = np.zeros((np.shape(rightEvectors)[0], np.shape(rightEvectors)[1]))

    for i in range(np.shape(rightEvectors)[0]):
        for j in range(np.shape(rightEvectors)[1]):
            busParticipation[i,j] = rightEvectors[i,j]*leftEvectors[j,i]
    
    return busParticipation


def printBusParticipation(busParticipation):

    for i in range(np.shape(busParticipation)[0]):
        print('')
        print(f'Modo {i+1}:', end=' ')
        for j in range(np.shape(busParticipation)[1]):
            print(f'{busParticipation[j,i]}', end=' ')