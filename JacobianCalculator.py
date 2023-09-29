import numpy as np
from scipy import linalg
import math
from classes import Measurement, StateData, busTypes

# ----------------------- REDUCED JACOBIAN --------------------------
def ReducedJacobian(jacobian: np.ndarray):

    jacobianPA = jacobian[0:29, 0:29]
    jacobianPV = jacobian[0:29, 29:59]
    jacobianQA = jacobian[29:59, 0:29]
    jacobianQV = jacobian[29:59, 29:59]
    
    reducedJacobian = jacobianQV-jacobianQA@linalg.inv(jacobianPA)@jacobianPV

    return reducedJacobian



# ------------------------ JacobianCalculator internal functions -----------------------------------
def PinjCalculatorAngleSame(hValues, hDataDict, Ybus, calcBus, gridTopology):   
    sum = 0.0
    Vi = hValues[hDataDict['v'+str(calcBus)].index]
    weight = (Vi*Vi*Ybus[calcBus-1, calcBus-1].imag)
    
    angles = np.zeros(2)
    angles[0] = hValues[hDataDict['a'+str(calcBus)].index]
    
    for bus in gridTopology:
        angles[1] = hValues[hDataDict['a'+str(bus)].index]
        Voltages = Vi*hValues[hDataDict['v'+str(bus)].index]
                
        rightSide = -Ybus[calcBus-1, bus-1].real*math.sin(angles[0]-angles[1]) + Ybus[calcBus-1, bus-1].imag*math.cos(angles[0]-angles[1])
            
        sum += Voltages*rightSide            
    
    return sum - weight

def PinjCalculatorAngleDiff(hValues, hDataDict, Ybus, calcBus, hData):
    angles = np.zeros(2)
    
    angles[0] = hValues[hDataDict['a'+str(calcBus)].index]
    angles[1] = hValues[hData.index]
    Voltages = hValues[hDataDict['v'+str(calcBus)].index]*hValues[hDataDict['v'+str(hData.bus)].index]
    
    rightSide = Ybus[calcBus-1, hData.bus-1].real*math.sin(angles[0]-angles[1]) - Ybus[calcBus-1, hData.bus-1].imag*math.cos(angles[0]-angles[1])  
  
    return Voltages*rightSide

def PinjCalculatorVoltageSame(hValues, hDataDict, Ybus, calcBus, gridTopology):
    sum = 0.0
    Vi = hValues[hDataDict['v'+str(calcBus)].index]
    weight = (Vi*Ybus[calcBus-1, calcBus-1].real)
    
    angles = np.zeros(2)
    angles[0] = hValues[hDataDict['a'+str(calcBus)].index]
    
    for bus in gridTopology:
        if bus == 1:
            angles[1] = 0.0
            Voltages = 1.06
        else:
            angles[1] = hValues[hDataDict['a'+str(bus)].index]
            Voltages = hValues[hDataDict['v'+str(bus)].index]
                
        rightSide = Ybus[calcBus-1, bus-1].real*math.cos(angles[0]-angles[1]) + Ybus[calcBus-1, bus-1].imag*math.sin(angles[0]-angles[1])
            
        sum += Voltages*rightSide            
    
    return sum + weight

def PinjCalculatorVoltageDiff(hValues, hDataDict, Ybus, calcBus, hData):
    angles = np.zeros(2)
    
    angles[0] = hValues[hDataDict['a'+str(calcBus)].index]
    angles[1] = hValues[hDataDict['a'+str(hData.bus)].index]
    Voltages = hValues[hDataDict['v'+str(calcBus)].index]
    
    rightSide = Ybus[calcBus-1, hData.bus-1].real*math.cos(angles[0]-angles[1]) + Ybus[calcBus-1, hData.bus-1].imag*math.sin(angles[0]-angles[1])  
    
    return Voltages*rightSide



def QinjCalculatorAngleSame(hValues, hDataDict, Ybus, calcBus, gridTopology):
    sum = 0.0
    Vi = hValues[hDataDict['v'+str(calcBus)].index]
    weight = (Vi*Vi*Ybus[calcBus-1, calcBus-1].real)
    
    angles = np.zeros(2)
    angles[0] = hValues[hDataDict['a'+str(calcBus)].index]
    
    for bus in gridTopology:
        if bus == 1:
            angles[1] = 0.0
            Voltages = 1.06
        else:
            angles[1] = hValues[hDataDict['a'+str(bus)].index]
            Voltages = Vi*hValues[hDataDict['v'+str(bus)].index]
                
        rightSide = Ybus[calcBus-1, bus-1].real*math.cos(angles[0]-angles[1]) + Ybus[calcBus-1, bus-1].imag*math.sin(angles[0]-angles[1])
            
        sum += Voltages*rightSide            
    
    return sum - weight

def QinjCalculatorAngleDiff(hValues, hDataDict, Ybus, calcBus, hData):
    angles = np.zeros(2)
    
    angles[0] = hValues[hDataDict['a'+str(calcBus)].index]
    angles[1] = hValues[hData.index]
    Voltages = hValues[hDataDict['v'+str(calcBus)].index]*hValues[hDataDict['v'+str(hData.bus)].index]
    
    rightSide = -Ybus[calcBus-1, hData.bus-1].real*math.cos(angles[0]-angles[1]) - Ybus[calcBus-1, hData.bus-1].imag*math.sin(angles[0]-angles[1])  
  
    return Voltages*rightSide

def QinjCalculatorVoltageSame(hValues, hDataDict, Ybus, calcBus, gridTopology):
    sum = 0.0
    Vi = hValues[hDataDict['v'+str(calcBus)].index]
    weight = (Vi*Ybus[calcBus-1, calcBus-1].imag)
    
    angles = np.zeros(2)
    angles[0] = hValues[hDataDict['a'+str(calcBus)].index]
    
    for bus in gridTopology:
        if bus == 1:
            angles[1] = 0.0
            Voltages = 1.06
        else:
            angles[1] = hValues[hDataDict['a'+str(bus)].index]
            Voltages = hValues[hDataDict['v'+str(bus)].index]
                
        rightSide = Ybus[calcBus-1, bus-1].real*math.sin(angles[0]-angles[1]) - Ybus[calcBus-1, bus-1].imag*math.cos(angles[0]-angles[1])
            
        sum += Voltages*rightSide            
    
    return sum - weight

def QinjCalculatorVoltageDiff(hValues, hDataDict, Ybus, calcBus, hData):
    angles = np.zeros(2)
    
    angles[0] = hValues[hDataDict['a'+str(calcBus)].index]
    angles[1] = hValues[hDataDict['a'+str(hData.bus)].index]
    Voltages = hValues[hDataDict['v'+str(calcBus)].index]
    
    rightSide = Ybus[calcBus-1, hData.bus-1].real*math.sin(angles[0]-angles[1]) - Ybus[calcBus-1, hData.bus-1].imag*math.cos(angles[0]-angles[1])  
    
    return Voltages*rightSide


# ------------------------ JacobianCalculator related functions ------------------------------------

def Pinj(measurement, hDataDict, Ybus, hValues, gridTopology, busConfiguration):
    jacobianRow = np.zeros(busConfiguration['PQ'].count*2 + busConfiguration['PV'].count)
    i = 0
    
    
    for hData in hDataDict.values():
        if hData.bus != busConfiguration['slack'].buses[0]:
            if hData.category == 0:
                if measurement.bus == hData.bus:
                    jacobianRow[i]=PinjCalculatorAngleSame(hValues, hDataDict, Ybus, measurement.bus, gridTopology)
                    i+=1
                
                else:
                    jacobianRow[i]=PinjCalculatorAngleDiff(hValues, hDataDict, Ybus, measurement.bus, hData)
                    i+=1
                
                
                    

            elif hData.category == 1:
                if all(PVbuses != hData.bus for PVbuses in busConfiguration['PV'].buses):
                    if measurement.bus == hData.bus:
                        jacobianRow[i]=PinjCalculatorVoltageSame(hValues, hDataDict, Ybus, measurement.bus, gridTopology)
                        i+=1
                        
                    else:
                        jacobianRow[i]=PinjCalculatorVoltageDiff(hValues, hDataDict, Ybus, measurement.bus, hData)
                        i+=1   
            
            
    
    return jacobianRow


def Qinj(measurement, hDataDict, Ybus, hValues, gridTopology, busConfiguration):
    jacobianRow = np.zeros(busConfiguration['PQ'].count*2 + busConfiguration['PV'].count)
    
    
    for i, (key, hData) in enumerate(hDataDict.items()):
        if hData.bus != busConfiguration['slack'].buses[0]:
            if hData.category == 0:
                if measurement.bus == hData.bus:
                    jacobianRow[i]=QinjCalculatorAngleSame(hValues, hDataDict, Ybus, measurement.bus, gridTopology)
                
                else:
                    jacobianRow[i]=QinjCalculatorAngleDiff(hValues, hDataDict, Ybus, measurement.bus, hData)
                    

            elif hData.category == 1:
                if all(PVbuses != hData.bus for PVbuses in busConfiguration['PV'].buses):
                    if measurement.bus == hData.bus:
                        jacobianRow[i]=QinjCalculatorVoltageSame(hValues, hDataDict, Ybus, measurement.bus, gridTopology)
                        
                    else:
                        jacobianRow[i]=QinjCalculatorVoltageDiff(hValues, hDataDict, Ybus, measurement.bus, hData)    
    
    return jacobianRow

   
# -------------------- JacobianCalculator Function -----------------------------------------------------
def JacobianCalculator(measurementList: list, hDataDict: dict, Ybus: np.ndarray, hValues: np.ndarray, gridTopology):
    
    busConfiguration = {
        'slack': busTypes(),
        'PV': busTypes(),
        'PQ': busTypes() 
    }
    
    for j in range(len(gridTopology)):
        if measurementList[j].busType == 0:
            busConfiguration['slack'].addBus(gridTopology[j])
        elif measurementList[j].busType == 1:
            busConfiguration['PV'].addBus(gridTopology[j])
        elif measurementList[j].busType == 2:
            busConfiguration['PQ'].addBus(gridTopology[j])

    jacobian = np.zeros((busConfiguration['PQ'].count*2 + busConfiguration['PV'].count, busConfiguration['PQ'].count*2 + busConfiguration['PV'].count), dtype=float)
    
    print(np.shape(jacobian))
    
    calcSelector = {
        0: Pinj,
        1: Qinj,
    }
    
    for j in range(len(measurementList)):
        if (measurementList[j].busType) == 2 or (measurementList[j].busType == 1 and measurementList[j].category == 0):
            jacobian[j]=calcSelector[measurementList[j].category](measurementList[j], hDataDict, Ybus, hValues, gridTopology, busConfiguration)
            
    
    return jacobian