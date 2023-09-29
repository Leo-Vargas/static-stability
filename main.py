import numpy as np
from scipy import linalg
import csv
from classes import Measurement, StateData
from functions import hSetup, is_symmetric, checkSensitivity, checkEquality, makeYbus
from JacobianCalculator import JacobianCalculator, ReducedJacobian

# --------------- entry data ---------------
measurementDict = {
    'p1': Measurement(0, 1, 0),
    'p2': Measurement(0, 2, 1),
    'p3': Measurement(0, 3, 1),
    'p4': Measurement(0, 4, 2),
    'p5': Measurement(0, 5, 2),
    'p6': Measurement(0, 6, 1),
    'p7': Measurement(0, 7, 2),
    'p8': Measurement(0, 8, 1),
    'p9': Measurement(0, 9, 2),
    'p10': Measurement(0, 10, 2),
    'p11': Measurement(0, 11, 2),
    'p12': Measurement(0, 12, 2),
    'p13': Measurement(0, 13, 2),
    'p14': Measurement(0, 14, 2),
    'q1': Measurement(1, 1, 0),
    'q2': Measurement(1, 2, 1),
    'q3': Measurement(1, 3, 1),
    'q4': Measurement(1, 4, 2),
    'q5': Measurement(1, 5, 2),
    'q6': Measurement(1, 6, 1),
    'q7': Measurement(1, 7, 2),
    'q8': Measurement(1, 8, 1),
    'q9': Measurement(1, 9, 2),
    'q10': Measurement(1, 10, 2),
    'q11': Measurement(1, 11, 2),
    'q12': Measurement(1, 12, 2),
    'q13': Measurement(1, 13, 2),
    'q14': Measurement(1, 14, 2),
}    
gridTopology = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])

linedata=np.array([
    (1,	2,	0.0193800000000000,	0.0591700000000000,	0.0264000000000000,	1),
    (1,	5,	0.0540300000000000,	0.223040000000000,	0.0246000000000000,	1),
    (2,	3,	0.0469900000000000,	0.197970000000000,	0.0219000000000000,	1),
    (2,	4,	0.0581100000000000,	0.176320000000000,	0.0170000000000000,	1),
    (2,	5,	0.0569500000000000,	0.173880000000000,	0.0173000000000000,	1),
    (3,	4,	0.0670100000000000,	0.171030000000000,	0.00640000000000000, 1),
    (4,	5,	0.0133500000000000,	0.0421100000000000,	0, 1),
    (4,	7,	0,	0.209120000000000,	0,	0.978000000000000),
    (4,	9,	0,	0.556180000000000,	0,	0.969000000000000),
    (5,	6,	0,	0.252020000000000,	0,	0.932000000000000),
    (6,	11,	0.0949800000000000,	0.198900000000000,	0,	1),
    (6,	12,	0.122910000000000,	0.255810000000000,	0, 1),
    (6,	13,	0.0661500000000000,	0.130270000000000,	0,	1),
    (7,	8,	0,	0.176150000000000,	0,	1),
    (7,	9,	0,	0.110010000000000,	0,	1),
    (9,	10,	0.0318100000000000,	0.0845000000000000,	0,	1),
    (9,	14,	0.127110000000000,	0.270380000000000,	0,	1),
    (10,	11,	0.0820500000000000,	0.192070000000000,	0,	1),
    (12,	13,	0.220920000000000,	0.199880000000000,	0,	1),
    (13,	14,	0.170930000000000,	0.348020000000000,	0,	1)
])



Ybus = makeYbus(linedata)
#print(Ybus)
#print('\n')
#is_symmetric(Ybus)



    
  
voltage = np.array([ 1.06, 1.045, 1.010, 1.018, 1.020, 1.070, 1.062, 1.090, 1.056, 1.051, 1.057, 1.055, 1.050, 1.036])
angle = np.array([-0.0, -5.0, -12.7, -10.3, -8.8, -14.2, -13.4, -13.4, -14.9, -15.1, -14.8, -15.1, -15.2, -16.])

[hValues, hDataDict] = hSetup(gridTopology, angle, voltage)

# --------------- Calculations ---------------
np.set_printoptions(precision=3)
jacobian = JacobianCalculator(list(measurementDict.values()), hDataDict, Ybus, hValues, gridTopology)
print(jacobian)
#print('\n')
 
#reducedJacobian = ReducedJacobian(jacobian)

#checkSensitivity(reducedJacobian)
#print(reducedJacobian)
#print('\n\n\n')
#print(reducedJacobian.shape)




#evalues, evectors = linalg.eig(reducedJacobian)

#evectorsInv = linalg.inv(evectors)
#evaluesMatrix = np.diag(evalues)
#print(evalues)
#print(evectors@D@evectorsInv)

