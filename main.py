import numpy as np
from scipy import linalg
import csv
from classes import Measurement, StateData
from functions import JacobianCalculator, hSetup, ReducedJacobian, is_symmetric, checkSensitivity, checkEquality, makeYbus

# --------------- entry data ---------------
measurementDict = {
    'p1': Measurement(0, (1,), 0),
    'p2': Measurement(0, (2,), 1),
    'p3': Measurement(0, (3,), 1),
    'p4': Measurement(0, (4,), 2),
    'p5': Measurement(0, (5,), 2),
    'p6': Measurement(0, (6,), 2),
    'p7': Measurement(0, (7,), 2),
    'p8': Measurement(0, (8,), 2),
    'p9': Measurement(0, (9,), 2),
    'q1': Measurement(1, (1,), 0),
    'q2': Measurement(1, (2,), 1),
    'q3': Measurement(1, (3,), 1),
    'q4': Measurement(1, (4,), 2),
    'q5': Measurement(1, (5,), 2),
    'q6': Measurement(1, (6,), 2),
    'q7': Measurement(1, (7,), 2),
    'q8': Measurement(1, (8,), 2),
    'q9': Measurement(1, (9,), 2)
}    
gridTopology = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])

linedata=np.array([
    (1, 4, 0, 0.0576, 0, 1),
    (4, 5, 0.017, 0.092, 0.158, 1),
    (5, 6, 0.039, 0.17, 0.358, 1),
    (3, 6, 0, 0.0586, 0, 1),
    (6, 7, 0.0119, 0.1008, 0.209, 1),
    (7, 8, 0.0085, 0.072, 0.149, 1),
    (8, 2, 0, 0.0625, 0, 1),
    (8, 9, 0.032, 0.161, 0.306, 1),
    (9, 4, 0.01, 0.085, 0.176, 1)
])

Ybus = makeYbus(linedata)
print(Ybus)
print('\n')
is_symmetric(Ybus)



    
"""    
voltage = np.array([ 1.043, 1.022, 1.013, 1.01, 1.012, 1.003, 1.01, 1.051, 1.044, 1.082, 1.057, 1.071, 1.042, 1.038, 1.045, 1.039, 1.028, 1.025, 1.029, 1.032, 1.033, 1.027, 1.022, 1.019, 1.001, 1.026, 1.011, 1.006, 0.995])
angle = np.array([-5.497, -8.004, -9.661, -14.381, -11.398, -13.158, -12.115, -14.434, -16.024, -14.434, -15.302, -15.302, -16.191, -16.278, -15.88, -16.188, -16.884, -17.052, -16.852, -16.468, -16.455, -16.662, -16.83, -16.424, -16.842, -15.912, -12.057, -17.136, -18.015])

[hValues, hDataDict] = hSetup(gridTopology, angle, voltage)

# --------------- Calculations ---------------
np.set_printoptions(precision=3)
jacobian = JacobianCalculator(list(measurementDict.values()), hDataDict, Ybus, hValues, gridTopology)
#print(jacobian)
#print('\n')
 
reducedJacobian = ReducedJacobian(jacobian)

#checkSensitivity(reducedJacobian)
#print(reducedJacobian)
#print('\n\n\n')
#print(reducedJacobian.shape)




evalues, evectors = linalg.eig(reducedJacobian)

evectorsInv = linalg.inv(evectors)
evaluesMatrix = np.diag(evalues)
#print(evalues)
#print(evectors@D@evectorsInv)

P, D, Q = np.linalg.svd(reducedJacobian, full_matrices=False)
print('----------------barreira ------------------')
print(D)

"""