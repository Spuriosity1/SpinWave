import numpy as np


# elementary rotation matrices

def Rx(theta):
    return np.array([
        [1,          0,            0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta),  np.cos(theta)]])

def Rz(theta):
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta),  np.cos(theta), 0],
        [          0,            0, 1]])

def Ry(theta):
    return np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [          0,            1,  0   ],
        [-np.sin(theta),  0, np.cos(theta)]])