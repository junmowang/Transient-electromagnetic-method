# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 20:28:14 2024

@author: zhu_y
"""

#### Program for calculating background resistivity
## Including: (1) numerical simulation of background resistivity; (2) Background resistivity is calculated from measured data

from scipy import interpolate
import numpy as np
import pandas as pd
#import sympy as sy
# import scipy as sc
from scipy import constants
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import math
#import cmath
import time


################################ （1）The background resistivity was calculated by numerical simulation ################################
def background_resistivity(S0,I,Re,Te_1,Te_2,T,epsilon_r,rou_log_ex,rou_log_vz):
    import TEM_correction
    H1=np.array([float('inf')])# Formation thickness
    S1=np.array([S0])## Background resistivity
    Kg_nor_ex1,Kg_log_ex1,Kg_nor_vz1,Kg_log_vz1=TEM_correction.rou_correction_1(I,Re,Te_1,Te_2,T,H1,S1,epsilon_r)
    rou_log_ex4=Kg_log_ex1*rou_log_ex##Corrected ex resistivity (logarithmic domain)
    rou_log_vz4=Kg_log_vz1*rou_log_vz##The resistivity of vz after correction
    S0_ex=rou_log_ex4[0,:]**-1
    S0_vz=rou_log_vz4[0,:]**-1
    return S0_ex,S0_vz## ex background conductivity value; vz background conductivity value

################################ （2）Background resistivity is calculated from measured data ################################
def background_resistivity_1(I,Re,Te_1,Te_2,T,epsilon_r,rou_log_ex):##Re_ex,T_ex;Re_vz,T_vz
    S0=abs(rou_log_ex[0,:])**-1#Take the first row as the initial conductivity
    H1=np.array([float('inf')])#formation thickness
    S1=np.array([S0])## Background conductivity (Note: S1 in this case is a one-dimensional array, not a scalar)
    import TEM_correction
    
    Kg_nor_ex1=np.zeros([rou_log_ex.shape[0],rou_log_ex.shape[1]])#Matrix to be assigned
    Kg_nor_vz1=np.zeros([rou_log_ex.shape[0],rou_log_ex.shape[1]])
    Kg_log_ex1=np.zeros([rou_log_ex.shape[0],rou_log_ex.shape[1]])
    Kg_log_vz1=np.zeros([rou_log_ex.shape[0],rou_log_ex.shape[1]])
    
    for j in np.arange(Re.shape[0]):##Calculate the resistivity background values at different measuring points
        Re1=Re[j,:].reshape(1,-1)
        S11=np.array([S1[0,j]])
        Kg_nor_ex1c,Kg_log_ex1c,Kg_nor_vz1c,Kg_log_vz1c=TEM_correction.rou_correction_1(I,Re1,Te_1,Te_2,T,H1,S11,epsilon_r)
        Kg_nor_ex1[:,j]=Kg_nor_ex1c.reshape(-1)
        Kg_log_ex1[:,j]=Kg_log_ex1c.reshape(-1)
        Kg_nor_vz1[:,j]=Kg_nor_vz1c.reshape(-1)
        Kg_log_vz1[:,j]=Kg_log_vz1c.reshape(-1)
    rou_log_ex4=Kg_log_ex1*rou_log_ex##Corrected ex resistivity (logarithmic domain)
    rou_log_vz4=Kg_log_vz1*rou_log_ex##The resistivity of vz after correction
    S0_ex=abs(rou_log_ex4[0,:])**-1##More accurate conductivity background value
    S0_vz=abs(rou_log_vz4[0,:])**-1
    return S0_ex,S0_vz## ex background conductivity value; vz background conductivity value