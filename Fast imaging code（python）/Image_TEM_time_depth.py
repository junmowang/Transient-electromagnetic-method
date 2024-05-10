# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 19:53:26 2024

@author: zhu_y
"""

#### Time depth conversion program
## It includes: (1) the time-depth transformation function of numerical simulation; (2) Time-depth transformation function of measured data

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

############################# (1) the time-depth transformation function of numerical simulation #############################
def T_H_1(rou0,TT,H,S,Te_1,Te_2):
    mu0=constants.mu_0#Vacuum permeability
    epsilon0=constants.epsilon_0#permittivity of vacuum
    dx=np.linalg.norm(Te_1-Te_2)#The length of the dipole
    #Mr. Chen. - Time-depth conversion
    def chen_t_h(rou,TT):
        rou=abs(rou)
        h=1.0*(2*TT*rou/mu0)**0.5#1.25/1.0/1.2
        return(h)

    #eaton- time-depth conversion
    def eaton(rou,TT,L):
        rou=abs(rou)
        P=(100*np.pi*TT*rou/(mu0*L**2))**0.25
        h=1000*(rou*TT*np.pi)**0.5*(1-np.exp(-P)/2)
        return(h)

    def h_model(H,S,rou0,h):#Standard curve of the model
        rou0=abs(rou0)
        H_h=np.arange(int(h.min()),math.ceil(h.max())+1)
        H_sum=np.cumsum(H)
        S_rou0=np.zeros(H_h.shape)
        n0=np.arange(H.shape[0])
        for i in n0:
            if i==n0[0]:
                S_rou0[0:np.argwhere(H_h==H[0])[0,0]]=1/S[i]
            elif i==n0[-1]:
                S_rou0[np.argwhere(H_h==sum(H[:-1]))[0,0]:np.argwhere(H_h==H_h[-1])[0,0]+1]=1/S[i]
            else:
                S_rou0[np.argwhere(H_h==H_sum[i-1])[0,0]:np.argwhere(H_h==H_sum[i])[0,0]]=1/S[i]
        return(H_h,S_rou0)

    chen_h1=chen_t_h(rou0,TT)#Mr. Chen. - Time Depth formula
    eaton_h1=eaton(rou0,TT,dx)#eaton- time depth formula

    H_hc,S_rouc=h_model(H,S,rou0,chen_h1)#Model standard curve -- Chen
    H_he,S_roue=h_model(H,S,rou0,eaton_h1)#Model standard curve -- eaton
    
    return chen_h1,H_hc,S_rouc,eaton_h1,H_he,S_roue#Chen depth, model depth, model resistivity; eaton depth, model depth, model resistivit


############################# (2) Time-depth transformation function of measured data #############################
def T_H_2(rou0,T,Te_1,Te_2):
    mu0=constants.mu_0#permeability of vacuum
    epsilon0=constants.epsilon_0#permittivity of vacuum
    dx=np.linalg.norm(Te_1-Te_2)#The length of the dipole
    TT=T.reshape(-1,1).repeat(rou0.shape[1],axis=1)
    #Mr. Chen. - Time-depth conversion
    def chen_t_h(rou,TT):
        rou=abs(rou)
        h=1.0*(2*TT*rou/mu0)**0.5#1.25/1.0/1.2
        return(h)

    #eaton- time-depth conversion
    def eaton(rou,TT,L):
        rou=abs(rou)
        P=(100*np.pi*TT*rou/(mu0*L**2))**0.25
        h=1000*(rou*TT*np.pi)**0.5*(1-np.exp(-P)/2)
        return(h)

    chen_h1=chen_t_h(rou0,TT)#Mr. Chen. - Time Depth formula
    eaton_h1=eaton(rou0,TT,dx)#eaton- time depth formula
    return chen_h1,eaton_h1#Chen Shin; eaton depth


