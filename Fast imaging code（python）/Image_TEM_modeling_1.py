# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:00:45 2024

@author: zhu_y
"""

#### An example model of transient electromagnetic forward modeling
### It includes: 
### (1) independent variables and model parameters; 
### (2) one-dimensional, uniform forward modeling; 
### (3) Calculate the resistivity of the whole period; 
### (4) resistivity correction; 
### (5) Time-depth transformation



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

##################### (1) independent variables and model parameters #####################
mu0=constants.mu_0#permeability of vacuum
epsilon0=constants.epsilon_0#permittivity of vacuum

##One-dimensional layered medium model parameters
H=np.array([500,300,float('inf')])#Three-layer thickness
S=np.array([0.01,0.1,0.01])#Electrical conductivity of each layer
epsilon_r=np.array([1.0,1.0,1.0])#The relative dielectric constant of each layer is 1.0

I=1#Power supply strength

Te_1=np.array([[-250,0,0]])#Long wire left end point
Te_2=np.array([[250,0,0]])#Long wire right end point


Re=np.array([[0,500,0]])#Different offset
# Re=np.array([[0,350,0],[0,500,0]])

T=np.logspace(-5,0,50)#observation time

plt.rcParams['font.sans-serif']=['Microsoft YaHei'] #It is used to display Chinese labels properly
start=time.time()

##################### (2) one-dimensional, uniform forward modeling #####################
judge=2### 1 —— dipole uniform half-space forward modeling; 2 —— one-dimensional forward modeling of long wires

import TEM_E_B_forwardmodeling_1
ex,vz=TEM_E_B_forwardmodeling_1.forwardmodeling(I,Re,Te_1,Te_2,T,H,S,epsilon_r,judge)

plt.figure(dpi=250)
plt.loglog(T,abs(ex))
plt.title('Electric field Ex')
plt.legend(['[0,350]','[0,500]','[0,1000]','[0,2000]'])
plt.figure(dpi=250)
plt.loglog(T,abs(vz))
plt.title('Magnetic field Vz')
plt.legend(['[0,350]','[0,500]','[0,1000]','[0,2000]'])

end1=time.time()

##################### (3) Calculate the resistivity of the whole period #####################
import TEM_Ex_Vz_rou
##Based on ex peak time imaging
rou_log_ex,rou_nor_ex,TT=TEM_Ex_Vz_rou.rou_ex(ex,I,Re,Te_1,Te_2,T)

plt.figure(dpi=250)
plt.loglog(TT,rou_nor_ex)
plt.loglog(TT,rou_log_ex,'-.')
plt.title('Electric field Ex - total apparent resistivity')
plt.legend(['[0,350]','[0,500]','[0,1000]','[0,2000]'])
# plt.xlim(0,1000)
# plt.ylim(10,1000)

##Based on vz peak time imaging
rou_log_vz,rou_nor_vz,TT=TEM_Ex_Vz_rou.rou_vz(vz,I,Re,Te_1,Te_2,T)

plt.figure(dpi=250)
plt.loglog(TT,rou_nor_vz)
plt.loglog(TT,rou_log_vz,'-.')
plt.title('Magnetic field Vz - total apparent resistivity')
plt.legend(['[0,350]','[0,500]','[0,1000]','[0,2000]'])

end2=time.time()

##################### (4) resistivity correction #####################

S0=0.01##Arbitrary background initial resistivity

### Calculate accurate background resistivity
import TEM_background_resistivity
S0_ex,S0_vz=TEM_background_resistivity.background_resistivity(S0,I,Re,Te_1,Te_2,T,epsilon_r,rou_log_ex,rou_log_vz)## ex 背景场值； vz 背景场值

S0=S0_ex[0]#Take the first resistivity as the final background resistivity

### Make corrections
import TEM_correction

### The correction coefficient uses one-dimensional uniformity
S2=S/S*S0##Background resistivity
Kg_nor_ex1,Kg_log_ex1,Kg_nor_vz1,Kg_log_vz1=TEM_correction.rou_correction_2(I,Re,Te_1,Te_2,T,H,S2,epsilon_r)

rou_log_ex4=Kg_log_ex1*rou_log_ex##Corrected ex resistivity (logarithmic domain)
rou_nor_ex4=Kg_nor_ex1*rou_nor_ex#(linear domain)
rou_log_vz4=Kg_log_vz1*rou_log_vz##Corrected vz resistivity (logarithmic domain)
rou_nor_vz4=Kg_nor_vz1*rou_nor_vz#(linear domain)

plt.figure(dpi=250)
plt.loglog(TT,rou_nor_ex4)
plt.loglog(TT,rou_log_ex4,'-.')
plt.title('Electric field Ex -- lifetime apparent resistivity (after correction)')
plt.legend(['[0,350]','[0,500]','[0,1000]','[0,2000]'])
plt.figure(dpi=250)
plt.loglog(TT,rou_nor_vz4)
plt.loglog(TT,rou_log_vz4,'-.')
plt.title('Magnetic field Vz - Lifetime apparent resistivity (corrected)')
plt.legend(['[0,350]','[0,500]','[0,1000]','[0,2000]'])

end3=time.time()

##################### (5) Time-depth transformation #####################

import TEM_time_depth

# depth(chen), model depth, model resistivity; eaton depth, model depth, model resistivity
chen_h1_ex,H_hc_ex,S_rouc_ex,eaton_h1_ex,H_he_ex,S_roue_ex=TEM_time_depth.T_H_1(rou_log_ex4,TT,H,S,Te_1,Te_2)#反演（2）校正后数据
chen_h1_vz,H_hc_vz,S_rouc_vz,eaton_h1_vz,H_he_vz,S_roue_vz=TEM_time_depth.T_H_1(rou_log_vz4,TT,H,S,Te_1,Te_2)#反演（2）


plt.figure(dpi=250)
plt.loglog(chen_h1_ex,rou_log_ex4)
plt.title('Total apparent resistivity (ex)- Depth - chen')
plt.loglog(H_hc_ex,S_rouc_ex,color='black')

plt.figure(dpi=250)
plt.loglog(eaton_h1_ex,rou_log_ex4)
plt.title('Total apparent resistivity (ex)- Depth -eaton')
plt.loglog(H_he_ex,S_roue_ex,color='black')

plt.figure(dpi=250)
plt.loglog(chen_h1_vz,rou_log_vz4)
plt.title('Total apparent resistivity (vz)- Depth - chen')
plt.loglog(H_hc_vz,S_rouc_vz,color='black')

plt.figure(dpi=250)
plt.loglog(eaton_h1_vz,rou_log_vz4)
plt.title('Total apparent resistivity (vz)- Depth -eaton')
plt.loglog(H_he_vz,S_roue_vz,color='black')


print('Forward time=',end1-start)
print('Imaging time=',end2-end1)
print('Correction time=',end3-end2)
print('Total time spent=',end3-start)