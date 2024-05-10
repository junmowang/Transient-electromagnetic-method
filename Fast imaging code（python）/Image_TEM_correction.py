# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 09:43:10 2024

@author: zhu_y
"""

### Calculation of correction factor program
## Eliminate the geometric effects of fast imaging resistivity
## It includes: (1) long wire uniform half-space time domain response program; (2) Correction coefficient program based on uniform half space; 
## (3) Calibration coefficient program based on one-dimensional layered media; (4) Calibration coefficient program of measured data

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


import TEM_E_B_forwardmodeling_1
import TEM_Ex_Vz_rou

############# (1) long wire uniform half-space time domain response program #############
def L_halfspace_T(I,Re,Te_1,Te_2,T,H,S,NL):
    mu0=constants.mu_0#permeability of vacuum
    epsilon0=constants.epsilon_0#permittivity of vacuum
    
    ##Long conductor parameter
    Te_0=(Te_1+Te_2)/2#Midpoint coordinates of long wire
    R=np.zeros(len(Re))#Wire midpoint, transmit-to-receive distance to be assigned matrix (vector)
    for i in np.arange(len(Re)):#Field point and conductor center point, transmit-receive distance to be assigned matrix (vector)
        R[i]=np.linalg.norm((Re-Te_0)[i])#Transmit-receive distance
    L=np.linalg.norm(Te_1-Te_2)#wire length
    ds=L/NL#The length of the dipole
    L_sin=(Te_2-Te_1)[0,1]/L#Y-axis change rate of a long wire
    L_cos=(Te_2-Te_1)[0,0]/L#Rate of change on the X-axis of a long wire
    Te_ou=np.zeros([NL,Te_1.shape[1]])
    for i in np.arange(NL):#Calculate the midpoint coordinates of each dipole
        Te_ou[i,0]=Te_1[0,0]+ds/2*L_cos+ds*i*L_cos
        Te_ou[i,1]=Te_1[0,1]+ds/2*L_sin+ds*i*L_sin
    R0=np.zeros([NL,Re.shape[0]])#Line source, transmit-receive distance to be assigned matrix
    for j in np.arange(Re.shape[0]):#Calculate the transmit-receive distance matrix
        for l in np.arange(NL):
            R0[l,j]=np.linalg.norm(Re[j,:]-Te_ou[l,:])#Transmit-receive distance matrix
    
    ##### vz, ex long wire - uniform half space - time domain
    vec_erf=np.frompyfunc(math.erf,1,1)#Error function erf vectorization in math
    ex0=np.zeros([len(T),Re.shape[0]])
    vz0=np.zeros([len(T),Re.shape[0]])
    for j in np.arange(Re.shape[0]):
        for k in np.arange(len(T)):
            K0_ex=0
            K0_vz=0
            for l in np.arange(NL):
                u=0.5*(mu0*S[0]/T[k])**0.5*R0[l,j]
                erf1=vec_erf(u)#.astype('float64')
                K1_ex=1/R0[l,j]**3*(erf1-2/np.pi**0.5*u*np.exp(-u**2))
                K0_ex=K0_ex+K1_ex
                K1_vz=(Re[j,1]-Te_ou[l,1])/R0[l,j]**5*(3*erf1-2/np.pi**0.5*u*(3+2*u**2)*np.exp(-u**2))
                K0_vz=K0_vz+K1_vz
            ex0[k,j]=K0_ex
            vz0[k,j]=K0_vz
    ex=-I/(2*np.pi*S[0])*ex0*ds
    vz=I/(2*np.pi*mu0*S[0])*vz0*ds
    
    return ex,vz


############# (2) Correction coefficient program based on uniform half space #############
def rou_correction_1(I,Re,Te_1,Te_2,T,H1,S1,epsilon_r):
    ex_ou,vz_ou=TEM_E_B_forwardmodeling_1.forwardmodeling(I,Re,Te_1,Te_2,T,H1,S1,epsilon_r,1)#Calculate the uniform half-space-dipole response
    ex_oL,vz_oL=L_halfspace_T(I,Re,Te_1,Te_2,T,H1,S1,10)#Calculate uniform half space - long conductor response
    
    rou_log_ex2,rou_nor_ex2,TT=TEM_Ex_Vz_rou.rou_ex(ex_ou,I,Re,Te_1,Te_2,T)#Calculate uniform half-space - dipole response-resistivity
    rou_log_vz2,rou_nor_vz2,TT=TEM_Ex_Vz_rou.rou_vz(vz_ou,I,Re,Te_1,Te_2,T)
    #
    rou_log_ex3,rou_nor_ex3,TT=TEM_Ex_Vz_rou.rou_ex(ex_oL,I,Re,Te_1,Te_2,T)#Calculate uniform half-space - long conductor response - resistivity
    rou_log_vz3,rou_nor_vz3,TT=TEM_Ex_Vz_rou.rou_vz(vz_oL,I,Re,Te_1,Te_2,T)
    #
    Kg_log_ex1=rou_log_ex2/rou_log_ex3#ex resistivity response - correction factor
    Kg_nor_ex1=rou_nor_ex2/rou_nor_ex3
    Kg_log_vz1=rou_log_vz2/rou_log_vz3#vz resistivity response - correction factor
    Kg_nor_vz1=rou_nor_vz2/rou_nor_vz3
    
    return Kg_nor_ex1,Kg_log_ex1,Kg_nor_vz1,Kg_log_vz1

############# (3) Calibration coefficient program based on one-dimensional layered media #############
def rou_correction_2(I,Re,Te_1,Te_2,T,H,S2,epsilon_r):
    ex_ou,vz_ou=TEM_E_B_forwardmodeling_1.forwardmodeling(I,Re,Te_1,Te_2,T,H,S2,epsilon_r,1)#Calculate the uniform half-space-dipole response
    ex_L,vz_L=TEM_E_B_forwardmodeling_1.forwardmodeling(I,Re,Te_1,Te_2,T,H,S2,epsilon_r,2)#Calculate uniform half-space-long conductor response
    
    rou_log_ex2,rou_nor_ex2,TT=TEM_Ex_Vz_rou.rou_ex(ex_ou,I,Re,Te_1,Te_2,T)#Calculate uniform half-space - dipole response-resistivity
    rou_log_vz2,rou_nor_vz2,TT=TEM_Ex_Vz_rou.rou_vz(vz_ou,I,Re,Te_1,Te_2,T)
    #
    rou_log_ex3,rou_nor_ex3,TT=TEM_Ex_Vz_rou.rou_ex(ex_L,I,Re,Te_1,Te_2,T)#Calculate uniform half-space - long conductor response - resistivity
    rou_log_vz3,rou_nor_vz3,TT=TEM_Ex_Vz_rou.rou_vz(vz_L,I,Re,Te_1,Te_2,T)
    #
    Kg_log_ex1=rou_log_ex2/rou_log_ex3#ex resistivity response - correction factor
    Kg_nor_ex1=rou_nor_ex2/rou_nor_ex3
    Kg_log_vz1=rou_log_vz2/rou_log_vz3#vz resistivity response - correction factor
    Kg_nor_vz1=rou_nor_vz2/rou_nor_vz3
    
    return Kg_nor_ex1,Kg_log_ex1,Kg_nor_vz1,Kg_log_vz1

############# (4) Calibration coefficient program of measured data #############
def rou_correction_3(I,Re,Te_1,Te_2,T,epsilon_r,S0_ex):##Re_ex,T_ex;Re_vz,T_vz

    H1=np.array([float('inf')])#formation thickness
    S1=np.array([S0_ex])## Background conductivity (Note: S1 in this case is a one-dimensional array, not a scalar)
    import TEM_correction
    
    Kg_nor_ex1=np.zeros([len(T),Re.shape[0]])#Matrix to be assigned
    Kg_nor_vz1=np.zeros([len(T),Re.shape[0]])
    Kg_log_ex1=np.zeros([len(T),Re.shape[0]])
    Kg_log_vz1=np.zeros([len(T),Re.shape[0]])
    
    for j in np.arange(Re.shape[0]):##Calculate the resistivity background values at different measuring points
        Re1=Re[j,:].reshape(1,-1)
        S11=np.array([S1[0,j]])
        Kg_nor_ex1c,Kg_log_ex1c,Kg_nor_vz1c,Kg_log_vz1c=TEM_correction.rou_correction_1(I,Re1,Te_1,Te_2,T,H1,S11,epsilon_r)
        Kg_nor_ex1[:,j]=Kg_nor_ex1c.reshape(-1)
        Kg_log_ex1[:,j]=Kg_log_ex1c.reshape(-1)
        Kg_nor_vz1[:,j]=Kg_nor_vz1c.reshape(-1)
        Kg_log_vz1[:,j]=Kg_log_vz1c.reshape(-1)
    
    return Kg_nor_ex1,Kg_log_ex1,Kg_nor_vz1,Kg_log_vz1##Return correction factor