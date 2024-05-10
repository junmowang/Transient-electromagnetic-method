# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 19:56:53 2024

@author: zhu_y
"""

### Rapid imaging - full time apparent resistivity calculation program
## Three methods: ex Peak MomentVZ Peak moment
## Includes: (1) two derivative methods of the program; (2) Two methods of resistivity procedures

# import xlrd
# from scipy.integrate import quad
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
# import numba as nb
#import functools

############################### (1) two derivative methods of the program ###############################
## Normal derivation
def nor_diff(ex,T,p):
    ex=-1*abs(ex)##absolute value
    def diff1(Tnew,ex1):#Discrete differentiation-forward difference formula
        K=Tnew.shape[0]#Calculate the number of arguments
        ext1=np.zeros((K-1,Tnew.shape[1]))#（Electric field）first derivative to be assigned matrix
        ext1=np.diff(ex1,axis=0)/np.diff(Tnew,axis=0)#forward difference
        Tnew1=Tnew[:-1,:]#Lost last row
        return ext1,Tnew1
    N=ex.shape[1]#Get the number of points
    TT=T.reshape(-1,1).repeat(N,axis=1)
    Tnew=np.logspace(np.log10(T.min()),np.log10(T.max()),p).reshape(-1,1).repeat(N,axis=1)#Encrypt observation time points
    ex1=np.zeros(Tnew.shape)#Create an interpolated crypto electric field
    for i in np.arange(N):#interpolation
        f=interpolate.interp1d(TT[:,i],ex[:,i],kind='cubic',fill_value="extrapolate")# linear; quadratic ; cubic(third order B-spline)
        ex1[:,i]=f(Tnew[:,i])#The electric field response after encryption

    ext1,Tnew1=diff1(Tnew,ex1)#Find the first derivative of the electric field
    
    loc_tpeak=np.argmax(ext1,axis=0)#The initial position of the maximum derivative value
    tpeak=np.zeros((1,N))#Matrix to be assigned at peak time
    dl=-5
    dr=6
    for i in np.arange(N):
        loc_t=Tnew1[loc_tpeak[i]+dl:loc_tpeak[i]+dr,i]
        tnew=np.linspace(loc_t.min(),loc_t.max(),p)#Time encryption
        loc_ext=ext1[loc_tpeak[i]+dl:loc_tpeak[i]+dr,i]
        h=interpolate.interp1d(loc_t,loc_ext,kind='cubic')# linear; quadratic ; cubic(third order B-spline)
        ext_new=h(tnew)
        loc=np.argmax(ext_new,axis=0)#The exact location of the maximum derivative
        tpeak[:,i]=tnew[loc]
    return ext1,Tnew1,tpeak#Encryption field, encryption time, peak time

## Take the derivative of the logarithm domain
def log_diff(vz,T,p):
    N=vz.shape[1]#Get the number of points
    TT=T.reshape(-1,1).repeat(N,axis=1)
    Tnew=np.logspace(np.log10(T.min()),np.log10(T.max()),p).reshape(-1,1).repeat(N,axis=1)#Encrypt observation time points
    
    log_vz=np.log10(abs(vz))#Logarithm of the initial electric field data (absolute value)
    log_TT=np.log10(TT)#Log the initial time data
    log_Tnew=np.log10(Tnew)#Logarithm of the encryption time
    
    log_vz1=np.zeros(Tnew.shape)#Create an interpolated crypto electric field
    log_vzt1=np.zeros([Tnew.shape[0],Tnew.shape[1]])#Create an encrypted first-order electric field after interpolation
    log_vzt2=np.zeros([Tnew.shape[0],Tnew.shape[1]])#Create an interpolated encrypted second-order electric field
    for i in np.arange(N):#interpolation
        f=CubicSpline(log_TT[:,i],log_vz[:,i])#The logarithmic domain is interpolated
        log_vz1[:,i]=f(log_Tnew[:,i])#The electric field response after encryption
        log_vzt1[:,i]=f(log_Tnew[:,i],1)#First order response of the electric field after encryption
        log_vzt2[:,i]=f(log_Tnew[:,i],2)#The second order response of the encrypted electric field
    vz1=-1*10**(log_vz1)#Take the logarithm and you get the electric field
    vzt1=vz1*log_vzt1/Tnew#We get a first-order electric field
    
    loc_tpeak=np.argmax(vzt1,axis=0)#The initial position of the maximum derivative value
    tpeak=np.zeros((1,N))#Matrix to be assigned at peak time
    dl=-5
    dr=6
    for i in np.arange(N):
        loc_t=Tnew[loc_tpeak[i]+dl:loc_tpeak[i]+dr,i]
        tnew=np.logspace(np.log10(loc_t.min()),np.log10(loc_t.max()),p)#Time encryption  
        loc_vzt=vzt1[loc_tpeak[i]+dl:loc_tpeak[i]+dr,i]
        h=CubicSpline(loc_t,loc_vzt)
        vzt_new=h(tnew)
        loc=np.argmax(vzt_new,axis=0)#The exact location of the maximum derivative
        tpeak[:,i]=tnew[loc]
    return vzt1,Tnew,tpeak#Encryption magnetic field, encryption time, peak time

############################### (2) Two methods of resistivity procedures ###############################
## Method 1: Based on ex peak moment
def rou_ex(ex,I,Re,Te_1,Te_2,T):
    TT=T.reshape(-1,1).repeat(len(Re),axis=1)
    dx=np.linalg.norm(Te_1-Te_2)#wire length
    Te_0=(Te_1+Te_2)/2#Midpoint coordinates of long wire
    R=np.zeros(len(Re))#Range to be assigned matrix (vector)
    for i in np.arange(len(Re)):#Calculated transmit-receive distance
        R[i]=np.linalg.norm((Re-Te_0)[i])#Transmit-receive distance
    RR=R.reshape(1,-1).repeat(len(T),axis=0)
    
    #Total apparent resistivity
    p=500#Cryptographic number
    
    ext1_nor,Tnew_nor,tpeak_nor=nor_diff(ex,T,p)
    tao_nor=TT/tpeak_nor#Calculated τ
    csmath_erf=np.frompyfunc(math.erf,1,1)#Error function erf vectorization in math
    u1=0.5*(10/tao_nor)**0.5
    erf1=csmath_erf(u1).astype('float64')
    rou_nor_ex=-np.pi**1.5/(I*dx)*RR[0,:]**3*ex/(np.pi**0.5/2*erf1-u1*np.exp(-u1**2))#Under normal circumstances, ex calculations are used
    
    ext1_log,Tnew_log,tpeak_log=log_diff(ex,T,p)
    tao_log=TT/tpeak_log#Calculated τ
    u2=0.5*(10/tao_log)**0.5
    erf2=csmath_erf(u2).astype('float64')
    rou_log_ex=-np.pi**1.5/(I*dx)*RR[0,:]**3*ex/(np.pi**0.5/2*erf2-u2*np.exp(-u2**2))#In the logarithm domain, ex is used
    
    return rou_log_ex,rou_nor_ex,TT

## Method 2: Peak time based on vz
def rou_vz(vz,I,Re,Te_1,Te_2,T):
    TT=T.reshape(-1,1).repeat(len(Re),axis=1)
    dx=np.linalg.norm(Te_1-Te_2)#wire length
    Te_0=(Te_1+Te_2)/2#Midpoint coordinates of long wire
    R=np.zeros(len(Re))#Range to be assigned matrix (vector)
    for i in np.arange(len(Re)):#Calculated transmit-receive distance
        R[i]=np.linalg.norm((Re-Te_0)[i])#Transmit-receive distance
    RR=R.reshape(1,-1).repeat(len(T),axis=0)
    YY=Re[:,1].reshape(1,-1).repeat(len(T),axis=0)
    
    #Total apparent resistivity
    p=500#Cryptographic number
    mu0=constants.mu_0
    
    vzt1_nor,Tnew_nor,tpeak_nor=nor_diff(vz,T,p)
    tao_nor=TT/tpeak_nor#Calculated τ
    csmath_erf=np.frompyfunc(math.erf,1,1)#Error function erf vectorization in math
    u1=0.5*(14/tao_nor)**0.5
    erf1=csmath_erf(u1).astype('float64')
    rou_nor_vz=np.pi**1.5*mu0/(I*dx)*RR[0,:]**5/YY[0,:]*abs(vz)/(3*np.pi**0.5/2*erf1-u1*np.exp(-1*u1**2)*(3+2*u1**2))#Under normal circumstances, vz calculations are used
    
    
    vzt1_log,Tnew_log,tpeak_log=log_diff(vz,T,p)
    tao_log=TT/tpeak_log#Calculated τ
    u2=0.5*(14/tao_log)**0.5
    erf2=csmath_erf(u2).astype('float64')
    rou_log_vz=np.pi**1.5*mu0/(I*dx)*RR[0,:]**5/YY[0,:]*abs(vz)/(3*np.pi**0.5/2*erf2-u2*np.exp(-1*u2**2)*(3+2*u2**2))#In the logarithm domain, vz is used
    
    return rou_log_vz,rou_nor_vz,TT