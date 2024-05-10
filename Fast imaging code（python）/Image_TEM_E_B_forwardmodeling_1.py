# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 10:14:29 2024

@author: zhu_y
"""

#### Forward transient electromagnetic program
### Supply current: negative step response
### Underground medium: uniform half-space, one-dimensional layered medium
### The field source is: dipole, long wire
### Observed components: Ex, Vz
### Includes: (1) One-dimensional, uniform forward programming


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

############################################  (1) One-dimensional, uniform forward programming ############################################
def forwardmodeling(I,Re,Te_1,Te_2,T,H,S,epsilon_r,judge):
    mu0=constants.mu_0#permeability of vacuum
    epsilon0=constants.epsilon_0#permittivity of vacuum
  
    if judge==1:##After the judgment is successful, the electric dipoles in the uniform half space are calculated
        #Field and source parameters
        ds=np.linalg.norm(Te_1-Te_2)#The length of the dipole
        Te_0=(Te_1+Te_2)/2#Midpoint coordinates of long wire
        R=np.zeros(len(Re))#Range to be assigned matrix (vector)
        for i in np.arange(len(Re)):#Calculated transmit-receive distance
            R[i]=np.linalg.norm((Re-Te_0)[i])#Transmit-receive distance
        #Uniform half space vz component
        TT=T.reshape(-1,1).repeat(len(Re),axis=1)
        RR=R.reshape(1,-1).repeat(len(T),axis=0)
        YY=Re[:,1].reshape(1,-1).repeat(len(T),axis=0)
        u=0.5*(mu0*S[0]/TT)**0.5*RR
        csmath=np.frompyfunc(math.erf,1,1)#Error function erf vectorization in math
        vz=I*ds/(np.pi**1.5*S[0]*mu0)*YY/RR**5*(3*np.pi**0.5/2*csmath(u).astype('float64')-u*np.exp(-1*u**2)*(3+2*u**2))
        #Uniform half space ex component
        ex=-1*I*ds/(np.pi**1.5*S[0])*RR**-3*(np.pi**0.5/2*csmath(u).astype('float64')-u*np.exp(-1*u**2))
        
    elif judge==2:##After successful judgment, the long wire of the layered medium is calculated
        ##One-dimensional layered medium model parameters
        H1=H#formation thickness
        S1=S#Electrical conductivity of each layer

        epsilon1=epsilon_r*epsilon0#Dielectric constant of each layer 
        
        ##Numerical filter parameter
        #The 120-point J0 linear digital filtering parameters of D. Geuptasarma and B. Shingh
        n1=1
        N1=120
        a1=-8.38850000000e+00
        s1= 9.04226468670e-02
        #The 140-point J1 linear digital filtering parameters of D.Guptasarma and B.Shingh
        n2=1
        N2=140
        a2=-7.91001919000e+00
        s2=8.79671439570e-02
        #Wang Huajun's 250 point cosine (sine) filtering parameter
        nc=-149
        NC=100
        delta=math.log(10)/20
        
        #Import filter coefficient
        ## The read path is adjusted according to the location of the parameter file ！！！！！！！！！
        original_data0=pd.read_excel("C:\\Users\\zhu_y\\Desktop\\filter coefficient\\120-J0.xlsx")#Read the 120-point J0 filter coefficient
        original_data1=pd.read_excel("C:\\Users\\zhu_y\\Desktop\\filter coefficient\\140-J1.xlsx")#Read the 140-point -J1 filter coefficient
        original_data2=pd.read_excel("C:\\Users\\zhu_y\\Desktop\\filter coefficient\\250-Cos.xlsx")#Read the 250 point-cosine filter coefficient
        original_data3=pd.read_excel("C:\\Users\\zhu_y\\Desktop\\filter coefficient\\250-Sin.xlsx")#Read 250 points - sinusoidal filter coefficient
        
        np_data0=original_data0.to_numpy()#Format conversion, converting DataFrame to ndarray
        np_data1=original_data1.to_numpy()
        np_data2=original_data2.to_numpy()
        np_data3=original_data3.to_numpy()
        
        re_data0=np_data0.reshape(-1)#Array reorganization, turning a two-dimensional array into a one-dimensional array (dimensionality reduction operation)
        re_data1=np_data1.T.reshape(-1)
        re_data2=np_data2.T.reshape(-1)
        re_data3=np_data3.T.reshape(-1)
        
        W1=re_data0.reshape(1,-1)#120 points -J0 filter coefficient
        W2=re_data1.reshape(1,-1)#140 points -J1 filter factor
        C=re_data2.reshape(1,-1)#250 point-cosine filter coefficient
        C1=re_data3.reshape(1,-1)#250 points - sinusoidal filtering coefficient
        
        ##Long conductor parameter
        Te_0=(Te_1+Te_2)/2#Midpoint coordinates of long wire
        R=np.zeros(len(Re))#Range to be assigned matrix (vector)
        for i in np.arange(len(Re)):#Field point and conductor center point, transmit-receive distance to be assigned matrix (vector
            R[i]=np.linalg.norm((Re-Te_0)[i])#Transmit-receive distance
        L=np.linalg.norm(Te_1-Te_2)#wire length
        NL=10#Number of dipoles divided by long wire
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
        R1=np.zeros(Re.shape[0])#Left endpoint entry, transmit-receive distance to be assigned matrix
        for i in np.arange(Re.shape[0]):
            R1[i]=np.linalg.norm(Re[i,:]-Te_1)
        R2=np.zeros(Re.shape[0])#Right endpoint entry, transmit-receive distance to be assigned matrix
        for i in np.arange(Re.shape[0]):
            R2[i]=np.linalg.norm(Re[i,:]-Te_2)
        
        #One-dimensional layered medium ex component
        def th(x):#hyperbolic tangent function
            y=np.e**(-2*x)
            z=(1-y)/(1+y)
            return z
        def Y1(H,S,epsilon,mu0,lamda,wou):
            N=len(H)#Recorded stratigraphic number
            for i in np.arange(N,0,-1):
                if i==N:
                    Y=(lamda**2+1j*wou*mu0*S[i-1]-wou**2*mu0*epsilon[i-1])**0.5/(1j*wou*mu0)
                else:
                    u=(lamda**2+1j*wou*mu0*S[i-1]-wou**2*mu0*epsilon[i-1])**0.5
                    z=1j*wou*mu0
                    Y0=u/z
                    Y=(Y0)*((Y+Y0*th(u*H[i-1]))/(Y0+Y*th(u*H[i-1])))
            return (Y)
        def Z1(H,S,epsilon,mu0,lamda,wou):
            N=len(H)
            for i in np.arange(N,0,-1):
                if i==N:
                    Z=(lamda**2+1j*wou*mu0*S[i-1]-wou**2*mu0*epsilon[i-1])**0.5/(1j*wou*epsilon[i-1]+S[i-1])
                else:
                    u=(lamda**2+1j*wou*mu0*S[i-1]-wou**2*mu0*epsilon[i-1])**0.5
                    y=1j*wou*epsilon[i-1]+S[i-1]
                    Z0=u/y
                    Z=(Z0)*((Z+Z0*th(u*H[i-1]))/(Z0+Z*th(u*H[i-1])))
            return (Z)
        
        ######### Calculated electric field Ex
        #Line source part
        RR0=R0.T.reshape(-1,1).repeat(W1.shape[1],axis=1).reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)
        ii0=np.arange(n1,N1+1).reshape(-1,1).repeat(R0.shape[0]*R0.shape[1],axis=1).T.reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)
        WW0=W1.reshape(-1,1).repeat(R0.shape[0]*R0.shape[1],axis=1).T.reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)
        lamda0=10**(a1+(ii0-1)*s1)/RR0
        TT0=T.reshape(-1,1).repeat(C.shape[1],axis=1).reshape(1,-1).repeat(W1.shape[1]*R0.shape[1]*R0.shape[0],axis=0)
        nn0=np.arange(nc,NC+1).reshape(-1,1).T.repeat(len(T),axis=0).reshape(1,-1).repeat(W1.shape[1]*R0.shape[1]*R0.shape[0],axis=0)
        CC0=C.reshape(-1,1).T.repeat(len(T),axis=0).reshape(1,-1).repeat(W1.shape[1]*R0.shape[1]*R0.shape[0],axis=0)
        wou0=np.exp(nn0*delta)/TT0
        Y00=Y1(H1,S1,epsilon1,mu0,lamda0,wou0)####### Time consuming 60%
        K0=1J*lamda0/(lamda0+1j*wou0*mu0*Y00)
        Ex0_0=np.imag(K0)/RR0/TT0*WW0*CC0*ds*mu0##Line source term - time-frequency transformation
        Ex0_1=np.zeros([R0.shape[1],len(T)])#Line source electric field EX response
        for k in np.arange(len(T)):
            for j in np.arange(R0.shape[1]):
                Ex0_1[j,k]=np.sum(Ex0_0[j*W1.shape[1]*R0.shape[0]:(j+1)*W1.shape[1]*R0.shape[0],k*C.shape[1]:(k+1)*C.shape[1]])
        Ex0=Ex0_1.T#Line source electric field
        
        #Left and right end part
        RR1=R1.reshape(-1,1).repeat(W2.shape[1],axis=1).reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)#left end point
        RR2=R2.reshape(-1,1).repeat(W2.shape[1],axis=1).reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)#right end point
        X1=(Re[:,0]-Te_1[:,0]).reshape(-1,1).repeat(W2.shape[1],axis=1).reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)
        X2=(Re[:,0]-Te_2[:,0]).reshape(-1,1).repeat(W2.shape[1],axis=1).reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)
        ii12=np.arange(n2,N2+1).reshape(-1,1).repeat(R0.shape[1],axis=1).T.reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)
        WW12=W2.reshape(-1,1).repeat(R0.shape[1],axis=1).T.reshape(-1,1).repeat(C.shape[1]*len(T),axis=1)
        lamda1=10**(a2+(ii12-1)*s2)/RR1
        lamda2=10**(a2+(ii12-1)*s2)/RR2
        TT12=T.reshape(-1,1).repeat(C.shape[1],axis=1).reshape(1,-1).repeat(W2.shape[1]*R0.shape[1],axis=0)
        nn12=np.arange(nc,NC+1).reshape(-1,1).T.repeat(len(T),axis=0).reshape(1,-1).repeat(W2.shape[1]*R0.shape[1],axis=0)
        en12=np.exp(nn12*delta)

        CC12=C.reshape(-1,1).T.repeat(len(T),axis=0).reshape(1,-1).repeat(W2.shape[1]*R0.shape[1],axis=0)#Cosine filtering coefficient
        CC23=C1.reshape(-1,1).T.repeat(len(T),axis=0).reshape(1,-1).repeat(W2.shape[1]*R0.shape[1],axis=0)#Sinusoidal filtering coefficient
        
        wou12=np.exp(nn12*delta)/TT12
        Y11=Y1(H1,S1,epsilon1,mu0,lamda1,wou12)####### Time consuming 7%
        Z11=Z1(H1,S1,epsilon1,mu0,lamda1,wou12)####### Time consuming 7%
        Y22=Y1(H1,S1,epsilon1,mu0,lamda2,wou12)####### Time consuming 7%
        Z22=Z1(H1,S1,epsilon1,mu0,lamda2,wou12)####### Time consuming 7%
        K1=lamda1*Z11/(lamda1+1j*wou12*epsilon0*Z11)-1/(lamda1/(1j*wou12*mu0)+Y11)
        K2=lamda2*Z22/(lamda2+1j*wou12*epsilon0*Z22)-1/(lamda2/(1j*wou12*mu0)+Y22)
        
        ##Endpoint term - time-frequency transformation
        # Ex1_0=np.imag(K1)*X1*WW12*CC12/en12/RR1**2#cosine transform
        # Ex2_0=np.imag(K2)*X2*WW12*CC12/en12/RR2**2#
        Ex1_0=np.real(K1)*X1*WW12*CC23/en12/RR1**2#sine transform
        Ex2_0=np.real(K2)*X2*WW12*CC23/en12/RR2**2#
        
        Ex1_1=np.zeros([R0.shape[1],len(T)])
        Ex2_1=np.zeros([R0.shape[1],len(T)])
        for k in np.arange(len(T)):
            for j in np.arange(R0.shape[1]):
                Ex1_1[j,k]=np.sum(Ex1_0[j*W2.shape[1]:(j+1)*W2.shape[1],k*C.shape[1]:(k+1)*C.shape[1]])
                Ex2_1[j,k]=np.sum(Ex2_0[j*W2.shape[1]:(j+1)*W2.shape[1],k*C.shape[1]:(k+1)*C.shape[1]])
        Ex1=Ex1_1.T#Left terminal electric field
        Ex2=Ex2_1.T#Right end electric field
        
        ## Total electric field
        # Ex_long=-I/(2**0.5*np.pi**1.5)*(Ex0+Ex1-Ex2)#cosine transform
        Ex_long=-I/(2**0.5*np.pi**1.5)*(Ex0+Ex1+Ex2)#sine transform
        ex=Ex_long
        
        ######### Calculate induced electromotive force Vz
        #Line source part
        RR0=R0.T.reshape(-1,1).repeat(W2.shape[1],axis=1).reshape(-1,1).repeat(C1.shape[1]*len(T),axis=1)
        ii0=np.arange(n2,N2+1).reshape(-1,1).repeat(R0.shape[0]*R0.shape[1],axis=1).T.reshape(-1,1).repeat(C1.shape[1]*len(T),axis=1)
        WW0=W2.reshape(-1,1).repeat(R0.shape[0]*R0.shape[1],axis=1).T.reshape(-1,1).repeat(C1.shape[1]*len(T),axis=1)
        lamda0=10**(a2+(ii0-1)*s2)/RR0
        TT0=T.reshape(-1,1).repeat(C1.shape[1],axis=1).reshape(1,-1).repeat(W2.shape[1]*R0.shape[1]*R0.shape[0],axis=0)
        nn0=np.arange(nc,NC+1).reshape(-1,1).T.repeat(len(T),axis=0).reshape(1,-1).repeat(W2.shape[1]*R0.shape[1]*R0.shape[0],axis=0)
        CC0=C1.reshape(-1,1).T.repeat(len(T),axis=0).reshape(1,-1).repeat(W2.shape[1]*R0.shape[1]*R0.shape[0],axis=0)
        wou0=np.exp(nn0*delta)/TT0
        ZZj=Re[:,2].reshape(-1,1).repeat(R0.shape[0]*W2.shape[1],axis=1).reshape(-1,1).repeat(C1.shape[1]*len(T),axis=1)
        YYj=Re[:,1].reshape(-1,1).repeat(R0.shape[0]*W2.shape[1],axis=1).reshape(-1,1).repeat(C1.shape[1]*len(T),axis=1)
        YYl=Te_ou[:,1].reshape(-1,1).repeat(W2.shape[1],axis=1).reshape(1,-1).repeat(R0.shape[1],axis=0).reshape(-1,1).repeat(C1.shape[1]*len(T),axis=1)
        Y00=Y1(H1,S1,epsilon1,mu0,lamda0,wou0)####### Time consuming 60%
        K0=2*lamda0**2*np.exp(lamda0*ZZj)/(lamda0+1j*wou0*mu0*Y00)## Kernel of Vz
        
        # Vz0_0=np.real(K0)/RR0**2*(YYj-YYl)*WW0*CC0*ds/TT0#cosine transform
        Vz0_0=np.imag(K0)/RR0**2*(YYj-YYl)*WW0*CC0*ds/TT0#sine transform
        
        Vz0_1=np.zeros([R0.shape[1],len(T)])#Line source electric field EX response
        for k in np.arange(len(T)):
            for j in np.arange(R0.shape[1]):
                Vz0_1[j,k]=np.sum(Vz0_0[j*W2.shape[1]*R0.shape[0]:(j+1)*W2.shape[1]*R0.shape[0],k*C1.shape[1]:(k+1)*C1.shape[1]])
        Vz0=Vz0_1.T#Line source electric field
        Vz_long=-I/(2*np.pi)**1.5*Vz0
        vz=Vz_long
        
    else:
        ex=np.nan
        vz=np.nan
        print ('Warning:  The judgment condition is incorrect')
    return ex,vz