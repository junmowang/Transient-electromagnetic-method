# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:01:36 2024

@author: zhu_y
"""

#### Reading program of measured electromagnetic field data
### The program contains: (1) Read data from txt file; (2) Calculation of apparent resistivity; (3) Calibration; (4) time-depth transformation

from scipy import interpolate
import numpy as np
import pandas as pd
#import sympy as sy
# import scipy as sc
from scipy import constants
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import math
# import cmath
import time

plt.rcParams['font.sans-serif']=['Microsoft YaHei'] #It is used to display Chinese labels properly

## Function that reads the contents of a file such as txt and converts it to a two-dimensional list with the list element as a string
def read_txt_file_to_2d_array(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            # Use newlines to split the text content to generate a two-dimensional array
            data_lines = content.split('\n')
            data_2d_array = [line.split() for line in data_lines]
            return data_2d_array
    except FileNotFoundError:
        return f" '{file_path}' not found。"

## Function to convert a string from a list into data, and the input can be a one-dimensional or two-dimensional list
def convert_to_float_list(input_list):
    try:
        # Initialize an empty list to store the float values
        float_list = []

        # Check if the input_list is a 2D list (list of lists)
        if all(isinstance(sublist, list) for sublist in input_list):
            # Convert each item in each sublist to float and append to float_list
            for sublist in input_list:
                float_sublist = [float(item) for item in sublist]
                float_list.append(float_sublist)
        else:
            # Convert each item in the input_list to float and append to float_list
            float_list = [float(item) for item in input_list]

        return float_list  # Return the final float_list

    except ValueError:
        return None  # If any ValueError occurs, return None

## Function to extract geometry and electromagnetic field data read into the data
def read_parameter(data_B1_VZ,p):# p is the judgment condition to determine whether the input data is VZ or EX
    N=int(convert_to_float_list(data_B1_VZ[12])[0])#Number of measured points
    I=convert_to_float_list(data_B1_VZ[18])[-1]#Read supply current
    Te_1=np.array(convert_to_float_list(data_B1_VZ[21])).reshape(1,-1)#The left endpoint term coordinates of the source
    Te_2=np.array(convert_to_float_list(data_B1_VZ[22])).reshape(1,-1)#The right endpoint term coordinates of the source
    Re=np.array(convert_to_float_list(data_B1_VZ[26:26+N]))[:,:-1]#Observation point coordinates
    NT=int(convert_to_float_list(data_B1_VZ[65])[2])#Observation time number
    T=np.array(convert_to_float_list(data_B1_VZ[68:68+NT]))[:,1]*10**-3#Observation time

    if p=='EX':
        loc=2
        E=10**-6
    elif p=='VZ':
        loc=3
        E=10**-9/constants.mu_0
    else:
        print("Condition p is incorrectly entered !!!")
    vz=np.zeros([NT,N])
    for j in np.arange(N):
        vz[:,j]=np.array(convert_to_float_list(data_B1_VZ[68+j*(NT+5):68+j*(NT+5)+NT]))[:,loc]
    vz=vz*E
    return vz,I,Re,Te_1,Te_2,T

## Coordinate transformation function, the source is placed on the X-axis, and the measurement point is transformed accordingly
def geometrical_parameter(Re,Te_1,Te_2):#The coordinates of the observation point, the left endpoint of the source, and the right endpoint of the source
    Re[:,2]=0#Let z=0, that is, both the source point and the field point are on the ground
    Te_1[:,2]=0
    Te_2[:,2]=0
      
    di=np.array([1,0,0])#Unit vector along the X-axis
    dL=Te_2-Te_1#A vector in the positive direction of the source
    cos_a=np.dot(dL,di)[0]/np.linalg.norm(dL)#Cosine of the source vector
    sin_a=np.linalg.norm(np.cross(dL,di))/np.linalg.norm(dL)#Sine of the source vector
    A=np.array([[cos_a,-sin_a,0.0],[sin_a,cos_a,0.0],[0.0,0.0,0.0]])#Rotation matrix
    
    ## Translation source point, measuring point
    Te_ou=(Te_1+Te_2)/2
    Te_1_a=Te_1-Te_ou#Left endpoint translation
    Te_2_a=Te_2-Te_ou#Right endpoint translation
    Re_a=Re-Te_ou#Translation of observation point
    
    ## Rotating source point, measuring point
    Te_1_b=np.dot(A,Te_1_a.T).reshape(1,-1)#Left end rotation
    Te_2_b=np.dot(A,Te_2_a.T).reshape(1,-1)#Right end rotation
    Re_b=np.zeros(Re_a.shape)
    for j in np.arange(Re.shape[0]):
        Re_b[j,:]=np.dot(A,Re_a[j,:])#Rotation of observation point
    
    return Re_b,Te_1_b,Te_2_b


############################### （1）Read measured data ###############################
## The functions of the effective data ex, I, Re, Te_1, Te_2 and T are extracted from the measured data

on=2## Criteria: on=1, first-line data; on=2, second-line data

if on==1:
    ##### Pay attention to file path ！！！
    data_B1_VZ=read_txt_file_to_2d_array("C:\\Users\\zhu_y\\Desktop\\快速成像\\某地实测数据\\B1-VZ.txt")#First line vz data
    data_L1_EX=read_txt_file_to_2d_array("C:\\Users\\zhu_y\\Desktop\\快速成像\\某地实测数据\\L1-EX.txt")#First line ex data
    
    vz,I_vz,Re_vz,Te_1_vz,Te_2_vz,T_vz=read_parameter(data_B1_VZ,'VZ')#First line vz data
    ex,I_ex,Re_ex,Te_1_ex,Te_2_ex,T_ex=read_parameter(data_L1_EX,'EX')#First line ex data
elif on==2:
    data_B2_VZ=read_txt_file_to_2d_array("C:\\Users\\zhu_y\\Desktop\\快速成像\\某地实测数据\\B2-VZ.txt")#Second line vz data
    data_L2_EX=read_txt_file_to_2d_array("C:\\Users\\zhu_y\\Desktop\\快速成像\\某地实测数据\\L2-EX.txt")#Second line ex data

    vz,I_vz,Re_vz,Te_1_vz,Te_2_vz,T_vz=read_parameter(data_B2_VZ,'VZ')#Second line vz data
    ex,I_ex,Re_ex,Te_1_ex,Te_2_ex,T_ex=read_parameter(data_L2_EX,'EX')#Second line ex data


## When the source line is placed on the X-axis, the coordinates of the measuring point change accordingly
Re_vz,Te_1_vz,Te_2_vz=geometrical_parameter(Re_vz,Te_1_vz,Te_2_vz)
Re_ex,Te_1_ex,Te_2_ex=geometrical_parameter(Re_ex,Te_1_ex,Te_2_ex)

## Define some general parameters
I=I_ex
Te_1=Te_1_ex
Te_2=Te_2_ex
T=T_ex

############################### （2）Calculate apparent resistivity ###############################
import TEM_Ex_Vz_rou#Import the program for calculating apparent resistivity
##ex peak time imaging
rou_log_ex,rou_nor_ex,TT_ex=TEM_Ex_Vz_rou.rou_ex(ex,I,Re_ex,Te_1,Te_2,T)

plt.figure(dpi=250)
plt.loglog(TT_ex,abs(rou_nor_ex))
plt.loglog(TT_ex,abs(rou_log_ex),'-.')
plt.title('Electric field Ex - total apparent resistivity')

##vz peak time imaging
rou_log_vz,rou_nor_vz,TT_vz=TEM_Ex_Vz_rou.rou_vz(vz,I,Re_vz,Te_1,Te_2,T)

plt.figure(dpi=250)
plt.loglog(TT_vz,abs(rou_nor_vz))
plt.loglog(TT_vz,abs(rou_log_vz),'-.')
plt.title('Magnetic field Vz - total apparent resistivity')

############################### （3）Perform resistivity correction ###############################
epsilon_r=np.array([1.0,1.0,1.0])#The relative dielectric constant of each layer is 1.0

## Calculate the exact background resistivity ex
import TEM_background_resistivity#Import the program for calculating background apparent resistivity
## Electric field ex background conductivity value (select S0_ex)
S0_ex,__=TEM_background_resistivity.background_resistivity_1(I,Re_ex,Te_1,Te_2,T_ex,epsilon_r,rou_log_ex)## ex background field value; vz background field value

import TEM_correction#Import the program for calculating the apparent resistivity correction factor

Kg_nor_ex1,Kg_log_ex1,__,__=TEM_correction.rou_correction_3(I,Re_ex,Te_1,Te_2,T_ex,epsilon_r,S0_ex)

rou_log_ex4=Kg_log_ex1*rou_log_ex##Corrected ex resistivity (logarithmic domain)
rou_nor_ex4=Kg_nor_ex1*rou_nor_ex#(linear domain)

## Calculate accurate background resistivity vz
import TEM_background_resistivity
## Electric field ex background conductivity value (select S0_vz)
__,S0_vz=TEM_background_resistivity.background_resistivity_1(I,Re_vz,Te_1,Te_2,T_vz,epsilon_r,rou_log_vz)## ex background field value; vz background field value

import TEM_correction

__,__,Kg_nor_vz1,Kg_log_vz1=TEM_correction.rou_correction_3(I,Re_vz,Te_1,Te_2,T_vz,epsilon_r,S0_vz)

rou_log_vz4=Kg_log_vz1*rou_log_vz##The resistivity of vz after correction
rou_nor_vz4=Kg_nor_vz1*rou_nor_vz#(linear domain)

## Take the absolute value of the resistivity
rou_log_ex4=abs(rou_log_ex4)
rou_nor_ex4=abs(rou_nor_ex4)
rou_log_vz4=abs(rou_log_vz4)
rou_nor_vz4=abs(rou_nor_vz4)


plt.figure(dpi=250)
plt.loglog(T,abs(rou_nor_ex4))
plt.loglog(T,abs(rou_log_ex4),'-.')
plt.title('Electric field Ex -- lifetime apparent resistivity (after correction)')
plt.legend(['[0,350]','[0,500]','[0,1000]','[0,2000]'])
plt.figure(dpi=250)
plt.loglog(T,abs(rou_nor_vz4))
plt.loglog(T,abs(rou_log_vz4),'-.')
plt.title('Magnetic field Vz - Lifetime apparent resistivity (corrected)')
plt.legend(['[0,350]','[0,500]','[0,1000]','[0,2000]'])

############################### （4）Time-depth transformation  ###############################
import TEM_time_depth#Import the calculating time depth conversion program

## Calculate depth (Chen), depth (eaton)
chen_h1_ex,eaton_h1_ex=TEM_time_depth.T_H_2(rou_log_ex4,T,Te_1,Te_2)#Calculated depth of corrected data
chen_h1_vz,eaton_h1_vz=TEM_time_depth.T_H_2(rou_log_vz4,T,Te_1,Te_2)#Calculated depth of corrected data


plt.figure(dpi=250)
plt.loglog(chen_h1_ex,rou_log_ex4)
plt.title('Total apparent resistivity (ex)- Depth - chen')

plt.figure(dpi=250)
plt.loglog(eaton_h1_ex,rou_log_ex4)
plt.title('Total apparent resistivity (ex)- Depth -eaton')

plt.figure(dpi=250)
plt.loglog(chen_h1_vz,rou_log_vz4)
plt.title('Total apparent resistivity (vz)- Depth - chen')

plt.figure(dpi=250)
plt.loglog(eaton_h1_vz,rou_log_vz4)
plt.title('Total apparent resistivity (vz)- Depth -eaton')
