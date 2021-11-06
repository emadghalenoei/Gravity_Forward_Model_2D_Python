# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 12:11:30 2020

@author: emadg
"""
import numpy as np
import math

def ModelSpace(xs,ys,Z0,ZEND,dZ,Pad_Length):

    Azimuth = math.atan2(xs[-1]-xs[0],ys[-1]-ys[0])
    dis = np.sqrt((xs-xs[0])**2 + (ys-ys[0])**2)
    z = np.arange(Z0+dZ/2, ZEND, dZ)
    N = math.ceil(len(z)/len(dis))
    deltaDis = abs((dis[1]-dis[0])/(N-1))
    N = math.ceil((dis[-1]-dis[0])/deltaDis)
    pad_num = math.ceil(Pad_Length/deltaDis);
    
    x = np.zeros(N+(2*pad_num))
    y = np.zeros(N+(2*pad_num))
    dismodel = np.zeros(N+(2*pad_num))
    
    x[0] = xs[0] - pad_num*deltaDis*math.sin(Azimuth)
    y[0] = ys[0] - pad_num*deltaDis*math.cos(Azimuth) 
    dismodel[0] = dis[0] - pad_num*deltaDis
    
    for i in range(N+(2*pad_num)):
        x[i] = x[0]+ i*deltaDis*math.sin(Azimuth)
        y[i] = y[0]+ i*deltaDis*math.cos(Azimuth)
        dismodel[i] = dismodel[0] +(i)*deltaDis

    X, Z = np.meshgrid(x, z)
    Y, Z = np.meshgrid(y, z)
    DISMODEL, Z = np.meshgrid(dismodel, z)
    return DISMODEL,X,Y,Z