# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 16:14:42 2020

@author: emad ghalenoei
"""

import numpy as np
import math
import os
import shutil
from ModelSpace import ModelSpace
from Gravity_Kernel_Expanded import Gravity_Kernel_Expanded
import matplotlib.pyplot as plt


plt.close('all')

fpath = os.getcwd()+'\Image'
if os.path.exists(fpath) and os.path.isdir(fpath):
    shutil.rmtree(fpath)
os.mkdir(fpath) 

fpath_loaddesk = os.getcwd()

Ndatapoints = 30;  # number of data points
xs = np.linspace(0,60000, Ndatapoints)  # define x of data points
ys = np.linspace(0,50000, Ndatapoints)  # define y of data points
dis_s = np.sqrt((xs-xs[0])**2 + (ys-ys[0])**2);  # compute distance from the first data point

# model space
Z0 = 0              # the shallower depth in model space  
ZEND = 10000        # the deepest depth in model space 
dZ = 100            # the width of prisms in z direction
Pad_Length = 6000   # padding length to minimize the edge effect

CX = 100            # number of prisms in x axis
CZ = 100            # number of prisms in z axis

Azimuth = math.atan2(xs[-1]-xs[0],ys[-1]-ys[0])
xmodel = np.linspace(xs[0]-Pad_Length*math.sin(Azimuth),xs[-1]+Pad_Length*math.sin(Azimuth),CX)
ymodel = np.linspace(ys[0]-Pad_Length*math.cos(Azimuth) ,ys[-1]+Pad_Length*math.cos(Azimuth) ,CX)
dismodel = np.linspace(dis_s[0]-Pad_Length,dis_s[-1]+Pad_Length,CX)
zmodel = np.linspace(Z0,ZEND,CZ)

X, Z = np.meshgrid(xmodel, zmodel)
Y, Z = np.meshgrid(ymodel, zmodel)
DISMODEL, Z = np.meshgrid(dismodel, zmodel)
    
    
dx=abs(X[0,1]-X[0,0]) 
dy=abs(Y[0,1]-Y[0,0])
dz = abs(Z[1,0]-Z[0,0])
dDis = abs(DISMODEL[0,1]-DISMODEL[0,0])

TrueDensityModel = np.load(fpath_loaddesk+'//'+'TrueDensityModel.npy')  #load or define your true model

Kernel_Grv = Gravity_Kernel_Expanded(DISMODEL,Z,dis_s)  # 2D gravity forward model kernel
Kernel_Grv = Kernel_Grv*1e8

   
dg_true = Kernel_Grv @ TrueDensityModel.flatten('F') # Unit(mGal)   # true gravity data without noise  

# Adding noise
noise_g_level = 0.05  #noise level
sigma_g_original=noise_g_level*max(abs(dg_true))  # standard deviation of noise
noise_g_original = sigma_g_original*np.random.randn(Ndatapoints) 
dg_simulated = dg_true + noise_g_original  # simulated and noisy gravity data

fig, axe = plt.subplots()

axe.plot(dis_s,dg_simulated, 'k-.',linewidth=2) #row=0, col=0
axe.plot(dis_s,dg_true, 'r--',linewidth=2) #row=0, col=0
axe.set(xlabel='X Profile (km)', ylabel='Gravity (mGal)')
plt.show()
figname = 'Gravity'
fignum = ''
fig.savefig(fpath+'/'+figname+str(fignum)+'.pdf')
plt.close(fig)    # close the figure window
    