#!/usr/bin/env python
# -*- coding: utf-8 -*-
#This is the depth only file,the Factor is most important result he
import os
import csv
import sys, os, shutil, re, math, time
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import smtplib

def CsvName(rootDir):
      list_dirs=os.walk(rootDir)
      filenames=[]
      for root,dirs,files in list_dirs:
          for f in files:
              filenames.append (os.path.join(rootDir,f))
      return filenames

def format(name):
    b=name
    a=r''
    c=a.join(b)
    return c
      
def Topasdata(filenames):   
    i=0
    A=[0]
    TopasEnergy=[] 
    with open(filenames) as idd:
        result= csv.reader(idd)
        for line in result:
            A.append(line)
    idd.closed
    n=0
    X=[]
    del A[0:9]
    while (n<len(A)):                  
          TopasEnergy.append(float(A[n][3])*10**10)
          X.append(n)          
          n=n+1
    Dose=TopasEnergy     
    return X,Dose

def AllTopasName(Energy,Snot):
    TopasPath='F:/Monte Carlo/SingleTopas/All/'+str(Snot)+'cm/'
    EnergyList=[151.0, 153.2, 155.3, 157.4, 159.5, 161.6, 163.9, 166.2, 168.8, 171.3, 173.7, 176.2, 178.6, 181.1, 183.4, 185.8, 188.2, 190.5, 193.0, 195.6, 198.3, 201.0, 203.7, 206.3, 209.0, 211.0, 211.6, 214.2, 216.7, 219.3, 221.8]
    if float(Energy) in EnergyList:
       position=EnergyList.index(float(Energy))
    else:
        print 'not in list'
        print Energy
        position=0
    if position>9:
        S1filename=TopasPath+'151.0_S1_Run_00'+str(position)+'.csv'
    else:
        S1filename=TopasPath+'151.0_S1_Run_000'+str(position)+'.csv'  
    return S1filename

        
def SinLossRate(Energy,Snot,Depth):
    if Energy<151 or Snot>38:
        print 'no Absorber data'   
    SnotList=[0,23,28,33,38]
    n=0
    while Snot-SnotList[n]>=0:
          n=n+1 
    Snotmax=SnotList[n]
    Snotmin=SnotList[n-1]    
    [X0,Y0max]=Topasdata(AllTopasName(Energy,0))
    [X0,Y0min]=Topasdata(AllTopasName(Energy,0))    
    Peak=Y0max.index(max(Y0max))   
    Y0max=np.array(Y0max[67:Peak])
    Y0min=np.array(Y0min[67:Peak])    
    [X,Ymax]=Topasdata(AllTopasName(Energy,Snotmax))  
    [X,Ymin]=Topasdata(AllTopasName(Energy,Snotmin))    
    Peak=Ymax.index(max(Ymax))    
    X=X[0:Peak]   
    Ymax=np.array(Ymax[0:Peak])
    Ymin=np.array(Ymin[0:Peak])                                
    
    Lossmax=(Y0max-Ymax)/Y0max
    Lossmin=(Y0min-Ymin)/Y0min   
    if Depth*10<len(X):
          position=X.index(Depth*10)
          LossRate=np.linspace(Lossmin[position],Lossmax[position],5)              
          return LossRate[int(Snot-Snotmin)]
    else:
          return 0
#Main
s_time = time.time()
PathTopas="/home/globaloncology/Desktop/topas/work/SingleTopas/All/"
bin=[1,1,350]
Energy=[151.0, 153.2, 155.3, 157.4, 159.5, 161.6, 163.9, 166.2, 168.8, 171.3, 173.7, 176.2, 178.6, 181.1, 183.4, 185.8, 188.2, 190.5, 193.0, 195.6, 198.3, 201.0, 203.7, 206.3, 209.0, 211.0, 211.6, 214.2, 216.7, 219.3, 221.8]
#MakeTopasFile(PathTopas,Energy,bin)
#FigureAll(Energy)
#LossRate(Energy)
SinLossRate(221.8,25,10)
#Email()
ttime = (time.time() - s_time)/60
print ttime,'mins'