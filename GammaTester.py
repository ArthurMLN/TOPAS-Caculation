# -*- coding: utf-8 -*-
###############################################################################
# CONSTANTS
###############################################################################
DICOM_FILE = "y:\\DICOM\\Rodriguez_ALA_3.9_TPS.dcm"
COUNTER_VAR = 0
OUTPUT_DIR = "/scratch"
OUTPUT_FOLDER='/home/globaloncology/Desktop/CPfolder/'
###############################################################################
# STANDARD IMPORTS
###############################################################################
import sys, os, shutil
import subprocess
import matplotlib.pyplot as plt
from scipy import interpolate 
import numpy as np
from matplotlib.pyplot import figure, show, savefig
from matplotlib import cm, colors
import os
import csv
import sys, os, shutil, re, math, time
import math

###############################################################################
# CUSTOM IMPORTS
###############################################################################
import GammaTester
import gamma_test_agent
import OPIMRT_ASCII_reader
import db_utilities
import HPLUSQA_Constants as HPC
import dose_plane_base
import OPIMRT_ASCII_reader
import utilities
import LossRate
###############################################################################
# Function
###############################################################################
def Topasdata(filenames,Energy,Depth):       
    #adjust depth to bin number
    DepthinCm=Depth 
    IDD=0   
    DepthinBin=int(Depth/0.1)    
    Range=250
    i=0
    Array=np.array([[ 0 for m in range(Range) ] for n in range(Range) ])
    with open(filenames) as idd:
        result= csv.reader(idd)
        for line in result:                           
                if i>8 and int(line[2])==DepthinBin:                   
                   Array[int(line[0])][int(line[1])] =float(line[3])*1000000
                   IDD=IDD+float(line[3])*1000000
                i=i+1               
                                         
    idd.closed     
    print IDD
    return Array,IDD
#----------------------------------------------------------------------------------------------
def spotposition(CPfileName,SpotsfileName,ShiftX,ShiftY):
    Result=[]    
    EnergyTable=[]
    i=0
    with open(CPfileName) as idd:
            result= csv.reader(idd)            
            for line in result:
                if i>1:                
                   EnergyTable=EnergyTable+[float(line[2]) for i in range(int(line[1]))]
                i=i+1               
    idd.closed
    i=0   
    n=0    
    with open(SpotsfileName) as idd:
            spotfile= csv.reader(idd)
            for line in spotfile:
                if i>1:
                    #line 1=x,line 2=y,line 3=weight,line 4=MU
                    try:
                        
                         Result.append([int(float(line[1])+ShiftX)*10,int(float(line[2])+ShiftY)*10,EnergyTable[n],float(line[4])])
                         n=n+1  
                    except:
                         Result.append([float(line[1])+ShiftX,float(line[2])+ShiftY,EnergyTable[n-1],float(line[4])])
                         print 'Spot file problem'                                    
                i=i+1
    idd.closed  
  
    return Result
        

def searchTopasFile(Energy):
    EnergyData=[72.5, 73.4, 74.3, 75.2, 76.1, 77.0, 77.9, 78.8, 79.7, 80.5, 81.4, 82.2, 83.1, 83.9, 84.7, 85.5, 85.6, 86.4, 87.2, 88.0, 88.8, 89.6, 90.4, 91.1, 91.9, 92.7, 93.5, 94.2, 95.0, 95.7, 96.9, 98.0, 99.5, 100.9, 102.4, 103.8, 105.2, 106.6, 108.0, 109.4, 110.7, 112.1, 113.4, 114.7, 116.0, 117.3, 118.6, 119.9, 121.2, 122.5, 124.0, 125.6, 127.4, 129.2, 131.0, 132.8, 134.6, 136.4, 138.1, 139.8, 141.6, 143.2, 144.9, 146.9, 148.8, 151.0, 153.2, 155.3, 157.4, 159.5, 161.6, 163.9, 166.2, 168.8, 171.3, 173.7, 176.2, 178.6, 181.1, 183.4, 185.8, 188.2, 190.5, 193.0, 195.6, 198.3, 201.0, 203.7, 206.3, 209.0, 211.0, 211.6, 214.2, 216.7, 219.3, 221.8]
            
    TopasPath='G:/7.11 double/'
    if float(Energy) in EnergyData:
       position=EnergyData.index(float(Energy))
    else:
        print 'not in list'
        print Energy
        position=0
    if position>9:
        S1filename=TopasPath+'96_S1_Run_00'+str(position)+'.csv'
    else:
        S1filename=TopasPath+'96_S1_Run_000'+str(position)+'.csv'        
    if position>9:
        S2filename=TopasPath+'96_S2_Run_00'+str(position)+'.csv'
    else:
        S2filename=TopasPath+'96_S2_Run_000'+str(position)+'.csv'
    return S1filename,S2filename
def SingleDose(Energy,DepthinCm,Snot):    
    Weight=0.2  
    if Energy<80:
        Weight=0.12
    if 80<=Energy<=150:
        Weight=0.08
    if  Energy>150:
        Weight=0.00  
    [S1,S2]=searchTopasFile(Energy)
    #Path='G:/7.11 double/'
    #S1=Path+'S1/'+str(Energy)+'.csv'
    #S2=Path+'S2/'+str(Energy)+'.csv'
    S1Dose,MontroIDDS1,=Topasdata(S1,Energy,DepthinCm)   
    
    if Energy<=150:
        #print 'double source'
        S2Dose,MontroIDDS2=Topasdata(S2,Energy,DepthinCm)       
        MontroIDD=MontroIDDS1*(1-Weight)+MontroIDDS2*Weight
        MeasurementIDD=SingleIDDDose(Energy,DepthinCm,Snot)
        Factor=(MeasurementIDD/MontroIDD)        
        Dose=(S1Dose*(1-Weight)+S2Dose*Weight)*Factor
        print Factor
                  
    else:
        MontroIDD=MontroIDDS1
        MeasurementIDD=SingleIDDDose(Energy,DepthinCm,Snot)
        try:
            Factor=(MeasurementIDD/MontroIDD)
        except:
            Factor=0
            print 'no Dose in plane'
        Dose=S1Dose*Factor      
    return Dose*1000
    
def SumDose(EnergyList,Depth,Snot):
    Range=1000
    i=0
    j=0
    SumDose=np.array([[ 0 for i in range(Range) ] for j in range(Range) ])
    A=[]
    #Thr=[[ 0 for i in range(Range) ] for j in range(Range) ]
    LastEnergy=EnergyList[0][2]    
    LastSingleDose=SingleDose(EnergyList[0][2],Depth,Snot) 
    i=0
    for SpotPosition in EnergyList:        
        if SpotPosition[2]!=LastEnergy:
            LastSingleDose=SingleDose(SpotPosition[2],Depth,Snot) 
        LastEnergy=SpotPosition[2] 
        #print LastEnergy 
        x=SpotPosition[0] 
        y=SpotPosition[1] 
        print np.min(LastSingleDose),'minum dose'
        SumDose[375+x:625+x,375+y:625+y]=SumDose[375+x:625+x,375+y:625+y]+LastSingleDose
        print np.min(SumDose),str(i),'negative numbner',np.max(SumDose)
        i=i+1   
    A=SumDose[375:625,375:625] *0.0001*0.3  
    SumDose=np.transpose(A)
    
    return SumDose
    
def SearhFolder(rootDir):
      list_dirs=os.walk(rootDir)
      filenames=[]
      for root,dirs,files in list_dirs:
          for f in files:
              filenames.append (os.path.join(rootDir,f))
      return filenames
def SingleIDDDose(Energy,Depth,Snot):
     path='F:/Monte Carlo/IDDMeasu/'
     NameList=SearhFolder(path)
     b=[]
     EnergyIDDList=[]
     position=0
     IDD=[]
     IDDList=[]
     X=[]#IDD depth
     Y=[]#IDD dose
     Energy=round(float(Energy),1)  
     
     for item in NameList:
          b.append(item[len(path):-4])
    
     for item in b:
           EnergyIDDList.append(round(float(item),1))
     if float(Energy) in EnergyIDDList:
         position=EnergyIDDList.index(Energy)         
     else:
         print ('no IDD file')
         print Energy        
     #print sorted(EnergyIDDList)
     file=open(NameList[position])
     for item in file:
          IDD.append([item[6:-3]])
     for item in IDD:
         b=item[0].find(',')
         if b>2:
               c=float(item[0][0:b-1])/10   #change to cm
               d=item[0][(b+3):-1]
               if 'E' not in d:                   
                   X.append(c)
                   Y.append(float(d))  
     dev=[]
     for item in X:
          dev.append(abs(item-Depth))
     if  min(dev)<0.1:                       #cm
         position=dev.index(min(dev))         
         Result=Y[position]
         #print 'original depth dose'
     else:                  
         position=dev.index(min(dev))
         if len(Y)-5>position>5:            
            Fitdose=Y[position-6:position]
            FitDepth=X[position-6:position]         
            Fitline=np.polyfit(FitDepth,Fitdose,3)
            FitFunction=np.poly1d(Fitline)           
            Result=FitFunction(Depth)
            if Result<0.1*sum(Fitdose)/len(Fitdose) or Result>2*sum(Fitdose)/len(Fitdose):
                    Result=(Y[position]+Y[position-1])/2        
            
         if position<=5:            
            Fitdose=Y[0:6]
            FitDepth=X[0:6]         
            Fitline=np.polyfit(FitDepth,Fitdose,3)
            FitFunction=np.poly1d(Fitline)
            Result=FitFunction(Depth)
            if Result<0.1*sum(Fitdose)/len(Fitdose)or Result>2*sum(Fitdose)/len(Fitdose):
                    Result=(Y[position]+Y[position-1])/2       
            
         if position>=len(Y)-5:          
            Result=0
     #print Result
     if Snot>23 and Energy>150:
           #SLossRate=1-LossRate.SinLossRate(Energy,Snot,Depth)
           #return Result*SLossRate
           return Result*0.94
     else:
           return Result
     
         
#---------------------------------------------------------------------------------------  
def OpgReader(filename,Factor):
    a = OPIMRT_ASCII_reader.OPIMRT_ASCII_Reader(filename,Factor,True)
    b=a.get_dose_array(filename)   
    CaculateUse=[[ 0 for i in range(len(b)) ] for j in range(len(b))]
    
    m=0
    n=0
    
    while m<len(b):
        n=0
        while n<len(b):
             CaculateUse[m][n]=b[-m-1][n]
             n=n+1
        m=m+1   
        
    return CaculateUse
    
def OpgfileWriter(DoseArray,FileName,ImageName):
    z=[]
    Header=['<opimrtascii>','\n'
    '<asciiheader>','\n'
    'File Version:       3','\n'
    'Separator:          ","','\n'
    'Workspace Name:','/n'     
    'File Name:          C:\Documents and Settings\user\Desktop\G3\Patient','\n'
    'Image Name:         1178377_ARL_20.0','\n'
    'Radiation Type:     Proton','\n'
    'Energy:             0.0 MV','\n'
    'SSD:                2430.0 mm','\n'
    'SID:                1000.0 mm','\n'
    'Field Size Cr:      250.0 mm','\n'
    'Field Size In:      250.0 mm','\n'
    'Data Type:          Abs. Dose','\n'
    'Data Factor:        1','\n'
    'Data Unit:          Gy','\n'
    'Length Unit:        mm','\n'
    'Plane:              XY','\n'
    'No. of Columns:     250','\n'
    'No. of Rows:        250','\n'
    'Number of Bodies:   1','\n'
    'Operators Note:     Corr.: Unif.Calib.,Bkgnd,Calib Output,Temp.,Press.,','\n'
    '</asciiheader>','\n'
    '<asciibody>','\n',
    'Plane Position:     0.0 mm','\n']
    #Header[6]='Image Name:      '+ImageName +','+'\n'
    Endbody=['</asciibody>','\n'
    '</opimrtascii>','\n']
    Xindex=np.linspace(-125,124,len(DoseArray))
    for item in Xindex:
        z.append(','+str(item))
    Xtitle='X[mm]'
    Ytitle='Y[mm]\n'
    Yindex=np.linspace(-125,124,len(DoseArray))
    DoseValue=DoseArray
    j=str()
    i=[]
    for item in DoseValue:
        for Digit in item:
            j=j+str(Digit)+','
        i.append(j)
        j=str()

    n=0
    opgfile=FileName
    with open(opgfile,'w+')as f:
         for item in Header:
             f.write('%s' % item)    
         f.write('%s'% Xtitle)
         for item in z:
             f.write('%s' % item)
         f.write('\n')    
         f.write('%s' % Ytitle)
         for item in Yindex:
             m=str(item)+','+i[n]
             f.write('%s\n' % m)
             n=n+1
         for item in Endbody:
             f.write('%s' % item)
         f.close()
    text=open(opgfile).read()
    text=text.replace('\n','\r\n')
    open(opgfile,'wb').write(text)
#-------------------------------------------------------------------------------------------
def Interpolate(image,Rate):
    Range=len(image)
    Rate=float(Rate)
    x=[]
    y=[]
    xnew=[]
    ynew=[]
    x=np.arange(0,Range,1)
    y=np.arange(0,Range,1)    

    f=interpolate.interp2d(x, y, image, kind='cubic') 
    xnew=np.arange(0,Range,1/Rate)
    ynew=np.arange(0,Range,1/Rate)
    Result=f(xnew,ynew)
    return Result
def Normalization(inputArray):
    Result=[[ 0 for i in range(len(inputArray)) ] for j in range(len(inputArray)) ]
    minNum=np.amin(inputArray)
    maxNum=np.amax(inputArray)
    x=0
    y=0
    for row in inputArray:
        y=0
        for num in row:
             
             num=float(num)
             Result[x][y]=(num-minNum)/(maxNum-minNum)
             y=y+1
        x=x+1
    return Result
def MatrixDataTransfer(Matrixx):
    Measurment=Interpolate(Matrixx,10)    
    Range=250
    BinRate=10 #1cm=10bin 
    MatrixxImage=[[ 0 for i in range(Range) ] for j in range(Range) ]
    LeftCenter=(len(Measurment)/2)-10
    m=0
    n=0
    X=0
    Y=0   
    thresh=np.amax(Measurment)*0.01    
    while m<len(Measurment):
        n=0
        while n<len(Measurment):                       
                  X=int(round((m-LeftCenter)*0.762-0.381)+Range/2)
                  Y=int(round((n-LeftCenter)*0.762-0.381)+Range/2)  
                  if X>249:
                      X=249
                  if Y>249:
                     Y=249                          
                  try:             
                     if Measurment[m][n]>thresh :         
                          MatrixxImage[X][Y]=Measurment[m][n]              
                  except:
                      print 'out of range',X,Y,m,n
                      MatrixxImage[X][Y]=Measurment[m-1][n-1]  
                  n=n+1
        m=m+1
    
    return MatrixxImage

def Allplot(Topas,Measurment,name):
   # Topas=Normalization(Topas)
    #Measurment=MatrixDataTransfer(Measurment)
    #Measurment=Normalization(Measurment)
    #Toaps=Interpolate(Topas,10)
    
    left=0
    right=250
    fig=plt.figure(figsize=[30,30])
    plt.subplot(2,2,1,aspect='equal')
    plt.grid(True)
    plt.imshow(Measurment)
    plt.xlim(left,right)
    plt.ylim(left,right)
    plt.colorbar()
    plt.subplot(2,2,2,aspect='equal')                       
    plt.imshow(Topas)
    plt.xlim(left,right)
    plt.ylim(left,right)
    cbar = plt.colorbar()
    plt.grid(True)  
    #cbar.set_ticks(np.linspace(0,0.9,6))
    plt.subplot(2,2,3,aspect='equal')
    plt.contour(Measurment)  
    plt.xlim(left,right)
    plt.ylim(left,right)
    cbar = plt.colorbar()
    plt.grid(True)  
    #cbar.set_ticks(np.linspace(0,0.9,6))  
    plt.subplot(2,2,4,aspect='equal')  
    #plt.contour(Topas)
    plt.contour(Measurment,linestyles='dashed')  
    plt.xlim(left,right)
    plt.ylim(left,right)    
    cbar = plt.colorbar() 
    plt.grid(True)
    plt.show()
    plt.savefig(name,dpi=128)
    plt.close('all')
#--------------------------------------------------------------------------------------------
def CsvName(rootDir):
      list_dirs=os.walk(rootDir)
      #print list_dirs
      filenames=[]
      for root,dirs,files in list_dirs:
          for f in files:
              filenames.append (os.path.join(rootDir,root,f))
      return filenames
def SearchFile(FilePath):
   Filenames=CsvName(FilePath)
   CP=[]
   Spot=[]
   field=[]
   OPG=[]
   for item in Filenames:
       if 'opg' in item:
           OPG.append(item)
       if 'Spots' in item:          
          Spot.append(item)
          a=item.find('-')
          b=item.find('-',a+1)
          #field.append(item[a+1:b])
       if 'CPs' in item:
           CP.append(item)
   #print Spot,CP,OPG
   Result=[[],[],[],[],[],[]]
   for item in OPG:                
          a=item.find('.opg')
          b=item.find('_',a-10)
          c=item.find('_',b+1) 
          d=item.find('/',b-10)         
          field.append(item[b+1:c])
          depth=item[c+1:a]
          if depth=='265':
                Result[4].append(5)
          else:
                Result[4].append(float(depth))
          Result[3].append(item)
          Result[5].append(item[d+1:b]) 
          
   n=0
   while n<len(field):
         for item in Spot:
             if field[n] in item:
                 Result[0].append(item)
                 Result[2].append(field[n])
         for item in CP:
             if field[n] in item:
                 Result[1].append(item)
         n=n+1
   return Result
def MatchCenter(Cal,Meas,Calname):
    from scipy import ndimage
    CenterMX,CenterMY=ndimage.measurements.center_of_mass(np.array(Meas))
    CenterCX,CenterCY=ndimage.measurements.center_of_mass(Cal)  
    Xshift=int(CenterMX-CenterCX)
    Yshift=int(CenterMY-CenterCY)
    print Xshift,Yshift
    plt.figure
    CaculateUse=[[ 0 for i in range(len(Cal)) ] for j in range(len(Cal))]
    m=0
    n=0
    try:
        while m<len(Cal):
              n=0
              while n<len(Cal):
                     if m+Yshift<len(Cal) and n+Xshift<len(Cal): 
                        CaculateUse[m+Yshift][n+Xshift]=Cal[m][n]
                     n=n+1
              m=m+1
    except:
             CaculateUse=Cal
    plt.imshow(Meas)
    plt.scatter(CenterMY,CenterMX,marker='x')
    plt.scatter((CenterCY+Yshift),(CenterCX+Xshift),color='r')
    plt.scatter(CenterCY,CenterCX,color='g')    
    plt.savefig(Calname[0:len(Calname)-4]+'.png',dpi=128)
    return CaculateUse
###############################################################################
# CONSTANTS
###############################################################################
DICOM_FILE = "y:\\DICOM\\Rodriguez_ALA_3.9_TPS.dcm"
COUNTER_VAR = 0
OUTPUT_DIR = "F:/"
INDEX_PAGE = '/home/hplusqa/public_html/mar07/index.php'
SPLIT_PAGE = '/home/hplusqa/public_html/mar07/split.php'

###############################################################################
# STANDARD IMPORTS
###############################################################################
import sys, os, re, math, time
import numpy as np
import multiprocessing
import pylab
import uuid

from matplotlib import rc
rc('text', usetex=False)



###############################################################################
# CUSTOM IMPORTS
###############################################################################
import OPIMRT_ASCII_reader

import utilities
import HPLUSQA_Constants as HPC



###############################################################################
# FUNCTIONS
###############################################################################
def gamma_for_position(gamma_normalization_dose, gamma_search_length, gamma_dose_requirement, gamma_percentage,\
    meas, block, tps, x_test_points, y_test_points, grid_x, grid_y, distance_arr, results_list, local_normalization):

    dose_normalization = 1.0/(gamma_normalization_dose*gamma_percentage) #Be careful "percentage is fraction, not precentage
    records = []

#    measured_dose_array = meas.interpolate(x_test_points,y_test_points)

    for i in range(len(x_test_points)):
        for j in range(len(y_test_points)):

            #check to see if point is inside block for this thread,
            # if not skip it
            counter = i*j + j
            if counter < block[0] or counter > block[1]:
                continue

            measured_dose = meas.interpolate(x_test_points[i],y_test_points[j])
            tps_dose      = tps.interpolate(x_test_points[i],y_test_points[j])
            delta = 1 - tps_dose/(measured_dose + 1E-10)
#            measured_dose = measured_dose_array(x_test_points[i],y_test_points[j])

            #skip point of low dose, i.e. dose < GAMMA_DOSE_REQUIREMENT
            if measured_dose/gamma_normalization_dose < gamma_dose_requirement:
                continue

            if local_normalization:
                dose_normalization = 1.0/(measured_dose*gamma_percentage)

            x_values = x_test_points[i] + grid_x - gamma_search_length
            y_values = y_test_points[j] + grid_y - gamma_search_length

            #calculate dose difference grid
            dose_matrix = tps.interpolate(x_values,y_values)
            dose_matrix -= measured_dose
            dose_matrix *= dose_normalization
            dose_matrix = dose_matrix**2

            #Add distance and dose in quadrature
            gamma_matrix = dose_matrix + distance_arr

            index = dose_matrix.argmin()
            n, m = dose_matrix.shape
            dose_x, dose_y = index/m, index%n
#            print "%s, %s, %   s, %s, %s" % (index, n, m, dose_x, dose_y)

            gamma_min = math.sqrt(gamma_matrix.min())
            index = gamma_matrix.argmin()
            n, m = gamma_matrix.shape
            index_x, index_y = index/m, index%n

            rec = (i, j, gamma_min, index_x, index_y, dose_x, dose_y, delta)
            records .append(rec)

    #only get the lock one time
    results_list += records
    return results_list



###############################################################################
# CLASS DEFINITION
###############################################################################

class GammaTester:
    normalization_dose=1
    RESULTS = []
    COUNTER = 1

    def get_distance_matrix(self, m, d):
        '''
        returns a x_points by y_points array of distances to the center.

        This array pre-asigns a distance to every point in the gamma index
        search grid and allows for the gamma to calculated for an entire matrix
        of points.
        '''
        n = m
        x_arr = np.array([np.linspace(-d,d, m)**2,]*n)
        y_arr = x_arr.transpose()

        d_arr = (y_arr + x_arr)

        return d_arr



    def get_ranges_for_threads(self, N, T):

        #figure out how many iterations for each thread
        block_sizes = [N/T]*T
        for i in range(N%T):
            block_sizes[i] += 1

        ranges = []
        index = 0
        for iBlock in range(T):
            block =[index, index+block_sizes[iBlock] - 1]
            index += block_sizes[iBlock]
            ranges.append(block)
#            print ranges[-1], (ranges[-1][1] - ranges[-1][0])
        return ranges


    def process_results(self, results_list, test_point_size, distance_arr ):
        #min_arr keeps track of where the minimum is located
        gamma_array = np.zeros((test_point_size,test_point_size))
        gamma_values = []
        number_passing = 0
        min_arr =  distance_arr*0.0
        dose_diff_arr = distance_arr*0.0

        average_diff = 0.0

        for i, result in enumerate(results_list):

            gamma_i, gamma_j, gamma, min_i, min_j, dose_i, dose_j, delta = result
            average_diff += delta
#            print "%d, %d, %d, %.3f" % (i, gamma_i, gamma_j, gamma)
            gamma_values.append(min(gamma,2.0))
            gamma_array[gamma_i, gamma_j] = gamma
#            gamma_array[gamma_i, gamma_j] =min(gamma,1.0)

            #don't track
            if distance_arr[min_i, min_j] > 0.01:
                min_arr[min_i, min_j] += 1

            m = dose_diff_arr.shape[0]
            if dose_i != 0 and dose_i !=m -1 and dose_j !=0 and dose_j !=m -1:
                dose_diff_arr[dose_i, dose_j] += 1
            if gamma < 1:
                number_passing += 1

        gamma_array = np.fliplr(gamma_array)
        gamma_array = np.flipud(gamma_array)

        return gamma_array, gamma_values, number_passing, min_arr, dose_diff_arr, average_diff/float(len(results_list))




    def run_gamma(self, meas, tps, gamma_normalization_dose, number_threads, local_normalization=0):

        print "\n\n--------------\nrun_gamma: extent %s\n-------------\n\n" % meas.get_extent()
        x_test_points = np.linspace(meas.get_extent()[1], meas.get_extent()[0], self.gamma_test_points )
        y_test_points = np.linspace(meas.get_extent()[3], meas.get_extent()[2], self.gamma_test_points )


        # the "+ 1 - (G3.GAMMA_SEARCH_GRID_SIZE%2)" part ensures that there are an odd number of points in the grid
        # which in turn ensures that there is a test point at distance = 0.
        n_search_grid = self.gamma_search_grid_size + 1 - (self.gamma_search_grid_size%2)

        grid_x = np.linspace(0,2.0*self.gamma_search_length,n_search_grid)
        grid_y = np.linspace(0,2.0*self.gamma_search_length,n_search_grid)

        #precalculate the distance terms because they are constant for each point in grid
        distance_arr = self.get_distance_matrix(n_search_grid,self.gamma_search_length)/(self.gamma_distance*self.gamma_distance + 1E-12)


        print "Order: %sx%sx%sx%s = %s" % \
            (len(grid_x),len(grid_y),len(x_test_points),len(y_test_points),
             len(grid_x)*len(grid_y)*len(x_test_points)*len(y_test_points))

        #interpolation is lazy, so call interpolate to get it set up before running threads
        #meas.interpolate(x_test_points[0],y_test_points[0])

        mgr = multiprocessing.Manager()
        results_list = mgr.list()
        procs = []
        blocks = self.get_ranges_for_threads(len(x_test_points)*len(y_test_points), number_threads)
        for block in blocks:

            proc = multiprocessing.Process(target=gamma_for_position, \
                            args=(gamma_normalization_dose, self.gamma_search_length, \
                            self.gamma_dose_requirement, self.gamma_percentage, \
                            meas, block, tps, x_test_points, y_test_points, \
                            grid_x, grid_y, distance_arr, results_list, local_normalization))

            procs.append(proc)
            proc.start()
            #print "started proc for block %s" % block

        all_done = False
        while not all_done:
            all_done = True
            for proc in procs:
                if proc.exitcode == None:
                    all_done = False
#                    print "waiting for processess to finish"
                    time.sleep(0.01)

        for proc in procs:
            proc.join()
            proc.terminate()
#            print "proc %s, exit codes %s" % (proc.name, proc.exitcode)


        results = []
        for r in results_list: results.append(r)

        return self.process_results(results, self.gamma_test_points, distance_arr)



    def get_difference(self, meas, tps, gamma_normalization_dose):

        ### THIS IS WIP AND DOES NOT YET WORK

        x_test_points = np.linspace(meas.get_extent()[1], meas.get_extent()[0], self.gamma_test_points )
        y_test_points = np.linspace(meas.get_extent()[3], meas.get_extent()[2], self.gamma_test_points )
        #tps_dose_matrix = tps.interpolate(x_test_points,y_test_points)
        #meas_dose_matrix = meas.interpolate(x_test_points,y_test_points)

        diff_matrix = meas - tps

        return pylab.imshow(diff_matrix, cmap=cm.spectral,  vmin=vmin, vmax=vmax, extent=meas.get_image_extent())


    def draw_contours(self, measurement, calc, normalization):

        meas_cont = measurement.get_contour(normalization,'solid', self.contour_lines)
        calc_cont = calc.get_contour(normalization,'dotted', self.contour_lines)
        #format_str = "%s%s" % (2,'%%')
        format_str = '%.0f%%'
        cb = pylab.colorbar(meas_cont, format=format_str)
        cb.ax.set_ylabel("100%% is %.2f Gy" % normalization, fontsize=20)
        pylab.grid(True)

        pylab.xlabel("x [cm]", fontsize=20)
        pylab.ylabel("y [cm]", fontsize=20)



    def strip_tag(self, tag):
        return tag.replace(",","").replace("(","").replace(")","").replace(" ","_")


    def setup_output_folder(self, tag):

        name = tag.split("_")
        name = "%s_%s" % (name[0], name[1])
        beam = tag.split("_")
        try:
            beam = beam[2] + "_" + beam[3]
        except IndexError:
            beam = "no_beam_name"

        output_folder = "%s/%s" % (OUTPUT_DIR, name)
        command = "mkdir %s " % output_folder
        print command
        os.popen(command)

        output_folder = "%s/%s" % (output_folder, beam)
        command = "mkdir %s" % output_folder
        print command
        os.popen(command)

        command = "cp %s %s %s" % (INDEX_PAGE, SPLIT_PAGE, output_folder)
        os.popen(command)

        return output_folder



    def generate_comparison_plots(self, measurement, tps,Calname,MeasName ):
        tag='test point'
        
        output_folder = 'F:/'
        pylab.close('all')
        pylab.clf()
        meas_plot = pylab.subplot(1,1,1, aspect='equal')





        #------------------------------- RUN GAMMA ---------------------------------------------
        normalization_dose=1
        local_normalization=False
        g_array, g_values, number_passing, min_arr, dose_diff_arr, average_diff = self.run_gamma(measurement, tps, normalization_dose
, self.number_threads, local_normalization)

        #------------------------------- MAKE GAMMA PLOT ----------------------------------------
        extent=[-self.gamma_search_length, self.gamma_search_length, -self.gamma_search_length, self.gamma_search_length]
        print "\n\nGAMMA SEARCH LENGTH, %s" % self.gamma_search_length
        pylab.clf()
        pylab.subplot(1,1,1, aspect='equal')
        min_image = pylab.imshow(min_arr, extent=extent, cmap=self.gamma_color_scheme)
        pylab.colorbar(min_image)
        pylab.grid(True)
        hist_title = "Gamma Position: %s" % (tag)
        pylab.suptitle(hist_title, fontsize=22)
        pylab.xlabel("search x [cm]", fontsize=20)
        pylab.ylabel("search y [cm]", fontsize=20)
        filename =Calname+'_Gamma.png'    
        pylab.savefig(filename)
        pylab.clf()


        #------------------------------- MAKE GAMMA HISTOGRAM ----------------------------------------
        pylab.close('all')
        pylab.subplot(1,1,1, aspect='equal')
        pylab.clf()
        hist = pylab.hist(g_values,self.gamma_hist_bins,range=[0,2])
        #hist_title = "Gamma Index: %s" % (tag)
        hist_title = "Gamma values for test positions"
        pylab.suptitle(hist_title, fontsize=22)
        # draw a default vline at x=1 that spans the yrange

        gamma_data_dict = { 'DTA': self.gamma_distance*10, \
                            'tolerance':self.gamma_percentage*100,
                            'normalization': normalization_dose*100, \
                            'dose_min': normalization_dose*self.gamma_dose_requirement*100, \
                            'test_points': len(g_values), \
                            'percent_passing' : 100*number_passing/float(len(g_values))  }

        line_1 = pylab.axvline(x=1, linewidth=2, color='r')
        text_line_1 = "DTA: %.1f mm" % (self.gamma_distance*10)
        text_line_2 = "Dose tolerance: %2.0f%%" %  (self.gamma_percentage*100)
        text_line_3 = "Normalization: %.1f cGy" % (normalization_dose*100)
        text_line_4 = "Min dose: %.1f cGy" % (normalization_dose*self.gamma_dose_requirement*100)
        text_line_5 = "Test points: %d" % (len(g_values))
        text_line_6 = "Points passed: %d" % (number_passing)
        text_line_7 = "Pass rate: %.1f%%" % (100*number_passing/float(len(g_values)) )
        text_line_8 = "Avg. Diff.: %.2f%%" % (100*average_diff)

        pylab.text(0.86, 1.12, text_line_1, transform = meas_plot.transAxes, fontsize=16)
        pylab.text(0.86, 1.05, text_line_2, transform = meas_plot.transAxes, fontsize=16)
        pylab.text(0.86, 0.98, text_line_3, transform = meas_plot.transAxes, fontsize=16)
        pylab.text(0.86, 0.91, text_line_4, transform = meas_plot.transAxes, fontsize=16)
        pylab.text(0.86, 0.84, text_line_5, transform = meas_plot.transAxes, fontsize=16)
        pylab.text(0.86, 0.77, text_line_6, transform = meas_plot.transAxes, fontsize=16)
#        pylab.text(0.86, 0.70, text_line_8, transform = meas_plot.transAxes, fontsize=16)

        pylab.text(0.86, 0.58, text_line_7, transform = meas_plot.transAxes, color='red', fontsize=19)
        pylab.xlabel("gamma value", fontsize=20)
        pylab.ylabel("counts", fontsize=20)
        filename =Calname+'_hist.png'        
        pylab.savefig(filename)
        pylab.clf()

        #------------------------------- MAKE GAMMA MAP ----------------------------------------
        g_image = pylab.imshow(g_array, extent=measurement.get_image_extent(), cmap=self.gamma_color_scheme, vmin=0, vmax=1.25)
        #title = "Gamma(%.0f mm, %.0f%%), %s" % (self.gamma_distance*10, self.gamma_percentage*100, tag)
        title = "Gamma(%.0f mm, %.0f%%)" % (self.gamma_distance*10, self.gamma_percentage*100)
        pylab.suptitle(title, fontsize=22)
        pylab.grid(True)
        pylab.xlabel("x [cm]", fontsize=20)
        pylab.ylabel("y [cm]", fontsize=20)
        cb = pylab.colorbar(g_image)
        cb.set_ticks(ticks = [0.0,0.25, 0.5, 0.75, 1.0, 1.25])
        cb.ax.set_ylabel('Gamma score', fontsize=22)
        filename = Calname+'_map.png'       
        pylab.savefig(filename)
        pylab.clf()

#        #------------------------------- DIFFERENCE PLOT ---------------------------------------------
        #diff_image = self.get_difference(measurement, tps, normalization_dose)
        #title = "Measured dose - calculated dose"
        #pylab.suptitle(title, fontsize=22)
        #pylab.grid(True)
        #pylab.xlabel("x [cm]", fontsize=20)
        #pylab.ylabel("y [cm]", fontsize=20)
        #filename = '/home/globaloncology/Desktop/Case/42.png'
        #filenames.append(filename)
        #pylab.savefig(filename)
        #pylab.clf()
        #diff_image = pylab.imshow(dose_diff_arr, extent=extent, cmap=self.gamma_color_scheme)
        #pylab.colorbar(diff_image)
        #pylab.grid(True)
        #pylab.clf()

        #return filenames, gamma_data_dict


    def __init__(self, gamma_normalization_method="measured_max"):
#        print "initializing %s" % self
        self.gamma_normalization_method = gamma_normalization_method

        self.gamma_search_length = 0.5 #"cm"
        self.gamma_dose_requirement = 0.1
        self.gamma_distance = 0.3
        self.gamma_percentage = 0.03
        self.gamma_color_scheme = HPC.GAMMA_COLOR_SCHEME
        self.gamma_search_grid_size = HPC.GAMMA_SEARCH_GRID_SIZE
        self.gamma_hist_bins = HPC.GAMMA_HIST_BINS
        self.gamma_test_points = HPC.GAMMA_TEST_POINTS
        self.contour_lines = HPC.CONTOUR_LINES
        self.number_threads = HPC.NUMBER_THREADS


###############################################################################
# MAIN
###############################################################################
def run_gamma_comp(Calname,MeasName):    
    meas_plane = OPIMRT_ASCII_reader.OPIMRT_ASCII_Reader(MeasName, 1, True)   
    tps_plane = OPIMRT_ASCII_reader.OPIMRT_ASCII_Reader(Calname,1, True)
    b=tps_plane.get_dose_array(Calname) 
    c=meas_plane.get_dose_array(MeasName)
    pylab.imshow(b)    
    pylab.savefig(Calname+'.png')
    pylab.imshow(c)
    pylab.show()
    pylab.savefig(MeasName+'.png')
    #meas_plane.auto_crop()
    #meas_plane.shift_xy(5,5)
    #tps_plane.interpolate_myself( meas_plane.get_extent(), HPC.INTERPOLATION_DISTANCE, shift_x=0.0, shift_y=0.0)
    pylab.subplot(2,1,1, aspect='equal')
    _img1 = meas_plane.get_image()
    pylab.colorbar(_img1)
##
    pylab.subplot(2,1,2, aspect='equal')
    img2 = tps_plane.get_image()
    pylab.colorbar(img2)

    print "running GammaTester . . ."
    mygt = GammaTester()
    
    print "generating comparison plots . . ."
    mygt.generate_comparison_plots(meas_plane, tps_plane,Calname,MeasName)

    print "completed run_gamma_comp() . . ."

if __name__ == '__main__':
   import cProfile    
   s_time = time.time()
   Depth=2+6.7
   Snot=0
   FilePath='F:\Case\Armstrong_1087226/'
   names=SearchFile(FilePath)
   #print names
   #0=spot,1=cp,2=field,3=opg,4=depth,5=id   
   n=0
   while n<len(names):
   #while  n<1:
       if names[2][n][0]=='A':          
          Snot=25
          ShiftX=0
          ShiftY=0
       if names[2][n][0]=='B':          
          Snot=30
          ShiftX=0
          ShiftY=0
       if names[2][n][0]=='C':          
          Snot=25
          ShiftX=0
          ShiftY=0
       if names[2][n][0]=='D':          
          Snot=29
          ShiftX=0-5
          ShiftY=0-5
       if names[2][n][0]=='E':          
          Snot=29
          ShiftX=0-5
          ShiftY=0-5
       if names[2][n][0]=='F':          
          Snot=28
          ShiftX=0
          ShiftY=0       
       #print Snot
       #print 'first Shift',ShiftX,ShiftY
       depth=names[4][n]+6.7 #if absorber     
       
       ImageName=names[2][n]+'_'+names[2][n]+'_'+names[2][n] 
       Meas=OpgReader(names[3][n],1.027)
       #print names[3][n]       
       MeasResultName=names[5][n]+'_'+names[2][n]+'_'+str(depth)+'_'+'meas.opg'
       #print  MeasResultName,'new measu'
       FinalMatrixx=MatrixDataTransfer(Meas)
       OpgfileWriter(FinalMatrixx,MeasResultName,ImageName+'meas') 
       
       Array=SumDose(spotposition(names[1][n],names[0][n],ShiftX,ShiftY),depth,Snot)
       TopasResultName=names[5][n]+'_'+names[2][n]+'_'+str(depth)+'_'+'Topas.opg'
       Calname= names[5][n]+'_'+names[2][n]+'_'+str(depth)+'_'+'Topas.opg'
       print Calname
       Array=MatchCenter(Array,FinalMatrixx,Calname)
       OpgfileWriter(Array,TopasResultName,ImageName)  
       Calname= names[5][n]+'_'+names[2][n]+'_'+str(depth)+'_'+'Topas.opg'
       Cal=OpgReader(names[5][n]+'_'+names[2][n]+'_'+str(depth)+'_'+'Topas.opg',1)
       
       Allplot(Cal,FinalMatrixx,names[5][n]+'_'+names[2][n]+str(depth)+'.jpg')
       n=n+1
#-----------------------------------------------------------------------------------------
       
       cProfile.run('run_gamma_comp(Calname,MeasResultName)','none',sort='cumulative')
   ttime = (time.time() - s_time)/60
   print ttime
 
