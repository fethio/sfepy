import csv
import numpy as nm
from math import sin, exp, pi, pow, sqrt
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
import os


"""

. Sfepy kodu calistirilir. Daha sonra vtk dosyalari Paraview atilir,
. Paraview icinden csv olarak kaydedilir. Her bir time-step icin 1 dosya cikar. Paraview icinde element iD'ler 0'dan baslar.

. Daha sonra sonuclari Plot yapmak icin bu Script Calistirilir.



Script Calisma Mantigi: csv Dosya icindeki her deger array'lere atanir.
                        createArray function icinde csv dosyasi okunurken, icindeki text yazilari (Bulk stress, mat id...) atlanir. 
                        Bu sebeple 1. elemanin index'i artik 0'dir. Paraview'de de 1. elemanin ID'si 0'dir.
                        Ek olarak value array'inde 0. index'te matId yazar, bu yuzden value[1]'den baslayarak olusturulmus stress-strain array'lerine atanir.

"""
###############################################################################

class  PlotData:
    def __init__(self,path,SelectionID):
      self.path = path
      os.chdir(path) # change current dictionary
      self.ordered_files = sorted(os.listdir(path), key=lambda x: (int(re.sub('\D','',x)),x)) # Array icindeki elemanlari dosya isimlerine gore(0'dan 20'ye) siraladi.
      self.SelectionID = SelectionID
      self.Time = nm.linspace(0,1,len(self.ordered_files))
      self.node_groups = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.u_x         = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.u_y         = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.u_z         = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Points_x    = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Points_y    = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Points_z    = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      
      self.bulk_stress_0 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.bulk_stress_1 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.bulk_stress_2 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.bulk_stress_3 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.bulk_stress_4 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.bulk_stress_5 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.bulk_stress_6 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.bulk_stress_7 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.bulk_stress_8 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
    
      self.neohook_stress_0 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.neohook_stress_1 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.neohook_stress_2 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.neohook_stress_3 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.neohook_stress_4 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.neohook_stress_5 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.neohook_stress_6 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.neohook_stress_7 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.neohook_stress_8 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      
    
      self.green_strain_0 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.green_strain_1 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.green_strain_2 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.green_strain_3 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.green_strain_4 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.green_strain_5 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.green_strain_6 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.green_strain_7 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.green_strain_8 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))

    def createArray(self,AnalysisMode):      
        i = 0       
        for csvFilename in self.ordered_files: # her bir step icin doner. self.ordered_files'lar tum step dosyalaridir.
          f = open(csvFilename, 'rb')
          f.next() # skip the header. index bu sayede 0'dan baslar.
          readFile = csv.reader ((f), delimiter=",",quoting=csv.QUOTE_NONNUMERIC)
          j = 0
          for index,value in enumerate(readFile): # 1 step'deki dosya icindeki tum degerler icin doner.
          # 7 eleman varsa index 0,1,2,3,4,5,6 doner. Bu yuzden elementID 0'dan vermeye basla.  
           for x in self.SelectionID:
              if x == index:
                  if AnalysisMode == 'Nodal':
                    self.u_x[i][j]         = value[1]
                    self.node_groups[i][j] = value[0]
                    self.u_y[i][j]         = value[2]
                    self.u_z[i][j]         = value[3]
                    self.Points_x[i][j]    = value[4]
                    self.Points_y[i][j]    = value[5]
                    self.Points_z[i][j]    = value[6]
                    
                    
                  if  AnalysisMode == 'Cell':                 
                      
                      self.neohook_stress_0[i][j]          = value[10]
                      self.neohook_stress_1[i][j]          = value[11]
                      self.neohook_stress_2[i][j]          = value[12]
                      self.neohook_stress_3[i][j]          = value[13]
                      self.neohook_stress_4[i][j]          = value[14]
                      self.neohook_stress_5[i][j]          = value[15]
                      self.neohook_stress_6[i][j]          = value[16]
                      self.neohook_stress_7[i][j]          = value[17]
                      self.neohook_stress_8[i][j]          = value[18]
                      
                      self.green_strain_0[i][j]            = value[19]
                      self.green_strain_1[i][j]            = value[20]
                      self.green_strain_2[i][j]            = value[21]
                      self.green_strain_3[i][j]            = value[22]
                      self.green_strain_4[i][j]            = value[23]
                      self.green_strain_5[i][j]            = value[24]
                      self.green_strain_6[i][j]            = value[25]
                      self.green_strain_7[i][j]            = value[26]
                      self.green_strain_8[i][j]            = value[27]
                      
                      self.bulk_stress_0[i][j]             = value[1]
                      self.bulk_stress_1[i][j]             = value[2]
                      self.bulk_stress_2[i][j]             = value[3]
                      self.bulk_stress_3[i][j]             = value[4]
                      self.bulk_stress_4[i][j]             = value[5]
                      self.bulk_stress_5[i][j]             = value[6]
                      self.bulk_stress_6[i][j]             = value[7]
                      self.bulk_stress_7[i][j]             = value[8]
                      self.bulk_stress_8[i][j]             = value[9]
                      #print("bulk_0: ",self.neohook_stress_0)


                      
                  j = j +1    
          i = i + 1
 
    
    def plotData(self):
        plt.figure(figsize=(12, 10))
        mpl.rcParams['lines.linewidth'] = 1.5

        
#        plt.subplot(321)
#        plt.plot(self.Time,self.total_stress_0,'k')
#        plt.xlabel('Time(s), $t$',size = 9)
#        plt.ylabel('Total Stress(MPa), $\sigma_{11}$',size = 9)
#        #plt.ylim(0,8.2)
        
        plt.subplot(221)
        plt.plot( self.Time,self.green_strain_0,'k')
        plt.xlabel('Time(s), $t$',size = 9)
        plt.ylabel('Deformation (mm/mm), $\epsilon_{11}$',size = 9)
        
        plt.subplot(222)
        plt.plot(self.Time,self.neohook_stress_0,'k')
        plt.xlabel('Time(s), $t$',size = 9)
        plt.ylabel('NeoHook Stress (MPa), $\sigma_{11}$',size = 9)
        #plt.ylim(0,0.62)
        

        
        plt.subplot(223)
        plt.plot(self.Time,self.bulk_stress_0,'k')
        plt.xlabel('Time(s), $t$',size = 9)
        plt.ylabel('Bulk Stress (MPa), $\sigma_{11}$',size = 9)
        

        
#        plt.subplot(131)
#        plt.plot( Time,self.u_x)
#        plt.xlabel('Time (s)')
#        plt.ylabel('Displacement (m),$u_x$')
        
        plt.tight_layout()
        plt.show()
    
    def printValue(self):
        print "   neohook        bulk       Defor"
        print ' %10.3e ' % self.neohook_stress_0[80],  # burdaki sayi step'i gosterir.
        print ' %10.3e ' % self.bulk_stress_0[80],
        print ' %10.3e ' % self.green_strain_0[80]
        

##############################################################################
pathPoint ='/home/gkdnz/Desktop/PlotOutput/Point'
pathCell = '/home/gkdnz/Desktop/PlotOutput/Cell'
pointID = [5]
cellID = [0] #739
AnalysisMode = ['Nodal','Cell']

#plotPoint = PlotData(pathPoint,pointID)
#plotPoint.createArray(AnalysisMode[0])
#plotPoint.plotData(Time)

plotCell = PlotData(pathCell,cellID)
plotCell.createArray(AnalysisMode[1])
plotCell.plotData()
plotCell.printValue()

