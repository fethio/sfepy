import csv
import numpy as nm
from math import sin, exp, pi, pow, sqrt
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
import os

"""

Sfepy kodu çalıştırılır. Daha sonra Paraview atılır, oradan csv olarak kaydedilir. Plot için bu Script Çalıştırılır.
Paraview içinde element ID'ler 0'dan başlar.

Script Çalışma Mantığı: csv Dosya içindeki her değer array'lere atanır.
                        createArray function içinde csv dosyasu okunurken, içindeki text yazıları (Bulk stress, mat id...) atlanır. 
                        Bu sebeple 1. elemanın index'i artık 0'dır. Paraview'de de 1. elemanın ID'si 0'dır.

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
      
      self.active_stress_0 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.active_stress_1 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.active_stress_2 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.active_stress_3 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.active_stress_4 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.active_stress_5 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.active_stress_6 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.active_stress_7 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.active_stress_8 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      
      self.Piola_Kirchhoff_stress_0 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Piola_Kirchhoff_stress_1 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Piola_Kirchhoff_stress_2 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Piola_Kirchhoff_stress_3 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Piola_Kirchhoff_stress_4 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Piola_Kirchhoff_stress_5 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Piola_Kirchhoff_stress_6 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Piola_Kirchhoff_stress_7 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.Piola_Kirchhoff_stress_8 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      
      self.total_stress_0 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.total_stress_1 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.total_stress_2 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.total_stress_3 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.total_stress_4 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.total_stress_5 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.total_stress_6 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.total_stress_7 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
      self.total_stress_8 = nm.zeros(shape=(len(self.ordered_files),len(SelectionID)))
    
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
        for csvFilename in self.ordered_files: # her bir time-step icin doner
          f = open(csvFilename, 'rb')
          f.next() # skip the header. index bu sayede 0'dan baslar.
          readFile = csv.reader ((f), delimiter=",",quoting=csv.QUOTE_NONNUMERIC)
          j = 0
          for index,value in enumerate(readFile): # 1 time-step'deki dosya icindeki tum degerler icin doner.
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
                      self.active_stress_0[i][j]           = value[1]
                      self.active_stress_1[i][j]           = value[2]
                      self.active_stress_2[i][j]           = value[3]
                      self.active_stress_3[i][j]           = value[4]
                      self.active_stress_4[i][j]           = value[5]
                      self.active_stress_5[i][j]           = value[6]
                      self.active_stress_6[i][j]           = value[7]
                      self.active_stress_7[i][j]           = value[8]
                      self.active_stress_8[i][j]           = value[9]
                      
                      self.total_stress_0[i][j]            = value[10]
                      self.total_stress_1[i][j]            = value[11]
                      self.total_stress_2[i][j]            = value[12]
                      self.total_stress_3[i][j]            = value[13]
                      self.total_stress_4[i][j]            = value[14]
                      self.total_stress_5[i][j]            = value[15]
                      self.total_stress_6[i][j]            = value[16]
                      self.total_stress_7[i][j]            = value[17]
                      self.total_stress_8[i][j]            = value[18]                      
                      
                      self.neohook_stress_0[i][j]          = value[19]
                      self.neohook_stress_1[i][j]          = value[20]
                      self.neohook_stress_2[i][j]          = value[21]
                      self.neohook_stress_3[i][j]          = value[22]
                      self.neohook_stress_4[i][j]          = value[23]
                      self.neohook_stress_5[i][j]          = value[24]
                      self.neohook_stress_6[i][j]          = value[25]
                      self.neohook_stress_7[i][j]          = value[26]
                      self.neohook_stress_8[i][j]          = value[27]
                      
                      self.green_strain_0[i][j]            = value[28]
                      self.green_strain_1[i][j]            = value[29]
                      self.green_strain_2[i][j]            = value[30]
                      self.green_strain_3[i][j]            = value[31]
                      self.green_strain_4[i][j]            = value[32]
                      self.green_strain_5[i][j]            = value[33]
                      self.green_strain_6[i][j]            = value[34]
                      self.green_strain_7[i][j]            = value[35]
                      self.green_strain_8[i][j]            = value[36]
                      
                      self.bulk_stress_0[i][j]             = value[37]
                      self.bulk_stress_1[i][j]             = value[38]
                      self.bulk_stress_2[i][j]             = value[39]
                      self.bulk_stress_3[i][j]             = value[40]
                      self.bulk_stress_4[i][j]             = value[41]
                      self.bulk_stress_5[i][j]             = value[42]
                      self.bulk_stress_6[i][j]             = value[43]
                      self.bulk_stress_7[i][j]             = value[44]
                      self.bulk_stress_8[i][j]             = value[45]
                      
                      self.Piola_Kirchhoff_stress_0[i][j]  = value[46]
                      self.Piola_Kirchhoff_stress_1[i][j]  = value[47]
                      self.Piola_Kirchhoff_stress_2[i][j]  = value[48]
                      self.Piola_Kirchhoff_stress_3[i][j]  = value[49]
                      self.Piola_Kirchhoff_stress_4[i][j]  = value[50]
                      self.Piola_Kirchhoff_stress_5[i][j]  = value[51]
                      self.Piola_Kirchhoff_stress_6[i][j]  = value[52]
                      self.Piola_Kirchhoff_stress_7[i][j]  = value[53]
                      self.Piola_Kirchhoff_stress_8[i][j]  = value[54]


                      
                  j = j +1    
          i = i + 1

    def Factivation(self,t):
      t_d = 0.0
      t_a = 0.4
      q = 3.5
      t_p = 0.8
      S = 200
      fa = (1+sin(pi*(t-t_d)/t_a - pi/2))/2
      if t <= t_d:
         f_act = 0.0
        
      if t_d < t <= t_d + t_a:
         f_act = pow(fa,q) 
           
      if t_d + t_a < t <= t_p:
         f_act = 1.0
           
      if t_p < t:
         f_act = exp(-S*(t-t_d-t_p))
      return f_act    
    
    def plotData(self):
        plt.figure(figsize=(12, 10))
        mpl.rcParams['lines.linewidth'] = 1.5
        plt.subplot(321)
        plt.ylim(0,1.2)
        plt.plot([t for t in self.Time] , [self.Factivation(t) for t in self.Time],'r')
        plt.xlabel('Time(s), $t$',size = 9)
        plt.ylabel('Activation, $f_a$',size = 9)
        
#        plt.subplot(321)
#        plt.plot(self.Time,self.total_stress_0,'k')
#        plt.xlabel('Time(s), $t$',size = 9)
#        plt.ylabel('Total Stress(MPa), $\sigma_{11}$',size = 9)
#        #plt.ylim(0,8.2)
        
        plt.subplot(322)
        plt.plot( self.Time,self.green_strain_0,'k')
        plt.xlabel('Time(s), $t$',size = 9)
        plt.ylabel('Deformation (mm/mm), $\epsilon_{11}$',size = 9)
        
        plt.subplot(323)
        plt.plot(self.Time,self.neohook_stress_0,'k')
        plt.xlabel('Time(s), $t$',size = 9)
        plt.ylabel('NeoHook Stress (MPa), $\sigma_{11}$',size = 9)
        #plt.ylim(0,0.62)
        
        plt.subplot(324)
        plt.plot(self.Time,self.active_stress_0,'k')
        plt.xlabel('Time(s), $t$',size = 9)
        plt.ylabel('Active Stress (MPa), $\sigma_{11}$',size = 9)
        #plt.ylim(0,0.9)
        
        plt.subplot(325)
        plt.plot(self.Time,self.Piola_Kirchhoff_stress_0,'k')
        plt.xlabel('Time(s), $t$',size = 9)
        plt.ylabel('$2^{nd}$ Piola Kirchhoff Stress (MPa), $\sigma_{11}$',size = 9)
        #plt.ylim(0,8.2)
        
        plt.subplot(326)
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
        print "   Piolo        active      neohook        bulk         Total        Defor"
        print ' %10.3e ' % self.Piola_Kirchhoff_stress_0[80],
        print ' %10.3e ' % self.active_stress_0[80],  
        print ' %10.3e ' % self.neohook_stress_0[80],  
        print ' %10.3e ' % self.bulk_stress_0[80],
        print ' %10.3e ' % self.total_stress_0[80],
        print ' %10.3e ' % self.green_strain_0[80]
        

##############################################################################
pathPoint ='/home/gkdnz/Desktop/PlotOutput/Point'
pathCell = '/home/gkdnz/Desktop/PlotOutput/Cell'
pointID = [5]
cellID = [739] #739
AnalysisMode = ['Nodal','Cell']

#plotPoint = PlotData(pathPoint,pointID)
#plotPoint.createArray(AnalysisMode[0])
#plotPoint.plotData(Time)

plotCell = PlotData(pathCell,cellID)
plotCell.createArray(AnalysisMode[1])
plotCell.plotData()
plotCell.printValue()

