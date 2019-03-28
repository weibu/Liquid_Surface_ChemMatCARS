from PyQt4 import uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtWebKit import *
from PyQt4 import QtCore, QtGui
from mplwidget import MplWidget
from matplotlib.widgets import MultiCursor
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from scipy import special
from scipy.special import *
from scipy.integrate import quad
import pylab as pl
import numpy as np
import time
import os
import math
import cmath
import glob
import sys
import xr_ref as xr
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
import matplotlib as mpl
import periodictable 
from periodictable import *


(Ui_MainWindow, QMainWindow) = uic.loadUiType('mainwindow.ui')

class MainWindow (QMainWindow):
    """MainWindow inherits QMainWindow"""

    def __init__ (self, parent = None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.tabWidget.setCurrentIndex(0)
        self.directory=os.getcwd()
        self.reffiles=[]
        self.reffitfiles=[]
        self.refedfiles=[]
        self.selectedreffiles_rows=[]
        self.selectedreffitfiles_rows=[]
        self.selectedrefedfiles_rows=[]
        self.refcal=[]
        self.sldcal=[]
        self.halftab='      '
        self.refsavefitindex=0
        self.refsavedataindex=0
        self.rodsavefitindex=0
        self.flusavefitindex=0
        self.gixsavefitindex=0
        self.rodfiles=[]
        self.rodfitfiles=[]
        self.selectedrodfiles_rows=[]
        self.selectedrodfitfiles_rows=[]
        self.flufiles=[]
        self.flufitfiles=[]
        self.selectedflufiles_rows=[]
        self.selectedflufitfiles_rows=[]
        self.gixfiles=[]
        self.gixfitfiles=[]
        self.gixedfiles=[]
        self.selectedgixfiles_rows=[]
        self.selectedgixfitfiles_rows=[]
        self.selectedgixedfiles_rows=[]
        self.eleradius=periodictable.constants.electron_radius*1e10
        self.avoganum=periodictable.constants.avogadro_number
        self.boltzmann=1.38065e-23
        mpl.rc('axes',color_cycle=['b','r','g','c','m','y','k'])
        self.ui.refPW.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.ui.refPW.canvas.ax.set_ylabel('Normalized Reflectivity')
        self.ui.refsldPW.canvas.ax.set_xlabel('Z'+' '+r'$[\AA]$')
        self.ui.refsldPW.canvas.ax.set_ylabel('Electron Density Profile'+' '+r'$[e/\AA^{3}]$')
        self.ui.refqofflabel.setText(u'\u212b'+u'\u207b'+u'\u00b9')
        self.ui.refqreslabel.setText(u'\u212b'+u'\u207b'+u'\u00b9')
        self.ui.rodrhounitlabel.setText('e/'+u'\u212b'+u'\u00b3')
        self.ui.rodalphalabel.setText('angle('+u'\u03B1'+')')
        self.ui.roddthlabel.setText('2'+u'\u03B8')
        self.ui.rodqofflabel.setText(u'\u212b'+u'\u207b'+u'\u00b9')
        self.ui.rodsizelabel.setText(u'\u212b')
        self.ui.rodsizereslabel.setText(u'\u212b')
        self.ui.rodroughlabel.setText(u'\u212b')
        self.ui.rodrholabel.setText(u'\u03c1'+'_sub')
        self.ui.rodPW.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.ui.rodPW.canvas.ax.set_ylabel('Intensity [a.u.]')
        self.ui.flusurunitlabel.setText(u'\u212b'+u'\u207b'+u'\u00b2')
        self.ui.flurhotopunitlabel.setText('e/'+u'\u212b'+u'\u00b3')
        self.ui.flurhobotunitlabel.setText('e/'+u'\u212b'+u'\u00b3')
        self.ui.flurelrhobotunitlabel.setText('e/' + u'\u212b' + u'\u00b3')
        self.ui.fluqcunitlabel.setText(u'\u212b' + u'\u207b' + u'\u00b9')
        self.ui.fluqofflabel.setText(u'\u212b'+u'\u207b'+u'\u00b9')  
        self.ui.flurhotoplabel.setText(u'\u03c1'+'_top')
        self.ui.flurhobotlabel.setText(u'\u03c1'+'_bot')
        self.ui.flubetatoplabel.setText(u'\u03B2'+'_top')
        self.ui.flubetabotlabel.setText(u'\u03B2'+'_bot(inc)')
        self.ui.flubetabot2label.setText(u'\u03B2'+'_bot(flu)')
        self.ui.gixPW.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.ui.gixPW.canvas.ax.set_ylabel('Intensity [a.u.]')
        self.ui.gixsldPW.canvas.ax.set_xlabel('Z'+' '+r'$[\AA]$')
        self.ui.gixsldPW.canvas.ax.set_ylabel('Electron Density Profile'+' '+r'$[e/\AA^{3}]$')
        self.ui.gixqofflabel.setText(u'\u212b'+u'\u207b'+u'\u00b9')
        self.ui.gixqmaxlabel.setText(u'\u212b'+u'\u207b'+u'\u00b9')
        self.ui.gixalphalabel.setText('angle('+u'\u03B1'+')')
        self.ui.gixdthlabel.setText('2'+u'\u03B8')
        self.errorlist=np.array([[1, 1.074], [2, 1.204], [3, 1.222], [4, 1.220], [5, 1.213], [6, 1.205], [7, 1.198], [8, 1.191], [9, 1.184], [10, 1.178], [11, 1.173], [12, 1.168], [13, 1.163], [14, 1.159], [15, 1.155], [16, 1.151], [17, 1.148], [18, 1.145], [19, 1.142], [20, 1.139], [22, 1.134], [24, 1.129], [26, 1.125], [28, 1.121], [30, 1.118], [32, 1.115], [34, 1.112], [36, 1.109], [38, 1.106], [40, 1.104], [42, 1.102], [44, 1.100], [46, 1.098], [48, 1.096], [50, 1.094], [60, 1.087], [70, 1.081], [80, 1.076], [90, 1.072], [100, 1.069], [120, 1.063], [140, 1.059], [160, 1.055], [180, 1.052]]) #, [3000, 1.050]])
        self.initRefPar()
        self.initRodPar()
        self.initFluPar()
        self.initGixPar()
        self.connect(self.ui.action_About, SIGNAL('triggered()'),self.showAbout)
        self.connect(self.ui.action_Open_REF_file, SIGNAL('triggered()'),self.openRefFile)
        self.connect(self.ui.addreffilePB, SIGNAL('clicked()'),self.addRefFile)
        self.connect(self.ui.reffileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedRefFile)
        self.connect(self.ui.reflegendCB,SIGNAL('stateChanged(int)'), self.updateRefPlot)
        self.connect(self.ui.reflegendlocCoB,SIGNAL('currentIndexChanged(int)'), self.updateRefPlot)
        self.connect(self.ui.reflogyCB, SIGNAL('stateChanged(int)'), self.updateRefPlot)
        self.connect(self.ui.rmreffilePB, SIGNAL('clicked()'), self.removeRefFile)
        self.connect(self.ui.addreffitfilePB, SIGNAL('clicked()'),self.addRefFitFile)
        self.connect(self.ui.reffitfileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedRefFitFile)
        self.connect(self.ui.rmreffitfilePB, SIGNAL('clicked()'), self.removeRefFitFile)
        self.connect(self.ui.refscalePB, SIGNAL('clicked()'),self.setRefPlotScale)
        self.connect(self.ui.addrefedfilePB, SIGNAL('clicked()'),self.addRefEDFile)
        self.connect(self.ui.refedfileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedRefEDFile)
        self.connect(self.ui.rmrefedfilePB, SIGNAL('clicked()'), self.removeRefEDFile)
        self.connect(self.ui.refedlegendCB,SIGNAL('stateChanged(int)'), self.updateRefEDPlot)
        self.connect(self.ui.refedlegendlocCoB,SIGNAL('currentIndexChanged(int)'), self.updateRefEDPlot)
        self.connect(self.ui.refedscalePB, SIGNAL('clicked()'), self.setEDPlotScale)
        self.connect(self.ui.insrefslabPB, SIGNAL('clicked()'), self.insRefSlab)
        self.connect(self.ui.rmrefslabPB, SIGNAL('clicked()'), self.rmRefSlab)
        self.connect(self.ui.refnumslabSB, SIGNAL('valueChanged(int)'), self.modRefSlab)
        self.connect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self.updateRefParaVal)
        self.connect(self.ui.refparTW,SIGNAL('cellDoubleClicked(int,int)'),self.setupRefPara)
        self.connect(self.ui.calrefCB, SIGNAL('stateChanged(int)'),self.updateRefCal)
        self.connect(self.ui.calsldCB, SIGNAL('stateChanged(int)'),self.updateRefCal)
        self.connect(self.ui.refsaveCB, SIGNAL('activated(int)'), self.saveRef)
        self.connect(self.ui.refloadCB, SIGNAL('activated(int)'), self.loadRef)
        self.connect(self.ui.fitrefPB,SIGNAL('clicked()'),self.fitRef)
        self.connect(self.ui.refqoffLE,SIGNAL('returnPressed()'),self.updateRefCal)
        self.connect(self.ui.refyscaleLE,SIGNAL('returnPressed()'),self.updateRefCal)
        self.connect(self.ui.refqresLE,SIGNAL('returnPressed()'),self.updateRefCal)
        self.connect(self.ui.refrrfCB, SIGNAL('stateChanged(int)'),self.updateRefCal)
        self.connect(self.ui.clerefconPB, SIGNAL('clicked()'),self.cleRefCon)
        self.connect(self.ui.refsysconPB, SIGNAL('clicked()'),self.updateRefSysPara)
        self.connect(self.ui.errcalPB, SIGNAL('clicked()'), self.errorCal)
        ##rod analysis
        self.connect(self.ui.action_Open_ROD_file, SIGNAL('triggered()'),self.openRodFile)
        self.connect(self.ui.addrodfilePB, SIGNAL('clicked()'),self.addRodFile)
        self.connect(self.ui.rodfileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedRodFile)
        self.connect(self.ui.rodlegendCB,SIGNAL('stateChanged(int)'), self.updateRodPlot)
        self.connect(self.ui.rodlegendlocCoB,SIGNAL('currentIndexChanged(int)'), self.updateRodPlot)
        self.connect(self.ui.rodlogyCB, SIGNAL('stateChanged(int)'), self.updateRodPlot)
        self.connect(self.ui.rmrodfilePB, SIGNAL('clicked()'), self.removeRodFile)
        self.connect(self.ui.addrodfitfilePB, SIGNAL('clicked()'),self.addRodFitFile)
        self.connect(self.ui.rodfitfileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedRodFitFile)
        self.connect(self.ui.rmrodfitfilePB, SIGNAL('clicked()'), self.removeRodFitFile)
        self.connect(self.ui.calrodCB, SIGNAL('stateChanged(int)'),self.updateRodCal)
        self.connect(self.ui.rodqoffLE,SIGNAL('returnPressed()'),self.updateRodCal)
        self.connect(self.ui.rodyscaleLE,SIGNAL('returnPressed()'),self.updateRodCal)
        self.connect(self.ui.rodsizeLE,SIGNAL('returnPressed()'),self.updateRodCal)
        self.connect(self.ui.rodsizeresLE,SIGNAL('returnPressed()'),self.updateRodCal)
        self.connect(self.ui.rodconLE,SIGNAL('returnPressed()'),self.updateRodCal)
        self.connect(self.ui.rodlinLE,SIGNAL('returnPressed()'),self.updateRodCal)
        self.connect(self.ui.rodroughLE,SIGNAL('returnPressed()'),self.updateRodCal)
        self.connect(self.ui.fitrodPB,SIGNAL('clicked()'),self.fitRod)
        self.connect(self.ui.rodconsPB, SIGNAL('clicked()'),self.updateRodPara)
        self.connect(self.ui.rodsaveCB, SIGNAL('activated(int)'), self.saveRod)
        self.connect(self.ui.rodloadCB, SIGNAL('activated(int)'), self.loadRod)
        self.connect(self.ui.rodscalePB, SIGNAL('clicked()'),self.setRodPlotScale)
        self.connect(self.ui.rodffPB, SIGNAL('clicked()'),self.formfactorShow)
        #flu analysis
        self.connect(self.ui.action_Open_FLU_file, SIGNAL('triggered()'),self.openFluFile)
        self.connect(self.ui.addflufilePB, SIGNAL('clicked()'),self.addFluFile)
        self.connect(self.ui.flufileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedFluFile)
        self.connect(self.ui.flulegendCB,SIGNAL('stateChanged(int)'), self.updateFluPlot)
        self.connect(self.ui.flulegendlocCoB,SIGNAL('currentIndexChanged(int)'), self.updateFluPlot)
        self.connect(self.ui.flulogyCB, SIGNAL('stateChanged(int)'), self.updateFluPlot)
        self.connect(self.ui.rmflufilePB, SIGNAL('clicked()'), self.removeFluFile)
        self.connect(self.ui.addflufitfilePB, SIGNAL('clicked()'),self.addFluFitFile)
        self.connect(self.ui.flufitfileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedFluFitFile)
        self.connect(self.ui.rmflufitfilePB, SIGNAL('clicked()'), self.removeFluFitFile)
        self.connect(self.ui.insflusubPB, SIGNAL('clicked()'), self.insFluIon)
        self.connect(self.ui.rmflusubPB, SIGNAL('clicked()'), self.rmFluIon)
        self.connect(self.ui.calfluCB, SIGNAL('stateChanged(int)'),self.updateFluCal)
        self.connect(self.ui.flushoqcCB, SIGNAL('stateChanged(int)'),self.updateFluPlot)
        self.connect(self.ui.flusurLE,SIGNAL('returnPressed()'),self.updateFluCal)
        self.connect(self.ui.fluqoffLE,SIGNAL('returnPressed()'),self.updateFluCal)
        self.connect(self.ui.fluyscaleLE,SIGNAL('returnPressed()'),self.updateFluCal)
        self.connect(self.ui.fluconLE,SIGNAL('returnPressed()'),self.updateFluCal)
        self.connect(self.ui.flulinLE,SIGNAL('returnPressed()'),self.updateFluCal)
        self.connect(self.ui.flusurcurLE,SIGNAL('returnPressed()'),self.updateFluCal)
        self.connect(self.ui.flufitPB,SIGNAL('clicked()'),self.fitFlu)
        self.connect(self.ui.fluconsPB, SIGNAL('clicked()'),self.updateFluPara)
        self.connect(self.ui.flusaveCB, SIGNAL('activated(int)'), self.saveFlu)
        self.connect(self.ui.fluloadCB, SIGNAL('activated(int)'), self.loadFlu) 
        self.connect(self.ui.fluscalePB, SIGNAL('clicked()'),self.setFluPlotScale)
        self.connect(self.ui.flubulLE, SIGNAL('returnPressed()'), self.updateFluQc)
        self.connect(self.ui.flurhobotLE, SIGNAL('returnPressed()'), self.updateFluQc)
        self.connect(self.ui.flurhotopLE, SIGNAL('returnPressed()'), self.updateFluQc)
        self.connect(self.ui.flusubTW,SIGNAL('cellChanged(int,int)'), self.updateFluQc)
        #gixos analysis
        self.connect(self.ui.action_Open_GIX_file, SIGNAL('triggered()'),self.openGixFile)
        self.connect(self.ui.addgixfilePB, SIGNAL('clicked()'),self.addGixFile)
        self.connect(self.ui.gixfileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedGixFile)
        self.connect(self.ui.gixlegendCB,SIGNAL('stateChanged(int)'), self.updateGixPlot)
        self.connect(self.ui.gixlegendlocCoB,SIGNAL('currentIndexChanged(int)'), self.updateGixPlot)
        self.connect(self.ui.gixlogyCB, SIGNAL('stateChanged(int)'), self.updateGixPlot)
        self.connect(self.ui.rmgixfilePB, SIGNAL('clicked()'), self.removeGixFile)
        self.connect(self.ui.gixscalePB, SIGNAL('clicked()'),self.setGixPlotScale)
        self.connect(self.ui.addgixfitfilePB, SIGNAL('clicked()'),self.addGixFitFile)
        self.connect(self.ui.gixfitfileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedGixFitFile)
        self.connect(self.ui.rmgixfitfilePB, SIGNAL('clicked()'), self.removeGixFitFile)
        self.connect(self.ui.addgixedfilePB, SIGNAL('clicked()'),self.addGixEDFile)
        self.connect(self.ui.gixedfileLW, SIGNAL('itemSelectionChanged()'),self.updateSelectedGixEDFile)
        self.connect(self.ui.rmgixedfilePB, SIGNAL('clicked()'), self.removeGixEDFile)
        self.connect(self.ui.gixedlegendCB,SIGNAL('stateChanged(int)'), self.updateGixEDPlot)
        self.connect(self.ui.gixedlegendlocCoB,SIGNAL('currentIndexChanged(int)'), self.updateGixEDPlot)
        self.connect(self.ui.gixedscalePB, SIGNAL('clicked()'), self.setGixEDPlotScale)
        self.connect(self.ui.insgixslabPB, SIGNAL('clicked()'), self.insGixSlab)
        self.connect(self.ui.rmgixslabPB, SIGNAL('clicked()'), self.rmGixSlab)
        self.connect(self.ui.gixnumslabSB, SIGNAL('valueChanged(int)'), self.modGixSlab)
        self.connect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self.updateGixParaVal)
        self.connect(self.ui.gixparTW,SIGNAL('cellDoubleClicked(int,int)'),self.setupGixPara)
        self.connect(self.ui.calgixCB, SIGNAL('stateChanged(int)'),self.updateGixCal)
        self.connect(self.ui.calgixsldCB, SIGNAL('stateChanged(int)'),self.updateGixCal)
        self.connect(self.ui.gixqoffLE,SIGNAL('returnPressed()'),self.updateGixCal)
        self.connect(self.ui.gixyscaleLE,SIGNAL('returnPressed()'),self.updateGixCal)
        self.connect(self.ui.gixqmaxLE,SIGNAL('returnPressed()'),self.updateGixCal)
        self.connect(self.ui.clegixconPB, SIGNAL('clicked()'),self.cleGixCon)
        self.connect(self.ui.gixsysconPB, SIGNAL('clicked()'),self.updateGixSysPara)
        self.connect(self.ui.fitgixPB,SIGNAL('clicked()'),self.fitGix)
        self.connect(self.ui.gixsaveCB, SIGNAL('activated(int)'), self.saveGix)
        self.connect(self.ui.gixloadCB, SIGNAL('activated(int)'), self.loadGix)
        
        
    def openRefFile(self):  #open ref files and also remove all current ref files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple REF Files to import', directory=self.directory, filter='RRF Files (*.rrf*;*_rrf.txt;*_rrf0.txt)')
        self.ui.tabWidget.setCurrentIndex(0)
        self.reffiles=map(str, f)
        self.directory=str(QFileInfo(self.reffiles[0]).absolutePath())
        self.updateRefFile()
            
    def updateRefFile(self): #update ref files in the listwidget
        self.ui.reffileLW.clear()
        for i in range(len(self.reffiles)):
            try:
                self.ui.reffileLW.addItem('#'+str(i+1)+self.halftab+str(self.reffiles[i].split('\\')[-2])+'\\'+str(self.reffiles[i].split('\\')[-1]))
            except:
                self.ui.reffileLW.addItem('#'+str(i+1)+self.halftab+str(self.reffiles[i].split('/')[-2])+'/'+str(self.reffiles[i].split('/')[-1]))

            
    def addRefFile(self): #add ref files into the listwidget and deselect all ref files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple REF Files to import', directory=self.directory, filter='RRF Files (*.rrf*;*_rrf.txt;*_rrf0.txt)')
        self.reffiles=self.reffiles+map(str, f)
        self.directory=str(QFileInfo(self.reffiles[0]).absolutePath())
        self.updateRefFile()
        
    def updateSelectedRefFile(self): #update the selected ref files in the listwidget
        selectedreffiles=self.ui.reffileLW.selectedItems()
        self.selectedreffiles_rows=[]
        for item in selectedreffiles:
            self.selectedreffiles_rows.append(self.ui.reffileLW.row(item))
        self.selectedreffiles_rows.sort()
        self.refscale=[[1,0,1,0] for i in range(len(self.selectedreffiles_rows))]
        self.updateRefPlot()
        
    def removeRefFile(self): #remove ref files in the listwidget and deselect all ref files in the listwidget
        items=self.ui.reffileLW.selectedItems()
        for item in items:
            self.reffiles.pop(self.ui.reffileLW.row(item))
        self.ui.reffileLW.clear()
        self.updateRefFile()
           
    def updateRefFitFile(self): #update ref fit files in the listwidget
        self.ui.reffitfileLW.clear()
        for i in range(len(self.reffitfiles)):
            try:
                self.ui.reffitfileLW.addItem('#'+str(i+1)+self.halftab+str(self.reffitfiles[i].split('\\')[-2])+'\\'+str(self.reffitfiles[i].split('\\')[-1]))
            except:
                self.ui.reffitfileLW.addItem('#'+str(i+1)+self.halftab+str(self.reffitfiles[i].split('/')[-2])+'/'+str(self.reffitfiles[i].split('/')[-1]))

    def addRefFitFile(self): #add ref fit files into the listwidget and deselect ref fit files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple REF Fit Files to import', directory=self.directory, filter='FIT Files (*.fit*;*_fit.txt;*_fit0.txt)')
        self.reffitfiles=self.reffitfiles+map(str, f)
        self.directory=str(QFileInfo(self.reffitfiles[0]).absolutePath())
        self.updateRefFitFile()
        
    def updateSelectedRefFitFile(self): #update the selected ref fit files in the listwidget
        selectedreffitfiles=self.ui.reffitfileLW.selectedItems()
        self.selectedreffitfiles_rows=[]
        for item in selectedreffitfiles:
            self.selectedreffitfiles_rows.append(self.ui.reffitfileLW.row(item))
        self.selectedreffitfiles_rows.sort()
        self.reffitscale=[[1,0,1,0] for i in range(len(self.selectedreffitfiles_rows))]
        self.updateRefPlot()
        
    def removeRefFitFile(self):  #remove ref fit files in the listwidget and deselect all ref fit files in the listwidget
        items=self.ui.reffitfileLW.selectedItems()
        for item in items:
            self.reffitfiles.pop(self.ui.reffitfileLW.row(item))
        self.ui.reffitfileLW.clear()
        self.updateRefFitFile()
        
    def updateRefPlot(self): #update the plot in the ref plotwidget
        self.ui.refPW.canvas.ax.clear()
        self.ui.refPW.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.ui.refPW.canvas.ax.set_ylabel('Normalized Reflectivity')
        if  len(self.selectedreffiles_rows)!=0: #plot ref files
            for i in range(len(self.selectedreffiles_rows)):
                data1=np.loadtxt(str(self.reffiles[self.selectedreffiles_rows[i]]), comments='#')
                self.ui.refPW.canvas.ax.errorbar(data1[:,0]*self.refscale[i][0]+self.refscale[i][1],data1[:,1]*self.refscale[i][2]+self.refscale[i][3],data1[:,2]*self.refscale[i][2],fmt='o',label='#'+str(self.selectedreffiles_rows[i]+1))
        if  len(self.selectedreffitfiles_rows)!=0: #plot ref fit files
            for i in range(len(self.selectedreffitfiles_rows)):
                data1=np.loadtxt(str(self.reffitfiles[self.selectedreffitfiles_rows[i]]), comments='#')
                self.ui.refPW.canvas.ax.errorbar(data1[:,0]*self.reffitscale[i][0]+self.reffitscale[i][1],data1[:,1]*self.reffitscale[i][2]+self.reffitscale[i][3],fmt='-',label='#'+str(self.selectedreffitfiles_rows[i]+1))
        if  self.ui.calrefCB.checkState()!=0:
                self.ui.refPW.canvas.ax.errorbar(np.array(self.refcal)[:,0],np.array(self.refcal)[:,1],fmt='-', label='cal')
        if self.ui.reflegendCB.checkState()!=0:
            self.ui.refPW.canvas.ax.legend(loc=self.ui.reflegendlocCoB.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        if self.ui.reflogyCB.checkState()!=0:
            self.ui.refPW.canvas.ax.set_yscale('log')
        else:
            self.ui.refPW.canvas.ax.set_yscale('linear')
        self.ui.refPW.canvas.draw()
        
    def setRefPlotScale(self): #set the scale of each data in the ref plot
        if len(self.selectedreffiles_rows)+len(self.selectedreffitfiles_rows)==0:
            self.messageBox('Warning:: No Ref or Fit files selected!')
        else:
            row_ref=len(self.selectedreffiles_rows)
            row_fit=len(self.selectedreffitfiles_rows)
            row=row_ref+row_fit
            Dialog=QDialog(self)
            self.uiplotscale=uic.loadUi('plotscale.ui', Dialog)
            self.uiplotscale.scaleTW.setRowCount(row) #set the table size; 4 column is fixed
            self.uiplotscale.show()
            self.uiplotscale.scaleLabel.setText('Reflectvity Plot Scale Setup: X=X*Factor+Offset')
            self.uiplotscale.scaleTW.setHorizontalHeaderLabels(QStringList()<<"X Factor"<<"X Offset"<<"Y Factor"<<"Y Offset") #set the horizontal header
            vlabel=QStringList() #set the vertical header 
            for i in range(row_ref):
                vlabel.append("Ref #"+str(self.selectedreffiles_rows[i]+1))
            for i in range(row_fit):
                vlabel.append("Fit #"+str(self.selectedreffitfiles_rows[i]+1))
            self.uiplotscale.scaleTW.setVerticalHeaderLabels(vlabel)
            for i in range(row_ref):  #set the initial values
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i,j,QTableWidgetItem(str(self.refscale[i][j])))
                    self.uiplotscale.scaleTW.item(i,j).setTextAlignment(Qt.AlignCenter)
            for i in range(row_fit):
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i+row_ref,j,QTableWidgetItem(str(self.reffitscale[i][j])))
                    self.uiplotscale.scaleTW.item(i+row_ref,j).setTextAlignment(Qt.AlignCenter)
            self.connect(self.uiplotscale.scaleTW, SIGNAL('cellChanged(int,int)'), self.updateRefPlotScale) #update the ref scale and plot
            self.connect(self.uiplotscale.closePB,SIGNAL('clicked()'), self.closePlotScale) #close the scale setup window
           
                       
    def updateRefPlotScale(self): #update the scale of each data in the ref plot
        row_ref=len(self.selectedreffiles_rows)
        row_fit=len(self.selectedreffitfiles_rows)
        self.refscale=[[float(str(self.uiplotscale.scaleTW.item(i,j).text())) for j in range(4)] for i in range(row_ref)]
        self.reffitscale=[[float(str(self.uiplotscale.scaleTW.item(i+row_ref,j).text())) for j in range(4)] for i in range(row_fit)]
        self.updateRefPlot()
        
    def closePlotScale(self): #close the plot scale window'
        self.uiplotscale.close()

    def updateRefEDFile(self): #update ed files in the listwidget
        self.ui.refedfileLW.clear()
        for i in range(len(self.refedfiles)):
            try:
                self.ui.refedfileLW.addItem('#'+str(i+1)+self.halftab+str(self.refedfiles[i].split('\\')[-2])+'\\'+str(self.refedfiles[i].split('\\')[-1]))
            except:
                self.ui.refedfileLW.addItem('#'+str(i+1)+self.halftab+str(self.refedfiles[i].split('/')[-2])+'/'+str(self.refedfiles[i].split('/')[-1]))

    
    def addRefEDFile(self): #add ref ed files into the listwidget and deselect all ed files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple ED Files to import', directory=self.directory, filter='FIT Files (*.sld*; *.ed*;*_sld.txt;*_ed.txt)')
        self.refedfiles=self.refedfiles+map(str, f)
        self.directory=str(QFileInfo(self.refedfiles[0]).absolutePath())
        self.updateRefEDFile()
        
    def updateSelectedRefEDFile(self): #update the selected ed files in the listwidget
        selectedrefedfiles=self.ui.refedfileLW.selectedItems()
        self.selectedrefedfiles_rows=[]
        for item in selectedrefedfiles:
            self.selectedrefedfiles_rows.append(self.ui.refedfileLW.row(item))
        self.selectedrefedfiles_rows.sort()
        self.refedscale=[[1,0,1,0] for i in range(len(self.selectedrefedfiles_rows))]
        self.updateRefEDPlot()
        
    def removeRefEDFile(self):  #remove ed files in the listwidget and deselect all ref fit files in the listwidget
        items=self.ui.refedfileLW.selectedItems()
        for item in items:
            self.refedfiles.pop(self.ui.refedfileLW.row(item))
        self.ui.refedfileLW.clear()
        self.updateRefEDFile()
        
    def updateRefEDPlot(self): #update the plot in the ref ed plotwidget
        self.ui.refsldPW.canvas.ax.clear()
        self.ui.refsldPW.canvas.ax.set_xlabel('Z'+' '+r'$[\AA]$')
        self.ui.refsldPW.canvas.ax.set_ylabel('Electron Density Profile'+' '+r'$[e/\AA^{3}]$')
        if  len(self.selectedrefedfiles_rows)!=0: #plot ref ed files
            for i in range(len(self.selectedrefedfiles_rows)):
                data1=np.loadtxt(str(self.refedfiles[self.selectedrefedfiles_rows[i]]), comments='#')
                self.ui.refsldPW.canvas.ax.errorbar(data1[:,0]*self.refedscale[i][0]+self.refedscale[i][1],data1[:,1]*self.refedscale[i][2]+self.refedscale[i][3],fmt='-',label='#'+str(self.selectedrefedfiles_rows[i]+1))
        if  self.ui.calsldCB.checkState()!=0:
                self.ui.refsldPW.canvas.ax.errorbar(np.array(self.sldcal)[:,0],np.array(self.sldcal)[:,1],fmt='-', label='cal')
        if self.ui.refedlegendCB.checkState()!=0:
            self.ui.refsldPW.canvas.ax.legend(loc=self.ui.refedlegendlocCoB.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        self.ui.refsldPW.canvas.draw()
        
    def setEDPlotScale(self): #set the scale of each data in the ed plot
        if len(self.selectedrefedfiles_rows)==0:
            self.messageBox('Warning:: No electron density files selected!')
        else:
            row_ed=len(self.selectedrefedfiles_rows)
            Dialog=QDialog(self)
            self.uiplotscale=uic.loadUi('plotscale.ui', Dialog)
            self.uiplotscale.scaleTW.setRowCount(row_ed) #set the table size; 4 column is fixed
            self.uiplotscale.show()
            self.uiplotscale.scaleLabel.setText('Electron Density Plot Scale Setup: X=X*Factor+Offset')
            self.uiplotscale.scaleTW.setHorizontalHeaderLabels(QStringList()<<"X Factor"<<"X Offset"<<"Y Factor"<<"Y Offset") #set the horizontal header
            vlabel=QStringList() #set the vertical header 
            for i in range(row_ed):
                vlabel.append("ED #"+str(self.selectedrefedfiles_rows[i]+1))
            self.uiplotscale.scaleTW.setVerticalHeaderLabels(vlabel)
            for i in range(row_ed):  #set the initial values
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i,j,QTableWidgetItem(str(self.refedscale[i][j])))
                    self.uiplotscale.scaleTW.item(i,j).setTextAlignment(Qt.AlignCenter)
            self.connect(self.uiplotscale.scaleTW,SIGNAL('cellChanged(int,int)'), self.updateEDPlotScale)
            self.connect(self.uiplotscale.closePB,SIGNAL('clicked()'), self.closePlotScale)
           
                       
    def updateEDPlotScale(self): #update the scale of each data in the ed plot
        row_ed=len(self.selectedrefedfiles_rows)
        self.refedscale=[[float(str(self.uiplotscale.scaleTW.item(i,j).text())) for j in range(4)] for i in range(row_ed)]
        self.updateRefEDPlot()
        
    def initRefPar(self): #initialize the refpar table
        self.ui.refparTW.horizontalHeader().setVisible(True)
        self.ui.refparTW.verticalHeader().setVisible(True)
        self.ui.refparTW.setHorizontalHeaderLabels(QStringList()<<'d ('+u'\u212b'+')'<<u'\u03c1'+' (e/'+u'\u212b'+u'\u00b3'+')'<<u'\u03bc'+' (cm'+u'\u207b'+u'\u00b9'+')'<<u'\u03c3'+' ('+u'\u212b'+')')
        top='top/0/0/3'
        bottom='bottom/0.333/0/NA'
        for i in range(4):
            self.ui.refparTW.setItem(0,i,QTableWidgetItem(top.split('/')[i]))
            self.ui.refparTW.setItem(1,i,QTableWidgetItem(bottom.split('/')[i]))
        self.ui.refnumslabSB.setValue(0)
        self.refpara={}  #initialize the parameter dictionary
        self.refpara[0]=[0,False, None,None]
        self.refpara[1]=[0,False, None,None]
        self.refpara[2]=[3,False, None,None]
        self.refpara[3]=[0.333,False, None,None]
        self.refpara[4]=[0,False, None,None]
        self.refsyspara={}  #initialize the ref system parameter dictonary
        self.refsyspara[0]=[float(self.ui.refqoffLE.text()), False, None, None]
        self.refsyspara[1]=[float(self.ui.refyscaleLE.text()), False, None, None]
        self.refsyspara[2]=[float(self.ui.refqresLE.text()), False, None, None]
        self.refsysCB=[self.ui.refqoffCB, self.ui.refyscaleCB, self.ui.refqresCB]
        self.updateRefParaName()
        
    def updateRefParaName(self):
        top=['rho_t','mu_t','sigma0']
        middle=[]
        bottom=['rho_b','mu_b']
        for i in range(self.ui.refparTW.rowCount()-2):
            layer=str(i+1)
            middle.extend(['d'+layer,'rho'+layer,'mu'+layer,'sigma'+layer])
        self.refparaname=top+middle+bottom
        self.refsysparaname=['q_off', 'y_scale', 'q_res']
        
    def insRefSlab(self, selrows=None): #add a slab in refpar table
        if selrows==None:
            insrows=self.ui.refparTW.selectionModel().selectedRows()
            insrows=[self.ui.refparTW.row(self.ui.refparTW.itemFromIndex(insrows[i])) for i in range(len(insrows))]
        else:
            insrows=[selrows]
        if len(insrows)!=1:
            self.messageBox('Warning:: Only one row can be selected!')
        elif insrows[0]==0:
            self.messageBox('Warning:: Cannot insert a layer above the top phase!')
        else:
            self.disconnect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self.updateRefParaVal)
            insrow=insrows[0]
            self.ui.refparTW.insertRow(insrow)
            for i in range(4):
                self.ui.refparTW.setItem(insrow,i,QTableWidgetItem('10/0.3/0/3'.split('/')[i]))
            self.connect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self.updateRefParaVal)
            self.ui.refnumslabSB.setValue(self.ui.refparTW.rowCount()-2)
            for i in list(reversed(range((insrow-1)*4+3,4*(self.ui.refparTW.rowCount()-3)+5))):   #update the parameter dictionary    
                self.refpara[i+4]=self.refpara[i]
            self.addRefParaDic(insrow)
            self.updateRefParaName()  #update the paramenter name list
          #  print self.refparaname
           # print self.refpara
            self.updateRefCal()            
            
    def addRefParaDic(self,row):
        self.refpara[(row-1)*4+3]=[10,False, None,None]
        self.refpara[(row-1)*4+4]=[0.3,False, None,None]
        self.refpara[(row-1)*4+5]=[0,False, None,None]
        self.refpara[(row-1)*4+6]=[3,False, None,None]
                    
                
            
    def rmRefSlab(self, selrows=None): #remove multiple slabs in refpar table
        row=self.ui.refparTW.rowCount()
        rmrows=self.ui.refparTW.selectionModel().selectedRows()   
        removerows=[]
        if selrows==None:
            for rmrow in rmrows:
                removerows.append(self.ui.refparTW.row(self.ui.refparTW.itemFromIndex(rmrow)))
                removerows.sort(reverse=True)
        else:
                removerows=selrows
        if len(removerows)==0:
            self.messageBox('Warning:: No layer is selected')
        else:
            for i in range(len(removerows)):
                if removerows[i] == 0:
                    self.messageBox('Warning:: Cannot remove the top phase!')
                elif removerows[i] == row-1:
                    self.messageBox('Warning:: Cannot remove the bottom phase!')
                else:       
                    self.ui.refparTW.removeRow(removerows[i])
                    for i in range(removerows[i]*4+3,len(self.refpara)):  #update the parameter dictionary
                        self.refpara[i-4]=self.refpara[i]                   #shift the parameters below the deleting row up
                    for key in range(len(self.refpara)-4,len(self.refpara)): #delete the last four parameters
                        self.refpara.pop(key)
          #  print self.refpara
            self.updateRefParaName()  #update the paramenter name list
            self.ui.refnumslabSB.setValue(self.ui.refparTW.rowCount()-2)
            self.updateRefCal()

            
    def modRefSlab(self): #modify refpar table based on the change of spin box
        diff=self.ui.refparTW.rowCount()-self.ui.refnumslabSB.value()-2
        row=self.ui.refparTW.rowCount()
        if diff>0:
            selrows=[]
            for i in range(diff):
                selrows.append(row-2-i)
               # self.ui.refparTW.removeRow(row-2-i)
           # print selrows
            self.rmRefSlab(selrows)
        elif diff<0:
            for i in range(-diff):
                self.insRefSlab(row-1)
#                self.ui.refparTW.insertRow(row-1)
#                for j in range(4):
#                    self.ui.refparTW.setItem(row-1,j,QTableWidgetItem('10/0.3/0/3'.split('/')[j]))
           # self.updateRefCal()
        
    def updateRefParaVal(self):
        #selrow=self.ui.refparTW.currentRow()
        selcol=self.ui.refparTW.currentColumn()
        #print selrow, selcol
#        if selrow==self.ui.refparTW.rowCount()-1:
#            paranum=selrow*4+selcol-2
#        else:
#            paranum=selrow*4+selcol-1
#        self.refpara[paranum][0]= float(str(self.ui.refparTW.item(selrow,selcol).text()))    # update the current value for selected cell in the ref parameter dictionary
        #print self.refpara
        if self.ui.refroughCB.checkState()!=0 and selcol==3:
            self.sameRough()  #fix all roughness
        self.updateRefCal()
        
    def sameRough(self):
        row=self.ui.refparTW.rowCount()
        samerough=float(str(self.ui.refparTW.item(0,3).text()))
        self.disconnect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self. updateRefParaVal)
        for i in range(1,row-1):
            self.ui.refparTW.setItem(i,3,QTableWidgetItem(str(samerough)))
            self.refpara[i*4+2][0]=samerough
        self.connect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self. updateRefParaVal)

            
    
    def updateRefCal(self): # caluate the Ref and Sld based on current parameters.
        row=self.ui.refparTW.rowCount()       
        d=[float(str(self.ui.refparTW.item(i+1,0).text())) for i in range(row-2)]
        rho=[float(str(self.ui.refparTW.item(i,1).text())) for i in range(row)]
        mu=[float(str(self.ui.refparTW.item(i,2).text())) for i in range(row)]
        sigma=[float(str(self.ui.refparTW.item(i,3).text())) for i in range(row-1)]
        if self.refsavedataindex==1:  # save the data after qoff correction & fit 
            syspara=[0,float(self.ui.refyscaleLE.text()),float(self.ui.refqresLE.text())]  
            data=np.loadtxt(str(self.reffiles[self.selectedreffiles_rows[0]]), comments='#')
            x=data[:,0]
            y=data[:,1]
            yerr=data[:,2]
            qoff=float(self.ui.refqoffLE.text())
            xnew=x+qoff
            xref=np.linspace(xnew[0],xnew[-1],800)
            self.refcal=np.vstack((xref,self.refCalFun(d,rho,mu,sigma,syspara,xref))).T  # fit data
            lamda=12.3984/float(self.ui.xenLE.text())  #correction for ref data
            frsnll,frsnl1=xr.parratt(x,lamda,[0,1],[rho[0],rho[-1]],[mu[0]/4/np.pi/1e8*lamda,mu[0]/4/np.pi/1e8*lamda]) 
            frsnllnew,frsnl1new=xr.parratt(xnew,lamda,[0,1],[rho[0],rho[-1]],[mu[0]/4/np.pi/1e8*lamda,mu[0]/4/np.pi/1e8*lamda])
            ynew=y*frsnll/frsnllnew
            yerrnew=yerr*frsnll/frsnllnew
            self.refdata=np.vstack((xnew,ynew,yerrnew)).T
        else:
            syspara=[float(self.ui.refqoffLE.text()),float(self.ui.refyscaleLE.text()),float(self.ui.refqresLE.text())] 
           # self.refedxmin=-4*sigma[0]
          #  self.refedxmax=np.sum(d)+4*sigma[-1]
            if self.refsavefitindex==1:
                xref=np.linspace(self.refxmin,self.refxmax,self.refnp)
                self.refcal=np.vstack((xref,self.refCalFun(d,rho,mu,sigma,syspara,xref))).T
            elif self.refsavefitindex==2:
                xsld=np.linspace(self.refedxmin,self.refedxmax,self.refnp)
                self.sldcal=np.vstack((xsld,self.sldCalFun(d,rho,sigma,xsld))).T
            else:
                if  self.ui.calrefCB.checkState()!=0:
                    if  len(self.selectedreffiles_rows)!=0: 
                        data=np.loadtxt(str(self.reffiles[self.selectedreffiles_rows[0]]), comments='#')
                        self.refxmax=np.max(data[:,0])
                        self.refxmin=np.min(data[:,0])
                    else:
                        self.refxmax=0.7
                        self.refxmin=0
                    xref=np.linspace(self.refxmin,self.refxmax,800)
                   # self.refcal=[[xref[i],self.refCalFun(d,rho,mu,sigma,xref[i])] for i in range(len(xref))]
                    self.refcal=np.vstack((xref,self.refCalFun(d,rho,mu,sigma,syspara,xref))).T
                self.updateRefPlot()
                if  self.ui.calsldCB.checkState()!=0:
                    if sigma[0]!=0 and sigma[-1]!=0:
                        xsld=np.linspace(-4*sigma[0],np.sum(d)+4*sigma[-1],800)
                    else:
                        xsld=np.linspace(-10,np.sum(d)+10,800)
                    #self.sldcal=[[xsld[i],self.sldCalFun(d,rho,sigma,xsld[i])] for i in range(len(xsld))]
                    self.sldcal=np.vstack((xsld,self.sldCalFun(d,rho,sigma,xsld))).T
                self.updateRefEDPlot()

        
    def refCalFun(self,d,rho,mu,sigma,syspara,x):
        qoff=syspara[0] 
        yscale=syspara[1]
        qres=syspara[2] 
        d=[abs(d[i]) for i in range(len(d))]
        rho=[abs(rho[i]) for i in range(len(rho))]
        mu=[abs(mu[i]) for i in range(len(mu))]
        sigma=[abs(sigma[i]) for i in range(len(sigma))]
        erad=self.eleradius   # classic electron radius
        slab=0.25
        k0=2*np.pi*float(self.ui.xenLE.text())/12.3984 # wave vector
        theta=x/2/k0   # convert q to theta
        length=np.sum(d)+4*(sigma[0]+sigma[-1]) # total length of inner slabs plus 4 times rougness for both sides
        steps=int(length/slab) # each sliced box has thickness of ~ 0.25 \AA
        xsld=np.linspace(-4*sigma[0],np.sum(d)+4*sigma[-1],steps) # get the x-axis for sld
        intrho=self.sldCalFun(d,rho,sigma,xsld)
        intmu=self.sldCalFun(d,mu,sigma,xsld)
        sd=length/steps # thickness for each slab
        sdel=[]
        sbet=[] 
       # t0=time.time()
        sdel.append(erad*2.0*np.pi/k0/k0*rho[0]) # delta for the top phase
        sbet.append(mu[0]/2/k0/1e8)        # beta for the top phase
        sdel=sdel+[intrho[i]*erad*2.0*np.pi/k0/k0 for i in range(len(intrho))] # add delta for the interface
        sbet=sbet+[intmu[i]/2/k0/1e8 for i in range(len(intmu))] # add beta for the interface
        sdel.append(erad*2.0*np.pi/k0/k0*rho[-1])  # delta for the bottom phase
        sbet.append(mu[-1]/2/k0/1e8)    # beta for the bottom phase         
        d=slab*np.ones_like(sdel)
        lamda=2*np.pi/k0
        fdel=erad*2.0*np.pi/k0/k0
        sdelf=np.array(sdel)/fdel
        #t1=time.time()
        #print 't1-t0=',t1-t0
        ref,refr=xr.parratt(x+qoff,lamda,d,sdelf,sbet)
        frsnll,frsnl1=xr.parratt(x,lamda,[0,1],[sdelf[0],sdelf[-1]],[sbet[0],sbet[-1]])        
        #t2=time.time()
        #print 't2-t1=', t2-t1
        # the following is to use python to calculate the ref, seems at least 10 times slower than f90. 
#        rrf=[]        
#        for i in range(len(theta)):
#            f=[]
#            an=[]
#            rn=[]
#            rn.append(0)
#            f.append(cmath.sqrt(complex(theta[i]*theta[i],-2*sbet[0]))) # f for the top phase
#            for j in range(steps):
#                f.append(cmath.sqrt(complex(theta[i]*theta[i]-2*(sdel[j+1]-sdel[0]), -2*sbet[j+1])))
#            f.append(cmath.sqrt(complex(theta[i]*theta[i]-2*(sdel[-1]-sdel[0]), -2*sbet[-1])))  # f for the bottom phase
#            for j in range(steps+1):
#                an.append((f[steps+1-j-1]-f[steps+1-j])/(f[steps+1-j-1]+f[steps+1-j]))
#                rn.append(cmath.exp(f[steps+1-j-1]*complex(0,-2*k0*sd))*(rn[j]+an[j])/(rn[j]*an[j]+1))
#            ref=abs(rn[-1])*abs(rn[-1])
#            rn1=(f[0]-f[-1])/(f[0]+f[-1])
#            frsnll=abs(rn1)*abs(rn1)
#            rrf.append(ref/frsnll)
        #t3=time.time()
        #print 't3-t2=', t3-t2                        
        #return ref/frsnll
        if self.ui.refrrfCB.checkState()!=0:
            return yscale*ref/frsnll
        else:
            return yscale*ref
        
    def sldCalFun(self,d,y,sigma,x):
        wholesld=[]
        if self.ui.refroughCB.checkState()!=0: #for fixed roughness
            for i in range(1,len(sigma)):
                sigma[i]=sigma[0]
        for i in range(len(sigma)):
            if sigma[i]<=0:
                sigma[i]=1e-5
        for j in range(len(x)):
            pos=[]
            erfx=[]
            pos.append(0)
            erfx.append(x[j]/sigma[0]/math.sqrt(2))
            for i in range(len(d)):
                pos.append(pos[i]+d[i])
                erfx.append((x[j]-pos[i+1])/sigma[i+1]/math.sqrt(2))
            sld=0
            for i in range(len(sigma)):
                sld=sld+math.erf(erfx[i])*(y[i+1]-y[i])
            if (sld+y[0]+y[-1])/2<0:
                wholesld.append(0)
            else:
                wholesld.append((sld+y[0]+y[-1])/2)
        return wholesld
        
    def saveRef(self):
        if str(self.ui.refsaveCB.currentText())=='Save Fit':
            self.refsavefitindex=1
            self.saveRefFitDig()
        elif str(self.ui.refsaveCB.currentText())=='Save ED':
            self.refsavefitindex=2
            self.saveRefFitDig()
        elif str(self.ui.refsaveCB.currentText())=='Save Para':
            self.saveRefPara()
        elif str(self.ui.refsaveCB.currentText())=='Save Data':
            self.saveRefData()
            
    def saveRefData(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Ref Data & Fit after the Q Offset Correction',directory=self.directory))
        fname1=self.saveFileName+'_rrf0.txt'
        fname2=self.saveFileName+'_fit0.txt'
        self.refsavedataindex=1
        try:
            self.updateRefCal()
            np.savetxt(fname1,self.refdata,fmt='%.4f\t%.4e\t%.4e')
            np.savetxt(fname2,self.refcal,fmt='%.4f\t%.4e')
            self.refsavedataindex=0
        except:
            self.refsavedataindex=0
        
    def saveRefFitDig(self):
        Dialog=QDialog(self)
        self.uirefsavefit=uic.loadUi('refsave.ui', Dialog)
        if self.refsavefitindex==1:
            self.uirefsavefit.label.setText('Save Reflectvity Fit/Calcualtion!')
            try:
                self.uirefsavefit.xminLE.setText(str(self.refxmin))
                self.uirefsavefit.xmaxLE.setText(str(self.refxmax))
            except:
                pass
        elif self.refsavefitindex==2: 
            row=self.ui.refparTW.rowCount()       
            d=[float(str(self.ui.refparTW.item(i+1,0).text())) for i in range(row-2)]
            sigma=[float(str(self.ui.refparTW.item(i,3).text())) for i in range(row-1)]
            self.uirefsavefit.label.setText('Save Electron Density Profile!')
            self.uirefsavefit.xminLE.setText(str(-4*sigma[0]))
            self.uirefsavefit.xmaxLE.setText(str(np.sum(d)+4*sigma[-1]))
        self.uirefsavefit.numpointLE.setText(str(400))
        self.uirefsavefit.show()
        self.connect(self.uirefsavefit.cancelPB, SIGNAL('clicked()'), self.cancelSaveRefFit)
        self.connect(self.uirefsavefit.okPB, SIGNAL('clicked()'), self.saveRefFit)
        
    def cancelSaveRefFit(self):
        self.uirefsavefit.close()
        self.refsavefitindex=0
        
    def saveRefFit(self):
        self.refnp=float(self.uirefsavefit.numpointLE.text())
        if float(self.uirefsavefit.xminLE.text())>=float(self.uirefsavefit.xmaxLE.text()) or float(self.uirefsavefit.numpointLE.text())<=0:
            self.messageBox('Warning::Starting points must be lower than ending points \n and numer of points must be large than 0!!')
        else:
            if self.refsavefitindex==1:       
                self.refxmin=float(self.uirefsavefit.xminLE.text())
                self.refxmax=float(self.uirefsavefit.xmaxLE.text())
                self.updateRefCal()
                self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Ref Fit Data',directory=self.directory))
                fname=self.saveFileName+'_fit.txt'
                np.savetxt(fname,self.refcal,fmt='%.4f\t%.4e')
            elif self.refsavefitindex==2:       
                self.refedxmin=float(self.uirefsavefit.xminLE.text())
                self.refedxmax=float(self.uirefsavefit.xmaxLE.text())
                self.updateRefCal()
                self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Electron Density Data',directory=self.directory))
                fname=self.saveFileName+'_ed.txt'
                np.savetxt(fname,self.sldcal,fmt='%.4e\t%.4e')
            self.refsavefitindex=0
            self.uirefsavefit.close()
    
        
    def setupRefPara(self): # constrains setup for ref parameters
        Dialog=QDialog(self)
        self.uirefpara=uic.loadUi('refpara.ui', Dialog)
        selrow=self.ui.refparTW.currentRow()
        selcol=self.ui.refparTW.currentColumn()
        if selrow==self.ui.refparTW.rowCount()-1:
            self.paranum=selrow*4+selcol-2
        else:
            self.paranum=selrow*4+selcol-1
        self.uirefpara.label.setText('Limits Setup of Parameter:'+self.refparaname[self.paranum])
        if self.refpara[self.paranum][2]!=None: 
            self.uirefpara.minCB.setCheckState(2)
            self.uirefpara.minLE.setText(str(self.refpara[self.paranum][2]))
        if self.refpara[self.paranum][3]!=None: 
            self.uirefpara.maxCB.setCheckState(2)
            self.uirefpara.maxLE.setText(str(self.refpara[self.paranum][3]))
        self.uirefpara.show()
        self.connect(self.uirefpara.cancelPB, SIGNAL('clicked()'), self.cancelRefPara)
        self.connect(self.uirefpara.okPB, SIGNAL('clicked()'), self.takeRefPara)
        
    def cancelRefPara(self):
        self.uirefpara.close()
        
    def takeRefPara(self):
        if self.uirefpara.minCB.checkState()!=0 and self.uirefpara.maxCB.checkState()!=0 and float(self.uirefpara.minLE.text())>float(self.uirefpara.maxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain!!!")
        else:
            if self.uirefpara.minCB.checkState()!=0:
                self.refpara[self.paranum][2]=float(self.uirefpara.minLE.text())
            else:
                self.refpara[self.paranum][2]=None
            if self.uirefpara.maxCB.checkState()!=0:
                self.refpara[self.paranum][3]=float(self.uirefpara.maxLE.text())
            else:
                self.refpara[self.paranum][3]=None
            self.uirefpara.close()
    
    def cleRefCon(self):
        for i in range(len(self.refpara)):
            self.refpara[i][2]=None
            self.refpara[i][3]=None
        for i in range(len(self.refsyspara)):
            self.refsyspara[i][2]=None
            self.refsyspara[i][3]=None
            
    def updateRefSysPara(self):
        Dialog=QDialog(self)
        self.uirefsyspara=uic.loadUi('refsyspara.ui', Dialog)
        if self.refsyspara[0][2]!=None:  #set up the current value
            self.uirefsyspara.qoffminCB.setCheckState(2)
            self.uirefsyspara.qoffminLE.setText(str(self.refsyspara[0][2]))
        if self.refsyspara[0][3]!=None:  
            self.uirefsyspara.qoffmaxCB.setCheckState(2)
            self.uirefsyspara.qoffmaxLE.setText(str(self.refsyspara[0][3]))
        if self.refsyspara[1][2]!=None:  
            self.uirefsyspara.yscaleminCB.setCheckState(2)
            self.uirefsyspara.yscaleminLE.setText(str(self.refsyspara[1][2]))
        if self.refsyspara[1][3]!=None:  
            self.uirefsyspara.yscalemaxCB.setCheckState(2)
            self.uirefsyspara.yscalemaxLE.setText(str(self.refsyspara[1][3]))
        if self.refsyspara[2][2]!=None:  
            self.uirefsyspara.qresminCB.setCheckState(2)
            self.uirefsyspara.qresminLE.setText(str(self.refsyspara[2][2]))
        if self.refsyspara[2][3]!=None:  
            self.uirefsyspara.qresmaxCB.setCheckState(2)
            self.uirefsyspara.qresmaxLE.setText(str(self.refsyspara[2][3]))
        self.uirefsyspara.show()
        self.connect(self.uirefsyspara.cancelPB, SIGNAL('clicked()'), self.cancelRefSysPara)
        self.connect(self.uirefsyspara.okPB, SIGNAL('clicked()'), self.takeRefSysPara)
        
    def cancelRefSysPara(self):
        self.uirefsyspara.close()
                
    def takeRefSysPara(self):
        if self.uirefsyspara.qoffminCB.checkState()!=0 and self.uirefsyspara.qoffmaxCB.checkState()!=0 and float(self.uirefsyspara.qoffminLE.text())>float(self.uirefsyspara.qoffmaxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain for Q offset!!!")
        elif self.uirefsyspara.yscaleminCB.checkState()!=0 and self.uirefsyspara.yscalemaxCB.checkState()!=0 and float(self.uirefsyspara.yscaleminLE.text())>float(self.uirefsyspara.yscalemaxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain for Y scale!!!")
        elif self.uirefsyspara.qresminCB.checkState()!=0 and self.uirefsyspara.qresmaxCB.checkState()!=0 and float(self.uirefsyspara.qresminLE.text())>float(self.uirefsyspara.qresmaxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain for Q resolution!!!")
        else:
            if self.uirefsyspara.qoffminCB.checkState()!=0:
                self.refsyspara[0][2]=float(self.uirefsyspara.qoffminLE.text())
            else:
                self.refsyspara[0][2]=None
            if self.uirefsyspara.qoffmaxCB.checkState()!=0:
                self.refsyspara[0][3]=float(self.uirefsyspara.qoffmaxLE.text())
            else:
                self.refsyspara[0][3]=None
            if self.uirefsyspara.yscaleminCB.checkState()!=0:
                self.refsyspara[1][2]=float(self.uirefsyspara.yscaleminLE.text())
            else:
                self.refsyspara[1][2]=None
            if self.uirefsyspara.yscalemaxCB.checkState()!=0:
                self.refsyspara[1][3]=float(self.uirefsyspara.yscalemaxLE.text())
            else:
                self.refsyspara[1][3]=None
            if self.uirefsyspara.qresminCB.checkState()!=0:
                self.refsyspara[2][2]=float(self.uirefsyspara.qresminLE.text())
            else:
                self.refsyspara[2][2]=None
            if self.uirefsyspara.qresmaxCB.checkState()!=0:
                self.refsyspara[2][3]=float(self.uirefsyspara.qresmaxLE.text())
            else:
                self.refsyspara[2][3]=None
            self.uirefsyspara.close()
         
    def getRefParaVal(self):
        for i in range(len(self.refparaname)-2): #get the current values except the bottom phase in the table 
                cell=divmod(i+1,4)  #get the cell index for each parameter
                self.refpara[i][0]=float(str(self.ui.refparTW.item(cell[0],cell[1]).text()))
        self.refpara[len(self.refparaname)-2][0]=float(str(self.ui.refparTW.item(cell[0]+1,1).text()))   #last row
        self.refpara[len(self.refparaname)-1][0]=float(str(self.ui.refparTW.item(cell[0]+1,2).text()))
        self.refsyspara[0][0]=float(self.ui.refqoffLE.text())  #system parameters
        self.refsyspara[1][0]=float(self.ui.refyscaleLE.text())
        self.refsyspara[2][0]=float(self.ui.refqresLE.text())
    
    
    def fitRef(self):
        self.getRefParaVal()
        index=self.ui.refparTW.selectionModel().selectedIndexes()
        row=self.ui.refparTW.rowCount()
        selrows=[self.ui.refparTW.row(self.ui.refparTW.itemFromIndex(index[i])) for i in range(len(index))]
        selcols=[self.ui.refparTW.column(self.ui.refparTW.itemFromIndex(index[i])) for i in range(len(index))]
        selparas=[]
        for i in range(len(selrows)):  #get selected parameters
            if selrows[i]!=row-1:
                selparas.append(selrows[i]*4+selcols[i]-1)
            else:
                selparas.append(selrows[i]*4+selcols[i]-2)
       # print selparas
        for i in range(len(self.refpara)):  #set selected parameters to be varied
            if i in selparas:
                self.refpara[i][1]=True
            else:
                self.refpara[i][1]=False
        #print self.refpara
        if self.ui.refqoffCB.checkState()!=0:   #set selected system paramenters to be varied
            self.refsyspara[0][1]=True
        else:
            self.refsyspara[0][1]=False
        if self.ui.refyscaleCB.checkState()!=0:  
            self.refsyspara[1][1]=True
        else:
            self.refsyspara[1][1]=False
        if self.ui.refqresCB.checkState()!=0:  
            self.refsyspara[2][1]=True
        else:
            self.refsyspara[2][1]=False
       # print self.refsyspara
        self.refparameter=Parameters()
        for i in range(len(self.refpara)):
            self.refparameter.add(self.refparaname[i], value=self.refpara[i][0],vary=self.refpara[i][1],min=self.refpara[i][2],max=self.refpara[i][3])
        for i in range(len(self.refsysparaname)):
            self.refparameter.add(self.refsysparaname[i], value=self.refsyspara[i][0],vary=self.refsyspara[i][1],min=self.refsyspara[i][2],max=self.refsyspara[i][3])
        #print self.refparameter
       
#        d=[self.refparameter[self.refparaname[i*4+3]].value for i in range(row-2)]
#        rho=[self.refparameter[self.refparaname[i*4]].value for i in range(row-1)]
#        mu=[self.refparameter[self.refparaname[i*4+1]].value for i in range(row-1)]
#        sigma=[self.refparameter[self.refparaname[i*4+2]].value for i in range(row-1)]
#        rho.append(self.refparameter[self.refparaname[-2]].value)  #add bottom phase
#        mu.append(self.refparameter[self.refparaname[-1]].value)  #add bottom phase
#        print d, rho, mu, sigma
       
        if  len(self.selectedreffiles_rows)!=1: #plot ref files
            self.messageBox("Please select only one set of data for fitting!")
        else:
            data=np.loadtxt(str(self.reffiles[self.selectedreffiles_rows[0]]), comments='#')
            ini=max(float(str(self.ui.reffitranLE.text()).split(':')[0]),data[0][0])
            fin=min(float(str(self.ui.reffitranLE.text()).split(':')[1]),data[-1][0])
            data1=data[np.where(np.logical_and(data[:,0]>=ini,data[:,0]<=fin))]
            x=data1[:,0]
            y=data1[:,1]
            if self.ui.referrCB.currentIndex()==0:
                yerr=data1[:,2]
            elif self.ui.referrCB.currentIndex()==1:
                yerr=np.sqrt(y)
            elif self.ui.referrCB.currentIndex()==2:
                yerr=y
            else:
                yerr=np.ones_like(x)
           # print yerr
            self.refresult=minimize(self.ref2min, self.refparameter, args=(x,y,yerr))
#            tempchi=self.refresult.redchi
#            print tempchi
#            self.refparameter=self.refresult.params
#            self.refresult=minimize(self.ref2min, self.refparameter, args=(x,y,yerr))
#            print self.refresult.redchi
#            while(np.abs(self.refresult.redchi-tempchi)/tempchi>0.02):
#                tempchi=self.refresult.redchi
#                self.refparameter=self.refresult.params
#                self.refresult=minimize(self.ref2min, self.refparameter,args=(x,y,yerr))
#                print self.refresult.redchi
            print(fit_report(self.refresult))
            residual=np.vstack((x,self.refresult.residual)).T
           # print residual
            self.disconnect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self.updateRefParaVal)
            if self.ui.refroughCB.checkState()!=0: #enforce the roughness to be same if set
                for i in range(1,row-1):
                    self.refresult.params[self.refparaname[4*i+2]].value=self.refresult.params[self.refparaname[2]].value
            for i in range(len(self.refparaname)-2): #put the best values except the bottom phase in the table 
                cell=divmod(i+1,4)  #get the cell index for each parameter
               # print str(result.params[self.refparaname[i]].value)
                self.ui.refparTW.setItem(cell[0],cell[1],QTableWidgetItem(format(self.refresult.params[self.refparaname[i]].value,'.4f')))
            self.ui.refparTW.setItem(row-1,1,QTableWidgetItem(format(self.refresult.params[self.refparaname[-2]].value, '.4f'))) # put the best values for the bottom phase
            self.ui.refparTW.setItem(row-1,2,QTableWidgetItem(format(self.refresult.params[self.refparaname[-1]].value, '.4f')))
            self.ui.refqoffLE.setText(format(self.refresult.params[self.refsysparaname[0]].value, '.6f'))  #put the best sys parameter values 
            self.ui.refyscaleLE.setText(format(self.refresult.params[self.refsysparaname[1]].value, '.3f'))
            self.ui.refqresLE.setText(format(self.refresult.params[self.refsysparaname[2]].value, '.6f'))
            self.connect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self.updateRefParaVal)
            self.ui.calrefCB.setCheckState(2)
            self.updateRefCal()
            self.ui.refChiLE.setText(format(self.refresult.redchi, '.3f'))
            self.ui.refparaTB.clear()
            fitinfo='Fitting Paramenters:\n'
            fitinfo=fitinfo+'Name\tStderr\tMin\tMax\n'
            selparas.sort()
            for i in selparas:
                fitinfo=fitinfo+self.refparaname[i]+'\t'+format(self.refresult.params[self.refparaname[i]].stderr, '.4f')+'\t'+str(self.refpara[i][2])+'\t'+str(self.refpara[i][3])+'\n'
            for i in range(3):
                if self.refsyspara[i][1]==True:
                   fitinfo=fitinfo+self.refsysparaname[i]+'\t'+format(self.refresult.params[self.refsysparaname[i]].stderr, '.4f')+'\t'+str(self.refsyspara[i][2])+'\t'+str(self.refsyspara[i][3])+'\n' 
            fitinfo=fitinfo+'********************************\n'
            fitinfo=fitinfo+'Fitting Residual:\n'
            for i in range(len(residual)):
                fitinfo=fitinfo+format(residual[i][0], '.3f')+'\t'+format(residual[i][1], '.4f')+'\n'
            self.ui.refparaTB.append(fitinfo)
            cursor=self.ui.refparaTB.textCursor()
            cursor.setPosition(0)
            self.ui.refparaTB.setTextCursor(cursor)
           # print self.refresult.params
            #print self.refresult.params[self.refparaname[-1]].stderr
            #print self.refparameter[self.refparaname[-1]]
        
    def ref2min(self, params, x, y, yerr): #residuel for ref fitting
        row=self.ui.refparTW.rowCount()
        d=[params[self.refparaname[i*4+3]].value for i in range(row-2)]
        rho=[params[self.refparaname[i*4]].value for i in range(row-1)]
        mu=[params[self.refparaname[i*4+1]].value for i in range(row-1)]
        sigma=[params[self.refparaname[i*4+2]].value for i in range(row-1)]
        rho.append(params[self.refparaname[-2]].value)  #add bottom phase
        mu.append(params[self.refparaname[-1]].value)  #add bottom phase
        syspara=[params[self.refsysparaname[i]].value for i in range(3)]
        model=self.refCalFun(d,rho,mu,sigma,syspara,x)
        return (model-y)/yerr
    
    def saveRefPara(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Reflectivity Fitting Parameters',directory=self.directory))
        fid=open(self.saveFileName+'_par.txt','w')
        try:
            fid.write('Chi_Square\t'+format(self.refresult.redchi, '.3f')+'\n')  #chisquare
        except:
            fid.write('Chi_Square\tNA\n') 
        fid.write('Error_Type\t'+str(self.ui.referrCB.currentText()).split()[0]+'\n')
        fid.write('Num_of_Layer\t'+str(self.ui.refparTW.rowCount()-2)+'\n')    #number of layers
        if self.ui.refroughCB.checkState()!=0:
            fid.write('Roughness\tFixed\n')
        else:
            fid.write('Roughness\tVary\n')
        fid.write('Para_Name\tValue\t\tVary\tStderr\t\tMin\tMax\n')
        for i in range((self.ui.refparTW.rowCount()-2)*4+5):
            try:
                fid.write(self.refparaname[i]+'\t\t'+format(self.refresult.params[self.refparaname[i]].value,'.3e')+'\t'+str(self.refpara[i][1])+'\t'+format(self.refresult.params[self.refparaname[i]].stderr,'.3e')+'\t'+str(self.refpara[i][2])+'\t'+str(self.refpara[i][3])+'\n')
            except:
                if i <=(self.ui.refparTW.rowCount()-2)*4+2:
                    cell=divmod(i+1,4)
                else:
                    cell=divmod(i+2,4)
                fid.write(self.refparaname[i]+'\t\t'+str(self.ui.refparTW.item(cell[0],cell[1]).text())+'\tNA\tNA\tNA\tNA\n')
        for i in range(3):
            try:
                fid.write(self.refsysparaname[i]+'\t\t'+format(self.refresult.params[self.refsysparaname[i]].value,'.3e')+'\t'+str(self.refsyspara[i][1])+'\t'+format(self.refresult.params[self.refsysparaname[i]].stderr,'.3e')+'\t'+str(self.refsyspara[i][2])+'\t'+str(self.refsyspara[i][3])+'\n')
            except:
                temp=[float(self.ui.refqoffLE.text()),float(self.ui.refyscaleLE.text()),float(self.ui.refqresLE.text())]   
                fid.write(self.refsysparaname[i]+'\t\t'+str(temp[i])+'\tNA\tNA\tNA\tNA\n')
        fid.close()
                
    def loadRef(self):
        if str(self.ui.refloadCB.currentText())=='Load Para':
            self.loadRefPara()
            
    def loadRefPara(self):
        filename=QFileDialog.getOpenFileName(caption='Select Parameter File to read', directory=self.directory, filter='Par Files (*.par*;*_par.txt)')
        self.directory=str(QFileInfo(filename).absolutePath())
        fid=open(filename)
        fdata=fid.readlines()
        fid.close()
        nlayer=eval(fdata[2].split('\t')[1])  #get number of layers
        roughness=fdata[3][:-1].split('\t')[1]  # get the roughness
        para=[]
        for i in range(5,len(fdata)):
            para.append(eval(fdata[i].split('\t')[2]))
        self.ui.calsldCB.setCheckState(0)
        self.ui.calrefCB.setCheckState(0)
        if roughness=='Fixed':
            self.ui.refroughCB.setCheckState(2)
        else:
            self.ui.refroughCB.setCheckState(0)
        self.ui.refnumslabSB.setValue(nlayer)
        self.disconnect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self. updateRefParaVal)
        for i in range(nlayer*4+5):
            if i<=nlayer*4+2:
                cell=divmod(i+1,4)
            else:
                cell=divmod(i+2,4)
            self.ui.refparTW.setItem(cell[0],cell[1],QTableWidgetItem(str(para[i])))
            self.refpara[i][0]=para[i]
        self.ui.refqoffLE.setText(format(para[-3],'.6f'))
        self.ui.refyscaleLE.setText(format(para[-2],'.3f'))
        self.ui.refqresLE.setText(format(para[-1],'.6f'))
        self.connect(self.ui.refparTW,SIGNAL('cellChanged(int,int)'), self. updateRefParaVal)
        self.ui.calsldCB.setCheckState(2)
        self.ui.calrefCB.setCheckState(2)
       # print self.refpara
    
    def errorCal(self):
        index=self.ui.refparTW.selectionModel().selectedIndexes()
        self.sysselparas=[]
        for i in range(len(self.refsyspara)):
            if self.refsysCB[i].checkState()!=0:
                self.sysselparas.append(i)
        #self.messageBox('Please have the best fitting value!')
        if  len(self.selectedreffiles_rows)!=1: 
            self.messageBox("Please select only one set of data for uncertainty calculation!")
        elif len(index)+len(self.sysselparas)!=1:
            self.messageBox("Please select only one parameter for uncertainty calculation!")
        else:
            self.getRefParaVal()   #get best fitting vaules and assige them to the new globle name
            self.refsysbestpara=[self.refsyspara[i][0] for i in range(len(self.refsyspara))]
            self.refbestpara=[self.refpara[i][0] for i in range(len(self.refpara))]
           # print self.refsysbestpara
           # print self.refbestpara
            Dialog=QDialog(self)
            self.uireferr1=uic.loadUi('referr1.ui',Dialog)
            if len(self.sysselparas)==1:
                self.uireferr1.label.setText('Uncertainty Calculation for Parameter:'+self.refsysparaname[self.sysselparas[0]])
                self.uireferr1.bestvalLE.setText(format(self.refsysbestpara[self.sysselparas[0]], '.2e'))
                self.uireferr1.leftranLE.setText(format(self.refsysbestpara[self.sysselparas[0]]*0.1, '.2e'))
                self.uireferr1.rightranLE.setText(format(self.refsysbestpara[self.sysselparas[0]]*0.1, '.2e'))
                self.refselparaname=self.refsysparaname[self.sysselparas[0]]
                self.paranum=-1
            else:
                selrow=self.ui.refparTW.currentRow()
                selcol=self.ui.refparTW.currentColumn()
                if selrow==self.ui.refparTW.rowCount()-1:
                    self.paranum=selrow*4+selcol-2
                else:
                    self.paranum=selrow*4+selcol-1
                self.uireferr1.label.setText('Uncertainty Calculation for Parameter:'+self.refparaname[self.paranum])
                self.uireferr1.bestvalLE.setText(format(self.refbestpara[self.paranum], '.4f'))
                self.uireferr1.leftranLE.setText(format(self.refbestpara[self.paranum]*0.1, '.4f'))
                self.uireferr1.rightranLE.setText(format(self.refbestpara[self.paranum]*0.1, '.4f'))
                self.refselparaname=self.refparaname[self.paranum]
                self.sysselparas=[-1]
            self.uireferr1.show()
            self.connect(self.uireferr1.cancelPB, SIGNAL('clicked()'), self.cancelerrorCal)
            self.connect(self.uireferr1.nextPB, SIGNAL('clicked()'), self.errorcalPara)
    
    def cancelerrorCal(self):
        self.uireferr1.close()
        
    def errorcalPara(self):
        self.referrbestval=float(self.uireferr1.bestvalLE.text())
        self.referrnumpoint=float(self.uireferr1.numpointLE.text())
        self.referrleftran=float(self.uireferr1.leftranLE.text())
        self.referrrightran=float(self.uireferr1.rightranLE.text())
        self.uireferr1.close()
        Dialog=QDialog(self)
        self.uireferr2=uic.loadUi('referr2.ui',Dialog)
        self.uireferr2.label.setText('Please select other parameters to fit when\ncalculation the uncertainty for '+self.refselparaname)
        self.uireferr2.show()
        self.connect(self.uireferr2.cancelPB, SIGNAL('clicked()'), self.cancelerrorcalPara)
        self.connect(self.uireferr2.nextPB, SIGNAL('clicked()'), self.errorcalCal)
        
    def cancelerrorcalPara(self):
        self.uireferr2.close()
        
    def errorcalCal(self):
        self.uireferr2.close()
        index=self.ui.refparTW.selectionModel().selectedIndexes()
        row=self.ui.refparTW.rowCount()
        selrows=[self.ui.refparTW.row(self.ui.refparTW.itemFromIndex(index[i])) for i in range(len(index))]
        selcols=[self.ui.refparTW.column(self.ui.refparTW.itemFromIndex(index[i])) for i in range(len(index))]
        selparas=[]
        for i in range(len(selrows)):  #get selected parameters
            if selrows[i]!=row-1:
                selparas.append(selrows[i]*4+selcols[i]-1)
            else:
                selparas.append(selrows[i]*4+selcols[i]-2)
        sysselparas=[]
        for i in range(len(self.refsyspara)):
            if self.refsysCB[i].checkState()!=0:
                sysselparas.append(i)
        try:
            selparas.remove(self.paranum)
        except:
            pass
        for i in range(len(self.refpara)):  #set selected parameters to be varied
            if i in selparas:
                self.refpara[i][1]=True
            else:
                self.refpara[i][1]=False
        for i in range(len(self.refsyspara)):  #set selected system parameters to be varied
            if i in sysselparas:
                self.refsyspara[i][1]=True
            else:
                self.refsyspara[i][1]=False
        data=np.loadtxt(str(self.reffiles[self.selectedreffiles_rows[0]]), comments='#')    #get data
        ini=max(float(str(self.ui.reffitranLE.text()).split(':')[0]),data[0][0])
        fin=min(float(str(self.ui.reffitranLE.text()).split(':')[1]),data[-1][0])
        data1=data[np.where(np.logical_and(data[:,0]>=ini,data[:,0]<=fin))]
        x=data1[:,0]
        y=data1[:,1]
        if self.ui.referrCB.currentIndex()==0:
            yerr=data1[:,2]
        elif self.ui.referrCB.currentIndex()==1:
            yerr=np.sqrt(y)
        elif self.ui.referrCB.currentIndex()==2:
            yerr=y
        else:
            yerr=np.ones_like(x)
        self.referrx=np.linspace(self.referrbestval-self.referrleftran, self.referrbestval+self.referrrightran,int(self.referrnumpoint+1)) #get x range for this parameter
        print self.referrx
        self.referr=[]
        self.referrx1=self.referrx[np.where(self.referrx>=self.referrbestval)]
        self.referrx2=self.referrx[np.where(self.referrx<self.referrbestval)][::-1]
        print self.referrx1, self.referrx2
        
        self.progressDialog=QProgressDialog('Calculating Chi-square','Abort',0,100)
        self.progressDialog.setWindowModality(Qt.WindowModal)
        self.progressDialog.setWindowTitle('Wait')
        self.progressDialog.setAutoClose(True)
        self.progressDialog.setAutoReset(True)
        self.progressDialog.setMinimum(1)
        self.progressDialog.setMaximum(len(self.referrx))
        self.progressDialog.show()   
        
        for i in range(len(self.referrx2)):
            self.referrpara=Parameters()
            for j in range(len(self.refpara)):
                if j==self.paranum:
                        self.referrpara.add(self.refparaname[j], value=self.referrx2[i],vary=self.refpara[j][1],min=self.refpara[j][2],max=self.refpara[j][3])
                else:
                    if i==0:
                        self.referrpara.add(self.refparaname[j], value=self.refbestpara[j],vary=self.refpara[j][1],min=self.refpara[j][2],max=self.refpara[j][3])
                    else:
                        self.referrpara.add(self.refparaname[j], value=self.reftemppara[j],vary=self.refpara[j][1],min=self.refpara[j][2],max=self.refpara[j][3])
            for j in range(len(self.refsysparaname)):
                if j==self.sysselparas[0]:
                    self.referrpara.add(self.refsysparaname[j], value=self.referrx2[i],vary=self.refsyspara[j][1],min=self.refsyspara[j][2],max=self.refsyspara[j][3])
                else:
                    if i==0:
                        self.referrpara.add(self.refsysparaname[j], value=self.refsysbestpara[j],vary=self.refsyspara[j][1],min=self.refsyspara[j][2],max=self.refsyspara[j][3])
                    else:
                        self.referrpara.add(self.refsysparaname[j], value=self.refsystemppara[j],vary=self.refsyspara[j][1],min=self.refsyspara[j][2],max=self.refsyspara[j][3])
            print self.referrpara
            self.referrresult=minimize(self.ref2min, self.referrpara,args=(x,y,yerr))
            self.referr.append(np.array([self.referrx2[i],self.referrresult.redchi]))
            print self.referrx2[i], self.referrresult.redchi
            self.refsystemppara=[self.referrresult.params[self.refsysparaname[j]].value for j in range(len(self.refsyspara))]
            self.reftemppara=[self.referrresult.params[self.refparaname[j]].value for j in range(len(self.refpara))]
            self.progressDialog.setLabelText('Calculating Chi-square for '+self.refselparaname+' at '+format(self.referrx2[i],'.4f'))    #update the progress dialog 
            self.updateProgress()
            if self.progressDialog.wasCanceled()==True:
                break
        
        for i in range(len(self.referrx1)):
            self.referrpara=Parameters()
            for j in range(len(self.refpara)):
                if j==self.paranum:
                        self.referrpara.add(self.refparaname[j], value=self.referrx1[i],vary=self.refpara[j][1],min=self.refpara[j][2],max=self.refpara[j][3])
                else:
                    if i==0:
                        self.referrpara.add(self.refparaname[j], value=self.refbestpara[j],vary=self.refpara[j][1],min=self.refpara[j][2],max=self.refpara[j][3])
                    else:
                        self.referrpara.add(self.refparaname[j], value=self.reftemppara[j],vary=self.refpara[j][1],min=self.refpara[j][2],max=self.refpara[j][3])
            for j in range(len(self.refsysparaname)):
                if j==self.sysselparas[0]:
                    self.referrpara.add(self.refsysparaname[j], value=self.referrx1[i],vary=self.refsyspara[j][1],min=self.refsyspara[j][2],max=self.refsyspara[j][3])
                else:
                    if i==0:
                        self.referrpara.add(self.refsysparaname[j], value=self.refsysbestpara[j],vary=self.refsyspara[j][1],min=self.refsyspara[j][2],max=self.refsyspara[j][3])
                    else:
                        self.referrpara.add(self.refsysparaname[j], value=self.refsystemppara[j],vary=self.refsyspara[j][1],min=self.refsyspara[j][2],max=self.refsyspara[j][3])
            print self.referrpara
            self.referrresult=minimize(self.ref2min, self.referrpara,args=(x,y,yerr))
            self.referr.append(np.array([self.referrx1[i],self.referrresult.redchi]))
            print self.referrx1[i], self.referrresult.redchi
            self.refsystemppara=[self.referrresult.params[self.refsysparaname[j]].value for j in range(len(self.refsyspara))]
            self.reftemppara=[self.referrresult.params[self.refparaname[j]].value for j in range(len(self.refpara))]
            self.progressDialog.setLabelText('Calculating Chi-square for '+self.refselparaname+' at '+format(self.referrx1[i],'.4f'))    #update the progress dialog 
            self.updateProgress()
            if self.progressDialog.wasCanceled()==True:
                break
       
        self.progressDialog.hide() 
        self.referr=np.array(self.referr)
        self.referr=self.referr[np.argsort(self.referr[:,0])]
        f=interp1d(self.errorlist[:,0],self.errorlist[:,1],kind='cubic')
        chisqfac=f(self.referrresult.nfree)
        f1=interp1d(self.referr[:,0],self.referr[:,1],kind='linear')
        x=np.linspace(self.referrx[0],self.referrx[-1],200)
        self.referrex=np.vstack((x,f1(x))).T
        minchi=np.min(self.referr[:,1])
        tarchi=minchi*chisqfac
        self.referrtar=np.array([[self.referrx[0],tarchi],[self.referrx[-1],tarchi]])
        minchipos=np.where(self.referr[:,1]==np.min(self.referr[:,1]))[0][0]
        try:
            f2=interp1d(self.referr[:,1][:minchipos+1],self.referr[:,0][:minchipos+1],kind='linear')
            lefterr=format(self.referrbestval-f2(tarchi),'.4f')
        except:
            lefterr='not found'
        try:
            f3=interp1d(self.referr[:,1][minchipos:],self.referr[:,0][minchipos:],kind='linear')
            righterr=format(f3(tarchi)-self.referrbestval,'.4f')      
        except:
            righterr='not found'
            
        
        Dialog=QDialog(self)
        self.uireferr3=uic.loadUi('referr3.ui',Dialog)
        self.uireferr3.label.setText('Plot for Chi-square vs Parameter: '+self.refselparaname)
        self.uireferr3.show()
        self.uireferr3.minchiLE.setText(format(minchi,'.2f'))
        self.uireferr3.tarchiLE.setText(format(tarchi,'.2f'))
        self.uireferr3.lefterrLE.setText(lefterr)
        self.uireferr3.righterrLE.setText(righterr)
        self.connect(self.uireferr3.logyCB,SIGNAL('stateChanged(int)'), self.referrcalPlot)
        self.connect(self.uireferr3.closePB, SIGNAL('clicked()'), self.closeerrorcalPlot)
        self.connect(self.uireferr3.savePB, SIGNAL('clicked()'), self.errorcalSave)        
        
        self.referrcalPlot()
        
    def referrcalPlot(self): 
        self.uireferr3.plotWidget.canvas.ax.clear()
        self.uireferr3.plotWidget.canvas.ax.set_xlabel(self.refselparaname)
        self.uireferr3.plotWidget.canvas.ax.set_ylabel('Chi-square')
        self.uireferr3.plotWidget.canvas.ax.errorbar(self.referr[:,0], self.referr[:,1], fmt='o')
        self.uireferr3.plotWidget.canvas.ax.errorbar(self.referrex[:,0], self.referrex[:,1], fmt='-')
        self.uireferr3.plotWidget.canvas.ax.errorbar(self.referrtar[:,0], self.referrtar[:,1], fmt='-')
        if self.uireferr3.logyCB.checkState()!=0:
            self.uireferr3.plotWidget.canvas.ax.set_yscale('log')
        else:
            self.uireferr3.plotWidget.canvas.ax.set_yscale('linear')
        self.uireferr3.plotWidget.canvas.draw()
        
    def closeerrorcalPlot(self):
        self.uireferr3.close()
        
    def errorcalSave(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Chi-sqaure Data for Parameter'+self.refselparaname,directory=self.directory))
        fid=open(self.saveFileName+'_chi.txt','w')
        fid.write('#target chi-square:\t'+self.uireferr3.tarchiLE.text()+'\n')
        fid.write('#left error:\t'+self.uireferr3.lefterrLE.text()+'\n')
        fid.write('#right error:\t'+self.uireferr3.righterrLE.text()+'\n')
        for i in range(len(self.referrex)):
            fid.write(format(self.referrex[i][0], '.4f')+'\t'+format(self.referrex[i][1],'.4f')+'\n')
        fid.close()
        
        
    def updateProgress(self):
        self.progressDialog.setValue(self.progressDialog.value()+1)  
    
        ################################################        
        #state the rod analysis section. 
        ################################################
    
    def openRodFile(self):  #open ref files and also remove all current ref files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple ROD Files to import', directory=self.directory, filter='ROD Files (*.rod*;*.cut*;*_cut.txt)')
        self.ui.tabWidget.setCurrentIndex(1)
        self.rodfiles=map(str, f)
        self.directory=str(QFileInfo(self.rodfiles[0]).absolutePath())
        self.updateRodFile()
    
    def updateRodFile(self): #update rod files in the listwidget
        self.ui.rodfileLW.clear()
        for i in range(len(self.rodfiles)):
            try:
                self.ui.rodfileLW.addItem('#'+str(i+1)+self.halftab+str(self.rodfiles[i].split('\\')[-2])+'\\'+str(self.rodfiles[i].split('\\')[-1]))
            except:
                self.ui.rodfileLW.addItem('#'+str(i+1)+self.halftab+str(self.rodfiles[i].split('/')[-2])+'/'+str(self.rodfiles[i].split('/')[-1]))

    def addRodFile(self): #add rod files into the listwidget and deselect all rod files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Rod Files to import', directory=self.directory, filter='Rod Files (*.rod*;*.cut*;*_cut.txt)')
        self.rodfiles=self.rodfiles+map(str, f)
        self.directory=str(QFileInfo(self.rodfiles[0]).absolutePath())
        self.updateRodFile()
        
    def updateSelectedRodFile(self): #update the selected rod files in the listwidget
        selectedrodfiles=self.ui.rodfileLW.selectedItems()
        self.selectedrodfiles_rows=[]
        for item in selectedrodfiles:
            self.selectedrodfiles_rows.append(self.ui.rodfileLW.row(item))
        self.selectedrodfiles_rows.sort()
        self.rodscale=[[1,0,1,0] for i in range(len(self.selectedrodfiles_rows))]
        self.updateRodPlot()
        
    def removeRodFile(self): #remove rod files in the listwidget and deselect all rod files in the listwidget
        items=self.ui.rodfileLW.selectedItems()
        for item in items:
            self.rodfiles.pop(self.ui.rodfileLW.row(item))
        self.ui.rodfileLW.clear()
        self.updateRodFile()
           
    def updateRodFitFile(self): #update rod fit files in the listwidget
        self.ui.rodfitfileLW.clear()
        for i in range(len(self.rodfitfiles)):
            try:
                self.ui.rodfitfileLW.addItem('#'+str(i+1)+self.halftab+str(self.rodfitfiles[i].split('\\')[-2])+'\\'+str(self.rodfitfiles[i].split('\\')[-1]))
            except:
                self.ui.rodfitfileLW.addItem('#'+str(i+1)+self.halftab+str(self.rodfitfiles[i].split('/')[-2])+'/'+str(self.rodfitfiles[i].split('/')[-1]))

    def addRodFitFile(self): #add rod fit files into the listwidget and deselect rod fit files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Rod Fit Files to import', directory=self.directory, filter='FIT Files (*.fit*; *_fit.txt)')
        self.rodfitfiles=self.rodfitfiles+map(str, f)
        self.directory=str(QFileInfo(self.rodfitfiles[0]).absolutePath())
        self.updateRodFitFile()
        
    def updateSelectedRodFitFile(self): #update the selected rod fit files in the listwidget
        selectedrodfitfiles=self.ui.rodfitfileLW.selectedItems()
        self.selectedrodfitfiles_rows=[]
        for item in selectedrodfitfiles:
            self.selectedrodfitfiles_rows.append(self.ui.rodfitfileLW.row(item))
        self.selectedrodfitfiles_rows.sort()
        self.rodfitscale=[[1,0,1,0] for i in range(len(self.selectedrodfitfiles_rows))]
        self.updateRodPlot()
        
    def removeRodFitFile(self):  #remove rod fit files in the listwidget and deselect all rod fit files in the listwidget
        items=self.ui.rodfitfileLW.selectedItems()
        for item in items:
            self.rodfitfiles.pop(self.ui.rodfitfileLW.row(item))
        self.ui.rodfitfileLW.clear()
        self.updateRodFitFile()
        
    def updateRodPlot(self): #update the plot in the rod plotwidget
        self.ui.rodPW.canvas.ax.clear()
        self.ui.rodPW.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.ui.rodPW.canvas.ax.set_ylabel('Intensity [a.u.]')
        if  len(self.selectedrodfiles_rows)!=0: #plot rod files
            for i in range(len(self.selectedrodfiles_rows)):
                data1=np.loadtxt(str(self.rodfiles[self.selectedrodfiles_rows[i]]), comments='#')
                self.ui.rodPW.canvas.ax.errorbar(data1[:,0]*self.rodscale[i][0]+self.rodscale[i][1],data1[:,1]*self.rodscale[i][2]+self.rodscale[i][3],data1[:,2]*self.rodscale[i][2],fmt='o',label='#'+str(self.selectedrodfiles_rows[i]+1))
        if  len(self.selectedrodfitfiles_rows)!=0: #plot rod fit files
            for i in range(len(self.selectedrodfitfiles_rows)):
                data1=np.loadtxt(str(self.rodfitfiles[self.selectedrodfitfiles_rows[i]]), comments='#')
                self.ui.rodPW.canvas.ax.errorbar(data1[:,0]*self.rodfitscale[i][0]+self.rodfitscale[i][1],data1[:,1]*self.rodfitscale[i][2]+self.rodfitscale[i][3],fmt='-',label='#'+str(self.selectedrodfitfiles_rows[i]+1))
        if  self.ui.calrodCB.checkState()!=0:
                self.ui.rodPW.canvas.ax.errorbar(np.array(self.rodcal)[:,0],np.array(self.rodcal)[:,1],fmt='-', label='cal')
        if self.ui.rodlegendCB.checkState()!=0:
            self.ui.rodPW.canvas.ax.legend(loc=self.ui.rodlegendlocCoB.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        if self.ui.rodlogyCB.checkState()!=0:
            self.ui.rodPW.canvas.ax.set_yscale('log')
        else:
            self.ui.rodPW.canvas.ax.set_yscale('linear')
        self.ui.rodPW.canvas.draw()
    
    def setRodPlotScale(self): #set the scale of each data in the rod plot
        if len(self.selectedrodfiles_rows)+len(self.selectedrodfitfiles_rows)==0:
            self.messageBox('Warning:: No Rod or Fit files selected!')
        else:
            row_rod=len(self.selectedrodfiles_rows)
            row_fit=len(self.selectedrodfitfiles_rows)
            row=row_rod+row_fit
            Dialog=QDialog(self)
            self.uiplotscale=uic.loadUi('plotscale.ui', Dialog)
            self.uiplotscale.scaleTW.setRowCount(row) #set the table size; 4 column is fixed
            self.uiplotscale.show()
            self.uiplotscale.scaleLabel.setText('Rod Plot Scale Setup: X=X*Factor+Offset')
            self.uiplotscale.scaleTW.setHorizontalHeaderLabels(QStringList()<<"X Factor"<<"X Offset"<<"Y Factor"<<"Y Offset") #set the horizontal header
            vlabel=QStringList() #set the vertical header 
            for i in range(row_rod):
                vlabel.append("Rod #"+str(self.selectedrodfiles_rows[i]+1))
            for i in range(row_fit):
                vlabel.append("Fit #"+str(self.selectedrodfitfiles_rows[i]+1))
            self.uiplotscale.scaleTW.setVerticalHeaderLabels(vlabel)
            for i in range(row_rod):  #set the initial values
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i,j,QTableWidgetItem(str(self.rodscale[i][j])))
                    self.uiplotscale.scaleTW.item(i,j).setTextAlignment(Qt.AlignCenter)
            for i in range(row_fit):
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i+row_rod,j,QTableWidgetItem(str(self.rodfitscale[i][j])))
                    self.uiplotscale.scaleTW.item(i+row_rod,j).setTextAlignment(Qt.AlignCenter)
            self.connect(self.uiplotscale.scaleTW, SIGNAL('cellChanged(int,int)'), self.updateRodPlotScale) #update the rod scale and plot
            self.connect(self.uiplotscale.closePB,SIGNAL('clicked()'), self.closePlotScale) #close the scale setup window
           
                       
    def updateRodPlotScale(self): #update the scale of each data in the rod plot
        row_rod=len(self.selectedrodfiles_rows)
        row_fit=len(self.selectedrodfitfiles_rows)
        self.rodscale=[[float(str(self.uiplotscale.scaleTW.item(i,j).text())) for j in range(4)] for i in range(row_rod)]
        self.rodfitscale=[[float(str(self.uiplotscale.scaleTW.item(i+row_rod,j).text())) for j in range(4)] for i in range(row_fit)]
        self.updateRodPlot()
        
  
    
    def initRodPar(self): #initialize the rod parameters
        self.rodparaname=['q_off','y_scale','size_z','size_res','bg_con','bg_lin','roughness'] 
        self.uirodLE=[self.ui.rodqoffLE,self.ui.rodyscaleLE,self.ui.rodsizeLE,self.ui.rodsizeresLE,self.ui.rodconLE,self.ui.rodlinLE,self.ui.rodroughLE]
        self.uirodCB=[self.ui.rodqoffCB,self.ui.rodyscaleCB,self.ui.rodsizeCB,self.ui.rodsizeresCB,self.ui.rodconCB,self.ui.rodlinCB,self.ui.rodroughCB]       
        self.rodpara={}  #initialize the rod  parameter dictonary
        for i in range(len(self.rodparaname)):
            self.rodpara[i]=[float(self.uirodLE[i].text()), False, None, None]
      
    def updateRodParVal(self): #update the rod parameters value
        for i in range(len(self.rodparaname)):
            self.rodpara[i][0]=float(self.uirodLE[i].text())   
    
    def updateRodCal(self): # caluate the Rod  based on current parameters.
        if self.ui.rodlipidCB.checkState()+self.ui.rodNPCB.checkState()!=2:
            self.messageBox('Warning: Please select either lipids or NPs for rod calculation!')
        else:
            self.updateRodParVal()
            rodpara=[self.rodpara[i][0] for i in range(7)]
            if self.rodsavefitindex==1:
                xrod=np.linspace(max(0,self.rodxmin),self.rodxmax,self.rodnp)
                self.rodcal=np.vstack((xrod,self.rodCalFun(rodpara,xrod))).T
            else:
                if  self.ui.calrodCB.checkState()!=0:
                    if  len(self.selectedrodfiles_rows)!=0: 
                        data=np.loadtxt(str(self.rodfiles[self.selectedrodfiles_rows[0]]), comments='#')
                        self.rodxmax=np.max(data[:,0])
                        self.rodxmin=max(0,np.min(data[:,0]))
                    else:
                        self.rodxmax=1
                        self.rodxmin=0
                    xrod=np.linspace(self.rodxmin,self.rodxmax,1600)
                    self.rodcal=np.vstack((xrod,self.rodCalFun(rodpara,xrod))).T
                self.updateRodPlot()
            
        
    def rodCalFun(self,rodpara,x):
        qoff=rodpara[0] 
        yscale=rodpara[1]
        size=abs(rodpara[2]) 
        sizeres=abs(rodpara[3])
        bgcon=rodpara[4]
        bglin=rodpara[5]
        roughness=rodpara[6]
        dth=float(self.ui.roddthLE.text())/180*np.pi  #out-of-plane angle in rad
        alpha=float(self.ui.rodalphaLE.text())/180*np.pi  #incident angle in rad
        k0=2*np.pi*float(self.ui.rodxenLE.text())/12.3984 # wave vector
        q1=x+qoff  #used for formfactor and roughness calculation   
        q2=np.array(x+qoff-k0*np.sin(alpha)) #used tansmission calculation
       # beta=np.arcsin(q2/k0)  #outgoing angle in rad
        beta=[np.arcsin(q2[i]/k0) for i in range(len(q2))]
        q3=k0*np.sqrt(2+2*np.sin(alpha)*np.sin(beta)-2*np.cos(alpha)*np.cos(beta)*np.cos(dth))  #total q for NP formfactor calculation
        erad=self.eleradius   # classic electron radius
        qc=2*np.sqrt(np.pi*erad*float(self.ui.rodrhoLE.text())) # half critical q for tansmission calculation
        rod=[] #return value
        for i in range(len(q1)):
            if q2[i]<0:
                rod.append(0) # return 0 for beta less zero
            else:
                fmax=cmath.sqrt(complex(q2[i]*q2[i]-qc*qc,0))
                f1=complex(q2[i],0)
                trans=4*abs(f1/(f1+fmax))*abs(f1/(f1+fmax))
                if self.ui.rodlipidCB.checkState()!=0:
                    formfac=special.sph_jn(0,q1[i]*size/2)[0][0]*special.sph_jn(0,q1[i]*size/2)[0][0] #use Bessel function here.
                else:
                    if sizeres==0:
                        formfac=9*(special.sph_jn(1,q3[i]*size)[0][1]/(q3[i]*size))**2
                    else:
                        formfac=quad(lambda t: 9*(special.sph_jn(1,q3[i]*t)[0][1]/(q3[i]*t))**2*np.exp(-(t-size)**2/2/sizeres**2), size-2.82*sizeres, size+2.82*sizeres)[0]/np.sqrt(2*np.pi)/sizeres
                rod.append(trans*formfac*yscale*np.exp(-q1[i]**2*roughness**2)+bgcon+bglin*q1[i])
        return rod
                
    def fitRod(self):
        if self.ui.rodlipidCB.checkState()+self.ui.rodNPCB.checkState()!=2:
            self.messageBox('Warning: Please select either lipids or NPs for rod calculation!')
        else: 
            self.updateRodParVal()
            for i in range(len(self.rodpara)):
                if self.uirodCB[i].checkState()!=0:  #set selected paramenters to be varied
                    self.rodpara[i][1]=True
                else:
                    self.rodpara[i][1]=False
            parastatus=np.array([self.rodpara[i][1] for i in range(len(self.rodpara))])
            selparas=np.where(parastatus==True)
            self.rodparameter=Parameters()
            for i in range(len(self.rodpara)):
                self.rodparameter.add(self.rodparaname[i], value=self.rodpara[i][0],vary=self.rodpara[i][1],min=self.rodpara[i][2],max=self.rodpara[i][3])
            if  len(self.selectedrodfiles_rows)!=1: #plot ref files
                self.messageBox("Please select only one set of data for fitting!")
            else:
                data=np.loadtxt(str(self.rodfiles[self.selectedrodfiles_rows[0]]), comments='#')
                ini=max(float(str(self.ui.rodfitranLE.text()).split(':')[0]),data[0][0])
                fin=min(float(str(self.ui.rodfitranLE.text()).split(':')[1]),data[-1][0])
                data1=data[np.where(np.logical_and(data[:,0]>=ini,data[:,0]<=fin))]
                x=data1[:,0]
                y=data1[:,1]
                if self.ui.roderrCB.currentIndex()==0:
                    yerr=data1[:,2]
                elif self.ui.roderrCB.currentIndex()==1:
                    yerr=np.sqrt(y)
                elif self.ui.roderrCB.currentIndex()==2:
                    yerr=y
                else:
                    yerr=np.ones_like(x)
               # print yerr
                self.rodresult=minimize(self.rod2min, self.rodparameter, args=(x,y,yerr))
                print(fit_report(self.rodresult))
                residual=np.vstack((x,self.rodresult.residual)).T
               # print residual
                self.ui.rodqoffLE.setText(format(self.rodresult.params[self.rodparaname[0]].value, '.6f'))
                self.ui.rodyscaleLE.setText(format(self.rodresult.params[self.rodparaname[1]].value, '.2e'))
                self.ui.rodsizeLE.setText(format(self.rodresult.params[self.rodparaname[2]].value, '.2f'))
                self.ui.rodsizeresLE.setText(format(self.rodresult.params[self.rodparaname[3]].value, '.2f'))
                self.ui.rodconLE.setText(format(self.rodresult.params[self.rodparaname[4]].value, '.2e'))
                self.ui.rodlinLE.setText(format(self.rodresult.params[self.rodparaname[5]].value, '.2e'))
                self.ui.rodroughLE.setText(format(self.rodresult.params[self.rodparaname[6]].value, '.2e'))
                self.ui.calrodCB.setCheckState(2)
                self.updateRodCal()
                self.ui.rodChiLE.setText(format(self.rodresult.redchi, '.3f'))
                self.ui.rodparaTB.clear()
                fitinfo='Fitting Paramenters:\n'
                fitinfo=fitinfo+'Name\tStderr\tMin\tMax\n'
                for i in selparas[0]:
                    fitinfo=fitinfo+self.rodparaname[i]+'\t'+format(self.rodresult.params[self.rodparaname[i]].stderr, '.4f')+'\t'+str(self.rodpara[i][2])+'\t'+str(self.rodpara[i][3])+'\n'
                fitinfo=fitinfo+'********************************\n'
                fitinfo=fitinfo+'Fitting Residual:\n'
                for i in range(len(residual)):
                    fitinfo=fitinfo+format(residual[i][0], '.3f')+'\t'+format(residual[i][1], '.4f')+'\n'
                self.ui.rodparaTB.append(fitinfo)
                cursor=self.ui.rodparaTB.textCursor()
                cursor.setPosition(0)
                self.ui.rodparaTB.setTextCursor(cursor)
    
    def rod2min(self, params, x, y, yerr): #residuel for rod fitting  
        rodpara=[params[self.rodparaname[i]] for i in range(len(params))]
        model=self.rodCalFun(rodpara,x)
        return (model-y)/yerr
        
    def updateRodPara(self):
        Dialog=QDialog(self)
        self.uirodpara=uic.loadUi('rodpara.ui', Dialog)
        self.uirodparaminCB=[self.uirodpara.qoffminCB,self.uirodpara.yscaleminCB,self.uirodpara.sizeminCB,self.uirodpara.sizeresminCB,self.uirodpara.bgconminCB,self.uirodpara.bglinminCB,self.uirodpara.roughminCB]
        self.uirodparamaxCB=[self.uirodpara.qoffmaxCB,self.uirodpara.yscalemaxCB,self.uirodpara.sizemaxCB,self.uirodpara.sizeresmaxCB,self.uirodpara.bgconmaxCB,self.uirodpara.bglinmaxCB,self.uirodpara.roughmaxCB]
        self.uirodparaminLE=[self.uirodpara.qoffminLE,self.uirodpara.yscaleminLE,self.uirodpara.sizeminLE,self.uirodpara.sizeresminLE,self.uirodpara.bgconminLE,self.uirodpara.bglinminLE,self.uirodpara.roughminLE]
        self.uirodparamaxLE=[self.uirodpara.qoffmaxLE,self.uirodpara.yscalemaxLE,self.uirodpara.sizemaxLE,self.uirodpara.sizeresmaxLE,self.uirodpara.bgconmaxLE,self.uirodpara.bglinmaxLE,self.uirodpara.roughmaxLE]
        for i in range(len(self.rodpara)):
            if self.rodpara[i][2]!=None:
                self.uirodparaminCB[i].setCheckState(2)
                self.uirodparaminLE[i].setText(str(self.rodpara[i][2]))
            if self.rodpara[i][3]!=None:
                self.uirodparamaxCB[i].setCheckState(2)
                self.uirodparamaxLE[i].setText(str(self.rodpara[i][3]))
        self.uirodpara.show()
        self.connect(self.uirodpara.cancelPB, SIGNAL('clicked()'), self.cancelRodPara)
        self.connect(self.uirodpara.okPB, SIGNAL('clicked()'), self.takeRodPara)
        
    def cancelRodPara(self):
        self.uirodpara.close()
                
    def takeRodPara(self):
        for i in range(len(self.rodpara)):
            if self.uirodparaminCB[i].checkState()!=0 and self.uirodparamaxCB[i].checkState()!=0 and float(self.uirodparaminLE[i].text())>float(self.uirodparamaxLE[i].text()):
                self.messageBox("Error:: Low constrain must be smaller than high constrain for "+ str(self.rodparaname[i])+"!!!")
                index=1                
                break
            else:
                index=0
                if self.uirodparaminCB[i].checkState()!=0:
                    self.rodpara[i][2]=float(self.uirodparaminLE[i].text())
                else:
                    self.rodpara[i][2]=None
                if self.uirodparamaxCB[i].checkState()!=0:
                    self.rodpara[i][3]=float(self.uirodparamaxLE[i].text())
                else: 
                    self.rodpara[i][3]=None
        if index==0:       
            self.uirodpara.close()        
    
    def saveRod(self):
        if str(self.ui.rodsaveCB.currentText())=='Save Fit':
            self.saveRodFitDig()
        elif str(self.ui.rodsaveCB.currentText())=='Save Para':
            self.saveRodPara()
    
    def saveRodPara(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Rod Scan Fitting Parameters',directory=self.directory))
        fid=open(self.saveFileName+'_par.txt','w')
        try:
            fid.write('Chi_Square\t'+format(self.rodresult.redchi, '.3f')+'\n')  #chisquare
        except:
            fid.write('Chi_Square\tNA\n') 
        if self.ui.rodlipidCB.checkState()!=0:
            fid.write('System_Type\tLipids\n')
        else:
            fid.write('System_Type\tNPs\n')
        fid.write('Error_Type\t'+str(self.ui.roderrCB.currentText()).split()[0]+'\n')
        fid.write('Para_Name\tValue\t\tVary\tStderr\t\tMin\tMax\n')
        for i in range(len(self.rodpara)):
            try:
                fid.write(self.rodparaname[i]+'   \t'+format(self.rodresult.params[self.rodparaname[i]].value,'.3e')+'\t'+str(self.rodpara[i][1])+'\t'+format(self.rodresult.params[self.rodparaname[i]].stderr,'.3e')+'\t'+str(self.rodpara[i][2])+'\t'+str(self.rodpara[i][3])+'\n')
            except:
                fid.write(self.rodparaname[i]+'   \t'+format(float(self.uirodLE[i].text()),'.3e')+'\tNA\tNA\t\tNA\tNA\n')
        fid.write('Constants:\n')
        fid.write('Xray_energy\t'+format(float(self.ui.rodxenLE.text()),'.3f')+'\n')
        fid.write('Rho_subphase\t'+format(float(self.ui.rodrhoLE.text()),'.3f')+'\n')
        fid.write('Angle_alpha\t'+format(float(self.ui.rodalphaLE.text()),'.4f')+'\n')
        fid.write('Angle_dth\t'+format(float(self.ui.roddthLE.text()),'.3f')+'\n')
        fid.close()    
    
    def loadRod(self):
        if str(self.ui.rodloadCB.currentText())=='Load Para':
            self.loadRodPara()
            
    def loadRodPara(self):
        filename=QFileDialog.getOpenFileName(caption='Select Parameter File to read', directory=self.directory, filter='Par Files (*.par*;*_par.txt)')
        self.directory=str(QFileInfo(filename).absolutePath())
        fid=open(filename)
        fdata=fid.readlines()
        fid.close()
        system=fdata[1][:-1].split('\t')[1]
        if system=='Lipids':
            self.ui.rodlipidCB.setCheckState(2)
            self.ui.rodNPCB.setCheckState(0)
        else:
            self.ui.rodlipidCB.setCheckState(0)
            self.ui.rodNPCB.setCheckState(2)
        para=[]
        for i in range(4,4+len(self.rodpara)):
            para.append(eval(fdata[i].split('\t')[1]))
        for i in range(len(self.rodpara)):
            self.uirodLE[i].setText(format(para[i],'.2e'))
        cons=[]
        for i in range(5+len(self.rodpara),len(fdata)):
            cons.append(eval(fdata[i].split('\t')[1]))
        self.ui.rodxenLE.setText(str(cons[0]))
        self.ui.rodrhoLE.setText(str(cons[1]))
        self.ui.rodalphaLE.setText(str(cons[2]))
        self.ui.roddthLE.setText(str(cons[3]))
        self.ui.calrodCB.setCheckState(2)
        self.updateRodCal()
    
    def saveRodFitDig(self):
        Dialog=QDialog(self)
        self.uirodsavefit=uic.loadUi('refsave.ui', Dialog)
        self.uirodsavefit.label.setText('Save Rod Scan Fit/Calcualtion!')
        try:
            self.uirodsavefit.xminLE.setText(str(self.rodxmin))
            self.uirodsavefit.xmaxLE.setText(str(self.rodxmax))
        except:
            pass
        self.uirodsavefit.numpointLE.setText(str(1600))
        self.uirodsavefit.show()
        self.connect(self.uirodsavefit.cancelPB, SIGNAL('clicked()'), self.cancelSaveRodFit)
        self.connect(self.uirodsavefit.okPB, SIGNAL('clicked()'), self.saveRodFit)
        
    def cancelSaveRodFit(self):
        self.uirodsavefit.close()
        self.rodsavefitindex=0
        
    def saveRodFit(self):
        if float(self.uirodsavefit.xminLE.text())>=float(self.uirodsavefit.xmaxLE.text()) or float(self.uirodsavefit.numpointLE.text())<=0:
            self.messageBox('Warning::Starting points must be lower than ending points \n and numer of points must be positive!!')
        else:
            self.rodsavefitindex=1  
            self.rodnp=float(self.uirodsavefit.numpointLE.text())
            self.rodxmin=float(self.uirodsavefit.xminLE.text())
            self.rodxmax=float(self.uirodsavefit.xmaxLE.text())
            self.updateRodCal()
            self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Rod Fit Data',directory=self.directory))
            fname=self.saveFileName+'_fit.txt'
            np.savetxt(fname,self.rodcal,fmt='%.4e\t%.4e')
            self.rodsavefitindex=0
            self.uirodsavefit.close()    
            
    def formfactorShow(self):
        Dialog=QDialog(self)
        ui=uic.loadUi('formfactorDialog.ui',Dialog)
        ui.show()
        cylffPixmap=QtGui.QPixmap('cylff.png')
        cylffscaledPixmap=cylffPixmap.scaled(ui.cylLabel.size(),Qt.KeepAspectRatio)
        ui.cylLabel.setPixmap(cylffscaledPixmap)
        sphffPixmap=QtGui.QPixmap('sphff.png')
        sphffscaledPixmap=sphffPixmap.scaled(ui.sphLabel.size(),Qt.KeepAspectRatio)
        ui.sphLabel.setPixmap(sphffscaledPixmap)
        Dialog.exec_()  
#        app = QApplication(sys.argv)
#        web = QWebView()
#        web.settings().setAttribute(QWebSettings.PluginsEnabled, True)
#        web.show()
#        web.load(QUrl('formfactor.pdf')) 
#        sys.exit(app.exec_())    
    
        ################################################        
        #start the fluorescence analysis section. 
        ################################################

    def initFluPar(self):  #initialize the flu parameters
        self.ui.flusubTW.horizontalHeader().setVisible(True)
        self.ui.flusubTW.verticalHeader().setVisible(True)
        self.ui.flusubTW.setHorizontalHeaderLabels(QStringList()<<'Element'<<'Composition'<<'Ionic Radius'+' ('+u'\u212b'+')')
        self.fluparaname=['sur_den','q_off','y_scale','bg_con','bg_lin','sur_cur']
        self.uifluLE=[self.ui.flusurLE, self.ui.fluqoffLE, self.ui.fluyscaleLE, self.ui.fluconLE, self.ui.flulinLE, self.ui.flusurcurLE]
        self.uifluCB=[self.ui.flusurCB, self.ui.fluqoffCB, self.ui.fluyscaleCB, self.ui.fluconCB, self.ui.flulinCB, self.ui.flusurcurCB]
        self.fluconsname=['Con_bulk', 'Xray_energy', 'Fluo_energy','Rho_topphase', 'Rho_botphase','Beta_topphase','Beta_bot(inc)','Beta_bot(flu)','Slit_vertical', 'Detector_len']  
        self.uifluconLE=[self.ui.flubulLE, self.ui.fluxenLE, self.ui.flufluenLE, self.ui.flurhotopLE, self.ui.flurhobotLE, self.ui.flubetatopLE, self.ui.flubetabotLE, self.ui.flubetabot2LE, self.ui.flusliLE, self.ui.fludetLE]
        self.flupara={}  #initialize the flu  parameter dictonary
        for i in range(len(self.fluparaname)):
            self.flupara[i]=[float(self.uifluLE[i].text()), False, None, None]  
        self.updateFluElement()
        self.updateFluQc()
            
    def openFluFile(self):  #open flu files and also remove all current ref files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Fluorescence Files to import', directory=self.directory, filter='Flu Files (*.flu*;*_flu.txt)')
        self.ui.tabWidget.setCurrentIndex(2)
        self.flufiles=map(str, f)
        self.directory=str(QFileInfo(self.flufiles[0]).absolutePath())
        self.updateFluFile()
    
    def updateFluFile(self): #update flu files in the listwidget
        self.ui.flufileLW.clear()
        for i in range(len(self.flufiles)):
            try:
                self.ui.flufileLW.addItem('#'+str(i+1)+self.halftab+str(self.flufiles[i].split('\\')[-2])+'\\'+str(self.flufiles[i].split('\\')[-1]))
            except:
                self.ui.flufileLW.addItem('#'+str(i+1)+self.halftab+str(self.flufiles[i].split('/')[-2])+'/'+str(self.flufiles[i].split('/')[-1]))

    def addFluFile(self): #add flu files into the listwidget and deselect all flu files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Fluorescence Files to import', directory=self.directory, filter='Flu Files (*.flu*;*_flu.txt)')
        self.flufiles=self.flufiles+map(str, f)
        self.directory=str(QFileInfo(self.flufiles[0]).absolutePath())
        self.updateFluFile()
        
    def updateSelectedFluFile(self): #update the selected flu files in the listwidget
        selectedflufiles=self.ui.flufileLW.selectedItems()
        self.selectedflufiles_rows=[]
        for item in selectedflufiles:
            self.selectedflufiles_rows.append(self.ui.flufileLW.row(item))
        self.selectedflufiles_rows.sort()
        self.fluscale=[[1,0,1,0] for i in range(len(self.selectedflufiles_rows))]
        self.updateFluPlot()
        
    def removeFluFile(self): #remove flu files in the listwidget and deselect all flu files in the listwidget
        items=self.ui.flufileLW.selectedItems()
        for item in items:
            self.flufiles.pop(self.ui.flufileLW.row(item))
        self.ui.flufileLW.clear()
        self.updateFluFile()
           
    def updateFluFitFile(self): #update flu fit files in the listwidget
        self.ui.flufitfileLW.clear()
        for i in range(len(self.flufitfiles)):
            try:
                self.ui.flufitfileLW.addItem('#'+str(i+1)+self.halftab+str(self.flufitfiles[i].split('\\')[-2])+'\\'+str(self.flufitfiles[i].split('\\')[-1]))
            except:
                self.ui.flufitfileLW.addItem('#'+str(i+1)+self.halftab+str(self.flufitfiles[i].split('/')[-2])+'/'+str(self.flufitfiles[i].split('/')[-1]))

    def addFluFitFile(self): #add flu fit files into the listwidget and deselect flu fit files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Fluorescence Fit Files to import', directory=self.directory, filter='FIT Files (*.fit*; *_fit.txt)')
        self.flufitfiles=self.flufitfiles+map(str, f)
        self.directory=str(QFileInfo(self.flufitfiles[0]).absolutePath())
        self.updateFluFitFile()
        
    def updateSelectedFluFitFile(self): #update the selected flu fit files in the listwidget
        selectedflufitfiles=self.ui.flufitfileLW.selectedItems()
        self.selectedflufitfiles_rows=[]
        for item in selectedflufitfiles:
            self.selectedflufitfiles_rows.append(self.ui.flufitfileLW.row(item))
        self.selectedflufitfiles_rows.sort()
        self.flufitscale=[[1,0,1,0] for i in range(len(self.selectedflufitfiles_rows))]
        self.updateFluPlot()
        
    def removeFluFitFile(self):  #remove flu fit files in the listwidget and deselect all flu fit files in the listwidget
        items=self.ui.flufitfileLW.selectedItems()
        for item in items:
            self.flufitfiles.pop(self.ui.flufitfileLW.row(item))
        self.ui.flufitfileLW.clear()
        self.updateFluFitFile()
        
    
    def updateFluPlot(self): #update the plot in the flu plotwidget
        self.ui.fluPW.canvas.ax.clear()
        self.ui.fluPW.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.ui.fluPW.canvas.ax.set_ylabel('Intensity [a.u.]')
        if  len(self.selectedflufiles_rows)!=0: #plot flu files
            for i in range(len(self.selectedflufiles_rows)):
                data1=np.loadtxt(str(self.flufiles[self.selectedflufiles_rows[i]]), comments='#')
                self.ui.fluPW.canvas.ax.errorbar(data1[:,0]*self.fluscale[i][0]+self.fluscale[i][1],data1[:,1]*self.fluscale[i][2]+self.fluscale[i][3],data1[:,2]*self.fluscale[i][2],fmt='o',label='#'+str(self.selectedflufiles_rows[i]+1))
        if  len(self.selectedflufitfiles_rows)!=0: #plot flu fit files
            for i in range(len(self.selectedflufitfiles_rows)):
                data1=np.loadtxt(str(self.flufitfiles[self.selectedflufitfiles_rows[i]]), comments='#')
                self.ui.fluPW.canvas.ax.errorbar(data1[:,0]*self.flufitscale[i][0]+self.flufitscale[i][1],data1[:,1]*self.flufitscale[i][2]+self.flufitscale[i][3],fmt='-',label='#'+str(self.selectedflufitfiles_rows[i]+1))
        if  self.ui.calfluCB.checkState()!=0:
                self.ui.fluPW.canvas.ax.errorbar(np.array(self.flucal)[:,0],np.array(self.flucal)[:,1],fmt='-', label='cal')
        if self.ui.flushoqcCB.checkState()!=0:
            self.ui.fluPW.canvas.ax.axvline(x=self.fluqc-float(self.ui.fluqoffLE.text()), linestyle='dashed', label='qc', color='k')
        if self.ui.flulegendCB.checkState()!=0:
            self.ui.fluPW.canvas.ax.legend(loc=self.ui.flulegendlocCoB.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        if self.ui.flulogyCB.checkState()!=0:
            self.ui.fluPW.canvas.ax.set_yscale('log')
        else:
            self.ui.fluPW.canvas.ax.set_yscale('linear')
        self.ui.fluPW.canvas.draw()
        
    def setFluPlotScale(self): #set the scale of each data in the flu plot 
        if len(self.selectedflufiles_rows)+len(self.selectedflufitfiles_rows)==0:
            self.messageBox('Warning:: No Fluorescence or Fit files selected!')
        else:
            row_flu=len(self.selectedflufiles_rows)
            row_fit=len(self.selectedflufitfiles_rows)
            row=row_flu+row_fit
            Dialog=QDialog(self)
            self.uiplotscale=uic.loadUi('plotscale.ui', Dialog)
            self.uiplotscale.scaleTW.setRowCount(row) #set the table size; 4 column is fixed
            self.uiplotscale.show()
            self.uiplotscale.scaleLabel.setText('Fluorescence Plot Scale Setup: X=X*Factor+Offset')
            self.uiplotscale.scaleTW.setHorizontalHeaderLabels(QStringList()<<"X Factor"<<"X Offset"<<"Y Factor"<<"Y Offset") #set the horizontal header
            vlabel=QStringList() #set the vertical header 
            for i in range(row_flu):
                vlabel.append("Flu #"+str(self.selectedflufiles_rows[i]+1))
            for i in range(row_fit):
                vlabel.append("Fit #"+str(self.selectedflufitfiles_rows[i]+1))
            self.uiplotscale.scaleTW.setVerticalHeaderLabels(vlabel)
            for i in range(row_flu):  #set the initial values
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i,j,QTableWidgetItem(str(self.fluscale[i][j])))
                    self.uiplotscale.scaleTW.item(i,j).setTextAlignment(Qt.AlignCenter)
            for i in range(row_fit):
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i+row_flu,j,QTableWidgetItem(str(self.flufitscale[i][j])))
                    self.uiplotscale.scaleTW.item(i+row_flu,j).setTextAlignment(Qt.AlignCenter)
            self.connect(self.uiplotscale.scaleTW, SIGNAL('cellChanged(int,int)'), self.updateFluPlotScale) #update the flu scale and plot
            self.connect(self.uiplotscale.closePB,SIGNAL('clicked()'), self.closePlotScale) #close the scale setup window
           
                       
    def updateFluPlotScale(self): #update the scale of each data in the flu plot
        row_flu=len(self.selectedflufiles_rows)
        row_fit=len(self.selectedflufitfiles_rows)
        self.fluscale=[[float(str(self.uiplotscale.scaleTW.item(i,j).text())) for j in range(4)] for i in range(row_flu)]
        self.flufitscale=[[float(str(self.uiplotscale.scaleTW.item(i+row_flu,j).text())) for j in range(4)] for i in range(row_fit)]
        self.updateFluPlot()

    def updateFluElement(self): #update the subphase ion info including xray propeties
        xrayen=float(self.ui.fluxenLE.text())
        k0=2*np.pi*xrayen/12.3984 # wave vector
        fluen=float(self.ui.flufluenLE.text())  
        k1=2*np.pi*fluen/12.3984 # wave vector for emission line
        row=self.ui.flusubTW.rowCount()
        con_bulk=float(self.ui.flubulLE.text())   #get the bulk concentration
        volume=0  #volume for ions for 1 L subphase
        numele=0  # total number of electrons from 1 L subphase
        beta=0  #beta for the subphase at incident x-ray energy
        beta1=0 #beta for the subphase at emission x-ray energy
        self.fluelepara={}  # setup dict for all elements in the subphase  
        for i in range(row):
            element=str(self.ui.flusubTW.item(i,0).text())  #get the element for this row
            try:
                self.fluelepara[i]=[element,float(str(self.ui.flusubTW.item(i,1).text())),float(str(self.ui.flusubTW.item(i,2).text())),elements.symbol(element).number, elements.symbol(element).xray.scattering_factors(energy=xrayen)[1], elements.symbol(element).xray.scattering_factors(energy=fluen)[1]]
                volume=volume+con_bulk*self.fluelepara[i][1]*4/3*np.pi*self.fluelepara[i][2]**3*self.avoganum/1e27
                numele=numele+con_bulk*self.fluelepara[i][1]*self.fluelepara[i][3]*self.avoganum
                beta=beta+2*np.pi*self.eleradius*con_bulk*self.fluelepara[i][1]*self.avoganum*self.fluelepara[i][4]/k0/k0/1e27
                beta1=beta1+2*np.pi*self.eleradius*con_bulk*self.fluelepara[i][1]*self.avoganum*self.fluelepara[i][5]/k1/k1/1e27
            except:
                self.messageBox('Error: unknown element ' + element+'!')
                break
        self.flutopdel=self.eleradius*2*np.pi/k0/k0*float(self.ui.flurhotopLE.text())  
        self.flutopbet=float(self.ui.flubetatopLE.text()) 
        flubotbeta=float(self.ui.flubetabotLE.text())   # beta of water for incident beam 
        flubotbeta2=float(self.ui.flubetabot2LE.text())   # beta of water for flurescenct beam 
        self.botrho=(1-volume)*float(self.ui.flurhobotLE.text())+numele/1e27
        self.flubotdel=self.eleradius*2*np.pi/k0/k0*self.botrho
      #  print beta, beta1
        self.flubotbeta=beta+(1-volume)*flubotbeta  #beta =3.462e-10 for water at 20keV
        self.flubotmu1=2*k1*(beta1+(1-volume)*flubotbeta2)  #beta= 1.24492e-9 for water at 14.148keV; mu for the emission line
        self.fluqc=2*np.sqrt(2)*k0*np.sqrt(self.flubotdel-self.flutopdel)  #get qc
    #    print self.fluqc, self.toprho
      #  print volume, toprho            
      #  print self.fluelepara

    def updateFluQc(self):
        self.updateFluElement()
        self.ui.flurelrhobotLE.setText(format(self.botrho, '.4f'))
        self.ui.fluqcLE.setText(format(self.fluqc, '.4f'))


    def insFluIon(self):  # add one ion in the subphase
        insrows=self.ui.flusubTW.selectionModel().selectedRows()
        insrows=[self.ui.flusubTW.row(self.ui.flusubTW.itemFromIndex(insrows[i])) for i in range(len(insrows))]
        if len(insrows)!=1:
            self.messageBox('Warning:: Only one row can be seleted!')
        else:
            self.ui.flusubTW.insertRow(insrows[0])
        for i in range(3):
            self.ui.flusubTW.setItem(insrows[0],i,QTableWidgetItem('Cl/2/1.80'.split('/')[i]))
            
    def rmFluIon(self): #remove one ion in the subphase
        rmrows=self.ui.flusubTW.selectionModel().selectedRows()   
        removerows=[]
        for rmrow in rmrows:
            removerows.append(self.ui.flusubTW.row(self.ui.flusubTW.itemFromIndex(rmrow)))
            removerows.sort(reverse=True)
        if len(removerows)==0:
            self.messageBox('Warning:: No ion is selected!!')
        else:
            for i in range(len(removerows)):
                self.ui.flusubTW.removeRow(removerows[i])
    
    
    def updateFluParVal(self): #update the flu parameters value
        for i in range(len(self.fluparaname)):
            self.flupara[i][0]=float(self.uifluLE[i].text())   
    
    def updateFluCal(self): # caluate the flu  based on current parameters.
        self.updateFluParVal()
        self.updateFluElement()
        flupara=[self.flupara[i][0] for i in range(len(self.flupara))]
        if self.flusavefitindex==1:
            xflu=np.linspace(max(0,self.fluxmin),self.fluxmax,self.flunp)
            self.flucal=np.vstack((xflu,self.fluCalFun(flupara,xflu))).T
        else:
            if  self.ui.calfluCB.checkState()!=0:
                if  len(self.selectedflufiles_rows)!=0: 
                    data=np.loadtxt(str(self.flufiles[self.selectedflufiles_rows[0]]), comments='#')
                    self.fluxmax=np.max(data[:,0])
                    self.fluxmin=max(0,np.min(data[:,0]))
                else:
                    self.fluxmax=self.fluqc+0.006  #only calculate the flu around qc (+/- 0.006)
                    self.fluxmin=self.fluqc-0.006
                xflu=np.linspace(self.fluxmin,self.fluxmax,200)
                self.flucal=np.vstack((xflu,self.fluCalFun(flupara,xflu))).T
            self.updateFluPlot()
            
    def fluCalFun(self,flupara,x):
        surden=flupara[0]
        qoff=flupara[1]
        yscale=flupara[2]
        bgcon=flupara[3]
        bglin=flupara[4]
        surcur=flupara[5]*1e10   # in unit of /AA
        conbulk=float(self.ui.flubulLE.text())   #get the bulk concentration
        k0=2*np.pi*float(self.ui.fluxenLE.text())/12.3984 # wave vector
        slit=float(self.ui.flusliLE.text())  #get slits size
        detlen=float(self.ui.fludetLE.text())*1e7  #get slits size in unit of /AA 
        topd=1/(self.flutopbet*2*k0) #get the absorption length in top phase
        qz=x+qoff
        alpha=qz/2/k0  #get incident angle 
        fprint=slit/alpha*1e7 #get the footprint in unit of /AA
        if surcur==0:  #no surface curvature  
            flu=[]
            #p_d=[]
            for i in range(len(alpha)):
                z1=(fprint[i]-detlen)/2*alpha[i]
                z2=(fprint[i]+detlen)/2*alpha[i]
                effd,trans=self.frsnllCal(self.flutopdel,self.flutopbet,self.flubotdel,self.flubotbeta,self.flubotmu1,k0,alpha[i]) 
                effv=effd*topd*np.exp(-detlen/2/topd)*(detlen*effd*np.exp(z2/alpha[i]/topd)*(np.exp(-z1/effd)-np.exp(-z2/effd))+topd*(np.exp(detlen/topd)-1)*(z1-z2))/(detlen*effd+topd*(z1-z2))
                int_sur=surden*topd*(np.exp(detlen/2/topd)-np.exp(-detlen/2/topd))   #surface intensity
                int_bulk=effv*self.avoganum*conbulk*self.fluelepara[0][1]/1e27  #bluk intensity; the element in the first row is the target element
                int_tot=yscale*trans*(int_sur+int_bulk)+bgcon+bglin*alpha[i]
                flu.append(int_tot)
             #   p_d.append(effd)
        else:  #with surface curvature 
            flu=[]
            for i in range(len(alpha)):
                bsum=0
                ssum=0
                steps=int((detlen+fprint[i])/2/1e6)  # use 0.1 mm as the step size
                stepsize=(detlen+fprint[i])/2/steps
                x=np.linspace(-fprint[i]/2, detlen/2, steps) # get the position fo single ran hitting the surface relative to the center of detector area with the step size "steps"
                for j in range(len(x)):
                    alphanew=alpha[i]-x[j]/surcur  # the incident angle at position x[j]
                    y1=-detlen/2-x[j]
                    y2=detlen/2-x[j]
                    effd,trans=self.frsnllCal(self.flutopdel,self.flutopbet,self.flubotdel,self.flubotbeta,self.flubotmu1,k0,alphanew) 
                    if x[j]>-detlen/2:
                        bsum=bsum+np.exp(-x[j]/topd)*trans*effd*(1.0-np.exp(-y2*alpha[i]/effd))
                        ssum=ssum+np.exp(-x[j]/topd)*trans
                    else:
                        bsum=bsum+np.exp(-x[j]/topd)*trans*effd*(np.exp(-y1*alpha[i]/effd)-np.exp(-y2*alpha[i]/effd))  #surface has no contribution at this region
                int_bulk=bsum*stepsize*self.avoganum*conbulk*self.fluelepara[0][1]/1e27
                int_sur=ssum*stepsize*surden
                int_tot=yscale*(int_bulk+int_sur)+bgcon+bglin*alpha[i]
                flu.append(int_tot)
        return flu
        
    def frsnllCal(self, dett,bett,detb,betb,mub,k0,alpha):
        f1=cmath.sqrt(complex(alpha*alpha,2*bett))
        fmax=cmath.sqrt(complex(alpha*alpha-2*(detb-dett),2*betb))
        length1=1/mub
        length2=1/(2*k0*fmax.imag)
        eff_d=length1*length2/(length1+length2)
        trans=4*abs(f1/(f1+fmax))*abs(f1/(f1+fmax))
       # frsnll=abs((f1-fmax)/(f1+fmax))*abs((f1-fmax)/(f1+fmax))
        return eff_d,trans
        
    def fitFlu(self):
        self.updateFluParVal()
        self.updateFluElement()
        for i in range(len(self.flupara)):
            if self.uifluCB[i].checkState()!=0:  #set the selected parameters to be varied
                self.flupara[i][1]=True
            else:
                self.flupara[i][1]=False
        parastatus=np.array([self.flupara[i][1] for i in range(len(self.flupara))])
        selparas=np.where(parastatus==True)
        self.fluparameter=Parameters()
        for i in range(len(self.flupara)):
            self.fluparameter.add(self.fluparaname[i], value=self.flupara[i][0],vary=self.flupara[i][1],min=self.flupara[i][2],max=self.flupara[i][3])
        if len(self.selectedflufiles_rows)!=1:
            self.messageBox('Please select only one set of data for fitting!')
        else:
            data=np.loadtxt(str(self.flufiles[self.selectedflufiles_rows[0]]), comments='#')
            ini=max(float(str(self.ui.flufitranLE.text()).split(':')[0]),np.min(data[:,0]))
            fin=min(float(str(self.ui.flufitranLE.text()).split(':')[1]),np.max(data[:,0]))
            data1=data[np.where(np.logical_and(data[:,0]>=ini,data[:,0]<=fin))]
            x=data1[:,0]
            y=data1[:,1]
            if self.ui.fluerrCB.currentIndex()==0:
                yerr=data1[:,2]
            elif self.ui.fluerrCB.currentIndex()==1:
                yerr=np.sqrt(y)
            elif self.ui.fluerrCB.currentIndex()==2:
                yerr=y
            else:
                yerr=np.ones_like(x)
            self.fluresult=minimize(self.flu2min, self.fluparameter, args=(x,y,yerr))
            print(fit_report(self.fluresult))
            residual=np.vstack((x,self.fluresult.residual)).T
            self.ui.flusurLE.setText(format(self.fluresult.params[self.fluparaname[0]].value, '.5f'))
            self.ui.fluqoffLE.setText(format(self.fluresult.params[self.fluparaname[1]].value, '.6f'))
            self.ui.fluyscaleLE.setText(format(self.fluresult.params[self.fluparaname[2]].value, '.2e'))
            self.ui.fluconLE.setText(format(self.fluresult.params[self.fluparaname[3]].value, '.2e'))
            self.ui.flulinLE.setText(format(self.fluresult.params[self.fluparaname[4]].value, '.2e'))
            self.ui.flusurcurLE.setText(format(self.fluresult.params[self.fluparaname[5]].value, '.2f'))
            self.ui.calfluCB.setCheckState(2)
            self.updateFluCal()
            self.ui.fluChiLE.setText(format(self.fluresult.redchi, '.3f'))
            self.ui.fluparaTB.clear()
            fitinfo='Fitting Paramenters:\n'
            fitinfo=fitinfo+'Name\tStderr\tMin\tMax\n'
            for i in selparas[0]:
                fitinfo=fitinfo+self.fluparaname[i]+'\t'+format(self.fluresult.params[self.fluparaname[i]].stderr, '.3e')+'\t'+str(self.flupara[i][2])+'\t'+str(self.flupara[i][3])+'\n'
            fitinfo=fitinfo+'********************************\n'
            fitinfo=fitinfo+'Fitting Residual:\n'
            for i in range(len(residual)):
                fitinfo=fitinfo+format(residual[i][0], '.3f')+'\t'+format(residual[i][1], '.4f')+'\n'
            self.ui.fluparaTB.append(fitinfo)
            cursor=self.ui.fluparaTB.textCursor()
            cursor.setPosition(0)
            self.ui.fluparaTB.setTextCursor(cursor)
            
    
    def flu2min(self,params,x,y,yerr): # residuel for flu fitting
        flupara=[params[self.fluparaname[i]] for i in range(len(params))]
        model=self.fluCalFun(flupara,x)
        return (model-y)/yerr
        
    def updateFluPara(self):
        Dialog=QDialog(self)
        self.uiflupara=uic.loadUi('flupara.ui', Dialog)
        self.uifluparaminCB=[self.uiflupara.surdenminCB,self.uiflupara.qoffminCB,self.uiflupara.yscaleminCB,self.uiflupara.bgconminCB,self.uiflupara.bglinminCB,self.uiflupara.surcurminCB]
        self.uifluparamaxCB=[self.uiflupara.surdenmaxCB,self.uiflupara.qoffmaxCB,self.uiflupara.yscalemaxCB,self.uiflupara.bgconmaxCB,self.uiflupara.bglinmaxCB,self.uiflupara.surcurmaxCB]
        self.uifluparaminLE=[self.uiflupara.surdenminLE,self.uiflupara.qoffminLE,self.uiflupara.yscaleminLE,self.uiflupara.bgconminLE,self.uiflupara.bglinminLE,self.uiflupara.surcurminLE]
        self.uifluparamaxLE=[self.uiflupara.surdenmaxLE,self.uiflupara.qoffmaxLE,self.uiflupara.yscalemaxLE,self.uiflupara.bgconmaxLE,self.uiflupara.bglinmaxLE,self.uiflupara.surcurmaxLE]
        for i in range(len(self.flupara)):
            if self.flupara[i][2]!=None:
                self.uifluparaminCB[i].setCheckState(2)
                self.uifluparaminLE[i].setText(str(self.flupara[i][2]))
            if self.flupara[i][3]!=None:
                self.uifluparamaxCB[i].setCheckState(2)
                self.uifluparamaxLE[i].setText(str(self.flupara[i][3]))
        self.uiflupara.show()
        self.connect(self.uiflupara.cancelPB, SIGNAL('clicked()'), self.cancelFluPara)
        self.connect(self.uiflupara.okPB, SIGNAL('clicked()'), self.takeFluPara)
        
    def cancelFluPara(self):
        self.uiflupara.close()
                
    def takeFluPara(self):
        for i in range(len(self.flupara)):
            if self.uifluparaminCB[i].checkState()!=0 and self.uifluparamaxCB[i].checkState()!=0 and float(self.uifluparaminLE[i].text())>float(self.uifluparamaxLE[i].text()):
                self.messageBox("Error:: Low constrain must be smaller than high constrain for "+ str(self.fluparaname[i])+"!!!")
                index=1                
                break
            else:
                index=0
                if self.uifluparaminCB[i].checkState()!=0:
                    self.flupara[i][2]=float(self.uifluparaminLE[i].text())
                else:
                    self.flupara[i][2]=None
                if self.uifluparamaxCB[i].checkState()!=0:
                    self.flupara[i][3]=float(self.uifluparamaxLE[i].text())
                else: 
                    self.flupara[i][3]=None
        if index==0:       
            self.uiflupara.close()   
            
    def saveFlu(self):
        if str(self.ui.flusaveCB.currentText())=='Save Fit':
            self.saveFluFitDig()
        elif str(self.ui.flusaveCB.currentText())=='Save Para':
            self.saveFluPara()
            
    def saveFluPara(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Fluorescence Fitting Parameters',directory=self.directory))
        fid=open(self.saveFileName+'_par.txt','w')
        try:
            fid.write('Chi_Square\t'+format(self.fluresult.redchi, '.3f')+'\n')  #chisquare
        except:
            fid.write('Chi_Square\tNA\n') 
        fid.write('Error_Type\t'+str(self.ui.fluerrCB.currentText()).split()[0]+'\n')
        fid.write('Para_Name\tValue\t\tVary\tStderr\t\tMin\tMax\n')
        for i in range(len(self.flupara)):
            try:
                fid.write(self.fluparaname[i]+'   \t'+format(self.fluresult.params[self.fluparaname[i]].value,'.3e')+'\t'+str(self.flupara[i][1])+'\t'+format(self.fluresult.params[self.fluparaname[i]].stderr,'.3e')+'\t'+str(self.flupara[i][2])+'\t'+str(self.flupara[i][3])+'\n')
            except:
                fid.write(self.fluparaname[i]+'   \t'+format(float(self.uifluLE[i].text()),'.3e')+'\tNA\tNA\t\tNA\tNA\n')
        fid.write('Constants:\n')
        for i in range(len(self.fluconsname)):
            fid.write(str(self.fluconsname[i])+'\t'+format(float(self.uifluconLE[i].text()),'.3e')+'\n')
        fid.write('Ele.\tComp.\tRad.\n')
        row=self.ui.flusubTW.rowCount()
        for i in range(row):
            fid.write(str(self.ui.flusubTW.item(i,0).text())+'\t'+str(self.ui.flusubTW.item(i,1).text())+'\t'+str(self.ui.flusubTW.item(i,2).text())+'\n')
        fid.close()  
    
    def loadFlu(self):
        if str(self.ui.fluloadCB.currentText())=='Load Para':
            self.loadFluPara()
            
    def loadFluPara(self):
        filename=QFileDialog.getOpenFileName(caption='Select Parameter File to read', directory=self.directory, filter='Par Files (*.par*;*_par.txt)')
        self.directory=str(QFileInfo(filename).absolutePath())
        fid=open(filename)
        fdata=fid.readlines()
        fid.close()
        self.ui.calfluCB.setCheckState(0)
        para=[]
        for i in range(3,3+len(self.flupara)):
            para.append(eval(fdata[i].split('\t')[1]))
        for i in range(len(self.flupara)):
            self.uifluLE[i].setText(format(para[i],'.3e'))
        cons=[]
        for i in range(4+len(self.flupara),4+len(self.flupara)+len(self.fluconsname)):
            cons.append(eval(fdata[i].split('\t')[1]))
        for i in range(len(self.fluconsname)):
            self.uifluconLE[i].setText(format(cons[i],'.3e'))
        elements={}
        for i in range(5+len(self.flupara)+len(self.fluconsname), len(fdata)):
            elements[i-(5+len(self.flupara)+len(self.fluconsname))]=fdata[i].split('\n')[0].split('\t')
        self.ui.flusubTW.setRowCount(len(elements))
        for i in range(len(elements)):
            for j in range(3):
                self.ui.flusubTW.setItem(i,j,QTableWidgetItem(elements[i][j]))
        self.ui.calfluCB.setCheckState(2)
        self.updateFluCal()
    
    def saveFluFitDig(self):
        Dialog=QDialog(self)
        self.uiflusavefit=uic.loadUi('refsave.ui', Dialog)
        self.uiflusavefit.label.setText('Save Fluorescence Fit/Calcualtion!')
        try:
            self.uiflusavefit.xminLE.setText(str(self.fluxmin))
            self.uiflusavefit.xmaxLE.setText(str(self.fluxmax))
        except:
            pass
        self.uiflusavefit.numpointLE.setText(str(200))
        self.uiflusavefit.show()
        self.connect(self.uiflusavefit.cancelPB, SIGNAL('clicked()'), self.cancelSaveFluFit)
        self.connect(self.uiflusavefit.okPB, SIGNAL('clicked()'), self.saveFluFit)
        
    def cancelSaveFluFit(self):
        self.uiflusavefit.close()
        self.flusavefitindex=0
        
    def saveFluFit(self):
        if float(self.uiflusavefit.xminLE.text())>=float(self.uiflusavefit.xmaxLE.text()) or float(self.uiflusavefit.numpointLE.text())<=0:
            self.messageBox('Warning::Starting points must be lower than ending points \n and numer of points must be large than 0!!')
        else:
            self.flusavefitindex=1  
            self.flunp=float(self.uiflusavefit.numpointLE.text())
            self.fluxmin=float(self.uiflusavefit.xminLE.text())
            self.fluxmax=float(self.uiflusavefit.xmaxLE.text())
            self.updateFluCal()
            self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Fluorescence Fit Data',directory=self.directory))
            fname=self.saveFileName+'_fit.txt'
            np.savetxt(fname,self.flucal,fmt='%.4e\t%.4e')
            self.flusavefitindex=0
            self.uiflusavefit.close()    
            
            
    
        ################################################        
        #start the GIXOS analysis section. 
        ################################################ 
   
    
    def initGixPar(self): #initialize the gixpar table
        self.ui.gixparTW.horizontalHeader().setVisible(True)
        self.ui.gixparTW.verticalHeader().setVisible(True)
        self.ui.gixparTW.setHorizontalHeaderLabels(QStringList()<<'d ('+u'\u212b'+')'<<u'\u03c1'+' (e/'+u'\u212b'+u'\u00b3'+')'<<u'\u03bc'+' (cm'+u'\u207b'+u'\u00b9'+')'<<u'\u03c3'+' ('+u'\u212b'+')')
        top='top/0/0/3'
        bottom='bottom/0.333/0/NA'
        for i in range(4):
            self.ui.gixparTW.setItem(0,i,QTableWidgetItem(top.split('/')[i]))
            self.ui.gixparTW.setItem(1,i,QTableWidgetItem(bottom.split('/')[i]))
        self.ui.gixnumslabSB.setValue(0)
        self.gixpara={}  #initialize the parameter dictionary
        self.gixpara[0]=[0,False, None,None]
        self.gixpara[1]=[0,False, None,None]
        self.gixpara[2]=[3,False, None,None]
        self.gixpara[3]=[0.333,False, None,None]
        self.gixpara[4]=[0,False, None,None]
        self.gixsyspara={}  #initialize the gix system parameter dictonary
        self.gixsyspara[0]=[float(self.ui.gixqoffLE.text()), False, None, None]
        self.gixsyspara[1]=[float(self.ui.gixyscaleLE.text()), False, None, None]
        self.gixsyspara[2]=[float(self.ui.gixqmaxLE.text()), False, None, None]
        self.gixsysCB=[self.ui.gixqoffCB, self.ui.gixyscaleCB, self.ui.gixqmaxCB]
        self.gixsysLE=[self.ui.gixqoffLE, self.ui.gixyscaleLE, self.ui.gixqmaxLE]
        self.gixconLE=[self.ui.gixxenLE, self.ui.gixalphaLE, self.ui.gixdthLE, self.ui.gixtemLE, self.ui.gixtenLE]
        self.gixconsname=['Xray_energy','Angle_alpha', 'Angle_dth', 'Temperature', 'Surface_Tension']
        self.updateGixParaName()
        
    def updateGixParaName(self):
        top=['rho_t','mu_t','sigma_0']
        middle=[]
        bottom=['rho_b','mu_b']
        for i in range(self.ui.gixparTW.rowCount()-2):
            layer=str(i+1)
            middle.extend(['d'+layer,'rho'+layer,'mu'+layer,'sigma'+layer])
        self.gixparaname=top+middle+bottom
        self.gixsysparaname=['q_off', 'y_scale', 'q_max']
    
    def openGixFile(self):  #open gix files and also remove all current ref files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple GIXOS Files to import', directory=self.directory, filter='GIXOS Files (*gix*;*gio*)')
        self.ui.tabWidget.setCurrentIndex(3)
        self.gixfiles=map(str, f)
        self.directory=str(QFileInfo(self.gixfiles[0]).absolutePath())
        self.updateGixFile()
            
    def updateGixFile(self): #update gixos files in the listwidget
        self.ui.gixfileLW.clear()
        for i in range(len(self.gixfiles)):
            try:
                self.ui.gixfileLW.addItem('#'+str(i+1)+self.halftab+str(self.gixfiles[i].split('\\')[-2])+'\\'+str(self.gixfiles[i].split('\\')[-1]))
            except:
                self.ui.gixfileLW.addItem('#'+str(i+1)+self.halftab+str(self.gixfiles[i].split('/')[-2])+'/'+str(self.gixfiles[i].split('/')[-1]))

            
    def addGixFile(self): #add gixos files into the listwidget and deselect all gixos files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple GIXOS Files to import', directory=self.directory, filter='GIXOS Files (*gix*;*gio*)')
        self.gixfiles=self.gixfiles+map(str, f)
        self.directory=str(QFileInfo(self.gixfiles[0]).absolutePath())
        self.updateGixFile()
    
    def updateSelectedGixFile(self): #update the selected gixos files in the listwidget
        selectedgixfiles=self.ui.gixfileLW.selectedItems()
        self.selectedgixfiles_rows=[]
        for item in selectedgixfiles:
            self.selectedgixfiles_rows.append(self.ui.gixfileLW.row(item))
        self.selectedgixfiles_rows.sort()
        self.gixscale=[[1,0,1,0] for i in range(len(self.selectedgixfiles_rows))]
        self.updateGixPlot()
        
    def removeGixFile(self): #remove gixos files in the listwidget and deselect all gixos files in the listwidget
        items=self.ui.gixfileLW.selectedItems()
        for item in items:
            self.gixfiles.pop(self.ui.gixfileLW.row(item))
        self.ui.gixfileLW.clear()
        self.updateGixFile()  
        
    def updateGixPlot(self): #update the plot in the gixos plotwidget
        self.ui.gixPW.canvas.ax.clear()
        self.ui.gixPW.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.ui.gixPW.canvas.ax.set_ylabel('Intensity [a.u.]')
        if  len(self.selectedgixfiles_rows)!=0: #plot ref files
            for i in range(len(self.selectedgixfiles_rows)):
                data1=np.loadtxt(str(self.gixfiles[self.selectedgixfiles_rows[i]]), comments='#')
                self.ui.gixPW.canvas.ax.errorbar(data1[:,0]*self.gixscale[i][0]+self.gixscale[i][1],data1[:,1]*self.gixscale[i][2]+self.gixscale[i][3],data1[:,2]*self.gixscale[i][2],fmt='o',label='#'+str(self.selectedgixfiles_rows[i]+1))
        if  len(self.selectedgixfitfiles_rows)!=0: #plot gixos fit files
            for i in range(len(self.selectedgixfitfiles_rows)):
                data1=np.loadtxt(str(self.gixfitfiles[self.selectedgixfitfiles_rows[i]]), comments='#')
                self.ui.gixPW.canvas.ax.errorbar(data1[:,0]*self.gixfitscale[i][0]+self.gixfitscale[i][1],data1[:,1]*self.gixfitscale[i][2]+self.gixfitscale[i][3],fmt='-',label='#'+str(self.selectedgixfitfiles_rows[i]+1))
        if  self.ui.calgixCB.checkState()!=0:
                self.ui.gixPW.canvas.ax.errorbar(np.array(self.gixcal)[:,0],np.array(self.gixcal)[:,1],fmt='-', label='cal')
        if self.ui.gixlegendCB.checkState()!=0:
            self.ui.gixPW.canvas.ax.legend(loc=self.ui.gixlegendlocCoB.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        if self.ui.gixlogyCB.checkState()!=0:
            self.ui.gixPW.canvas.ax.set_yscale('log')
        else:
            self.ui.gixPW.canvas.ax.set_yscale('linear')
        self.ui.gixPW.canvas.draw()
    
    def setGixPlotScale(self): #set the scale of each data in the gixos plot
        if len(self.selectedgixfiles_rows)+len(self.selectedgixfitfiles_rows)==0:
            self.messageBox('Warning:: No Ref or Fit files selected!')
        else:
            row_gix=len(self.selectedgixfiles_rows)
            row_fit=len(self.selectedgixfitfiles_rows)
            row=row_gix+row_fit
            Dialog=QDialog(self)
            self.uiplotscale=uic.loadUi('plotscale.ui', Dialog)
            self.uiplotscale.scaleTW.setRowCount(row) #set the table size; 4 column is fixed
            self.uiplotscale.show()
            self.uiplotscale.scaleLabel.setText('GIXOS Plot Scale Setup: X=X*Factor+Offset')
            self.uiplotscale.scaleTW.setHorizontalHeaderLabels(QStringList()<<"X Factor"<<"X Offset"<<"Y Factor"<<"Y Offset") #set the horizontal header
            vlabel=QStringList() #set the vertical header 
            for i in range(row_gix):
                vlabel.append("Gix #"+str(self.selectedgixfiles_rows[i]+1))
            for i in range(row_fit):
                vlabel.append("Fit #"+str(self.selectedgixfitfiles_rows[i]+1))
            self.uiplotscale.scaleTW.setVerticalHeaderLabels(vlabel)
            for i in range(row_gix):  #set the initial values
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i,j,QTableWidgetItem(str(self.gixscale[i][j])))
                    self.uiplotscale.scaleTW.item(i,j).setTextAlignment(Qt.AlignCenter)
            for i in range(row_fit):
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i+row_gix,j,QTableWidgetItem(str(self.gixfitscale[i][j])))
                    self.uiplotscale.scaleTW.item(i+row_gix,j).setTextAlignment(Qt.AlignCenter)
            self.connect(self.uiplotscale.scaleTW, SIGNAL('cellChanged(int,int)'), self.updateGixPlotScale) #update the ref scale and plot
            self.connect(self.uiplotscale.closePB,SIGNAL('clicked()'), self.closePlotScale) #close the scale setup window    
    
                           
    def updateGixPlotScale(self): #update the scale of each data in the gixos plot
        row_gix=len(self.selectedgixfiles_rows)
        row_fit=len(self.selectedgixfitfiles_rows)
        self.gixscale=[[float(str(self.uiplotscale.scaleTW.item(i,j).text())) for j in range(4)] for i in range(row_gix)]
        self.gixfitscale=[[float(str(self.uiplotscale.scaleTW.item(i+row_gix,j).text())) for j in range(4)] for i in range(row_fit)]
        self.updateGixPlot()
        
    def updateGixFitFile(self): #update gixos fit files in the listwidget
        self.ui.gixfitfileLW.clear()
        for i in range(len(self.gixfitfiles)):
            try:
                self.ui.gixfitfileLW.addItem('#'+str(i+1)+self.halftab+str(self.gixfitfiles[i].split('\\')[-2])+'\\'+str(self.gixfitfiles[i].split('\\')[-1]))
            except:
                self.ui.gixfitfileLW.addItem('#'+str(i+1)+self.halftab+str(self.gixfitfiles[i].split('/')[-2])+'/'+str(self.gixfitfiles[i].split('/')[-1]))

    def addGixFitFile(self): #add gixos fit files into the listwidget and deselect gixos fit files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple GIXOS Fit Files to import', directory=self.directory, filter='FIT Files (*fit*)')
        self.gixfitfiles=self.gixfitfiles+map(str, f)
        self.directory=str(QFileInfo(self.gixfitfiles[0]).absolutePath())
        self.updateGixFitFile()
        
    def updateSelectedGixFitFile(self): #update the selected gixos fit files in the listwidget
        selectedgixfitfiles=self.ui.gixfitfileLW.selectedItems()
        self.selectedgixfitfiles_rows=[]
        for item in selectedgixfitfiles:
            self.selectedgixfitfiles_rows.append(self.ui.gixfitfileLW.row(item))
        self.selectedgixfitfiles_rows.sort()
        self.gixfitscale=[[1,0,1,0] for i in range(len(self.selectedgixfitfiles_rows))]
        self.updateGixPlot()
        
    def removeGixFitFile(self):  #remove gix fit files in the listwidget and deselect all gix fit files in the listwidget
        items=self.ui.gixfitfileLW.selectedItems()
        for item in items:
            self.gixfitfiles.pop(self.ui.gixfitfileLW.row(item))
        self.ui.gixfitfileLW.clear()
        self.updateGixFitFile()
        
        
        
    def updateGixEDFile(self): #update ed files in the listwidget
        self.ui.gixedfileLW.clear()
        for i in range(len(self.gixedfiles)):
            try:
                self.ui.gixedfileLW.addItem('#'+str(i+1)+self.halftab+str(self.gixedfiles[i].split('\\')[-2])+'\\'+str(self.gixedfiles[i].split('\\')[-1]))
            except:
                self.ui.gixedfileLW.addItem('#'+str(i+1)+self.halftab+str(self.gixedfiles[i].split('/')[-2])+'/'+str(self.gixedfiles[i].split('/')[-1]))

    
    def addGixEDFile(self): #add gixos ed files into the listwidget and deselect all ed files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple ED Files to import', directory=self.directory, filter='FIT Files (*.sld*; *.ed*;*_sld.txt;*_ed.txt)')
        self.gixedfiles=self.gixedfiles+map(str, f)
        self.directory=str(QFileInfo(self.gixedfiles[0]).absolutePath())
        self.updateGixEDFile()
        
    def updateSelectedGixEDFile(self): #update the selected ed files in the listwidget
        selectedgixedfiles=self.ui.gixedfileLW.selectedItems()
        self.selectedgixedfiles_rows=[]
        for item in selectedgixedfiles:
            self.selectedgixedfiles_rows.append(self.ui.gixedfileLW.row(item))
        self.selectedgixedfiles_rows.sort()
        self.gixedscale=[[1,0,1,0] for i in range(len(self.selectedgixedfiles_rows))]
        self.updateGixEDPlot()
        
    def removeGixEDFile(self):  #remove ed files in the listwidget and deselect all ed files in the listwidget
        items=self.ui.gixedfileLW.selectedItems()
        for item in items:
            self.gixedfiles.pop(self.ui.gixedfileLW.row(item))
        self.ui.gixedfileLW.clear()
        self.updateGixEDFile()
        
    def updateGixEDPlot(self): #update the plot in the gix ed plotwidget
        self.ui.gixsldPW.canvas.ax.clear()
        self.ui.gixsldPW.canvas.ax.set_xlabel('Z'+' '+r'$[\AA]$')
        self.ui.gixsldPW.canvas.ax.set_ylabel('Electron Density Profile'+' '+r'$[e/\AA^{3}]$')
        if  len(self.selectedgixedfiles_rows)!=0: #plot gixos ed files
            for i in range(len(self.selectedgixedfiles_rows)):
                data1=np.loadtxt(str(self.gixedfiles[self.selectedgixedfiles_rows[i]]), comments='#')
                self.ui.gixsldPW.canvas.ax.errorbar(data1[:,0]*self.gixedscale[i][0]+self.gixedscale[i][1],data1[:,1]*self.gixedscale[i][2]+self.gixedscale[i][3],fmt='-',label='#'+str(self.selectedgixedfiles_rows[i]+1))
        if  self.ui.calgixsldCB.checkState()!=0:
                self.ui.gixsldPW.canvas.ax.errorbar(np.array(self.gixsldcal)[:,0],np.array(self.gixsldcal)[:,1],fmt='-', label='cal')
        if self.ui.gixedlegendCB.checkState()!=0:
            self.ui.gixsldPW.canvas.ax.legend(loc=self.ui.gixedlegendlocCoB.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        self.ui.gixsldPW.canvas.draw()
        
    def setGixEDPlotScale(self): #set the scale of each data in the ed plot
        if len(self.selectedgixedfiles_rows)==0:
            self.messageBox('Warning:: No electron density files selected!')
        else:
            row_ed=len(self.selectedgixedfiles_rows)
            Dialog=QDialog(self)
            self.uiplotscale=uic.loadUi('plotscale.ui', Dialog)
            self.uiplotscale.scaleTW.setRowCount(row_ed) #set the table size; 4 column is fixed
            self.uiplotscale.show()
            self.uiplotscale.scaleLabel.setText('Electron Density Plot Scale Setup: X=X*Factor+Offset')
            self.uiplotscale.scaleTW.setHorizontalHeaderLabels(QStringList()<<"X Factor"<<"X Offset"<<"Y Factor"<<"Y Offset") #set the horizontal header
            vlabel=QStringList() #set the vertical header 
            for i in range(row_ed):
                vlabel.append("ED #"+str(self.selectedgixedfiles_rows[i]+1))
            self.uiplotscale.scaleTW.setVerticalHeaderLabels(vlabel)
            for i in range(row_ed):  #set the initial values
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i,j,QTableWidgetItem(str(self.gixedscale[i][j])))
                    self.uiplotscale.scaleTW.item(i,j).setTextAlignment(Qt.AlignCenter)
            self.connect(self.uiplotscale.scaleTW,SIGNAL('cellChanged(int,int)'), self.updateGixEDPlotScale)
            self.connect(self.uiplotscale.closePB,SIGNAL('clicked()'), self.closePlotScale)
           
                       
    def updateGixEDPlotScale(self): #update the scale of each data in the ed plot
        row_ed=len(self.selectedgixedfiles_rows)
        self.gixedscale=[[float(str(self.uiplotscale.scaleTW.item(i,j).text())) for j in range(4)] for i in range(row_ed)]
        self.updateGixEDPlot()
        
    def insGixSlab(self, selrows=None): #add a slab in gix par table
        if selrows==None:
            insrows=self.ui.gixparTW.selectionModel().selectedRows()
            insrows=[self.ui.gixparTW.row(self.ui.gixparTW.itemFromIndex(insrows[i])) for i in range(len(insrows))]
        else:
            insrows=[selrows]
        if len(insrows)!=1:
            self.messageBox('Warning:: Only one row can be selected!')
        elif insrows[0]==0:
            self.messageBox('Warning:: Cannot insert a layer above the top phase!')
        else:
            self.disconnect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self.updateGixParaVal)
            insrow=insrows[0]
            self.ui.gixparTW.insertRow(insrow)
            for i in range(4):
                self.ui.gixparTW.setItem(insrow,i,QTableWidgetItem('10/0.3/0/3'.split('/')[i]))
            self.connect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self.updateGixParaVal)
            self.ui.gixnumslabSB.setValue(self.ui.gixparTW.rowCount()-2)
            for i in list(reversed(range((insrow-1)*4+3,4*(self.ui.gixparTW.rowCount()-3)+5))):   #update the parameter dictionary    
                self.gixpara[i+4]=self.gixpara[i]
            self.addGixParaDic(insrow)
            self.updateGixParaName()  #update the paramenter name list
          #  print self.refparaname
           # print self.refpara
            self.updateGixCal()     
            
    def addGixParaDic(self,row):
        self.gixpara[(row-1)*4+3]=[10,False, None,None]
        self.gixpara[(row-1)*4+4]=[0.3,False, None,None]
        self.gixpara[(row-1)*4+5]=[0,False, None,None]
        self.gixpara[(row-1)*4+6]=[3,False, None,None]
                    
    def updateGixParaVal(self):
        selcol=self.ui.gixparTW.currentColumn()
        if self.ui.gixroughCB.checkState()!=0 and selcol==3:
            self.sameGixRough()  #fix all roughness
        self.updateGixCal()   

    def sameGixRough(self):
        row=self.ui.gixparTW.rowCount()
        samerough=float(str(self.ui.gixparTW.item(0,3).text()))
        self.disconnect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self. updateGixParaVal)
        for i in range(1,row-1):
            self.ui.gixparTW.setItem(i,3,QTableWidgetItem(str(samerough)))
            self.gixpara[i*4+2][0]=samerough
        self.connect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self. updateGixParaVal)
            
    def rmGixSlab(self, selrows=None): #remove multiple slabs in gix par table
        row=self.ui.gixparTW.rowCount()
        rmrows=self.ui.gixparTW.selectionModel().selectedRows()   
        removerows=[]
        if selrows==None:
            for rmrow in rmrows:
                removerows.append(self.ui.gixparTW.row(self.ui.gixparTW.itemFromIndex(rmrow)))
                removerows.sort(reverse=True)   #remove the lower layer first
        else:
            removerows=selrows
        if len(removerows)==0:
            self.messageBox('Warning:: No layer is selected')
        else:
            for i in range(len(removerows)):
                if removerows[i] == 0:
                    self.messageBox('Warning:: Cannot remove the top phase!')
                elif removerows[i] == row-1:
                    self.messageBox('Warning:: Cannot remove the bottom phase!')
                else:       
                    self.ui.gixparTW.removeRow(removerows[i])
                    for i in range(removerows[i]*4+3,len(self.gixpara)):  #update the parameter dictionary
                        self.gixpara[i-4]=self.gixpara[i]                   #shift the parameters below the deleting row up
                    for key in range(len(self.gixpara)-4,len(self.gixpara)): #delete the last four parameters
                        self.gixpara.pop(key)
          #  print self.refpara
            self.updateGixParaName()  #update the paramenter name list
            self.ui.gixnumslabSB.setValue(self.ui.gixparTW.rowCount()-2)
            self.updateGixCal()

            
    def modGixSlab(self): #modify Gix par table based on the change of spin box
        diff=self.ui.gixparTW.rowCount()-self.ui.gixnumslabSB.value()-2
        row=self.ui.gixparTW.rowCount()
        if diff>0:   #remove 
            selrows=[]
            for i in range(diff):
                selrows.append(row-2-i)
            self.rmGixSlab(selrows)
        elif diff<0:   #insert
            for i in range(-diff):
                self.insGixSlab(row-1)

    def updateGixCal(self):  #calculate the GIXOS and SLD based on current parameters
        row=self.ui.gixparTW.rowCount()       
        d=[float(str(self.ui.gixparTW.item(i+1,0).text())) for i in range(row-2)]
        rho=[float(str(self.ui.gixparTW.item(i,1).text())) for i in range(row)]
        mu=[float(str(self.ui.gixparTW.item(i,2).text())) for i in range(row)]
        sigma=[float(str(self.ui.gixparTW.item(i,3).text())) for i in range(row-1)]
        syspara=[float(self.ui.gixqoffLE.text()),float(self.ui.gixyscaleLE.text()),float(self.ui.gixqmaxLE.text())] 
        if self.gixsavefitindex==1:
            xgix=np.linspace(self.gixxmin,self.gixxmax,self.gixnp)
            self.gixcal=np.vstack((xgix,self.gixCalFun(d,rho,mu,sigma,syspara,xgix))).T
        elif self.gixsavefitindex==2:
            xsld=np.linspace(self.gixedxmin,self.gixedxmax,self.gixnp)
            self.gixsldcal=np.vstack((xsld,self.sldCalFun(d,rho,sigma,xsld))).T  
        else:
            if  self.ui.calgixCB.checkState()!=0:
                if  len(self.selectedgixfiles_rows)!=0: 
                    data=np.loadtxt(str(self.gixfiles[self.selectedgixfiles_rows[0]]), comments='#')
                    self.gixxmax=np.max(data[:,0])
                    self.gixxmin=np.min(data[:,0])
                else:
                    self.gixxmax=0.7
                    self.gixxmin=0
                xgix=np.linspace(self.gixxmin,self.gixxmax,800)
                self.gixcal=np.vstack((xgix,self.gixCalFun(d,rho,mu, sigma,syspara,xgix))).T
            self.updateGixPlot()
            
            if  self.ui.calgixsldCB.checkState()!=0:
                if sigma[0]!=0 and sigma[-1]!=0:
                    xsld=np.linspace(-4*sigma[0],np.sum(d)+4*sigma[-1],800)
                else:
                    xsld=np.linspace(-10,np.sum(d)+10,800)
                self.gixsldcal=np.vstack((xsld,self.sldCalFun(d,rho,sigma,xsld))).T
            self.updateGixEDPlot()

    def gixCalFun(self,d,rho,mu,sigma,syspara,x):
        qoff=syspara[0] 
        x=qoff+x
        yscale=syspara[1]
        qmax=syspara[2]  
        temperature=float(self.ui.gixtemLE.text())+273.15 # temperature in Kelvin
        tension=float(self.ui.gixtenLE.text())/1000  # tension in N/m
        d=[abs(d[i]) for i in range(len(d))]
        rho=[abs(rho[i]) for i in range(len(rho))]
        mu=[abs(mu[i]) for i in range(len(mu))]
        sigma=[abs(sigma[i]) for i in range(len(sigma))]
        erad=self.eleradius   # classic electron radius
        k0=2*np.pi*float(self.ui.gixxenLE.text())/12.3984 # wave vector
        slab=0.25
        del_alpha=np.arcsin(qoff/2/k0)  # get alpha offset from qoffset
        alpha=float(self.ui.gixalphaLE.text())/180*np.pi+del_alpha  # get alpha
        beta=np.arcsin(x/k0-np.sin(alpha))  #get beta
        twoth=float(self.ui.gixdthLE.text())/180*np.pi #get two theta
        qxy=k0*np.sqrt(np.cos(alpha)**2+np.cos(beta)**2-2*np.cos(alpha)*np.cos(beta)*np.cos(twoth))  #get qxy
       # ftprint=float(self.ui.gixincsliLE.text())/np.sin(alpha)  #get the footprint 
      #  delbeta=np.sqrt((ftprint*np.sin(beta))**2+float(self.ui.gixoutsliLE.text()))  #get the delta beta
        length=np.sum(d)+4*(sigma[0]+sigma[-1]) # total length of inner slabs plus 4 times rougness for both sides
        steps=int(length/slab) # each sliced box has thickness of ~ 0.25 \AA
        xsld=np.linspace(-4*sigma[0],np.sum(d)+4*sigma[-1],steps) # get the x-axis for sld
        intrho=self.sldCalFun(d,rho,sigma,xsld)
        intmu=self.sldCalFun(d,mu,sigma,xsld)
        sd=length/steps # thickness for each slab
        sdel=[]
        sbet=[] 
        sdel.append(erad*2.0*np.pi/k0/k0*rho[0]) # delta for the top phase
        sbet.append(mu[0]/2/k0/1e8)        # beta for the top phase
        sdel=sdel+[intrho[i]*erad*2.0*np.pi/k0/k0 for i in range(len(intrho))] # add delta for the interface
        sbet=sbet+[intmu[i]/2/k0/1e8 for i in range(len(intmu))] # add beta for the interface
        sdel.append(erad*2.0*np.pi/k0/k0*rho[-1])  # delta for the bottom phase
        sbet.append(mu[-1]/2/k0/1e8)    # beta for the bottom phase         
        d=slab*np.ones_like(sdel)
        lamda=2*np.pi/k0
        fdel=erad*2.0*np.pi/k0/k0
        sdelf=np.array(sdel)/fdel
        ref,refr=xr.parratt(x,lamda,d,sdelf,sbet)
        frsnll,frsnl1=xr.parratt(x,lamda,[0,1],[sdelf[0],sdelf[-1]],[sbet[0],sbet[-1]])    
        eta=self.boltzmann*temperature*x**2*1e20/2/np.pi/tension
        ##get Fresnel transmission 
        trans_bet=[]
        for i in range(len(beta)):
            if beta[i]>0:
                f1=cmath.sqrt(complex(beta[i]*beta[i],2*sbet[0]))
                fmax=cmath.sqrt(complex(beta[i]*beta[i]-2*(sdel[-1]-sdel[0]),2*sbet[-1]))
                trans=4*abs(f1/(f1+fmax))*abs(f1/(f1+fmax))
            else:
                trans=0
            trans_bet.append(trans)
        return yscale*ref/frsnll*np.array(trans_bet)*eta/x**2/qmax**eta*qxy**(eta-2)*np.cos(beta)#*delbeta
        
        
    def setupGixPara(self): # constrains setup for gix parameters
        Dialog=QDialog(self)
        self.uigixpara=uic.loadUi('refpara.ui', Dialog)
        selrow=self.ui.gixparTW.currentRow()
        selcol=self.ui.gixparTW.currentColumn()
        if selrow==self.ui.gixparTW.rowCount()-1:
            self.paranum=selrow*4+selcol-2
        else:
            self.paranum=selrow*4+selcol-1
        self.uigixpara.label.setText('Limits Setup of Parameter:'+self.gixparaname[self.paranum])
        if self.gixpara[self.paranum][2]!=None: 
            self.uigixpara.minCB.setCheckState(2)
            self.uigixpara.minLE.setText(str(self.gixpara[self.paranum][2]))
        if self.gixpara[self.paranum][3]!=None: 
            self.uigixpara.maxCB.setCheckState(2)
            self.uigixpara.maxLE.setText(str(self.gixpara[self.paranum][3]))
        self.uigixpara.show()
        self.connect(self.uigixpara.cancelPB, SIGNAL('clicked()'), self.cancelGixPara)
        self.connect(self.uigixpara.okPB, SIGNAL('clicked()'), self.takeGixPara)
        
    def cancelGixPara(self):
        self.uigixpara.close()
        
    def takeGixPara(self):
        if self.uigixpara.minCB.checkState()!=0 and self.uigixpara.maxCB.checkState()!=0 and float(self.uigixpara.minLE.text())>float(self.uigixpara.maxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain!!!")
        else:
            if self.uigixpara.minCB.checkState()!=0:
                self.gixpara[self.paranum][2]=float(self.uigixpara.minLE.text())
            else:
                self.gixpara[self.paranum][2]=None
            if self.uigixpara.maxCB.checkState()!=0:
                self.gixpara[self.paranum][3]=float(self.uigixpara.maxLE.text())
            else:
                self.gixpara[self.paranum][3]=None
            self.uigixpara.close()
    
    def cleGixCon(self):
        for i in range(len(self.gixpara)):
            self.gixpara[i][2]=None
            self.gixpara[i][3]=None
        for i in range(len(self.gixsyspara)):
            self.gixsyspara[i][2]=None
            self.gixsyspara[i][3]=None
            
            
            
    def updateGixSysPara(self):
        Dialog=QDialog(self)
        self.uigixsyspara=uic.loadUi('gixsyspara.ui', Dialog)
        if self.gixsyspara[0][2]!=None:  #set up the current value
            self.uigixsyspara.qoffminCB.setCheckState(2)
            self.uigixsyspara.qoffminLE.setText(str(self.gixsyspara[0][2]))
        if self.gixsyspara[0][3]!=None:  
            self.uigixsyspara.qoffmaxCB.setCheckState(2)
            self.uigixsyspara.qoffmaxLE.setText(str(self.gixsyspara[0][3]))
        if self.gixsyspara[1][2]!=None:  
            self.uigixsyspara.yscaleminCB.setCheckState(2)
            self.uigixsyspara.yscaleminLE.setText(str(self.gixsyspara[1][2]))
        if self.gixsyspara[1][3]!=None:  
            self.uigixsyspara.yscalemaxCB.setCheckState(2)
            self.uigixsyspara.yscalemaxLE.setText(str(self.gixsyspara[1][3]))
        if self.gixsyspara[2][2]!=None:  
            self.uigixsyspara.qmaxminCB.setCheckState(2)
            self.uigixsyspara.qmaxminLE.setText(str(self.gixsyspara[2][2]))
        if self.gixsyspara[2][3]!=None:  
            self.uigixsyspara.qmaxmaxCB.setCheckState(2)
            self.uigixsyspara.qmaxmaxLE.setText(str(self.gixsyspara[2][3]))
        self.uigixsyspara.show()
        self.connect(self.uigixsyspara.cancelPB, SIGNAL('clicked()'), self.cancelGixSysPara)
        self.connect(self.uigixsyspara.okPB, SIGNAL('clicked()'), self.takeGixSysPara)
        
    def cancelGixSysPara(self):
        self.uigixsyspara.close()
                
    def takeGixSysPara(self):
        if self.uigixsyspara.qoffminCB.checkState()!=0 and self.uigixsyspara.qoffmaxCB.checkState()!=0 and float(self.uigixsyspara.qoffminLE.text())>float(self.uigixsyspara.qoffmaxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain for Q offset!!!")
        elif self.uigixsyspara.yscaleminCB.checkState()!=0 and self.uigixsyspara.yscalemaxCB.checkState()!=0 and float(self.uigixsyspara.yscaleminLE.text())>float(self.uigixsyspara.yscalemaxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain for Y scale!!!")
        elif self.uigixsyspara.qmaxminCB.checkState()!=0 and self.uigixsyspara.qmaxmaxCB.checkState()!=0 and float(self.uigixsyspara.qmaxminLE.text())>float(self.uigixsyspara.qmaxmaxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain for Q_max!!!")
        else:
            if self.uigixsyspara.qoffminCB.checkState()!=0:
                self.gixsyspara[0][2]=float(self.uigixsyspara.qoffminLE.text())
            else:
                self.gixsyspara[0][2]=None
            if self.uigixsyspara.qoffmaxCB.checkState()!=0:
                self.gixsyspara[0][3]=float(self.uigixsyspara.qoffmaxLE.text())
            else:
                self.gixsyspara[0][3]=None
            if self.uigixsyspara.yscaleminCB.checkState()!=0:
                self.gixsyspara[1][2]=float(self.uigixsyspara.yscaleminLE.text())
            else:
                self.gixsyspara[1][2]=None
            if self.uigixsyspara.yscalemaxCB.checkState()!=0:
                self.gixsyspara[1][3]=float(self.uigixsyspara.yscalemaxLE.text())
            else:
                self.gixsyspara[1][3]=None
            if self.uigixsyspara.qmaxminCB.checkState()!=0:
                self.gixsyspara[2][2]=float(self.uigixsyspara.qmaxminLE.text())
            else:
                self.gixsyspara[2][2]=None
            if self.uigixsyspara.qmaxmaxCB.checkState()!=0:
                self.gixsyspara[2][3]=float(self.uigixsyspara.qmaxmaxLE.text())
            else:
                self.gixsyspara[2][3]=None
            self.uigixsyspara.close()
         
    def getGixParaVal(self):
        for i in range(len(self.gixparaname)-2): #get the current values except the bottom phase in the table 
                cell=divmod(i+1,4)  #get the cell index for each parameter
                self.gixpara[i][0]=float(str(self.ui.gixparTW.item(cell[0],cell[1]).text()))
        self.gixpara[len(self.gixparaname)-2][0]=float(str(self.ui.gixparTW.item(cell[0]+1,1).text()))   #last row
        self.gixpara[len(self.gixparaname)-1][0]=float(str(self.ui.gixparTW.item(cell[0]+1,2).text()))
        self.gixsyspara[0][0]=float(self.ui.gixqoffLE.text())  #system parameters
        self.gixsyspara[1][0]=float(self.ui.gixyscaleLE.text())
        self.gixsyspara[2][0]=float(self.ui.gixqmaxLE.text())
        
        
    def fitGix(self):
        self.getGixParaVal()
        index=self.ui.gixparTW.selectionModel().selectedIndexes()
        row=self.ui.gixparTW.rowCount()
        selrows=[self.ui.gixparTW.row(self.ui.gixparTW.itemFromIndex(index[i])) for i in range(len(index))]
        selcols=[self.ui.gixparTW.column(self.ui.gixparTW.itemFromIndex(index[i])) for i in range(len(index))]
        selparas=[]
        for i in range(len(selrows)):  #get selected parameters
            if selrows[i]!=row-1:
                selparas.append(selrows[i]*4+selcols[i]-1)
            else:
                selparas.append(selrows[i]*4+selcols[i]-2)
       # print selparas
        for i in range(len(self.gixpara)):  #set selected parameters to be varied
            if i in selparas:
                self.gixpara[i][1]=True
            else:
                self.gixpara[i][1]=False
       # print self.gixpara
        for i in range(len(self.gixsysCB)):
            if self.gixsysCB[i].checkState()!=0:
                self.gixsyspara[i][1]=True
            else:
                self.gixsyspara[i][1]=False
      #  print self.gixsyspara
        self.gixparameter=Parameters()
        for i in range(len(self.gixpara)):
            self.gixparameter.add(self.gixparaname[i], value=self.gixpara[i][0],vary=self.gixpara[i][1],min=self.gixpara[i][2],max=self.gixpara[i][3])
        for i in range(len(self.gixsysparaname)):
            self.gixparameter.add(self.gixsysparaname[i], value=self.gixsyspara[i][0],vary=self.gixsyspara[i][1],min=self.gixsyspara[i][2],max=self.gixsyspara[i][3])
        #print self.refparameter
       
        if  len(self.selectedgixfiles_rows)!=1: 
            self.messageBox("Please select only one set of data for fitting!")
        else:
            data=np.loadtxt(str(self.gixfiles[self.selectedgixfiles_rows[0]]), comments='#')
            ini=max(float(str(self.ui.gixfitranLE.text()).split(':')[0]),data[0][0])
            fin=min(float(str(self.ui.gixfitranLE.text()).split(':')[1]),data[-1][0])
            data1=data[np.where(np.logical_and(data[:,0]>=ini,data[:,0]<=fin))]
            x=data1[:,0]
            y=data1[:,1]
            if self.ui.gixerrCB.currentIndex()==0:
                yerr=data1[:,2]
            elif self.ui.gixerrCB.currentIndex()==1:
                yerr=np.sqrt(y)
            elif self.ui.gixerrCB.currentIndex()==2:
                yerr=y
            else:
                yerr=np.ones_like(x)
           # print yerr
           # print self.refparameter
            self.gixresult=minimize(self.gix2min, self.gixparameter, args=(x,y,yerr))
            tempchi=self.gixresult.redchi
            self.gixparameter=self.gixresult.params
            self.gixresult=minimize(self.gix2min, self.gixparameter, args=(x,y,yerr))
            while(np.abs(self.gixresult.redchi-tempchi)/tempchi>0.02):
                tempchi=self.gixresult.redchi
                self.gixparameter=self.gixresult.params
                self.gixresult=minimize(self.gix2min, self.gixparameter, args=(x,y,yerr))
            print(fit_report(self.gixresult))
            residual=np.vstack((x,self.gixresult.residual)).T
           # print residual
            self.disconnect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self.updateGixParaVal)
            if self.ui.gixroughCB.checkState()!=0: #enforce the roughness to be same if set
                for i in range(1,row-1):
                    self.gixresult.params[self.gixparaname[4*i+2]].value=self.gixresult.params[self.gixparaname[2]].value
            for i in range(len(self.gixparaname)-2): #put the best values except the bottom phase in the table 
                cell=divmod(i+1,4)  #get the cell index for each parameter
               # print str(result.params[self.refparaname[i]].value)
                self.ui.gixparTW.setItem(cell[0],cell[1],QTableWidgetItem(format(self.gixresult.params[self.gixparaname[i]].value,'.4f')))
            self.ui.gixparTW.setItem(row-1,1,QTableWidgetItem(format(self.gixresult.params[self.gixparaname[-2]].value, '.4f'))) # put the best values for the bottom phase
            self.ui.gixparTW.setItem(row-1,2,QTableWidgetItem(format(self.gixresult.params[self.gixparaname[-1]].value, '.4f')))
            self.ui.gixqoffLE.setText(format(self.gixresult.params[self.gixsysparaname[0]].value, '.6f'))  #put the best sys parameter values 
            self.ui.gixyscaleLE.setText(format(self.gixresult.params[self.gixsysparaname[1]].value, '.3e'))
            self.ui.gixqmaxLE.setText(format(self.gixresult.params[self.gixsysparaname[2]].value, '.3f'))
            self.connect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self.updateGixParaVal)
            self.ui.calgixCB.setCheckState(2)
            self.updateGixCal()
            self.ui.gixChiLE.setText(format(self.gixresult.redchi, '.3f'))
            self.ui.gixparaTB.clear()
            fitinfo='Fitting Paramenters:\n'
            fitinfo=fitinfo+'Name\tStderr\tMin\tMax\n'
            selparas.sort()
            for i in selparas:
                fitinfo=fitinfo+self.gixparaname[i]+'\t'+format(self.gixresult.params[self.gixparaname[i]].stderr, '.4f')+'\t'+str(self.gixpara[i][2])+'\t'+str(self.gixpara[i][3])+'\n'
            for i in range(len(self.gixsysCB)):
                if self.gixsyspara[i][1]==True:
                   fitinfo=fitinfo+self.gixsysparaname[i]+'\t'+format(self.gixresult.params[self.gixsysparaname[i]].stderr, '.4f')+'\t'+str(self.gixsyspara[i][2])+'\t'+str(self.gixsyspara[i][3])+'\n' 
            fitinfo=fitinfo+'********************************\n'
            fitinfo=fitinfo+'Fitting Residual:\n'
            for i in range(len(residual)):
                fitinfo=fitinfo+format(residual[i][0], '.3f')+'\t'+format(residual[i][1], '.4f')+'\n'
            self.ui.gixparaTB.append(fitinfo)
            cursor=self.ui.gixparaTB.textCursor()
            cursor.setPosition(0)
            self.ui.gixparaTB.setTextCursor(cursor)
            
    def gix2min(self, params, x, y, yerr): #residuel for gixos fitting
        row=self.ui.gixparTW.rowCount()
        d=[params[self.gixparaname[i*4+3]].value for i in range(row-2)]
        rho=[params[self.gixparaname[i*4]].value for i in range(row-1)]
        mu=[params[self.gixparaname[i*4+1]].value for i in range(row-1)]
        sigma=[params[self.gixparaname[i*4+2]].value for i in range(row-1)]
        rho.append(params[self.gixparaname[-2]].value)  #add bottom phase
        mu.append(params[self.gixparaname[-1]].value)  #add bottom phase
        syspara=[params[self.gixsysparaname[i]].value for i in range(3)]
        model=self.gixCalFun(d,rho,mu,sigma,syspara,x)
        return (model-y)/yerr
        
    def saveGix(self):
        if str(self.ui.gixsaveCB.currentText())=='Save Fit':
            self.gixsavefitindex=1
            self.saveGixFitDig()
        elif str(self.ui.gixsaveCB.currentText())=='Save ED':
            self.gixsavefitindex=2
            self.saveGixFitDig()
        elif str(self.ui.gixsaveCB.currentText())=='Save Para':
            self.saveGixPara()
        elif str(self.ui.gixsaveCB.currentText())=='Save Data':
            self.saveGixData()
            
    def saveGixFitDig(self):
        Dialog=QDialog(self)
        self.uigixsavefit=uic.loadUi('refsave.ui', Dialog)
        if self.gixsavefitindex==1:
            self.uigixsavefit.label.setText('Save GIXOS Fit/Calcualtion!')
            try:
                self.uigixsavefit.xminLE.setText(str(self.gixxmin))
                self.uigixsavefit.xmaxLE.setText(str(self.gixxmax))
            except:
                pass
        elif self.gixsavefitindex==2: 
            row=self.ui.gixparTW.rowCount()       
            d=[float(str(self.ui.gixparTW.item(i+1,0).text())) for i in range(row-2)]
            sigma=[float(str(self.ui.gixparTW.item(i,3).text())) for i in range(row-1)]
            self.uigixsavefit.label.setText('Save Electron Density Profile!')
            self.uigixsavefit.xminLE.setText(str(-4*sigma[0]))
            self.uigixsavefit.xmaxLE.setText(str(np.sum(d)+4*sigma[-1]))
        self.uigixsavefit.numpointLE.setText(str(400))
        self.uigixsavefit.show()
        self.connect(self.uigixsavefit.cancelPB, SIGNAL('clicked()'), self.cancelSaveGixFit)
        self.connect(self.uigixsavefit.okPB, SIGNAL('clicked()'), self.saveGixFit)
        
    def cancelSaveGixFit(self):
        self.uigixsavefit.close()
        self.gixsavefitindex=0
        
    def saveGixFit(self):
        self.gixnp=float(self.uigixsavefit.numpointLE.text())
        if float(self.uigixsavefit.xminLE.text())>=float(self.uigixsavefit.xmaxLE.text()) or float(self.uigixsavefit.numpointLE.text())<=0:
            self.messageBox('Warning::Starting points must be lower than ending points \n and numer of points must be large than 0!!')
        else:
            if self.gixsavefitindex==1:       
                self.gixxmin=float(self.uigixsavefit.xminLE.text())
                self.gixxmax=float(self.uigixsavefit.xmaxLE.text())
                self.updateGixCal()
                self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Gixos Fit Data',directory=self.directory))
                fname=self.saveFileName+'_fit.txt'
                np.savetxt(fname,self.gixcal,fmt='%.4f\t%.4e')
            elif self.gixsavefitindex==2:       
                self.gixedxmin=float(self.uigixsavefit.xminLE.text())
                self.gixedxmax=float(self.uigixsavefit.xmaxLE.text())
                self.updateGixCal()
                self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Electron Density Data',directory=self.directory))
                fname=self.saveFileName+'_ed.txt'
                np.savetxt(fname,self.gixsldcal,fmt='%.4e\t%.4e')
            self.gixsavefitindex=0
            self.uigixsavefit.close()   
        
        
    
    def saveGixPara(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save GIXOS Fitting Parameters',directory=self.directory))
        fid=open(self.saveFileName+'_par.txt','w')
        try:
            fid.write('Chi_Square\t'+format(self.gixresult.redchi, '.3f')+'\n')  #chisquare
        except:
            fid.write('Chi_Square\tNA\n') 
        fid.write('Error_Type\t'+str(self.ui.gixerrCB.currentText()).split()[0]+'\n')
        fid.write('Num_of_Layer\t'+str(self.ui.gixparTW.rowCount()-2)+'\n')    #number of layers
        if self.ui.gixroughCB.checkState()!=0:
            fid.write('Roughness\tFixed\n')
        else:
            fid.write('Roughness\tVary\n')
        fid.write('Para_Name\tValue\t\tVary\tStderr\t\tMin\tMax\n')
        for i in range((self.ui.gixparTW.rowCount()-2)*4+5):
            try:
                fid.write(self.gixparaname[i]+'\t\t'+format(self.gixresult.params[self.gixparaname[i]].value,'.3e')+'\t'+str(self.gixpara[i][1])+'\t'+format(self.gixresult.params[self.gixparaname[i]].stderr,'.3e')+'\t'+str(self.gixpara[i][2])+'\t'+str(self.gixpara[i][3])+'\n')
            except:
                if i <=(self.ui.gixparTW.rowCount()-2)*4+2:
                    cell=divmod(i+1,4)
                else:
                    cell=divmod(i+2,4)
                fid.write(self.gixparaname[i]+'\t\t'+str(self.ui.gixparTW.item(cell[0],cell[1]).text())+'\tNA\tNA\tNA\tNA\n')
        for i in range(len(self.gixsyspara)):
            try:
                fid.write(self.gixsysparaname[i]+'\t\t'+format(self.gixresult.params[self.gixsysparaname[i]].value,'.3e')+'\t'+str(self.gixsyspara[i][1])+'\t'+format(self.gixresult.params[self.gixsysparaname[i]].stderr,'.3e')+'\t'+str(self.gixsyspara[i][2])+'\t'+str(self.gixsyspara[i][3])+'\n')
            except:
                temp=[float(self.gixsysLE[i].text()) for i in range(len(self.gixsysLE))]   
                fid.write(self.gixsysparaname[i]+'\t\t'+str(temp[i])+'\tNA\tNA\tNA\tNA\n')
        fid.write('Constants:\n')
        for i in range(len(self.gixconsname)):
            fid.write(str(self.gixconsname[i])+'\t'+format(float(self.gixconLE[i].text()),'.3e')+'\n')                
        fid.close()
   
    def loadGix(self):
        if str(self.ui.gixloadCB.currentText())=='Load Para':
            self.loadGixPara()   
   
    def loadGixPara(self):
        filename=QFileDialog.getOpenFileName(caption='Select Parameter File to read', directory=self.directory, filter='Par Files (*.par*;*_par.txt)')
        self.directory=str(QFileInfo(filename).absolutePath())
        fid=open(filename)
        fdata=fid.readlines()
        fid.close()
        self.ui.calgixCB.setCheckState(0)
        self.ui.calgixsldCB.setCheckState(0)
        nlayer=eval(fdata[2].split('\t')[1])  #get number of layers
        roughness=fdata[3][:-1].split('\t')[1]  # get the roughness
        if roughness=='Fixed':
            self.ui.gixroughCB.setCheckState(2)
        else:
            self.ui.gixroughCB.setCheckState(0)
        self.ui.gixnumslabSB.setValue(nlayer)
        self.disconnect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self. updateGixParaVal)
        para=[]
        for i in range(5,5+nlayer*4+5+len(self.gixsysparaname)):
            para.append(eval(fdata[i].split('\t')[2]))
        for i in range(nlayer*4+5):   #for structure parameter
            if i<=nlayer*4+2:
                cell=divmod(i+1,4)
            else:
                cell=divmod(i+2,4)
            self.ui.gixparTW.setItem(cell[0],cell[1],QTableWidgetItem(str(para[i])))
            self.gixpara[i][0]=para[i]            
        for i in range(len(self.gixsysparaname)):   #for system parameter
            self.gixsysLE[i].setText(format(para[i+nlayer*4+5],'.3e'))
        cons=[]
        for i in range(5+nlayer*4+5+len(self.gixsysparaname)+1, len(fdata)):
            cons.append(eval(fdata[i].split('\t')[1]))
        for i in range(len(self.gixconLE)):
            self.gixconLE[i].setText(format(cons[i],'.3e'))
        self.connect(self.ui.gixparTW,SIGNAL('cellChanged(int,int)'), self. updateGixParaVal)     
        self.ui.calgixCB.setCheckState(2)
        self.ui.calgixsldCB.setCheckState(2)
        self.updateGixCal()    
    
    
    
    def showAbout(self):
        cwd=os.getcwd()
        files=['mainwindow.py','mainwindow.ui','main.py','mplwidget.py']
        fname=[cwd+'/'+fname for fname in files]
        updateTime=max([os.path.getmtime(fn) for fn in fname])
        self.messageBox('Surface X-ray Scattering Data Analyzer\nVersion: 17.01\nLast Update: '+time.strftime("%m/%d/%Y %I:%M:%S %p",time.localtime(updateTime))+'\nCopyright belongs to:\n\tWei Bu <bu@cars.uchicago.edu>',title='About')
   
    def messageBox(self,text,title='Warning'):
        mesgbox=QMessageBox()
        mesgbox.setText(text)
        mesgbox.setWindowTitle(title)
        mesgbox.exec_()    
        
    def test(self):
        print 'I am here'
    