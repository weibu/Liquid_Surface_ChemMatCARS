from PyQt4 import uic, QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from spec_routines import specread
from mca_routines import mcaread
from mplwidget import MplWidget
from matplotlib.widgets import MultiCursor
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from TwoDDetector import TwoDDetector
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
from sympy.solvers import solve 
from sympy import Symbol
import scipy as sy
import pylab as pl
import numpy as np
import time
import os
import glob
import sys
import matplotlib as mpl
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
from functools import partial
#from TwoD_Integrate import integrate_2d


(Ui_MainWindow, QMainWindow) = uic.loadUiType('mainwindow.ui')

class MainWindow (QMainWindow):
    """MainWindow inherits QMainWindow"""

    def __init__ (self, parent = None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.pilatus=TwoDDetector('Pilatus')
        self.bruker=TwoDDetector('Bruker')
        np.seterr(invalid='ignore',divide='ignore')
        self.mcaDirName='vortex'
        self.ccdDirName='apex'
        self.pilDirName='pilatus'
        self.beamline='APS-15IDC'
        self.specData={}
        self.specPar={}
        self.startLineNum=0
        self.endScanNum=0
        self.refData=[]
        self.refInfo=[]
        self.mcaIntData=[]
        self.refQc=float(self.ui.refQcLineEdit.text())
        self.absfac=float(self.ui.AbsFacLineEdit.text())
        self.reffiles=[]
        self.selectedBgFrameNums=[]
        self.selectedreffiles_rows=[]
        self.mcafiles=[]
        self.selectedmcafiles_rows=[]
        self.selectedplotfiles_rows=[]
        self.ui.PlotWidget.setCurrentIndex(0)
        self.directory=os.getcwd()
        self.pilGIDshow=0
        self.pilGISAXSshow=0
        self.pilmovieindex=0
        self.refpatchindex=0
        self.formfactor_r=0
        self.formfactor_rxy=0
        self.formfactor_rz=0
        self.mcafitstatus=0
        self.mcafitallstatus=0
        self.cutDirItems=['H Cut', 'V Cut', 'Qz Cut', 'Qxy Cut']
        self.ui.pilCutDirComboBox.clear()
        self.ui.pilCutDirComboBox.addItems(self.cutDirItems)
        self.plotfiles=[]
        self.twodplotfiles=[]
        self.halftab='      '
        self.xyzformat='x=%.3f,y=%.3f,z=%.2e'
        # mpl.rc('axes',color_cycle=['b','r','g','c','m','y','k'])
        self.ui.plotMoveUpPushButton.setIcon(QIcon('arrow_up.png'))
        self.ui.plotMoveDownPushButton.setIcon(QIcon('arrow_down.png'))        
        self.connect(self.ui.actionAbout, SIGNAL('triggered()'),self.showAbout)
        self.connect(self.ui.actionOpen_Spec_File, SIGNAL('triggered()'),self.openSpecFile)
        self.connect(self.ui.speUpSpeFilePushButton, SIGNAL('clicked()'),self.readSpecFile)
        self.connect(self.ui.actionAPS_15IDC, SIGNAL('triggered()'),self.selectAPS_15IDC)
        self.connect(self.ui.actionAPS_9IDC, SIGNAL('triggered()'),self.selectAPS_9IDC)
        self.connect(self.ui.actionG_l2, SIGNAL('triggered()'), self.calg_l2)
        self.connect(self.ui.actionG_l3, SIGNAL('triggered()'), self.calg_l3)
        self.connect(self.ui.actionAbs_ratio, SIGNAL('triggered()'), self.calabsrat)
        self.connect(self.ui.scanListWidget, SIGNAL('itemSelectionChanged()'),self.scanListChanged)
       # self.connect(self.ui.scanListWidget, SIGNAL('itemClicked(QListWidgetItem*)'),self.scanListChanged)
        self.connect(self.ui.spLogXCheckBox, SIGNAL('stateChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spLogYCheckBox, SIGNAL('stateChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spGridCheckBox, SIGNAL('stateChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spLegendCheckBox, SIGNAL('stateChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spNCheckBox, SIGNAL('stateChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spLegendLocComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spXComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spYComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spY2ComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spNComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spExportPushButton, SIGNAL('clicked()'),self.specSaveData)
        self.connect(self.ui.scansLineEdit, SIGNAL('returnPressed()'),self.scanListInputChanged)
        self.connect(self.ui.imageListWidget, SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)        
        self.connect(self.ui.mcaLogXCheckBox, SIGNAL('stateChanged(int)'),self.updateMcaPlotData)
        self.connect(self.ui.mcaLogYCheckBox, SIGNAL('stateChanged(int)'),self.updateMcaPlotData)
        self.connect(self.ui.mcaGridCheckBox, SIGNAL('stateChanged(int)'),self.updateMcaPlotData)
        self.connect(self.ui.mcaLegendCheckBox, SIGNAL('stateChanged(int)'),self.updateMcaPlotData)
        self.connect(self.ui.mcaNormComboBox, SIGNAL('currentIndexChanged(int)'),self.updateMcaPlotData)
        self.connect(self.ui.mcaLegendLocComboBox, SIGNAL('currentIndexChanged(int)'),self.updateMcaPlotData)
        self.connect(self.ui.mcaCalibCheckBox, SIGNAL('stateChanged(int)'),self.updateMcaPlotData)
        self.connect(self.ui.mcaCalibConLineEdit, SIGNAL('returnPressed()'),self.updateMcaPlotData)
        self.connect(self.ui.mcaCalibLinLineEdit, SIGNAL('returnPressed()'),self.updateMcaPlotData)
        self.connect(self.ui.mcaCalibQuaLineEdit, SIGNAL('returnPressed()'),self.updateMcaPlotData)
        self.connect(self.ui.gixSumCheckBox,SIGNAL('stateChanged(int)'),self.update2dPlots)
        self.connect(self.ui.gixBPCfacLineEdit,SIGNAL('returnPressed()'),self.imageSelectedScanChanged)
        self.connect(self.ui.pilBPCfacLineEdit,SIGNAL('returnPressed()'),self.imageSelectedScanChanged)
        self.connect(self.ui.pilPinholeCheckBox,SIGNAL('stateChanged(int)'),self.update2dPlots)
        self.connect(self.ui.gixLogIntCheckBox, SIGNAL('stateChanged(int)'),self.update2dPlots)
        self.connect(self.ui.pilLogIntCheckBox, SIGNAL('stateChanged(int)'),self.update2dPlots)
        self.connect(self.ui.imageSelectAllCheckBox, SIGNAL('clicked()'),self.selectAllImages)
        self.connect(self.ui.imageListWidget, SIGNAL('clicked(QModelIndex)'),self.unSelectedAllImages)
        self.connect(self.ui.gixSpecCheckBox,SIGNAL('stateChanged(int)'),self.update2dPlots)
        self.connect(self.ui.gixCcd_OffLineEdit,SIGNAL('returnPressed()'),self.update2dPlots)
        self.connect(self.ui.gixCMapComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
        self.connect(self.ui.gixMinLineEdit, SIGNAL('returnPressed()'),self.update2dPlots)
        self.connect(self.ui.gixMaxLineEdit, SIGNAL('returnPressed()'),self.update2dPlots)
        self.connect(self.ui.gixAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
        self.connect(self.ui.gixMinHorizontalSlider, SIGNAL('sliderReleased()'),self.updateMinSlider)
        self.connect(self.ui.gixMaxHorizontalSlider, SIGNAL('sliderReleased()'),self.updateMaxSlider)
        self.connect(self.ui.pilSpecCheckBox,SIGNAL('stateChanged(int)'),self.update2dPlots)
        self.connect(self.ui.pilCMapComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
        self.connect(self.ui.pilMinLineEdit, SIGNAL('returnPressed()'),self.update2dPlots)
        self.connect(self.ui.pilMaxLineEdit, SIGNAL('returnPressed()'),self.update2dPlots)
        self.connect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
        self.connect(self.ui.pilMinHorizontalSlider, SIGNAL('sliderReleased()'),self.updatePilMinSlider)
        self.connect(self.ui.pilMaxHorizontalSlider, SIGNAL('sliderReleased()'),self.updatePilMaxSlider)
        self.connect(self.ui.gixRefPushButton,SIGNAL('clicked()'),self.refPlotWin)
        self.connect(self.ui.pilRefPushButton,SIGNAL('clicked()'),self.refPlotWin)
        self.connect(self.ui.pilRefPatchPushButton,SIGNAL('clicked()'),self.pilRefPatchData)
        self.connect(self.ui.gixDBPosLineEdit,SIGNAL('returnPressed()'),self.update2dPlots)
        self.connect(self.ui.gixSDDistLineEdit,SIGNAL('returnPressed()'),self.update2dPlots)
        self.connect(self.ui.refSlitLineEdit,SIGNAL('returnPressed()'),self.refQzListSelectionChanged)
        self.connect(self.ui.refBGOffLineEdit,SIGNAL('returnPressed()'),self.refQzListSelectionChanged)
        self.connect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
        self.connect(self.ui.refBGDirComboBox,SIGNAL('currentIndexChanged(int)'),self.refQzListSelectionChanged)
        self.connect(self.ui.refAnalyzePushButton,SIGNAL('clicked()'),self.refAnalyze)
        self.connect(self.ui.refBPCFacLineEdit, SIGNAL('returnPressed()'),self.refAnalyze)
        #self.connect(self.ui.AbsFacLineEdit, SIGNAL('returnPressed()'),self.refDoBPC)
        self.connect(self.ui.refAcceptPushButton, SIGNAL('clicked()'), self.updateRefDataList)
        self.connect(self.ui.refAnalyzeAllPushButton, SIGNAL('clicked()'),self.refAnalyzeAll)
        #self.connect(self.ui.refDataListWidget,SIGNAL('itemChanged(QListWidgetItem*)'), self.updateRefPlotData)
        self.connect(self.ui.refQcLineEdit,SIGNAL('returnPressed()'),self.updateRefPlotData)
        self.connect(self.ui.refQoffLineEdit,SIGNAL('returnPressed()'),self.updateRefPlotData)
        self.connect(self.ui.refRRFCheckBox,SIGNAL('stateChanged(int)'),self.updateRefPlotData)
        self.connect(self.ui.refLogYCheckBox,SIGNAL('stateChanged(int)'),self.updateRefPlotData)
        self.connect(self.ui.refQzSqrCheckBox,SIGNAL('stateChanged(int)'),self.updateRefPlotData)
        self.connect(self.ui.refDelPushButton,SIGNAL('clicked()'),self.delRefDataList)
        self.connect(self.ui.refSeaWinLineEdit, SIGNAL('returnPressed()'),self.refAnalyze)
        self.connect(self.ui.refExportPushButton,SIGNAL('clicked()'), self.saveRefData)
        self.connect(self.ui.refAddFilePushButton, SIGNAL('clicked()'), self.addRefFiles)
        self.connect(self.ui.refRemoveFilePushButton, SIGNAL('clicked()'), self.removeRefFiles)
        self.connect(self.ui.refRefFileListWidget, SIGNAL('itemSelectionChanged()'), self.updateSelectedRefFiles)
        self.connect(self.ui.refLegendCheckBox,SIGNAL('stateChanged(int)'), self.updateRefPlotData)
        self.connect(self.ui.refLegendLocComboBox,SIGNAL('currentIndexChanged(int)'), self.updateRefPlotData)
        #self.connect(self.ui.refCenCheckBox, SIGNAL('stateChanged(int)'),)
        self.connect(self.ui.refCenLineEdit, SIGNAL('returnPressed()'), self.refAnalyze)
        self.connect(self.ui.gixLineAreaPushButton, SIGNAL('clicked()'), self.updateCutData)
        self.connect(self.ui.gixIntRangeLineEdit, SIGNAL('returnPressed()'), self.updateCutData)
        self.connect(self.ui.pilLineAreaPushButton, SIGNAL('clicked()'), self.updateCutData)
        self.connect(self.ui.pilIntRangeLineEdit, SIGNAL('returnPressed()'), self.updateCutData)
        #self.connect(self.ui.gixCutDirComboBox, SIGNAL('currentIndexChanged(int)'), self.updateCutData)
        self.connect(self.ui.cutLogXCheckBox, SIGNAL('stateChanged(int)'),self.updateCutPlotData)
        self.connect(self.ui.cutLogYCheckBox, SIGNAL('stateChanged(int)'),self.updateCutPlotData)
        self.connect(self.ui.cutGridCheckBox, SIGNAL('stateChanged(int)'),self.updateCutPlotData)
        self.connect(self.ui.cutLegendCheckBox, SIGNAL('stateChanged(int)'),self.updateCutPlotData)
        self.connect(self.ui.cutLegendLocComboBox, SIGNAL('currentIndexChanged(int)'),self.updateCutPlotData)
        self.connect(self.ui.cutErrorbarCheckBox, SIGNAL('stateChanged(int)'),self.updateCutData)
        self.connect(self.ui.cutExportPushButton, SIGNAL('clicked()'),self.saveCutData)
        self.connect(self.ui.cutOffsetCheckBox, SIGNAL('stateChanged(int)'),self.updateCutData)
        self.connect(self.ui.cutOffsetLineEdit, SIGNAL('returnPressed()'),self.updateCutData)
        self.connect(self.ui.mcaExportPushButton, SIGNAL('clicked()'),self.saveMcaData)
        self.connect(self.ui.mcaOffsetCheckBox, SIGNAL('stateChanged(int)'),self.updateMcaPlotData)
        self.connect(self.ui.mcaOffsetLineEdit, SIGNAL('returnPressed()'),self.updateMcaPlotData)
        self.connect(self.ui.mcaFitPushButton, SIGNAL('clicked()'), self.mcaPeakFit)
        self.connect(self.ui.mcaFitAllPushButton, SIGNAL('clicked()'), self.mcaPeakFitAll)
     #   self.connect(self.ui.mcaSumPushButton, SIGNAL('clicked()'), self.mcaSum)
        self.connect(self.ui.mcaSumAllPushButton, SIGNAL('clicked()'), self.mcaSumAll)
        self.connect(self.ui.mcaAcceptPushButton, SIGNAL('clicked()'), self.mcaAcceptPeak)
     #   self.connect(self.ui.mcaFrameCheckBox, SIGNAL('stateChanged(int)'), self.updateMcaInt)
     #   self.connect(self.ui.mcaQzCheckBox, SIGNAL('stateChanged(int)'), self.updateMcaInt)
     #   self.connect(self.ui.mcaEnergyCheckBox, SIGNAL('stateChanged(int)'), self.updateMcaInt)
        self.connect(self.ui.mcaXAxisComboBox,SIGNAL('currentIndexChanged(int)'), self.updateMcaInt)
        self.connect(self.ui.mcaDelPushButton, SIGNAL('clicked()'), self.mcaDelIntData)
        self.connect(self.ui.mcaMergePushButton, SIGNAL('clicked()'), self.mcaIntMergeDia)
        self.connect(self.ui.mcaExportIntPushButton, SIGNAL('clicked()'), self.mcaSaveIntData)
        self.connect(self.ui.mcaAddFilePushButton, SIGNAL('clicked()'),self.addMcaFiles)
        self.connect(self.ui.mcaRemoveFilePushButton, SIGNAL('clicked()'),self.removeMcaFiles)
        self.connect(self.ui.mcaFileListWidget, SIGNAL('itemSelectionChanged()'),self.updateSelectedMcaFiles)
        self.connect(self.ui.cutClearGraphPushButton, SIGNAL('clicked()'),self.clearCutGraph)
        self.connect(self.ui.backgroundPushButton, SIGNAL('clicked()'),self.addBGImagestoList)
        self.connect(self.ui.removePushButton, SIGNAL('clicked()'), self.removeBGImages)
        self.connect(self.ui.backgroundListWidget, SIGNAL('itemSelectionChanged()'),self.bgSelectionChanged)
        self.connect(self.ui.bgSelectAllCheckBox, SIGNAL('stateChanged(int)'),self.bgSelectAll)
        self.connect(self.ui.bgFacLineEdit, SIGNAL('returnPressed()'),self.ccdSelectedScanChanged)
        self.connect(self.ui.scanListWidget, SIGNAL('doubleClicked(QModelIndex)'),self.displayScanInfo)
       # self.connect(self.ui.pilGIDPushButton, SIGNAL('clicked()'),self.pilGIDData)
       # self.connect(self.ui.pilGIDPatchPushButton, SIGNAL('clicked()'),self.pilGIDpatchData)
        self.connect(self.ui.gidComboBox, SIGNAL('activated(int)'),self.pilGID)
        self.connect(self.ui.gisaxsComboBox, SIGNAL('activated(int)'),self.pilGISAXS)
        self.connect(self.ui.pilMoviePushButton, SIGNAL('clicked()'),self.pilMovieShow)
        self.connect(self.ui.pilBurnPushButton, SIGNAL('clicked()'),self.pilBurnTest)
        self.connect(self.ui.pilGioxsPushButton,SIGNAL('clicked()'),self.gioxsDisplay)
        self.connect(self.ui.refNormFacCheckBox, SIGNAL('stateChanged(int)'),self.updateRefPlotData)
        self.connect(self.ui.pilAspectComboBox, SIGNAL('currentIndexChanged(int)'),self.update2dPlots)
        self.connect(self.ui.plotAddFilePushButton, SIGNAL('clicked()'),self.addPlotFile)
        self.connect(self.ui.plotRemoveFilePushButton, SIGNAL('clicked()'), self.removePlotFile)
        self.connect(self.ui.plotFileListWidget, SIGNAL('itemSelectionChanged()'),self.updateSelectedPlotFile)
        self.connect(self.ui.plotLogXCheckBox, SIGNAL('stateChanged(int)'),self.updatePlotPlot)
        self.connect(self.ui.plotLogYCheckBox, SIGNAL('stateChanged(int)'),self.updatePlotPlot)
        self.connect(self.ui.plotGridCheckBox, SIGNAL('stateChanged(int)'),self.updatePlotPlot)
        self.connect(self.ui.plotLegendCheckBox, SIGNAL('stateChanged(int)'),self.updatePlotPlot)
        self.connect(self.ui.plotLegendLocComboBox, SIGNAL('currentIndexChanged(int)'),self.updatePlotPlot)
        self.connect(self.ui.onedLabelSizeSpinBox,SIGNAL('valueChanged(int)'),self.updatePlotPlot)
        self.connect(self.ui.onedTickSizeSpinBox,SIGNAL('valueChanged(int)'),self.updatePlotPlot)
        self.connect(self.ui.onedStyComboBox,SIGNAL('currentIndexChanged(int)'),self.updatePlotPlot)
        self.connect(self.ui.plotRSPushButton, SIGNAL('clicked()'), self.setPlotScale)
        self.connect(self.ui.plotPFPushButton, SIGNAL('clicked()'), self.dataPeakFit)
        self.connect(self.ui.plotMoveUpPushButton, SIGNAL('clicked()'), self.upPlotFile)
        self.connect(self.ui.plotMoveDownPushButton, SIGNAL('clicked()'), self.downPlotFile)
        self.connect(self.ui.graPushButton, SIGNAL('clicked()'),self.graPlotPlot)
        self.connect(self.ui.savegraPushButton, SIGNAL('clicked()'),self.saveGraPlotData)
        self.connect(self.ui.twodPlotAddFilePushButton, SIGNAL('clicked()'),self.add2dPlotFile)
        self.connect(self.ui.twodPlotRemoveFilePushButton, SIGNAL('clicked()'), self.remove2dPlotFile)
        self.connect(self.ui.twodPlotFileListWidget, SIGNAL('itemSelectionChanged()'),self.updateSelected2dPlotFile)
        self.connect(self.ui.twodCMapComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlotPlot)
        self.connect(self.ui.twodLogIntCheckBox, SIGNAL('stateChanged(int)'),self.update2dPlotPlot)
        self.connect(self.ui.twodMinHorizontalSlider, SIGNAL('sliderReleased()'),self.update2dPlotPlot)
        self.connect(self.ui.twodMaxHorizontalSlider, SIGNAL('sliderReleased()'),self.update2dPlotPlot)
        self.connect(self.ui.twodMinLineEdit, SIGNAL('returnPressed()'),self.update2dPlotPlot)
        self.connect(self.ui.twodMaxLineEdit, SIGNAL('returnPressed()'),self.update2dPlotPlot)
        self.connect(self.ui.twodUpdatePlotPushButton, SIGNAL('clicked()'),self.update2dIntPolData)
        self.connect(self.ui.twodIntPolComboBox, SIGNAL('currentIndexChanged(int)'),self.update2dPlotPlot)
        self.connect(self.ui.twodLabelSizeSpinBox,SIGNAL('valueChanged(int)'),self.update2dPlotPlot)
        self.connect(self.ui.twodTickSizeSpinBox,SIGNAL('valueChanged(int)'),self.update2dPlotPlot)
        self.connect(self.ui.twodColBarCheckBox,SIGNAL('stateChanged(int)'),self.update2dPlotPlot)
        self.connect(self.ui.twodMinHorizontalSlider, SIGNAL('sliderReleased()'), self.update2dMinSlider)
        self.connect(self.ui.twodMaxHorizontalSlider, SIGNAL('sliderReleased()'), self.update2dMaxSlider)
        
    def selectAPS_15IDC(self):
        self.ui.statusBar.clearMessage()
        self.beamline='APS-15IDC'
        self.ui.actionAPS_9IDC.setChecked(False)
        self.ui.statusBar.showMessage('APS-15IDC Selected')
    
    def selectAPS_9IDC(self):
        self.ui.statusBar.clearMessage()
        self.beamline='APS-9IDC'
        self.ui.actionAPS_15IDC.setChecked(False)
        self.ui.statusBar.showMessage('APS-9IDC Selected')
        
    def showAbout(self):
        cwd=os.getcwd()
        files=['mainwindow.py','mainwindow.ui','main.py','mca_routines.py','mpl2dwidget.py','mplwidget.py','spec_routines.py','TwoDDetector.py','logo.png']
        fname=[cwd+'/'+fname for fname in files]
        updateTime=max([os.path.getmtime(fn) for fn in fname])
        self.messageBox('LSS-Reader\n Version: 17.01\nLast Update: '+time.strftime("%m/%d/%Y %I:%M:%S %p",time.localtime(updateTime))+'\nCopyright belongs to:\n\tWei Bu <weibu1977@gmail.com>\n\tMrinal K Bera <nayanbera@gmail.com>',title='About')
        
    def openSpecFile(self):
        self.ui.statusBar.clearMessage()       
        self.ui.imagesLabel.setText('Detector:')
        self.specFileName=QFileDialog.getOpenFileName(caption='Open Spec File')
        if self.specFileName!='':
            self.clearAll()
            self.directory=str(QFileInfo(self.specFileName).absolutePath())
            self.directory_old=self.directory
        else:
            self.directory=self.directory_old
        specData_old=self.specData
        specPar_old=self.specPar
        startLineNum_old=self.startLineNum
        try:
            self.specData={}
            self.specPar={}
            self.startLineNum=0
            self.endScanNum=0
            self.readSpecFile()
            self.ui.specFileLabel.setText('SpecFile: '+str(self.specFileName.split('/')[-1]))
        except:
            self.specData=specData_old
            self.specPar=specPar_old
            self.startLineNum=startLineNum_old
            print "Either you have just loaded the old file or there is something wrong with Specfile reading"

        
    def clearAll(self):
        self.disconnect(self.ui.scanListWidget, SIGNAL('itemSelectionChanged()'),self.scanListChanged)
        self.disconnect(self.ui.imageListWidget, SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged) 
        self.disconnect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
        self.disconnect(self.ui.spXComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.disconnect(self.ui.spYComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.disconnect(self.ui.spY2ComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.disconnect(self.ui.spNComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.ui.imageListWidget.clear()
        self.ui.statusBar.clearMessage()
        self.ui.scanListWidget.clear()
        self.ui.specPlotMplWidget.canvas.ax.clear()
        self.ui.spXComboBox.clear()
        self.ui.spYComboBox.clear()
        self.ui.spY2ComboBox.clear()
        self.ui.spNComboBox.clear()
        self.ui.mcaPlotMplWidget.canvas.ax.clear()
        self.ui.mcaPlotMplWidget.canvas.draw()
        self.ui.refBadPixListWidget.clear()
        self.ui.cutPlotMplWidget.canvas.ax.clear()
        self.ui.cutPlotMplWidget.canvas.draw()
        self.ui.gixMplWidget.canvas.ax.clear()
        self.ui.gixMplWidget.canvas.draw()
        self.ui.refADDataPlotWidget.canvas.ax.clear()
        self.ui.refADDataPlotWidget.canvas.draw()
        self.ui.refQzListWidget.clear() 
        self.connect(self.ui.scanListWidget, SIGNAL('itemSelectionChanged()'),self.scanListChanged)
        self.connect(self.ui.imageListWidget, SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged) 
        self.connect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
        self.connect(self.ui.spXComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spYComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spY2ComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.connect(self.ui.spNComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
         
        
    def readSpecFile(self):
        '''
        Reading spec File
        '''
        self.specRead=specread(self.specFileName,startLineNum=self.startLineNum, endScanNum=self.endScanNum, beamline=self.beamline,data=self.specData,par=self.specPar)
        self.clearAll()
        if self.startLineNum==0:
            inipil=0
        else:
            inipil=1
        self.startLineNum=self.specRead.endLineNum
        self.specData=self.specRead.Data
        self.specPar=self.specRead.Par
        self.endScanNum=self.specData['NumOfScans']
       # print self.endScanNum
        if self.specData['NumOfScans']==0:
            self.ui.statusBar.showMessage(self.SpecData['Message'])
        else:
            if self.specData['Error']:
                self.messageBox('Warning:: The spec file has identical scan numbers!!')
            self.scanlines=[self.specData[i]['ScanLine'] for i in list(set(self.specData['ScanNums']))]
            self.ui.scanListWidget.addItems(self.scanlines)
            if np.any([self.specPar[i]['Detector']=='Vortex' for i in range(1,self.specData['NumOfScans']+1)]):
                if QDir(self.directory+'/'+self.mcaDirName).exists()==True:
                    self.mcaDir=self.directory+'/'+self.mcaDirName
                else:
                    self.messageBox('Warning:: The default MCA directory not found, please choose!!')
                    self.mcaDir=QFileDialog.getExistingDirectory(caption='Open MCA Directory')
                self.mcafhead=str(self.mcaDir)+'/'+self.specFileName.split('/')[-1]+'_'
                self.mcaftail='_mca' 
                print "MCA data Found!!"
            if np.any([self.specPar[i]['Detector']=='Bruker' for i in range(1,self.specData['NumOfScans']+1)]):
                try:
                    self.ccdfhead=str(self.ccdDir)+'/'+self.specFileName.split('/')[-1]+'_SCAN-'
                    self.ccdftail='.sfrm'
                except:
                    if QDir(self.directory+'/'+self.ccdDirName).exists()==True:
                        self.ccdDir=self.directory+'/'+self.ccdDirName
                    else:
                        self.messageBox('Warning:: The default CCD directory not found, please choose!!')
                        self.ccdDir=QFileDialog.getExistingDirectory(caption='Open CCD Directory')
                    self.ccdfhead=str(self.ccdDir)+'/'+self.specFileName.split('/')[-1]+'_SCAN-'
                    self.ccdftail='.sfrm'
                print "Bruker Images Found!!"
            if np.any([self.specPar[i]['Detector']=='Pilatus' for i in range(1,self.specData['NumOfScans']+1)]):
                if QDir(self.directory+'/'+self.pilDirName).exists()==True:
                    self.pilDir=self.directory+'/'+self.pilDirName
                else:
                    self.messageBox('Warning:: The default Piltus directory not found, please choose!!')
                    self.pilDir=QFileDialog.getExistingDirectory(caption='Open Pilatus Directory')
                self.pilfhead=str(self.pilDir)+'/'+self.specFileName.split('/')[-1]+'_SCAN-'
                self.pilftail='.tif'
                print "Pilatus Images Found!!"
                if inipil==0:
                    self.ui.refBGOffLineEdit.setText('0.5')   #preset the background    
                    self.ui.refSlitLineEdit.setText('21,11')  #preset the size of ROI
                    self.ui.refQcLineEdit.setText('0.0217') #preset the critical q
                    self.ui.AbsFacLineEdit.setText('1.42') #for Al 50um at 10keV
        self.ui.statusBar.showMessage('Done')
        self.scanCenter={}


    def scanListChanged(self):
        self.ui.statusBar.clearMessage()
        self.ui.imageSelectAllCheckBox.setCheckState(0)
        self.ui.imageListWidget.clear()
        self.disconnect(self.ui.spXComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.disconnect(self.ui.spYComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.disconnect(self.ui.spY2ComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.disconnect(self.ui.spNComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
        self.selectedScans=self.ui.scanListWidget.selectedItems()
        #self.selectedScanNums=[self.ui.scanListWidget.row(items) for items in self.selectedScans]
        self.selectedScanNums=[int(str(items.text()).split()[1]) for items in self.selectedScans]
        self.ui.scansLineEdit.setText(str([item for item in self.selectedScanNums])[1:-1])
        if self.checkSameScans()==False or self.selectedScanNums==[]:
            self.ui.statusBar.showMessage('Error:: The scans are not identical or some scans have no data!!')
            self.disconnect(self.ui.scanListWidget, SIGNAL('itemSelectionChanged()'),self.scanListChanged)
            for item in self.selectedScans:
                self.ui.scanListWidget.setItemSelected(item,False)
            self.connect(self.ui.scanListWidget, SIGNAL('itemSelectionChanged()'),self.scanListChanged)
            self.ui.specPlotMplWidget.canvas.ax.clear()
            self.ui.specPlotMplWidget.canvas.draw()
        else:
            self.ui.spXComboBox.clear()
            self.ui.spYComboBox.clear()
            self.ui.spY2ComboBox.clear()
            self.ui.spNComboBox.clear()
            self.ui.spXComboBox.addItems(self.specData[self.selectedScanNums[0]]['ScanVar'])
            self.ui.spYComboBox.addItems(self.specData[self.selectedScanNums[0]]['ScanVar'])
            self.ui.spY2ComboBox.addItems(['None'])
            self.ui.spY2ComboBox.addItems(self.specData[self.selectedScanNums[0]]['ScanVar'])
         #   ycol=self.ui.spYComboBox.findText(self.specData['YCol'])
         #   if ycol>0:
         #       self.ui.spYComboBox.setCurrentIndex(self.ui.spYComboBox.findText(self.specData['YCol']))
         #   else:
         #       self.ui.spYComboBox.setCurrentIndex(len(self.specData[self.selectedScanNums[0]]['ScanVar'])-1)
            self.ui.spNComboBox.addItems(self.specData[self.selectedScanNums[0]]['ScanVar'])
         #   ncol=self.ui.spYComboBox.findText(self.specData['NCol'])
         #   if ncol>0:
         #       self.ui.spNComboBox.setCurrentIndex(self.ui.spNComboBox.findText(self.specData['NCol']))
         #   else:
         #       self.ui.spNComboBox.setCurrentIndex(len(self.specData[self.selectedScanNums[0]]['ScanVar'])-2)   
            if self.specData[self.selectedScanNums[0]]['ScanVar'][0]=='H':
                self.ui.spXComboBox.setCurrentIndex(2)
            else: 
                self.ui.spXComboBox.setCurrentIndex(0)
            self.ui.spYComboBox.setCurrentIndex(len(self.specData[self.selectedScanNums[0]]['ScanVar'])-1)
            self.ui.spY2ComboBox.setCurrentIndex(0)
            self.ui.spNComboBox.setCurrentIndex(len(self.specData[self.selectedScanNums[0]]['ScanVar'])-2)
            self.connect(self.ui.spXComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
            self.connect(self.ui.spYComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
            self.connect(self.ui.spY2ComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
            self.connect(self.ui.spNComboBox, SIGNAL('currentIndexChanged(int)'),self.updateSpecPlotData)
            try:
                self.updateSpecPlotData()
            except:
                self.ui.statusBar.showMessage('Warning:: The scan(s) are not regular Spec scan')
            if self.specPar[self.selectedScanNums[0]]['Detector']=='Vortex':
                self.updateMcaImageList()
                self.ui.imagesLabel.setText('Vortex:')
                self.ui.PlotWidget.setCurrentIndex(0)
            if self.specPar[self.selectedScanNums[0]]['Detector']=='Bruker':
                self.updateCcdImageList()
                self.ui.PlotWidget.setCurrentIndex(0)
                self.ui.imagesLabel.setText('Bruker:')
            if self.specPar[self.selectedScanNums[0]]['Detector']=='Pilatus':
                self.updatePilImageList()
                self.ui.PlotWidget.setCurrentIndex(0)
                self.ui.imagesLabel.setText('Pilatus:')
        self.ui.statusBar.showMessage('Done')
            
    def scanListInputChanged(self):
        inputscannumbers=str(self.ui.scansLineEdit.text()).split(',')
        self.disconnect(self.ui.scanListWidget, SIGNAL('itemSelectionChanged()'),self.scanListChanged)
        self.selectedScans=self.ui.scanListWidget.selectedItems()
        self.selectedScanNums=[int(str(items.text()).split()[1]) for items in self.selectedScans]
        try:    
            for i in self.selectedScanNums:
                self.ui.scanListWidget.setItemSelected(self.ui.scanListWidget.item(i-1),False)
            self.selectedScanNums=[]
        except:
            self.selectedScanNums=[]
        for item in inputscannumbers:
            scan=map(int, item.split('-'))
            if len(scan)>1:
                self.selectedScanNums=self.selectedScanNums+range(scan[0],scan[1]+1)
            else:
                self.selectedScanNums=self.selectedScanNums+scan
        #self.selectedScanNums=[item for item in self.selectedScanNums]
        self.connect(self.ui.scanListWidget, SIGNAL('itemSelectionChanged()'),self.scanListChanged)
        snum=[int(item.split()[1]) for item in self.scanlines]
        for i in self.selectedScanNums:
            self.ui.scanListWidget.setItemSelected(self.ui.scanListWidget.item(snum.index(i)),True)
        self.ui.scanListWidget.setCurrentRow(self.selectedScanNums[-1]-1)
        self.ui.statusBar.showMessage('Done')
        

    def updateSpecPlotData(self):
        self.ui.statusBar.clearMessage()
        self.ui.specPlotMplWidget.canvas.ax.clear()
        self.specData['YCol']=str(self.ui.spYComboBox.currentText())
        self.specData['NCol']=str(self.ui.spNComboBox.currentText())
        title='File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1]
        if self.checkSameScans()==False:
            self.ui.statusBar.showMessage('Error:: The scans are not identical or some scans have no data!!')
        else:
            for i in self.selectedScanNums:
                x=self.specData[i][str(self.ui.spXComboBox.currentText())]
                y=self.specData[i][str(self.ui.spYComboBox.currentText())]
                n=self.specData[i][str(self.ui.spNComboBox.currentText())]
                if str(self.ui.spY2ComboBox.currentText())=='None':
                    self.label='S '+str(i)
                    self.specPlot(x,y,n)
                else: 
                    self.label='S '+str(i)+str(self.ui.spYComboBox.currentText())
                    self.specPlot(x,y,n)
                    y2=self.specData[i][str(self.ui.spY2ComboBox.currentText())]
                    self.label='S '+str(i)+str(self.ui.spY2ComboBox.currentText())
                    self.specPlot(x,y2,n)
            yerr=pl.sqrt(np.abs(y))
            if self.ui.spNCheckBox.checkState()!=0:
                yerr=pl.sqrt(y/n**2+y**2/n**3)
                y=y/n
            try:
                self.specPeakFit(x,y,yerr)
                self.peakPos='%.4f'%x[np.argmax(y)]
                if (self.specFitPar[0]+self.specFitPar[1]/2.0)>x[0] and (self.specFitPar[0]+self.specFitPar[1]/2.0)>x[-1]:
                    self.peakCenter='%.4f'%(self.specFitPar[0]-self.specFitPar[1]/2.0)
                    self.peakFWHM='%.4f'%self.specFitPar[5]
                    #print '1', self.specFitPar[0]
                    self.ui.specPlotMplWidget.canvas.ax.set_title(title+'\n'+'Peak= '+self.peakPos+', Edge= '+self.peakCenter+', Edge Width= '+self.peakFWHM)
                elif (self.specFitPar[0]-self.specFitPar[1]/2.0)<x[0] and (self.specFitPar[0]-self.specFitPar[1]/2.0)<x[-1]:
                    self.peakCenter='%.4f'%(self.specFitPar[0]+self.specFitPar[1]/2.0)
                    self.peakFWHM='%.4f'%self.specFitPar[5]
                    #print '2', self.specFitPar[0]
                    self.ui.specPlotMplWidget.canvas.ax.set_title(title+'\n'+'Peak= '+self.peakPos+', Edge= '+self.peakCenter+', Edge Width= '+self.peakFWHM)
                else:    
                    self.peakCenter='%.4f'%self.specFitPar[0]
                    self.peakFWHM='%.4f'%self.specFitPar[1]
                    #print '3', self.specFitPar[0]
                    self.ui.specPlotMplWidget.canvas.ax.set_title(title+'\n'+'Peak= '+self.peakPos+', Cen= '+self.peakCenter+', FWHM= '+self.peakFWHM)
                if self.specFitPar[7]<1:
                    self.ui.specPlotMplWidget.canvas.ax.plot(x,self.specPeakFun1(self.specFitPar, x),'g--',label='Fit')
                else:
                    self.ui.specPlotMplWidget.canvas.ax.plot(x,self.specPeakFun2(self.specFitPar, x),'g--',label='Fit')
                self.scanCenter[self.selectedScanNums[-1]]=float(format(self.specFitPar[0],'.3f'))
            except:
                self.peakPos='%.4f'%x[np.argmax(y)]
                self.peakCenter=self.peakPos
                self.peakFWHM='0'
        self.specPlotSettings()
        self.ui.statusBar.showMessage('Done')
        
    def specSaveData(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save SPEC scans', directory=self.directory))
        for i in self.selectedScanNums:
            x=self.specData[i][str(self.ui.spXComboBox.currentText())]
            y=self.specData[i][str(self.ui.spYComboBox.currentText())]
            n=self.specData[i][str(self.ui.spNComboBox.currentText())]
            yerr=pl.sqrt(y)
            self.fname=self.saveFileName+'_'+str(i)+'_'+str(self.ui.spXComboBox.currentText())+'_'+str(self.ui.spYComboBox.currentText())+'.txt'
            if self.ui.spNCheckBox.checkState()!=0:
                yerr=pl.sqrt(y/n**2+y**2/n**3)
                y=y/n
                self.fname=self.saveFileName+'_'+str(i)+'_'+str(self.ui.spXComboBox.currentText())+'_'+str(self.ui.spYComboBox.currentText())+'_'+str(self.ui.spNComboBox.currentText())+'.txt'
            specscandata=[[x[j],y[j],yerr[j]] for j in range(len(x))]
            np.savetxt(self.fname,specscandata,fmt='%.4f\t%.4e\t%.4e')
            if str(self.ui.spY2ComboBox.currentText())!='None':
                y2=self.specData[i][str(self.ui.spY2ComboBox.currentText())]
                yer2r=pl.sqrt(y2)
                self.fname=self.saveFileName+'_'+str(i)+'_'+str(self.ui.spXComboBox.currentText())+'_'+str(self.ui.spY2ComboBox.currentText())+'.txt'
                if self.ui.spNCheckBox.checkState()!=0:
                    yerr2=pl.sqrt(y2/n**2+y2**2/n**3)
                    y2=y2/n
                    self.fname=self.saveFileName+'_'+str(i)+'_'+str(self.ui.spXComboBox.currentText())+'_'+str(self.ui.spY2ComboBox.currentText())+'_'+str(self.ui.spNComboBox.currentText())+'.txt'
                specscandata=[[x[j],y2[j],yerr[j]] for j in range(len(x))]
                np.savetxt(self.fname,specscandata,fmt='%.4f\t%.4e\t%.4e')
            
            
    def specPlotSettings(self):
        self.spLogX=self.ui.spLogXCheckBox.checkState()
        self.spLogY=self.ui.spLogYCheckBox.checkState()
        self.spGrid=self.ui.spGridCheckBox.checkState()
        self.ui.specPlotMplWidget.canvas.ax.set_xscale('linear')
        self.ui.specPlotMplWidget.canvas.ax.set_yscale('linear')
        self.ui.specPlotMplWidget.canvas.ax.set_xlabel(str(self.ui.spXComboBox.currentText()))
        self.ui.specPlotMplWidget.canvas.ax.set_ylabel(str(self.ui.spYComboBox.currentText()))
        if str(self.ui.spY2ComboBox.currentText())=='None':
            self.ui.specPlotMplWidget.canvas.ax.set_ylabel(str(self.ui.spYComboBox.currentText()))
        else:
            self.ui.specPlotMplWidget.canvas.ax.set_ylabel(str(self.ui.spYComboBox.currentText())+' and '+str(self.ui.spY2ComboBox.currentText()))
        self.ui.specPlotMplWidget.canvas.ax.grid(b=False)
        if self.spLogX!=0:
            self.ui.specPlotMplWidget.canvas.ax.set_xscale('log')
        if self.spLogY!=0:
            self.ui.specPlotMplWidget.canvas.ax.set_yscale('log')
        if self.spGrid!=0:
            self.ui.specPlotMplWidget.canvas.ax.grid(b=True,color='r',linestyle='--')
        if self.ui.spLegendCheckBox.checkState()!=0:
            self.ui.specPlotMplWidget.canvas.ax.legend(loc=self.ui.spLegendLocComboBox.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        self.ui.specPlotMplWidget.canvas.draw()
            
        
    def specPlot(self,x,y,n):
        yerr=pl.sqrt(y)
        if self.ui.spNCheckBox.checkState()!=0:
            yerr=pl.sqrt(y/n**2+y**2/n**3)
            y=y/n
        self.ui.specPlotMplWidget.canvas.ax.errorbar(x,y,yerr,label=self.label,fmt='o-')
        for i in range(len(x)):
            self.ui.specPlotMplWidget.canvas.ax.text(x[i],y[i],str(i),color='r',fontsize=16)
        
    def specPeakFit(self,x,y,yerr):
        par=[x[np.argmax(y)],np.abs(x[-1]-x[0])/2.0,np.max(y),0.0,0.0,np.abs(x[0]-x[1]),0.0,1]
        p=leastsq(self.specPeakRes1,par,args=(x,y,yerr),maxfev=5000)
        if p[0][1]<np.abs(x[0]-x[1]):
            p=leastsq(self.specPeakRes2,par,args=(x,y,yerr),maxfev=5000)
        self.specFitPar=p[0]
        
    
    def specPeakRes1(self,par,x,y,yerr):
        if np.all(np.where(yerr>0,1,0))==True:
            return (y-self.specPeakFun1(par,x))/yerr
        else:
            return (y-self.specPeakFun1(par,x))
            
    def specPeakRes2(self,par,x,y,yerr):
        if np.all(np.where(yerr>0,1,0))==True:
            return (y-self.specPeakFun2(par,x))/yerr
        else:
            return (y-self.specPeakFun2(par,x))
            
        
    def specPeakFun1(self, par, x):
        par[1]=np.abs(par[1])
        par[5]=np.abs(par[5])
        par[7]=-1
        return par[2]*((1+np.tanh((x-par[0]+par[1]/2)/par[5]))/2+(np.tanh(-(x-par[0]-par[1]/2)/par[5])+1)/2)+par[3]+par[4]*(x-par[0])+par[6]*(x-par[0])**2
    
    def specPeakFun2(self, par, x):
        par[1]=np.abs(par[1])
        sig=par[1]/2.0/np.sqrt(2.0*np.log(2.0))
        par[7]=1
        return par[2]*np.exp(-(x-par[0])**2/2.0/sig**2)+par[3]+par[4]*x
        

        
    def checkSameScans(self):
        scanvar=[self.specData[i]['ScanVar'] for i in self.selectedScanNums]
        return all(x==scanvar[0] for  x in scanvar)
    
    def checkSameQzs(self,qzvalues):
        return all(np.abs(x-qzvalues[0])<0.0001 for x in qzvalues)
        
    def checkSameArrays(self,anyarray):
        return all(x==anyarray[0] for x in anyarray)
        
    def updateMcaImageList(self):
        self.disconnect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
        self.ui.statusBar.clearMessage()
        for i in self.selectedScanNums:
            self.ui.scanListWidget.setItemSelected(self.ui.scanListWidget.item(i-1),True)
        self.mcaFileNames=[self.mcafhead+str(i)+self.mcaftail for i in self.selectedScanNums]
        self.mcaData={}
        self.mcaPar={}
        start=0 
        for i in range(len(self.mcaFileNames)):  
            self.mcaread=mcaread(self.mcaFileNames[i],beamline=self.beamline)
            data=self.mcaread.Data
            par=self.mcaread.Par
            for j in range(data['NumOfScans']):                
                self.mcaData[start]=data[j]
                self.mcaPar[start]=par[j]
                self.ui.imageListWidget.addItem('S# '+str(self.selectedScanNums[i])+'\tmS# '+str(j+1)+'\tQz='+str(self.mcaPar[start]['Q'][2]))
                start=start+1           
        self.connect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
        
    def mcaSelectedScanChanged(self):
        self.ui.statusBar.clearMessage()
        self.selectedMcaScans=self.ui.imageListWidget.selectedItems()
        self.selectedMcaScanNums=[self.ui.imageListWidget.row(items) for items in self.selectedMcaScans]
        self.mcafitplotstatus=0
        try:
            self.ui.PlotWidget.setCurrentIndex(5)
            self.updateMcaPlotData()
        except:
            self.ui.statusBar.showMessage('Warning:: Some scans maybe missing')
            
            
            
    def updateMcaPlotData(self):
        self.ui.statusBar.clearMessage()
        self.ui.mcaPlotMplWidget.canvas.ax.clear()
        self.nomMcaData={}
        fact=1.0

        if len(self.selectedMcaScanNums) > 10:
            self.progressDialog = QProgressDialog('Reading MCA Frames', 'Abort', 0, 100)
            self.progressDialog.setWindowModality(Qt.WindowModal)
            self.progressDialog.setWindowTitle('Wait')
            self.progressDialog.setAutoClose(True)
            self.progressDialog.setAutoReset(True)
            self.progressDialog.setMinimum(1)
            self.progressDialog.setMaximum(len(self.selectedMcaScanNums))
            self.progressDialog.show()

        for i in self.selectedMcaScanNums:
            y=pl.array(self.mcaData[i]['Vortex'],dtype='float')
            x=pl.arange(1,len(y)+1)
            if self.ui.mcaCalibCheckBox.checkState()!=0:
                self.ui.mcaCalibConLineEdit.setText(str(self.mcaPar[i]['Calib'][0]))
                self.ui.mcaCalibLinLineEdit.setText(str(self.mcaPar[i]['Calib'][1]))
                self.ui.mcaCalibQuaLineEdit.setText(str(self.mcaPar[i]['Calib'][2]))
            con=float(self.ui.mcaCalibConLineEdit.text())
            lin=float(self.ui.mcaCalibLinLineEdit.text())
            qua=float(self.ui.mcaCalibQuaLineEdit.text())
            x=con+lin*x+qua*x**2
            if str(self.ui.mcaNormComboBox.currentText())!='None':
                monitor=str(self.ui.mcaNormComboBox.currentText())
                n=self.mcaData[i][monitor]
                tc=self.mcaPar[i]['Time'][1]/self.mcaPar[i]['Time'][2]  #real count time correction
                self.nomMcaData[i]=np.vstack((x,y/n*tc,pl.sqrt(y+y**2/n)/n*tc)).transpose()
            else:
                n=1
                tc=1
                self.nomMcaData[i]=np.vstack((x,y,pl.sqrt(y))).transpose()
            self.mcaLabel=str(i+1)+' Qz='+'%.4f'%self.mcaPar[i]['Q'][2]
            if self.ui.mcaOffsetCheckBox.checkState()!=0:
                fact=fact*float(self.ui.mcaOffsetLineEdit.text())
            self.mcaPlot(x,y,n,tc,fact)
            if len(self.selectedMcaScanNums) > 10:
                self.progressDialog.setLabelText('Reading MCA Frame #' + str(i))
                self.updateProgress()
                if self.progressDialog.wasCanceled() == True:
                    break
        if len(self.selectedMcaScanNums) > 10:
            self.progressDialog.hide()
        self.mcaPlotSettings()
        self.ui.statusBar.showMessage('Done')
        
    def mcaPlot(self,x,y,n,tc,fact):
        if str(self.ui.mcaNormComboBox.currentText())!='None':
            yerr=pl.sqrt(y+y**2/n)/n*tc
        else:
            yerr=pl.sqrt(y)
        y=y/n*tc
        self.ui.mcaPlotMplWidget.canvas.ax.errorbar(x,fact*y,fact*yerr,label=self.mcaLabel,fmt='o-')

        
    def mcaPlotSettings(self):
        self.mcaLogX=self.ui.mcaLogXCheckBox.checkState()
        self.mcaLogY=self.ui.mcaLogYCheckBox.checkState()
        self.mcaGrid=self.ui.mcaGridCheckBox.checkState()
        self.ui.mcaPlotMplWidget.canvas.ax.set_xscale('linear')
        self.ui.mcaPlotMplWidget.canvas.ax.set_yscale('linear')
        self.ui.mcaPlotMplWidget.canvas.ax.set_xlabel('Energy')
        self.ui.mcaPlotMplWidget.canvas.ax.set_ylabel('Intensity')
        self.ui.mcaPlotMplWidget.canvas.ax.set_title('Fluorescence Spectrum')
        #self.ui.specPlotMplWidget.canvas.ax.legend('_nolegend_')
        self.ui.mcaPlotMplWidget.canvas.ax.grid(b=False)
        if self.mcaLogX!=0:
            self.ui.mcaPlotMplWidget.canvas.ax.set_xscale('log')
        if self.mcaLogY!=0:
            self.ui.mcaPlotMplWidget.canvas.ax.set_yscale('log')
        if self.mcaGrid!=0:
            self.ui.mcaPlotMplWidget.canvas.ax.grid(b=True,color='r',linestyle='--')
        if self.ui.mcaLegendCheckBox.checkState()!=0:
            self.ui.mcaPlotMplWidget.canvas.ax.legend(loc=self.ui.mcaLegendLocComboBox.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        if self.mcafitplotstatus==1:
             fit=np.array(self.peakfitdata)
             xran=np.abs(fit[:,0][-1]-fit[:,0][0])
             yran=np.max(fit[:,1])-np.min(fit[:,1])
             self.ui.mcaPlotMplWidget.canvas.ax.plot(fit[:,0],fit[:,1],'r-')
             self.ui.mcaPlotMplWidget.canvas.ax.set_xlim(fit[:,0][0]-0.25*xran,fit[:,0][-1]+0.25*xran)
             self.ui.mcaPlotMplWidget.canvas.ax.set_ylim(np.min(fit[:,1])-0.2*yran,np.max(fit[:,1])+0.2*yran)
        self.ui.mcaPlotMplWidget.canvas.draw()
            
    def saveMcaData(self):
        self.saveMcaFileName=str(QFileDialog.getSaveFileName(caption='Save Mca data',directory=self.directory))
        for i in self.selectedMcaScanNums:
            self.fmcaName=self.saveMcaFileName+str(self.ui.imageListWidget.item(i).text().split('\t')[0])+'_'+'%.4f'%self.mcaPar[i]['Q'][2]+'_mca.txt'
            np.savetxt(self.fmcaName,self.nomMcaData[i],fmt='%.4f\t%.4e\t%.4e')
            
    def mcaPeakFit(self):  #fit the mca spectrum
        self.mcafitstatus=1 #let self.peakfit know it's mca fit
        self.mcafitallstatus=0
        self.peakFit()
        
    def mcaPeakFitAll(self): #fit all selected spectrum
        self.mcafitallstatus=1
        self.peakFit()
 
    # def mcaSum(self):  # sum one selected spectrum with given range
    #     data=self.nomMcaData[self.selectedMcaScanNums[0]]
    #     ini=max(float(str(self.ui.mcaRanLineEdit.text()).split(':')[0]),data[0][0])
    #     fin=min(float(str(self.ui.mcaRanLineEdit.text()).split(':')[1]),data[-1][0])
    #     print ini, fin
    #     dataran=np.where((np.logical_and(data[0:]<=fin,data[0:]>=ini)))
    #     print dataran[0]
    #     newdata=data[dataran[0][0]:dataran[0][-1]+1]
    #     print newdata
    #     print np.sum(newdata[:,1]), np.sqrt(np.sum(newdata[:,2]**2))
    
    def mcaSumAll(self):  #sum all selected spectrum with given range
        self.mcaIntData=[]
        self.ui.mcaDataListWidget.clear()
        for i in range(len(self.selectedMcaScanNums)):
            data=self.nomMcaData[self.selectedMcaScanNums[i]]
            ini=max(float(str(self.ui.mcaRanLineEdit.text()).split(':')[0]),data[0][0])
            fin=min(float(str(self.ui.mcaRanLineEdit.text()).split(':')[1]),data[-1][0])     
            dataran=np.where((np.logical_and(data[:,0]<=fin,data[:,0]>=ini)))
            newdata=data[dataran[0][0]:dataran[0][-1]+1]
            frame=self.selectedMcaScanNums[i]+1
            qz=self.mcaPar[self.selectedMcaScanNums[i]]['Q'][2]
            try:
                energy=self.mcaPar[self.selectedMcaScanNums[i]]['Energy']
            except:
                energy=10.0
            inten=np.sum(newdata[:,1])
            error=np.sqrt(np.sum(newdata[:,2]**2))
            #print newdata, inten, error
            self.mcaIntData.append([frame, qz, energy, inten, error])
            if self.ui.mcaXAxisComboBox.currentIndex()==0:
                string=str(int(frame))+'\t'+'%.4e'%inten+'\t'+'%.4e'%error 
            elif self.ui.mcaXAxisComboBox.currentIndex()==1:
                string='%.4f'%qz+'\t'+'%.4e'%inten+'\t'+'%.4e'%error
            else:
                string='%.4f'%energy+'\t'+'%.4e'%inten+'\t'+'%.4e'%error
            self.ui.mcaDataListWidget.addItem(string)
        self.updateMcaIntPlot()
        self.command='Flu Sum, scans=['+str([item for item in np.sort(self.selectedScanNums)])[1:-1]+'], energy range=['+str(ini)+':'+str(fin)+']'
        self.ui.commandLineEdit.setText(self.command)
    
    
    def mcaAcceptPeak(self,num=0): #accecpt the fitting result, save the data, and update the mca integral plot 
        frame=self.selectedMcaScanNums[num]+1
        qz=self.mcaPar[self.selectedMcaScanNums[num]]['Q'][2]
        try:
            energy=self.mcaPar[self.selectedMcaScanNums[num]]['Energy']
        except:
            energy=10
       # print self.fitpeakpara
        inten=self.fitpeakpara[0][1]*self.fitpeakpara[0][2]
        error=np.sqrt(self.fitpeakpara[0][1]**2*self.fiterror[2]**2+self.fitpeakpara[0][2]**2*self.fiterror[1]**2)
       # print frame, qc, energy, inten, error
        self.mcaIntData.append([frame, qz, energy, inten, error])
        if self.ui.mcaXAxisComboBox.currentIndex()==0:
            string=str(int(frame))+'\t'+'%.4e'%inten+'\t'+'%.4e'%error 
        elif self.ui.mcaXAxisComboBox.currentIndex()==1:
            string='%.4f'%qz+'\t'+'%.4e'%inten+'\t'+'%.4e'%error
        else:
            string='%.4f'%energy+'\t'+'%.4e'%inten+'\t'+'%.4e'%error
        self.ui.mcaDataListWidget.addItem(string)
        self.updateMcaIntPlot()
        
    def updateMcaIntPlot(self):
        self.ui.mcaIntPlotMplWidget.canvas.ax.clear()
        self.ui.mcaIntPlotMplWidget.canvas.ax.set_ylabel('Integral Intensity')
        try:
            data=np.array(self.mcaIntData)
            y=data[:,3]
            yerr=data[:,4]
            if self.ui.mcaXAxisComboBox.currentIndex()==0:
                x=data[:,0]
                self.ui.mcaIntPlotMplWidget.canvas.ax.set_xlabel('Frame')
            elif self.ui.mcaXAxisComboBox.currentIndex()==1:
                x=data[:,1]
                self.ui.mcaIntPlotMplWidget.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
            else:
                x=data[:,2]
                self.ui.mcaIntPlotMplWidget.canvas.ax.set_xlabel('Energy (keV)')
            self.ui.mcaIntPlotMplWidget.canvas.ax.errorbar(x,y,yerr,fmt='o',label='#0')
        except:
            self.messageBox('no mca integral data to process.')
        if len(self.selectedmcafiles_rows)!=0: #plot the loaded mca int data
            for i in range(len(self.selectedmcafiles_rows)):
                data1=np.loadtxt(str(self.mcafiles[self.selectedmcafiles_rows[i]]),comments='#')
                self.ui.mcaIntPlotMplWidget.canvas.ax.errorbar(data1[:,0],data1[:,1],data1[:,2],fmt='o',label='#'+str(self.selectedmcafiles_rows[i]+1))
        self.ui.mcaIntPlotMplWidget.canvas.ax.legend(loc=2,frameon=False,scatterpoints=0,numpoints=1)
        self.ui.mcaIntPlotMplWidget.canvas.draw()
    
    def updateMcaInt(self):
        try:
            self.ui.mcaDataListWidget.clear()
            for i in range(len(self.mcaIntData)):
                if self.ui.mcaXAxisComboBox.currentIndex()==0:
                    string=str(self.mcaIntData[i][0])+'\t'+'%.4e'%self.mcaIntData[i][3]+'\t'+'%.4e'%self.mcaIntData[i][4] 
                elif self.ui.mcaXAxisComboBox.currentIndex()==1:
                    string='%.4f'%self.mcaIntData[i][1]+'\t'+'%.4e'%self.mcaIntData[i][3]+'\t'+'%.4e'%self.mcaIntData[i][4]
                else:
                    string='%.4f'%self.mcaIntData[i][2]+'\t'+'%.4e'%self.mcaIntData[i][3]+'\t'+'%.4e'%self.mcaIntData[i][4]
                self.ui.mcaDataListWidget.addItem(string)
            self.updateMcaIntPlot()
        except:
            self.messageBox('no mca integral data to process.')
                
    def mcaIntMergeDia(self):
        Dialog=QDialog(self)
        self.uimcamerge=uic.loadUi('mcamerge.ui', Dialog)
        self.uimcamerge.label.setText('Please provide the resolution of x')
        self.uimcamerge.show()
        self.connect(self.uimcamerge.cancelPushButton, SIGNAL('clicked()'), self.mcaMergeDiaClose)
        self.connect(self.uimcamerge.okPushButton, SIGNAL('clicked()'), self.mcaIntMerge)
            
    def mcaMergeDiaClose(self):
        self.uimcamerge.close()
        
    
    def mcaIntMerge(self):
        resolution=float(self.uimcamerge.resLineEdit.text())
        try:
            if self.ui.mcaXAxisComboBox.currentIndex()==0:
                index=0
            elif self.ui.mcaXAxisComboBox.currentIndex()==1:
                index=1
            else:
                index=2
            self.mcaIntData=np.array(self.mcaIntData)
            self.mcaIntData=self.mcaIntData[np.argsort(self.mcaIntData[:,index])]
           # print self.mcaIntData
            data=[self.mcaIntData[0]]
           # print data
           # print data[-1]
            for i in range(1,len(self.mcaIntData)):
                if np.abs(self.mcaIntData[i][index]-data[-1][index])<=resolution:
                    a = data[-1][3] ** 2 / data[-1][4] ** 2
                    b = self.mcaIntData[i][3] ** 2 / self.mcaIntData[i][4] ** 2
                    c = data[-1][3] / data[-1][4] ** 2
                    d = self.mcaIntData[i][3] / self.mcaIntData[i][4] ** 2
                    data[-1][3] = (a + b) / (c + d)
                    data[-1][4] = np.sqrt((a + b) / (c + d) ** 2)
                else:
                    data=np.vstack((data,self.mcaIntData[i]))
            self.mcaIntData=data
            self.uimcamerge.close()
            self.updateMcaInt()
        except:
            self.messageBox('no mca integral data to process.')
            
                
                    
    
     
    def mcaDelIntData(self):
        items=self.ui.mcaDataListWidget.selectedItems()
        self.mcaIntData=list(self.mcaIntData)
        for item in items:
            self.mcaIntData.pop(self.ui.mcaDataListWidget.row(item))
            self.ui.mcaDataListWidget.takeItem(self.ui.mcaDataListWidget.row(item))
        self.updateMcaIntPlot()
        
    def mcaSaveIntData(self):
        try:
            self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Fluorescence Integral Data',directory=self.directory))
            fid=open(self.saveFileName+'_flu.txt','w')
            try:
                fid.write('#'+str(self.command)+'\n')
            except:
                pass
            for i in range(len(self.mcaIntData)):
                if self.ui.mcaXAxisComboBox.currentIndex()==0:
                    fid.write(str(self.mcaIntData[i][0])+'\t'+str(self.mcaIntData[i][3])+'\t'+str(self.mcaIntData[i][4])+'\n')
                elif self.ui.mcaXAxisComboBox.currentIndex()==1:
                    fid.write(str(self.mcaIntData[i][1])+'\t'+str(self.mcaIntData[i][3])+'\t'+str(self.mcaIntData[i][4])+'\n')
                else:
                    fid.write(str(self.mcaIntData[i][2])+'\t'+str(self.mcaIntData[i][3])+'\t'+str(self.mcaIntData[i][4])+'\n')
            fid.close()
        except:
                print 'no mca integral data to save.'
                
    def addMcaFiles(self): #add mca files with integral intensity in the list
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Files to import', directory=self.directory)
        self.mcafiles=self.mcafiles+map(str, f)
        self.mcafnames=[]
        self.ui.mcaFileListWidget.clear()
        for i in range(len(self.mcafiles)):
            s=str(self.mcafiles[i])
            self.mcafnames.append(s[s.rfind('/')+1:])
            self.ui.mcaFileListWidget.addItem('#'+str(i+1)+'\t'+self.mcafnames[i])

    def removeMcaFiles(self):
        items=self.ui.mcaFileListWidget.selectedItems()
        for item in items:
            self.mcafnames.pop(self.ui.mcaFileListWidget.row(item))
            self.mcafiles.pop(self.ui.mcaFileListWidget.row(item))
        self.ui.mcaFileListWidget.clear()
        for i in range(len(self.reffnames)):
            self.ui.mcaFileListWidget.addItem('#'+str(i+1)+'\t'+self.mcafnames[i])    
    
    def updateSelectedMcaFiles(self):
        selectedmcafiles=self.ui.mcaFileListWidget.selectedItems()
        self.selectedmcafiles_rows=[]
        for item in selectedmcafiles:
            self.selectedmcafiles_rows.append(self.ui.mcaFileListWidget.row(item))
        self.updateMcaIntPlot()

    def messageBox(self,text,title='Warning'):
        mesgbox=QMessageBox()
        mesgbox.setText(text)
        mesgbox.setWindowTitle(title)
        mesgbox.exec_()
    
    def updateProgress(self):
        self.progressDialog.setValue(self.progressDialog.value()+1)
        
    def selectAllImages(self):
        if self.ui.imageSelectAllCheckBox.checkState()!=0:
            self.disconnect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
            for i in range(self.ui.imageListWidget.count()):
                self.ui.imageListWidget.setItemSelected(self.ui.imageListWidget.item(i),True)
            self.connect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
            self.imageSelectedScanChanged()
        else:
            self.unSelectedAllImages()
                
    def unSelectedAllImages(self):
        if self.ui.imageSelectAllCheckBox.checkState()!=0 and self.specPar[self.selectedScanNums[0]]['Detector']=='Bruker':
            self.disconnect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
            for items in self.selectedCcdFrames:
                self.ui.imageListWidget.setItemSelected(items,False)
            self.connect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
            self.ui.imageSelectAllCheckBox.setCheckState(0)
        elif self.ui.imageSelectAllCheckBox.checkState()!=0 and self.specPar[self.selectedScanNums[0]]['Detector']=='Pilatus':
            self.disconnect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
            for items in self.selectedPilFrames:
                self.ui.imageListWidget.setItemSelected(items,False)
            self.connect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
            self.ui.imageSelectAllCheckBox.setCheckState(0)
        elif self.ui.imageSelectAllCheckBox.checkState()!=0 and self.specPar[self.selectedScanNums[0]]['Detector']=='Vortex':
            self.ui.imageSelectAllCheckBox.setCheckState(0)
        self.ui.gixSumCheckBox.setCheckState(0)


    def imageSelectedScanChanged(self):
        #self.ui.imageSelectAllCheckBox.setCheckState(0)
        self.pilGIDshow=0
        self.pilGISAXSshow=0
        self.ui.statusBar.showMessage('Plotting....Wait!!')
        if  self.specPar[self.selectedScanNums[0]]['Detector']=='Vortex':
            self.mcaSelectedScanChanged()
            self.ui.statusBar.showMessage('Done')
        elif self.specPar[self.selectedScanNums[0]]['Detector']=='Bruker':
            self.det='Bruker'            
            self.ccdSelectedScanChanged()
            self.ui.statusBar.showMessage('Done')
        elif self.specPar[self.selectedScanNums[0]]['Detector']=='Pilatus':
            self.det='Pilatus'
            self.pilSelectedScanChanged()
            self.ui.statusBar.showMessage('Done')
 
        
    def updateCcdImageList(self):
        self.disconnect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
        self.ui.statusBar.clearMessage()
        self.numFrames={}
        self.ccdFileNames=[]
        self.ccdFileQzs=[]
        self.ccdMonc=np.array([])
        self.ccdX=[]
        self.ccdY=[]
        self.ccd_Dist=[]
        self.ccd_PDist=[]
        self.ccdAlpha=[]
        self.ccd_Sh=[]
        self.ccd_Wavelength=[]
        self.ccd_Tth=[]
        self.ccd_Chi=[]
        self.ccd_Gl2=[]
        self.ccd_AbsNum=[]
        self.ccd_Phi=[]
        self.ccdX_off=[]
        self.ccdY_off=[]
        for i in self.selectedScanNums:
            #self.ui.scanListWidget.setItemSelected(self.ui.scanListWidget.item(i),True)
            try:
                self.numFrames[i]=len(self.specData[i][self.specData[i]['ScanVar'][0]])
                self.ccdMonc=np.append(self.ccdMonc, self.specData[i]['Monc'])
                fname=self.ccdfhead+str(i)+'_'
                filelist=sorted(glob.glob(str(fname+'*')))
                for j in range(self.numFrames[i]):
                    self.ccdFileNames.append(filelist[j])
                    fnum=filelist[j].split('_')[-1].split('.')[0].lstrip('0')
                    if fnum=='':
                        fnum='0'
                    #self.ccdFileNames=self.ccdFileNames+[self.ccdfhead+str(i)+'_%04d'%(j,)+self.ccdftail]
                    self.ui.imageListWidget.addItem('S# '+str(i)+'\tF# '+fnum+'\tQz='+'%.4f'%self.specData[i]['L'][j]+'\tQx='+'%.4f'%self.specData[i]['H'][j]+'\tAbs=%d'%self.specPar[i]['Absorber'])
                    self.ccdFileQzs.append(self.specData[i]['L'][j])
                    self.ccdAlpha.append(np.arcsin(self.ccdFileQzs[-1]*self.specPar[i]['Wavelength']/4.0/np.pi))
                    self.ccdX.append(self.specPar[i]['DBPos'][0])
                    self.ccd_Wavelength.append(self.specPar[i]['Wavelength'])
                    self.ccdY.append(self.specPar[i]['DBPos'][1])
                    self.ccdX_off.append(self.specPar[i]['CCD_X'])
                    self.ccdY_off.append(self.specPar[i]['CCD_Y'])
                    self.ccd_Dist.append(self.specPar[i]['S2D_Dist'])
                    self.ccd_PDist.append(self.specPar[i]['S7D_Dist'])
                    self.ccd_Sh.append(self.specPar[i]['an_Sam_H'])
                    self.ccd_Tth.append(self.specPar[i]['Two_Thet'])
                    self.ccd_Chi.append(self.specPar[i]['Chi'])
                    self.ccd_Gl2.append(self.specPar[i]['g_l2'])
                    self.ccd_AbsNum.append(self.specPar[i]['Absorber'])
                    self.ccd_Phi.append(self.specPar[i]['Phi'])
            except:
                self.messageBox('Warning:: #S '+str(i+1)+' has been canceled or no data found.')
        self.ccdData={}
        self.ccdLogData={}
        self.ccdErrorData={}
        self.connect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)

    def updatePilImageList(self):
        self.disconnect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)
        self.ui.statusBar.clearMessage()
        self.numFrames={}
        self.pilFileNames=[]
        self.pilFileQzs=[]
        self.pilMonc=np.array([])
        self.pilX=[]
        self.pilY=[]
        self.pil_Dist=[]
        self.pil_PDist=[]
        self.pilAlpha=[]
        self.pilTrueAlpha=[]
        self.pil_Sh=[]
        self.pil_Wavelength=[]
        self.pil_Tth=[]
        self.pil_Chi=[]
        self.pil_Gl2=[]
        self.pil_AbsNum=[]
        self.pil_Phi=[]
        self.pil_Dth=[]
        self.pilX_off=[]
        self.pilY_off=[]
        self.pilFileQxs=[]
        self.pilFileQys=[]
        self.pilScanNum=[]
        self.pilFrameNum=[]
        self.pil_S1v=[]
        self.pil_S1h=[]
        self.pil_S4h=[]
        self.pil_S5h=[]
        for i in self.selectedScanNums:
            #self.ui.scanListWidget.setItemSelected(self.ui.scanListWidget.item(i),True)
#            try:
            self.numFrames[i]=len(self.specData[i][self.specData[i]['ScanVar'][0]])  #get frame number for each scan
            try:
                self.pilMonc=np.append(self.pilMonc, self.specData[i]['Monc'])
            except:
                self.pilMonc=np.append(self.pilMonc, self.specData[i]['monc'])
            fname=self.pilfhead+str(i)+'_'
            filelist=sorted(glob.glob(str(fname+'*')))
            for j in range(self.numFrames[i]):
                #self.pilFileNames=self.pilFileNames+[self.pilfhead+str(i)+('_%0'+str(self.specPar[i]['ImageNumber'])+'d')%(j,)+self.pilftail]
                try:
                    self.pilFileNames.append(filelist[j])
                    fnum=filelist[j].split('_')[-1].split('.')[0].lstrip('0')
                    if fnum=='':
                        fnum='0'
                    self.ui.imageListWidget.addItem('S# '+str(i)+'\tF# '+fnum+'\tQz='+'%.4f'%self.specData[i]['L'][j]+'\tQx='+'%.4f'%self.specData[i]['H'][j]+'\tAbs=%d'%self.specPar[i]['Absorber'])
                    self.pilScanNum.append(i)  
                    self.pilFrameNum.append(j)
                    self.pilFileQzs.append(self.specData[i]['L'][j])
                    self.pilFileQxs.append(self.specData[i]['H'][j])
                    self.pilFileQys.append(self.specData[i]['K'][j])
                    self.pilAlpha.append(np.arcsin(self.pilFileQzs[-1]*self.specPar[i]['Wavelength']/4.0/np.pi))
                    self.pilTrueAlpha.append(self.specPar[i]['In_Rot']*np.pi/180.0)
                    self.pilX.append(self.specPar[i]['DBPos'][1])
                    self.pil_Wavelength.append(self.specPar[i]['Wavelength'])
                    self.pilY.append(self.specPar[i]['DBPos'][0])
                    self.pil_Dist.append(self.specPar[i]['S2D_Dist'])
                    self.pil_PDist.append(self.specPar[i]['S7D_Dist'])
                    self.pil_Sh.append(self.specPar[i]['Sample_H'])
                    self.pil_S1v.append(self.specPar[i]['S1T']+self.specPar[i]['S1T'])
                    self.pil_S1h.append(self.specPar[i]['S1L']+self.specPar[i]['S1R'])
                    self.pil_S4h.append(self.specPar[i]['S4L']+self.specPar[i]['S4R'])
                    self.pil_S5h.append(self.specPar[i]['S5L']+self.specPar[i]['S5R'])
                    try:
                        self.pil_Dth.append(self.specData[i]['Det_Th'][j]*np.pi/180.0)
                    except:
                        self.pil_Dth.append(2.0*np.arcsin(np.sqrt(self.pilFileQxs[-1]**2+self.pilFileQys[-1]**2)*self.specPar[i]['Wavelength']/4.0/np.pi))
                    self.pil_Tth.append(self.specPar[i]['Two_Thet'])
                    self.pil_Chi.append(self.specPar[i]['Chi'])
                    self.pil_Gl2.append(self.specPar[i]['g_l2'])
                    self.pil_AbsNum.append(self.specPar[i]['Absorber'])
                    self.pil_Phi.append(self.specPar[i]['Phi'])
                except:
                    print 'Warning:: Frame '+str(j)+' doesnot exist!!'
#            except:
#                self.messageBox('Warning:: #S '+str(i+1)+' has been canceled or no data found.')
        self.pilData={}
        self.pilLogData={}
        self.pilErrorData={}
        self.pilGIOXSData={}
        self.pilGIOXSErrorData={}
        self.connect(self.ui.imageListWidget,SIGNAL('itemSelectionChanged()'),self.imageSelectedScanChanged)

    def ccdSelectedScanChanged(self):
        self.ui.statusBar.clearMessage()
        self.selectedCcdFrames=self.ui.imageListWidget.selectedItems()
        bad_pix=int(self.ui.gixBPCfacLineEdit.text())
        self.selectedCcdFramesNums=[self.ui.imageListWidget.row(items) for items in self.selectedCcdFrames]
        self.ccdSelectedQzs=[self.ccdFileQzs[i] for i in self.selectedCcdFramesNums]
        self.ccdSelectedMonc=[self.ccdMonc[i] for i in self.selectedCcdFramesNums]
        self.ccdSelectedX=[self.ccdX[i] for i in self.selectedCcdFramesNums]   #this is the direct beam x pixel
        self.ccdSelectedY=[self.ccdY[i] for i in self.selectedCcdFramesNums]   # this is the direct beam y pixel
        self.ccdSelected_Dist=[self.ccd_Dist[i] for i in self.selectedCcdFramesNums]   # this is the distance btw sample and CCD
        self.ccdSelected_PDist=[self.ccd_PDist[i] for i in self.selectedCcdFramesNums] # this is the distance between the two slits before the detector. Useful only for pinhole kind of setup
        self.ccdSelected_Sh=[-self.ccd_Gl2[i]*np.tan(self.ccdAlpha[i]) for i in self.selectedCcdFramesNums]#[self.ccd_Sh[i] for i in self.selectedCcdFramesNums]
        self.ccdSelectedXoff=[self.ccdX_off[i] for i in self.selectedCcdFramesNums]
        self.ccdSelectedYoff=[self.ccdY_off[i] for i in self.selectedCcdFramesNums]
        self.ccdSelected_Alpha=[self.ccdAlpha[i] for i in self.selectedCcdFramesNums]
        self.ccdSelected_Wavelength=[self.ccd_Wavelength[i] for i in self.selectedCcdFramesNums]
        self.ccdSelected_Tth=[self.ccd_Tth[i] for i in self.selectedCcdFramesNums]
        self.ccdSelected_Chi=[self.ccd_Chi[i] for i in self.selectedCcdFramesNums]
        self.ccdSelected_Gl2=[self.ccd_Gl2[i] for i in self.selectedCcdFramesNums]
        self.ccdSelected_AbsNum=[self.ccd_AbsNum[i] for i in self.selectedCcdFramesNums]
        self.ccdSelected_Phi=[self.ccd_Phi[i] for i in self.selectedCcdFramesNums]
        fac=float(self.ui.bgFacLineEdit.text())
        j=0
        self.progressDialog=QProgressDialog('Reading CCD Images','Abort',0,100)
        self.progressDialog.setWindowModality(Qt.WindowModal)
        self.progressDialog.setWindowTitle('Wait')
        self.progressDialog.setAutoClose(True)
        self.progressDialog.setAutoReset(True)
        self.progressDialog.setMinimum(1)
        self.progressDialog.setMaximum(len(self.selectedCcdFramesNums))
        self.progressDialog.show()
        
        for i in self.selectedCcdFramesNums:
            if len(self.ui.backgroundListWidget.selectedItems())==0:
                self.bruker.openFile(self.ccdFileNames[i],bad_pix=bad_pix)
                self.ccdData[i]=self.bruker.imageData
                self.ccdLogData[i]=np.log10(self.ccdData[i]+1)
                self.ccdErrorData[i]=self.bruker.errorData
                self.xyzformat='x=%.3f,y=%.3f,z=%d'
            else:
                self.bruker.openFile(self.ccdFileNames[i],bad_pix=bad_pix)
                self.ccdData[i]=self.bruker.imageData
                self.ccdErrorData[i]=self.bruker.errorData
                self.bruker.sumFiles({i:self.ccdData[i]},{i:self.ccdErrorData[i]},absfac=self.absfac,absnum=[self.ccdSelected_AbsNum[j]],mon=[self.ccdSelectedMonc[j]]) # For the Normalization Calculation Only!!
                self.ccdData[i]=self.bruker.imageData-fac*self.bgData
                if self.ui.bgIgnoreNegCheckBox.checkState()!=0:
                    self.ccdData[i]=np.where(self.ccdData[i]<0,0,self.ccdData[i])
                self.ccdLogData[i]=np.log10(np.where(self.ccdData[i]<=0,1e-10,self.ccdData[i]))
                self.ccdErrorData[i]=np.sqrt(self.bruker.errorData**2+fac**2*self.bgError**2)
                j=j+1
                self.xyzformat='x=%.3f,y=%.3f,z=%.2e'
            self.progressDialog.setLabelText('Reading CCD Frame #'+str(i))     
            self.updateProgress()
            if self.progressDialog.wasCanceled()==True:
                break
        self.progressDialog.hide()
        self.ccdFrameNumsWithSameQz={}
        self.ccdFrameAbsWithSameQz={}
        for i in self.selectedCcdFramesNums:
            if float(format(self.ccdFileQzs[i],'.4f')) in self.ccdFrameNumsWithSameQz:                
                self.ccdFrameNumsWithSameQz[float(format(self.ccdFileQzs[i],'.4f'))].append(i)
                self.ccdFrameAbsWithSameQz[float(format(self.ccdFileQzs[i],'.4f'))].append(self.ccd_AbsNum[i])
            else:
                self.ccdFrameNumsWithSameQz[float(format(self.ccdFileQzs[i],'.4f'))]=[i]
                self.ccdFrameAbsWithSameQz[float(format(self.ccdFileQzs[i],'.4f'))]=[self.ccd_AbsNum[i]]
        try:
            self.imageMax=np.max(self.ccdData[self.selectedCcdFramesNums[0]])
            #self.imageMin=np.max([1e-10,np.min(self.ccdData[self.selectedCcdFramesNums[0]])])
            self.imageMin=np.min(self.ccdData[self.selectedCcdFramesNums[0]][np.nonzero(self.ccdData[self.selectedCcdFramesNums[0]])])
        except:
            self.imageMax=100
            self.imageMin=0
        self.ui.gixMaxLineEdit.setText(str(self.imageMax))
        self.ui.gixMinLineEdit.setText(str(self.imageMin))
        self.ui.gixMinHorizontalSlider.setRange(0, 100)
        self.ui.gixMaxHorizontalSlider.setRange(0, 100)
        self.ui.gixMinHorizontalSlider.setValue(0)
        self.ui.gixMaxHorizontalSlider.setValue(100)
#        try:
        self.ui.PlotWidget.setCurrentIndex(1)
        self.updateCcdPlotData()
#        except:
#            self.messageBox('Warning:: Some scans maybe missing')
            
    def pilSelectedScanChanged(self):
        self.pilmovieindex=0
        #        self.disconnect=(self.ui.gidComboBox, SIGNAL('currentIndexChanged(int)'), self.pilGID)
#        self.ui.gidComboBox.setCurrentIndex(0)
#        self.connect=(self.ui.gidComboBox, SIGNAL('currentIndexChanged(int)'), self.pilGID)
        bad_pix=int(self.ui.pilBPCfacLineEdit.text())
        self.ui.statusBar.clearMessage()
        self.selectedPilFrames=self.ui.imageListWidget.selectedItems()
        self.selectedPilFramesNums=np.sort([self.ui.imageListWidget.row(items) for items in self.selectedPilFrames])
        self.pilSelectedScanNum=[self.pilScanNum[i] for i in self.selectedPilFramesNums]
        self.pilSelectedQzs=[self.pilFileQzs[i] for i in self.selectedPilFramesNums]
        self.pilSelectedMonc=[self.pilMonc[i] for i in self.selectedPilFramesNums]
        self.pilSelectedX=[self.pilX[i] for i in self.selectedPilFramesNums]   #this is the direct beam x pixel
        self.pilSelectedY=[self.pilY[i] for i in self.selectedPilFramesNums]   # this is the direct beam y pixel
        self.pilSelected_Dist=[self.pil_Dist[i] for i in self.selectedPilFramesNums]   # this is the distance btw sample and Pilatus
        self.pilSelected_PDist=[self.pil_PDist[i] for i in self.selectedPilFramesNums] # this is the distance between the two slits before the detector. Useful only for pinhole kind of setup
        self.pilSelected_Sh=[-self.pil_Gl2[i]*np.tan(self.pilAlpha[i]) for i in self.selectedPilFramesNums]#[self.ccd_Sh[i] for i in self.selectedCcdFramesNums]
#        self.ccdSelectedXoff=[self.ccdX_off[i] for i in self.selectedCcdFramesNums]
#        self.ccdSelectedYoff=[self.ccdY_off[i] for i in self.selectedCcdFramesNums]
        self.pilSelected_Alpha=[self.pilAlpha[i] for i in self.selectedPilFramesNums]
        self.pilSelected_TrueAlpha=[self.pilTrueAlpha[i] for i in self.selectedPilFramesNums]
        self.pilSelected_Wavelength=[self.pil_Wavelength[i] for i in self.selectedPilFramesNums]
        self.pilSelected_Dth=[self.pil_Dth[i] for i in self.selectedPilFramesNums]
        self.pilSelected_Tth=[self.pil_Tth[i] for i in self.selectedPilFramesNums]
        self.pilSelected_Chi=[self.pil_Chi[i] for i in self.selectedPilFramesNums]
        self.pilSelected_Gl2=[self.pil_Gl2[i] for i in self.selectedPilFramesNums]
        self.pilSelected_AbsNum=[self.pil_AbsNum[i] for i in self.selectedPilFramesNums]
        self.pilSelected_Phi=[self.pil_Phi[i] for i in self.selectedPilFramesNums]
        fac=float(self.ui.bgFacLineEdit.text())
        j=0
        self.progressDialog=QProgressDialog('Reading Pilatus Images','Abort',0,100)
        self.progressDialog.setWindowModality(Qt.WindowModal)
        self.progressDialog.setWindowTitle('Wait')
        self.progressDialog.setAutoClose(True)
        self.progressDialog.setAutoReset(True)
        self.progressDialog.setMinimum(1)
        self.progressDialog.setMaximum(len(self.selectedPilFramesNums))
        self.progressDialog.show()
        for i in self.selectedPilFramesNums:
            if len(self.ui.backgroundListWidget.selectedItems())==0:
                self.pilatus.openFile(self.pilFileNames[i],bad_pix=bad_pix)
                self.pilData[i]=self.pilatus.imageData
                self.pilLogData[i]=np.log10(self.pilData[i]+1)
                self.pilErrorData[i]=self.pilatus.errorData
                self.xyzformat='x=%.3f,y=%.3f,z=%.2e'
            else:
                self.pilatus.openFile(self.pilFileNames[i],bad_pix=bad_pix)
                self.pilData[i]=self.pilatus.imageData
                self.pilErrorData[i]=self.pilatus.errorData
                self.pilatus.sumFiles({i:self.pilData[i]},{i:self.pilErrorData[i]},absfac=self.absfac,absnum=[self.pilSelected_AbsNum[j]],mon=[self.pilSelectedMonc[j]]) # For the Normalization Calculation Only!!
                self.pilGIOXSData[i]=self.pilatus.imageData   #for GIOXS data, bg substraction has to be done in 1D format. So keep it without bg image substraction
                self.pilGIOXSErrorData[i]=self.pilatus.errorData
                self.pilData[i]=self.pilatus.imageData-fac*self.bgData
                if self.ui.bgIgnoreNegCheckBox.checkState()!=0:
                    self.pilData[i]=np.where(self.pilData[i]<0,0,self.pilData[i])
                self.pilLogData[i]=np.log10(np.where(self.pilData[i]<=0,1e-10,self.pilData[i]))
                self.pilErrorData[i]=np.sqrt(self.pilatus.errorData**2+fac**2*self.bgError**2)
                j=j+1
                self.xyzformat='x=%.3f,y=%.3f,z=%.2e'
            self.progressDialog.setLabelText('Reading Pilatus Frame #'+str(i))     
            self.updateProgress()
            if self.progressDialog.wasCanceled()==True:
                break
        self.progressDialog.hide()
        self.pilFrameNumsWithSameQz={}
        self.pilFrameAbsWithSameQz={} 
        for i in self.selectedPilFramesNums:
            if float(format(self.pilFileQzs[i],'.4f')) in self.pilFrameNumsWithSameQz:                
                self.pilFrameNumsWithSameQz[float(format(self.pilFileQzs[i],'.4f'))].append(i)
                self.pilFrameAbsWithSameQz[float(format(self.pilFileQzs[i],'.4f'))].append(self.pil_AbsNum[i])
            else:
                self.pilFrameNumsWithSameQz[float(format(self.pilFileQzs[i],'.4f'))]=[i]
                self.pilFrameAbsWithSameQz[float(format(self.pilFileQzs[i],'.4f'))]=[self.pil_AbsNum[i]]
        try:
            self.imageMax=np.max(self.pilData[self.selectedPilFramesNums[0]])
            if self.imageMax<1:
                self.imageMax=1
            #self.imageMin=np.max([1e-10,np.min(self.pilData[self.selectedPilFramesNums[0]])])
            self.imageMin=np.min(self.pilData[self.selectedPilFramesNums[0]][np.nonzero(self.pilData[self.selectedPilFramesNums[0]])])
        except:
            self.imageMax=100
            self.imageMin=0
        self.ui.pilMaxLineEdit.setText(str(self.imageMax))
        self.ui.pilMinLineEdit.setText(str(self.imageMin))
        self.ui.pilMinHorizontalSlider.setRange(0, 100)
        self.ui.pilMaxHorizontalSlider.setRange(0, 100)
        self.ui.pilMinHorizontalSlider.setValue(0)
        self.ui.pilMaxHorizontalSlider.setValue(100)
#        try:
        self.ui.PlotWidget.setCurrentIndex(2)
        self.updatePilPlotData()
#        except:
#            self.messageBox('Warning:: Some scans maybe missing')
            
    def updateMaxSlider(self):
        self.ui.gixMaxLineEdit.setText(str(self.ui.gixMaxHorizontalSlider.value()*(self.imageMax-self.imageMin)/100.0))
        self.ui.gixMinHorizontalSlider.setMaximum(int(float(self.ui.gixMaxLineEdit.text())*100.0/(self.imageMax-self.imageMin)))
        self.updateCcdPlotData()
        
    def updatePilMaxSlider(self):
        self.ui.pilMaxLineEdit.setText(str(self.ui.pilMaxHorizontalSlider.value()*(self.imageMax-self.imageMin)/100.0))
        self.ui.pilMinHorizontalSlider.setMaximum(int(float(self.ui.pilMaxLineEdit.text())*100.0/(self.imageMax-self.imageMin)))
        self.updatePilPlotData()
        
    def updateMinSlider(self):
        self.ui.gixMinLineEdit.setText(str(self.ui.gixMinHorizontalSlider.value()*(self.imageMax-self.imageMin)/100.0))
        self.updateCcdPlotData()
        
    def updatePilMinSlider(self):
        self.ui.pilMinLineEdit.setText(str(self.ui.pilMinHorizontalSlider.value()*(self.imageMax-self.imageMin)/100.0))
        self.updatePilPlotData()
        
    def addBGImagestoList(self):        
        self.ui.backgroundListWidget.clear()
        for item in self.ui.imageListWidget.selectedItems():
            self.ui.backgroundListWidget.addItem(item.text())
        self.bgccdData={}
        self.bgccdErrorData={}
        self.bgpilData={}
        self.bgpilErrorData={}
        self.bgMonc={}
        self.bgAbsNum={}
        self.bgScanNum={}
        self.bgWavelength={}
        self.bgDist={}
        self.bgAlpha={}
        self.bgTrueAlpha={}
        j=0
        if self.det=='Bruker':
            for i in self.selectedCcdFramesNums:
                self.bgccdData[j]=self.ccdData[i]
                self.bgccdErrorData[j]=self.ccdErrorData[i]
                self.bgMonc[j]=self.ccdSelectedMonc[j]
                self.bgAbsNum[j]=self.ccdSelected_AbsNum[j]
                j=j+1
        elif self.det=='Pilatus':
            for i in self.selectedPilFramesNums:
                self.bgpilData[j]=self.pilData[i]
                self.bgpilErrorData[j]=self.pilErrorData[i]
                self.bgMonc[j]=self.pilSelectedMonc[j]
                self.bgWavelength[j]=self.pilSelected_Wavelength[j]
                self.bgDist[j]=self.pilSelected_Dist[j]
                self.bgAlpha[j]=self.pilSelected_Alpha[j]
                self.bgTrueAlpha[j]=self.pilSelected_TrueAlpha[j]
                self.bgAbsNum[j]=self.pilSelected_AbsNum[j]
                self.bgScanNum[j]=self.pilSelectedScanNum[j]
                j=j+1
            
    def bgSelectAll(self):
        self.disconnect(self.ui.backgroundListWidget, SIGNAL('itemSelectionChanged()'),self.bgSelectionChanged)
        if self.ui.bgSelectAllCheckBox.checkState()!=0:
            for i in range(self.ui.backgroundListWidget.count()):
                self.ui.backgroundListWidget.setItemSelected(self.ui.backgroundListWidget.item(i),True)
        else:
            for i in range(self.ui.backgroundListWidget.count()):
                self.ui.backgroundListWidget.setItemSelected(self.ui.backgroundListWidget.item(i),False)
        self.connect(self.ui.backgroundListWidget, SIGNAL('itemSelectionChanged()'),self.bgSelectionChanged)
        self.bgSelectionChanged()
        
            
    def bgSelectionChanged(self):
        self.selectedBgFrameNums=[self.ui.backgroundListWidget.row(item) for item in self.ui.backgroundListWidget.selectedItems()]
        #print self.selectedBgFrameNums
        self.selectedBgData={}
        self.selectedBgError={}
        self.selectedBgAbsNum=[]
        self.selectedBgScanNum=[]
        self.selectedBgMonc=[]
        self.selectedBgWavelength=[]
        self.selectedBgDist=[]
        self.selectedBgAlpha=[]
        self.selectedBgTrueAlpha=[]
        
        if self.det=='Bruker':
            for i in self.selectedBgFrameNums:
                self.selectedBgData[i]=self.bgccdData[i]
                self.selectedBgError[i]=self.bgccdErrorData[i]
                self.selectedBgMonc.append(self.bgMonc[i])
                self.selectedBgAbsNum.append(self.bgAbsNum[i])
                self.bruker.sumFiles(self.selectedBgData,self.selectedBgError,absfac=self.absfac,absnum=self.selectedBgAbsNum,mon=self.selectedBgMonc)
                self.bgData=self.bruker.imageData
                self.bgError=self.bruker.errorData
        elif self.det=='Pilatus':
            for i in self.selectedBgFrameNums:
                self.selectedBgData[i]=self.bgpilData[i]
                self.selectedBgError[i]=self.bgpilErrorData[i]
                self.selectedBgMonc.append(self.bgMonc[i])
                self.selectedBgWavelength.append(self.bgWavelength[i])
                self.selectedBgDist.append(self.bgDist[i])
                self.selectedBgAlpha.append(self.bgAlpha[i])
                self.selectedBgTrueAlpha.append(self.bgTrueAlpha[i])
                self.selectedBgAbsNum.append(self.bgAbsNum[i])
                self.selectedBgScanNum.append(self.bgScanNum[i])
            self.pilatus.sumFiles(self.selectedBgData,self.selectedBgError,absfac=self.absfac,absnum=self.selectedBgAbsNum,mon=self.selectedBgMonc)
            self.bgData=self.pilatus.imageData
            self.bgError=self.pilatus.errorData
        self.pilSelectedScanChanged()
#        if len(self.ui.backgroundListWidget.selectedItems())>0:
#            self.ui.gixSumCheckBox.setCheckState(2)
        
        
    def removeBGImages(self):
        self.disconnect(self.ui.backgroundListWidget, SIGNAL('itemSelectionChanged()'),self.bgSelectionChanged)
        for item in self.ui.backgroundListWidget.selectedItems():
            self.ui.backgroundListWidget.takeItem(self.ui.backgroundListWidget.row(item))  
        self.connect(self.ui.backgroundListWidget, SIGNAL('itemSelectionChanged()'),self.bgSelectionChanged)
        self.ui.bgSelectAllCheckBox.setCheckState(0)
        
    def update2dPlots(self):
        if self.det=='Bruker':
            self.updateCcdPlotData()
        elif self.det=='Pilatus':
            self.updatePilPlotData()

            
    def updateCcdPlotData(self):
        if len(self.selectedCcdFramesNums)!=0:
            if len(self.selectedCcdFramesNums)>10 and self.ui.gixSumCheckBox.checkState()<1:
                self.messageBox("The number of selected frames are more than 10. So plotting only first 10 images")
                N=10
            else:
                N=len(self.selectedCcdFramesNums)   
            self.ui.gixMplWidget.canvas.fig.clf()
            cmap=str(self.ui.gixCMapComboBox.currentText())
            self.ui.gixMaxHorizontalSlider.setValue(int(float(self.ui.gixMaxLineEdit.text())*100/(self.imageMax-self.imageMin)))
            self.ui.gixMinHorizontalSlider.setValue(int(float(self.ui.gixMinLineEdit.text())*100/(self.imageMax-self.imageMin)))
            vmax=float(self.ui.gixMaxLineEdit.text())
            vmin=float(self.ui.gixMinLineEdit.text())
            ccdXMin=np.ones_like(self.selectedCcdFramesNums)
            ccdYMin=np.ones_like(self.selectedCcdFramesNums)
            ccdXMax=float(self.bruker.NCOLS)*np.ones_like(self.selectedCcdFramesNums)
            ccdYMax=float(self.bruker.NROWS)*np.ones_like(self.selectedCcdFramesNums)
            xlabel='PixX'
            ylabel='PixY'
            self.wavelength=self.ccdSelected_Wavelength
            self.alpha=self.ccdSelected_Alpha
            if self.ui.gixSpecCheckBox.checkState()!=0:
                if self.checkSameArrays(self.ccdSelectedX)==False or self.checkSameArrays(self.ccdSelectedY)==False or self.checkSameArrays(self.ccdSelected_Dist)==False:
                    self.messageBox('Warning:: The seleted frames have different CCD parameters; the parameters of the first frames are used here!!')
                self.ui.gixSDDistLineEdit.setText(str(self.ccdSelected_Dist[0]))
                self.ui.gixDBPosLineEdit.setText(str(self.ccdSelectedX[0])+','+str(self.ccdSelectedY[0]))
                self.ui.gixCcd_OffLineEdit.setText('%0.2f'%self.ccdSelectedXoff[0]+','+'%0.2f'%self.ccdSelectedYoff[0])
                self.xoff=np.array(map(int, np.array(self.ccdSelectedXoff)/0.06))
                self.yoff=np.array(map(int, np.array(self.ccdSelectedYoff)/0.06))
                self.distance=np.array(self.ccdSelected_Dist)
                self.xcenter=np.array(self.ccdSelectedX)
                self.ycenter=np.array(self.ccdSelectedY)
            else:
                self.xoff=np.array(map(int,float(self.ui.gixCcd_OffLineEdit.text().split(',')[0])*np.ones_like(self.selectedCcdFramesNums)/0.06))
                self.yoff=np.array(map(int,float(self.ui.gixCcd_OffLineEdit.text().split(',')[-1])*np.ones_like(self.selectedCcdFramesNums)/0.06))
                self.distance=float(self.ui.gixSDDistLineEdit.text())*np.ones_like(self.selectedCcdFramesNums)
                self.xcenter=float(self.ui.gixDBPosLineEdit.text().split(',')[0])*np.ones_like(self.selectedCcdFramesNums)
                self.ycenter=float(self.ui.gixDBPosLineEdit.text().split(',')[-1])*np.ones_like(self.selectedCcdFramesNums)
                
            if str(self.ui.gixAxesComboBox.currentText())=='Angles':
                xlabel='Psi [Degrees]'
                ylabel='Theta [Degrees]'
                ccdXMin=(ccdXMin-(self.xcenter-self.xoff))*0.06*180.0/np.pi/self.distance
                ccdXMax=(ccdXMax-(self.xcenter-self.xoff))*0.06*180.0/np.pi/self.distance
                ccdYMin=(((self.ycenter-self.yoff)-ccdYMin)*0.06-self.ccdSelected_Sh)*180.0/np.pi/self.distance
                ccdYMax=(((self.ycenter-self.yoff)-ccdYMax)*0.06-self.ccdSelected_Sh)*180.0/np.pi/self.distance
            elif str(self.ui.gixAxesComboBox.currentText())=='Q':
                xlabel='Qxy [1/Angs]'
                ylabel='Qz [1/Angs]'
                ccdXMin=4.0*np.pi*np.sin((ccdXMin-(self.xcenter-self.xoff))*0.06/self.distance/2.0)/self.wavelength
                ccdXMax=4.0*np.pi*np.sin((ccdXMax-(self.xcenter-self.xoff))*0.06/self.distance/2.0)/self.wavelength
                ccdYMin=2.0*np.pi*(np.sin((((self.ycenter-self.yoff)-ccdYMin)*0.06-self.ccdSelected_Sh)/self.distance)+np.sin(self.alpha))/self.wavelength
                ccdYMax=2.0*np.pi*(np.sin((((self.ycenter-self.yoff)-ccdYMax)*0.06-self.ccdSelected_Sh)/self.distance)+np.sin(self.alpha))/self.wavelength
                
            if self.ui.gixSumCheckBox.checkState()!=0:
                if self.checkSameQzs(self.ccdSelectedQzs)==True:
                    self.sumData=self.ccdData[self.selectedCcdFramesNums[0]]
                    for i in self.selectedCcdFramesNums[1:]:
                        self.sumData=self.sumData+self.ccdData[i]
                    self.sumData=self.sumData/float(len(self.selectedCcdFramesNums))
                    self.logData=np.log10(np.where(self.sumData<=0,1e-10,self.sumData))
                    ax=self.ui.gixMplWidget.canvas.fig.add_subplot(1,1,1)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)
                    ax.set_title('Sum'+str(self.selectedCcdFramesNums))
                    self.extent=[ccdXMin[0],ccdXMax[0],ccdYMax[0],ccdYMin[0]]
                    if self.ui.gixLogIntCheckBox.checkState()!=0:
                        p=ax.imshow(self.logData,interpolation='nearest',extent=self.extent,vmax=np.log10(vmax),vmin=np.log10(vmin),cmap=cmap,aspect='equal') 
                        self.Zdata=self.logData
                    else:
                        p=ax.imshow(self.sumData,interpolation='nearest',extent=self.extent,vmax=vmax,vmin=vmin,cmap=cmap,aspect='equal')
                        self.Zdata=self.sumData 
                    ax.format_coord=self.format_coord
                    self.ui.gixMplWidget.canvas.fig.colorbar(p)
                else:
                    self.messageBox('Error:: The selected frames have different qz values!!')
            else:
                row=2
                while N/row+pl.where(pl.mod(N,row)>0,1,0)>=row:
                    row=row+1
                row=row-1
                col=N/row+pl.where(pl.mod(N,row)>0,1,0)
                num=1
                ax={}
                for i in self.selectedCcdFramesNums[:N]:
                    ax[i]=self.ui.gixMplWidget.canvas.fig.add_subplot(row,col,num)
                    ax[i].set_xlabel(xlabel)
                    ax[i].set_ylabel(ylabel)
                    ax[i].set_title('F #'+str(i))
                    self.extent=[ccdXMin[num-1],ccdXMax[num-1],ccdYMax[num-1],ccdYMin[num-1]]
                    if self.ui.gixLogIntCheckBox.checkState()!=0:
                        p=ax[i].imshow(self.ccdLogData[i],interpolation='nearest',extent=self.extent,vmax=np.log10(vmax),vmin=np.log10(vmin),cmap=cmap,aspect='equal')
                        self.Zdata=self.ccdLogData[i] 
                    else:
                        p=ax[i].imshow(self.ccdData[i],interpolation='nearest',extent=self.extent,vmax=vmax,vmin=vmin,cmap=cmap,aspect='equal')
                        self.Zdata=self.ccdData[i]
                    ax[i].format_coord=self.format_coord
                    self.ui.gixMplWidget.canvas.fig.colorbar(p)
                    num=num+1
            self.ui.gixMplWidget.canvas.draw()
        
        
    def updatePilPlotData(self):
        aspect=str(self.ui.pilAspectComboBox.currentText())
        if self.pilGIDshow==0 and self.pilGISAXSshow==0 and len(self.selectedPilFramesNums)!=0:
            if len(self.selectedPilFramesNums)>10 and self.ui.gixSumCheckBox.checkState()<1:
                #self.messageBox("The number of selected frames are more than 10. So plotting only first 10 images")
                N=10
            else:
                N=len(self.selectedPilFramesNums)    
            self.ui.pilMplWidget.canvas.fig.clf()
            cmap=str(self.ui.pilCMapComboBox.currentText())
            self.ui.pilMaxHorizontalSlider.setValue(int(float(self.ui.pilMaxLineEdit.text())*100/(self.imageMax-self.imageMin)))
            self.ui.pilMinHorizontalSlider.setValue(int(float(self.ui.pilMinLineEdit.text())*100/(self.imageMax-self.imageMin)))
            vmax=float(self.ui.pilMaxLineEdit.text())
            vmin=float(self.ui.pilMinLineEdit.text())
            pilXMin=np.ones_like(self.selectedPilFramesNums)-1
            pilYMin=np.ones_like(self.selectedPilFramesNums)-1
            pilXMax=float(self.pilatus.NCOLS)*np.ones_like(self.selectedPilFramesNums)-1   #[194,...,194]
            pilYMax=float(self.pilatus.NROWS)*np.ones_like(self.selectedPilFramesNums)-1
            xlabel='PixX'
            ylabel='PixY'
            self.wavelength=np.array(self.pilSelected_Wavelength)
            self.alpha=np.array(self.pilSelected_Alpha)
            self.truealpha=np.array(self.pilSelected_TrueAlpha)
            self.dth=np.array(self.pilSelected_Dth)
            if self.ui.pilSpecCheckBox.checkState()!=0: #if this CB is checked, use spec value
                if self.checkSameArrays(self.pilSelectedX)==False or self.checkSameArrays(self.pilSelectedY)==False or self.checkSameArrays(self.pilSelected_Dist)==False:
                    self.messageBox('Warning:: The selected frames have different Pilatus parameters; the parameters of the first frames are used here!!')
                try:
                    self.ui.pilSDDistLineEdit.setText(str(self.pilSelected_Dist[0]))
                    self.ui.pilDBPosLineEdit.setText(str(int(pilXMax[0])-self.pilSelectedX[0])+','+str(self.pilSelectedY[0]))
                   # if self.ui.pilPinholeCheckBox.checkState()!=0:
                   #     self.distance=np.array(self.pilSelected_PDist)
                   # else:
                   #     self.distance=np.array(self.pilSelected_Dist)
                    self.distance=np.array(self.pilSelected_Dist)
                    self.xcenter=np.array(pilXMax[0]-self.pilSelectedX) # image was flipped from left to right, i.e., x=194
                    self.ycenter=np.array(self.pilSelectedY)
                except:
                    pass
            else:   #if this CB is unchecked, use input value   
                self.distance=float(self.ui.pilSDDistLineEdit.text())*np.ones_like(self.selectedPilFramesNums)
                self.xcenter=float(self.ui.pilDBPosLineEdit.text().split(',')[0])*np.ones_like(self.selectedPilFramesNums)
                self.ycenter=float(self.ui.pilDBPosLineEdit.text().split(',')[-1])*np.ones_like(self.selectedPilFramesNums)
                #self.ui.pilSDDistLineEdit.setText(str(self.pilSelected_Dist[0]))
                #self.ui.pilDBPosLineEdit.setText(str(self.pilSelectedX[0])+','+str(self.pilSelectedY[0]))
                #if self.ui.pilPinholeCheckBox.checkState()!=0:
                #    self.distance=np.array(self.pilSelected_PDist)
                #else:
                #    self.distance=self.distance
                #self.xcenter=np.array(self.pilSelectedX)
                #self.ycenter=np.array(self.pilSelectedY).ones_like(self.selectedPilFramesNums)
            if str(self.ui.pilAxesComboBox.currentText())=='Angles':
                xlabel=r'$2\theta$'+' '+r'$[\deg]$'
                ylabel=r'$\beta$'+' '+r'$[\deg]$'
                pilYMin=((pilYMin-self.ycenter)*0.172/self.distance+np.arcsin(2*np.sin(self.alpha)-np.sin(self.truealpha)))*180/np.pi
                pilYMax=((pilYMax-self.ycenter)*0.172/self.distance+np.arcsin(2*np.sin(self.alpha)-np.sin(self.truealpha)))*180/np.pi
                if self.ui.pilPinholeCheckBox.checkState()!=0:
                    self.distance=np.array(self.pilSelected_PDist)     #distance between S4 & S5 for pinhole setup
                #else:
                #    self.distance=np.array(self.pilSelected_Dist)
                pilXMin=(self.dth+(pilXMin-self.xcenter)*0.172/self.distance)*180.0/np.pi
                pilXMax=(self.dth+(pilXMax-self.xcenter)*0.172/self.distance)*180.0/np.pi
                pilXCen=self.dth*180.0/np.pi
                
            elif str(self.ui.pilAxesComboBox.currentText())=='Q':
                xlabel=r'$Q_{xy}$'+' '+r'$[\AA^{-1}]$'
               # ylabel='Qz [1/Angs]'
                ylabel=r'$Q_z$'+' '+r'$[\AA^{-1}]$'
                pilYMin=2.0*np.pi*(np.sin((pilYMin-self.ycenter)*0.172/self.distance+np.arcsin(2*np.sin(self.alpha)-np.sin(self.truealpha)))+np.sin(self.truealpha))/self.wavelength
                pilYMax=2.0*np.pi*(np.sin((pilYMax-self.ycenter)*0.172/self.distance+np.arcsin(2*np.sin(self.alpha)-np.sin(self.truealpha)))+np.sin(self.truealpha))/self.wavelength
                if self.ui.pilPinholeCheckBox.checkState()!=0:
                    self.distance=np.array(self.pilSelected_PDist)     #distance between S4 & S5 for pinhole setup
                #else:
                #    self.distance=np.array(self.pilSelected_Dist)
                pilXMin=4.0*np.pi*np.sin((self.dth+(pilXMin-self.xcenter)*0.172/self.distance)/2.0)/self.wavelength
                pilXMax=4.0*np.pi*np.sin((self.dth+(pilXMax-self.xcenter)*0.172/self.distance)/2.0)/self.wavelength
                pilXCen=4.0*np.pi*np.sin(self.dth/2.0)/self.wavelength
                
            if self.pilmovieindex==1:
                for i in self.selectedPilFramesNums:
                    self.ui.pilMplWidget.canvas.fig.clf()
                    ax=self.ui.pilMplWidget.canvas.fig.add_subplot(1,1,1)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel) 
                    ax.set_title('F #'+str(i),fontsize=16)
                    self.extent=[pilXMin[i],pilXMax[i],pilYMin[i],pilYMax[i]]
                    if self.ui.pilLogIntCheckBox.checkState()!=0:
                        p=ax.imshow(self.pilLogData[i],interpolation='bicubic',extent=self.extent,vmax=np.log10(vmax),vmin=np.log10(vmin),cmap=cmap,aspect=aspect,origin='lower')
                        self.Zdata=self.pilLogData[i] 
                    else:
                        p=ax.imshow(self.pilData[i],interpolation='bicubic',extent=self.extent,vmax=vmax,vmin=vmin,cmap=cmap,aspect=aspect,origin='lower')
                        self.Zdata=self.pilData[i]
                    ax.format_coord=self.format_pil_coord
                    ax.autoscale(False)
                    ax.plot([pilXCen[i],pilXCen[i]], [pilYMin[i],pilYMax[i]], 'w--', linewidth=2.0)
                    self.ui.pilMplWidget.canvas.fig.colorbar(p)
                    self.ui.pilMplWidget.canvas.fig.suptitle('File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1])
                    self.ui.pilMplWidget.canvas.fig.tight_layout()
                    self.ui.pilMplWidget.canvas.fig.subplots_adjust(top=0.9)
                    self.ui.pilMplWidget.canvas.draw()
                    self.ui.pilMplWidget.canvas.flush_events()
                    time.sleep(0.5)
                    
            else:
                if self.ui.gixSumCheckBox.checkState()!=0:   #if "Sum All" is checked; merge the frames having the same qz
                    if self.checkSameQzs(self.pilSelectedQzs)==True: # check if qz are same in selected frames
                        self.sumData=self.pilData[self.selectedPilFramesNums[0]]   #self.pilData only has pixel and intensity;  
                        for i in self.selectedPilFramesNums[1:]:
                            self.sumData=self.sumData+self.pilData[i]
                        self.sumData=self.sumData/float(len(self.selectedPilFramesNums)) # average the intensity
                        self.logData=np.log10(np.where(self.sumData<=0,1e-10,self.sumData))
                        ax=self.ui.pilMplWidget.canvas.fig.add_subplot(1,1,1)
                        ax.set_xlabel(xlabel)
                        ax.set_ylabel(ylabel)
                        ax.set_title('Sum'+str(self.selectedPilFramesNums))
                        self.extent=[pilXMin[0],pilXMax[0],pilYMin[0],pilYMax[0]]
                        if self.ui.pilLogIntCheckBox.checkState()!=0:
                            p=ax.imshow(self.logData,interpolation='bicubic',extent=self.extent,vmax=np.log10(vmax),vmin=np.log10(vmin),cmap=cmap,aspect=aspect,origin='lower') 
                            self.Zdata=self.logData
                        else:
                            p=ax.imshow(self.sumData,interpolation='bicubic',extent=self.extent,vmax=vmax,vmin=vmin,cmap=cmap,aspect=aspect,origin='lower')
                            self.Zdata=self.sumData 
                        ax.format_coord=self.format_pil_coord
                        self.ui.pilMplWidget.canvas.fig.colorbar(p)
                    else:
                        self.messageBox('Error:: The selected frames have different qz values!!')
                else:
                    row=2
                    while N/row+pl.where(pl.mod(N,row)>0,1,0)>=row:  #setup subplot format for n of selected frames. 
                        row=row+1
                    row=row-1
                    col=N/row+pl.where(pl.mod(N,row)>0,1,0)
                    num=1
                    ax={}
                    for i in self.selectedPilFramesNums[:N]:
                        ax[i]=self.ui.pilMplWidget.canvas.fig.add_subplot(row,col,num)
                        ax[i].set_xlabel(xlabel,fontsize=16)
                        ax[i].set_ylabel(ylabel,fontsize=16)
                        ax[i].set_title('F #'+str(i),fontsize=16)
                        self.extent=[pilXMin[num-1],pilXMax[num-1],pilYMin[num-1],pilYMax[num-1]]
                        if self.ui.pilLogIntCheckBox.checkState()!=0:
                            p=ax[i].imshow(self.pilLogData[i],interpolation='bicubic',extent=self.extent,vmax=np.log10(vmax),vmin=np.log10(vmin),cmap=cmap,aspect=aspect,origin='lower')
                            self.Zdata=self.pilLogData[i] 
                        else:
                            p=ax[i].imshow(self.pilData[i],interpolation='bicubic',extent=self.extent,vmax=vmax,vmin=vmin,cmap=cmap,aspect=aspect,origin='lower')
                            self.Zdata=self.pilData[i]
                        ax[i].format_coord=self.format_pil_coord
                        self.ui.pilMplWidget.canvas.fig.colorbar(p,format="%.1e")
                        num=num+1
                self.ui.pilMplWidget.canvas.fig.suptitle('File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1])
                self.ui.pilMplWidget.canvas.fig.tight_layout()
                self.ui.pilMplWidget.canvas.fig.subplots_adjust(top=0.9)
              # self.ui.pilMplWidget.canvas.fig.subplots_adjust(left=0)
              #  self.ui.pilMplWidget.canvas.fig.subplots_adjust(right=1)
                self.ui.pilMplWidget.canvas.draw()
        elif self.pilGIDshow!=0 and len(self.selectedPilFramesNums)!=0:            
            self.pilGIDPlot()
        elif self.pilGISAXSshow!=0 and len(self.selectedPilFramesNums)!=0:
            self.showPilGISAXS()
            
    def pilMovieShow(self):
        self.disconnect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
        self.ui.pilAxesComboBox.setCurrentIndex(2)
        self.connect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
        self.pilmovieindex=1
        self.disconnect(self.ui.pilPinholeCheckBox,SIGNAL('stateChanged(int)'),self.update2dPlots)
        self.ui.pilPinholeCheckBox.setCheckState(2)
        self.updatePilPlotData()
        self.pilmovieindex=0
        self.ui.pilPinholeCheckBox.setCheckState(0)
        self.connect(self.ui.pilPinholeCheckBox,SIGNAL('stateChanged(int)'),self.update2dPlots)
        
    def gioxsDisplay(self):
        if len(self.selectedBgFrameNums)==0 or len(self.selectedPilFramesNums)==0:
            self.messageBox('Warning: No GIOXS data or background selected!!')
        else:
            Dialog=QDialog(self)
            self.uipdgioxs=uic.loadUi('pdgioxs.ui', Dialog)
            self.uipdgioxs.frameLineEdit.setText(str(len(self.selectedPilFramesNums)))
            self.uipdgioxs.show()
            self.connect(self.uipdgioxs.closePushButton, SIGNAL('clicked()'), self.pilGioxsClose)
            self.connect(self.uipdgioxs.runPushButton, SIGNAL('clicked()'), self.pilGioxsRun)
            self.connect(self.uipdgioxs.savePushButton, SIGNAL('clicked()'), self.pilGioxsSave)
            self.pilGioxsRun()
            
    def pilGioxsClose(self):
        self.uipdgioxs.close()
        
    def pilGioxsSave(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save gioxs data',directory=self.directory))
        for i in range(len(self.gioxsdata)):
            fname=self.saveFileName+self.gioxsname[i]
            np.savetxt(fname,self.gioxsdata[i],fmt='%.4f\t%.4e\t%.4e')
        fname=self.saveFileName+'_bg_gio.txt'
        np.savetxt(fname,self.gioxsbgdata,fmt='%.4f\t%.4e\t%.4e')


        
    def pilGioxsRun(self):
        self.uipdgioxs.plotWidget.canvas.ax.clear()
        aveframe=int(self.uipdgioxs.frameLineEdit.text())
        if aveframe<1:
            aveframe=max(np.abs(aveframe),1)
            self.messageBox('Warning: Average frames should be a positive interger!\n'+str(aveframe)+' is used as ave frames in the calculation')
        avepoint=int(self.uipdgioxs.pointLineEdit.text())
        if avepoint<1:
            avepoint=max(np.abs(avepoint),1)
            self.messageBox('Warning: Average points should be a positive interger!\n'+str(avepoints)+' is used as ave points in the calculation')
        ini=self.xcenter[0]-4
        fin=self.xcenter[0]+4
        cen=[self.xcenter[0], self.ycenter[0]]
        loop=[divmod(len(self.selectedPilFramesNums),aveframe)[0],divmod(len(self.selectedPilFramesNums),aveframe)[1]]
        if loop[1]!=0:
            loop[0]=loop[0]+1
        self.gioxsdata={}
        self.gioxsname={}
        self.pilatus.plotHint(self.bgData,self.bgError,hroi=[ini,fin],cen=cen,ax_type='Q', wavelength=self.selectedBgWavelength[0],s2d_dist=self.selectedBgDist[0],alpha=self.selectedBgAlpha[0],truealpha=self.selectedBgTrueAlpha[0])
        self.gioxsbgdata=self.pilatus.hintData
        for i in range(loop[0]):
            data=self.pilGIOXSData[self.selectedPilFramesNums[aveframe*i]]
            errordata=self.pilGIOXSErrorData[self.selectedPilFramesNums[aveframe*i]]**2
            if i==loop[0]-1 and loop[1]!=0: # last average 
                for j in range(aveframe*i+1,aveframe*i+loop[1]):
                    data=data+self.pilGIOXSData[self.selectedPilFramesNums[j]]
                    errordata=errordata+self.pilGIOXSErrorData[self.selectedPilFramesNums[j]]**2
            else:
                for j in range(aveframe*i+1,aveframe*(i+1)):
                    data=data+self.pilGIOXSData[self.selectedPilFramesNums[j]]
                    errordata=errordata+self.pilGIOXSErrorData[self.selectedPilFramesNums[j]]**2
            if i==loop[0]-1 and loop[1]!=0: # last average 
                data=data/loop[1]
                errordata=np.sqrt(errordata)/loop[1]
                self.pilatus.plotHint(data,errordata,hroi=[ini,fin],cen=cen,ax_type='Q', wavelength=self.pilSelected_Wavelength[0],s2d_dist=self.pilSelected_Dist[0],alpha=self.pilSelected_Alpha[0],truealpha=self.pilSelected_TrueAlpha[0])
                datasum=self.pilatus.hintData
                datagioxs=self.pilGioxsSub(datagioxs,self.gioxsbgdata)
               # print datasum
                self.gioxsdata[i]=self.aveData(datagioxs,avepoint)
                self.gioxsname[i]='_frame'+str(self.selectedPilFramesNums[aveframe*i])+'-'+str(self.selectedPilFramesNums[aveframe*i+loop[1]-1])+'_gio.txt'
                self.uipdgioxs.plotWidget.canvas.ax.errorbar(self.gioxsdata[i][:,0],self.gioxsdata[i][:,1],self.gioxsdata[i][:,2],fmt='o-',label='F #'+str(self.selectedPilFramesNums[aveframe*i])+'-'+str(self.selectedPilFramesNums[aveframe*i+loop[1]-1]))
            else:
                data=data/aveframe
                errordata=np.sqrt(errordata)/aveframe
                self.pilatus.plotHint(data,errordata,hroi=[ini,fin],cen=cen,ax_type='Q', wavelength=self.pilSelected_Wavelength[0],s2d_dist=self.pilSelected_Dist[0],alpha=self.pilSelected_Alpha[0],truealpha=self.pilSelected_TrueAlpha[0])
                datasum=self.pilatus.hintData
                datagioxs=self.pilGioxsSub(datasum,self.gioxsbgdata)
               # print datasum
                self.gioxsdata[i]=self.aveData(datagioxs,avepoint)
                #print datasum[:,0]
                #print self.gioxsdata[i][:,0]
                self.gioxsname[i]='_frame'+str(self.selectedPilFramesNums[aveframe*i])+'-'+str(self.selectedPilFramesNums[aveframe*(i+1)-1])+'_gio.txt'
                self.uipdgioxs.plotWidget.canvas.ax.errorbar(self.gioxsdata[i][:,0],self.gioxsdata[i][:,1],self.gioxsdata[i][:,2],fmt='o-',label='F #'+str(self.selectedPilFramesNums[aveframe*i])+'-'+str(self.selectedPilFramesNums[aveframe*(i+1)-1]))
       # self.uipdgioxs.plotWidget.canvas.ax.clear()
       # self.uipdgioxs.plotWidget.canvas.ax.errorbar(bgdata[:,0],bgdata[:,1],bgdata[:,2],fmt='o-')
        self.uipdgioxs.plotWidget.canvas.ax.errorbar(self.gioxsbgdata[:,0],self.gioxsbgdata[:,1],self.gioxsbgdata[:,2],fmt='o-',label='background')
        self.uipdgioxs.plotWidget.canvas.ax.legend(frameon=False,scatterpoints=0,numpoints=1)
        self.uipdgioxs.plotWidget.canvas.ax.set_yscale('log')
        self.uipdgioxs.plotWidget.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.uipdgioxs.plotWidget.canvas.ax.set_ylabel('Intensity')
        self.uipdgioxs.plotWidget.canvas.ax.set_title('File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1]+'\nbg S# '+str([item for item in np.unique(np.sort(self.selectedBgScanNum))])[1:-1])
        self.uipdgioxs.plotWidget.canvas.draw()
        self.command='GIOXS, singal scans=['+str([item for item in np.sort(self.selectedScanNums)])[1:-1]+'], bg scans=['+str([item for item in np.unique(np.sort(self.selectedBgScanNum))])[1:-1]+'], ave_frames='+str(aveframe)+', ave_points='+str(avepoint)
        self.ui.commandLineEdit.setText(self.command)                
        
    def pilGioxsSub(self, data1, data2): #for Gioxs 1d data substration
        data1x=data1[:,0]
        data1y=data1[:,1]
        data1yerr=data1[:,2]
        data2x=data2[:,0]
        data2y=data2[:,1]
        data2yerr=data2[:,2]
        overlap1=np.where(np.logical_and(data1x<=data2x[-1],data1x>=data2x[0]))
        f=interp1d(data2x,data2y,kind='cubic')
        ferr=interp1d(data2x,data2yerr,kind='cubic')
        data3y=f(data1x[overlap1])
        data3yerr=ferr(data1x[overlap1])
        datanewx=data1x[overlap1]
        datanewy=data1y[overlap1]-data3y
        datanewyerr=np.sqrt(data1yerr[overlap1]**2+data3yerr**2)
        data=np.vstack((datanewx,datanewy,datanewyerr)).T
        #print data
        #print data1x
        #print data2x
        return data
    
    def aveData(self, data, avepoint):
        if avepoint==1:
            return data
        else:
            loop=[divmod(len(data),avepoint)[0],divmod(len(data),avepoint)[1]]
            newdata=np.array([0,0,0])
            for i in range(loop[0]):
                newdata=np.vstack((newdata,np.array([np.average(data[:,0][i*avepoint:(i+1)*avepoint]),np.average(data[:,1][i*avepoint:(i+1)*avepoint]),np.sqrt(np.sum(data[:,2][i*avepoint:(i+1)*avepoint]**2))/avepoint])))
            if loop[1]!=0:
                newdata=np.vstack((newdata,np.array([np.average(data[:,0][loop[0]*avepoint:loop[0]*avepoint+loop[1]]),np.average(data[:,1][loop[0]*avepoint:loop[0]*avepoint+loop[1]]),np.sqrt(np.sum(data[:,2][loop[0]*avepoint:loop[0]*avepoint+loop[1]]**2))/loop[1]])))
            return newdata[1:]
    
    def pilBurnTest(self):
        self.pilSelectedQxs=[self.pilFileQxs[i] for i in self.selectedPilFramesNums]
        if self.pilSelectedQxs.count(self.pilSelectedQxs[0])!=len(self.pilSelectedQxs) or self.pilSelectedQzs.count(self.pilSelectedQzs[0])!=len(self.pilSelectedQzs) or self.pilSelectedScanNum.count(self.pilSelectedScanNum[0])!=len(self.pilSelectedScanNum):
            self.messageBox('Warning: It is not a burn test scan or multiple scans are selected!')
        else:
            Dialog=QDialog(self)                
            self.uipdburn=uic.loadUi('pdburn.ui', Dialog)
            self.uipdburn.show()
            self.connect(self.uipdburn.gidComboBox, SIGNAL('activated(int)'),self.pilGidBurn)  #do GID burn analysis
            #self.connect(self.uipdburn.rangeLineEdit, SIGNAL('returnPressed()'),self.pilGidBurn)
            self.connect(self.uipdburn.refPushButton, SIGNAL('clicked()'),self.pilRefBurn) #do REF burn analysis
            self.connect(self.uipdburn.closePushButton,SIGNAL('clicked()'),self.pilBurnClose) #close the window 
        
    def pilBurnClose(self):
        self.uipdburn.close()
        
    def expdecay(self, x, a,b):  #define exponential decay
        return a*np.exp(-x/b)
        
        
    def pilRefBurn(self):
        if self.pilSelected_Dth[0] > 0.001:
            self.messageBox('Warning: It it not a REF burn test. Please use GID burn test button!')
        else:
            self.uipdburn.plotWidget.canvas.ax.clear()
            self.uipdburn.plotWidget.canvas.ax.set_xlabel('Time (sec.)')
            self.uipdburn.plotWidget.canvas.ax.set_ylabel('Intensity')
            ctime=float(self.scanlines[self.selectedScanNums[0]-1].split(' ')[-1])  #get the counting time
            scan=self.scanlines[self.selectedScanNums[0]-1].split(' ')[1] #get scan number
            self.uipdburn.Label.setText('Burn test for Scan #'+scan)
            slitx=int(self.ui.refSlitLineEdit.text().split(',')[0]) #get the size of ROI
            slity=int(self.ui.refSlitLineEdit.text().split(',')[1])
            self.slit=[slitx,slity]        
            self.dir=str(self.ui.refBGDirComboBox.currentText()) #get the background offset direction
            bg=float(self.ui.refBGOffLineEdit.text())  #get the background offset
            data=[]
            for i in self.selectedPilFramesNums:
                cenx=self.pilatus.NCOLS-self.pilX[i]-1 #image was flipped left to right
                ceny=self.pilY[i]
                cenx,ceny=self.pilatus.peakLocate(self.pilData[i],self.pilErrorData[i],cen=[cenx,ceny],mon=self.pilMonc[i])
                if self.dir=='H':
                    self.bg=bg*self.pil_Dist[i]*np.pi/180/0.172/slitx
                else:
                    self.bg=bg*self.pil_Dist[i]*np.pi/180/0.172/slity
                sig,sigerr,lbg,lbgerr,rbg,rbgerr=self.pilatus.sumROI(slit=self.slit,cen=[cenx,ceny],dir=self.dir,bg=self.bg)
                data.append([ctime*i,sig-(lbg+rbg)/2.0, np.sqrt(sigerr**2+(lbgerr**2+rbgerr**2)/4)])
            data=np.array(data)
            x=data[:,0]
            y=data[:,1]
            yerr=data[:,2]
            pfit,pcov=curve_fit(self.expdecay,x,y,p0=[y[0],(x[-1]-x[0])/np.log(y[0]/y[-1])],sigma=yerr,maxfev=5000)
            tauerr=int(pcov[1][1]**0.5)
            tau=int(pfit[1])
            self.uipdburn.plotWidget.canvas.ax.errorbar(x,y,yerr,fmt='o-')
            self.uipdburn.plotWidget.canvas.ax.plot(x,self.expdecay(x,pfit[0],pfit[1]),'r--',label=r'$\tau$'+' = '+str(format(tau,'.1e'))+' ('+str(format(tauerr,'.1e'))+') sec.')
            self.uipdburn.plotWidget.canvas.ax.legend(loc=1,frameon=False,scatterpoints=0,numpoints=1) 
            self.uipdburn.plotWidget.canvas.ax.set_xlim(x[0],x[-1])
            self.uipdburn.plotWidget.canvas.ax.set_title('File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1])
            self.uipdburn.plotWidget.canvas.draw()
            self.command='REF Burn, scan='+str(self.selectedScanNums)
            self.ui.commandLineEdit.setText(self.command)  
            
    def pilGidBurn(self):
        if self.pilSelected_Dth[0]==0:
            self.messageBox('Warning: It is not a GID burn test. Please click REF burn test button!')
        else:
            ctime=float(self.scanlines[self.selectedScanNums[0]-1].split(' ')[-1])  #get the counting time
            scan=self.scanlines[self.selectedScanNums[0]-1].split(' ')[1] #get scan number
            self.uipdburn.Label.setText('Burn test for Scan #'+scan)
            ini=float(self.uipdburn.rangeLineEdit.text().split(':')[0]) #get the  ini qz value from lineedit
            fin=float(self.uipdburn.rangeLineEdit.text().split(':')[1]) #get the final qz value from lineedit
            cen=self.pilatus.NCOLS-self.pilSelectedX[0]-1
            ini=max(int((2.0*np.arcsin(ini*self.pilSelected_Wavelength[0]/4.0/np.pi)-self.pilSelected_Dth[0])*self.pilSelected_PDist[0]/0.172)+cen,0) # conver it to pixel and check with the boundary 
            fin=min(int((2.0*np.arcsin(fin*self.pilSelected_Wavelength[0]/4.0/np.pi)-self.pilSelected_Dth[0])*self.pilSelected_PDist[0]/0.172)+cen,self.pilatus.NCOLS-1)
            row=np.arange(self.pilatus.NROWS)
            col=np.arange(ini,fin+1)
            col,row=np.meshgrid(col,row)
            data=[]
            self.uipdburn.plotWidget.canvas.ax.clear()
            self.uipdburn.plotWidget.canvas.ax.set_xlabel('Time (sec.)')
            for i in self.selectedPilFramesNums:
                self.pilatus.plotVint(self.pilData[i],self.pilErrorData[i],hroi=[ini,fin])
                datai=self.pilatus.vintData
                meanleft=np.mean(datai[0:3,1])  #pass a straight line bewteen first and last three points in the selected range, then subtracte it from the data, keep the errorbar unchange
                meanright=np.mean(datai[-3:,1])
                slope=(meanright-meanleft)/(datai[-2,0]-datai[1,0])
                datai[:,1]=datai[:,1]-slope*(datai[:,0]-datai[1,0])-meanleft
                #datai[:,2]=np.sqrt(np.abs(datai[:,1]))   #need to recheck the error bar
                if self.uipdburn.gidComboBox.currentText()=='Intensity':
                    data.append([ctime*i,np.sum(datai[:,1])/self.pilMonc[i],np.sqrt(np.sum(datai[:,1])**2/self.pilMonc[i]**3+np.sum(datai[:,2]**2)/self.pilMonc[i]**2)])  
                    self.uipdburn.plotWidget.canvas.ax.set_ylabel('Intensity')
                elif self.uipdburn.gidComboBox.currentText()=='Location':
                    data.append([ctime*i,np.sum(datai[:,1]*datai[:,0])/np.sum(datai[:,1]),np.sqrt(np.sum(datai[:,0]**2*datai[:,2]**2)/np.sum(datai[:,1])**2+np.sum(datai[:,1]*datai[:,0])**2*np.sum(datai[:,2]**2)/np.sum(datai[:,1])**4)])  
                elif self.uipdburn.gidComboBox.currentText()=='Width':
                    location=np.sum(datai[:,1]*datai[:,0])/np.sum(datai[:,1])
                    datai[:,0]=(datai[:,0]-location)**2
                    data.append([ctime*i,np.sum(datai[:,1]*datai[:,0])/np.sum(datai[:,1]),np.sqrt(np.sum(datai[:,0]**2*datai[:,2]**2)/np.sum(datai[:,1])**2+np.sum(datai[:,1]*datai[:,0])**2*np.sum(datai[:,2]**2)/np.sum(datai[:,1])**4)])  
            data=np.array(data)
            x=data[:,0]
            y=data[:,1]
            yerr=data[:,2]
            if self.uipdburn.gidComboBox.currentText()=='Location': 
                self.uipdburn.plotWidget.canvas.ax.set_ylabel('Relative Location '+r'$[\AA^{-1}]$')
                y=4.0*np.pi*np.sin((self.pilSelected_Dth[0]+(data[:,1]-cen)*0.172/self.pilSelected_PDist[0])/2.0)/self.pilSelected_Wavelength[0]
                yerr=4.0*np.pi*np.cos((self.pilSelected_Dth[0]+(data[:,1]-cen)*0.172/self.pilSelected_PDist[0])/2.0)/self.pilSelected_Wavelength[0]*data[:,2]*0.172/2/self.pilSelected_PDist[0]
                y=y-y[0]  
            if self.uipdburn.gidComboBox.currentText()=='Width': 
                self.uipdburn.plotWidget.canvas.ax.set_ylabel('FWHM '+r'$[\AA^{-1}]$')
                y=np.sqrt(y)
                yerr=0.5/y*yerr
                y=4.0*np.pi*np.sin((data[:,1]*0.172/self.pilSelected_PDist[0])/2.0)/self.pilSelected_Wavelength[0]
                yerr=4.0*np.pi*np.cos((data[:,1]*0.172/self.pilSelected_PDist[0])/2.0)/self.pilSelected_Wavelength[0]*data[:,2]*0.172/2/self.pilSelected_PDist[0]
            self.uipdburn.plotWidget.canvas.ax.errorbar(x,y,yerr,fmt='o-')
            if self.uipdburn.gidComboBox.currentText()=='Intensity':
                self.uipdburn.plotWidget.canvas.ax.set_ylabel('Intensity')
                pfit,pcov=curve_fit(self.expdecay,x,y,p0=[y[0],(x[-1]-x[0])/np.log(y[0]/y[-1])],sigma=yerr,maxfev=5000)
                tauerr=int(pcov[1][1]**0.5)
                tau=int(pfit[1])
                self.uipdburn.plotWidget.canvas.ax.plot(x,self.expdecay(x,pfit[0],pfit[1]),'r--',label=r'$\tau$'+' = '+str(format(tau,'.1e'))+' ('+str(format(tauerr,'.1e'))+') sec.')
                self.uipdburn.plotWidget.canvas.ax.legend(loc=1,frameon=False,scatterpoints=0,numpoints=1) 
            self.uipdburn.plotWidget.canvas.ax.set_xlim(x[0],x[-1])
            self.uipdburn.plotWidget.canvas.ax.set_title('File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1])
            self.uipdburn.plotWidget.canvas.draw()
            self.command='GID Burn, scan='+str(self.selectedScanNums)+', Qz range=['+str(self.uipdburn.rangeLineEdit.text())+'], mode='+str(self.uipdburn.gidComboBox.currentText())
            self.ui.commandLineEdit.setText(self.command)             
                 
    def pilGISAXS(self):
        if str(self.ui.gisaxsComboBox.currentText())=='Show GISAXS':
            self.showPilGISAXS()
        elif str(self.ui.gisaxsComboBox.currentText())=='Save GISAXS':
            self.savePilGISAXS()

    def showPilGISAXS(self):
        self.pilGISAXSshow=1
        self.pilGIDshow=0
        if len(self.selectedPilFramesNums)!=1 and self.ui.gixSumCheckBox.checkState()==0:
            self.messageBox('Warning: Please select only one frame or use summed frames')
        else:
            if self.ui.gixSumCheckBox.checkState() != 0:  # if "Sum All" is checked; merge the frames having the same qz
                self.sumData = self.pilData[self.selectedPilFramesNums[0]]  # self.pilData only has pixel and intensity;
                self.sumErrData = self.pilErrorData[self.selectedPilFramesNums[0]] ** 2
                for i in self.selectedPilFramesNums[1:]:
                    self.sumData = self.sumData + self.pilData[i]
                    self.sumErrData = self.sumErrData + self.pilErrorData[i] ** 2
                self.sumData = self.sumData / float(len(self.selectedPilFramesNums))  # average the intensity
                self.sumErrData = np.sqrt(self.sumErrData) / float(len(self.selectedPilFramesNums))  # average the error bar
            if str(self.ui.pilAxesComboBox.currentText())!='Q':
                self.disconnect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
                self.ui.pilAxesComboBox.setCurrentIndex(2)
                self.connect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
            if self.ui.pilSpecCheckBox.checkState() != 0:  # if this CB is checked, use spec value
                if self.checkSameArrays(self.pilSelectedX) == False or self.checkSameArrays(self.pilSelectedY) == False or self.checkSameArrays(self.pilSelected_Dist) == False:
                    self.messageBox('Warning:: The selected frames have different Pilatus parameters; the parameters of the first frames are used here!!')
                try:
                    self.ui.pilSDDistLineEdit.setText(str(self.pilSelected_Dist[0]))
                    self.ui.pilDBPosLineEdit.setText(str(int(self.pilatus.NCOLS-self.pilSelectedX[0])) + ',' + str(self.pilSelectedY[0]))
                    self.distance = np.array(self.pilSelected_Dist)
                    self.xcenter = np.array(self.pilatus.NCOLS- self.pilSelectedX)  # image was flipped from left to right, i.e., x=194
                    self.ycenter = np.array(self.pilSelectedY)
                except:
                    pass
            else:  # if this CB is unchecked, use input value
                self.distance = float(self.ui.pilSDDistLineEdit.text()) * np.ones_like(self.selectedPilFramesNums)
                self.xcenter = float(self.ui.pilDBPosLineEdit.text().split(',')[0]) * np.ones_like(self.selectedPilFramesNums)
                self.ycenter = float(self.ui.pilDBPosLineEdit.text().split(',')[-1]) * np.ones_like(self.selectedPilFramesNums)
            self.ui.pilMplWidget.canvas.fig.clf()
            cmap = str(self.ui.pilCMapComboBox.currentText())
            vmax = float(self.ui.pilMaxLineEdit.text())
            vmin = float(self.ui.pilMinLineEdit.text())
            pilX = np.arange(0, self.pilatus.NCOLS)
            pilY = np.arange(0, self.pilatus.NROWS)
            aspect = str(self.ui.pilAspectComboBox.currentText())
            xlabel = r'$Q_{xy}$' + ' ' + r'$[\AA^{-1}]$'
            ylabel = r'$Q_z$' + ' ' + r'$[\AA^{-1}]$'
            ax = self.ui.pilMplWidget.canvas.fig.add_subplot(1, 1, 1)
            ax.set_xlabel(xlabel, fontsize=16)
            ax.set_ylabel(ylabel, fontsize=16)
            pilY_deg = np.arctan(-(self.ycenter[0] - pilY) * 0.172 / self.distance[0]) + np.arcsin(2 * np.sin(self.alpha[0]) - np.sin(self.truealpha[0]))  # beta in rad
            if self.ui.pilPinholeCheckBox.checkState()!= 0:
                self.distance = np.array(self.pilSelected_PDist)  # distance between S4 & S5 for pinhole setup
            pilX_deg = self.dth[0] + np.arctan(-(self.xcenter[0] - pilX) * 0.172 / self.distance[0])  # dth in rad
            pilX_deg, pilY_deg = np.meshgrid(pilX_deg, pilY_deg)  # mesh x-y corrdinates
            self.pilXbin = 2.0 * np.pi / self.wavelength[0] * np.sqrt(np.cos(self.truealpha[0]) ** 2 + np.cos(pilY_deg) ** 2 - 2 * np.cos(self.truealpha[0]) * np.cos(pilY_deg) * np.cos(pilX_deg))
            self.pilYbin = 2.0 * np.pi / self.wavelength[0] * (np.sin(self.truealpha[0]) + np.sin(pilY_deg))
            qxyrange = np.linspace(np.min(self.pilXbin), np.max(self.pilXbin), len(self.pilXbin[0]) * (np.max(self.pilXbin) - np.min(self.pilXbin)) / (
                                           self.pilXbin[0][-1] - self.pilXbin[0][0]))
            qxy, qz = np.meshgrid(qxyrange, self.pilYbin[:, 0])
            if len(self.selectedPilFramesNums)==1:
                 index = self.selectedPilFramesNums[0]
                 self.sumData=self.pilData[index]
                 self.sumErrData=self.pilErrorData[index]
            self.pilDataBin = griddata((self.pilXbin.ravel(), self.pilYbin.ravel()), self.sumData.ravel(), (qxy, qz), method='linear', fill_value=1e-10)
            self.pilErrorDataBin = griddata((self.pilXbin.ravel(), self.pilYbin.ravel()), self.sumErrData.ravel(), (qxy, qz), method='linear', fill_value=1e-10)
            self.pilXbin=qxy
            self.pilYbin=qz
            self.extent=[self.pilXbin[0][0],self.pilXbin[0][-1],self.pilYbin[0][0],self.pilYbin[-1][0]]
            if self.ui.pilLogIntCheckBox.checkState() != 0:
                self.logData = np.log10(np.where(self.pilDataBin <= 0, 1e-10, self.pilDataBin))
                p = ax.imshow(self.logData, interpolation='bicubic', extent=self.extent, vmax=np.log10(vmax),vmin=np.log10(vmin), cmap=cmap, aspect=aspect, origin='lower')
                self.Zdata = self.logData
            else:
                p = ax.imshow(self.pilDataBin, interpolation='bicubic', extent=self.extent, vmax=vmax, vmin=vmin, cmap=cmap,aspect=aspect, origin='lower')
                self.Zdata = self.pilDataBin
            ax.set_aspect(aspect)
            ax.set_title("GISAXS")
            ax.format_coord=self.format_pil_coord
            self.ui.pilMplWidget.canvas.fig.colorbar(p)
            self.ui.pilMplWidget.canvas.fig.tight_layout()
            self.ui.pilMplWidget.canvas.fig.subplots_adjust(top=0.9)
            self.ui.pilMplWidget.canvas.draw()

    def savePilGISAXS(self):
        try:
            self.saveFileName = str(QFileDialog.getSaveFileName(caption='Save gioxs data', directory=self.directory))
            data2d = []
            fname_2d = self.saveFileName + '_gisaxs.txt'
            for i in range(len(self.pilDataBin)):
                for j in range(len(self.pilDataBin[0])):
                    data2d.append([self.pilXbin[i][j], self.pilYbin[i][j], self.pilDataBin[i][j], self.pilErrorDataBin[i][j]])
            np.savetxt(fname_2d, np.array(data2d), fmt='%.4e\t%.4e\t%.4e\t%.4e')
        except:
            self.messageBox("Warning: No GISAXS data available!")

    def pilGID(self):
        if str(self.ui.gidComboBox.currentText())=='GID Patch': 
            self.pilGIDRelax=0
            self.pilGIDpatchData()
        elif str(self.ui.gidComboBox.currentText())=='GID Patch Relax':
            self.pilGIDRelax=1
            self.pilGIDpatchData()
        elif str(self.ui.gidComboBox.currentText())=='Save GID':
            self.saveGIDData()
        elif str(self.ui.gidComboBox.currentText())=='Form Factor':
            self.setGIDFormFactor()
            
    def setGIDFormFactor(self):
        Dialog=QDialog(self)                
        self.uipdformfac=uic.loadUi('pdformfac.ui', Dialog)
        self.uipdformfac.show()
        self.connect(self.uipdformfac.updatePushButton, SIGNAL('clicked()'),self.updateFormFactor)  #update gid plot with formfactor
        self.connect(self.uipdformfac.okPushButton, SIGNAL('clicked()'),self.presetFormFactor) #preset formfactor
        
    def presetFormFactor(self):
        if self.uipdformfac.sphCheckBox.checkState()!=0 and self.uipdformfac.ellCheckBox.checkState()!=0:
            self.messageBox('Warning: Both "Sphere" and "Ellipsoid" are selected.\n Please uncheck one of them at least.')
        elif self.uipdformfac.sphCheckBox.checkState()!=0:
            self.formfactor_rxy=0
            self.formfactor_rz=0
            self.formfactor_r=float(self.uipdformfac.sphLineEdit.text())
            if self.formfactor_r<=0:
                self.messageBox('Warning: The radius of the sphere should be positive!')
            else:
                self.uipdformfac.close()
        elif self.uipdformfac.ellCheckBox.checkState()!=0:
            self.formfactor_r=0
            self.formfactor_rxy=float(self.uipdformfac.ellXYLineEdit.text())
            self.formfactor_rz=float(self.uipdformfac.ellZLineEdit.text())
            if self.formfactor_rxy<=0 or self.formfactor_rz<=0:
                self.messageBox('Warning: The radius of the ellipsoid should be positive!')
            else:
                self.uipdformfac.close()
        else:
            self.formfactor_r=0
            self.formfactor_rxy=0
            self.formfactor_rz=0
            self.uipdformfac.close()
            
    def updateFormFactor(self):
        try:
            temp=len(self.pilGIDData)
            if self.uipdformfac.sphCheckBox.checkState()!=0 and self.uipdformfac.ellCheckBox.checkState()!=0:
                self.messageBox('Warning: Both "Sphere" and "Ellipsoid" are selected.\n Please uncheck one of them at least.')
            elif self.uipdformfac.sphCheckBox.checkState()!=0:
                self.formfactor_rxy=0
                self.formfactor_rz=0
                self.formfactor_r=float(self.uipdformfac.sphLineEdit.text())
                if self.formfactor_r<=0:
                    self.messageBox('Warning: The radius of the sphere should be positive!')
                else:
                    self.pilGIDPlot()
            elif self.uipdformfac.ellCheckBox.checkState()!=0:
                self.formfactor_r=0
                self.formfactor_rxy=float(self.uipdformfac.ellXYLineEdit.text())
                self.formfactor_rz=float(self.uipdformfac.ellZLineEdit.text())
                if self.formfactor_rxy<=0 or self.formfactor_rz<=0:
                    self.messageBox('Warning: The radius of the ellipsoid should be positive!')
                else:
                    self.pilGIDPlot()
            else:
                self.formfactor_r=0
                self.formfactor_rxy=0
                self.formfactor_rz=0
                self.pilGIDPlot()
        except:
            self.messageBox('Waning: GID patch data does not exist!\n Please do the GID patch first')
            
    
    def pilGIDpatchData(self):
        if self.pilGIDshow==0:
            self.disconnect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
            if self.ui.pilAxesComboBox.currentIndex()==0:
                self.ui.pilAxesComboBox.setCurrentIndex(2)
            self.connect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
        self.pilGIDshow=1
        self.pilGISAXSshow=0
        self.ui.pilMplWidget.canvas.fig.clf()
        self.pilxleft=self.xcenter-(float(self.ui.pilHSlitLineEdit.text())-1)/2
        self.pilxright=self.xcenter+(float(self.ui.pilHSlitLineEdit.text())-1)/2
        self.pilhintLQDatadic={}  #create the dictionary for low qz 
        self.pilhintMQDatadic={}  #create the dictionary for medium qz
        self.pilhintHQDatadic={}  #create the dictionary for high qz
        areainfo={}  # create a dictinary for area cal info which are contants for each scan
        averagemonc={} #create the dicionary for average monc for each scan
        qz=np.sort(list(set([float(format(self.pilSelectedQzs[i], '.3f')) for i in range(len(self.pilSelectedQzs))]))) # get and sort the qz in the selected frames
        self.absfac=float(self.ui.AbsFacLineEdit.text())
        j=0   #for xcenter and ycenter, which are already defined in self.selectedPilFramesNums
        for i in self.selectedPilFramesNums:  # get average monc for each scan
            ckey=self.pilScanNum[i]
            if ckey in averagemonc:
                averagemonc[ckey].append(self.pilMonc[i])
            else: 
                averagemonc[ckey]=[self.pilMonc[i]]
        #print averagemonc, ckey
        #print np.average(averagemonc[ckey])
        
        for i in self.selectedPilFramesNums:
            self.pilatus.plotHint({i:self.pilData[i]},{i:self.pilErrorData[i]},absnum=[self.pil_AbsNum[i]],cen=[self.xcenter[j],self.ycenter[j]], hroi=[self.pilxleft[j], self.pilxright[j]], vroi=[5,486], ax_type='Angles', wavelength=[self.pil_Wavelength[i]],s2d_dist=[self.pil_Dist[i]], sh=[self.pil_Sh[i]],alpha=self.pilAlpha[i], truealpha=self.pilTrueAlpha[i], mon=[self.pilMonc[i]])
            j=j+1
            ckey=self.pilScanNum[i]
            if ckey in areainfo:
                pass
            else:
             #   print ckey, self.pilTrueAlpha[i]
                areainfo[ckey]=[self.pil_S1h[i],self.pil_S1v[i],self.pil_S4h[i],self.pil_S5h[i],self.pilTrueAlpha[i], self.pil_Dist[i], self.pil_PDist[i]]   #s1h, s1v, s4h, alpha, sample to detector, slit to detector
            dth=float(format(self.pil_Dth[i]*180/np.pi,'.3f'))**np.ones_like(self.pilatus.hintData[:,0])            
            if np.abs(self.pilFileQzs[i]-qz[0])<0.001:  #low qz
                if ckey in self.pilhintLQDatadic:
                    self.pilhintLQDatadic[ckey].append(np.vstack((dth,self.pilatus.hintData[:,0],self.pilatus.hintData[:,1]*np.average(averagemonc[ckey]),self.pilatus.hintData[:,2]*np.average(averagemonc[ckey]))).transpose())
                else:
                    self.pilhintLQDatadic[ckey]=[np.vstack((dth,self.pilatus.hintData[:,0],self.pilatus.hintData[:,1]*np.average(averagemonc[ckey]),self.pilatus.hintData[:,2]*np.average(averagemonc[ckey]))).transpose()]
            elif np.abs(self.pilFileQzs[i]-qz[1])<0.001: #medium qz
                if ckey in self.pilhintMQDatadic:
                    self.pilhintMQDatadic[ckey].append(np.vstack((dth,self.pilatus.hintData[:,0],self.pilatus.hintData[:,1]*np.average(averagemonc[ckey]),self.pilatus.hintData[:,2]*np.average(averagemonc[ckey]))).transpose())
                else:
                    self.pilhintMQDatadic[ckey]=[np.vstack((dth,self.pilatus.hintData[:,0],self.pilatus.hintData[:,1]*np.average(averagemonc[ckey]),self.pilatus.hintData[:,2]*np.average(averagemonc[ckey]))).transpose()]
            else:  #high qz
                if ckey in self.pilhintHQDatadic:
                    self.pilhintHQDatadic[ckey].append(np.vstack((dth,self.pilatus.hintData[:,0],self.pilatus.hintData[:,1]*np.average(averagemonc[ckey]),self.pilatus.hintData[:,2]*np.average(averagemonc[ckey]))).transpose())
                else:
                    self.pilhintHQDatadic[ckey]=[np.vstack((dth,self.pilatus.hintData[:,0],self.pilatus.hintData[:,1],self.pilatus.hintData[:,2])).transpose()]
        LQkeys=self.sortedScans(self.pilhintLQDatadic) #sorted scan numbers for low, middle, high qz. 
        MQkeys=self.sortedScans(self.pilhintMQDatadic)
        HQkeys=self.sortedScans(self.pilhintHQDatadic)
       # print areainfo, LQkeys, MQkeys
       # print self.pilhintLQDatadic
       # print self.pilhintLQDatadic[LQkeys[0]][0][0][0]
      
      # Area and polarization corrections;
        for i in range(len(LQkeys)):
            for j in range(len(self.pilhintLQDatadic[LQkeys[i]])):
                areacorr=self.calArea(areainfo[LQkeys[i]], self.pilhintLQDatadic[LQkeys[i]][j][:,0][0])
                self.pilhintLQDatadic[LQkeys[i]][j][:,2]=self.pilhintLQDatadic[LQkeys[i]][j][:,2]/areacorr/(1-(np.sin(self.pilhintLQDatadic[LQkeys[i]][j][:,0]*np.pi/180)*np.cos(self.pilhintLQDatadic[LQkeys[i]][j][:,1]*np.pi/180))**2)
                self.pilhintLQDatadic[LQkeys[i]][j][:,3]=self.pilhintLQDatadic[LQkeys[i]][j][:,3]/areacorr/(1-(np.sin(self.pilhintLQDatadic[LQkeys[i]][j][:,0]*np.pi/180)*np.cos(self.pilhintLQDatadic[LQkeys[i]][j][:,1]*np.pi/180))**2)
        for i in range(len(MQkeys)):
            for j in range(len(self.pilhintMQDatadic[MQkeys[i]])):
                areacorr=self.calArea(areainfo[MQkeys[i]], self.pilhintMQDatadic[MQkeys[i]][j][:,0][0])
                self.pilhintMQDatadic[MQkeys[i]][j][:,2]=self.pilhintMQDatadic[MQkeys[i]][j][:,2]/areacorr/(1-(np.sin(self.pilhintMQDatadic[MQkeys[i]][j][:,0]*np.pi/180)*np.cos(self.pilhintMQDatadic[MQkeys[i]][j][:,1]*np.pi/180))**2)
                self.pilhintMQDatadic[MQkeys[i]][j][:,3]=self.pilhintMQDatadic[MQkeys[i]][j][:,3]/areacorr/(1-(np.sin(self.pilhintMQDatadic[MQkeys[i]][j][:,0]*np.pi/180)*np.cos(self.pilhintMQDatadic[MQkeys[i]][j][:,1]*np.pi/180))**2)
        for i in range(len(HQkeys)):
            for j in range(len(self.pilhintHQDatadic[HQkeys[i]])):
                areacorr=self.calArea(areainfo[HQkeys[i]], self.pilhintHQDatadic[HQkeys[i]][j][:,0][0])
                self.pilhintHQDatadic[HQkeys[i]][j][:,2]=self.pilhintHQDatadic[HQkeys[i]][j][:,2]/areacorr/(1-(np.sin(self.pilhintHQDatadic[HQkeys[i]][j][:,0]*np.pi/180)*np.cos(self.pilhintHQDatadic[HQkeys[i]][j][:,1]*np.pi/180))**2)
                self.pilhintHQDatadic[HQkeys[i]][j][:,3]=self.pilhintHQDatadic[HQkeys[i]][j][:,3]/areacorr/(1-(np.sin(self.pilhintHQDatadic[HQkeys[i]][j][:,0]*np.pi/180)*np.cos(self.pilhintHQDatadic[HQkeys[i]][j][:,1]*np.pi/180))**2)
       # print self.pilhintLQDatadic
        if len(LQkeys)>1: #for low qz
            self.pilLQPatchData=self.imagehjoin(self.pilhintLQDatadic[LQkeys[0]],self.pilhintLQDatadic[LQkeys[1]]) #horizontally patch the first two images
            for i in range(2,len(LQkeys)):
                self.pilLQPatchData=self.imagehjoin(self.pilLQPatchData,self.pilhintLQDatadic[LQkeys[i]])  #horizontally patch other images at low qz if necessary 
        else:
            self.pilLQPatchData=self.pilhintLQDatadic[LQkeys[0]]
        if len(MQkeys)>1:  #for middle qz
            self.pilMQPatchData=self.imagehjoin(self.pilhintMQDatadic[MQkeys[0]],self.pilhintMQDatadic[MQkeys[1]]) #horizontally patch the first two images
            for i in range(2,len(MQkeys)):
                self.pilMQPatchData=self.imagehjoin(self.pilMQPatchData,self.pilhintMQDatadic[MQkeys[i]])  #horizontally patch other images at low qz if necessary 
        elif len(MQkeys)==1:
            self.pilMQPatchData=self.pilhintMQDatadic[MQkeys[0]]
        else:
            self.pilMQPatchData=[]
        if len(HQkeys)>1: #for middle qz
            self.pilHQPatchData=self.imagehjoin(self.pilhintHQDatadic[HQkeys[0]],self.pilhintHQDatadic[HQkeys[1]]) #horizontally patch the first two images
            for i in range(2,len(HQkeys)):
                self.pilHQPatchData=self.imagehjoin(self.pilHQPatchData,self.pilhintHQDatadic[HQkeys[i]])  #horizontally patch other images at low qz if necessary 
        elif len(HQkeys)==1:
            self.pilHQPatchData=self.pilhintHQDatadic[MQkeys[0]]
        else:
            self.pilHQPatchData=[]
        self.pilPatchData=self.pilLQPatchData # need vertical patch later on
        if len(self.pilMQPatchData)!=0:
            self.pilPatchData=self.imagevjoin(self.pilPatchData,self.pilMQPatchData)
        if len(self.pilHQPatchData)!=0:
            self.pilPatchData=self.imagevjoin(self.pilPatchData,self.pilHQPatchData)
        
        pilGIDXAxis=[self.pilPatchData[i][0][0] for i in range(len(self.pilPatchData))]  #get x-axis
        pilGIDYAxis=list(self.pilPatchData[0][:,1])  #get y-axis
        self.pilGIDXAxs, self.pilGIDYAxs=np.meshgrid(pilGIDXAxis,pilGIDYAxis) # mesh the x-y corrdinates
        self.pilGIDData=np.array([np.hstack(self.pilPatchData[i][:,2]) for i in range(len(self.pilPatchData))]).T  # get z-axis 
        self.pilGIDDataErr=np.array([np.hstack(self.pilPatchData[i][:,3]) for i in range(len(self.pilPatchData))]).T  # get error bars
        self.imageMax=np.max(self.pilGIDData)
        #self.imageMin=max(1e-10,np.min(self.pilGIDData))
        self.imageMin=np.min(self.pilGIDData[np.nonzero(self.pilGIDData)])
        self.ui.pilMaxLineEdit.setText(str(self.imageMax)) #set the vmax/vmin in the lineedit and slider 
        self.ui.pilMinLineEdit.setText(str(self.imageMin))
        self.ui.pilMaxHorizontalSlider.setValue(int(float(self.ui.pilMaxLineEdit.text())*100/(self.imageMax-self.imageMin)))
        self.ui.pilMinHorizontalSlider.setValue(int(float(self.ui.pilMinLineEdit.text())*100/(self.imageMax-self.imageMin)))
        self.pilGIDPlot()   # go to plot function

    def calArea(self, areainfo, dth):  #calculate the area correction for pilatum100 GID measurement
        s1y=areainfo[0] # footrprint size perpendicular to the beam, y
        s4h=areainfo[2]
        s5h=areainfo[3]
        alpha=areainfo[4]
        s1x=areainfo[1]/np.sin(alpha)  #footprint size along the beam, x        
        sam2d=areainfo[5]-32.0   # sample to detector distance; there is 32 mm shorter if we use hard slits. 
        sli2s=areainfo[5]-areainfo[6]   # s2 slits to sample distance
        dth=dth*np.pi/180.0   #convert to rad. 
        if s5h> float(self.ui.pilHSlitLineEdit.text())*0.172:    
            s5h=float(self.ui.pilHSlitLineEdit.text())*0.172
            sam2d=areainfo[5]
        phi=np.arcsin((s4h+s5h)/2.0/np.sqrt((s4h+s5h)**2/4.0+(sam2d-sli2s)**2))
        dist=(s4h/2.0*sam2d+s5h/2.0*sli2s)/np.sqrt((s4h+s5h)**2/4.0+(sam2d-sli2s)**2)
      #  print s1y, s1x, s4h, alpha, s5h, sam2d, sli2s, dth
        seg=[]
        seg.append([-s1x/2.0, s1y/2.0,-s1x/2.0,-s1y/2.0])    #segment 1 of the footprint, (x1, x2, y1, y2)
        seg.append([-s1x/2.0,-s1y/2.0, s1x/2.0,-s1y/2.0])    #segment 2 of the footprint, (x1, x2, y1, y2)
        seg.append([ s1x/2.0,-s1y/2.0, s1x/2.0, s1y/2.0])    #segment 3 of the footprint, (x1, x2, y1, y2)
        seg.append([ s1x/2.0, s1y/2.0,-s1x/2.0, s1y/2.0])    #segment 4 of the footprint, (x1, x2, y1, y2)
        
        conners=[[-s1x/2,-s1y/2], [s1x/2,-s1y/2],[s1x/2,s1y/2],[-s1x/2,s1y/2]]  #conners for the footprint
        pointb=[]         #cross point for the bottom line
        pointt=[]         #cross point for the top line

        
        for i in range(4):  # get the cross points for both top and bottom lines
            if seg[i][0]==seg[i][2]:    #vertical line  
                xtemp=seg[i][0]
                ybtemp=(xtemp*np.sin(dth-phi)-dist)/np.cos(dth-phi)
                yttemp=(xtemp*np.sin(dth+phi)+dist)/np.cos(dth+phi)
                if ybtemp>-s1y/2 and ybtemp<s1y/2:
                    pointb.append([xtemp, ybtemp])
                else:
                    pointb.append('None')
                if yttemp>-s1y/2 and yttemp<s1y/2:
                    pointt.append([xtemp, yttemp])
                else:
                    pointt.append('None')
            else:      #horizontal line
                ytemp=seg[i][1]
                xbtemp=(ytemp*np.cos(dth-phi)+dist)/np.sin(dth-phi)
                xttemp=(ytemp*np.cos(dth+phi)-dist)/np.sin(dth+phi)
                if xbtemp>-s1x/2 and xbtemp<s1x/2:
                    pointb.append([xbtemp, ytemp])
                else:
                    pointb.append('None')
                if xttemp>-s1x/2 and xttemp<s1x/2:
                    pointt.append([xttemp, ytemp])
                else:
                    pointt.append('None')
                  
        listb=[]
        listt=[] 
        for i in range(4):  #get segment numbers from top and bottom line
            if pointb[i]!='None':
                listb.append(i)
            if pointt[i]!='None':
                listt.append(i)
                # print listb, listt
        bpop=[]    # list for the footprint conners need to be poped out
        if len(listb)!=0:  #get conner points need to be poped out 
            bpop=range(listb[0],listb[1])
        if len(listt)!=0:
            tpop=range(listt[1],listt[0]+4)
            for i in range(len(tpop)):
                if tpop[i]>3:
                    bpop.append(tpop[i]-4)
                else:
                    bpop.append(tpop[i])
        bpop.sort(reverse=True)   
        for i in bpop:    #pop selected conner points
            conners.pop(i)
        for i in range(len(listb)):  #add points from bottom line    
            conners.append(pointb[listb[i]])
        for i in range(len(listt)):  #add points from top line
            conners.append(pointt[listt[i]])
        conangle=np.array([np.arctan2(conners[i][1],conners[i][0]) for i in range(len(conners))])  # get angles from the each point on the polygon 
        conorder=[conners[i] for i in np.argsort(conangle)]     # sort the points for the area calculation

        area=0  # area for the polygon
        for i in range(len(conorder)):
            if i==len(conorder)-1:
                area=area+conorder[i][0]*conorder[0][1]-conorder[i][1]*conorder[0][0]
            else:
                area=area+conorder[i][0]*conorder[i+1][1]-conorder[i][1]*conorder[i+1][0]
        #print np.abs(area)/2/s1x/s1y,dth
       # print dth, np.abs(area)/2.0
        return np.abs(area)/2/s1x/s1y  #return the normalized polygon area 

        
    def pilGIDPlot(self):
        self.ui.pilMplWidget.canvas.fig.clf()
        aspect=str(self.ui.pilAspectComboBox.currentText())
        vmax=float(self.ui.pilMaxLineEdit.text())
        vmin=float(self.ui.pilMinLineEdit.text())
        cmap=str(self.ui.pilCMapComboBox.currentText())
        #self.gidax=self.ui.pilMplWidget.canvas.fig.add_subplot(1,1,1)
        gs=gridspec.GridSpec(5,12)
        self.gidax=self.ui.pilMplWidget.canvas.fig.add_subplot(gs[0:3,3:7])  #for 2D data
        self.gidax_h=self.ui.pilMplWidget.canvas.fig.add_subplot(gs[3:5,3:7]) #for h-cut data
        self.gidax_v=self.ui.pilMplWidget.canvas.fig.add_subplot(gs[0:3,7:9])  #for v-cut data
        self.gidax_vc=self.ui.pilMplWidget.canvas.fig.add_subplot(gs[3:5,7:9])  #for vslide-cut data
        alpha=self.truealpha[0]  #get alpha
        k0=2*np.pi/self.wavelength[0] #get wavevector    
        
        if str(self.ui.pilAxesComboBox.currentText())=='Angles':   
            self.pilGIDxlabel=r'$2\theta$'+' '+r'$[\deg]$'
            self.pilGIDylabel=r'$\beta$'+' '+r'$[\deg]$'
            fint=interp1d(self.pilGIDXAxs[0], self.pilGIDData, kind='linear')   # interpolate the data along x-axis
            newpilGIDXAxis=np.linspace(self.pilGIDXAxs[0][0],self.pilGIDXAxs[0][-1],len(self.pilGIDXAxs[0])*10)
            inter_pilGIDData=fint(newpilGIDXAxis)
            self.Zdata=inter_pilGIDData
            self.extent=[self.pilGIDXAxs[0][0],self.pilGIDXAxs[0][-1],self.pilGIDYAxs[0][0],self.pilGIDYAxs[-1][-1]]
            if self.ui.pilLogIntCheckBox.checkState()!=0:
                self.gid_p=self.gidax.imshow(np.log10(inter_pilGIDData),interpolation='bicubic',extent=self.extent,vmax=np.log10(vmax),vmin=np.log10(vmin),cmap=cmap,aspect=aspect,origin='lower')            
            else:
                self.gid_p=self.gidax.imshow(inter_pilGIDData,interpolation='bicubic',extent=self.extent,vmax=vmax,vmin=vmin,cmap=cmap,aspect=aspect,origin='lower') 
           # self.ui.pilMplWidget.canvas.fig.colorbar(self.gid_p)
            #get the h-cut data and label
            self.gid_h=np.vstack((self.pilGIDXAxs[0,:],np.sum(self.pilGIDData,axis=0),np.sqrt(np.sum(self.pilGIDDataErr**2,axis=0)))).transpose()
            self.gidax_h.set_xlabel(r'$2\theta$'+' '+r'$[\deg]$')
            self.gidax_h.set_title('Integrated over '+r'$\beta$')
            #get the v-cut data
            self.gid_v=np.vstack((self.pilGIDYAxs[:,0],np.sum(self.pilGIDData,axis=1),np.sqrt(np.sum(self.pilGIDDataErr**2,axis=1)))).transpose()
            self.gidax_v.set_ylabel(r'$\beta$'+' '+r'$[\deg]$')
            self.gidax_v.set_title('Integrated over '+r'$2\theta$')
            #get the v-slice-cut data
            vc_ran=np.linspace(np.where(self.pilGIDYAxs[:,0]>0)[0][0],np.where(self.pilGIDYAxs[:,0]>0)[0][-1],5)
            vc_ran=map(int,vc_ran)
            self.gid_vc={}
            for i in range(4):
                self.gid_vc[i]=np.vstack((self.pilGIDXAxs[0,:],np.sum(self.pilGIDData[vc_ran[i]:vc_ran[i+1],:],axis=0),np.sqrt(np.sum(self.pilGIDDataErr[vc_ran[i]:vc_ran[i+1],:]**2,axis=0)))).transpose()
            self.gidax_vc.set_title(r'$\beta$'+'-binned')
            self.gidax_vc.set_xlabel(r'$2\theta$'+' '+r'$[\deg]$')
            
        
        if str(self.ui.pilAxesComboBox.currentText())=='Q':
            self.pilGIDxlabel=r'$Q_{xy}$'+' '+r'$[\AA^{-1}]$'
            self.pilGIDylabel=r'$Q_z$'+' '+r'$[\AA^{-1}]$'
            self.pilGIDXAxs_Q=k0*np.sqrt(np.cos(alpha)**2+np.cos(self.pilGIDYAxs/180*np.pi)**2-2*np.cos(alpha)*np.cos(self.pilGIDYAxs/180*np.pi)*np.cos(self.pilGIDXAxs/180*np.pi))
            self.pilGIDYAxs_Q=k0*(np.sin(alpha)+np.sin(self.pilGIDYAxs/180*np.pi))
            if self.pilGIDRelax==0:
                qxyrange=np.where(np.logical_and(self.pilGIDXAxs_Q[0]<self.pilGIDXAxs_Q[-1][-1],self.pilGIDXAxs_Q[0]>self.pilGIDXAxs_Q[-1][0]))  #find the qxy range in the smallest rectangular grids. 
                qxy,qz=np.meshgrid(self.pilGIDXAxs_Q[0][qxyrange],self.pilGIDYAxs_Q[:,0])
                self.pilGIDData_Q=griddata((self.pilGIDXAxs_Q.ravel(),self.pilGIDYAxs_Q.ravel()),self.pilGIDData.ravel(), (qxy,qz), method='linear')
                self.pilGIDDataErr_Q=griddata((self.pilGIDXAxs_Q.ravel(),self.pilGIDYAxs_Q.ravel()),self.pilGIDDataErr.ravel(), (qxy,qz), method='linear')
            else:
               # print  self.pilGIDXAxs_Q[0][0], self.pilGIDXAxs_Q[0][-1], self.pilGIDXAxs_Q[-1][0], self.pilGIDXAxs_Q[-1][-1] 
                #print np.min(self.pilGIDXAxs_Q), np.max(self.pilGIDXAxs_Q)
                qxyrange=np.linspace(np.min(self.pilGIDXAxs_Q), np.max(self.pilGIDXAxs_Q),len(self.pilGIDXAxs_Q[0])*(np.max(self.pilGIDXAxs_Q)-np.min(self.pilGIDXAxs_Q))/(self.pilGIDXAxs_Q[0][-1]-self.pilGIDXAxs_Q[0][0]))
               # print len(self.pilGIDXAxs_Q[0]),  len(self.pilGIDXAxs_Q[0])*(np.max(self.pilGIDXAxs_Q)-np.min(self.pilGIDXAxs_Q))/(self.pilGIDXAxs_Q[0][-1]-self.pilGIDXAxs_Q[0][0])
               # print qxyrange
                #qxyrange=np.where(np.logical_and(self.pilGIDXAxs_Q[0]<self.pilGIDXAxs_Q[-1][-1],self.pilGIDXAxs_Q[0]>self.pilGIDXAxs_Q[-1][0]))  #find the qxy range in the biggest rectangular grids. (relax mode)
                qxy,qz=np.meshgrid(qxyrange,self.pilGIDYAxs_Q[:,0])
                self.pilGIDData_Q=griddata((self.pilGIDXAxs_Q.ravel(),self.pilGIDYAxs_Q.ravel()),self.pilGIDData.ravel(), (qxy,qz), method='linear',fill_value=1e-10)
                self.pilGIDDataErr_Q=griddata((self.pilGIDXAxs_Q.ravel(),self.pilGIDYAxs_Q.ravel()),self.pilGIDDataErr.ravel(), (qxy,qz), method='linear',fill_value=1e-10)
            self.pilGIDXAxs_Q=qxy
            self.pilGIDYAxs_Q=qz
            fint=interp1d(self.pilGIDXAxs_Q[0], self.pilGIDData_Q, kind='linear')   # interpolate the data along x-axis
            newpilGIDXAxis=np.linspace(self.pilGIDXAxs_Q[0][0],self.pilGIDXAxs_Q[0][-1],len(self.pilGIDXAxs_Q[0])*10)
            inter_pilGIDData_Q=fint(newpilGIDXAxis)
            self.Zdata=inter_pilGIDData_Q
            self.extent=[self.pilGIDXAxs_Q[0][0],self.pilGIDXAxs_Q[0][-1],self.pilGIDYAxs_Q[0][0],self.pilGIDYAxs_Q[-1][-1]]
            if self.ui.pilLogIntCheckBox.checkState()!=0:
                self.gid_p=self.gidax.imshow(np.log10(inter_pilGIDData_Q),interpolation='bicubic',extent=self.extent,vmax=np.log10(vmax),vmin=np.log10(vmin),cmap=cmap,aspect=aspect,origin='lower')            
            else:
                self.gid_p=self.gidax.imshow(inter_pilGIDData_Q,interpolation='bicubic',extent=self.extent,vmax=vmax,vmin=vmin,cmap=cmap,aspect=aspect,origin='lower')            
           
            if self.formfactor_r>0:  #for the sphere ff
                self.gidax.autoscale(False)
                qxymin=self.pilGIDXAxs_Q[0][np.where(self.pilGIDXAxs_Q[0]>0)[0][0]]
                qzmin=self.pilGIDYAxs_Q[:,0][np.where(self.pilGIDYAxs_Q[:,0]>0)[0][0]]
                qmin=np.sqrt(qxymin**2+qzmin**2)  # get min q in the 2d plot
                qmax=np.sqrt(self.pilGIDXAxs_Q[0][-1]**2+self.pilGIDYAxs_Q[-1][0]**2) # get max q in the 2d plot
                x=Symbol('x')
                anmin=solve((x-qmin*self.formfactor_r)*(x**2-8/3)-2,x)
                anmax=solve((x-qmax*self.formfactor_r)*(x**2-8/3)-2,x)                
                nmin=anmin[-1]/np.pi-0.5     #get min n based on radius and min q
                nmax=anmax[-1]/np.pi-0.5     #get min n based on radius and min q
               # print int(nmin+1), int(nmax+1),anmin,anmax
                if int(nmin)==int(nmax):
                    self.messageBox('Warning: not form factor minimun in the plot range.')
                for i in range(int(nmin+1),int(nmax+1)):                 # get n range
                    totq=((i+0.5)*np.pi-2/((i+0.5)**2*np.pi**2-8/3))/self.formfactor_r #get total q for n
                    plotqxy=np.linspace(self.pilGIDXAxs_Q[0][0],totq,100)   # get the x-axis for formfactor for ith order
                    plotqz=np.sqrt(totq**2-plotqxy**2)       # get the y-axis for formfactor for ith order
                    self.gidax.plot(plotqxy, plotqz, 'w-', linewidth=1.5)
            if self.formfactor_rxy>0 and self.formfactor_rz>0: #for ellipsoid ff
                self.gidax.autoscale(False)
                qxymin=self.pilGIDXAxs_Q[0][np.where(self.pilGIDXAxs_Q[0]>0)[0][0]]
                qzmin=self.pilGIDYAxs_Q[:,0][np.where(self.pilGIDYAxs_Q[:,0]>0)[0][0]]
                qmin=np.sqrt(qxymin**2*self.formfactor_rxy**2+qzmin**2*self.formfactor_rz**2) #get min q for ellipsoid
                qmax=np.sqrt(self.pilGIDXAxs_Q[0][-1]**2*self.formfactor_rxy**2+self.pilGIDYAxs_Q[-1][0]**2*self.formfactor_rz**2) #get max q for ellipsoid
                x=Symbol('x')
                anmin=solve((x-qmin)*(x**2-8/3)-2,x)
                anmax=solve((x-qmax)*(x**2-8/3)-2,x)                
                nmin=anmin[-1]/np.pi-0.5     #get min n based on radius and min q
                nmax=anmax[-1]/np.pi-0.5     #get min n based on radius and min q
                #print int(nmin+1),int(nmax+1),anmin,anmax
                if int(nmin)==int(nmax):
                    self.messageBox('Warning: not form factor minimun in the plot range.')
                for i in range(int(nmin+1),int(nmax+1)):                 # get n range
                    totq=(i+0.5)*np.pi-2/((i+0.5)**2*np.pi**2-8/3) #get total qr for n
                    plotqxy=np.linspace(self.pilGIDXAxs_Q[0][0],totq/self.formfactor_rxy,100)   # get the x-axis for formfactor for ith order
                    plotqz=np.nan_to_num(np.sqrt((totq**2-plotqxy**2*self.formfactor_rxy**2)/self.formfactor_rz**2))       # get the y-axis for formfactor for ith order
                    self.gidax.plot(plotqxy, plotqz, 'w-', linewidth=1.5)
             #get the h-cut data and label
            self.gid_h=np.vstack((self.pilGIDXAxs_Q[0,:],np.sum(self.pilGIDData_Q,axis=0),np.sqrt(np.sum(self.pilGIDDataErr_Q**2,axis=0)))).transpose()
            self.gidax_h.set_xlabel(r'$Q_{xy}$'+' '+r'$[\AA^{-1}]$')
            self.gidax_h.set_title('Integrated over '+r'$Q_z$')
            #get the v-cut data
            self.gid_v=np.vstack((self.pilGIDYAxs_Q[:,0],np.sum(self.pilGIDData_Q,axis=1),np.sqrt(np.sum(self.pilGIDDataErr_Q**2,axis=1)))).transpose()
            self.gidax_v.set_ylabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
            self.gidax_v.set_title('Integrated over '+r'$Q_{xy}$')
            #get the v-slice-cut data
            vc_ran=np.linspace(np.where(self.pilGIDYAxs_Q[:,0]>0)[0][0],np.where(self.pilGIDYAxs_Q[:,0]>0)[0][-1],5)
            self.gid_vc={}
            for i in range(4):
                self.gid_vc[i]=np.vstack((self.pilGIDXAxs_Q[0,:],np.sum(self.pilGIDData_Q[int(vc_ran[i]):int(vc_ran[i+1]),:],axis=0),np.sqrt(np.sum(self.pilGIDDataErr_Q[int(vc_ran[i]):int(vc_ran[i+1]),:]**2,axis=0)))).transpose()
            self.gidax_vc.set_title(r'$Q_z$'+'-binned')
            self.gidax_vc.set_xlabel(r'$Q_{xy}$ '+r'$[\AA^{-1}]$')
        
        self.gidax.set_aspect(aspect)
        self.gidax.set_xlabel(self.pilGIDxlabel)
        self.gidax.set_ylabel(self.pilGIDylabel)
        self.gidax.set_title('Two-dimensional Data')
        self.gidax.format_coord=self.format_pil_coord
        vc_max=max(np.max(self.gid_vc[0][:,1]),np.max(self.gid_vc[1][:,1]),np.max(self.gid_vc[2][:,1]),np.max(self.gid_vc[3][:,1]))
        #self.gidax_v.errorbar(self.gid_v[:,1]/np.max(self.gid_v[:,1]),self.gid_v[:,0],xerr=self.gid_v[:,2]/np.max(self.gid_v[:,1]),fmt='b-')
        self.gidax_v.errorbar(self.gid_v[:,1],self.gid_v[:,0],xerr=self.gid_v[:,2],fmt='b-')
        #self.gidax_h.errorbar(self.gid_h[:,0],self.gid_h[:,1]/np.max(self.gid_h[:,1]),self.gid_h[:,2]/np.max(self.gid_h[:,1]),fmt='b-')
        self.gidax_h.errorbar(self.gid_h[:,0],self.gid_h[:,1],self.gid_h[:,2],fmt='b-')
        for i in range(4):
            #self.gidax_vc.errorbar(self.gid_vc[i][:,0], (self.gid_vc[i][:,1]+vc_max*i)/vc_max/4, self.gid_vc[i][:,2]/vc_max/4, fmt='b-')
            self.gidax_vc.errorbar(self.gid_vc[i][:,0], self.gid_vc[i][:,1]+vc_max*i, self.gid_vc[i][:,2], fmt='b-')
        for i in range(3):
            #self.gidax_vc.plot(self.gid_vc[i][:,0],(i+1)*0.25*np.ones_like(self.gid_vc[i]),'r--')
            self.gidax_vc.plot(self.gid_vc[i][:,0],(i+1)*vc_max*np.ones_like(self.gid_vc[i]),'r--')
        self.gidax_h.set_ylabel('Intensity (a.u.)')
        self.gidax_v.set_xlabel('Intensity (a.u.)')
        self.gidax_vc.set_ylabel('Intensity (a.u.)')
        self.gidax_v.yaxis.tick_right()
        self.gidax_v.yaxis.set_label_position("right")
        self.gidax_h.set_xlim(self.gid_h[:,0][0],self.gid_h[:,0][-1])
        self.gidax_vc.set_xlim(self.gid_h[:,0][0],self.gid_h[:,0][-1])
        self.gidax_v.set_ylim(self.gid_v[:,0][0],self.gid_v[:,0][-1])
        if self.ui.pilLogIntCheckBox.checkState()!=0:
            self.gidax_h.set_yscale('log')
            self.gidax_v.set_xscale('log')
        else:
            self.gidax_h.set_yscale('linear')
            self.gidax_v.set_xscale('linear')
        #    self.gidax_v.xaxis.set_ticks(np.linspace(0,1,3))
        #    self.gidax_h.yaxis.set_ticks(np.linspace(0,1,3))
        #self.gidax_vc.xaxis.set_ticks(np.linspace(float(format(self.gid_vc[0][:,0][0],'.2g')),float(format(self.gid_vc[0][:,0][-1],'.2g')),3))
        #self.gidax_v.locator_params(axis = 'x', nbins = 4)
        self.gidax_vc.locator_params(axis = 'x', nbins = 4)
        self.ui.pilMplWidget.canvas.fig.suptitle('File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1])
        temp2d=make_axes_locatable(self.gidax)
        cax2d=temp2d.append_axes("right", size="5%",pad=0.05,)
        #cax2d = self.ui.pilMplWidget.canvas.fig.add_axes([0.565, 0.475, 0.008, 0.425]) 
        self.ui.pilMplWidget.canvas.fig.colorbar(self.gid_p,cax=cax2d,format="%.1e")
        
        self.ui.pilMplWidget.canvas.fig.tight_layout()
        self.ui.pilMplWidget.canvas.fig.subplots_adjust(top=0.9)
        self.ui.pilMplWidget.canvas.draw()
        self.command=str(self.ui.gidComboBox.currentText())+', scans=['+str([item for item in np.sort(self.selectedScanNums)])[1:-1]+'], sli=['+self.ui.pilHSlitLineEdit.text()+'], /'+self.ui.pilAxesComboBox.currentText()
        self.ui.commandLineEdit.setText(self.command)
        
    def saveGIDData(self):
#        print self.pilGIDData_Q
#        print self.pilGIDXAxs_Q[0]
#        print self.pilGIDYAxs_Q[:,0]
#        print len(self.pilGIDData_Q), len(self.pilGIDData_Q[0]),  len(self.pilGIDXAxs_Q[0]), len(self.pilGIDYAxs_Q[:,0])
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save gioxs data', directory=self.directory))
        data2d=[]
        if str(self.ui.pilAxesComboBox.currentText())=='Q':
            fname_h=self.saveFileName+'_qxy_gid.txt'
            fname_v=self.saveFileName+'_qz_gid.txt'
            fname_2d=self.saveFileName+'_q2d_gid.txt'
            for i in range(len(self.pilGIDData_Q)):
                for j in range(len(self.pilGIDData_Q[0])):
                    data2d.append([self.pilGIDXAxs_Q[i][j],self.pilGIDYAxs_Q[i][j],self.pilGIDData_Q[i][j],self.pilGIDDataErr_Q[i][j]])
        elif str(self.ui.pilAxesComboBox.currentText())=='Angles':
            fname_h=self.saveFileName+'_dth_gid.txt'
            fname_v=self.saveFileName+'_beta_gid.txt'
            fname_2d=self.saveFileName+'_agnle2d_gid.txt'
            for i in range(len(self.pilGIDData)):
                for j in range(len(self.pilGIDData[0])):
                    data2d.append([self.pilGIDXAxs[i][j],self.pilGIDYAxs[i][j],self.pilGIDData[i][j],self.pilGIDDataErr[i][j]])
#        print np.array(data2d)
#        print len(np.array(data2d)), len(np.array(data2d[0]))
        np.savetxt(fname_h,self.gid_h,fmt='%.4e\t%.4e\t%.4e')
        np.savetxt(fname_v,self.gid_v,fmt='%.4e\t%.4e\t%.4e')
        np.savetxt(fname_2d,np.array(data2d),fmt='%.4e\t%.4e\t%.4e\t%.4e')
        
        
    def imagehjoin(self, data1, data2): #join two images horizontally 
        data1x=np.array([data1[i][0][0] for i in range(len(data1))])  #get dth/x range for data1
        data1z=np.array([np.sum(data1[i][:,2]) for i in range(len(data1))]) #get integrated (along beta) for each dth in data1
        data1zerr=np.array([np.sqrt(np.sum(data1[i][:,3]**2)) for i in range(len(data1))]) # get error bars
        data2x=np.array([data2[i][0][0] for i in range(len(data2))])  #get dth range for data2
        data2z=np.array([np.sum(data2[i][:,2]) for i in range(len(data2))])  #get integrated (along beta) for each dth in data2
        data2zerr=np.array([np.sqrt(np.sum(data2[i][:,3]**2)) for i in range(len(data2))]) #get error bars
        data1step=np.abs(data1x[0]-data1x[1])
        data2step=np.abs(data2x[0]-data2x[1])
        if data2x[0]>=data1x[-1]:
            self.messageBox('Warning: no overlap area between two scans!')
            return np.vstack((data1,data2))
        else:  #patching
            overlap2=np.where(data2x<=data1x[-1]) #create an array for where data2 having smaller dth than the point in data1
            overlap1=np.where(np.logical_and(data1x>=data2x[0],data1x<=data2x[overlap2[0][-1]])) #create an array for where data1 having larger dth than the point in data2
            f=interp1d(data2x[overlap2],data2z[overlap2],kind='linear')
            ferr=interp1d(data2x[overlap2],data2zerr[overlap2],kind='linear')
            nonzero=np.nonzero(f(data1x[overlap1])) #get rid of nonzero element in f(data1x[overlap1])
            overlap1z=data1z[overlap1][nonzero]
            overlap1zerr=data1zerr[overlap1][nonzero]
            overlap2z=f(data1x[overlap1])[nonzero]
            overlap2zerr=ferr(data1x[overlap1])[nonzero]
           # print overlap1y, overlap1yerr
            #print overlap2y, overlap2yerr
            x=data1x[overlap1][nonzero]   #get x, y, yerr for the ratio between two images in the overlap range
            y=overlap1z/overlap2z
            yerr=y*np.sqrt(((overlap1zerr/overlap1z)**2+(overlap2zerr/overlap2z)**2))
            pfit,pcov=curve_fit(self.polynomial,x,y,p0=[1],sigma=yerr) #fit the data with the scale factor 
            fac=pfit[0]
            facerr=pcov[0][0]**0.5
           # print fac, facerr
            for i in range(len(data2)):  # scale the data2 with the fac and facerr
                data2[i][:,2]=data2[i][:,2]*fac
                data2[i][:,3]=np.sqrt(data2[i][:,3]**2*fac**2+data2[i][:,2]**2*facerr**2/fac**2)
           # print overlap1, overlap2
           # print data1x[overlap1], data2x[overlap2]
            data=[]
            for i in range(overlap1[0][0]):  #part 1 for patching data 
                data.append(data1[i])
            if data1step<=data2step:  # if data1 has finer or equal step size; part 2 for patching
                data2zoverlap=np.array([np.hstack(data2[i][:,2].T) for i in overlap2[0]]).T
                data2zerroverlap=np.array([np.hstack(data2[i][:,3].T) for i in overlap2[0]]).T
               # print data2zoverlap,  data2zerroverlap, data2x[overlap2]
                f2d=interp2d(data2x[overlap2],data2[0][:,1],data2zoverlap, kind='linear')  #interpolate the data2 in the overlap2
                f2derr=interp2d(data2x[overlap2],data2[0][:,1],data2zerroverlap, kind='linear')
                data3=f2d(data1x[overlap1],data2[0][:,1])  #get data2 data at data1 x-axis from the interpolation 
                data3err=f2derr(data1x[overlap1],data2[0][:,1]) 
                for i in range(overlap1[0][0],overlap1[0][-1]+1):
                    #use avreage 
                    data.append(np.vstack((data1[i][:,0],data1[i][:,1],(data1[i][:,2]+data3[:,i-overlap1[0][0]])/2.0,np.sqrt(data1[i][:,3]**2+data3err[:,i-overlap1[0][0]]**2)/2.0)).transpose())
                    #use minimum error combination for the overlap area                    
                    #data.append(np.vstack((data1[i][:,0],data1[i][:,1],np.nan_to_num((data1[i][:,2]*data3err[:,i-overlap1[0][0]]**2+data1[i][:,3]**2*data3[:,i-overlap1[0][0]])/(data1[i][:,3]**2+data3err[:,i-overlap1[0][0]]**2)),np.nan_to_num(data1[i][:,3]*data3err[:,i-overlap1[0][0]]/np.sqrt(data1[i][:,3]**2+data3err[:,i-overlap1[0][0]]**2)))).transpose())
                    #use weighted average
                    #data.append(np.vstack((data1[i][:,0],data1[i][:,1],np.nan_to_num((data1[i][:,2]**2*data3err[:,i-overlap1[0][0]]**2+data1[i][:,3]**2*data3[:,i-overlap1[0][0]]**2)/(data1[i][:,3]**2*data3[:,i-overlap1[0][0]]+data3err[:,i-overlap1[0][0]]**2*data1[i][:,2])),np.nan_to_num(np.sqrt(data1[i][:,3]**2*data1[i][:,2]**2*data3err[:,i-overlap1[0][0]]**4+data1[i][:,3]**4*data3[:,i-overlap1[0][0]]**2*data3err[:,i-overlap1[0][0]]**2)/(data1[i][:,3]**2*data3[:,i-overlap1[0][0]]+data3err[:,i-overlap1[0][0]]**2*data1[i][:,2])))).transpose())
               # print data1[overlap1[0][-1]], data2[overlap2[0][-1]], data[-1]
            else: # if data2 has finer step size
                data1zoverlap=np.array([np.hstack(data1[i][:,2].T) for i in overlap1[0]]).T
                data1zerroverlap=np.array([np.hstack(data1[i][:,3].T) for i in overlap1[0]]).T
               # print data2zoverlap,  data2zerroverlap, data2x[overlap2]
                f2d=interp2d(data1x[overlap1],data1[0][:,1],data1zoverlap, kind='linear')  #interpolate the data1 in the overlap1
                f2derr=interp2d(data1x[overlap1],data1[0][:,1],data1zerroverlap, kind='linear')
                data3=f2d(data2x[overlap2],data1[0][:,1])  #get data1 data at data2 x-axis from the interpolation 
                data3err=f2derr(data2x[overlap2],data1[0][:,1]) 
                for i in range(overlap2[0][0],overlap2[0][-1]+1):
                    data.append(np.vstack((data2[i][:,0],data2[i][:,1],(data2[i][:,2]+data3[:,i-overlap2[0][0]])/2.0,np.sqrt(data2[i][:,3]**2+data3err[:,i-overlap2[0][0]]**2)/2.0)).transpose())
                    #data.append(np.vstack((data2[i][:,0],data2[i][:,1],np.nan_to_num((data2[i][:,2]*data3err[:,i-overlap2[0][0]]**2+data2[i][:,3]**2*data3[:,i-overlap2[0][0]])/(data2[i][:,3]**2+data3err[:,i-overlap2[0][0]]**2)),np.nan_to_num(data2[i][:,3]*data3err[:,i-overlap2[0][0]]/np.sqrt(data2[i][:,3]**2+data3err[:,i-overlap2[0][0]]**2)))).transpose())
                    #data.append(np.vstack((data2[i][:,0],data2[i][:,1],np.nan_to_num((data2[i][:,2]**2*data3err[:,i-overlap2[0][0]]**2+data2[i][:,3]**2*data3[:,i-overlap2[0][0]]**2)/(data2[i][:,3]**2*data3[:,i-overlap2[0][0]]+data3err[:,i-overlap2[0][0]]**2*data2[i][:,2])),np.nan_to_num(np.sqrt(data2[i][:,3]**2*data2[i][:,2]**2*data3err[:,i-overlap2[0][0]]**4+data2[i][:,3]**4*data3[:,i-overlap2[0][0]]**2*data3err[:,i-overlap2[0][0]]**2)/(data2[i][:,3]**2*data3[:,i-overlap2[0][0]]+data3err[:,i-overlap2[0][0]]**2*data2[i][:,2])))).transpose())
               # print data1[overlap1[0][0]][-10:], data1[overlap1[0][1]][-10:]
               # print data3[:,0][-10:], data3[:,1][-10:], data3err[:,1][-10:], data3[:,2][-10:]
              #  print data2[overlap2[0][1]][-10:], data[-20][-10:], data[-19][-10:]
            if data1x[-1]>=data2x[-1]: #part 3 for patching
                for i in range(overlap1[0][-1]+1,len(data1)):
                    data.append(data1[i])
            else: 
                for i in range(overlap2[0][-1]+1,len(data2)):
                    data.append(data2[i])
            #print data1[overlap1[0][-1]], data2[-1], data[-1]   
            #print data
            return data
    
    def imagevjoin(self, data1, data2):
        data1y=data1[0][:,1]  #get beta/y range for data1
        data2y=data2[0][:,1]  #get beta/y range for data2
        if data1y[-1]<data2y[0]:
            self.messageBox('Warning: no overlap area between two scans!')
            return np.vstack((data1,data2))
        else:
            data1x=np.array([data1[i][0][0] for i in range(len(data1))])  #get dth/x range for data1
            data2x=np.array([data2[i][0][0] for i in range(len(data2))]) #get dth/x range for data2
            data2z=np.array([np.hstack(data2[i][:,2]) for i in range(len(data2))]).T  #get z for data2
            data2zerr=np.array([np.hstack(data2[i][:,3]) for i in range(len(data2))]).T #get zerr for data2
            fdata2=interp2d(data2x,data2y,data2z, kind='linear')  #interpolate data2 to match data1 x-axis
            fdata2err=interp2d(data2x,data2y,data2zerr,kind='linear')
            data2znew=fdata2(data1x,data2y)
            data2zerrnew=fdata2err(data1x,data2y)
            data2new=[np.vstack((data1x[0]*np.ones_like(data2y),data2y,data2znew[:,0],data2zerrnew[:,0])).transpose()]
            for i in range(1,len(data1x)):
                data2new.append(np.vstack((data1x[i]*np.ones_like(data2y),data2y,data2znew[:,i],data2zerrnew[:,i])).transpose())
            data2=data2new
            overlap1=np.where(data1y>=data2y[0]) #create an array for where data1 having larger beta than the point in data2
            overlap2=np.where(data2y<=data1y[-1]) #create an array for where data2 having smaller beta than the point in data1
            #print overlap1, overlap2, len(overlap1[0]), len(overlap2[0])
            data1zint=np.array([np.sum([data1[i][j,2] for i in range(len(data1))]) for j in overlap1[0]]) #get integrated for each beat for data1 in overlap1
            data1zinterr=np.array([np.sqrt(np.sum([data1[i][j,3]**2 for i in range(len(data1))])) for j in overlap1[0]]) #get error bars
            data2zint=np.array([np.sum([data2[i][j,2] for i in range(len(data2))]) for j in overlap2[0]]) #get integrated for each beat for data2 in overlap2
            data2zinterr=np.array([np.sqrt(np.sum([data2[i][j,3]**2 for i in range(len(data2))])) for j in overlap2[0]]) #get error bars
           # print len(data2),data1zint,data2zint
           # print np.sum(data1zint)*len(data2zint)/(np.sum(data2zint)*len(data1zint))
            f=interp1d(data2y[overlap2],data2zint, kind='cubic')
            ferr=interp1d(data2y[overlap2],data2zinterr,kind='cubic')
            #print data1y[overlap1][:-1]
            nonzero=np.nonzero(f(data1y[overlap1][:-1]))
            overlap1z=data1zint[:-1][nonzero]
            overlap1zerr=data1zinterr[:-1][nonzero]
            overlap2z=f(data1y[overlap1][:-1])[nonzero]
            overlap2zerr=ferr(data1y[overlap1][:-1])[nonzero]
           # print len(overlap1z), len(overlap2z)
            x=data1y[overlap1][:-1][nonzero]   #get x, y, yerr for the ratio between two images in the overlap range
            y=overlap1z/overlap2z
            yerr=y*np.sqrt(((overlap1zerr/overlap1z)**2+(overlap2zerr/overlap2z)**2))
            pfit,pcov=curve_fit(self.polynomial,x,y,p0=[1],sigma=yerr) #fit the data with the scale factor 
            fac=pfit[0]
            facerr=pcov[0][0]**0.5
            #print np.sum(y)/len(y), fac, facerr
            for i in range(len(data2)):  # scale the data2 with the fac and facerr
                data2[i][:,2]=data2[i][:,2]*fac
                data2[i][:,3]=np.sqrt(data2[i][:,3]**2*fac**2+data2[i][:,2]**2*facerr**2/fac**2)
            #print data2[0][:,2][:overlap2[0][-1]+1]
            data2zoverlap=np.array([np.hstack(data2[i][:,2][:overlap2[0][-1]+1].T) for i in range(len(data2))]).T
            data2zerroverlap=np.array([np.hstack(data2[i][:,3][:overlap2[0][-1]+1].T) for i in range(len(data2))]).T
            f2d=interp2d(data1x,data2y[overlap2],data2zoverlap, kind='linear')  #interpolate the data2 in the overlap2
            f2derr=interp2d(data1x,data2y[overlap2],data2zerroverlap, kind='linear')
            data3=f2d(data1x,data1y[overlap1][:-1])  #get data2 data at data1 y-axis from the interpolation 
            data3err=f2derr(data1x,data1y[overlap1][:-1]) #get data2 error bars at data1 y-axis from the interpolation 
            #print  data3[:,:10]
            #print len(data3[:,:10][0]), data3[:,:10][0]
            #print data1[0][0:overlap1[0][0]]
            #print data2zoverlap[:,:10]
            #print data1y[overlap1][-10:], data2y[overlap2][-10:]
            #print data2[0][0,2],data2[0][0,3],data2[-1][0,2],data2[-1][0,3]
            #print data3[-1], data3err[-1]
           # data=np.vstack(data1[i][0:overlap1[0][0]]
            data1z=np.array([np.hstack(data1[i][:,2]) for i in range(len(data1))])  
            data1zerr=np.array([np.hstack(data1[i][:,3]) for i in range(len(data1))])
            data2z=np.array([np.hstack(data2[i][:,2]) for i in range(len(data2))])  
            data2zerr=np.array([np.hstack(data2[i][:,3]) for i in range(len(data2))])
            data=[np.vstack((data1y[0]*np.ones_like(data1x),data1x,data1z[:,0],data1zerr[:,0])).transpose()]
            for i  in range(1,overlap1[0][0]):
                data.append(np.vstack((data1y[i]*np.ones_like(data1x),data1x,data1z[:,i],data1zerr[:,i])).transpose())
            for i in range(overlap1[0][0],overlap1[0][-1]):
                data.append(np.vstack((data1y[i]*np.ones_like(data1x),data1x,(data1z[:,i]+data3.T[:,i-overlap1[0][0]])/2.0,np.sqrt(data1zerr[:,i]**2+data3err.T[:,i-overlap1[0][0]]**2)/2.0)).transpose())
            #print data1z[:,-1],data1z[:,-2],data1zerr[:,-2]
            #print data[-1]
            for i in range(overlap2[0][-1],len(data2[0])):
                data.append(np.vstack((data2y[i]*np.ones_like(data1x),data1x,data2z[:,i],data2zerr[:,i])).transpose())
            dataz=np.array([np.hstack(data[i][:,2]) for i in range(len(data))]) 
            datazerr=np.array([np.hstack(data[i][:,3]) for i in range(len(data))]) 
            datay=np.array([data[i][0][0] for i in range(len(data))])
            datanew=[np.vstack((data1x[0]*np.ones_like(datay),datay,dataz[:,0],datazerr[:,0])).transpose()]
            for i in range(1,len(data1x)):
                datanew.append(np.vstack((data1x[i]*np.ones_like(datay),datay,dataz[:,i],datazerr[:,i])).transpose())
            #print len(datanew)
        return datanew
    
    
    def sortedScans(self,data):  #reture the order of scans from low dth to higt dth, qlso works for sorting qz in the ref 
        if len(data)<2:
            return data.keys()
        else:
            x={}
            sorted_keys=[]
            for keys in data.keys():
                try:
                    x[keys]=sorted([data[keys][0][0][0],data[keys][-1][0][0]]) #get dth range for each key, i.e., scan number
                except:
                    x[keys]=sorted([data[keys][0][0],data[keys][-1][0]])  #get qz range for each key
            sorted_x=sorted(x.items(), key=lambda t:t[1])  #sort the dth range
            for i in range(len(sorted_x)):
                sorted_keys.append(sorted_x[i][0]) #make an array for keys (scan number) based on the dth range
            return sorted_keys
        
    def polynomial(self,x,*p):    #define polynomial function
        return sum([p[i]*x**i for i in range(len(p))])
        
    def pilGIDData(self):
        #print 'I am missing'
        if self.pilGIDshow==0:
            self.disconnect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
            self.ui.pilAxesComboBox.setCurrentIndex(1)
            self.connect(self.ui.pilAxesComboBox, SIGNAL('currentIndexChanged(int)'), self.update2dPlots)
        self.pilGIDshow=1
        aspect=str(self.ui.pilAspectComboBox.currentText())
        vmax=float(self.ui.pilMaxLineEdit.text())
        vmin=float(self.ui.pilMinLineEdit.text())
        self.ui.pilMplWidget.canvas.fig.clf()
        self.pilxleft=self.xcenter-(float(self.ui.pilHSlitLineEdit.text())-1)/2
        self.pilxright=self.xcenter+(float(self.ui.pilHSlitLineEdit.text())-1)/2
        self.pilhintDatadic={}
        cmap=str(self.ui.pilCMapComboBox.currentText())
        self.absfac=float(self.ui.AbsFacLineEdit.text())
        #sortedFrameNums=np.array(self.selectedPilFramesNums)[np.argsort(np.array(self.pil_Dth))]
        j=0   #for xcenter and ycenter, which are already defined in self.selectedPilFramesNums
        for i in self.selectedPilFramesNums:
            self.pilatus.plotHint({i:self.pilData[i]},{i:self.pilErrorData[i]},absfac=self.absfac,absnum=[self.pil_AbsNum[i]],cen=[self.pilX[i],self.pilY[i]], hroi=[self.pilxleft[j]-1, self.pilxright[j]-1], vroi=[5,486], ax_type='Angles', wavelength=[self.pil_Wavelength[i]],s2d_dist=[self.pil_Dist[i]], sh=[self.pil_Sh[i]],alpha=self.pilAlpha[i], truealpha=self.pilTrueAlpha[i], mon=[self.pilMonc[i]])
            #ang=((self.pilatus.hintData[:,0]-self.ycenter[j])*0.172/self.distance[0]+np.arcsin(2.0*np.sin(self.pilAlpha[i])-np.sin(self.pilTrueAlpha[i])))*180/np.pi #convert pixle to beta
            #print self.pilatus.hintData[0][0]             
            j=j+1
            ckey=float(format(self.pil_Dth[i]*180/np.pi,'.2f'))
            if ckey in self.pilhintDatadic:
 #               self.pilhintDatadic[ckey].append(np.vstack((ang,self.pilatus.hintData[:,1],self.pilatus.hintData[:,2])).transpose())  #same dth but may be different starting qz 
                self.pilhintDatadic[ckey].append(np.vstack((self.pilatus.hintData[:,0],self.pilatus.hintData[:,1],self.pilatus.hintData[:,2])).transpose())
            else:
 #               self.pilhintDatadic[ckey]=[np.vstack((ang,self.pilatus.hintData[:,1],self.pilatus.hintData[:,2])).transpose()]
                self.pilhintDatadic[ckey]=[np.vstack((self.pilatus.hintData[:,0],self.pilatus.hintData[:,1],self.pilatus.hintData[:,2])).transpose()]
        lengths=np.array([len(self.pilhintDatadic[key]) for key in sorted(self.pilhintDatadic.keys())])
#        if np.any(lengths!=lengths[0]):
#            self.messageBox('Error:: Couldnot patch the frames perhaps some frames are missing for doing #the patch.')
#        else:
        self.pilhintData=[]
        self.pilhintErrorData=[]
        pilXAxis=[]
        for keys in sorted(self.pilhintDatadic.keys()):
            data=sorted(self.pilhintDatadic[keys],key=lambda x:x[0,0])  #each data for same qxy 
            pdata={}
            for i in range(len(data)): #creating dictonary of same starting qz's with averaging the same starting qz arrays
                c1keys=float(format(data[i][0,0],'.3f'))
                if c1keys in pdata: #same starting qz
                    pdata[c1keys]=(pdata[c1keys]+data[i])/2.0   #averaging 
                else: #different starting qz
                    pdata[c1keys]=data[i]
            keys1=sorted(pdata.keys())  #starting qz for different images with same qxy (i.e., dth)
            pdata2=pdata[keys1[0]]  # data of low starting qz for this dth 
            for i in range(1,len(keys1)):
                overlap1=np.where(pdata2[:,0]>=pdata[keys1[i]][0,0]) #create a list of pixles in lower qz data, which has larger qz than the first point of higher qz data 
                overlap2=np.where(pdata[keys1[i]][:,0]<=pdata2[-1,0]) #create a list of pixles in higher qz data, which has smaller qz than the last point of smaller qz data 
                f=interp1d(pdata[keys1[i]][overlap2,0][0],pdata[keys1[i]][overlap2,1][0]) #create the interpreation function based on the higher qz data in the overlap area
                fac=np.sum(pdata2[overlap1[0][:-1],1])/np.sum(f(pdata2[overlap1[0][:-1],0])) #calculate the scale factor
                pdata[keys1[i]][:,1]=fac*pdata[keys1[i]][:,1]
                pdata[keys1[i]][:,2]=fac*pdata[keys1[i]][:,2]
                pdata2=np.vstack((pdata2[:overlap1[0][0],:],pdata[keys1[i]]))
            self.pilhintData.append(list(pdata2[:,1]))  #sige intensity for pilhintData for this dth
            self.pilhintErrorData.append(list(pdata2[:,2]))  #sign errorbar of intensity for pilhintData for this dth
            pilXAxis.append(float(keys))  #sign x-axis for this dth
        pilYAxis=list(pdata2[:,0])  #sign y-axis for all dth
        self.pilhintData=np.array(self.pilhintData).transpose()
        self.pilhintErrorData=np.array(self.pilhintErrorData).transpose()
        self.pilXAxs, self.pilYAxs=np.meshgrid(pilXAxis,pilYAxis)
        self.pil_xlabel=r'$2\theta$'+' '+r'$[\deg]$'
        self.pil_ylabel=r'$\beta$'+' '+r'$[\deg]$'
        vmax=np.max(self.pilhintData)

#            if str(self.ui.pilAxesComboBox.currentText())=='Angles':
#                self.pilYAxs=((self.pilYAxs-self.ycenter[0])*0.172/self.distance[0]+self.alpha[0])*180/np.pi
#                self.pil_ylabel= 'Theta [Degrees]'
        if str(self.ui.pilAxesComboBox.currentText())=='Q':
            self.pilYAxs==2.0*np.pi*(np.sin(self.pilYAxs*np.pi/180.0)+np.sin(self.pilTrueAlpha[0]))/self.wavelength[0] #this coversion is approximation. 
            self.pilXAxs=4.0*np.pi*np.sin(self.pilXAxs*np.pi/2.0/180.0)/self.wavelength[0]  #this coversion is approximation.
            self.pil_xlabel=r'$Q_{xy}$'+' '+r'$[\AA^{-1}]$'
            self.pil_ylabel=r'$Q_z$'+' '+r'$[\AA^{-1}]$'
        self.gidax=self.ui.pilMplWidget.canvas.fig.add_subplot(1,1,1)
        if self.ui.pilLogIntCheckBox.checkState()!=0:
            self.gid_p=self.gidax.pcolor(self.pilXAxs, self.pilYAxs, np.log10(self.pilhintData),cmap=pl.cm.get_cmap(cmap),vmin=np.log10(vmin), vmax=np.log10(vmax))
        else:
            fint=interp1d(self.pilXAxs[0], self.pilhintData, kind='cubic')
            #print self.pilXAxs[0], self.pilhintData, len(self.pilhintData), len(self.pilhintData[0])
            newpilXAxis=np.linspace(self.pilXAxs[0][0],self.pilXAxs[0][-1],len(self.pilXAxs[0])*10)
            inter_pilhitData=fint(newpilXAxis)
            self.extent=[self.pilXAxs[0][0],self.pilXAxs[0][-1],self.pilYAxs[0][0],self.pilYAxs[-1][-1]]
            self.gid_p=self.gidax.imshow(inter_pilhitData,interpolation='bicubic',extent=self.extent,vmax=vmax,vmin=vmin,cmap=cmap,aspect=aspect,origin='lower')            
            #self.gid_p=self.gidax.pcolor(self.pilXAxs, self.pilYAxs, self.pilhintData,cmap=pl.cm.get_cmap(cmap),vmin=vmin, vmax=vmax)           
        #self.gidax.set_xlim((np.min(self.pilXAxs),np.max(self.pilXAxs)))
        #self.gidax.set_ylim((np.min(self.pilYAxs),np.max(self.pilYAxs)))
        self.gidax.set_aspect(aspect)
        self.gidax.set_xlabel(self.pil_xlabel)
        self.gidax.set_ylabel(self.pil_ylabel)
        self.ui.pilMplWidget.canvas.fig.colorbar(self.gid_p)
        self.ui.pilMplWidget.canvas.draw()
       
    

    def format_coord(self,x,y):  # cursor interactive 
        numrows, numcols = self.Zdata.shape
        col = int(np.abs((x-self.extent[0]))*numrows/np.abs((self.extent[1]-self.extent[0])))
        row = int(np.abs((y-self.extent[3]))*numcols/np.abs((self.extent[3]-self.extent[2])))
        if col>=0 and col<numcols and row>=0 and row<numrows:
            z = self.Zdata[row,col]
            return self.xyzformat%(x, y, z)
        else:
            return 'x=%.2f, y=%.2f'%(x, y)
    
    def format_pil_coord(self,x,y):  # cursor interactive 
        numrows, numcols = self.Zdata.shape
        col = int(np.abs((x-self.extent[0]))*numcols/np.abs((self.extent[1]-self.extent[0])))
        row = int(np.abs((y-self.extent[2]))*numrows/np.abs((self.extent[3]-self.extent[2])))
        if col>=0 and col<numcols and row>=0 and row<numrows:
            z = self.Zdata[row,col]
            return self.xyzformat%(x, y, z)
        else:
            return 'x=%.2f, y=%.2f'%(x, y)
    
    def format_coord_ref(self,x,y):
        numrows, numcols = self.refZdata.shape
        col = int(np.abs((x-self.extent[0]))*numrows/np.abs((self.extent[1]-self.extent[0])))
        row = int(np.abs((y-self.extent[3]))*numcols/np.abs((self.extent[3]-self.extent[2])))
        if col>=0 and col<numcols and row>=0 and row<numrows:
            z = self.refZdata[row,col]
            return 'x=%d,y=%d,z=%.1e'%(x, y, z)
        else:
            return 'x=%d, y=%d'%(x, y)
        
    def refPlotWin(self):
#        if self.ui.gixSumCheckBox.checkState()==0:
#            self.ui.gixSumCheckBox.setCheckState(2)
#            self.messageBox('Warning:: Multiple frames selected, summing over all the frames for Reflectivity!!')
#            self.updateCcdPlotData()
        self.ui.PlotWidget.setCurrentIndex(3)
        self.ui.refCenCheckBox.setCheckState(0)
        self.ui.refComLineEdit.clear()
        if self.refpatchindex==1:
            self.ui.refDataListWidget.clear()
            self.ui.refNormRefListWidget.clear()
            self.refData=[]
        self.refpatchindex=0  #different index for pilatus patch model
        if self.det=='Bruker':
            self.updateRefCcdPlotData()
        elif self.det=='Pilatus':
            self.updateRefPilPlotData()
            
        
       
    def updateRefCcdPlotData(self):
        self.ui.refADDataPlotWidget.canvas.fig.clf()
#        dist=float(self.ui.gixSDDistLineEdit.text())
#        cmap=str(self.ui.gixCMapComboBox.currentText())
        self.disconnect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
        self.ui.refQzListWidget.clear()
        self.scanframe={}
        self.refSortedFrameNums={}
        if self.ui.gixMergeSameQzCheckBox.checkState()!=0:
            j=0
            for key in sorted(self.ccdFrameNumsWithSameQz.keys()):
                self.scanframe[key]={}
                for item in self.ccdFrameNumsWithSameQz[key]:
                    ekey=int(str(self.ui.imageListWidget.item(item).text()).split('\t')[0].split('#')[1])
                    if ekey in self.scanframe[key]:
                        self.scanframe[key][ekey].append(int(str(self.ui.imageListWidget.item(item).text()).split('\t')[1].split('#')[1]))
                    else:
                        self.scanframe[key][ekey]=[int(str(self.ui.imageListWidget.item(item).text()).split('\t')[1].split('#')[1])]
                if all(x==self.ccdFrameAbsWithSameQz[key][0] for x in self.ccdFrameAbsWithSameQz[key]):     #for frames with the same absorber
                    self.ui.refQzListWidget.addItem('Qz= '+str(key)+'\t {S:F}= '+str(self.scanframe[key])+'\t Abs= '+str([self.ccdFrameAbsWithSameQz[key][0]]))
                    self.refSortedFrameNums[j]=self.ccdFrameNumsWithSameQz[key]
                    j=j+1
                else:
                    for i in range(len(self.ccdFrameNumsWithSameQz[key])):
                        sftext=str(self.ui.imageListWidget.item(self.ccdFrameNumsWithSameQz[key][i]).text()).split('\t')
                        self.ui.refQzListWidget.addItem('Qz= '+str(key)+'\t {S:F}= '+str({int(sftext[0].split('#')[1]):[int(sftext[1].split('#')[1])]})+'\t Abs= '+str([self.ccdFrameAbsWithSameQz[key][i]]))
                        self.refSortedFrameNums[j]=[self.ccdFrameNumsWithSameQz[key][i]]
                        j=j+1
        else:
            j=0
            for i in self.selectedCcdFramesNums:
                sftext=str(self.ui.imageListWidget.item(i).text()).split('\t')
                self.ui.refQzListWidget.addItem('Qz= '+str(self.ccdFileQzs[i])+'\t {S:F}= '+str({int(sftext[0].split('#')[1]):[int(sftext[1].split('#')[1])]})+'\t Abs='+str(self.ccd_AbsNum[i]))
                self.refSortedFrameNums[j]=[i]
                j=j+1
        self.connect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
        self.ui.refQzListWidget.setItemSelected(self.ui.refQzListWidget.item(0),True)
        
    def updateRefPilPlotData(self):
        self.ui.refADDataPlotWidget.canvas.fig.clf()
        dist=float(self.ui.pilSDDistLineEdit.text())
        cmap=str(self.ui.pilCMapComboBox.currentText())
        self.disconnect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
        self.ui.refQzListWidget.clear()
        self.scanframe={}
        self.refSortedFrameNums={}
        if self.ui.pilMergeSameQzCheckBox.checkState()!=0:
            j=0
            for key in sorted(self.pilFrameNumsWithSameQz.keys()):
                self.scanframe[key]={}
                for item in self.pilFrameNumsWithSameQz[key]:
                    ekey=int(str(self.ui.imageListWidget.item(item).text()).split('\t')[0].split('#')[1])
                    if ekey in self.scanframe[key]:
                        self.scanframe[key][ekey].append(int(str(self.ui.imageListWidget.item(item).text()).split('\t')[1].split('#')[1]))
                    else:
                        self.scanframe[key][ekey]=[int(str(self.ui.imageListWidget.item(item).text()).split('\t')[1].split('#')[1])]
                if all(x==self.pilFrameAbsWithSameQz[key][0] for x in self.pilFrameAbsWithSameQz[key]):                   
                    self.ui.refQzListWidget.addItem('Qz= '+str(key)+'\t {S:F}= '+str(self.scanframe[key])+'\t Abs= '+str([self.pilFrameAbsWithSameQz[key][0]]))
                    self.refSortedFrameNums[j]=self.pilFrameNumsWithSameQz[key]
                    j=j+1
                else:
                    for i in range(len(self.pilFrameNumsWithSameQz[key])):
                        sftext=str(self.ui.imageListWidget.item(self.pilFrameNumsWithSameQz[key][i]).text()).split('\t')
                        self.ui.refQzListWidget.addItem('Qz= '+str(key)+'\t {S:F}= '+str({int(sftext[0].split('#')[1]):[int(sftext[1].split('#')[1])]})+'\t Abs= '+str([self.pilFrameAbsWithSameQz[key][i]]))
                        self.refSortedFrameNums[j]=[self.pilFrameNumsWithSameQz[key][i]]
                        j=j+1
        else:
            j=0
            for i in self.selectedPilFramesNums:
                sftext=str(self.ui.imageListWidget.item(i).text()).split('\t')
                self.ui.refQzListWidget.addItem('Qz= '+str(self.pilFileQzs[i])+'\t {S:F}= '+str({int(sftext[0].split('#')[1]):[int(sftext[1].split('#')[1])]})+'\t Abs='+str(self.pil_AbsNum[i]))
                self.refSortedFrameNums[j]=[i]
                j=j+1
        self.connect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
        self.ui.refQzListWidget.setItemSelected(self.ui.refQzListWidget.item(0),True)
        
    def refQzListSelectionChanged(self):
        self.ui.refADDataPlotWidget.canvas.fig.clf()
        cmap=str(self.ui.pilCMapComboBox.currentText())
        self.refQzSelectedItems=self.ui.refQzListWidget.selectedItems()        
        self.refQzSelected=float(str(self.refQzSelectedItems[0].text()).split()[1])
        self.absfac=float(self.ui.AbsFacLineEdit.text())       
        self.dir=str(self.ui.refBGDirComboBox.currentText())
        if self.det=='Bruker':
            slitx=int(self.ui.refSlitLineEdit.text().split(',')[0])
            slity=int(self.ui.refSlitLineEdit.text().split(',')[1])
            self.slit=[slitx,slity]
            self.bg=float(self.ui.refBGOffLineEdit.text())
            dist=float(self.ui.gixSDDistLineEdit.text())
            self.selFrameNums=self.refSortedFrameNums[self.ui.refQzListWidget.row(self.refQzSelectedItems[0])]
            #cenx=self.ccdX[self.selFrameNums[0]]-int((dist+self.ccd_Gl2[self.selFrameNums[0]])*(np.tan(self.ccd_Tth[self.selFrameNums[0]]*np.pi/180)-np.tan(2*self.ccd_Phi[self.selFrameNums[0]]*np.pi/180))/0.06)
            cenx=self.ccdX[self.selFrameNums[0]]-int((dist+self.ccd_Gl2[self.selFrameNums[0]])*(np.tan(self.ccd_Tth[self.selFrameNums[0]]*np.pi/180)-np.tan(2*np.arcsin(self.wavelength[self.selFrameNums[0]]/np.pi*1.9236/4)))/0.06)
            ceny=self.ccdY[self.selFrameNums[0]]-int((self.ccdAlpha[self.selFrameNums[0]]*dist+self.ccd_Sh[self.selFrameNums[0]])/0.06)
            self.mon=[self.ccdMonc[i] for i in self.selFrameNums]
            self.absnum=[self.ccd_AbsNum[i] for i in self.selFrameNums]
            self.vmax=float(self.ui.gixMaxLineEdit.text())
            self.vmin=float(self.ui.gixMinLineEdit.text())
            self.cen=[cenx,ceny]
            if self.ui.refCenCheckBox.checkState()==0:
                self.ui.refCenLineEdit.setText(str(self.cen[0])+','+str(self.cen[1]))
            else:
                self.cen=[int(self.ui.refCenLineEdit.text().split(',')[0]),int(self.ui.refCenLineEdit.text().split(',')[1])]  
            self.selData={}
            self.selErrorData={}
            if self.dir=='H':
                self.bg=self.bg*self.ccd_Dist[self.selFrameNums[0]]*np.pi/180/0.06/slitx
            else:
                self.bg=self.bg*self.ccd_Dist[self.selFrameNums[0]]*np.pi/180/0.06/slity   
            for i in self.selFrameNums:
                self.selData[i]=self.ccdData[i]
                self.selErrorData[i]=self.ccdErrorData[i]
            self.bruker.setROI(self.selData,self.selErrorData,absfac=self.absfac,absnum=self.absnum,slit=self.slit,cen=self.cen,bg=self.bg,dir=self.dir,mon=self.mon)
            slit=[self.slit[1], self.slit[0]]
            cen=[self.cen[1], self.cen[0]]
            self.ax={}
            if self.dir=='H':
                self.refextent=[self.bruker.bglcen-slit[1], self.bruker.bgrcen+slit[1]+1, cen[0]+slit[0]+1, cen[0]-slit[0]]
                self.ax[1]=self.ui.refADDataPlotWidget.canvas.fig.add_subplot(2,1,1)
                self.vmax=np.max(self.bruker.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.bruker.bglcen-slit[1]:self.bruker.bgrcen+slit[1]+1])
                self.vmin=np.max(self.bruker.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.bruker.bglcen-slit[1]:self.bruker.bgrcen+slit[1]+1])
                self.refZdata=self.bruker.imageData
                p=self.ax[1].imshow(self.bruker.imageROI[cen[0]-slit[0]:cen[0]+slit[0]+1, self.bruker.bglcen-slit[1]:self.bruker.bgrcen+slit[1]+1], origin='upper', aspect='auto', vmax=self.vmax, cmap=cmap, extent=self.refextent, interpolation='nearest')
                self.ax[1].set_ylabel('w/o Corr')
                self.ax[1].format_coord=self.format_coord_ref
            else:
                self.refextent=[cen[1]+slit[1]+1, cen[1]-slit[1], self.bruker.bgucen-slit[0], self.bruker.bgdcen+slit[0]+1]
                self.ax[1]=self.ui.refADDataPlotWidget.canvas.fig.add_subplot(1,2,1)            
                self.vmax=np.max(self.bruker.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.bruker.bglcen-slit[1]:self.bruker.bgrcen+slit[1]+1])
                self.vmin=np.max(self.bruker.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.bruker.bglcen-slit[1]:self.bruker.bgrcen+slit[1]+1])
                self.refZdata=self.bruker.imageData
                p=self.ax[1].imshow(self.bruker.imageROI[self.bruker.bgucen-slit[0]:self.bruker.bgdcen+slit[0]+1,cen[1]-slit[1]:cen[1]+slit[1]+1], origin='upper', aspect='auto', vmax=self.vmax, cmap=cmap, extent=self.refextent, interpolation='nearest')
                self.ax[1].set_ylabel('w/o Corr')
                self.ax[1].format_coord=self.format_coord_ref
            self.ui.refADDataPlotWidget.canvas.fig.colorbar(p)
            self.ui.refADDataPlotWidget.canvas.draw()
        elif self.det=='Pilatus': 
            slitx=int(self.ui.refSlitLineEdit.text().split(',')[0])
            slity=int(self.ui.refSlitLineEdit.text().split(',')[1])
            self.slit=[slitx,slity]
            self.bg=float(self.ui.refBGOffLineEdit.text())           
            self.selFrameNums=self.refSortedFrameNums[self.ui.refQzListWidget.row(self.refQzSelectedItems[0])]
            cenx=self.pilatus.NCOLS-self.pilX[self.selFrameNums[0]]-1 #image was flipped left to right
            ceny=self.pilY[self.selFrameNums[0]]
            self.mon=[self.pilMonc[i] for i in self.selFrameNums]
            self.absnum=[self.pil_AbsNum[i] for i in self.selFrameNums]
            self.vmax=float(self.ui.pilMaxLineEdit.text())
            self.vmin=float(self.ui.pilMinLineEdit.text())
            self.cen=[cenx,ceny]
            if self.ui.refCenCheckBox.checkState()==0:
                self.ui.refCenLineEdit.setText(str(self.cen[0])+','+str(self.cen[1]))
            else:
                self.cen=[int(self.ui.refCenLineEdit.text().split(',')[0]),int(self.ui.refCenLineEdit.text().split(',')[1])]  
            self.selData={}
            self.selErrorData={}
            if self.dir=='H':
                self.bg=self.bg*self.pil_Dist[self.selFrameNums[0]]*np.pi/180/0.172/slitx
            else:
                self.bg=self.bg*self.pil_Dist[self.selFrameNums[0]]*np.pi/180/0.172/slity         
            for i in self.selFrameNums:
                self.selData[i]=self.pilData[i]
                self.selErrorData[i]=self.pilErrorData[i]
            self.pilatus.setROI(self.selData,self.selErrorData,absfac=self.absfac,absnum=self.absnum,slit=self.slit,cen=self.cen,bg=self.bg,dir=self.dir,mon=self.mon)
            slit=[self.slit[1], self.slit[0]]
            cen=[self.cen[1], self.cen[0]]
            self.ax={}
            if self.dir=='H':
                self.refextent=[max(0,self.pilatus.bglcen-slit[1]), self.pilatus.bgrcen+slit[1]+1, cen[0]-slit[0], cen[0]+slit[0]+1]
                self.ax[1]=self.ui.refADDataPlotWidget.canvas.fig.add_subplot(2,1,1)
                self.vmax=np.max(self.pilatus.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, max(0,self.pilatus.bglcen-slit[1]):self.pilatus.bgrcen+slit[1]+1])
               # self.vmin=np.max(self.pilatus.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, max(0,self.pilatus.bglcen-slit[1]):self.pilatus.bgrcen+slit[1]+1])
                self.refZdata=self.pilatus.imageData
                p=self.ax[1].imshow(self.pilatus.imageROI[cen[0]-slit[0]:cen[0]+slit[0]+1, max(0,self.pilatus.bglcen-slit[1]):self.pilatus.bgrcen+slit[1]+1], origin='lower', aspect='auto', vmax=self.vmax, cmap=cmap, extent=self.refextent, interpolation='nearest')
                self.ax[1].set_ylabel('w/o Corr')
                self.ax[1].format_coord=self.format_coord_ref
            else:
                self.refextent=[cen[1]+slit[1]+1, cen[1]-slit[1], self.pilatus.bgucen-slit[0], self.pilatus.bgdcen+slit[0]+1]
                self.ax[1]=self.ui.refADDataPlotWidget.canvas.fig.add_subplot(1,2,1)            
                self.vmax=np.max(self.pilatus.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.pilatus.bglcen-slit[1]:self.pilatus.bgrcen+slit[1]+1])
                #self.vmin=np.max(self.pilatus.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.pilatus.bglcen-slit[1]:self.pilatus.bgrcen+slit[1]+1])
                self.refZdata=self.pilatus.imageData
                p=self.ax[1].imshow(self.pilatus.imageROI[self.pilatus.bgucen-slit[0]:self.pilatus.bgdcen+slit[0]+1,cen[1]-slit[1]:cen[1]+slit[1]+1], origin='lower', aspect='auto', vmax=self.vmax, cmap=cmap, extent=self.refextent, interpolation='nearest')
                self.ax[1].set_ylabel('w/o Corr')
                self.ax[1].format_coord=self.format_coord_ref
            self.ui.refADDataPlotWidget.canvas.fig.colorbar(p)
            self.ui.refADDataPlotWidget.canvas.draw()
         
        
        
    def refAnalyze(self):
        self.refQzListSelectionChanged()
        cmap=str(self.ui.gixCMapComboBox.currentText())
        if self.det=='Bruker':
            self.ui.refBadPixListWidget.clear()
            fac=float(self.ui.refSeaWinLineEdit.text())
            self.cenfit=0
            if self.ui.refCenCheckBox.checkState()!=0:
                self.cen=[int(self.ui.refCenLineEdit.text().split(',')[0]),int(self.ui.refCenLineEdit.text().split(',')[1])]
                self.cenfit=1
            self.par=self.bruker.peakFind(self.selData,self.selErrorData,slit=self.slit,cen=self.cen,absfac=self.absfac,absnum=self.absnum,fac=fac,mon=self.mon,cenfit=self.cenfit)
            self.cen=[int(np.floor(self.par[1])), int(np.floor(self.par[3]))]
            slit=[self.slit[1], self.slit[0]]
            cen=[self.cen[1], self.cen[0]]
            bfac=float(self.ui.refBPCFacLineEdit.text())
            try:
                sig,sigerr,lbg,lbgerr,rbg,rbgerr=self.bruker.badPix_corr(par=self.par,slit=self.slit,cen=self.cen,bad=1,bfac=bfac,dir=self.dir,bg=self.bg,plot=0)
            except:
                self.messageBox('Warning:: Please change to different search window by +/- 0.5.')
            if self.dir=='H':
                self.ui.refBadPixListWidget.addItem('------Signal------')
                for items in self.bruker.badSig:
                    self.ui.refBadPixListWidget.addItem(str(items[0])+'\t'+str(items[1]))
                self.ui.refBadPixListWidget.addItem('------Left BG------')
                for items in self.bruker.badLeft:
                    self.ui.refBadPixListWidget.addItem(str(items[0])+'\t'+str(items[1]))
                self.ui.refBadPixListWidget.addItem('------Right BG------')
                for items in self.bruker.badRight:
                    self.ui.refBadPixListWidget.addItem(str(items[0])+'\t'+str(items[1]))
            else:
                self.ui.refBadPixListWidget.addItem('------Signal------')
                for items in self.bruker.badSig:
                    self.ui.refBadPixListWidget.addItem(str(items[0])+'\t'+str(items[1]))
                self.ui.refBadPixListWidget.addItem('------Up BG------')
                for items in self.bruker.badUp:
                    self.ui.refBadPixListWidget.addItem(str(items[0])+'\t'+str(items[1]))
                self.ui.refBadPixListWidget.addItem('-----Down BG------')
                for items in self.bruker.badDown:
                    self.ui.refBadPixListWidget.addItem(str(items[0])+'\t'+str(items[1]))
            self.bruker.setROI(self.bruker.imageData,self.bruker.errorData,slit=self.slit,cen=self.cen,absfac=self.absfac,absnum=self.absnum,bg=self.bg,dir=self.dir,mon=None)  
            if self.dir=='H':            
                self.ax[2]=self.ui.refADDataPlotWidget.canvas.fig.add_subplot(2,1,2)
                self.refextent=[ self.bruker.bglcen-slit[1], self.bruker.bgrcen+slit[1]+1, cen[0]+slit[0]+1, cen[0]-slit[0]]
                self.vmax=np.max(self.bruker.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.bruker.bglcen-slit[1]:self.bruker.bgrcen+slit[1]+1])
                self.vmin=np.min(self.bruker.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.bruker.bglcen-slit[1]:self.bruker.bgrcen+slit[1]+1])
                p=self.ax[2].imshow(self.bruker.imageROI[cen[0]-slit[0]:cen[0]+slit[0]+1, self.bruker.bglcen-slit[1]:self.bruker.bgrcen+slit[1]+1], origin='upper', aspect='auto', vmax=self.vmax, cmap=cmap, extent=self.refextent, interpolation='nearest')
                self.ax[2].set_ylabel('w Corr')
                self.ax[2].format_coord=self.format_coord_ref
            else:
                self.ax[2]=self.ui.refADDataPlotWidget.canvas.fig.add_subplot(1,2,2)
                self.ax[2].clear()
                self.refextent=[cen[1]+slit[1]+1, cen[1]-slit[1], self.bruker.bgucen-slit[0], self.bruker.bgdcen+slit[0]+1]
                self.vmax=np.max(self.bruker.imageData[self.bruker.bgucen-slit[0]:self.bruker.bgdcen+slit[0]+1,cen[1]-slit[1]:cen[1]+slit[1]+1])
                self.vmin=np.min(self.bruker.imageData[self.bruker.bgucen-slit[0]:self.bruker.bgdcen+slit[0]+1,cen[1]-slit[1]:cen[1]+slit[1]+1])
                p=self.ax[2].imshow(self.bruker.imageROI[self.bruker.bgucen-slit[0]:self.bruker.bgdcen+slit[0]+1,cen[1]-slit[1]:cen[1]+slit[1]+1], origin='upper', aspect='auto', vmax=self.vmax, cmap=cmap, extent=self.refextent, interpolation='nearest')
                self.ax[2].set_ylabel('w Corr')
                self.ax[2].format_coord=self.format_coord_ref
        elif self.det=='Pilatus':
            self.ui.refBadPixListWidget.clear()
            fac=float(self.ui.refSeaWinLineEdit.text())
            self.cenfit=1
            #Peak searching is dummy it is not doing anything#
            if self.ui.refCenCheckBox.checkState()!=0:
                self.cen=[int(self.ui.refCenLineEdit.text().split(',')[0]),int(self.ui.refCenLineEdit.text().split(',')[1])]
            self.par=self.pilatus.peakFind(self.selData,self.selErrorData,slit=self.slit,cen=self.cen,absfac=self.absfac,absnum=self.absnum,fac=fac,mon=self.mon,cenfit=self.cenfit)
            self.cen=[int(np.floor(self.par[1])), int(np.floor(self.par[3]))]
            slit=[self.slit[1], self.slit[0]]
            cen=[self.cen[1], self.cen[0]]
            bfac=float(self.ui.refBPCFacLineEdit.text())
            sig,sigerr,lbg,lbgerr,rbg,rbgerr=self.pilatus.sumROI(slit=self.slit,cen=self.cen,dir=self.dir,bg=self.bg)
            #print sig,sigerr,lbg,lbgerr,rbg,rbgerr
            self.pilatus.setROI(self.pilatus.imageData,self.pilatus.errorData,slit=self.slit,cen=self.cen,absfac=self.absfac,absnum=self.absnum,bg=self.bg,dir=self.dir,mon=None)  
            if self.dir=='H':
                self.refextent=[max(0,self.pilatus.bglcen-slit[1]), self.pilatus.bgrcen+slit[1]+1, cen[0]-slit[0], cen[0]+slit[0]+1]
                self.ax[2]=self.ui.refADDataPlotWidget.canvas.fig.add_subplot(2,1,2)
                self.vmax=np.max(self.pilatus.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, max(0,self.pilatus.bglcen-slit[1]):self.pilatus.bgrcen+slit[1]+1])
               # self.vmin=np.max(self.pilatus.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.pilatus.bglcen-slit[1]:self.pilatus.bgrcen+slit[1]+1])
                self.refZdata=self.pilatus.imageData
                p=self.ax[2].imshow(self.pilatus.imageROI[cen[0]-slit[0]:cen[0]+slit[0]+1, max(0,self.pilatus.bglcen-slit[1]):self.pilatus.bgrcen+slit[1]+1], origin='lower', aspect='auto', vmax=self.vmax, cmap=cmap, extent=self.refextent, interpolation='nearest')
                self.ax[2].set_ylabel('w Corr')
                self.ax[2].format_coord=self.format_coord_ref
            else:
                self.refextent=[cen[1]+slit[1]+1, cen[1]-slit[1], self.pilatus.bgucen-slit[0], self.pilatus.bgdcen+slit[0]+1]
                self.ax[2]=self.ui.refADDataPlotWidget.canvas.fig.add_subplot(1,2,2)            
                self.vmax=np.max(self.pilatus.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.pilatus.bglcen-slit[1]:self.pilatus.bgrcen+slit[1]+1])
                #self.vmin=np.max(self.pilatus.imageData[cen[0]-slit[0]:cen[0]+slit[0]+1, self.pilatus.bglcen-slit[1]:self.pilatus.bgrcen+slit[1]+1])
                self.refZdata=self.pilatus.imageData
                p=self.ax[2].imshow(self.pilatus.imageROI[self.pilatus.bgucen-slit[0]:self.pilatus.bgdcen+slit[0]+1,cen[1]-slit[1]:cen[1]+slit[1]+1], origin='lower', aspect='auto', vmax=self.vmax, cmap=cmap, extent=self.refextent, interpolation='nearest')
                self.ax[2].set_ylabel('w Corr')
                self.ax[2].format_coord=self.format_coord_ref
        self.ui.refADDataPlotWidget.canvas.fig.colorbar(p)
        self.ui.refADDataPlotWidget.canvas.draw()
        self.fsig=sig-(lbg+rbg)/2.0
        self.fsigerr=np.sqrt(sigerr**2+(lbgerr**2+rbgerr**2)/4)
        self.fsigerrper=self.fsigerr/self.fsig*100.0
        self.ui.refQLineEdit.setText('%.4f'%float(str(self.refQzSelectedItems[0].text()).split()[1]))
        self.ui.refRLineEdit.setText('%.4e'%self.fsig)
        self.ui.refRErrLineEdit.setText('%.4e'%self.fsigerr)
        self.ui.refPerRErrLineEdit.setText('%.4f'%self.fsigerrper)
    
    def updateRefDataList(self):
        if self.refpatchindex==0: #regular ref
            comments=str(self.refQzSelectedItems[0].text()).split('\t')[1:]
            self.refData.append([float(str(self.refQzSelectedItems[0].text()).split()[1]),self.fsig,self.fsigerr])
            self.refInfo.append('#'+comments[0]+'\t'+comments[1]+'\t#Comments:'+str(self.ui.refComLineEdit.text()))
            self.ui.refDataListWidget.addItem('%.4f'%float(str(self.refQzSelectedItems[0].text()).split()[1])+'\t'+'%.4e'%self.fsig+'\t'+'%.4e'%self.fsigerr+'\t'+self.refInfo[-1])
        else:  #pilatus patch ref
            self.ui.refDataListWidget.clear()
            self.ui.refNormRefListWidget.clear()
            for i in range(len(self.refData)):
                self.ui.refDataListWidget.addItem('%.4f'%self.refData[i][0]+'\t'+'%.4e'%self.refData[i][1]+'\t'+'%.4e'%self.refData[i][2])
        self.updateRefPlotData()
        
    def refAnalyzeAll(self):
        self.refAnalyze()
        self.updateRefDataList()
        self.progressDialog=QProgressDialog('Analyzing Images for reflectivity calculations','Abort',1,self.ui.refQzListWidget.count())
        self.progressDialog.setWindowModality(Qt.WindowModal)
        self.progressDialog.setWindowTitle('Wait')
        self.progressDialog.setAutoClose(True)
        self.progressDialog.setAutoReset(True)
        self.progressDialog.setMinimum(1)
        self.progressDialog.setMaximum(self.ui.refQzListWidget.count())
        self.progressDialog.show()
        for i in range(1,self.ui.refQzListWidget.count()):
            self.disconnect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
            self.ui.refQzListWidget.setItemSelected(self.ui.refQzListWidget.item(i-1),False)
            self.connect(self.ui.refQzListWidget, SIGNAL('itemSelectionChanged()'),self.refQzListSelectionChanged)
            self.ui.refQzListWidget.setItemSelected(self.ui.refQzListWidget.item(i),True)
            self.refAnalyze()
            self.updateRefDataList()
            self.progressDialog.setLabelText('Reading Frames for Qz= '+str(self.refQzSelected))     
            self.updateProgress()
            if self.progressDialog.wasCanceled()==True:
                break
        self.progressDialog.hide()
            
            
        
    def delRefDataList(self):
        items=self.ui.refDataListWidget.selectedItems()
        for item in items:
            self.refData.pop(self.ui.refDataListWidget.row(item))
            if self.refpatchindex==0:
                self.refInfo.pop(self.ui.refDataListWidget.row(item))
            self.ui.refDataListWidget.takeItem(self.ui.refDataListWidget.row(item))
        self.updateRefPlotData()
        
    def fresnel(self,q,qc):
        fre=[]
        for qx in q:
            if qx<qc:
                fre.append(1)
            else:
                fre.append(((qx-np.sqrt(qx**2-qc**2))/(qx+np.sqrt(qx**2-qc**2)))**2)
        return np.array(fre)
        
    def pilRefPatchData(self):
        self.refpatchindex=1
        slitx=int(self.ui.refSlitLineEdit.text().split(',')[0]) #get the size of ROI
        slity=int(self.ui.refSlitLineEdit.text().split(',')[1])
        self.slit=[slitx,slity]        
        self.dir=str(self.ui.refBGDirComboBox.currentText()) #get the background offset direction
        bg=float(self.ui.refBGOffLineEdit.text())  #get the background offset
        self.command='Ref Patch, scans=['+str([item for item in np.sort(self.selectedScanNums)])[1:-1]+'], slits=['+str([item for item in self.slit])[1:-1]+'], bg_angle='+str(bg)+', bg_direction='+str(self.dir)
        self.ui.commandLineEdit.setText(self.command)         #update the command lineedit
        self.scanqzframe={}  #build a dictionary for all scans with their qz ranges and data for sorting the scans.
        self.pilRefPeakCen=[]
        for i in self.selectedPilFramesNums:
            qkey=self.pilScanNum[i]
            scenx=self.pilatus.NCOLS-self.pilX[i]-1 #image was flipped left to right
            sceny=self.pilY[i]
            if qkey in self.scanqzframe:
                selData=self.pilData[i]
                selErrorData=self.pilErrorData[i]
                cenx,ceny=self.pilatus.peakLocate(selData,selErrorData,cen=[scenx,sceny],mon=self.pilMonc[i])
                if np.abs(cenx-scenx)>2:
                    self.pilRefPeakCen.append([qkey,self.pilFrameNum[i],'x', np.abs(cenx-scenx)])                    
                    #self.messageBox('The reflected beam center is off by '+str(np.abs(cenx-scenx))+ ' pixels from the default center \nalong x directions for Scan #'+str(qkey)+' Frame #'+str(self.pilFrameNum[i])+'!')
                if np.abs(ceny-sceny)>2:
                    self.pilRefPeakCen.append([qkey,self.pilFrameNum[i],'y', np.abs(ceny-sceny)])  
                    #self.messageBox('The reflected beam center is off by '+str(np.abs(ceny-sceny))+ ' pixels from the default center \nalong y directions for Scan #'+str(qkey)+' Frame #'+str(self.pilFrameNum[i])+'!')
                if self.dir=='H':
                    self.bg=bg*self.pil_Dist[i]*np.pi/180/0.172/slitx
                else:
                    self.bg=bg*self.pil_Dist[i]*np.pi/180/0.172/slity
                if self.pilFileQzs[i]<0.03 and self.pilFileQzs[i]> 0.02 and self.ui.pilLowQzCorCheckBox.checkState()!=0:
                    self.slit=[slitx,5]
                else:
                    self.slit=[slitx,slity]    
                sig,sigerr,lbg,lbgerr,rbg,rbgerr=self.pilatus.sumROI(slit=self.slit,cen=[cenx,ceny],dir=self.dir,bg=self.bg)
                #self.scanqzframe[qkey].append(np.vstack((self.pilFileQzs[i],sig-(lbg+rbg)/2.0, np.sqrt(sigerr**2+(lbgerr**2+rbgerr**2)/4))).T)
                self.scanqzframe[qkey]=np.vstack((self.scanqzframe[qkey],(np.array([self.pilFileQzs[i],sig-(lbg+rbg)/2.0, np.sqrt(sigerr**2+(lbgerr**2+rbgerr**2)/4),sig,sigerr,lbg,rbg]))))
            else:
                selData=self.pilData[i]
                selErrorData=self.pilErrorData[i]
                cenx,ceny=self.pilatus.peakLocate(selData,selErrorData,cen=[scenx,sceny],mon=self.pilMonc[i])
                if np.abs(cenx-scenx)>2:
                    self.pilRefPeakCen.append([qkey,self.pilFrameNum[i],'x', np.abs(cenx-scenx)])    
                   # self.messageBox('The reflected beam center is off by more than 3 pixels from the default centen \nalong x directions for Scan #'+str(qkey)+' Frame #'+str(self.pilFrameNum[i])+'!')
                if np.abs(ceny-sceny)>2:
                    self.pilRefPeakCen.append([qkey,self.pilFrameNum[i],'y', np.abs(ceny-sceny)])  
                   # self.messageBox('The reflected beam center is off by more than 3 pixels from the default centen \nalong y directions for Scan #'+str(qkey)+' Frame #'+str(self.pilFrameNum[i])+'!')
                if self.dir=='H':
                    self.bg=bg*self.pil_Dist[i]*np.pi/180/0.172/slitx
                else:
                    self.bg=bg*self.pil_Dist[i]*np.pi/180/0.172/slity
                if self.pilFileQzs[i]<0.03 and self.pilFileQzs[i]> 0.02:
                    self.slit=[slitx,5]
                else:
                    self.slit=[slitx,slity]    
                sig,sigerr,lbg,lbgerr,rbg,rbgerr=self.pilatus.sumROI(slit=self.slit,cen=[cenx,ceny],dir=self.dir,bg=self.bg)  
                #self.scanqzframe[qkey]=[np.vstack((self.pilFileQzs[i],sig-(lbg+rbg)/2.0, np.sqrt(sigerr**2+(lbgerr**2+rbgerr**2)/4))).T]
                self.scanqzframe[qkey]=np.array([self.pilFileQzs[i],sig-(lbg+rbg)/2.0, np.sqrt(sigerr**2+(lbgerr**2+rbgerr**2)/4),sig,sigerr,lbg,rbg])
        self.pilrefsort=self.sortedScans(self.scanqzframe)
       # print self.scanqzframe
        #print self.pilrefsort
        if len(self.pilRefPeakCen)!=0:
            self.pilRefPeakCenDis()
        if len(self.pilrefsort)==0:
            self.messageBox('Warning: no frame selected!!!')
        elif len(self.pilrefsort)==1:
            self.refData=list(self.scanqzframe[self.pilrefsort[0]])
           #print self.refData
            self.ui.PlotWidget.setCurrentIndex(3)
            self.updateRefDataList()
        else:
            self.pilrefscanrag={}  #find the identical scans
            for key in self.scanqzframe:
                self.pilrefscanrag[key]=str(format(self.scanqzframe[key][0][0],'.3f'))+'_'+str(format(self.scanqzframe[key][-1][0],'.3f'))+'_'+str(len(self.scanqzframe[key]))
          #  print self.pilrefscanrag
            dd={}
            for key,value in self.pilrefscanrag.items():
                try:
                    dd[value].append(key)
                except:
                    dd[value]=[key]
          #  print dd
            self.pilrefsamescannum=[l for l in dd.values() if len(l)>1] #this is list with the identical scans numbers
          #  print self.pilrefsamescannum
            if len(self.pilrefsamescannum)!=0: #if identical scans found, give user choice to merge, patch, and drop scans. 
                self.messageBox('Identical scans were found! Please select your choices; regular patch is defult!')
                self.pilRefSameScan()
            else:
                self.pilRefPatch()
                
    def pilRefPeakCenDis(self):
        Dialog=QDialog(self)
        ui=uic.loadUi('pilPeakCenDialog.ui',Dialog)
        ui.show()
        line='Scan\tFrame\tAxis\t# of Pixels off\n'
        for i in range(len(self.pilRefPeakCen)):
            line=line+str(self.pilRefPeakCen[i][0])+'\t'+str(self.pilRefPeakCen[i][1])+'\t'+self.pilRefPeakCen[i][2]+'\t'+str(self.pilRefPeakCen[i][3])+'\n'
        ui.peakCenTextBrowser.append(line)
        
        Dialog.exec_()              
                
    def pilRefPatch(self):            
        Dialog=QDialog(self)                
        self.uirefpatch=uic.loadUi('refpatch.ui', Dialog)
        self.uirefpatch.show()
        self.pilreforder=0
        self.connect(self.uirefpatch.dropPushButton, SIGNAL('clicked()'),self.pilRefPatchDrop)  #drop the data
        self.connect(self.uirefpatch.nextPushButton,SIGNAL('clicked()'),self.pilRefPatchNext) #go to the next overlap
        self.pilRefPatchPlot()
        
    def pilRefSameScan(self):
        Dialog=QDialog(self)
        self.uirefpatchsamescan=uic.loadUi('refpatchsamescan.ui', Dialog)
        self.uirefpatchsamescan.show()
        self.pilrefsamescanorder=0
        self.connect(self.uirefpatchsamescan.dropPushButton, SIGNAL('clicked()'),self.pilRefSSDrop) #drop scans
        self.connect(self.uirefpatchsamescan.mergePushButton, SIGNAL('clicked()'),self.pilRefSSMerge) #merge scans
        self.connect(self.uirefpatchsamescan.nextPushButton, SIGNAL('clicked()'),self.pilRefSSNext) #go to the next overlap
        self.pilRefSameScanPlot()
    
    def pilRefSameScanPlot(self):  
        if self.pilrefsamescanorder==len(self.pilrefsamescannum):
            self.pilrefsort=self.sortedScans(self.scanqzframe)
            self.uirefpatchsamescan.close()
            self.pilRefPatch()
        else:
            self.uirefpatchsamescan.Label.setText('Indentical Reflectivity Scans: '+str([item for item in np.sort(self.pilrefsamescannum[self.pilrefsamescanorder])])[1:-1])
            if self.pilrefsamescanorder==len(self.pilrefsamescannum)-1:
                self.uirefpatchsamescan.nextPushButton.setText('Finish!') 
            self.uirefpatchsamescan.plotWidget.canvas.ax.clear()
            self.uirefpatchsamescan.plotWidget.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
            self.uirefpatchsamescan.plotWidget.canvas.ax.set_ylabel('Reflectivity')
            for i in range(len(self.pilrefsamescannum[self.pilrefsamescanorder])):
                self.uirefpatchsamescan.plotWidget.canvas.ax.errorbar(self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][i]][:,0],self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][i]][:,1], self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][i]][:,2], fmt='o-', label='#'+str(self.pilrefsamescannum[self.pilrefsamescanorder][i]))
            self.uirefpatchsamescan.plotWidget.canvas.ax.legend(frameon=False,scatterpoints=0,numpoints=1)
            self.uirefpatchsamescan.plotWidget.canvas.draw()

    def pilRefSSDrop(self):
        dropnumbers=str(self.uirefpatchsamescan.dropLineEdit.text()).split(',')
        dropnumbers=map(int,dropnumbers)
        for i in range(len(dropnumbers)):
            try:
                self.scanqzframe.pop(dropnumbers[i])
                self.pilrefsamescannum[self.pilrefsamescanorder].remove(dropnumbers[i])
            except:
                pass
        self.pilRefSameScanPlot()            
        #print  self.scanqzframe
        #print  self.pilrefsamescannum
        
    
    def pilRefSSMerge(self):
        datax=self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][0]][:,0]
        datay=self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][0]][:,1]
        datayerr=self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][0]][:,2]**2
        for i in range(1,len(self.pilrefsamescannum[self.pilrefsamescanorder])):
            datay=datay+self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][i]][:,1]
            datayerr=datayerr+self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][i]][:,2]**2
            self.scanqzframe.pop(self.pilrefsamescannum[self.pilrefsamescanorder][i])
        datay=datay/len(self.pilrefsamescannum[self.pilrefsamescanorder])
        datayerr=np.sqrt(datayerr)/len(self.pilrefsamescannum[self.pilrefsamescanorder])
        data=np.vstack((datax,datay,datayerr)).T
        self.scanqzframe[self.pilrefsamescannum[self.pilrefsamescanorder][0]]=data
        self.pilrefsamescannum[self.pilrefsamescanorder]=[self.pilrefsamescannum[self.pilrefsamescanorder][0]]
        #print data
        #print self.pilrefsamescannum
        #print self.scanqzframe
        self.pilRefSameScanPlot()
        
    def pilRefSSNext(self):
        self.pilrefsamescanorder=self.pilrefsamescanorder+1 
        self.pilRefSameScanPlot()
    
    def pilRefPatchPlot(self):
        if self.pilreforder==len(self.pilrefsort)-1:
            self.uirefpatch.close()
            self.refData=list(self.pilrefdata)
            self.ui.PlotWidget.setCurrentIndex(3)
          #  print self.pilrefdata
            self.updateRefDataList()
        elif self.pilreforder==0:
            if self.pilreforder==len(self.pilrefsort)-2:
                self.uirefpatch.nextPushButton.setText('Finish!')
            self.uirefpatch.Label.setText('Reflectivity Patch between Scans '+str(self.pilrefsort[self.pilreforder])+' and '+str(self.pilrefsort[self.pilreforder+1]))
            if self.scanqzframe[self.pilrefsort[self.pilreforder]][-1][0]<self.scanqzframe[self.pilrefsort[self.pilreforder+1]][0][0]:
                self.messageBox('Warning: no overlap qz between slected frames from scan '+str(self.pilrefsort[self.pilreforder])+' and '+str(self.pilrefsort[self.pilreforder+1])+'!!!')
            else:
               #print self.scanqzframe[self.pilrefsort[self.pilreforder]],self.scanqzframe[self.pilrefsort[self.pilreforder+1]]
                fac,xmin,xmax,ymin,ymax=self.oneDjoin(self.scanqzframe[self.pilrefsort[self.pilreforder]],self.scanqzframe[self.pilrefsort[self.pilreforder+1]],mtype='fac')
                #print self.scanqzframe[self.pilrefsort[self.pilreforder]],self.scanqzframe[self.pilrefsort[self.pilreforder+1]]                
                self.uirefpatch.plotWidget.canvas.ax.clear()
                self.uirefpatch.plotWidget.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
                self.uirefpatch.plotWidget.canvas.ax.set_ylabel('Reflectivity')
                self.uirefpatch.plotWidget.canvas.ax.errorbar(self.scanqzframe[self.pilrefsort[self.pilreforder]][:,0],self.scanqzframe[self.pilrefsort[self.pilreforder]][:,1],self.scanqzframe[self.pilrefsort[self.pilreforder]][:,2],fmt='o-',label='#'+str(self.pilrefsort[self.pilreforder]))
                self.uirefpatch.plotWidget.canvas.ax.errorbar(self.scanqzframe[self.pilrefsort[self.pilreforder+1]][:,0],self.scanqzframe[self.pilrefsort[self.pilreforder+1]][:,1]*fac,self.scanqzframe[self.pilrefsort[self.pilreforder+1]][:,2]*fac,fmt='o-',label='#'+str(self.pilrefsort[self.pilreforder+1]))
                for i in range(len(self.scanqzframe[self.pilrefsort[self.pilreforder]])):
                    if self.scanqzframe[self.pilrefsort[self.pilreforder]][i,0]>=xmin:
                        self.uirefpatch.plotWidget.canvas.ax.text(self.scanqzframe[self.pilrefsort[self.pilreforder]][i,0]+0.02*(xmax-xmin),self.scanqzframe[self.pilrefsort[self.pilreforder]][i,1]+0.02*(ymax-ymin),str(i-len(self.scanqzframe[self.pilrefsort[self.pilreforder]])),color='b',fontsize=16) 
                for i in range(len(self.scanqzframe[self.pilrefsort[self.pilreforder+1]])): 
                    if self.scanqzframe[self.pilrefsort[self.pilreforder+1]][i,0]<=xmax:
                        self.uirefpatch.plotWidget.canvas.ax.text(self.scanqzframe[self.pilrefsort[self.pilreforder+1]][i,0]-0.04*(xmax-xmin),self.scanqzframe[self.pilrefsort[self.pilreforder+1]][i,1]*fac-0.06*(ymax-ymin),str(i),color='r',fontsize=16)                
                self.uirefpatch.plotWidget.canvas.ax.legend(frameon=False,scatterpoints=0,numpoints=1)
                self.uirefpatch.plotWidget.canvas.ax.set_xlim(xmin,xmax)
                self.uirefpatch.plotWidget.canvas.ax.set_ylim(ymin,ymax)
              #  print xmin,xmax,ymin,ymax, fac
                self.uirefpatch.plotWidget.canvas.draw()
        else:
            if self.pilreforder==len(self.pilrefsort)-2:
                self.uirefpatch.nextPushButton.setText('Finish!')
            self.uirefpatch.Label.setText('Reflectivity Patch between Scans '+str(self.pilrefsort[self.pilreforder])+' and '+str(self.pilrefsort[self.pilreforder+1]))
            if self.scanqzframe[self.pilrefsort[self.pilreforder]][-1][0]<self.scanqzframe[self.pilrefsort[self.pilreforder+1]][0][0]:
                self.messageBox('Warning: no overlap qz between slected frames from scan '+str(self.pilrefsort[self.pilreforder])+' and '+str(self.pilrefsort[self.pilreforder+1])+'!!!')
            else:
                fac,xmin,xmax,ymin,ymax=self.oneDjoin(self.pilrefdata,self.scanqzframe[self.pilrefsort[self.pilreforder+1]],mtype='fac')
                self.uirefpatch.plotWidget.canvas.ax.clear()
                self.uirefpatch.plotWidget.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
                self.uirefpatch.plotWidget.canvas.ax.set_ylabel('Reflectivity')
                self.uirefpatch.plotWidget.canvas.ax.errorbar(self.pilrefdata[:,0],self.pilrefdata[:,1],self.pilrefdata[:,2],fmt='o-',label='#'+str(self.pilrefsort[self.pilreforder]))
                self.uirefpatch.plotWidget.canvas.ax.errorbar(self.scanqzframe[self.pilrefsort[self.pilreforder+1]][:,0],self.scanqzframe[self.pilrefsort[self.pilreforder+1]][:,1]*fac,self.scanqzframe[self.pilrefsort[self.pilreforder+1]][:,2]*fac,fmt='o-',label='#'+str(self.pilrefsort[self.pilreforder+1]))
                for i in range(len(self.pilrefdata)):
                    if self.pilrefdata[i,0]>=xmin:
                        self.uirefpatch.plotWidget.canvas.ax.text(self.pilrefdata[i,0]+0.02*(xmax-xmin),self.pilrefdata[i,1]+0.02*(ymax-ymin),str(i-len(self.pilrefdata)),color='b',fontsize=16) 
                for i in range(len(self.scanqzframe[self.pilrefsort[self.pilreforder+1]])): 
                    if self.scanqzframe[self.pilrefsort[self.pilreforder+1]][i,0]<=xmax:
                        self.uirefpatch.plotWidget.canvas.ax.text(self.scanqzframe[self.pilrefsort[self.pilreforder+1]][i,0]-0.04*(xmax-xmin),self.scanqzframe[self.pilrefsort[self.pilreforder+1]][i,1]*fac-0.06*(ymax-ymin),str(i),color='r',fontsize=16)                
                self.uirefpatch.plotWidget.canvas.ax.legend(frameon=False,scatterpoints=0,numpoints=1)
                self.uirefpatch.plotWidget.canvas.ax.set_xlim(xmin,xmax)
                self.uirefpatch.plotWidget.canvas.ax.set_ylim(ymin,ymax)
            #    print xmin,xmax,ymin,ymax, fac
                self.uirefpatch.plotWidget.canvas.draw()
    
    def pilRefPatchDrop(self):
        dropnumbers=str(self.uirefpatch.dropLineEdit.text()).split(',')
        dropnumbers=map(int,dropnumbers)
        dropnumbers.sort(key=abs,reverse=True) #sort the dropnumbers to make sure starting dropping the points away from the end.  
        for i in range(len(dropnumbers)):
            if dropnumbers[i]>=0:
                self.scanqzframe[self.pilrefsort[self.pilreforder+1]]=np.delete(self.scanqzframe[self.pilrefsort[self.pilreforder+1]],dropnumbers[i],axis=0)
            else:
                if self.pilreforder==0:
                    self.scanqzframe[self.pilrefsort[0]]=np.delete(self.scanqzframe[self.pilrefsort[0]],dropnumbers[i],axis=0)
                else:
                    self.pilrefdata=np.delete(self.pilrefdata,dropnumbers[i],axis=0)
        
      #  print self.pilrefdata
        self.pilRefPatchPlot()
        
        
    def pilRefPatchNext(self):
        if self.pilreforder==0:
            self.pilrefdata=self.oneDjoin(self.scanqzframe[self.pilrefsort[self.pilreforder]],self.scanqzframe[self.pilrefsort[self.pilreforder+1]],mtype='merge')
        else:
            self.pilrefdata=self.oneDjoin(self.pilrefdata,self.scanqzframe[self.pilrefsort[self.pilreforder+1]],mtype='merge')
        self.pilreforder=self.pilreforder+1
        self.pilRefPatchPlot()
        
     
    def oneDjoin(self,data1, data2, mtype='fac'): #join two 1D data
        if len(data1[0])>3:   #keep only first three columns
            for i in range(len(data1[0])-3):
                data1=np.delete(data1,-1,1)
        if len(data2[0])>3:   #keep only first three columns
            for i in range(len(data2[0])-3):
                data2=np.delete(data2,-1,1)
        data1x=data1[:,0] #get qz from data1
        data1y=data1[:,1] #get ref from data1
        data2x=data2[:,0] #get qz from data2
        data2y=data2[:,1] #get ref from data2
        
       # print data1
       # print data2
        overlap2=np.where(data2x<=data1x[-1]) #create an array for where data2 having smaller qz than points in data1
        #overlap1=np.where(np.logical_and(data1x>=data2x[0],data1x<=data2x[overlap2[0][-1]])) #create an array for where data1 having larger qz than the first point in data2 and smaller qz than the last point in data2[overlap2]
        overlap1=np.where(data1x>=data2x[0])  
       # print data1x[overlap1], data2x[overlap2]
        #f=interp1d(data2x[overlap2],data2y[overlap2], kind='linear')
        f=interp1d(data2x,data2y, kind='cubic')
        if mtype=='fac':
            fac=np.sum(data1y[overlap1])/np.sum(f(data1x[overlap1]))
           # print min(np.min(data1y[overlap1]),np.min(data2y[overlap2]*fac))
            if overlap1[0][0]==0:
                overlap1[0][0]=1
            if overlap2[0][-1]==len(data2x)-1:
                overlap2[0][-1]=len(data2x)-2
            return fac, data1x[overlap1[0][0]-1],data2x[overlap2[0][-1]+1],min(np.min(data1y[overlap1[0][0]-1:]),np.min(data2y[:overlap2[0][-1]+2]*fac)),max(np.max(data1y[overlap1[0][0]-1:]),np.max(data2y[:overlap2[0][-1]+2]*fac))
            
        elif mtype=='merge':
            data1yerr=data1[:,2] #get error from data1
            data2yerr=data2[:,2] #get error from data2
            #ferr=interp1d(data2x[overlap2],data2yerr[overlap2], kind='linear')
            ferr=interp1d(data2x,data2yerr, kind='cubic')
            x=data1x[overlap1]
            y=data1y[overlap1]/f(x)
            yerr=y*np.sqrt(((data1yerr[overlap1]/data1y[overlap1])**2+(ferr(x)/f(x))**2))
            pfit,pcov=curve_fit(self.polynomial,x,y,p0=[1],sigma=yerr) #fit the data with the scale factor 
            fac=pfit[0]
            print fac
            facerr=pcov[0][0]**0.5
            data2y=data2y*fac           #scale data2 with the fac and facerr
            data2yerr=np.sqrt(data2yerr**2*fac**2+data2y**2*facerr**2/fac**2)
            data2[:,1]=data2y  #this also changed the data self.scanqzframe[self.pilrefsort[self.pilreforder+1]]
            data2[:,2]=data2yerr
            data=data1[:overlap1[0][0],:] #get the first part of merged data from data1
            #f2=interp1d(data2x[overlap2],data2y[overlap2], kind='linear')
            #f2err=interp1d(data2x[overlap2],data2yerr[overlap2], kind='linear')
            f2=interp1d(data2x,data2y, kind='cubic')
            f2err=interp1d(data2x,data2yerr, kind='cubic')
            data3y=f2(data1x[overlap1])   #interpolate the data2 after the scaling
            data3yerr=f2err(data1x[overlap1])
            for i in range(overlap1[0][0],overlap1[0][-1]+1): #get the second part of merged data from data2
                data=np.vstack((data,(data1x[i], (data1y[i]+data3y[i-overlap1[0][0]])/2.0,np.sqrt(data1yerr[i]**2+data3yerr[i-overlap1[0][0]]**2)/2.0)))
            if data1x[-1]>=data2x[-1]:  #part 3 patching
                datalast=data1[overlap1[0][-1]+1:,:]
            else:
                datalast=data2[overlap2[0][-1]+1:,:]
            data=np.vstack((data,datalast))
        
            return data
         
    
    
    def addRefFiles(self):
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Files to import', directory=self.directory, filter='Ref Files (*_ref.txt;*.ref*;*_rrf.txt)')
        self.reffiles=self.reffiles+map(str, f)
        self.reffnames=[]
        self.ui.refRefFileListWidget.clear()
        for i in range(len(self.reffiles)):
            s=str(self.reffiles[i])
            self.reffnames.append(s[s.rfind('/')+1:])
            self.ui.refRefFileListWidget.addItem('#'+str(i+1)+'\t'+self.reffnames[i])

    def removeRefFiles(self):
        items=self.ui.refRefFileListWidget.selectedItems()
        #print items
        for item in items:
            #print self.reffiles, item
            #print self.ui.refRefFileListWidget.row(item)
            self.reffnames.pop(self.ui.refRefFileListWidget.row(item))
            self.reffiles.pop(self.ui.refRefFileListWidget.row(item))
            #self.ui.refRefFileListWidget.takeItem(self.ui.refRefFileListWidget.row(item))
        self.ui.refRefFileListWidget.clear()
        for i in range(len(self.reffnames)):
            self.ui.refRefFileListWidget.addItem('#'+str(i+1)+'\t'+self.reffnames[i])
    
    def updateSelectedRefFiles(self):
        selectedreffiles=self.ui.refRefFileListWidget.selectedItems()
        self.selectedreffiles_rows=[]
        for item in selectedreffiles:
            self.selectedreffiles_rows.append(self.ui.refRefFileListWidget.row(item))
        self.updateRefPlotData()

        
    def updateRefPlotData(self):
        self.ui.refRefDataPlotwidget.canvas.ax.clear()
        self.ui.refNormRefListWidget.clear()
        self.ui.refRefDataPlotwidget.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$')
        self.ui.refRefDataPlotwidget.canvas.ax.set_ylabel('Reflectivity')
        self.refQc=float(self.ui.refQcLineEdit.text())
        self.refQoff=float(self.ui.refQoffLineEdit.text())
        #print self.refpatchindex
        if len(self.refData)!=0:
            data=np.array(self.refData)
            if self.refpatchindex==0:
                self.refInfoNew=[self.refInfo[i] for i in np.argsort(data[:,0])]
                data=data[np.argsort(data[:,0])]
            if data[0,0]<self.refQc:
                norm=np.max(data[np.where(data[:,0]<self.refQc)][:,1]) #use the max intensity for the q below qc as the normalization factor. 
                
                #norm=data[0,1]
                #normerr=data[0,2]
                normerr=0
                data[:,1]=data[:,1]/norm
                data[:,2]=np.sqrt(data[:,2]**2/norm**2+normerr**2*data[:,1]**2/norm**2)
                perErr=data[:,2]*100.0/data[:,1]
                if self.refpatchindex==0:
                    self.ui.refNormRefListWidget.addItem('Q\tRef\tRef_Err\t%Err\t#S')
                    for i in range(len(data[:,0])):
                        self.ui.refNormRefListWidget.addItem('%.4f'%data[i,0]+'\t'+'%.4e'%data[i,1]+'\t'+'%.4e'%data[i,2]+'\t'+'%.4f'%perErr[i]+'\t'+self.refInfoNew[i])
                else:
                    self.ui.refNormRefListWidget.addItem('Q\tRef\tRef_Err\t%Err')
                    for i in range(len(data[:,0])):
                        self.ui.refNormRefListWidget.addItem('%.4f'%data[i,0]+'\t'+'%.4e'%data[i,1]+'\t'+'%.4e'%data[i,2]+'\t'+'%.4f'%perErr[i])
            if self.ui.refNormFacCheckBox.checkState()>0:
                norm=float(self.ui.refNormLineEdit.text())
                #normerr=np.sqrt(norm)
                data[:,1]=data[:,1]/norm
                data[:,2]=data[:,2]/norm                
                #data[:,2]=np.sqrt(data[:,2]**2/norm**2+normerr**2*data[:,1]**2/norm**4)
                perErr=data[:,2]*100.0/data[:,1]
                if self.refpatchindex==0:
                    self.ui.refNormRefListWidget.addItem('Q\tRef\tRef_Err\t%Err\t#S')
                    for i in range(len(data[:,0])):
                        self.ui.refNormRefListWidget.addItem('%.4f'%data[i,0]+'\t'+'%.4e'%data[i,1]+'\t'+'%.4e'%data[i,2]+'\t'+'%.4f'%perErr[i]+'\t'+self.refInfoNew[i])
                else:
                    self.ui.refNormRefListWidget.addItem('Q\tRef\tRef_Err\t%Err')
                    for i in range(len(data[:,0])):
                        self.ui.refNormRefListWidget.addItem('%.4f'%data[i,0]+'\t'+'%.4e'%data[i,1]+'\t'+'%.4e'%data[i,2]+'\t'+'%.4f'%perErr[i])        
            self.normRefData=data[np.argsort(data[:,0])]
            if self.ui.refRRFCheckBox.checkState()!=0:
                self.ui.refRefDataPlotwidget.canvas.ax.set_ylabel('R/RF')
                data[:,0]=data[:,0]+self.refQoff
                data[:,1]=data[:,1]/self.fresnel(data[:,0],self.refQc)
                data[:,2]=data[:,2]/self.fresnel(data[:,0],self.refQc)
                self.RrfData=data
            if self.ui.refQzSqrCheckBox.checkState()!=0:
                self.ui.refRefDataPlotwidget.canvas.ax.set_xlabel('Qz^2')
                data[:,0]=data[:,0]**2
            if len(data[0])!=3:  #for pilatus single scan data with left and right background
                self.ui.refRefDataPlotwidget.canvas.ax.errorbar(data[:,0],data[:,3],data[:,4],fmt='o',label='s')                    
                self.ui.refRefDataPlotwidget.canvas.ax.errorbar(data[:,0],data[:,5],fmt='-',label='lb')
                self.ui.refRefDataPlotwidget.canvas.ax.errorbar(data[:,0],data[:,6],fmt='-',label='rb')
            else: 
                self.ui.refRefDataPlotwidget.canvas.ax.errorbar(data[:,0],data[:,1],data[:,2],fmt='o',label='#0')
        else:
            print 'no current ref data to plot'
     #plot selected ref files
        if  len(self.selectedreffiles_rows)!=0:
            for i in range(len(self.selectedreffiles_rows)):
                data1=np.loadtxt(str(self.reffiles[self.selectedreffiles_rows[i]]), comments='#')
                if self.ui.refRRFCheckBox.checkState()!=0:
                    self.ui.refRefDataPlotwidget.canvas.ax.set_ylabel('R/RF')
                    if str(self.reffiles[self.selectedreffiles_rows[i]]).split('_')[-1]!='rrf.txt':
                        data1[:,0]=data1[:,0]+self.refQoff
                        data1[:,1]=data1[:,1]/self.fresnel(data1[:,0],self.refQc)
                        data1[:,2]=data1[:,2]/self.fresnel(data1[:,0],self.refQc)
                if self.ui.refQzSqrCheckBox.checkState()!=0:
                    self.ui.refRefDataPlotwidget.canvas.ax.set_xlabel('Qz^2')
                    data1[:,0]=data1[:,0]**2
                self.ui.refRefDataPlotwidget.canvas.ax.errorbar(data1[:,0],data1[:,1],data1[:,2],fmt='o',label='#'+str(self.selectedreffiles_rows[i]+1))
        # done
        if self.ui.refLogYCheckBox.checkState()!=0:
            self.ui.refRefDataPlotwidget.canvas.ax.set_yscale('log')
        else:
            self.ui.refRefDataPlotwidget.canvas.ax.set_yscale('linear')
        if self.ui.refLegendCheckBox.checkState()!=0:
            self.ui.refRefDataPlotwidget.canvas.ax.legend(loc=self.ui.refLegendLocComboBox.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        try:
            self.ui.refRefDataPlotwidget.canvas.ax.set_title('File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1])
        except:
            pass
        self.ui.refRefDataPlotwidget.canvas.draw()
        self.saveTemData()
        
    def saveRefData(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Reflectivity',directory=self.directory))
        if self.det=='Bruker':
            fid=open(self.saveFileName+'_ref.txt','w')
            for i in range(len(self.refData)):
                fid.write(str(self.normRefData[i,0])+'\t'+str(self.normRefData[i,1])+'\t'+str(self.normRefData[i,2])+'\t'+self.refInfoNew[i]+'\n')
            fid.close()
            try:
                fid=open(self.saveFileName+'_rrf.txt','w')
                for i in range(len(self.RrfData)):
                    fid.write(str(self.RrfData[i,0])+'\t'+str(self.RrfData[i,1])+'\t'+str(self.RrfData[i,2])+'\t'+self.refInfoNew[i]+'\n')
                fid.close()
            except:
                print 'No rrf data available!'
        elif self.det=='Pilatus':
            reffile=self.saveFileName+'_ref.txt'
            rrffile=self.saveFileName+'_rrf.txt'
            np.savetxt(reffile,self.normRefData,fmt='%.4e\t%.4e\t%.4e')
            try:
                np.savetxt(rrffile,self.RrfData,fmt='%.4e\t%.4e\t%.4e')
            except:
                print 'No rrf data available!'
    
    def saveTemData(self):
        fid=open(self.directory+'/temp_ref.txt','w')
        for i in range(len(self.refData)):
            if self.refpatchindex==0:
                fid.write(str(self.normRefData[i,0])+'\t'+str(self.normRefData[i,1])+'\t'+str(self.normRefData[i,2])+'\t'+self.refInfoNew[i]+'\n')
        fid.close()
    
    def clearCutGraph(self):
        self.ui.cutPlotMplWidget.canvas.ax.clear()
        self.ui.cutPlotMplWidget.canvas.draw()
        
    def updateCutData(self):
        if self.det=='Bruker':
            self.updateCcdCutData()
        elif self.det=='Pilatus':
           # print self.pilGIDshow, self.pilGISAXSshow
            if self.pilGIDshow==0 and self.pilGISAXSshow==0:
                self.updatePilCutData()
            else:
                self.updatePilGIDCutData()
        
    def updateCcdCutData(self):        
        self.selCutCcdData={}
        self.selCutCcdErrorData={}
        self.cutData={}
        self.cutLabel={}
        ini=float(str(self.ui.gixIntRangeLineEdit.text()).split(':')[0])
        fin=float(str(self.ui.gixIntRangeLineEdit.text()).split(':')[1])
        for i in self.selectedCcdFramesNums:
            self.selCutCcdData[i]=self.ccdData[i]#np.where(self.ccdData[i]<0,0,self.ccdData[i])#*self.absfac**self.ccd_AbsNum[i]/self.ccdMonc[i]
            self.selCutCcdErrorData[i]=self.ccdErrorData[i]#np.sqrt(self.ccdErrorData[i]**2/self.ccdMonc[i]**2+self.ccdData[i]**2/self.ccdMonc[i]**3)*self.absfac**self.ccd_AbsNum[i]
        if self.ui.gixSumCheckBox.checkState()!=0:
            self.selCutCcdData[-1]=self.selCutCcdData[self.selectedCcdFramesNums[0]]
            self.selCutCcdErrorData[-1]=self.selCutCcdErrorData[self.selectedCcdFramesNums[0]]**2
            for i in self.selectedCcdFramesNums[1:]:
                self.selCutCcdData[-1]=self.selCutCcdData[-1]+self.selCutCcdData[i]
                self.selCutCcdErrorData[-1]=self.selCutCcdErrorData[-1]+self.selCutCcdErrorData[i]**2
            self.selCutCcdData[-1]=self.selCutCcdData[-1]/len(self.selectedCcdFramesNums)
            self.selCutCcdErrorData[-1]=np.sqrt(self.selCutCcdErrorData[-1])/len(self.selectedCcdFramesNums)
            if self.ui.gixCutDirComboBox.currentText()=='H Cut':
                self.bruker.plotVint(self.selCutCcdData[-1],self.selCutCcdErrorData[-1],absfac=self.absfac,absnum=self.ccdSelected_AbsNum[0],cen=[self.xcenter[0]-self.xoff[0],self.ycenter[0]-self.yoff[0]],  hroi=None,  vroi=[max(int(ini)-1,0), min(int(fin)-1,1024)], ax_type=str(self.ui.gixAxesComboBox.currentText()), wavelength=self.wavelength[0],s2d_dist=self.distance[0],sh=self.ccdSelected_Sh[0],alpha=self.alpha[0], mon=None)
                self.cutData[-1]=self.bruker.vintData
            elif self.ui.gixCutDirComboBox.currentText()=='Qz Cut':
                ini=-int((self.distance[0]*np.arcsin(ini*self.wavelength[0]/2.0/np.pi-np.sin(self.alpha[0]))+self.ccdSelected_Sh[0])/0.06)+self.ycenter[0]-self.yoff[0]
                fin=-int((self.distance[0]*np.arcsin(fin*self.wavelength[0]/2.0/np.pi-np.sin(self.alpha[0]))+self.ccdSelected_Sh[0])/0.06)+self.ycenter[0]-self.yoff[0]                
               # print ini, fin
                self.bruker.plotVint(self.selCutCcdData[-1],self.selCutCcdErrorData[-1],absfac=self.absfac,absnum=self.ccdSelected_AbsNum[0], cen=[self.xcenter[0]-self.xoff[0],self.ycenter[0]-self.yoff[0]],  hroi=None,  vroi=[max(fin-1,0), min(ini-1,1024)], ax_type='Q', wavelength=self.wavelength[0],s2d_dist=self.distance[0],sh=self.ccdSelected_Sh[0],alpha=self.alpha[0], mon=None)
                self.cutData[-1]=self.bruker.vintData
            elif self.ui.gixCutDirComboBox.currentText()=='V Cut':
                self.bruker.plotHint(self.selCutCcdData[-1],self.selCutCcdErrorData[-1],absfac=self.absfac,absnum=self.ccdSelected_AbsNum[0],cen=[self.xcenter[0]-self.xoff[0],self.ycenter[0]],  hroi=[max(int(ini)-1,0), min(int(fin)-1,1024)],  vroi=None, ax_type=str(self.ui.gixAxesComboBox.currentText()), wavelength=self.wavelength[0],s2d_dist=self.distance[0],sh=self.ccdSelected_Sh[0],alpha=self.alpha[0], mon=None)
                self.cutData[-1]=self.bruker.hintData
            else:
                ini=int(2.0*self.distance[0]*np.arcsin(ini*self.wavelength[0]/4.0/np.pi)/0.06)+self.xcenter[0]-self.xoff[0]
                fin=int(2.0*self.distance[0]*np.arcsin(fin*self.wavelength[0]/4.0/np.pi)/0.06)+self.xcenter[0]-self.xoff[0]
                self.bruker.plotHint(self.selCutCcdData[-1],self.selCutCcdErrorData[-1],absfac=self.absfac,absnum=self.ccdSelected_AbsNum[0],cen=[self.xcenter[0]-self.xoff[0],self.ycenter[0]-self.yoff[0]],  hroi=[max(int(ini)-1,0), min(int(fin)-1,1024)],  vroi=None, ax_type='Q', wavelength=self.wavelength[0],s2d_dist=self.distance[0],sh=self.ccdSelected_Sh[0],alpha=self.alpha[0], mon=None)
                self.cutData[-1]=self.bruker.hintData
            self.cutLabel[-1]='summed cuts'
        else:
            j=0
            for i in self.selectedCcdFramesNums:
                if self.ui.gixCutDirComboBox.currentText()=='H Cut':
                    self.bruker.plotVint(self.selCutCcdData[i],self.selCutCcdErrorData[i],absfac=self.absfac,absnum=self.ccdSelected_AbsNum[j], cen=[self.xcenter[j]-self.xoff[j],self.ycenter[j]-self.yoff[j]], hroi=None, vroi=[int(ini)-1, int(fin)-1], ax_type=str(self.ui.gixAxesComboBox.currentText()), wavelength=self.wavelength[j],s2d_dist=self.distance[j],sh=self.ccdSelected_Sh[j],alpha=self.alpha[j], mon=None)
                    self.cutData[i]=self.bruker.vintData
                elif self.ui.gixCutDirComboBox.currentText()=='Qz Cut':
                    ini1=-int((self.distance[j]*np.arcsin(ini*self.wavelength[j]/2.0/np.pi-np.sin(self.alpha[j]))+self.ccdSelected_Sh[j])/0.06)+self.ycenter[j]
                    fin1=-int((self.distance[j]*np.arcsin(fin*self.wavelength[j]/2.0/np.pi-np.sin(self.alpha[j]))+self.ccdSelected_Sh[j])/0.06)+self.ycenter[j]                
                    self.bruker.plotVint(self.selCutCcdData[i],self.selCutCcdErrorData[i],absfac=self.absfac,absnum=self.ccdSelected_AbsNum[j],cen=[self.xcenter[j]-self.xoff[j],self.ycenter[j]-self.yoff[j]],  hroi=None,  vroi=[fin1-1, ini1-1], ax_type='Q', wavelength=self.wavelength[j],s2d_dist=self.distance[j],sh=self.ccdSelected_Sh[j],alpha=self.alpha[j], mon=None)
                    self.cutData[i]=self.bruker.vintData    
                elif self.ui.gixCutDirComboBox.currentText()=='V Cut':
                    self.bruker.plotHint(self.selCutCcdData[i],self.selCutCcdErrorData[i],absfac=self.absfac,absnum=self.ccdSelected_AbsNum[j],cen=[self.xcenter[j]-self.xoff[j],self.ycenter[j]-self.yoff[j]], hroi=[int(ini)-1, int(fin)-1],  vroi=None, ax_type=str(self.ui.gixAxesComboBox.currentText()), wavelength=self.wavelength[j],s2d_dist=self.distance[j],sh=self.ccdSelected_Sh[j],alpha=self.alpha[j], mon=None)
                    self.cutData[i]=self.bruker.hintData
                else:
                    ini1=int(2.0*self.distance[j]*np.arcsin(ini*self.wavelength[j]/4.0/np.pi)/0.06)+self.xcenter[j]
                    fin1=int(2.0*self.distance[j]*np.arcsin(fin*self.wavelength[j]/4.0/np.pi)/0.06)+self.xcenter[j]
                    self.bruker.plotHint(self.selCutCcdData[i],self.selCutCcdErrorData[i],absfac=self.absfac,absnum=self.ccdSelected_AbsNum[j], cen=[self.xcenter[j]-self.xoff[j],self.ycenter[j]-self.yoff[j]], hroi=[ini1-1, fin1-1],  vroi=None, ax_type='Q', wavelength=self.wavelength[j],s2d_dist=self.distance[j],sh=self.ccdSelected_Sh[j],alpha=self.alpha[j], mon=None)
                    self.cutData[i]=self.bruker.hintData
                self.cutLabel[i]=str(self.ui.imageListWidget.item(i).text().split('\t')[0])+' '+str(self.ui.imageListWidget.item(i).text().split('\t')[1])
                j=j+1
        if self.ui.gixSumCheckBox.checkState()!=0:
            if self.ui.cutErrorbarCheckBox.checkState()!=0:
                self.ui.cutPlotMplWidget.canvas.ax.errorbar(self.cutData[-1][:,0],self.cutData[-1][:,1],self.cutData[-1][:,2],fmt='o-',label=self.cutLabel[-1])
            else:
                self.ui.cutPlotMplWidget.canvas.ax.plot(self.cutData[-1][:,0],self.cutData[-1][:,1],'o',label=self.cutLabel[-1])
                
        else:
            fact=1
            if self.ui.cutErrorbarCheckBox.checkState()!=0:
                for i in self.selectedCcdFramesNums:
                    self.ui.cutPlotMplWidget.canvas.ax.errorbar(self.cutData[i][:,0],fact*self.cutData[i][:,1],fact*self.cutData[i][:,2],fmt='o-',label=self.cutLabel[i])
                    if self.ui.cutOffsetCheckBox.checkState()!=0:
                        fact=fact*float(self.ui.cutOffsetLineEdit.text())
            else:      
                for i in self.selectedCcdFramesNums:
                    self.ui.cutPlotMplWidget.canvas.ax.plot(self.cutData[i][:,0],fact*self.cutData[i][:,1],'o',label=self.cutLabel[i])
                    if self.ui.cutOffsetCheckBox.checkState()!=0:
                        fact=fact*float(self.ui.cutOffsetLineEdit.text())
        self.updateCutPlotData()
        
    def updatePilCutData(self):  
        self.ui.cutPlotMplWidget.canvas.ax.clear()
        self.selCutPilData={}
        self.selCutPilErrorData={}
        self.cutData={}
        self.cutLabel={}
        ini=float(str(self.ui.pilIntRangeLineEdit.text()).split(':')[0])
        fin=float(str(self.ui.pilIntRangeLineEdit.text()).split(':')[1])
        for i in self.selectedPilFramesNums:
            self.selCutPilData[i]=self.pilData[i]#np.where(self.ccdData[i]<0,0,self.ccdData[i])#*self.absfac**self.ccd_AbsNum[i]/self.ccdMonc[i]
            self.selCutPilErrorData[i]=self.pilErrorData[i]#np.sqrt(self.ccdErrorData[i]**2/self.ccdMonc[i]**2+self.ccdData[i]**2/self.ccdMonc[i]**3)*self.absfac**self.ccd_AbsNum[i]  
        if self.ui.gixSumCheckBox.checkState()!=0:
            self.selCutPilData[-1]=self.selCutPilData[self.selectedPilFramesNums[0]]
            self.selCutPilErrorData[-1]=self.selCutPilErrorData[self.selectedPilFramesNums[0]]**2
            for i in self.selectedPilFramesNums[1:]:
                self.selCutPilData[-1]=self.selCutPilData[-1]+self.selCutPilData[i]
                self.selCutPilErrorData[-1]=self.selCutPilErrorData[-1]+self.selCutPilErrorData[i]**2
            self.selCutPilData[-1]=self.selCutPilData[-1]/len(self.selectedPilFramesNums)
            self.selCutPilErrorData[-1]=np.sqrt(self.selCutPilErrorData[-1])/len(self.selectedPilFramesNums)
            if self.ui.pilCutDirComboBox.currentText()=='H Cut':  #vertical integration; provide PixY 
                self.pilatus.plotVint(self.selCutPilData[-1],self.selCutPilErrorData[-1],absfac=self.absfac,absnum=self.pilSelected_AbsNum[0],cen=[self.xcenter[0],self.ycenter[0]],  hroi=None,  vroi=[int(ini)-1, int(fin)-1], ax_type=str(self.ui.pilAxesComboBox.currentText()), wavelength=self.wavelength[0],s2d_dist=self.distance[0],sh=self.pilSelected_Sh[0],alpha=self.alpha[0], dth=self.dth[0],mon=None)
                self.cutData[-1]=self.pilatus.vintData
            elif self.ui.pilCutDirComboBox.currentText()=='Qz Cut':
                beta=np.arcsin(2*np.sin(self.alpha[0])-np.sin(self.truealpha[0]))  # get the beta value 
                ini=int(np.tan(np.arcsin(ini*self.wavelength[0]/2.0/np.pi-np.sin(self.alpha[0]))-beta)*self.pilSelected_Dist[0]/0.172)+self.ycenter[0]
                fin=int(np.tan(np.arcsin(fin*self.wavelength[0]/2.0/np.pi-np.sin(self.alpha[0]))-beta)*self.pilSelected_Dist[0]/0.172)+self.ycenter[0]
                #ini=-int((self.distance[0]*np.arcsin(ini*self.wavelength[0]/2.0/np.pi-np.sin(self.alpha[0]))+self.pilSelected_Sh[0])/0.172)+self.ycenter[0]
                #fin=-int((self.distance[0]*np.arcsin(fin*self.wavelength[0]/2.0/np.pi-np.sin(self.alpha[0]))+self.pilSelected_Sh[0])/0.172)+self.ycenter[0]  
                self.pilatus.plotVint(self.selCutPilData[-1],self.selCutPilErrorData[-1],absfac=self.absfac,absnum=self.pilSelected_AbsNum[0], cen=[self.xcenter[0],self.ycenter[0]],  hroi=None,  vroi=[int(ini)-1, int(fin)-1], ax_type='Q', wavelength=self.wavelength[0],s2d_dist=self.distance[0],sh=self.pilSelected_Sh[0],alpha=self.alpha[0], dth=self.dth[0], mon=None)
                self.cutData[-1]=self.pilatus.vintData
            elif self.ui.pilCutDirComboBox.currentText()=='V Cut':
                self.pilatus.plotHint(self.selCutPilData[-1],self.selCutPilErrorData[-1],absfac=self.absfac,absnum=self.pilSelected_AbsNum[0],cen=[self.xcenter[0],self.ycenter[0]],  hroi=[int(ini)-1, int(fin)-1],  vroi=None, ax_type=str(self.ui.pilAxesComboBox.currentText()), wavelength=self.wavelength[0],s2d_dist=self.distance[0],sh=self.pilSelected_Sh[0],alpha=self.alpha[0], truealpha=self.truealpha[0], mon=None)
                self.cutData[-1]=self.pilatus.hintData
            else:
                ini=int((2.0*np.arcsin(ini*self.wavelength[0]/4.0/np.pi)-self.dth[0])*self.distance[0]/0.172)+self.xcenter[0]
                fin=int((2.0*np.arcsin(fin*self.wavelength[0]/4.0/np.pi)-self.dth[0])*self.distance[0]/0.172)+self.xcenter[0]
                self.pilatus.plotHint(self.selCutPilData[-1],self.selCutPilErrorData[-1],absfac=self.absfac,absnum=self.pilSelected_AbsNum[0],cen=[self.xcenter[0],self.ycenter[0]],  hroi=[int(ini)-1, int(fin)-1],  vroi=None, ax_type='Q', wavelength=self.wavelength[0],s2d_dist=self.pilSelected_Dist[0],sh=self.pilSelected_Sh[0],alpha=self.alpha[0],truealpha=self.truealpha[0], mon=None)
                self.cutData[-1]=self.pilatus.hintData
            self.cutLabel[-1]='summed cuts'
        else:
            j=0
            for i in self.selectedPilFramesNums:
                if self.ui.pilCutDirComboBox.currentText()=='H Cut':
                    self.pilatus.plotVint(self.selCutPilData[i],self.selCutPilErrorData[i],absfac=self.absfac,absnum=self.pilSelected_AbsNum[j], cen=[self.xcenter[j],self.ycenter[j]], hroi=None, vroi=[int(ini)-1, int(fin)-1], ax_type=str(self.ui.pilAxesComboBox.currentText()), wavelength=self.wavelength[j],s2d_dist=self.distance[j],sh=self.pilSelected_Sh[j],alpha=self.alpha[j],dth=self.dth[j], mon=None)
                    self.cutData[i]=self.pilatus.vintData
                elif self.ui.pilCutDirComboBox.currentText()=='Qz Cut':
                    beta=np.arcsin(2*np.sin(self.alpha[j])-np.sin(self.truealpha[j]))  # get the beta value 
                    ini1=int(np.tan(np.arcsin(ini*self.wavelength[j]/2.0/np.pi-np.sin(self.alpha[j]))-beta)*self.distance[j]/0.172)+self.ycenter[j]
                    fin1=int(np.tan(np.arcsin(fin*self.wavelength[j]/2.0/np.pi-np.sin(self.alpha[j]))-beta)*self.distance[j]/0.172)+self.ycenter[j]
                    #ini1=-int((self.distance[j]*np.arcsin(ini*self.wavelength[j]/2.0/np.pi-np.sin(self.alpha[j]))+self.pilSelected_Sh[j])/0.172)+self.ycenter[j]
                    #fin1=-int((self.distance[j]*np.arcsin(fin*self.wavelength[j]/2.0/np.pi-np.sin(self.alpha[j]))+self.pilSelected_Sh[j])/0.172)+self.ycenter[j]                
                    self.pilatus.plotVint(self.selCutPilData[i],self.selCutPilErrorData[i],absfac=self.absfac,absnum=self.pilSelected_AbsNum[j],cen=[self.xcenter[j],self.ycenter[j]],  hroi=None,  vroi=[int(ini1)-1, int(fin1)-1], ax_type='Q', wavelength=self.wavelength[j],s2d_dist=self.distance[j],sh=self.pilSelected_Sh[j],alpha=self.alpha[j], dth=self.dth[j], mon=None)
                    self.cutData[i]=self.pilatus.vintData    
                elif self.ui.pilCutDirComboBox.currentText()=='V Cut':
                    self.pilatus.plotHint(self.selCutPilData[i],self.selCutPilErrorData[i],absfac=self.absfac,absnum=self.pilSelected_AbsNum[j],cen=[self.xcenter[j],self.ycenter[j]], hroi=[int(ini)-1, int(fin)-1],  vroi=None, ax_type=str(self.ui.pilAxesComboBox.currentText()), wavelength=self.wavelength[j],s2d_dist=self.distance[j],sh=self.pilSelected_Sh[j],alpha=self.alpha[j],truealpha=self.truealpha[j], mon=None)
                    self.cutData[i]=self.pilatus.hintData
                else:
                    ini1=int((2.0*np.arcsin(ini*self.wavelength[j]/4.0/np.pi)-self.dth[j])*self.distance[j]/0.172)+self.xcenter[j]
                    fin1=int((2.0*np.arcsin(fin*self.wavelength[j]/4.0/np.pi)-self.dth[j])*self.distance[j]/0.172)+self.xcenter[j]
                   # ini1=int(2.0*self.distance[j]*np.arcsin(ini*self.wavelength[j]/4.0/np.pi)/0.172)+self.xcenter[j]
                    #fin1=int(2.0*self.distance[j]*np.arcsin(fin*self.wavelength[j]/4.0/np.pi)/0.172)+self.xcenter[j]
                    self.pilatus.plotHint(self.selCutPilData[i],self.selCutPilErrorData[i],absfac=self.absfac,absnum=self.pilSelected_AbsNum[j], cen=[self.xcenter[j],self.ycenter[j]], hroi=[int(ini1-1), int(fin1-1)],  vroi=None, ax_type='Q', wavelength=self.wavelength[j],s2d_dist=self.distance[j],sh=self.pilSelected_Sh[j],alpha=self.alpha[j], truealpha=self.truealpha[j], mon=None)
                    self.cutData[i]=self.pilatus.hintData
                self.cutLabel[i]=str(self.ui.imageListWidget.item(i).text().split('\t')[0])+' '+str(self.ui.imageListWidget.item(i).text().split('\t')[1])
                j=j+1
        if self.ui.gixSumCheckBox.checkState()!=0:
            if self.ui.cutErrorbarCheckBox.checkState()!=0:
                self.ui.cutPlotMplWidget.canvas.ax.errorbar(self.cutData[-1][:,0],self.cutData[-1][:,1],self.cutData[-1][:,2],fmt='o-',label=self.cutLabel[-1])
            else:
                self.ui.cutPlotMplWidget.canvas.ax.plot(self.cutData[-1][:,0],self.cutData[-1][:,1],'o',label=self.cutLabel[-1])
                
        else:
            fact=1
            if self.ui.cutErrorbarCheckBox.checkState()!=0:
                for i in self.selectedPilFramesNums:
                    self.ui.cutPlotMplWidget.canvas.ax.errorbar(self.cutData[i][:,0],fact*self.cutData[i][:,1],fact*self.cutData[i][:,2],fmt='-',label=self.cutLabel[i])
                    if self.ui.cutOffsetCheckBox.checkState()!=0:
                        fact=fact*float(self.ui.cutOffsetLineEdit.text())
            else:      
                for i in self.selectedPilFramesNums:
                    self.ui.cutPlotMplWidget.canvas.ax.plot(self.cutData[i][:,0],fact*self.cutData[i][:,1],'-',label=self.cutLabel[i])
                    if self.ui.cutOffsetCheckBox.checkState()!=0:
                        fact=fact*float(self.ui.cutOffsetLineEdit.text())
        self.updateCutPlotData()
    
    def updatePilGIDCutData(self):
        self.ui.cutPlotMplWidget.canvas.ax.clear()
        start=float(self.ui.pilIntRangeLineEdit.text().split(':')[0])
        end=float(self.ui.pilIntRangeLineEdit.text().split(':')[1])
        if str(self.ui.pilAxesComboBox.currentText())=='Angles':   
            if str(self.ui.pilCutDirComboBox.currentText())=='H Cut':
                ini=np.where(self.pilGIDYAxs[:,0]>=start)[0][0]
                fin=np.where(self.pilGIDYAxs[:,0]<=end)[0][-1]
                xaxis=self.pilGIDXAxs[0,:]
                self.cutData=np.vstack((xaxis,np.sum(self.pilGIDData[ini:fin,:],axis=0),np.sqrt(np.sum(self.pilGIDDataErr[ini:fin,:]**2,axis=0)))).transpose()
            elif str(self.ui.pilCutDirComboBox.currentText())=='V Cut':
                ini=np.where(self.pilGIDXAxs[0,:]>=start)[0][0]
                fin=np.where(self.pilGIDXAxs[0,:]<=end)[0][-1]
                xaxis=self.pilGIDYAxs[:,0]
                self.cutData=np.vstack((xaxis,np.sum(self.pilGIDData[:,ini:fin],axis=1),np.sqrt(np.sum(self.pilGIDDataErr[:,ini:fin]**2,axis=1)))).transpose()
        elif str(self.ui.pilAxesComboBox.currentText())=='Q':
            if self.pilGISAXSshow==0:
                if str(self.ui.pilCutDirComboBox.currentText())=='H Cut':
                    ini=np.where(self.pilGIDYAxs_Q[:,0]>=start)[0][0]
                    fin=np.where(self.pilGIDYAxs_Q[:,0]<=end)[0][-1]
                    xaxis=self.pilGIDXAxs_Q[0,:]
                    self.cutData=np.vstack((xaxis,np.sum(self.pilGIDData_Q[ini:fin,:],axis=0),np.sqrt(np.sum(self.pilGIDDataErr_Q[ini:fin,:]**2,axis=0)))).transpose()
                elif str(self.ui.pilCutDirComboBox.currentText())=='V Cut':
                    ini=np.where(self.pilGIDXAxs_Q[0,:]>=start)[0][0]
                    fin=np.where(self.pilGIDXAxs_Q[0,:]<=end)[0][-1]
                    xaxis=self.pilGIDYAxs_Q[:,0]
                    self.cutData=np.vstack((xaxis,np.sum(self.pilGIDData_Q[:,ini:fin],axis=1),np.sqrt(np.sum(self.pilGIDDataErr_Q[:,ini:fin]**2,axis=1)))).transpose()
            else:
               # print 'I am here'
                if str(self.ui.pilCutDirComboBox.currentText())=='H Cut':
                    ini=np.where(self.pilYbin[:,0]>=start)[0][0]
                    fin=np.where(self.pilYbin[:,0]<=end)[0][-1]
                   # print ini, fin
                    xaxis=self.pilXbin[0,:]
                    self.cutData=np.vstack((xaxis,np.sum(self.pilDataBin[ini:fin,:],axis=0),np.sqrt(np.sum(self.pilErrorDataBin[ini:fin,:]**2,axis=0)))).transpose()
                elif str(self.ui.pilCutDirComboBox.currentText())=='V Cut':
                    ini=np.where(self.pilXbin[0,:]>=start)[0][0]
                    fin=np.where(self.pilXbin[0,:]<=end)[0][-1]
                    xaxis=self.pilYbin[:,0]
                    self.cutData=np.vstack((xaxis,np.sum(self.pilDataBin[:,ini:fin],axis=1),np.sqrt(np.sum(self.pilErrorDataBin[:,ini:fin]**2,axis=1)))).transpose()
        if self.ui.cutErrorbarCheckBox.checkState()!=0:
            self.ui.cutPlotMplWidget.canvas.ax.errorbar(self.cutData[:,0],self.cutData[:,1],self.cutData[:,2],fmt='b-',label='GID Cut')
        else:
            self.ui.cutPlotMplWidget.canvas.ax.plot(self.cutData[:,0],self.cutData[:,1],'b-',label='GID Cut')
        self.updateCutPlotData()
    
    def updateCutPlotData(self):
        self.ui.PlotWidget.setCurrentIndex(4) 
        self.cutLogX=self.ui.cutLogXCheckBox.checkState()
        self.cutLogY=self.ui.cutLogYCheckBox.checkState()
        self.cutGrid=self.ui.cutGridCheckBox.checkState()
        self.ui.cutPlotMplWidget.canvas.ax.set_xscale('linear')
        self.ui.cutPlotMplWidget.canvas.ax.set_yscale('linear')
        self.ui.cutPlotMplWidget.canvas.ax.set_ylabel('Intensity')
        self.ui.cutPlotMplWidget.canvas.ax.grid(b=False)
        title='File: '+self.specFileName+' S# '+str([item for item in np.sort(self.selectedScanNums)])[1:-1]
        if self.det=='Pilatus':
            currentcutindex=self.ui.pilCutDirComboBox.currentText()
            cutrange=str(self.ui.pilIntRangeLineEdit.text())
        elif self.det=='Bruker':
            currentcutindex=self.ui.gixCutDirComboBox.currentText()
            cutrange=str(self.ui.gixIntRangeLineEdit.text())
        if currentcutindex=='H Cut':
            self.ui.cutPlotMplWidget.canvas.ax.set_title(title+'\n'+'H Cut     Vint Range['+cutrange+']', fontsize=16)
            self.ui.cutPlotMplWidget.canvas.ax.set_xlabel('Pixels', fontsize=16)
        elif currentcutindex=='V Cut':
            self.ui.cutPlotMplWidget.canvas.ax.set_title(title+'\n'+'V Cut     Hint Range['+cutrange+']', fontsize=16)
            self.ui.cutPlotMplWidget.canvas.ax.set_xlabel('Pixels', fontsize=16)
        elif currentcutindex=='Qz Cut':
            self.ui.cutPlotMplWidget.canvas.ax.set_title(title+'\n'+'Qz Cut     Qz Range['+cutrange+']', fontsize=16)
            self.ui.cutPlotMplWidget.canvas.ax.set_xlabel(r'$Q_{xy}$'+' '+r'$[\AA^{-1}]$', fontsize=16)
        else:
            self.ui.cutPlotMplWidget.canvas.ax.set_title(title+'\n'+'Qxy Cut     Qxy Range['+cutrange+']', fontsize=16)
            self.ui.cutPlotMplWidget.canvas.ax.set_xlabel(r'$Q_z$'+' '+r'$[\AA^{-1}]$', fontsize=16)
        if self.cutLogX!=0:
            self.ui.cutPlotMplWidget.canvas.ax.set_xscale('log')
        if self.cutLogY!=0:
            self.ui.cutPlotMplWidget.canvas.ax.set_yscale('log')
        if self.cutGrid!=0:
            self.ui.cutPlotMplWidget.canvas.ax.grid(b=True,color='r',linestyle='--')
        if self.ui.cutLegendCheckBox.checkState()!=0:
            self.ui.cutPlotMplWidget.canvas.ax.legend(loc=self.ui.cutLegendLocComboBox.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        self.ui.cutPlotMplWidget.canvas.draw()
        
    def saveCutData(self):
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Cuts', directory=self.directory))
        if self.det=='Bruker':            
            if self.ui.gixSumCheckBox.checkState()!=0:
                self.fname=self.saveFileName+str(self.ui.imageListWidget.item(0).text().split('\t')[0])+'_sumcut.txt'
                np.savetxt(self.fname,self.cutData[-1],fmt='%.4f\t%.4e\t%.4e')
            else:
                for i in self.selectedCcdFramesNums:
                    self.fname=self.saveFileName+str(self.ui.imageListWidget.item(i).text().split('\t')[0])+'_'+str(self.ui.imageListWidget.item(i).text().split('\t')[1])+'_cut.txt'
                    np.savetxt(self.fname,self.cutData[i],fmt='%.4f\t%.4e\t%.4e')
        elif self.det=='Pilatus':
            if self.pilGIDshow==0 and self.pilGISAXSshow==0:
                if self.ui.gixSumCheckBox.checkState()!=0:
                    self.fname=self.saveFileName+str(self.ui.imageListWidget.item(0).text().split('\t')[0])+'_sumcut.txt'
                    np.savetxt(self.fname,self.cutData[-1],fmt='%.4f\t%.4e\t%.4e')
                else:
                    for i in self.selectedPilFramesNums:
                        self.fname=self.saveFileName+str(self.ui.imageListWidget.item(i).text().split('\t')[0])+'_'+str(self.ui.imageListWidget.item(i).text().split('\t')[1])+'_cut.txt'
                        np.savetxt(self.fname,self.cutData[i],fmt='%.4f\t%.4e\t%.4e')
            else:
                self.fname=self.saveFileName
                np.savetxt(self.fname,self.cutData,fmt='%.4f\t%.4e\t%.4e')
        
            
            
    def displayScanInfo(self):
        Dialog=QDialog(self)
        ui=uic.loadUi('scanInfoDialog.ui',Dialog)
        ui.show()
        snumtit=''
        impsnumtit=''
        for i in self.selectedScanNums:
            try:
                for j in range(self.specRead.Data[i]['StartScanLineNum'],self.specRead.Data[i+1]['StartScanLineNum']):
                    snumtit=snumtit+self.specRead.SpecFileFull[j]
            except:
                for j in range(self.specRead.Data[i]['StartScanLineNum'],len(self.specRead.SpecFileFull)):
                    snumtit=snumtit+self.specRead.SpecFileFull[j]
            try:
                impsnumtit=impsnumtit+self.specData[i]['ScanLine']+'\n'
                impsnumtit=impsnumtit+'Detector:\t'+self.specPar[i]['Detector']+'\n\n'
                impsnumtit=impsnumtit+'Energy:\t'+str(format(12.39842/self.specPar[i]['Wavelength'],'.4f'))+' keV\n\n'
                impsnumtit=impsnumtit+'Wavelength:\t'+str(self.specPar[i]['Wavelength'])+' '+u'\u212b'+'\n\n'  
                impsnumtit=impsnumtit+'Absorber\t'+str(self.specPar[i]['Absorber'])+'\n\n'
                impsnumtit=impsnumtit+'SX\t'+str(self.specPar[i]['Sample_X'])+'\n\n'
                impsnumtit=impsnumtit+'Qx\tQy\tQz\n'
                impsnumtit=impsnumtit+str(self.specPar[i]['Q'][0])+'\t'+str(self.specPar[i]['Q'][1])+'\t'+str(self.specPar[i]['Q'][2])+'\n\n'
                impsnumtit=impsnumtit+'S1 [H, V]\tS4 [H, V]\tS5 [H]\n'
                impsnumtit=impsnumtit+'['+str(format(self.specPar[i]['S1L']+self.specPar[i]['S1R'],'.3f'))+', '+str(format(self.specPar[i]['S1T']+self.specPar[i]['S1B'],'.3f'))+']\t'+'['+str(format(self.specPar[i]['S4L']+self.specPar[i]['S4R'],'.1f'))+', '+str(format(self.specPar[i]['S4T']+self.specPar[i]['S4B'],'.1f'))+']\t'+'['+str(format(self.specPar[i]['S5L']+self.specPar[i]['S5R'],'.1f'))+']\n\n'
                impsnumtit=impsnumtit+'S2 [H, V]\tS3 [H, V]\tS7 [H, V]\n'
                impsnumtit=impsnumtit+'['+str(format(self.specPar[i]['S2L']+self.specPar[i]['S2R'],'.1f'))+', '+str(format(self.specPar[i]['S2T']+self.specPar[i]['S2B'],'.1f'))+']\t'+'['+str(format(self.specPar[i]['S3L']+self.specPar[i]['S3R'],'.1f'))+', '+str(format(self.specPar[i]['S3T']+self.specPar[i]['S3B'],'.1f'))+']\t'+'['+str(format(self.specPar[i]['S7L']+self.specPar[i]['S7R'],'.1f')+', '+str(format(self.specPar[i]['S7T']+self.specPar[i]['S7B'],'.1f')))+']\n\n'
                if self.specPar[i]['Detector']=='Pilatus':
                    impsnumtit=impsnumtit+'pd[x, y]:\t['+str(self.specPar[i]['DBPos'][0])+', '+str(self.specPar[i]['DBPos'][1])+']\n'
                elif self.specPar[i]['Detector']=='Bruker':
                    impsnumtit=impsnumtit+'ad[x, y]:\t['+str(self.specPar[i]['DBPos'][0])+', '+str(self.specPar[i]['DBPos'][1])+']\n'
                impsnumtit=impsnumtit+'-----------------------------------------------\n'
            except:
                impsnumtit='No detailed info found in SPEC file!!!'
        ui.scanInfoTextBrowser.append(snumtit)
        ui.impInfoTextBrowser.append(impsnumtit)
        
        Dialog.exec_()  
        
            
        
    def g_lfunction(self,par,x,y):
        return y-par[0]+par[1]*np.tan(x*np.pi/180)
        
    def calg_l2(self):
        Dialog=QDialog(self)
        self.uig_l2=uic.loadUi('g_l2Dialog.ui', Dialog)
        self.uig_l2.show()
        scannum=[item for item in self.selectedScanNums]
        #self.g_l2alpha=[self.specPar[i]['In_Rot'] for i in self.selectedScanNums]
        self.g_l2alpha = [round(np.arctan(round(self.specPar[i]['Q'][2],3)*self.specPar[i]['Wavelength']/4/np.pi)*180/np.pi,5) for i in self.selectedScanNums]
        self.g_l2center=[self.scanCenter[i] for i in self.selectedScanNums]
        self.g_l2=[self.specPar[i]['g_l2'] for i in self.selectedScanNums]
        self.uig_l2.g_l2ScanLineEdit.setText(str(scannum)[1:-1])
        self.uig_l2.g_l2CenterLineEdit.setText(str(self.g_l2center)[1:-1])
        self.connect(self.uig_l2.g_l2CalPushButton,SIGNAL('clicked()'), self.g_l2plot)
        Dialog.exec_()
        
    def g_l2plot(self):
        self.uig_l2.g_l2TextBrowser.clear()
        if self.uig_l2.g_l2centercheckBox.checkState()!=0:
            self.g_l2center=map(float, str(self.uig_l2.g_l2CenterLineEdit.text()).split(','))
        par=[0,self.specPar[self.selectedScanNums[0]]['g_l2']]
        p=leastsq(self.g_lfunction,par,args=(np.array(self.g_l2alpha),np.array(self.g_l2center)),maxfev=5000)
        fp=p[0]
        self.uig_l2.g_l2plotWidget.canvas.ax.clear()
        self.uig_l2.g_l2plotWidget.canvas.ax.plot(self.g_l2alpha,self.g_l2center,'ro')
        self.uig_l2.g_l2plotWidget.canvas.ax.plot(np.array(self.g_l2alpha),fp[0]-fp[1]*np.tan(np.array(self.g_l2alpha)*np.pi/180),'-b',label=('g_l2= %.3f'%fp[1]))
        self.uig_l2.g_l2plotWidget.canvas.ax.set_xlabel('Alpha [Degrees]')
        self.uig_l2.g_l2plotWidget.canvas.ax.set_ylabel('Sh_Center [mm]')
        self.uig_l2.g_l2plotWidget.canvas.ax.legend(frameon=False,scatterpoints=0,numpoints=1)
        self.uig_l2.g_l2plotWidget.canvas.draw()
        self.uig_l2.g_l2TextBrowser.append(':Old Values:')
        self.uig_l2.g_l2TextBrowser.append('-----------------------------------------')
        self.uig_l2.g_l2TextBrowser.append('Scan \t Alpha \t Sh Nom \t Sh Cen \t N-C' )
        j=0
        for i in self.selectedScanNums:
            self.uig_l2.g_l2TextBrowser.append(str(i+1)+' \t %.5f '%self.g_l2alpha[j]+'\t %.4f'%(-np.tan(self.g_l2alpha[j]*np.pi/180)*self.g_l2[j])+'\t %.4f '%self.g_l2center[j]+'\t %.4f '%(-np.tan(self.g_l2alpha[j]*np.pi/180)*self.g_l2[j]-self.g_l2center[j]))
            j=j+1
        self.uig_l2.g_l2TextBrowser.append('\n')
        self.uig_l2.g_l2TextBrowser.append(':New Values:')
        self.uig_l2.g_l2TextBrowser.append('-----------------------------------------')
        self.uig_l2.g_l2TextBrowser.append('Scan \t Alpha \t Sh Nom \t Sh Cen \t N-C' )
        j=0
        for i in self.selectedScanNums:
            self.uig_l2.g_l2TextBrowser.append(str(i+1)+' \t %.5f '%self.g_l2alpha[j]+'\t %.4f'%(fp[0]-np.tan(self.g_l2alpha[j]*np.pi/180)*fp[1])+'\t %.4f '%self.g_l2center[j]+'\t %.4f '%(-np.tan(self.g_l2alpha[j]*np.pi/180)*fp[1]+fp[0]-self.g_l2center[j]))
            j=j+1
        self.uig_l2.g_l2TextBrowser.append('\n')
        self.uig_l2.g_l2TextBrowser.append('old g_l2\t= %.3f'%self.g_l2[0])
        self.uig_l2.g_l2TextBrowser.append('new g_l2\t= %.3f'%fp[1])
        self.uig_l2.g_l2TextBrowser.append('offset\t= %.3f'%fp[0])
        
    def calg_l3(self):
        Dialog=QDialog(self)
        self.uig_l3=uic.loadUi('g_l3Dialog.ui', Dialog)
        self.uig_l3.show()
        scannum=[item for item in self.selectedScanNums]
        #self.g_l3alpha=[self.specPar[i]['In_Rot'] for i in self.selectedScanNums]
        self.g_l3alpha = [round(np.arctan(round(self.specPar[i]['Q'][2], 3) * self.specPar[i]['Wavelength'] / 4 / np.pi) * 180 / np.pi, 5) for i in self.selectedScanNums]
        self.g_l3center=[self.scanCenter[i] for i in self.selectedScanNums]
        self.g_l3=[self.specPar[i]['g_l3'] for i in self.selectedScanNums]
        self.g_l2=[self.specPar[i]['g_l2'] for i in self.selectedScanNums]
        self.uig_l3.g_l3ScanLineEdit.setText(str(scannum)[1:-1])
        self.uig_l3.g_l3CenterLineEdit.setText(str(self.g_l3center)[1:-1])
        self.connect(self.uig_l3.g_l3CalPushButton,SIGNAL('clicked()'), self.g_l3plot)
        Dialog.exec_()
        
    def g_l3plot(self):
        self.uig_l3.g_l3TextBrowser.clear()
        if self.uig_l3.g_l3centercheckBox.checkState()!=0:
            self.g_l3center=map(float, str(self.uig_l3.g_l3CenterLineEdit.text()).split(','))
        par=[0,self.specPar[self.selectedScanNums[0]]['g_l2']-self.specPar[self.selectedScanNums[0]]['g_l3']]
        p=leastsq(self.g_lfunction,par,args=(np.array(self.g_l3alpha),np.array(self.g_l3center)),maxfev=5000)
        fp=p[0]
        self.uig_l3.g_l3plotWidget.canvas.ax.clear()
        self.uig_l3.g_l3plotWidget.canvas.ax.plot(self.g_l3alpha,self.g_l3center,'ro')
        self.uig_l3.g_l3plotWidget.canvas.ax.plot(np.array(self.g_l3alpha),fp[0]-fp[1]*np.tan(np.array(self.g_l3alpha)*np.pi/180),'-b',label=('g_l2-g_l3= %.3f'%fp[1]))
        self.uig_l3.g_l3plotWidget.canvas.ax.set_xlabel('Alpha [Degrees]')
        self.uig_l3.g_l3plotWidget.canvas.ax.set_ylabel('Oh_Center [mm]')
        self.uig_l3.g_l3plotWidget.canvas.ax.legend(frameon=False,scatterpoints=0,numpoints=1)
        self.uig_l3.g_l3plotWidget.canvas.draw()
        self.uig_l3.g_l3TextBrowser.append(':Old Values:')
        self.uig_l3.g_l3TextBrowser.append('-----------------------------------------')
        self.uig_l3.g_l3TextBrowser.append('Scan \t Alpha \t Oh Nom \t Oh Cen \t N-C' )
        j=0
        for i in self.selectedScanNums:
            self.uig_l3.g_l3TextBrowser.append(str(i+1)+' \t %.5f '%self.g_l3alpha[j]+'\t %.4f'%(np.tan(self.g_l3alpha[j]*np.pi/180)*(self.g_l3[j]-self.g_l2[j]))+'\t %.4f '%self.g_l3center[j]+'\t %.4f '%(np.tan(self.g_l3alpha[j]*np.pi/180)*(self.g_l3[j]-self.g_l2[j])-self.g_l3center[j]))
            j=j+1
        self.uig_l3.g_l3TextBrowser.append('\n')
        self.uig_l3.g_l3TextBrowser.append(':New Values:')
        self.uig_l3.g_l3TextBrowser.append('-----------------------------------------')
        self.uig_l3.g_l3TextBrowser.append('Scan \t Alpha \t Oh Nom \t Oh Cen \t N-C' )
        j=0
        for i in self.selectedScanNums:
            self.uig_l3.g_l3TextBrowser.append(str(i+1)+' \t %.5f '%self.g_l3alpha[j]+'\t %.4f'%(fp[0]-np.tan(self.g_l3alpha[j]*np.pi/180)*fp[1])+'\t %.4f '%self.g_l3center[j]+'\t %.4f '%(-np.tan(self.g_l3alpha[j]*np.pi/180)*fp[1]+fp[0]-self.g_l3center[j]))
            j=j+1
        self.uig_l3.g_l3TextBrowser.append('\n')
        self.uig_l3.g_l3TextBrowser.append('old (g_l2-g_l3)\t= %.3f'%(self.g_l2[0]-self.g_l3[0]))
        self.uig_l3.g_l3TextBrowser.append('new (g_l2_g_l3)\t= %.3f'%fp[1])
        self.uig_l3.g_l3TextBrowser.append('offset\t\t= %.3f'%fp[0])
    
    def absratiofunction(self,par,x,y):
        return y-np.log(par[0])+x*np.log(par[1])+par[2]*par[0]/par[1]**x
        #return y-par[0]/np.power(par[1],x)
        
    def calabsrat(self):
        scannum=self.selectedScanNums[0]
        if self.specData[scannum]['ScanVar'][0]!='Absorber':
            self.messageBox('Error:: This scan is NOT an absorber scan!!')
        else:
            Dialog=QDialog(self)
            uiabs=uic.loadUi('absratioDialog.ui', Dialog)
            uiabs.show()
            x=self.specData[scannum]['Absorber']
            z=self.specData[scannum]['Monc']
            ave=sum(z)/float(len(z))
            y=np.log(self.specData[scannum][str(self.ui.spYComboBox.currentText())]/z*ave/self.specData[scannum]['Seconds'])
            par=[1.0e+10,1.5,1.0e-8]
            #par=[1.0e+10,1.6]
            p=leastsq(self.absratiofunction,par,args=(x,y),maxfev=5000)
            fp=p[0]
            #print fp, y
            uiabs.absRatplotWidget.canvas.ax.clear()
            uiabs.absRatplotWidget.canvas.ax.plot(x,np.exp(y),'ro')
            uiabs.absRatplotWidget.canvas.ax.plot(x,fp[0]/fp[1]**x*np.exp(-fp[2]*par[0]/fp[1]**x),'-b')
            uiabs.absRatplotWidget.canvas.ax.set_xlabel('Absorber Number')
            uiabs.absRatplotWidget.canvas.ax.set_ylabel('Intensity')
            uiabs.absRatplotWidget.canvas.ax.set_yscale('log')
            uiabs.absRatplotWidget.canvas.draw()
            uiabs.absRatScanLineEdit.setText(str(scannum+1))
            uiabs.absRatAbsLineEdit.setText(str(fp[1]))
            uiabs.absRatSouLineEdit.setText('%.3e'%fp[0])
            uiabs.absRatDeaLineEdit.setText('%.3e'%fp[2])
            Dialog.exec_()
        
    def updatePlotFile(self): #update plot files in the listwidget
        self.ui.plotFileListWidget.clear()
        for i in range(len(self.plotfiles)):
            try:  #for pc
                self.ui.plotFileListWidget.addItem('#'+str(i+1)+self.halftab+str(self.plotfiles[i].split('\\')[-2])+'\\'+str(self.plotfiles[i].split('\\')[-1]))
            except: #for mac and linux
                self.ui.plotFileListWidget.addItem('#'+str(i+1)+self.halftab+str(self.plotfiles[i].split('/')[-2])+'/'+str(self.plotfiles[i].split('/')[-1]))

    def addPlotFile(self): #add plot files into the listwidget and deselect all ref files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Files to import', directory=self.directory, filter='Files (*.*)')
        self.plotfiles=self.plotfiles+map(str, f)
        self.updatePlotFile()
        
    def updateSelectedPlotFile(self): #update the selected ref files in the listwidget
        selectedplotfiles=self.ui.plotFileListWidget.selectedItems()
        self.fitplotstatus=0
        self.selectedplotfiles_rows=[]
        for item in selectedplotfiles:
            self.selectedplotfiles_rows.append(self.ui.plotFileListWidget.row(item))
        self.selectedplotfiles_rows.sort()
        self.plotscale=[[1,0,1,0] for i in range(len(self.selectedplotfiles_rows))]
        self.updatePlotPlot()
        
    def removePlotFile(self): #remove plot files in the listwidget and deselect all ref files in the listwidget
        items=self.ui.plotFileListWidget.selectedItems()
        items.reverse()
        for item in items:
            self.plotfiles.pop(self.ui.plotFileListWidget.row(item))
        self.ui.plotFileListWidget.clear()
        self.updatePlotFile()
                
        
    def upPlotFile(self):
        rows=self.selectedplotfiles_rows
        for i in rows:
            self.plotfiles.insert(max(0,i-1),self.plotfiles.pop(i))
        self.updatePlotFile()
        for i in rows:
            self.ui.plotFileListWidget.setItemSelected(self.ui.plotFileListWidget.item(i-1),True)
        
    def downPlotFile(self):
        rows=self.selectedplotfiles_rows
        rows.reverse()
        for i in self.selectedplotfiles_rows:
            self.plotfiles.insert(min(len(self.plotfiles),i+1),self.plotfiles.pop(i))   
        self.updatePlotFile()
        for i in rows:
            self.ui.plotFileListWidget.setItemSelected(self.ui.plotFileListWidget.item(i+1),True)
        
    def updatePlotPlot(self): #update the plot in the plot plotwidget
        self.ui.PlotMplWidget.canvas.ax.clear()
        labelsize=self.ui.onedLabelSizeSpinBox.value()
        ticksize=self.ui.onedTickSizeSpinBox.value()
        if self.ui.onedStyComboBox.currentIndex()==0:
            style='o'
        elif self.ui.onedStyComboBox.currentIndex()==1:
            style='-'
        elif self.ui.onedStyComboBox.currentIndex()==2:
            style='o-'
        if  len(self.selectedplotfiles_rows)!=0: #plot plot files
            for i in range(len(self.selectedplotfiles_rows)):
                data1=np.loadtxt(str(self.plotfiles[self.selectedplotfiles_rows[i]]), comments='#')
                data1=data1[np.argsort(data1[:,0])]
                try:  #for pc
                    datalabel=str(self.plotfiles[self.selectedplotfiles_rows[i]].split('\\')[-1])
                except: #for mac and linux
                    datalabel=str(self.plotfiles[self.selectedplotfiles_rows[i]].split('/')[-1])
                try:
                    self.ui.PlotMplWidget.canvas.ax.errorbar(data1[:,0]*self.plotscale[i][0]+self.plotscale[i][1],data1[:,1]*self.plotscale[i][2]+self.plotscale[i][3],data1[:,2]*self.plotscale[i][2],fmt=style,label='#'+str(self.selectedplotfiles_rows[i]+1)+' '+datalabel)
                except:
                    self.ui.PlotMplWidget.canvas.ax.plot(data1[:,0]*self.plotscale[i][0]+self.plotscale[i][1],data1[:,1]*self.plotscale[i][2]+self.plotscale[i][3],'-',label='#'+str(self.selectedplotfiles_rows[i]+1)+' '+datalabel)
        if self.fitplotstatus!=0:  #add fit result
            fitdata=np.array(self.peakfitdata)
            fitbgdata=np.array(self.peakbgfitdata)
            self.ui.PlotMplWidget.canvas.ax.plot(fitdata[:,0],fitdata[:,1],'r-')
            self.ui.PlotMplWidget.canvas.ax.plot(fitbgdata[:,0],fitbgdata[:,1],'g-')
        if self.ui.plotLogXCheckBox.checkState()!=0:
            self.ui.PlotMplWidget.canvas.ax.set_xscale('log')
        if self.ui.plotLogYCheckBox.checkState()!=0:
            self.ui.PlotMplWidget.canvas.ax.set_yscale('log')
        if self.ui.plotGridCheckBox.checkState()!=0:
            self.ui.PlotMplWidget.canvas.ax.grid(b=True,color='r',linestyle='--')    
        if self.ui.plotLegendCheckBox.checkState()!=0:
            self.ui.PlotMplWidget.canvas.ax.legend(loc=self.ui.plotLegendLocComboBox.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        self.ui.PlotMplWidget.canvas.ax.set_xlabel('X', fontsize=labelsize)
        self.ui.PlotMplWidget.canvas.ax.set_ylabel('Y', fontsize=labelsize)
        self.ui.PlotMplWidget.canvas.ax.tick_params(labelsize=ticksize)
        self.ui.PlotMplWidget.canvas.draw()
        
    def graPlotPlot(self): #deriative data wiht x-axis having equal spacing
        self.ui.PlotMplWidget.canvas.ax.clear()
        labelsize=self.ui.onedLabelSizeSpinBox.value()
        ticksize=self.ui.onedTickSizeSpinBox.value()
        if self.ui.onedStyComboBox.currentIndex()==0:
            style='o'
        elif self.ui.onedStyComboBox.currentIndex()==1:
            style='-'
        elif self.ui.onedStyComboBox.currentIndex()==2:
            style='o-'
        if  len(self.selectedplotfiles_rows)!=0: #plot plot files
            for i in range(len(self.selectedplotfiles_rows)):
                data1=np.loadtxt(str(self.plotfiles[self.selectedplotfiles_rows[i]]), comments='#')
                data1=data1[np.argsort(data1[:,0])]
                x=np.linspace(data1[:,0][0],data1[:,0][-1],len(data1))
                fint=interp1d(data1[:,0], data1[:,1], kind='linear')
                y=fint(x)
                y=np.gradient(y)
                self.gradata=np.vstack((x,y)).T
                try:  #for pc
                    datalabel=str(self.plotfiles[self.selectedplotfiles_rows[i]].split('\\')[-1])
                except: #for mac and linux
                    datalabel=str(self.plotfiles[self.selectedplotfiles_rows[i]].split('/')[-1])
                self.ui.PlotMplWidget.canvas.ax.plot(self.gradata[:,0],self.gradata[:,1],style,label='#'+str(self.selectedplotfiles_rows[i]+1)+' '+datalabel)
        if self.ui.plotLogXCheckBox.checkState()!=0:
            self.ui.PlotMplWidget.canvas.ax.set_xscale('log')
        if self.ui.plotLogYCheckBox.checkState()!=0:
            self.ui.PlotMplWidget.canvas.ax.set_yscale('log')
        if self.ui.plotGridCheckBox.checkState()!=0:
            self.ui.PlotMplWidget.canvas.ax.grid(b=True,color='r',linestyle='--')    
        if self.ui.plotLegendCheckBox.checkState()!=0:
            self.ui.PlotMplWidget.canvas.ax.legend(loc=self.ui.plotLegendLocComboBox.currentIndex()+1,frameon=False,scatterpoints=0,numpoints=1)
        self.ui.PlotMplWidget.canvas.ax.set_xlabel('X', fontsize=labelsize)
        self.ui.PlotMplWidget.canvas.ax.set_ylabel('Y', fontsize=labelsize)
        self.ui.PlotMplWidget.canvas.ax.set_title('Derivative', fontsize=labelsize)
        self.ui.PlotMplWidget.canvas.ax.tick_params(labelsize=ticksize)
        self.ui.PlotMplWidget.canvas.draw()
        
    def saveGraPlotData(self): #save the derivative data
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Derivative Data',directory=self.directory))
        fitfile=self.saveFileName+'_der.txt'
        np.savetxt(fitfile,self.gradata,fmt='%.4e\t%.4e')
    
    def setPlotScale(self): #set the scale of each data in the plot
        if len(self.selectedplotfiles_rows)==0:
            self.messageBox('Warning:: No files selected!')
        else:
            row=len(self.selectedplotfiles_rows)
            Dialog=QDialog(self)
            self.uiplotscale=uic.loadUi('plotscale.ui', Dialog)
            self.uiplotscale.scaleTW.setRowCount(row) #set the table size; 4 column is fixed
            self.uiplotscale.show()
            self.uiplotscale.scaleLabel.setText('Plot Scale Setup: X=X*Factor+Offset')
            self.uiplotscale.scaleTW.setHorizontalHeaderLabels(QStringList()<<"X Factor"<<"X Offset"<<"Y Factor"<<"Y Offset") #set the horizontal header
            vlabel=QStringList() #set the vertical header 
            for i in range(row):
                vlabel.append("#"+str(self.selectedplotfiles_rows[i]+1))
            self.uiplotscale.scaleTW.setVerticalHeaderLabels(vlabel)
            for i in range(row):  #set the initial values
                for j in range(4):
                    self.uiplotscale.scaleTW.setItem(i,j,QTableWidgetItem(str(self.plotscale[i][j])))
                    self.uiplotscale.scaleTW.item(i,j).setTextAlignment(Qt.AlignCenter)
            self.connect(self.uiplotscale.scaleTW, SIGNAL('cellChanged(int,int)'), self.updatePlotScale) #update the ref scale and plot
            self.connect(self.uiplotscale.closePB,SIGNAL('clicked()'), self.closePlotScale) #close the scale setup window
            
    def updatePlotScale(self): #update the scale of each data in the plot
        row=len(self.selectedplotfiles_rows)
        self.plotscale=[[float(str(self.uiplotscale.scaleTW.item(i,j).text())) for j in range(4)] for i in range(row)]
        self.updatePlotPlot()
        
    def closePlotScale(self): #close the plot scale window'
        self.uiplotscale.close()
        
        
    def dataPeakFit(self): #fit the data in the "Plot" work window
        self.mcafitstatus=0
        self.peakFit()
        
    def peakFit(self): #fit the data with lor and gua functions and the polynomial background
        if self.mcafitstatus==0 and len(self.selectedplotfiles_rows)!=1 and self.mcafitallstatus==0:
            self.messageBox('Warning:: Please select only one file!!')
        elif self.mcafitstatus==1 and len(self.nomMcaData)!=1 and self.mcafitallstatus==0:
            self.messageBox('Warning:: Please select only one mca frame')
        else:
            Dialog=QDialog(self)
            self.uipeakfit=uic.loadUi('peakfit.ui', Dialog)
            self.uipeakfit.show()
            try:
                self.uipeakfit.rangeLineEdit.setText(str(self.peakfitranini)+':'+str(self.peakfitranfin))
            except:
                pass
            peakrow=self.uipeakfit.numberOfPeakSpinBox.value()
            bgrow=self.uipeakfit.bgSpinBox.value()+1
            self.uipeakfit.peakTW.setRowCount(peakrow) #set the row for peak parameter; 3 colomn is fixed
            self.uipeakfit.bgTW.setRowCount(bgrow) #set the row for bg parameters
            self.uipeakfit.peakTW.setHorizontalHeaderLabels(QStringList()<<"Location"<<"Intensity"<<"Width")
            peakvlabel=QStringList()
            bgvlabel=QStringList()
            for i in range(peakrow):
                peakvlabel.append("#"+str(i+1))
            self.uipeakfit.peakTW.setVerticalHeaderLabels(peakvlabel)
            for i in range(bgrow):
                bgvlabel.append("C"+str(i))
            self.uipeakfit.bgTW.setVerticalHeaderLabels(bgvlabel)
            self.fitpeakpara=[[1,1,1] for i in range(peakrow)]
            self.fitbgpara=[0 for i in range(bgrow)]
            self.peakparadic={}
            self.peakbgparadic={}            
            for i in range(peakrow): # set values
                for j in range(3):
                    self.uipeakfit.peakTW.setItem(i,j,QTableWidgetItem(str(self.fitpeakpara[i][j])))
                    self.peakparadic[3*i+j]=[self.fitpeakpara[i][j],False,None,None]
            for i in range(bgrow): 
                self.uipeakfit.bgTW.setItem(i,0,QTableWidgetItem(str(self.fitbgpara[i])))
                self.peakbgparadic[i]=[self.fitbgpara[i],False,None,None]
            self.updatePeakParaName()
            self.connect(self.uipeakfit.closePushButton, SIGNAL('clicked()'), self.closePeakFit) #close the peak fit window
            self.connect(self.uipeakfit.fitPushButton, SIGNAL('clicked()'), self.fitPeakFit) #fit the selected data
            self.connect(self.uipeakfit.numberOfPeakSpinBox, SIGNAL('valueChanged(int)'), self.modPeakPara) #change parameters with number of peaks
            self.connect(self.uipeakfit.bgSpinBox, SIGNAL('valueChanged(int)'), self.modBGPara) #change parameters for the bg
            self.connect(self.uipeakfit.exportFitPushButton,SIGNAL('clicked()'), self.savePeakFit) #save peak fit
            self.connect(self.uipeakfit.exportParaPushButton,SIGNAL('clicked()'), self.savePeakPara) #save peak fit
            self.connect(self.uipeakfit.peakTW,SIGNAL('cellDoubleClicked(int,int)'),partial(self.setupPeakPara,'peak')) #setup the peak para limits. 
            self.connect(self.uipeakfit.bgTW,SIGNAL('cellDoubleClicked(int,int)'),partial(self.setupPeakPara,'bg')) #setup the peak bg para limits.
            
    
    def closePeakFit(self): # close the peak fit window
        self.uipeakfit.close()
        
    def updatePeakParaName(self): #update peak and bg para name
        self.peakparaname=[]
        self.peakbgparaname=[]
        for i in range(self.uipeakfit.peakTW.rowCount()):
            self.peakparaname=self.peakparaname+['Loc'+str(i+1),'Int'+str(i+1),'Wid'+str(i+1)]
        for i in range(self.uipeakfit.bgTW.rowCount()):
            self.peakbgparaname.append('C'+str(i))
            
    
    def modPeakPara(self): #modify peak parameter table
        diff=self.uipeakfit.peakTW.rowCount()-self.uipeakfit.numberOfPeakSpinBox.value()
        row=self.uipeakfit.peakTW.rowCount()
        if diff>0:
            for i in range(diff):
                self.uipeakfit.peakTW.removeRow(row-i-1)
                for j in range(3):
                    self.peakparadic.pop(3*(row-i-1)+j)
        else:
            for i in range(-diff):
                self.uipeakfit.peakTW.insertRow(row)
                for j in range(3):
                    self.uipeakfit.peakTW.setItem(row,j,QTableWidgetItem('1/1/1'.split('/')[j]))
                    self.peakparadic[3*row+j]=[1,False,None,None]
                row=row+1
        peakvlabel=QStringList()
        for i in range(row-diff):
                peakvlabel.append("#"+str(i+1))
        self.uipeakfit.peakTW.setVerticalHeaderLabels(peakvlabel)
        self.updatePeakParaName()
        
        
    def modBGPara(self): #modify bg parameter table
        diff=self.uipeakfit.bgTW.rowCount()-self.uipeakfit.bgSpinBox.value()-1
        row=self.uipeakfit.bgTW.rowCount()
        if diff>0:
            for i in range(diff):
                self.uipeakfit.bgTW.removeRow(row-i-1)
                self.peakbgparadic.pop(row-i-1)
        else:
            for i in range(-diff):
                self.uipeakfit.bgTW.insertRow(row)
                self.uipeakfit.bgTW.setItem(row,0,QTableWidgetItem('0'))
                self.peakbgparadic[row]=[0,False,None,None]
                row=row+1
        bgvlabel=QStringList()
        for i in range(row-diff):
                bgvlabel.append("C"+str(i))
        self.uipeakfit.bgTW.setVerticalHeaderLabels(bgvlabel)
        self.updatePeakParaName()
        
    def setupPeakPara(self, index):
        Dialog=QDialog(self)
        self.uipara=uic.loadUi('peakpara.ui', Dialog)
        if index=='peak':
            selrow=self.uipeakfit.peakTW.currentRow()
            selcol=self.uipeakfit.peakTW.currentColumn()
            self.paranum=selrow*3+selcol
            self.uipara.label.setText('Limits Setup of Parameter:'+self.peakparaname[self.paranum])
            if self.peakparadic[self.paranum][2]!=None: 
                self.uipara.minCB.setCheckState(2)
                self.uipara.minLE.setText(str(self.peakparadic[self.paranum][2]))
            if self.peakparadic[self.paranum][3]!=None: 
                self.uipara.maxCB.setCheckState(2)
                self.uipara.maxLE.setText(str(self.peakparadic[self.paranum][3]))
        else:
            self.paranum=self.uipeakfit.bgTW.currentRow()
            self.uipara.label.setText('Limits Setup of Parameter:'+self.peakbgparaname[self.paranum])
            if self.peakbgparadic[self.paranum][2]!=None: 
                self.uipara.minCB.setCheckState(2)
                self.uipara.minLE.setText(str(self.peakbgparadic[self.paranum][2]))
            if self.peakbgparadic[self.paranum][3]!=None: 
                self.uipara.maxCB.setCheckState(2)
                self.uipara.maxLE.setText(str(self.peakbgparadic[self.paranum][3]))
        self.uipara.show()
        self.connect(self.uipara.cancelPB, SIGNAL('clicked()'), self.cancelPeakPara)
        self.connect(self.uipara.okPB, SIGNAL('clicked()'), partial(self.takePeakPara,index))
        
    def cancelPeakPara(self):
        self.uipara.close()
        
    def takePeakPara(self, index):
        if self.uipara.minCB.checkState()!=0 and self.uipara.maxCB.checkState()!=0 and float(self.uipara.minLE.text())>float(self.uipara.maxLE.text()):
            self.messageBox("Error:: Low constrain must be smaller than high constrain!!!")
        else:
            if index=='peak':
                if self.uipara.minCB.checkState()!=0:
                    self.peakparadic[self.paranum][2]=float(self.uipara.minLE.text())
                else:
                    self.peakparadic[self.paranum][2]=None
                if self.uipara.maxCB.checkState()!=0:
                    self.peakparadic[self.paranum][3]=float(self.uipara.maxLE.text())
                else:
                    self.peakparadic[self.paranum][3]=None
            else:
                if self.uipara.minCB.checkState()!=0:
                    self.peakbgparadic[self.paranum][2]=float(self.uipara.minLE.text())
                else:
                    self.peakbgparadic[self.paranum][2]=None
                if self.uipara.maxCB.checkState()!=0:
                    self.peakbgparadic[self.paranum][3]=float(self.uipara.maxLE.text())
                else:
                    self.peakbgparadic[self.paranum][3]=None
            self.uipara.close()
        
    def fitPeakFit(self):
        if self.mcafitallstatus==0:
            if self.mcafitstatus==0:
                data=np.loadtxt(str(self.plotfiles[self.selectedplotfiles_rows[0]]), comments='#')
#                self.ui.PlotMplWidget.canvas.ax.clear()
#                self.ui.PlotMplWidget.canvas.ax.errorbar(data[:,0],data[:,1],data[:,2],fmt='o-')
            else:
                data=self.nomMcaData[self.selectedMcaScanNums[0]]
            self.fitPeak(data)
        else:  #for mca data fit all
           # print self.selectedMcaScanNums
           # print self.nomMcaData
            self.mcaIntData=[]
            self.ui.mcaDataListWidget.clear()
            for i in range(len(self.selectedMcaScanNums)):
                data=self.nomMcaData[self.selectedMcaScanNums[i]]
                self.mcafitstatus==1
                self.mcaLabelOne=str(i+1)
                self.fitPeak(data)
                self.mcaAcceptPeak(num=i)
            self.mcafitallstatus=0
            self.uipeakfit.close()
        
    
    def fitPeak(self,data): #fit the peak 
        data=data[np.argsort(data[:,0])]
        ini=max(float(str(self.uipeakfit.rangeLineEdit.text()).split(':')[0]),data[0][0])
        fin=min(float(str(self.uipeakfit.rangeLineEdit.text()).split(':')[1]),data[-1][0])
        self.peakfitranini=float(str(self.uipeakfit.rangeLineEdit.text()).split(':')[0])
        self.peakfitranfin=float(str(self.uipeakfit.rangeLineEdit.text()).split(':')[1])
        data1=data[np.where(np.logical_and(data[:,0]>=ini,data[:,0]<=fin))]
        peaknumber=self.uipeakfit.numberOfPeakSpinBox.value()
        peaktype=str(self.uipeakfit.peakTypeComboBox.currentText())
        bgorder=self.uipeakfit.bgSpinBox.value()
        x=data1[:,0]
        y=data1[:,1]
        yerr=data1[:,2]
        yerr[yerr==0]=np.min(yerr[np.nonzero(yerr)]) #replace zero with minumun nonzero value
        x0=np.linspace(x[0],x[-1],len(x)*10)
        self.fitpeakpara=[[float(str(self.uipeakfit.peakTW.item(i,j).text())) for j in range(3)] for i in range(peaknumber)]
        if peaknumber==1:
            inten=np.max(y)
            locat=np.sum(x*y)/np.sum(y)
            width=np.sqrt(np.sum(y*(x-locat)**2)/np.sum(y)/2)
            self.fitpeakpara=[[locat,inten,width]]
        self.fitbgpara=[float(str(self.uipeakfit.bgTW.item(i,0).text())) for i in range(bgorder+1)]
        if self.mcafitstatus!=0 or self.mcafitallstatus!=0:   # for fluorescence fitting 
            if peaktype=='Gaussian':        
                p=[[value for sublist in self.fitpeakpara for value in sublist],self.fitbgpara]
                flag=[0,peaknumber]
            else:
                p=[[value for sublist in self.fitpeakpara for value in sublist],self.fitbgpara]
                flag=[1,peaknumber]
            p0=[float(value) for sublist in p for value in sublist]
            pfit,pcov,info,errmsg,ier=leastsq(self.peakres,p0,args=(x,y,yerr,flag), full_output=1)
     #       print pcov
    #        print info,errmsg,ier
            residual=self.peakres(pfit,x,y,yerr,flag)
            self.chisquare=np.sum(residual**2)/(len(x)-len(p0))
            self.fiterror=[np.sqrt(self.chisquare*pcov[i,i]) for i in range(len(p0))]
            if flag[0]==0:
                self.peakfitdatay=sum(self.gaufun(x0,pfit[3*i:3*i+3]) for i in range(flag[1]))+self.polyfun(x0,pfit[3*flag[1]:])
            else:
                self.peakfitdatay=sum(self.lorfun(x0,pfit[3*i:3*i+3]) for i in range(flag[1]))+self.polyfun(x0,pfit[3*flag[1]:])
            self.peakfitdata=[[x0[i],self.peakfitdatay[i]] for i in range(len(x0))]
        else: # regluar fitting 
            self.fitpeakpara=[value for sublist in self.fitpeakpara for value in sublist]
            for i in range(3*peaknumber):
                self.peakparadic[i][0]=self.fitpeakpara[i]
            for i in range(bgorder+1):
                self.peakbgparadic[i][0]=float(str(self.uipeakfit.bgTW.item(i,0).text()))
            index=self.uipeakfit.peakTW.selectionModel().selectedIndexes()  #get selected cell for peak
            bgindex=self.uipeakfit.bgTW.selectionModel().selectedIndexes() #get selected cell for bg 
            selrows=[self.uipeakfit.peakTW.row(self.uipeakfit.peakTW.itemFromIndex(index[i])) for i in range(len(index))] #get selected cells for peak and bg
            selcols=[self.uipeakfit.peakTW.column(self.uipeakfit.peakTW.itemFromIndex(index[i])) for i in range(len(index))]
            selbgrows=[self.uipeakfit.bgTW.row(self.uipeakfit.bgTW.itemFromIndex(bgindex[i])) for i in range(len(bgindex))]
            #print selrows, selcols, selbgrows
            selparas=[3*selrows[i]+selcols[i] for i in range(len(selrows))]
            for i in range(len(self.peakparadic)):  #set selected peak parameters to be varied
                if i in selparas:
                    self.peakparadic[i][1]=True
                else:
                    self.peakparadic[i][1]=False
            for i in range(len(self.peakbgparadic)):  #set selected bg paramenter to be varied
                if i in selbgrows:
                    self.peakbgparadic[i][1]=True
                else:
                    self.peakbgparadic[i][1]=False
                
            self.peakparameter=Parameters()  #set up the parameter for lmfit
            for i in range(len(self.peakparadic)):
                self.peakparameter.add(self.peakparaname[i], value=self.peakparadic[i][0],vary=self.peakparadic[i][1],min=self.peakparadic[i][2],max=self.peakparadic[i][3])
            for i in range(len(self.peakbgparadic)):
                self.peakparameter.add(self.peakbgparaname[i], value=self.peakbgparadic[i][0],vary=self.peakbgparadic[i][1],min=self.peakbgparadic[i][2],max=self.peakbgparadic[i][3])
            
            self.peakresult=minimize(self.peak2min, self.peakparameter, args=(x,y,yerr))
            print(fit_report(self.peakresult))

            peakfitres=[self.peakresult.params[self.peakparaname[i]].value for i in range(peaknumber*3)]
            bgfitres=[self.peakresult.params[self.peakbgparaname[i]].value for i in range(bgorder+1)]
            if peaktype=='Gaussian':
                self.peakfitdatay=sum(self.gaufun(x0,peakfitres[3*i:3*i+3]) for i in range(peaknumber))+self.polyfun(x0,bgfitres)
            else:
                self.peakfitdatay=sum(self.lorfun(x0,peakfitres[3*i:3*i+3]) for i in range(peaknumber))+self.polyfun(x0,bgfitres)
            self.peakbgfitdata=[[x0[i],self.polyfun(x0[i],bgfitres)] for i in range(len(x0))]
            self.peakfitdata=[[x0[i],self.peakfitdatay[i]] for i in range(len(x0))]
        if self.mcafitstatus==0 and self.mcafitallstatus==0:
            self.fitplotstatus=1
            self.updatePlotPlot()
#            self.ui.PlotMplWidget.canvas.ax.plot(x0,self.peakfitdatay,'r-')
#            self.ui.PlotMplWidget.canvas.draw()
        elif self.mcafitstatus!=0 and self.mcafitallstatus==0:
            self.mcafitplotstatus=1
            self.updateMcaPlotData()
        else:
            self.ui.mcaPlotMplWidget.canvas.ax.clear()
            self.ui.mcaPlotMplWidget.canvas.ax.errorbar(data[:,0],data[:,1],data[:,2],fmt='o-', label=self.mcaLabelOne)
            self.ui.mcaPlotMplWidget.canvas.ax.plot(x0,self.peakfitdatay,'r-')
            fit=np.array(self.peakfitdata)
            xran=np.abs(fit[:,0][-1]-fit[:,0][0])
            yran=np.max(fit[:,1])-np.min(fit[:,1])
            self.ui.mcaPlotMplWidget.canvas.ax.set_xlim(fit[:,0][0]-0.25*xran,fit[:,0][-1]+0.25*xran)
            self.ui.mcaPlotMplWidget.canvas.ax.set_ylim(np.min(fit[:,1])-0.2*yran,np.max(fit[:,1])+0.2*yran)
            self.ui.mcaPlotMplWidget.canvas.ax.set_xlabel('Energy')
            self.ui.mcaPlotMplWidget.canvas.ax.set_ylabel('Intensity')
            self.ui.mcaPlotMplWidget.canvas.ax.set_title('Fluorescence Spectrum')
            self.ui.mcaPlotMplWidget.canvas.ax.legend(loc=2,frameon=False,scatterpoints=0,numpoints=1)
            self.ui.mcaPlotMplWidget.canvas.draw()
            self.ui.mcaPlotMplWidget.canvas.flush_events()
            time.sleep(0.25)
        #self.fitpeakpara=[[pfit[3*i],pfit[3*i+1],abs(pfit[3*i+2])] for i in range(peaknumber)]
        #self.fitbgpara=[pfit[3*peaknumber+i] for i in range(bgorder+1)]
        if self.mcafitstatus!=0 or self.mcafitallstatus!=0:   # for fluorescence fitting 
            self.fitpeakpara=[[pfit[3*i],pfit[3*i+1],abs(pfit[3*i+2])] for i in range(peaknumber)]
            self.fitbgpara=[pfit[3*peaknumber+i] for i in range(bgorder+1)]
            for i in range(peaknumber): # set values
                for j in range(3):
                    self.uipeakfit.peakTW.setItem(i,j,QTableWidgetItem(str(float('%.4g' % self.fitpeakpara[i][j]))))
            for i in range(bgorder+1): 
                    self.uipeakfit.bgTW.setItem(i,0,QTableWidgetItem(str(float('%.4g' % self.fitbgpara[i]))))
            self.uipeakfit.chiLineEdit.setText(str(float('%.3f' % self.chisquare)))
        else:    # for regular fitting
            for i in range(peaknumber): # set values
                for j in range(3):
                    self.uipeakfit.peakTW.setItem(i,j,QTableWidgetItem(str(float('%.4g' % peakfitres[3*i+j]))))
            for i in range(bgorder+1): 
                    self.uipeakfit.bgTW.setItem(i,0,QTableWidgetItem(str(float('%.4g' % bgfitres[i]))))
            self.uipeakfit.chiLineEdit.setText(str(float('%.3f' % self.peakresult.redchi)))
    
    def gaufun(self,x,p):
        return abs(p[1])*np.exp(-(x-p[0])**2/2/p[2]**2)
        
    def lorfun(self,x,p):
        return abs(p[1])/(1+(x-p[0])**2/p[2]**2)

    def polyfun(self,x,p):
        if self.uipeakfit.invCheckBox.checkState()!=0:
            return sum([p[i]/x**i for i in range(len(p))])
        else:
            return sum([p[i]*x**i for i in range(len(p))])
        
    def peak2min(self, params, x, y, yerr):
        peakrow=self.uipeakfit.numberOfPeakSpinBox.value()
        peaktype=str(self.uipeakfit.peakTypeComboBox.currentText())
        bgrow=self.uipeakfit.bgSpinBox.value()
        peakpara=[params[self.peakparaname[i]].value for i in range(peakrow*3)]
        bgpara=[params[self.peakbgparaname[i]].value for i in range(bgrow+1)]
       # print peakpara, bgpara
        sum=self.polyfun(x,bgpara)        
        for i in range(peakrow):
            if peaktype=='Gaussian':
                sum=sum+self.gaufun(x,peakpara[3*i:3*i+3])
            else:
                sum=sum+self.lorfun(x,peakpara[3*i:3*i+3])
        return (y-sum)/yerr
        
    def peakres(self,p,x,y,yerr,flag): #define gauss and lorent  
        sum=self.polyfun(x,p[3*flag[1]:])        
        if flag[0]==0:
            for i in range(flag[1]):
                sum=sum+self.gaufun(x,p[3*i:3*i+3])
        else:
            for i in range(flag[1]):
                sum=sum+self.lorfun(x,p[3*i:3*i+3])
        return (y-sum)/yerr
    
    def savePeakFit(self): # save fit data
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Fit Data',directory=self.directory))
        fitfile=self.saveFileName+'_fit.txt'
       # print self.peakfitdata
        np.savetxt(fitfile,self.peakfitdata,fmt='%.4e\t%.4e')

    def savePeakPara(self): # save fit parameters
        self.saveFileName=str(QFileDialog.getSaveFileName(caption='Save Fit Paremeters',directory=self.directory))
        fid=open(self.saveFileName+'_par.txt','w')
        try:  #for pc
            fid.write('Data File:\t'+str(self.plotfiles[self.selectedplotfiles_rows[0]].split('\\')[-1])+'\n')
        except:  #for mac and linux
            fid.write('Data File:\t'+str(self.plotfiles[self.selectedplotfiles_rows[0]].split('/')[-1])+'\n')
        fid.write('Peak Type:\t'+str(self.uipeakfit.peakTypeComboBox.currentText())+'\n')
        fid.write('# of Peaks:\t'+str(int(len(self.peakparaname)/3))+'\n')
        fid.write('BG Order:\t'+str(len(self.peakbgparaname)-1)+'\n')       
        fid.write('Chi_Square:\t'+str(float('%.4f' % self.peakresult.redchi))+'\n')
        fid.write('*********************\n')
        fid.write('\t\tValue\t\tVary\tStderr\t\tMin\tMax\n')
        for i in range(int(len(self.peakparaname)/3)):
            fid.write('Peak #'+str(i+1)+' pos:\t'+format(self.peakresult.params[self.peakparaname[i*3]].value, '.3e')+'\t'+str(self.peakparadic[i*3][1])+'\t'+format(self.peakresult.params[self.peakparaname[i*3]].stderr, '.3e')+'\t'+str(self.peakparadic[i*3][2])+'\t'+str(self.peakparadic[i*3][3])+'\n')
            fid.write('Peak #'+str(i+1)+' int:\t'+format(self.peakresult.params[self.peakparaname[i*3+1]].value, '.3e')+'\t'+str(self.peakparadic[i*3+1][1])+'\t'+format(self.peakresult.params[self.peakparaname[i*3+1]].stderr, '.3e')+'\t'+str(self.peakparadic[i*3+1][2])+'\t'+str(self.peakparadic[i*3+1][3])+'\n')
            fid.write('Peak #'+str(i+1)+' wid:\t'+format(self.peakresult.params[self.peakparaname[i*3+2]].value, '.3e')+'\t'+str(self.peakparadic[i*3+2][1])+'\t'+format(self.peakresult.params[self.peakparaname[i*3+2]].stderr, '.3e')+'\t'+str(self.peakparadic[i*3+2][2])+'\t'+str(self.peakparadic[i*3+2][3])+'\n')
        for i in range(len(self.fitbgpara)):
            fid.write('Background C'+str(i)+':\t'+format(self.peakresult.params[self.peakbgparaname[i]].value, '.3e')+'\t'+str(self.peakbgparadic[i][1])+'\t'+format(self.peakresult.params[self.peakbgparaname[i]].stderr, '.3e')+'\t'+str(self.peakbgparadic[i][2])+'\t'+str(self.peakbgparadic[i][3])+'\n')
        fid.close()
        
########################  for 2d plot######################

    def update2dPlotFile(self): #update 2dplot files in the listwidget
        self.ui.twodPlotFileListWidget.clear()
        for i in range(len(self.twodplotfiles)):
            try:  #for pc
                self.ui.twodPlotFileListWidget.addItem('#'+str(i+1)+self.halftab+str(self.twodplotfiles[i].split('\\')[-2])+'\\'+str(self.twodplotfiles[i].split('\\')[-1]))
            except: #for mac and linux
                self.ui.twodPlotFileListWidget.addItem('#'+str(i+1)+self.halftab+str(self.twodplotfiles[i].split('/')[-2])+'/'+str(self.twodplotfiles[i].split('/')[-1]))

    def add2dPlotFile(self): #add plot files into the listwidget and deselect all ref files in the listwidget
        f=QFileDialog.getOpenFileNames(caption='Select Multiple Files to import', directory=self.directory, filter='Files (*.*)')
        self.twodplotfiles=self.twodplotfiles+map(str, f)
        self.update2dPlotFile()
        
    def updateSelected2dPlotFile(self): #update the selected data files in the listwidget
        selected2dplotfiles=self.ui.twodPlotFileListWidget.selectedItems()
        self.selected2dplotfiles_rows=[]
        for item in selected2dplotfiles:
            self.selected2dplotfiles_rows.append(self.ui.twodPlotFileListWidget.row(item))
        self.selected2dplotfiles_rows.sort()
        if len(self.selected2dplotfiles_rows)!=0:
            self.update2dPlotData()
        
    def remove2dPlotFile(self): #remove plot files in the listwidget and deselect all ref files in the listwidget
        items=self.ui.twodPlotFileListWidget.selectedItems()
        items.reverse()
        for item in items:
            self.twodplotfiles.pop(self.ui.twodPlotFileListWidget.row(item))
        self.ui.twodPlotFileListWidget.clear()
        self.update2dPlotFile()        
     
    def update2dPlotData(self):
        data=np.loadtxt(str(self.twodplotfiles[self.selected2dplotfiles_rows[0]]), comments='#')
        x=np.unique(data[:,0])
        y=np.unique(data[:,1])
        if len(x)*len(y)!=len(data) or len(x)==1 or len(y)==1:
            self.messageBox('Warning:: The selected data file is not in 2d data format!')
        else:
            self.plot2dindex=0
            self.data2dint=np.reshape(data[:,2],(len(y),len(x)))
            self.extent=[x[0],x[-1],y[0],y[-1]]
            self.Zdata=self.data2dint
            self.data2dx, self.data2dy=np.meshgrid(x,y)
            self.imageMax=np.max(self.data2dint)
            self.imageMin=np.min(self.data2dint[np.where(self.data2dint>0)])
            self.ui.twodMaxLineEdit.setText(str(self.imageMax)) #set the vmax/vmin in the lineedit and slider 
            self.ui.twodMinLineEdit.setText(str(self.imageMin))
            self.ui.twodMaxHorizontalSlider.setValue(int(float(self.ui.twodMaxLineEdit.text())*100/(self.imageMax-self.imageMin)))
            self.ui.twodMinHorizontalSlider.setValue(int(float(self.ui.twodMinLineEdit.text())*100/(self.imageMax-self.imageMin)))
            self.update2dPlotPlot()
            
    def update2dMaxSlider(self):
        self.ui.twodMaxLineEdit.setText(str(self.ui.twodMaxHorizontalSlider.value()*(self.imageMax-self.imageMin)/100.0))
        self.ui.twodMinHorizontalSlider.setMaximum(int(float(self.ui.twodMaxLineEdit.text())*100.0/(self.imageMax-self.imageMin)))
       # print 'I ma here'
        self.update2dPlotPlot()
        
    def update2dMinSlider(self):
        self.ui.twodMinLineEdit.setText(str(self.ui.twodMinHorizontalSlider.value()*(self.imageMax-self.imageMin)/100.0))
        self.update2dPlotPlot()
        
    def update2dIntPolData(self):
        if len(self.selected2dplotfiles_rows)==0:
            self.messageBox('Please select one 2d data file first!')
        else:
            self.plot2dindex=1            
            pointden=self.ui.twodDenFacSpinBox.value()
            fint=interp1d(self.data2dx[0], self.data2dint, kind='linear')
            newdata2dx=np.linspace(self.data2dx[0][0],self.data2dx[0][-1],len(self.data2dx[0])*pointden)
            self.intpoldata2dint=fint(newdata2dx)
            self.Zdata = self.intpoldata2dint
           # print np.min(self.intpoldata2dint), np.min(self.data2dint)
            self.update2dPlotPlot()
    
    def update2dPlotPlot(self):
        self.ui.twodPlotMplWidget.canvas.fig.clf()
        cmap=str(self.ui.twodCMapComboBox.currentText())
        vmax=float(self.ui.twodMaxLineEdit.text())
        vmin=float(self.ui.twodMinLineEdit.text())
        labelsize=self.ui.twodLabelSizeSpinBox.value()
        ticksize=self.ui.twodTickSizeSpinBox.value()
        self.plot2dax=self.ui.twodPlotMplWidget.canvas.fig.add_subplot(1,1,1)
        self.plot2dax.set_xlabel(r'$Q_{xy}$'+' '+r'$[\AA^{-1}]$', fontsize=labelsize)
        self.plot2dax.set_ylabel(r'$Q_{z}$'+' '+r'$[\AA^{-1}]$', fontsize=labelsize)
        if self.plot2dindex==0:
            if self.ui.twodLogIntCheckBox.checkState()!=0:
                self.plot2d_p=self.plot2dax.pcolor(self.data2dx, self.data2dy, np.log10(self.data2dint),cmap=pl.cm.get_cmap(cmap),vmin=np.log10(vmin), vmax=np.log10(vmax))
            else:
                self.plot2d_p=self.plot2dax.pcolor(self.data2dx, self.data2dy, self.data2dint, cmap=pl.cm.get_cmap(cmap), vmin=vmin, vmax=vmax)
            self.plot2dax.set_xlim((np.min(self.data2dx),np.max(self.data2dx)))
            self.plot2dax.set_ylim((np.min(self.data2dy),np.max(self.data2dy)))
        else:
            interpol=str(self.ui.twodIntPolComboBox.currentText())
            if self.ui.twodLogIntCheckBox.checkState()!=0:
                self.plot2d_p=self.plot2dax.imshow(np.log10(self.intpoldata2dint),interpolation=interpol,extent=self.extent,cmap=cmap,vmin=np.log10(vmin), vmax=np.log10(vmax), aspect='auto', origin='lower')
            else:
                self.plot2d_p=self.plot2dax.imshow(self.intpoldata2dint,interpolation=interpol,extent=self.extent,cmap=cmap,vmin=vmin, vmax=vmax, aspect='auto', origin='lower')
        if self.ui.twodColBarCheckBox.checkState()!=0:
            cb=self.ui.twodPlotMplWidget.canvas.fig.colorbar(self.plot2d_p)
            cb.ax.tick_params(labelsize=ticksize)
        self.plot2dax.tick_params(labelsize=ticksize)
        self.plot2dax.format_coord=self.format_pil_coord
        self.ui.twodPlotMplWidget.canvas.draw()
    
    
    def questionDialog(self, question='Do you want to continue?'):
        reply=QMessageBox.question(self, 'Message',
            question, QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        return reply
        
        