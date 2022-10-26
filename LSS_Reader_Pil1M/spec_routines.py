#!/usr/bin/env python
import sys
import os
import numpy as np
from PyQt4.QtGui import QProgressDialog
from PyQt4.QtCore import Qt

class specread:
    def __init__(self, specfile, startLineNum=0, endScanNum=0, beamline='APS-15IDC',det='CCD',data={},par={}):
        self.Data=data
        self.Par=par
        self.specfile=specfile
        if beamline=='APS-15IDC':
            self.APS_15IDC(startLineNum=startLineNum, endScanNum=endScanNum, det=det)
        if beamline=='APS-9IDC':
            self.APS_9IDC(startLineNum=startLineNum,endScanNum=endScanNum, det=det)
        
    def updateProgress(self):
        self.progressDialog.setValue(self.progressDialog.value()+1)

    
    
    def APS_15IDC(self,startLineNum=0, endScanNum=0, det='CCD'):
        """
        Function to read a complete spec File collected at APS 15IDC 
        """
        self.progressDialog=QProgressDialog('Reading scans form SPEC File:','Abort',0,100)
        self.progressDialog.setWindowModality(Qt.WindowModal)
        self.progressDialog.setWindowTitle('Wait')
        self.progressDialog.setAutoClose(True)
        self.progressDialog.setAutoReset(True)
        self.progressDialog.setMinimum(1)
        self.Data['YCol']='Apex2'
        self.Data['NCol']='Monc'
        fid=open(self.specfile)
        fdata=fid.readlines()
        self.SpecFileFull=fdata
        fid.close()
        if fdata[0][:2]!='#F':
            self.Data['NumOfScans']=0
            self.Data['Message']='The file is not a valid specfile!!'
            print 'Error:: The file is not a valid specfile!!'
        else:
            startScanLineNums=[i for i in range(startLineNum,len(fdata)) if fdata[i][:2]=='#S']
            #print startScanLineNums
            #temp=[startScanLineNums[i]-startScanLineNums[i-1] for i in range(1,len(startScanLineNums))]
            #print temp
            
            #print np.where(np.diff(startScanLineNums)==2)
            self.progressDialog.setMaximum(len(startScanLineNums))
            self.progressDialog.show()
            if startLineNum>0:
                startScanLineNums=sorted(startScanLineNums)
            numOfScans=len(startScanLineNums)
            scanLines=[fdata[i] for i in startScanLineNums]
            templist=[scanLines[i].split()[1] for i in range(len(scanLines))]
           # print templist
            templist2=[int(templist[i])-int(templist[i-1]) for i in range(1,len(templist))]
           # print templist2
            
            if startLineNum==0:
                tmp=0
                self.Data['NumOfScans']=0#numOfScans
                self.Data['ScanNums']=[]
                offset=0
#                self.Data['ScanLines']=[]#scanLines
#                self.Data['StartScanLineNums']=[]#startScanLineNums
            else:
                tmp=self.Data['NumOfScans']
                self.Data['NumOfScans']=self.Data['NumOfScans']-1#+numOfScans
                offset=1
#                self.Data['ScanLines']=self.Data['ScanLines'][:-1]#+scanLines[1:]
#                self.Data['StartScanLineNums']=self.Data['StartScanLineNums'][:-1]#+startScanLineNums
            if int(scanLines[-1].split()[1])!=len(startScanLineNums)+endScanNum-offset:
                print len(startScanLineNums), scanLines[-1].split()[1], endScanNum
                print np.where(np.array(templist2)!=1)[0]+2
                self.Data['Error']=True
                self.Data['Message']='There are identical/missing scans in the file'+'\n'+'Identical/missing scan numbers:'+str(np.where(np.array(templist2)!=1)[0]+2)
            else:
                self.Data['Error']=False
            for i in range(numOfScans): 
                start=startScanLineNums[i]+1
                line=fdata[start]
                num=int(fdata[startScanLineNums[i]].split()[1])
                i=i+tmp
                self.Data[num]={}
                self.Par[num]={}
                if fdata[start-1].split()[2]=='getandsave_mca':
                    self.Par[num]['Detector']='Vortex'
                    self.Data[num]['ScanVar']='Empty'
                else:
                    self.Par[num]['Detector']='Monitor'
                tmpdata=[]
                while line[:2]!='\n' and line[:2]!='#C':
                    if line[:2]=='#P':
                        parName=line[4:].split()
                        start=start+1
                        parValue=map(eval,fdata[start][1:].split())
                        for j in range(len(parName)):
                            self.Par[num][parName[j]]=parValue[j]
                    if line[:2]=='#W':
                        tmppar=line[2:].split()
                        self.Par[num]['Wavelength']=eval(tmppar[1])
                    if line[:3]=='#G0':
                        self.Par[num]['g_l1']=float(line[4:].split()[5])
                        self.Par[num]['g_l2']=float(line[4:].split()[6])
                        self.Par[num]['g_l3']=float(line[4:].split()[7])
                    if line[:2]=='#A':
                        tmppar=line[2:].split()
                        self.Par[num]['Absorber']=eval(tmppar[1])
                    if line[:2]=='#Q':
                        tmppar=line[2:].split()
                        self.Par[num]['Q']=map(eval, tmppar)
                    if line[:2]=='#V':
                        self.Par[num]['Detector']='Vortex'
                    if line[:3]=='#B0':
                        tmppar=line[3:].split('.')
                        self.Par[num]['ImageNumber']=len(line[3:].split('_')[-1].split('.')[0])
#                        if tmppar[1]=='tif\n':
#                            self.Par[num]['Detector']='Pilatus'
#                        elif tmppar[1]=='sfrm\n':
#                            self.Par[num]['Detector']='Bruker'
                    if line[:3]=='#B1':
                        try:
                            tmppar=map(eval, line[3:].split())
                        except:
                            tmppar=map(eval, line[3:].split()[:-1])
                        self.Par[num]['Detector']='Pilatus'
                        self.Par[num]['DBPos']=tmppar[:2]
                        self.Par[num]['S2D_Dist']=tmppar[2]
                        self.Par[num]['S7D_Dist']=tmppar[3]
                    if line[:3]=='#B5':
                        tmppar=map(eval, line[3:].split())
                        self.Par[num]['Detector']='Pilatus1M'
                        self.Par[num]['DBPos']=tmppar[:2]
                        self.Par[num]['S2D_Dist']=tmppar[2]
                        self.Par[num]['S7D_Dist']=tmppar[3]
                    if line[:2]=='#L':
                        scanVar=line[3:-1].split()
                        self.Data[num]['ScanVar']=scanVar
                    if line[0]!='#':
                        try:
                            tmpdata.append(map(eval, line.split(  )))
                        except:
                            self.Data[num]['Message']='Something wrong with Scan Number %d',num,'.Please check the the scan in the specfile.'
                            print 'Something wrong with Scan Number %d',num
                    start=start+1
                    try:
                        line=fdata[start]
                    except:
                        break
                if self.Data[num]['ScanVar']!='Empty':
                    for j in range(len(scanVar)):
                        try:
                            self.Data[num][scanVar[j]]=np.array(tmpdata,dtype='float')[:,j]
                        except:
                            self.Data[num][scanVar[j]]=None
                if len(self.Par[num])==1:
                    self.Par[num]['Message']='No parameters!!'
                self.progressDialog.setLabelText('Reading Scan #'+str(num))         
                self.updateProgress()
                self.Data['NumOfScans']=num
#                self.Data['ScanLines']=self.Data['ScanLines']+[scanLines[num-tmp]]
                self.Data[num]['ScanLine']=fdata[startScanLineNums[i-tmp]]
                self.Data[num]['StartScanLineNum']=startScanLineNums[i-tmp]
                self.endLineNum=startScanLineNums[i-tmp]
                self.Data['ScanNums'].append(num)
                if self.progressDialog.wasCanceled()==True:
                    break
        self.progressDialog.hide()
            
                  
    
    def APS_9IDC(self,startLineNum=0,endScanNum=0, det='CCD'):
        """
        Function to read a complete spec File collected at APS 9IDC 
        """
        self.progressDialog=QProgressDialog('Reading scans form SPEC File:','Abort',0,100)
        self.progressDialog.setWindowModality(Qt.WindowModal)
        self.progressDialog.setWindowTitle('Wait')
        self.progressDialog.setAutoClose(True)
        self.progressDialog.setAutoReset(True)
        self.progressDialog.setMinimum(1)
        self.Data['YCol']='Bicron1'
        self.Data['NCol']='i2'
        fid=open(self.specfile)
        fdata=fid.readlines()
        self.SpecFileFull=fdata
        fid.close()
        if fdata[0][:2]!='#F':
            self.Data['NumOfScans']=0
            self.Data['Message']='The file is not a valid specfile!!'
            print 'Error:: The file is not a valid specfile!!'
        else:
            startScanLineNums=[i for i in range(startLineNum, len(fdata)) if fdata[i][:2]=='#S']
            self.progressDialog.setMaximum(len(startScanLineNums))
            self.progressDialog.show()
            self.endLineNum=startScanLineNums[-1]
            if startLineNum>0:
                startScanLineNums=sorted(startScanLineNums)
            self.Data['StartScanLineNums']=startScanLineNums
            numOfScans=len(self.Data['StartScanLineNums'])
            scanLines=[fdata[i] for i in startScanLineNums]
            if startLineNum==0:
                tmp=0
                self.Data['NumOfScans']=0#numOfScans
                self.Data['ScanLines']=[]#scanLines
                self.Data['StartScanLineNums']=[]#startScanLineNums
                self.Par['ParName']=[]
                for i in range(startScanLineNums[0]):
                    line=fdata[i].split()
                    if fdata[i][:2]=='#O':
                        self.Par['ParName']=self.Par['ParName']+line[1:]
            else:
                tmp=self.Data['NumOfScans']
                self.Data['NumOfScans']=self.Data['NumOfScans']-1#+numOfScans
                self.Data['ScanLines']=self.Data['ScanLines'][:-1]#+scanLines[1:]
                self.Data['StartScanLineNums']=self.Data['StartScanLineNums'][:-1]#+startScanLineNums
            for i in range(numOfScans):
                start=startScanLineNums[i]+1
                line=fdata[start]
                i=i+tmp
                self.Data[i]={}
                self.Par[i]={}
                if fdata[start-1].split()[2]=='getandsave_mca' or fdata[start-1].split()[2]=='MCAscanpt':
                    self.Par[i]['Mca']=1
                    self.Data[i]['ScanVar']='Empty'
                else:
                    self.Par[i]['Mca']=0
                self.Par[i]['CCD']=0
                tmpdata=[]
                pstart=0
                while line[:2]!='\n' and line[:2]!='#C':
                    if line[:2]=='#P':
                        parValue=map(eval,fdata[start].split()[1:])
                        for j in range(len(parValue)):
                            self.Par[i][self.Par['ParName'][pstart]]=parValue[j]
                            pstart=pstart+1
                    if line[:2]=='#Q':
                        tmppar=line[2:].split()
                        self.Par[i]['Q']=map(eval, tmppar)
                    if line[:2]=='#L':
                        scanVar=line[3:-1].split()
                        self.Data[i]['ScanVar']=scanVar
                    if line[0]!='#':
                        tmpdata.append(map(eval, line.split(  )))
                    start=start+1
                    line=fdata[start]
                for j in range(len(scanVar)):
                    try:
                        self.Data[i][scanVar[j]]=np.array(tmpdata)[:,j]
                    except:
                        self.Data[i][scanVar[j]]=None
                if len(self.Par[i])==1:
                    self.Par[i]['Message']='No parameters!!'
                self.progressDialog.setLabelText('Reading scans form SPEC File: '+str(i+1)) 
                self.updateProgress()
                self.Data['NumOfScans']=i
                self.Data['ScanLines']=self.Data['ScanLines']+[scanLines[i-tmp]]
                self.Data['StartScanLineNums']=self.Data['StartScanLineNums']+[startScanLineNums[i-tmp]]
                self.endLineNum=startScanLineNums[i-tmp]
                if self.progressDialog.wasCanceled()==True:
                    break
        self.progressDialog.hide()
                    
