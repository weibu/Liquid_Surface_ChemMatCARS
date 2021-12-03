#!/usr/bin/env python
import sys
import os
import numpy as np

class mcaread:
    def __init__(self, mcafile, beamline='APS-15IDC'):
        self.Data={}
        self.Par={}
        self.mcafile=mcafile
        if beamline=='APS-15IDC':
            self.APS_15IDC()
        if beamline=='APS-9IDC':
            self.APS_9IDC()
    
    
    def APS_15IDC(self):
        """
        Function to read a complete spec mca collected at APS 15IDC 
        """
        fid=open(self.mcafile)
        fdata=fid.readlines()
        fid.close()
        if fdata[0][:2]!='#F':
            self.Data['NumOfScans']=0
            self.Data['Message']='The file is not a valid specfile!!'
            print 'Error:: The file is not a valid specfile!!'
        else:
            startScanLineNums=[i for i in range(len(fdata)) if fdata[i][:2]=='#S']
            self.Data['StartScanLineNums']=startScanLineNums
            numOfScans=len(self.Data['StartScanLineNums'])
            self.Data['NumOfScans']=numOfScans
            scanLines=[fdata[i] for i in startScanLineNums]
            self.Data['ScanLines']=scanLines
            for i in range(numOfScans):
                start=startScanLineNums[i]+1
                line=fdata[start]
                self.Data[i]={}
                self.Par[i]={}
                tmpdata=[]
                while line[:2]!='@A':
                    if line[:7]=='#@CTIME':
                        tmppar=line[7:].split()
                        try:
                            self.Par[i]['Time']=map(eval, tmppar)
                        except:
                            self.Par[i]['Time']=[eval(tempar[0]), eval(tempar[0]), eval(tempar[0])]
                    if line[:5]=='#Monc':
                        tmppar=line[5:].split()
                        self.Data[i]['Monc']=eval(tmppar[0])
                        try: # When running for Hutch B and D
                            self.Data[i]['Monb']=eval(tmppar[1])
                            self.Data[i]['Mond']=eval(tmppar[2])
                        except:
                            pass
                    if line[:2]=='#Q':
                        tmppar=line[2:].split()
                        self.Par[i]['Q']=map(eval, tmppar)
                    if line[:7]=='#@CALIB':
                        tmppar=line[7:].split()
                        self.Par[i]['Calib']=map(eval, tmppar)
                    if line[:7]=='#Energy':
                        tempar=line[7:].split()
                        self.Par[i]['Energy']=eval(tempar[0])
                    start=start+1
                    line=fdata[start]
                self.Data[i]['Vortex']=map(eval, line[2:-2].split())
                start=start+1
                line=fdata[start]
                while line!='\n':
                    if line[-2]=='\\':
                        self.Data[i]['Vortex']=self.Data[i]['Vortex']+map(eval, line[:-2].split())
                    else:
                        self.Data[i]['Vortex']=self.Data[i]['Vortex']+map(eval, line[:-1].split())
                    start=start+1
                    line=fdata[start]
                if len(self.Par[i])==0:
                    self.Par[i]['Message']='No parameters!!'  
                  
    
    def APS_9IDC(self):
        """
        Function to read a complete spec mca collected at APS 15IDC 
        """
        fid=open(self.mcafile)
        fdata=fid.readlines()
        fid.close()
        if fdata[0][:2]!='#F':
            self.Data['NumOfScans']=0
            self.Data['Message']='The file is not a valid specfile!!'
            print 'Error:: The file is not a valid specfile!!'
        else:
            startScanLineNums=[i for i in range(len(fdata)) if fdata[i][:2]=='#S']
            self.Data['StartScanLineNums']=startScanLineNums
            numOfScans=len(self.Data['StartScanLineNums'])
            self.Data['NumOfScans']=numOfScans
            scanLines=[fdata[i] for i in startScanLineNums]
            self.Data['ScanLines']=scanLines
            for i in range(numOfScans):
                start=startScanLineNums[i]+1
                line=fdata[start]
                self.Data[i]={}
                self.Par[i]={}
                tmpdata=[]
                while line[:2]!='@A':
                    if line[:2]=='#T':
                        tmppar=line[2:].split()
                        self.Par[i]['Time']=eval(tmppar[0])
                    if line[:3]=='#i2':
                        tmppar=line[3:].split()
                        self.Data[i]['Monc']=eval(tmppar[0])
                    if line[:2]=='#Q':
                        tmppar=line[2:].split()
                        self.Par[i]['Q']=map(eval, tmppar)
                    if line[:7]=='#@CALIB':
                        tmppar=line[7:].split()
                        self.Par[i]['Calib']=map(eval, tmppar)
                    start=start+1
                    line=fdata[start]
                self.Data[i]['Vortex']=map(eval, line[2:-2].split())
                start=start+1
                line=fdata[start]
                while line!='\n':
                    if line[-2]=='\\':
                        self.Data[i]['Vortex']=self.Data[i]['Vortex']+map(eval, line[:-2].split())
                    else:
                        self.Data[i]['Vortex']=self.Data[i]['Vortex']+map(eval, line[:-1].split())
                    start=start+1
                    line=fdata[start]
                if len(self.Par[i])==0:
                    self.Par[i]['Message']='No parameters!!' 