import sqlite3
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import time
import os
from scipy.optimize import curve_fit
import shutil
import h5py
from scipy.interpolate import UnivariateSpline
import filecmp
from scipy import signal
import subprocess

def resFunc(x,a,b,c):
    return (a+b*x)**c

def gausswLine(x,a,b,c,mu,sigma): #guassian + line for fitting
    return (        
        #line
        a+c*(x-mu)+
        #photopeakd
        (b/(sigma*np.sqrt(2.0*np.pi)))*np.exp(-(x-mu)**2/(2.0*sigma**2))
        )

def linFunc(x,b):
    return b*x

def doEcal(k40PeakChan,energies,nonLinearity):

    slope=k40PeakChan/1460.8
    newChans=nonLinearity*linFunc(energies,slope) #units of channels
    #need to take the inverse so do a spline
    spl = UnivariateSpline(energies,newChans,k=2,s=0.0,ext=0)
    ens=np.linspace(0.0,3500.0,5E4)
    chanSpl=spl(ens)

    eCal=[]
    for i in range(len(newChans)):
        index=np.abs(chanSpl-i).argmin()
        eCal.append(ens[index])

    return np.array(eCal)

class database(object): #database class

    def __init__(self,DBfilename):
        self.conn = sqlite3.connect(DBfilename)
        self.conn.row_factory = sqlite3.Row
        #self.getTables()
        #self.tblname=tblname
        #self.getHeader(tblename)
        #self.getData()

    def getTables(self):
        c = self.conn.cursor()
        c.execute("select name from sqlite_master where type = 'table'")
        self.tables=c.fetchall()
        c.close()

    def getHeader(self,tblname):
        c=self.conn.cursor()
        c.execute("select rowid,* from %s" % (tblname))
        self.header=np.array(c.fetchone().keys())
        c.close()

    def getAllData(self,tblname):
        c=self.conn.cursor()
        #c.execute("select rowid,* from %s where tsm>=? and tsm<=?" %(self.tblname,),(tsmbegin,tsmend))
        c.execute("select rowid,* from %s"  % (tblname))
        self.data=self.tableToArray(c)
        c.close()

    def getSomeData(self,tblname,label,t0,tf):
        c=self.conn.cursor()
        c.execute("select rowid,* from %s where %s>=%.1f and %s<%.1f" % (tblname,label,t0,label,tf))
        self.data=self.tableToArray(c)
        return

    def getSomeDataIf(self,tblname,label,t0,tf,x):
        c=self.conn.cursor()
        c.execute("select rowid,* from %s where %s>=%.1f and %s<%.1f and Gate1=%d " % (tblname,label,t0,label,tf,x))
        self.data=self.tableToArray(c)
        return

    def getColumn(self,tbleName,column):
        c=self.conn.cursor()
        c.execute("select %s from %s" % (column,tbleName))
        self.data=c.fetchall()
        c.close()

    def getFirstAndLastRow(self,tblname):
        c=self.conn.cursor()
        c.execute("select rowid,* from %s order by rowid asc limit 1" % (tblname))
        firstRow=self.tableToArray(c)
        c.execute("select rowid,* from %s order by rowid desc limit 1" % (tblname))
        lastRow=self.tableToArray(c)
        return firstRow,lastRow

    def tableToArray(self,c):
        rows=c.fetchall()
        data=[]
        for i in range(len(rows)):
            data.append(list(rows[i]))
        #print data
        #data=np.array(data)
        return data        

    def closeConn(self):
        self.conn.close()
        self.data=[]

def rebin_Archer(x1, y1, x2):
    """
    Rebin histogram values y1 from old bin edges x1 to new edges x2.

    Input
    -----
     * x1 : m+1 array of old bin edges.
     * y1 : m array of old histogram values. This is the total number in 
              each bin, not an average.
     * x2 : n+1 array of new bin edges.

    Returns
    -------
     * y2 : n array of rebinned histogram values.

    The rebinning algorithm assumes that the counts in each old bin are
    uniformly distributed in that bin.

    Bins in x2 that are entirely outside the range of x1 are assigned 0.
    """

    x1 = np.asarray(x1)
    y1 = np.asarray(y1)
    x2 = np.asarray(x2)
    
    
    # allocating y2 vector
    m  = y1.size
    n  = x2.size - 1
    y2 = np.zeros(n, dtype=y1.dtype)

    i_place = np.searchsorted(x2, x1)

    # find out where x2 intersects with x1, this will determine which x2 bins
    # we need to consider
    start_pos = 0
    end_pos = n

    start_pos_test = np.where(i_place == 0)[0]
    if start_pos_test.size > 0:
        start_pos = start_pos_test[-1]

    end_pos_test = np.where(i_place == x1.size)[0]
    if end_pos_test.size > 0:
        end_pos = end_pos_test[0]

    # the first bin totally covers x1 range
    if (start_pos == end_pos - 1
        and i_place[start_pos] == 0
        and i_place[start_pos + 1] == x1.size):
 
        y2[start_pos] = y1[start_pos]
        return y2

    #print len(x1),len(y1),len(y2)
    #print len(i_place)
    #print i_place[980:len(i_place)]
    
    # first(0th) X1 bin
    ibin=0
    i_place_lower=i_place[ibin]
    i_place_upper=i_place[ibin+1]
    if i_place_lower==0 and i_place_upper>0:
        if i_place_upper==1:
            y2_index = i_place_upper - 1
            y2[y2_index] += y1[ibin]
            #print 'First bin %d counts go to y2 bin-%d' % (y1[ibin],y2_index)
            #print y2 
        if i_place_upper==2:
            print 'first bin 2 overlaps'

    # redistributing X1 bins with [1,m-1] indeces
    for ibin in range(1,m-1):
        if y1[ibin]==0:
            continue
        i_place_lower=i_place[ibin]
        i_place_upper=i_place[ibin+1]
        x1min = x1[ibin]
        x1max = x1[ibin+1]
        #Stop at the end
        if i_place_lower >= n-2:
            return y2
        
        # X1 bin fully inside X2 bin:
        if i_place_lower == i_place_upper:
            y2_index = i_place_upper - 1
            #print 'y2_index %d -- i_place_lower=%d -- ibin=%d' % (y2_index,i_place_lower, ibin)
            y2[y2_index] += y1[ibin]
            #print '%d-th bin %d counts go to y2 bin-%d: X1 overlaps with X2' % (ibin,y1[ibin],y2_index)
            #print y2
        # X1 bin overlaps w/ two X2 bins:
        if i_place_lower == i_place_upper - 1:
            # X2 value that "splits" the X1 bin
            x2_ov1 = x2[i_place_lower]
            # Roll the dice y1[ibin]-times
            for irand in range(0, int(y1[ibin])):
                probValue = np.random.uniform(x1min,x1max)
                #print 'rand-%d probV=%f x2_ov1=%f' %(irand,probValue,x2_ov1)
                if probValue < x2_ov1:
                    # send photon to lower bin :))
                    y2_index = i_place_upper - 2
                    y2[y2_index] += 1
                else:
                    y2_index = i_place_upper - 1
                    y2[y2_index] += 1
            #print '%d-th bin %d counts go to y2 bin-%d: X1 is completely in X2' % (ibin,y1[ibin],y2_index)
            #print y2
        # X1 bin overplaps w/ three X2 bins
        if i_place_lower == i_place_upper - 2:
            print '2 overlap bins'
        
    return y2


class algorithms:
    
    def __init__(self,foreTimeWindowSize,bkgTimeWindowSize,SPRTupperThresholds,SPRTlowerThresholds,kSigmaThresholds,isotopeList):
        
        # define ze variables
        self.foreTimes=[]
        self.bkgTimes=[]
        self.foreground=[]
        self.background=[]
        self.metrics=np.zeros(len(isotopeList))
        self.alarm=False
        self.alarms=[]
        self.ID=''
        self.IDs=[]
        self.aggSpec=np.array([])
        self.aggLiveTime=0.0
        
        # these would go in a configuration file
        self.foreTimeWindowSize=foreTimeWindowSize
        self.bkgTimeWindowSize=bkgTimeWindowSize
        self.isotopeList=isotopeList
        # self.a=np.log10((1-beta)/alpha)
        # self.b=np.log10(beta/(1-alpha))
        self.SPRTupperThresholds=SPRTupperThresholds
        self.SPRTlowerThresholds=SPRTlowerThresholds
        self.kSigmaThresholds=kSigmaThresholds
        
        self.startAlg=False

    def doTimeWindows(self,liveTime,spectrum):

        # print liveTime

        # initialize/generate foreground and background windows
        if len(self.foreTimes) == 0:
            self.foreground.append(spectrum)
            self.foreTimes.append(liveTime)
        else:
            self.foreground.append(spectrum)
            self.foreTimes.append(liveTime)
        if np.around(np.sum(self.foreTimes)) > self.foreTimeWindowSize:
            recordBkg=True
            if len(self.bkgTimes) == 0:
                self.background.append(self.foreground[0])
                self.bkgTimes.append(self.foreTimes[0])
                del self.foreTimes[0]
                del self.foreground[0]
            elif np.around(np.sum(self.bkgTimes)) < self.bkgTimeWindowSize:
                self.background.append(self.foreground[0])
                self.bkgTimes.append(self.foreTimes[0])
                del self.foreTimes[0]
                del self.foreground[0]
            else:
                # if not self.alarm:
                #    self.background.append(self.foreground[0])
                #    self.bkgTimes.append(self.foreTimes[0])
                #    del self.bkgTimes[0]
                #    del self.background[0]
                self.background.append(self.foreground[0])
                self.bkgTimes.append(self.foreTimes[0])
                del self.bkgTimes[0]
                del self.background[0]
                del self.foreTimes[0]
                del self.foreground[0]
                # print 'print fore comp',np.sum(self.foreTimes),self.foreTimeWindowSize
                # print 'print bkg comp',np.sum(self.bkgTimes),self.bkgTimeWindowSize
                self.startAlg=True
        # need this so our time windows don't get out of wack
        if np.around(np.sum(self.foreTimes)) > self.foreTimeWindowSize:
            while True:
                del self.foreground[0]
                del self.foreTimes[0]
                if np.around(np.sum(self.foreTimes)) <= self.foreTimeWindowSize:
                    break
        if np.around(np.sum(self.bkgTimes)) > self.bkgTimeWindowSize:
            while True:
                del self.background[0]
                del self.bkgTimes[0]
                if np.around(np.sum(self.bkgTimes)) <= self.bkgTimeWindowSize:
                    break
                
            return
        
    def doSPRT(self,liveTime,energies,spectrum):
        #pull out spectrum region of interest
        spectrum=np.array(spectrum)

        recordBkg=False

        self.doTimeWindows(liveTime,spectrum)
        bkgcnts=np.sum(self.background,axis=0)
        forecnts=np.sum(self.foreground,axis=0)

        self.IDs=[]
        self.alarms=[]
        #self.alarm=False
        #self.ID=''
        if self.startAlg:
            for j in range(len(self.isotopeList)):
                indexs=np.where((energies >= self.isotopeList[j]['ROI_Mins'][0]) & (energies <= self.isotopeList[j]['ROI_Maxs'][0]))[0]
                bkg=bkgcnts[indexs]
                fore=forecnts[indexs]
                bkgCR=np.sum(bkg)/np.sum(self.bkgTimes)
                foreCR=np.sum(fore)/np.sum(self.foreTimes)
                #liklihood ratio test for a Poisson variable
                thisMetric=(foreCR-bkgCR)-foreCR*np.log10(foreCR/bkgCR) # count rate is not quite Poisson... throwing caution to the wind
                if np.isnan(thisMetric):
                    thisMetric=0.0
                self.metrics[j]=self.metrics[j]+thisMetric
                if self.metrics[j] <= self.SPRTlowerThresholds[j]:
                    self.metrics[j]=0.0
                if self.metrics[j] >= self.SPRTupperThresholds[j]:
                    self.alarms.append(True)
                    self.IDs.append(self.isotopeList[j])
                else:
                    self.alarms.append(False)
                #print j,thisMetric,self.metrics[j]
                #print '\t',foreCR,bkgCR,np.sum(self.bkgTimes),np.sum(self.foreTimes)
        if True in self.alarms:
            #print alarms
            self.alarm=True
            #use ROI with maximum energy for ROI ID
            maxEnergy=0
            maxIndex=0
            for i in range(len(self.IDs)):
                if self.IDs[i]['Energies'][0] > maxEnergy:
                    maxEnergy=self.IDs[i]['Energies'][0]
                    maxIndex=i
            self.ID=self.IDs[maxIndex]
            if len(self.aggSpec) == 0:
                self.aggSpec=spectrum
                self.aggLiveTime=liveTime
            else:
                self.aggSpec=self.aggSpec+spectrum
                self.aggLiveTime=self.aggLiveTime+liveTime
        else:
            self.ID=''
            self.aggSpec=np.array([])
            self.aggLiveTime=0.0
            self.alarm=False

        return

    def dokSigma(self,liveTime,energies,spectrum):

        #pull out spectrum region of interest
        spectrum=np.array(spectrum)

        recordBkg=False

        self.doTimeWindows(liveTime,spectrum)
        bkgcnts=np.sum(self.background,axis=0)
        forecnts=np.sum(self.foreground,axis=0)

        self.alarms=[]
        if self.startAlg:
            for j in range(len(self.isotopeList)):
                indexs=np.where((energies >= self.isotopeList[j]['ROI_Mins'][0]) & (energies <= self.isotopeList[j]['ROI_Maxs'][0]))[0]
                bkg=bkgcnts[indexs]
                fore=forecnts[indexs]
                bkgCR=np.sum(bkg)/np.sum(self.bkgTimes)
                foreCR=np.sum(fore)/np.sum(self.foreTimes)
                #liklihood ratio test for a Poisson variable
                thisMetric=(foreCR-bkgCR)/(np.sqrt(np.sum(bkg))/np.sum(self.bkgTimes))
                if np.isnan(thisMetric):
                    thisMetric=0.0
                self.metrics[j]=thisMetric
                if self.metrics[j] >= self.kSigmaThresholds[j]:
                    self.alarms.append(True)
                else:
                    self.alarms.append(False)
                #print j,thisMetric,self.metrics[j]
                #print '\t',foreCR,bkgCR,np.sum(self.bkgTimes),np.sum(self.foreTimes)
        if True in self.alarms:
            #print alarms
            self.alarm=True
            if len(self.aggSpec) == 0:
                self.aggSpec=spectrum
                self.aggLiveTime=liveTime
            else:
                self.aggSpec=self.aggSpec+spectrum
                self.aggLiveTime=self.aggLiveTime+liveTime
        else:
            self.aggSpec=np.array([])
            self.aggLiveTime=0.0
            self.alarm=False

        return


### some definitions
xfmt = mdates.DateFormatter('%H:%M:%S')
dates=np.vectorize(dt.datetime.fromtimestamp)

energies=np.linspace(0,3000,1000)

plotIt=False

### look at declared transfer events
f=open('transfers_03022018.txt','r')
transferTimes=[]
icount=0
for line in f:
    if icount > 0:
        line=line.split('\t')
        dum=line[1].split('/')
        if len(dum[0]) == 1:
            dum[0]='0'+dum[0]
        if len(dum[1]) == 1:
            dum[1]+'0'+dum[1]
        dum='/'.join(dum)
        dum=dt.datetime.strptime(dum,'%m/%d/%y %H:%M')
        dum=time.mktime(dum.timetuple())
        transferTimes.append(dum)
    icount=icount+1

f.close()

transferTimes=np.array(transferTimes)

##############

#### read algorithm data
hdf5Path='hdf5Files/'
#hdf5Path='hdf5Files/08142018/'
hdf5Files=os.listdir(hdf5Path)
hdf5Files=np.array([x for x in hdf5Files if '.hdf5' in x])
hdf5Files=np.array([x for x in hdf5Files if '_ALG' in x])
dums=[x.split('MUSE01-')[-1].split('_ALG.hdf5')[0] for x in hdf5Files]
fmt='%Y-%m-%dT%H.%M.%S.%f'
hdf5TimeStamps=np.array([time.mktime(dt.datetime.strptime(x,fmt).timetuple()) for x in dums])
indexs=np.argsort(hdf5TimeStamps)
hdf5Files=hdf5Files[indexs]
hdf5TimeStamps=hdf5TimeStamps[indexs]

sprtMetrics=[]
kSigmaMetrics=[]
countRates=[]
allTimes=[]
alarmSpectra=[]
alarmLiveTimes=[]
alarmTimes=[]
alarmLiveTimes=[]
for hdf5File in hdf5Files:
    f=h5py.File(hdf5Path+hdf5File,"r")
    sprtMetrics=sprtMetrics+f['sprtMetrics'][:].tolist()
    kSigmaMetrics=kSigmaMetrics+f['kSigmaMetrics'][:].tolist()
    countRates=countRates+f['countRates'][:].tolist()
    allTimes=allTimes+f['times'][:].tolist()
    
    alarmSpectra=alarmSpectra+f['alarmSpectra'][:].tolist()
    alarmLiveTimes=alarmLiveTimes+f['alarmLiveTimes'][:].tolist()
    alarmTimes=alarmTimes+f['alarmTimes'][:].tolist()
    f.close()

sprtMetrics=np.array(sprtMetrics)
kSigmaMetrics=np.array(kSigmaMetrics)
countRates=np.array(countRates)
allTimes=np.array(allTimes)
alarmSpectra=np.array(alarmSpectra)
alarmLiveTimes=np.array(alarmLiveTimes)
alarmTimes=np.array(alarmTimes)
alarmLiveTimes=np.array(alarmLiveTimes)

#sum alarms close together in time
sumedAlarmTimes=[alarmTimes[0]]
sumedAlarmSpectra=[alarmSpectra[0]]
sumedAlarmLiveTimes=[alarmLiveTimes[0]]
sumedAlarmStart=alarmTimes[0][0]
sumedAlarmEnd=alarmTimes[0][1]
#sumedroiIDs=[roiIDs[0]]

icount=0
for i in range(1,len(alarmTimes)):
    #if np.abs(alarmTimes[i][0]-alarmTimes[i-1][0]) < 300.0 and roiIDs[i] != roiIDs[-1]:
    if np.abs(alarmTimes[i][0]-alarmTimes[i-1][1]) < 300.0:
        sumedAlarmSpectra[icount]=sumedAlarmSpectra[icount]+alarmSpectra[i]
        sumedAlarmLiveTimes[icount]=sumedAlarmLiveTimes[icount]+alarmLiveTimes[i]
        sumedAlarmTimes[icount][1]=alarmTimes[i][1]
    else:
        sumedAlarmTimes.append(alarmTimes[i])
        sumedAlarmSpectra.append(alarmSpectra[i])
        sumedAlarmLiveTimes.append(alarmLiveTimes[i])
        icount=icount+1
    #if sumedAlarmTimes[-1][1] < sumedAlarmTimes[-1][0]:
    #    print lknlkn

sumedAlarmSpectra=np.array(sumedAlarmSpectra)
sumedAlarmLiveTimes=np.array(sumedAlarmLiveTimes)
sumedAlarmTimes=np.array(sumedAlarmTimes)


#get backgrounds

hdf5Path='hdf5Files/'
#hdf5Path='hdf5Files/08142018/'
hdf5Files=os.listdir(hdf5Path)
hdf5Files=np.array([x for x in hdf5Files if '.hdf5' in x])
hdf5Files=np.array([x for x in hdf5Files if '_ALG' not in x])
dums=[x.split('GuardShack-')[-1].split('.hdf5')[0] for x in hdf5Files]
dums=[x.split('Supernode-')[-1].split('.hdf5')[0] for x in dums]
dums=[x.split('MUSE01-')[-1].split('.hdf5')[0] for x in dums]
fmt='%Y-%m-%dT%H.%M.%S.%f'
hdf5TimeStamps=np.array([time.mktime(dt.datetime.strptime(x,fmt).timetuple()) for x in dums])
indexs=np.argsort(hdf5TimeStamps)
hdf5Files=hdf5Files[indexs]
hdf5TimeStamps=hdf5TimeStamps[indexs]

intTime=600.0
deltaT=60.0

bkgSpecs=[]
bkgLTs=[]
bkgTimes=[]
for i in range(len(sumedAlarmTimes)):
    indexs=np.where( (hdf5TimeStamps < sumedAlarmTimes[i][0]) )[0]
    jndexs=np.where( (hdf5TimeStamps > sumedAlarmTimes[i][1]) ) [0]
    if len(jndexs) == 0:
        j=-1
    else:
        kndexs=range(indexs[-1],jndexs[0]+1)
        j=kndexs[0]
    try:
        f=h5py.File(hdf5Path+hdf5Files[j])
        times=f['2x4x16Times'][:]
        liveTimes=f['2x4x16LiveTimes'][:]
        spectra=f['2x4x16Spectra'][:]
        f.close()
    except:
        print 'file load abort!',hdf5Files[j],dates(sumedAlarmTimes[i])
        bkgTimes.append([sumedAlarmTimes[i][0]-deltaT-intTime,sumedAlarmTimes[i][0]-deltaT])
        bkgSpecs.append(np.zeros(len(energies)))
        bkgLTs.append(0.0)
        continue
    indexs=np.where( (times >= sumedAlarmTimes[i][0]-deltaT-intTime) & (times <= sumedAlarmTimes[i][0]-deltaT) )[0]
    if len(indexs) > 0:
        bkgTimes.append([sumedAlarmTimes[i][0]-deltaT-intTime,sumedAlarmTimes[i][0]-deltaT])
        bkgSpecs.append(np.sum(spectra[indexs],axis=0))
        bkgLTs.append(np.sum(liveTimes[indexs]))
    else:
        #try to go back one file
        j=j-1
        try:
            f=h5py.File(hdf5Path+hdf5Files[j])
            times=f['2x4x16Times'][:]
            liveTimes=f['2x4x16LiveTimes'][:]
            spectra=f['2x4x16Spectra'][:]
            f.close()
            indexs=np.where( (times >= sumedAlarmTimes[i][0]-deltaT-intTime) & (times <= sumedAlarmTimes[i][0]-deltaT) )[0]
            if len(indexs) > 0:
                bkgTimes.append([sumedAlarmTimes[i][0]-deltaT-intTime,sumedAlarmTimes[i][0]-deltaT])
                bkgSpecs.append(np.sum(spectra[indexs],axis=0))
                bkgLTs.append(np.sum(liveTimes[indexs]))
            else:
                print lknlkn # go to exception
        except:
            bkgTimes.append([sumedAlarmTimes[i][0]-deltaT-intTime,sumedAlarmTimes[i][0]-deltaT])
            bkgSpecs.append(np.zeros(len(energies)))
            bkgLTs.append(0.0)
            print 'zeroed spectrum',hdf5Files[j+1],dates(sumedAlarmTimes[i])

bkgSpecs1=np.array(bkgSpecs)
bkgLTs1=np.array(bkgLTs)
bkgTimes1=np.array(bkgTimes)

#

intTime=600.0
deltaT=960.0

bkgSpecs=[]
bkgLTs=[]
bkgTimes=[]
bkgEndTimes=[]
for i in range(len(sumedAlarmTimes)):
    indexs=np.where( (hdf5TimeStamps < sumedAlarmTimes[i][0]) )[0]
    jndexs=np.where( (hdf5TimeStamps > sumedAlarmTimes[i][1]) ) [0]
    if len(jndexs) == 0:
        j=-1
    else:
        kndexs=range(indexs[-1],jndexs[0]+1)
        j=kndexs[0]
    try:
        f=h5py.File(hdf5Path+hdf5Files[j])
        times=f['2x4x16Times'][:]
        liveTimes=f['2x4x16LiveTimes'][:]
        spectra=f['2x4x16Spectra'][:]
        f.close()
    except:
        print 'file load abort!',hdf5Files[j],dates(sumedAlarmTimes[i])
        bkgTimes.append([sumedAlarmTimes[i][0]-deltaT-intTime,sumedAlarmTimes[i][0]-deltaT])
        bkgSpecs.append(np.zeros(len(energies)))
        bkgLTs.append(0.0)
        continue
    indexs=np.where( (times >= sumedAlarmTimes[i][0]-deltaT-intTime) & (times <= sumedAlarmTimes[i][0]-deltaT) )[0]
    if len(indexs) > 0:
        bkgTimes.append([sumedAlarmTimes[i][0]-deltaT-intTime,sumedAlarmTimes[i][0]-deltaT])
        bkgSpecs.append(np.sum(spectra[indexs],axis=0))
        bkgLTs.append(np.sum(liveTimes[indexs]))
    else:
        #try to go back one file
        j=j-1
        try:
            f=h5py.File(hdf5Path+hdf5Files[j])
            times=f['2x4x16Times'][:]
            liveTimes=f['2x4x16LiveTimes'][:]
            spectra=f['2x4x16Spectra'][:]
            f.close()
            indexs=np.where( (times >= sumedAlarmTimes[i][0]-deltaT-intTime) & (times <= sumedAlarmTimes[i][0]-deltaT) )[0]
            if len(indexs) > 0:
                bkgTimes.append([sumedAlarmTimes[i][0]-deltaT-intTime,sumedAlarmTimes[i][0]-deltaT])
                bkgSpecs.append(np.sum(spectra[indexs],axis=0))
                bkgLTs.append(np.sum(liveTimes[indexs]))
            else:
                print lknlkn # go to exception
        except:
            bkgTimes.append([sumedAlarmTimes[i][0]-deltaT-intTime,sumedAlarmTimes[i][0]-deltaT])
            bkgSpecs.append(np.zeros(len(energies)))
            bkgLTs.append(0.0)
            print 'zeroed spectrum',hdf5Files[j+1],dates(sumedAlarmTimes[i])

bkgSpecs2=np.array(bkgSpecs)
bkgLTs2=np.array(bkgLTs)
bkgTimes2=np.array(bkgTimes)

f=h5py.File('alarmData.hdf5',"w")
dset=f.create_dataset('2x4x16AlarmSpectra',(len(sumedAlarmSpectra),len(sumedAlarmSpectra[0])),dtype='f', compression="gzip")
dset[...]=sumedAlarmSpectra
dset=f.create_dataset('2x4x16AlarmLiveTimes',(len(sumedAlarmLiveTimes),),dtype='f', compression="gzip")
dset[...]=sumedAlarmLiveTimes
dset=f.create_dataset('2x4x16AlarmTimes',(len(sumedAlarmTimes),len(sumedAlarmTimes[0])),dtype='float64',compression="gzip")
dset[...]=sumedAlarmTimes
#
dset=f.create_dataset('2x4x16BkgSpectra1',(len(bkgSpecs1),len(bkgSpecs1[0])),dtype='f', compression="gzip")
dset[...]=bkgSpecs1
dset=f.create_dataset('2x4x16BkgLiveTimes1',(len(bkgLTs1),),dtype='f', compression="gzip")
dset[...]=bkgLTs1
dset=f.create_dataset('2x4x16BkgTimes1',(len(bkgTimes1),len(bkgTimes1[0])),dtype='float64', compression="gzip")
dset[...]=bkgTimes1
#
dset=f.create_dataset('2x4x16BkgSpectra2',(len(bkgSpecs2),len(bkgSpecs2[0])),dtype='f', compression="gzip")
dset[...]=bkgSpecs2
dset=f.create_dataset('2x4x16BkgLiveTimes2',(len(bkgLTs2),),dtype='f', compression="gzip")
dset[...]=bkgLTs2
dset=f.create_dataset('2x4x16BkgTimes2',(len(bkgTimes2),len(bkgTimes2[0])),dtype='float64', compression="gzip")
dset[...]=bkgTimes2
f.close()


for i in range(len(sumedAlarmSpectra)):
    ttimes=dates(sumedAlarmTimes[i])
    start=ttimes[0]
    end=ttimes[1]
    startString='%s/%d %s:%s:%s - ' % (start.month,start.day,start.hour,start.minute,start.second)
    endString='%s/%d %s:%s:%s - ' % (end.month,end.day,end.hour,end.minute,end.second)
    print '%s - %s' % (start,end)
    #make a plot
    plt.figure()
    plt.step(energies,sumedAlarmSpectra[i]/sumedAlarmLiveTimes[i],lw=2,where='mid',label='Fore.')
    plt.step(energies,bkgSpecs1[i]/bkgLTs1[i],lw=2,where='mid',label='Bkg')
    plt.step(energies,(sumedAlarmSpectra[i]/sumedAlarmLiveTimes[i]-bkgSpecs1[i]/bkgLTs1[i]).clip(min=0),lw=2,where='mid',label='Diff')
    plt.title('%s - %s' % (start,end))
    plt.grid(True)
    plt.xlabel('Energy (keV)', size=14)
    plt.ylabel('Count Rate',size=14)
    plt.legend()
    plt.yscale('log')
    plt.savefig('figs/alarmFigs/%s-%s-%s_%s-%s-%s.png' % (start.strftime('%m'),start.strftime('%d'),start.strftime('%Y'),start.strftime('%H'),start.strftime('%M'),start.strftime('%S')),format='png')
    plt.show(False)
    #plt.close()
    #write an N42
    templateFile='n42Examplev7.n42'
    f=open(templateFile,'r')
    contents=f.read()
    f.close()
    #
    # Foreground
    #
    dateString='%s-%s-%sT%s:%s:%s.000Z' % (start.strftime('%Y'),start.strftime('%m'),start.strftime('%d'),start.strftime('%H'),start.strftime('%M'),start.strftime('%S'))
    contents=contents.replace('$$$$foreDate$$$$',dateString)
    realTime=sumedAlarmTimes[i][-1]-sumedAlarmTimes[i][0]
    contents=contents.replace('$$$$foreRealTime$$$$','PT%.3fS' % realTime)
    liveTime=sumedAlarmLiveTimes[i]
    contents=contents.replace('$$$$foreLiveTime$$$$','PT%.3fS' % liveTime)
    specString=' '.join([str(int(x)) for x in sumedAlarmSpectra[i].tolist()])
    contents=contents.replace('$$$$foreSpectrum$$$$',specString)
    #
    # Background 1
    #
    ttimes=dates(bkgTimes1[i])
    start=ttimes[0]
    dateString='%s-%s-%sT%s:%s:%s.000Z' % (start.strftime('%Y'),start.strftime('%m'),start.strftime('%d'),start.strftime('%H'),start.strftime('%M'),start.strftime('%S'))
    contents=contents.replace('$$$$bkgDate1$$$$',dateString)
    realTime=bkgTimes1[i][1]-bkgTimes1[i][0]
    contents=contents.replace('$$$$bkgRealTime1$$$$','PT%.3fS' % realTime)
    liveTime=bkgLTs1[i]
    contents=contents.replace('$$$$bkgLiveTime1$$$$','PT%.3fS' % liveTime)
    specString=' '.join([str(int(x)) for x in bkgSpecs1[i].tolist()])
    contents=contents.replace('$$$$bkgSpectrum1$$$$',specString)
    #
    # Background 2
    #
    ttimes=dates(bkgTimes2[i])
    start=ttimes[0]
    dateString='%s-%s-%sT%s:%s:%s.000Z' % (start.strftime('%Y'),start.strftime('%m'),start.strftime('%d'),start.strftime('%H'),start.strftime('%M'),start.strftime('%S'))
    contents=contents.replace('$$$$bkgDate2$$$$',dateString)
    realTime=bkgTimes2[i][1]-bkgTimes2[i][0]
    contents=contents.replace('$$$$bkgRealTime2$$$$','PT%.3fS' % realTime)
    liveTime=bkgLTs2[i]
    contents=contents.replace('$$$$bkgLiveTime2$$$$','PT%.3fS' % liveTime)
    specString=' '.join([str(int(x)) for x in bkgSpecs2[i].tolist()])
    contents=contents.replace('$$$$bkgSpectrum2$$$$',specString)
    #
    ttimes=dates(sumedAlarmTimes[i])
    start=ttimes[0]
    outFile='figs/alarmFigs/%s-%s-%s_%s-%s-%s.n42' % (start.strftime('%m'),start.strftime('%d'),start.strftime('%Y'),start.strftime('%H'),start.strftime('%M'),start.strftime('%S'))
    shutil.copyfile(templateFile, outFile)
    f=open(outFile,'w')
    f.write(contents)
    f.close()


#attempt to identify any peaks
idEns=np.array([609.0,662.0,1294.0,1173.0])
ids=['Rain','Cs-137','Ar-41','Co-60']

#load res func params
resParams=np.load('resParams_superNode_2x4x16.npy')

f=open('alarmInfo.csv','w')
f.write('Date,Start Time,End Time,Start Time (posix),End Time (posix),Count Rate, IDs\n')
for i in range(len(sumedAlarmTimes)):
    ttimes=dates(sumedAlarmTimes[i])
    start=ttimes[0]
    end=ttimes[1]
    date='%s/%s/%s' % (start.strftime('%m'),start.strftime('%d'),start.strftime('%Y'))
    startTime='%s:%s:%s' % (start.strftime('%H'),start.strftime('%M'),start.strftime('%S'))
    endTime='%s:%s:%s' % (end.strftime('%H'),end.strftime('%M'),end.strftime('%S'))
    plt.figure()
    plt.step(energies,sumedAlarmSpectra[i]/sumedAlarmLiveTimes[i],lw=2,label='Alarm')
    plt.step(energies,bkgSpecs1[i]/bkgLTs1[i],lw=2,label='Bkg 1')
    plt.step(energies,bkgSpecs2[i]/bkgLTs2[i],lw=2,label='Bkg 2')
    if np.sum(bkgLTs1) > 0:
        diff=sumedAlarmSpectra[i]/sumedAlarmLiveTimes[i]-bkgSpecs1[i]/bkgLTs1[i]
    else:
        diff=np.zeros(1000)
    plt.step(energies,diff,lw=2,label='Diff')
    sumedBins=np.sum(diff.reshape(-1, 4), axis=1)
    sumedEns=12.0*(np.arange(len(sumedBins))+1)
    plt.step(sumedEns,sumedBins,lw=2)
    peakind = signal.find_peaks_cwt(sumedBins, [5,10,2])
    #look for good peaks
    peakEns=[]
    for j in range(len(peakind)):
        plt.plot(sumedEns[peakind[j]],sumedBins[peakind[j]],marker='*',ls='',color='k')
        spec=diff
        sigma=resFunc(sumedEns[peakind[j]],resParams[0],resParams[1],resParams[2])
        eMin=sumedEns[peakind[j]]-5*sigma
        eMax=sumedEns[peakind[j]]+5*sigma
        indexs=np.where( (energies >= eMin) & (energies <= eMax) )[0]
        xs=energies[indexs]
        ys=spec[indexs]
        a=np.min(ys)
        b=np.max(ys)/(sigma*np.sqrt(2.0*np.pi))
        c=(ys[-1]-ys[0])/(eMax-eMin)
        mu=sumedEns[peakind[j]]
        goodPeak=False
        try:
            popt,pcov=curve_fit(gausswLine,xs,ys,p0=[a,b,c,mu,sigma],bounds=([-0.3,0.0,-1.0,eMin,sigma-sigma*.2],[1E9,1E9,1.0,eMax,sigma+sigma*0.2]))
            if np.sqrt(pcov[-1][-1])/(popt[-1]) < 0.15 and np.sqrt(pcov[1][1])/(popt[1]) < 2.0:
                goodPeak=True
                print (pcov[-1][-1])/(popt[-1]),'>',0.15
            if goodPeak:
                print
                print sumedEns[peakind[j]]
                print '\t',np.sqrt(pcov[-1][-1])/(popt[-1])
                peakEns.append(popt[3]-4.0)
                plt.plot(xs,gausswLine(xs,*popt),lw=2,marker='',ls='-',color='k')
                for k in range(len(popt)):
                    print '\t',popt[k],np.sqrt(np.diag(pcov)[k])
        except:
            #print '\t',a,b,c,mu,sigma
            pass
    plt.legend()
    plt.xlabel('Energy (keV)',size=14)
    plt.xlabel('Count Rate',size=14)
    plt.grid(True)
    plt.yscale('log')
    plt.savefig('figs/alarmFigs/%s-%s-%s_%s-%s-%s_peaks.png' % (start.strftime('%m'),start.strftime('%d'),start.strftime('%Y'),start.strftime('%H'),start.strftime('%M'),start.strftime('%S')),format='png')
    plt.show(False)
    #plt.close()
    #assign ids to peaks
    alarmIDs=[]
    for j in range(len(peakEns)):
        eDiff=np.abs((peakEns[j]-idEns))/peakEns[j]
        index=np.argmin(eDiff)
        if eDiff[index] < 0.05:
            alarmIDs.append(ids[index])
            break
    string='%s,%s,%s,%d,%d,%.2f,%s\n' % (date,startTime,endTime,int(time.mktime(start.timetuple())),int(time.mktime(end.timetuple())),np.sum(diff),' '.join(alarmIDs))
    f.write(string)
    #print lknlkn

f.close()

########### get video data
    
#camera facing out of node
videoPath='/rawdata/MINOS/MINOS_Video/superNode1/'
videoFiles=[]
videoStartTimes=[]
videoEndTimes=[]
for root, subdirs, files in os.walk(videoPath):
    #print root,files
    files=[x for x in files if '.mp4' in x]
    if len(files) > 0:
        videoFiles=videoFiles+[root+'/'+x for x in files]
        videoStartTimes=videoStartTimes+[x.split('/')[-1].split('_')[0] for x in files]
        videoEndTimes=videoEndTimes+[x.split('/')[-1].split('_')[1] for x in files]
        
videoFiles=np.array(videoFiles)
videoStartTimes=np.array(videoStartTimes)
videoEndTimes=np.array(videoEndTimes)
indexs=np.argsort(videoStartTimes)
videoFiles=videoFiles[indexs]
videoStartTimes=videoStartTimes[indexs].astype('int')
videoEndTimes=videoEndTimes[indexs].astype('int')
	

#load SPRT alarms
f=h5py.File('alarmData.hdf5',"r")
alarmSpectra=f['2x4x16AlarmSpectra'][:]
alarmLiveTimes=f['2x4x16AlarmLiveTimes'][:]
alarmTimes=f['2x4x16AlarmTimes'][:]
f.close()

#between September 9th and December 14th, video was 55 seconds ahead of radiation
startError=1504929600*1000
endError=1513314000*1000

basePath='/rawdata/MINOS/MINOS_Video/superNode1/'
os.chdir(basePath)
for i in range(len(alarmTimes)):
    start=int(alarmTimes[i][0]*1000)-90*1000 # get 90 seconds before the alarm start
    end=int(alarmTimes[i][1]*1000)+90*1000 # get 90 seconds seconds after the alarm end
    if (end-start)/1000.0 > 600: #limit to 10 minute videos
	end = start+600*1000
    indexs=np.where( (videoStartTimes >= start) & (videoEndTimes <= end) )[0]
    print start,end,len(indexs)
    if len(indexs) > 0: #we have a video
	string='concat:'
	for j in range(len(indexs)):
	    index=indexs[j]
	    if j > 0:
		string=string+'|' + videoFiles[index].split(basePath)[1]
	    else:
		string=string+videoFiles[index].split(basePath)[1]
	start=dt.datetime.fromtimestamp(start/1000.0)
	end=dt.datetime.fromtimestamp(end/1000.0)
	dateString='%s-%s-%s_%s-%s-%s_%s-%s-%s_%s-%s-%s' % (start.strftime('%m'),start.strftime('%d'),start.strftime('%Y'),start.strftime('%H'),start.strftime('%M'),start.strftime('%S'),end.strftime('%m'),end.strftime('%d'),end.strftime('%Y'),end.strftime('%H'),end.strftime('%M'),end.strftime('%S'))
	command=["ffmpeg","-i",string,"-c","copy","/home/xae/figs/alarmFigs/%s_Camera1.mp4" % dateString]
	subprocess.call(command)
        print dateString

os.chdir('/home/xae/')
#camera facing towards of node
videoPath='/rawdata/MINOS/MINOS_Video/superNode2/'
videoFiles=[]
videoStartTimes=[]
videoEndTimes=[]
for root, subdirs, files in os.walk(videoPath):
    #print root,files
    files=[x for x in files if '.mp4' in x]
    if len(files) > 0:
        videoFiles=videoFiles+[root+'/'+x for x in files]
        videoStartTimes=videoStartTimes+[x.split('/')[-1].split('_')[0] for x in files]
        videoEndTimes=videoEndTimes+[x.split('/')[-1].split('_')[1] for x in files]
videoFiles=np.array(videoFiles)
videoStartTimes=np.array(videoStartTimes)
videoEndTimes=np.array(videoEndTimes)
indexs=np.argsort(videoStartTimes)
videoFiles=videoFiles[indexs]
videoStartTimes=videoStartTimes[indexs].astype('int')
videoEndTimes=videoEndTimes[indexs].astype('int')

	
basePath='/rawdata/MINOS/MINOS_Video/superNode2/'
os.chdir(basePath)
for i in range(len(alarmTimes)):
    start=int(alarmTimes[i][0]*1000)-90*1000 # get 300 seconds before the alarm start
    end=int(alarmTimes[i][1]*1000)+90*1000 # get 90 seconds seconds after the alarm end
    if (end-start)/1000.0 > 600: #limit to 10 minute videos
	end = start+600*1000
    indexs=np.where( (videoStartTimes >= start) & (videoEndTimes <= end) )[0]
    print start,end,len(indexs)
    if len(indexs) > 0: #we have a video
	string='concat:'
	for j in range(len(indexs)):
	    index=indexs[j]
	    if j > 0:
		string=string+'|' + videoFiles[index].split(basePath)[1]
	    else:
		string=string+videoFiles[index].split(basePath)[1]
	start=dt.datetime.fromtimestamp(start/1000.0)
	end=dt.datetime.fromtimestamp(end/1000.0)
	dateString='%s-%s-%s_%s-%s-%s_%s-%s-%s_%s-%s-%s' % (start.strftime('%m'),start.strftime('%d'),start.strftime('%Y'),start.strftime('%H'),start.strftime('%M'),start.strftime('%S'),end.strftime('%m'),end.strftime('%d'),end.strftime('%Y'),end.strftime('%H'),end.strftime('%M'),end.strftime('%S'))
	command=["ffmpeg","-i",string,"-c","copy","/home/xae/figs/alarmFigs/%s_Camera2.mp4" % dateString]
	subprocess.call(command)

