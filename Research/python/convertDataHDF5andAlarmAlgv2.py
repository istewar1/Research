import sqlite3
import numpy as np
import matplotlib
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

energies=np.linspace(0,3000,1024)

plotIt=False
node = 'MUSE10'
inPath='/Volumes/Ian External HD/Node Data/sqlitefiles/'+node+'/'
eCal_path = '/Volumes/Ian External HD/Node Data/energy_pairs/'+node+'/'
hdf5outPath=inPath

#list all database files
dbFiles=os.listdir(inPath)
dbFiles=[x for x in dbFiles if 'sqlite3' in x]


binEdges=np.arange(0,3006,6)

energies3x3=np.load(eCal_path+'eCal_3x3.npy')
k40PeakChan=np.abs(energies3x3-1460.8).argmin()
slope=k40PeakChan/1460.8
nonLinearity=range(len(energies3x3))/linFunc(energies3x3,slope)

count = 0  # Use
count1 = 0 # Use
eCal3x3=np.linspace(0,3000,1024)

for dbFile in dbFiles:
    ### intialize algorithms
    isotopeList = [
        {'Isotope': 'Total Spectrum', 'Energies': [0.0], 'ROI_Mins': [0.0], 'ROI_Maxs': [3000.0], 'Fit_Mins': [0.0],
         'Fit_Maxs': [3000.0], 'Sigmas': [0.0]}]

    foreTimeWindowSize = 2.0
    bkgTimeWindowSize = 30.0
    upperThresholds2x4x16 = [1000.0]
    lowerThresholds2x4x16 = np.zeros(len(isotopeList))
    kSigmaThresholds = [20.0]

    sprt = algorithms(foreTimeWindowSize, bkgTimeWindowSize, upperThresholds2x4x16, lowerThresholds2x4x16,
                      kSigmaThresholds, isotopeList)
    kSigma = algorithms(foreTimeWindowSize, bkgTimeWindowSize, upperThresholds2x4x16, lowerThresholds2x4x16,
                        kSigmaThresholds, isotopeList)

    icount = 0
    startPeakLocs = []
    initSpectra = []

    # algorithm arrays
    alarms = []
    sprtMetrics = []
    kSigmaMetrics = []

    alarmSpec = []
    alarmSpectra = []
    alarmLiveTimes = []
    alarmTimes = []
    alarmIndexs = []
    alarmStartTime = 0

    countRates = []
    dumIDs = []
    allTimes = []

    try:
        print 'Reading',dbFile
        outFile=dbFile.replace('.sqlite3','.hdf5')
        outFile=node+'-'+'-'.join(outFile.split('-')[1:])
        hdf5write=False
    except:
        continue

    detDB=database(inPath+dbFile)
    detDB.getTables()
    tables=[str(x[0]) for x in detDB.tables]

    dataTable='Det_3x3_data'
    if dataTable in tables:
        detDB.getHeader(dataTable)
        header=np.array(detDB.header)
        detDB.getColumn(dataTable,'Live_Time')
        liveTimes=np.array(detDB.data)[:,0]
        length=10*(len(liveTimes)/10)
        detDB.getColumn(dataTable,'Time')
        times=np.array(detDB.data)[:,0]
        detDB.getColumn(dataTable,'Spectrum__IntArray')
        spectra=[np.fromstring(x[0],sep=',',dtype='int') for x in detDB.data]
        indexs=np.where(liveTimes < 0.05)[0]
        if len(indexs) > 0:
            f4 = open('liveTimeIssues.txt','a+')
            for index in indexs:
                f4.write('3x3,%s,%.2f,%.6e,%d\n' % (dbFile,times[index],liveTimes[index],np.sum(spectra,axis=1)[index]))
            f4.close()
        #sum every 10 elements
        liveTimes=np.add.reduceat(liveTimes[:length],np.arange(0,length,10))
        spectra=np.add.reduceat(spectra[:length],np.arange(0,length,10))
        times=np.array([times[i+5] for i in range(0,length,10)])
        ###
        # Run algorithms, calibrate and rebin spectra
        ###
        windowSize=900.0
        windowSpec=[]
        windowLiveTime=[]
        rebinedSpectra=[]
        peakLocs=[]
        alarm=False
        plt.figure()
        summed = np.sum(spectra[800:950],axis=0)
        plt.plot(summed);plt.grid(alpha=0.5)
        plt.yscale('log')
        plt.show(False); plt.close('all')
        for i in range(len(spectra)):
            liveTime=liveTimes[i]
            spectrum=spectra[i]
            time=times[i]
            #print i,len(allTimes),len(countRates),len(sprtMetrics)
            allTimes.append(time)
            ###
            # Calibrate and Rebin
            ###
            if alarm:
                windowSpec=[]
                windowLiveTime=[]
            else:
                windowSpec.append(spectrum)
                windowLiveTime.append(liveTime)
            #if we have enough data, then recalibrate energy
            if np.around(np.sum(windowLiveTime) > windowSize):
                ### recalibrate
                spec=np.sum(windowSpec,axis=0)/np.sum(windowLiveTime)
                chanMin=425
                chanMax=545
                xs=np.arange(chanMin,chanMax)
                ys=spec[chanMin:chanMax]
                sigma=30.0
                a=np.min(ys)
                b=np.max(ys)/(sigma*np.sqrt(2.0*np.pi))
                c=(ys[-1]-ys[0])/(chanMax-chanMin)
                mu=xs[np.argmax(ys)]
                #sometimes the fit doesn't work, so use a try statement
                try:
                    popt,pcov=curve_fit(gausswLine,xs,ys,p0=[a,b,c,mu,sigma],bounds=([0.0,0.0,-1.0,chanMin,5.0],[1E9,1E9,0.0,chanMax,50.0]))
                    if plotIt:
                        print dates([times[i]])[0]
                        print '%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (popt[0],popt[1],popt[2],popt[3],popt[4])
                        print '%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (np.sqrt(pcov[0][0]),np.sqrt(pcov[1][1]),np.sqrt(pcov[2][2]),np.sqrt(pcov[3][3]),np.sqrt(pcov[4][4]))
                        print
                except:
                    print i
                    popt=peakLocs[-1][1]
                    print '\t peak fit failed'
                peakLocs.append([times[i],popt[3],np.sqrt(pcov[3][3])])
                eCal=doEcal(popt[3],energies3x3,nonLinearity)
                windowSpec=[]
                windowLiveTime=[]
                if icount == 0:
                    icount=icount+1
                    eCalEdges=[eCal[0]-0.5*(eCal[1]-eCal[0])]+[eCal[j]+0.5*(eCal[j+1]-eCal[j]) for j in range(len(eCal)-1)]+[eCal[-1]+0.5*(eCal[-1]-eCal[-2])]
                    for k in range(len(initSpectra)):
                        rebinedSpectrum=rebin_Archer(eCalEdges, initSpectra[k], binEdges)
                        rebinedSpectra.append(rebinedSpectrum)
            
            #energy calibration is bin centers, make an array of edges
            if icount == 0:
                initSpectra.append(spectrum)
                dum=[]

                ###
                # start algorithms
                ###
                #sprt.doSPRT(liveTime,energies,spectrum)
                #kSigma.dokSigma(liveTime,energies,spectrum)
            else:
                eCalEdges=[eCal[0]-0.5*(eCal[1]-eCal[0])]+[eCal[j]+0.5*(eCal[j+1]-eCal[j]) for j in range(len(eCal)-1)]+[eCal[-1]+0.5*(eCal[-1]-eCal[-2])]
                rebinedSpectrum=rebin_Archer(eCalEdges, spectra[i], binEdges)
                rebinedSpectra.append(rebinedSpectrum)
                ###
                # Run algorithms
                ###


        if not hdf5write:
            f=h5py.File(hdf5outPath+outFile,"w")
            hdf5write=True

        dset=f.create_dataset('3x3Spectra',(len(rebinedSpectra),len(rebinedSpectra[0])),dtype='f', compression="gzip")
        dset[...]=rebinedSpectra
        dset=f.create_dataset('3x3LiveTimes',(len(liveTimes),),dtype='f', compression="gzip")
        dset[...]=liveTimes
        dset=f.create_dataset('3x3Times',(len(times),),dtype='float64',compression="gzip")
        dset[...]=times
        fileName=dbFile.replace('.sqlite3','_peakLocs.npy')
        np.save(inPath+'logs/'+fileName,peakLocs)

    f.close()
    

