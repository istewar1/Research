import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import time

def writeSPE(label,liveTime,realTime,time,spectrum,eCal=[0.0,0.0,0.0]):
    contents='$SPEC_ID:\r\n%s\r\n' % label
    contents=contents+'$MEAS_TIM:\r\n%.4f %.4f\r\n' % (liveTime,realTime)
    contents=contents+'$DATE_MEA:\r\n%s\r\n' % (time.strftime('%d-%b-%Y %I:%M:%S %p'))
    contents=contents+'$DATA:\r\n'
    for i in range(len(spectrum)):
        if i == 0:
            contents=contents+'%d %d\r\n' % (spectrum[i],len(spectrum))
        else:
            contents=contents+'%d\r\n' % (spectrum[i])
    contents=contents+'$ENER_FIT:\r\n%.4f %.4f %.4f\r\n' % (eCal[0],eCal[1],eCal[2])
    contents=contents+'$MCA_CAL:\r\n%.4f %.4f %.4f\r\n' % (eCal[0],eCal[1],eCal[2])
    contents=contents+'$ENDRECORD'
    f=open(label+'.SPE','w')
    f.write(contents)
    f.close()
    return
    

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

# some definitions
xfmt = mdates.DateFormatter('%H:%M:%S')
dates=np.vectorize(dt.datetime.fromtimestamp)

####### lets grab all the data from the database
dbFile='MUSE04-2018-10-25T14.19.34.827.sqlite3'
detDB=database(dbFile)
dataTable='Det_2x4x16_data'
detDB.getHeader(dataTable)
header=np.array(detDB.header)
detDB.getColumn(dataTable,'Live_Time')
liveTimes=np.array(detDB.data)[:,0]
length=10*(len(liveTimes)/10)
detDB.getColumn(dataTable,'Time')
times=np.array(detDB.data)[:,0]
detDB.getColumn(dataTable,'Spectrum__IntArray')
spectra=np.array([np.fromstring(x[0],sep=',',dtype='int') for x in detDB.data])

#look at the total count rate to select time windows
fig=plt.figure()
ttimes=dates(times)
plt.step(ttimes,np.sum(spectra,axis=1)/liveTimes,lw=2)
plt.grid(True)
formatter = mdates.DateFormatter('%m/%d %H:%M')
plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
fig.autofmt_xdate()
plt.show(False)

####### lets select the background data
start='10/25/2018 14:20'
end='10/25/2018 14:50'

start=time.mktime(dt.datetime.strptime(start,'%m/%d/%Y %H:%M').timetuple())
end=time.mktime(dt.datetime.strptime(end,'%m/%d/%Y %H:%M').timetuple())

indexs=np.where( (times >=start) & (times <= end) )[0]

liveTime=np.sum(liveTimes[indexs])
spectrum=np.sum(spectra[indexs],axis=0)

#make a plot
plt.figure()
plt.plot(range(len(spectrum)),spectrum/liveTime,lw=2)
plt.grid(True)
plt.yscale('log')
plt.show(False)

#save to SPE
label='Background'
realTime=times[indexs[-1]]-times[indexs[0]]
timeStamp=dt.datetime.fromtimestamp(start)
writeSPE(label,liveTime,realTime,timeStamp,spectrum,eCal=[0.0,0.0,0.0])

####### lets select the Ba-133
start='10/25/2018 15:00'
end='10/25/2018 15:05'

start=time.mktime(dt.datetime.strptime(start,'%m/%d/%Y %H:%M').timetuple())
end=time.mktime(dt.datetime.strptime(end,'%m/%d/%Y %H:%M').timetuple())

indexs=np.where( (times >=start) & (times <= end) )[0]

liveTime=np.sum(liveTimes[indexs])
spectrum=np.sum(spectra[indexs],axis=0)

plt.figure()
plt.plot(range(len(spectrum)),spectrum/liveTime,lw=2)
plt.grid(True)
plt.yscale('log')
plt.show(False)

label='Ba-133'
realTime=times[indexs[-1]]-times[indexs[0]]
timeStamp=dt.datetime.fromtimestamp(start)
writeSPE(label,liveTime,realTime,timeStamp,spectrum,eCal=[0.0,0.0,0.0])

#Cs-137
start='10/25/2018 15:07'
end='10/25/2018 15:12'

start=time.mktime(dt.datetime.strptime(start,'%m/%d/%Y %H:%M').timetuple())
end=time.mktime(dt.datetime.strptime(end,'%m/%d/%Y %H:%M').timetuple())

indexs=np.where( (times >=start) & (times <= end) )[0]

liveTime=np.sum(liveTimes[indexs])
spectrum=np.sum(spectra[indexs],axis=0)

plt.figure()
plt.plot(range(len(spectrum)),spectrum/liveTime,lw=2)
plt.grid(True)
plt.yscale('log')
plt.show(False)

label='Cs-137'
realTime=times[indexs[-1]]-times[indexs[0]]
timeStamp=dt.datetime.fromtimestamp(start)
writeSPE(label,liveTime,realTime,timeStamp,spectrum,eCal=[0.0,0.0,0.0])

####### lets select the Eu-152
start='10/25/2018 15:15'
end='10/25/2018 15:25'

start=time.mktime(dt.datetime.strptime(start,'%m/%d/%Y %H:%M').timetuple())
end=time.mktime(dt.datetime.strptime(end,'%m/%d/%Y %H:%M').timetuple())

indexs=np.where( (times >=start) & (times <= end) )[0]

liveTime=np.sum(liveTimes[indexs])
spectrum=np.sum(spectra[indexs],axis=0)

plt.figure()
plt.plot(range(len(spectrum)),spectrum/liveTime,lw=2)
plt.grid(True)
plt.yscale('log')
plt.show(False)

label='Eu-152'
realTime=times[indexs[-1]]-times[indexs[0]]
timeStamp=dt.datetime.fromtimestamp(start)
writeSPE(label,liveTime,realTime,timeStamp,spectrum,eCal=[0.0,0.0,0.0])

####### lets select the Co-60
start='10/25/2018 15:35'
end='10/25/2018 15:40'

start=time.mktime(dt.datetime.strptime(start,'%m/%d/%Y %H:%M').timetuple())
end=time.mktime(dt.datetime.strptime(end,'%m/%d/%Y %H:%M').timetuple())

indexs=np.where( (times >=start) & (times <= end) )[0]

liveTime=np.sum(liveTimes[indexs])
spectrum=np.sum(spectra[indexs],axis=0)

plt.figure()
plt.plot(range(len(spectrum)),spectrum/liveTime,lw=2)
plt.grid(True)
plt.yscale('log')
plt.show(False)

label='Co-60'
realTime=times[indexs[-1]]-times[indexs[0]]
timeStamp=dt.datetime.fromtimestamp(start)
writeSPE(label,liveTime,realTime,timeStamp,spectrum,eCal=[0.0,0.0,0.0])

