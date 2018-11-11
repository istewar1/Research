import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import time
from scipy.optimize import curve_fit

def writeSPE(path,label, liveTime, realTime, time, spectrum, eCal=[0.0, 0.0, 0.0]):
    contents = '$SPEC_ID:\r\n%s\r\n' % label
    contents = contents + '$MEAS_TIM:\r\n%.4f %.4f\r\n' % (liveTime, realTime)
    contents = contents + '$DATE_MEA:\r\n%s\r\n' % (time.strftime('%d-%b-%Y %I:%M:%S %p'))
    contents = contents + '$DATA:\r\n'
    for i in range(len(spectrum)):
        if i == 0:
            contents = contents + '%d %d\r\n' % (spectrum[i], len(spectrum))
        else:
            contents = contents + '%d\r\n' % (spectrum[i])
    contents = contents + '$ENER_FIT:\r\n%.4f %.4f %.4f\r\n' % (eCal[0], eCal[1], eCal[2])
    contents = contents + '$MCA_CAL:\r\n%.4f %.4f %.4f\r\n' % (eCal[0], eCal[1], eCal[2])
    contents = contents + '$ENDRECORD:'
    f = open(path+label + '.SPE', 'w')
    f.write(contents)
    f.close()
    return

def quadratic(x,a,b,c):
    return a + b * x + c * x * x

class database(object): #database class

    def __init__(self,DBfilename):
        self.conn = sqlite3.connect(DBfilename)
        self.conn.row_factory = sqlite3.Row

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

def quadratic(x,a,b,c):
    return a + b * x + c * x * x

xfmt = mdates.DateFormatter('%H:%M:%S')
dates = np.vectorize(dt.datetime.fromtimestamp)

inPath = '/Volumes/IAN USB/WIND/GADRAS dat Creation/'
dbFile = inPath+'archerair2012-WIND_Ba133-2018-10-30T16.30.27.689.sqlite3'
detDB = database(dbFile)
tables =['Det1_BR_4in_data','Det3_BL_4in_data','Det5_FC_2in_data']
energy_pairs = ['eCal_Det1.npy','eCal_Det3.npy','eCal_Det5.npy']

for i in range(len(tables)):
    dataTable = tables[i]
    energy = np.load(inPath+'/energy_pairs/'+energy_pairs[i])
    detDB.getColumn(dataTable, 'Spectrum__IntArray')
    spectra = np.array([np.fromstring(x[0], sep=',', dtype='int') for x in detDB.data])
    detDB.getColumn(dataTable, 'Time')
    times = np.array(detDB.data)
    times = np.array([float(x[0]) for x in times])
    detDB.getColumn(dataTable, 'Live_Time')
    liveTimes = np.array(detDB.data)
    liveTimes = np.array([float(x[0]) for x in liveTimes])

    # Creating energy parameters from energy cal previously created
    energy = np.load(inPath+'energy_pairs/'+energy_pairs[i])
    popt, popc = curve_fit(quadratic, range(0, len(energy)), energy)
    values = quadratic(np.arange(0,len(energy)), *popt)

    liveTime = np.sum(liveTimes)
    spectrum = np.sum(spectra, axis=0)

    fig,ax = plt.subplots()
    ax.plot(values,spectrum,'r',label='quadratic fit');ax.grid(alpha=0.5)
    ax.plot(energy, spectrum, 'b', label='quadratic fit');ax.set_yscale('log');ax.legend()
    ax.set_title(dataTable)
    ax.set_xlabel('Energy (keV)')
    plt.show(False)

    # save to SPE
    label = dataTable+'_'+'Ba133'
    realTime = times[-1] - times[0]
    timeStamp = dt.datetime.fromtimestamp(times[0])
    writeSPE(inPath,label, liveTime, realTime, timeStamp, spectrum, eCal=[float(x) for x in popt])

create_energy_parameters = False
if create_energy_parameters:
    inPath = '/Volumes/IAN USB/WIND/GADRAS dat Creation/energy_pairs/'
    e_file = 'eCal_Det1.npy'

    energy_pairs = np.load(inPath + e_file)
    popt, popc = curve_fit(quadratic, energy_pairs, range(1, 2049))
    values = quadratic(np.array(range(0, 2048)), *popt)
    print popt
    plt.figure();
    plt.plot(values, 'r--')
    plt.plot(energy_pairs, 'b')
    plt.grid()
    plt.show()
