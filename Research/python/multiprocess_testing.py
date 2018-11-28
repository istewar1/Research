import numpy as np
import time
import sqlite3
from multiprocessing import Pool

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

def string_to_int(data):
    return [np.fromstring(x[0],sep=',',dtype='int') for x in data]

node = 'MUSE01'
inPath='/Volumes/Ian External HD/Node Data/sqlitefiles/'+node+'/'
dbFile = 'MUSE01-2018-09-28T23.49.40.590.sqlite3'

dataTable = 'Det_3x3_data'
detDB = database(inPath + dbFile)
detDB.getColumn(dataTable, 'Spectrum__IntArray')

start = time.time()
spectra=[np.fromstring(x[0],sep=',',dtype='int') for x in detDB.data]
end = time.time()
print np.shape(spectra)
print end-start

start = time.time()
strarray = np.array([np.array(x) for x in detDB.data])
print
end = time.time()
print strarray[0]
print np.shape(s)
print end-start




