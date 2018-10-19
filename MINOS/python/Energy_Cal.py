import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as md
import os
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
from scipy.interpolate import UnivariateSpline
import sqlite3

def generate_calibration_axis(bins1, energies1):
    """
    Perform a Akima spline interpolation between the given bins and energies, and a linear polynomial
    extrapolation beyond the known peaks, which is fit between the last two known peaks
    :param bins: Bins at which each energy is located
    :param energies: Energies being used for calibration
    :param display: set this to (characterizerinstance).NONE, (future)
    :return: the new x_axis energy data from the spline fit
    """
    bins = np.array(bins1)
    energies = np.array(energies1)
    bins = np.append(bins,0)
    energies = np.append(energies,0.0)
    # sort the bins and energies
    bins, energies = zip(*sorted(zip(bins, energies)))
    print("Using Bins: " + str(bins))
    print("Using Energies: " + str(energies))
    akima_fit = interpolate.Akima1DInterpolator(x=bins, y=energies)
    ext = np.polyfit([bins[-2], bins[-1]], [energies[-2], energies[-1]], 1)
    interpolated_energy_data = akima_fit(np.arange(0, bins[-1]))
    extrapolated_energy_data = np.polyval(ext, np.arange(bins[-1], 2048))
    energy_data = np.concatenate([interpolated_energy_data, extrapolated_energy_data])
    print(len(energy_data))
    return energy_data

global database

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
        c.execute("select rowid,* from %s order by Time asc"  % (tblname))
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
        return data        

    def closeConn(self):
        self.conn.close()
        self.data=[]

#Gaussian + line
def peakFunc(x,a,b,c,mu,sigma):
    '''
    Peak fit function with Gaussian on sloped line to account
    for down scatter contribution from higher E photons.
    '''
    return (        
        #line
        a+c*x+
        #photopeak
        b*np.exp(-(x-mu)**2/(2.0*sigma**2))
        )
        
def getIsotopeEnergies_large(isotope):
    '''
    Description: Provides energy and ROI for isotope.
    Param: isotope = string (e.g. 'Cs-137')
    '''
    peaks = []; peak_windows_lower = []; peak_windows_upper = []
    if isotope == "Cs-137":
        peaks.append(662.0)
        peak_windows_lower.append(570)
        peak_windows_upper.append(780)
    elif isotope == "Co-60":
        peaks.append(1173.0)
        peak_windows_lower.append(1050)
        peak_windows_upper.append(1240)
        peaks.append(1332.0)
        peak_windows_lower.append(1250)
        peak_windows_upper.append(1420)
    elif isotope == "Ba-133":
        peaks.append(80.0)
        peak_windows_lower.append(60)
        peak_windows_upper.append(100)
        peaks.append(356.0)
        peak_windows_lower.append(320)
        peak_windows_upper.append(450)
    elif isotope == "Co-57":
        peaks.append(122.0)
        peak_windows_lower.append(110)
        peak_windows_upper.append(135)
    elif isotope == "K-40":
        peaks.append(1461.0)
        peak_windows_lower.append(1280)
        peak_windows_upper.append(1600)
    elif isotope == "Background":
        peaks.append(1461.0)
        peak_windows_lower.append(1280)
        peak_windows_upper.append(1600)
    elif isotope == "Th-232":
        peaks.append(2614.0)
        peak_windows_lower.append(2430)
        peak_windows_upper.append(2800)
    elif isotope == "Eu-152":
        peaks.append(40.118)
        peak_windows_lower.append(25)
        peak_windows_upper.append(55)
        peaks.append(121.78)
        peak_windows_lower.append(100)
        peak_windows_upper.append(140)
        peaks.append(344.29)
        peak_windows_lower.append(300)
        peak_windows_upper.append(390)
        peaks.append(778.904)
        peak_windows_lower.append(710)
        peak_windows_upper.append(850)
        peaks.append(1408.006)
        peak_windows_lower.append(1300)
        peak_windows_upper.append(1550)
    elif isotope == "Am-241":
        peaks.append(59.5)
        peak_windows_lower.append(42)
        peak_windows_upper.append(75)
    elif isotope == "Ho-166":
        peaks.append(184.41)
        peak_windows_lower.append(160)
        peak_windows_upper.append(220)
        peaks.append(280.46)
        peak_windows_lower.append(250)
        peak_windows_upper.append(320)
        peaks.append(410.95)
        peak_windows_lower.append(375)
        peak_windows_upper.append(450)
        peaks.append(711.68)
        peak_windows_lower.append(680)
        peak_windows_upper.append(730)
        peaks.append(810.29)
        peak_windows_lower.append(785)
        peak_windows_upper.append(850)
    return peaks,peak_windows_lower,peak_windows_upper
   
inPath = '/Users/i6o/Downloads/'
energy = np.genfromtxt(inPath+'Example_Energy_Pairs.json',delimiter=',')
energy = np.array(energy)
dbFiles = [x for x in os.listdir(inPath) if '.sqlite3' in x]
dbFiles = ['MUSE11-Energy_cal-2012-01-08T15.00.43.988.sqlite3']

for dbFile in dbFiles:
    print dbFile
    detDB=database(inPath+dbFile)

    dataTable='Det_2x4x16_215_data'
    detDB.getColumn(dataTable,'Spectrum__IntArray')
    spectra_2x4x16=np.array([np.fromstring(x[0],sep=',',dtype='int') for x in detDB.data])
    detDB.getColumn(dataTable,'Time')
    Time_2x4x16=np.array(detDB.data)[:,0]
    dateTime_2x4x16 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_2x4x16])
    detDB.getColumn(dataTable,'Live_Time')
    Live_Time_2x4x16=np.array(detDB.data)[:,0]
    
    cpms_2x4x16 = np.array([np.sum(x) for x in spectra_2x4x16])
    cps_2x4x16 = np.add.reduceat(cpms_2x4x16, np.arange(0, len(cpms_2x4x16), 10))
    dates_cps_2x4x16 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_2x4x16[::10]])
    spectra_seconds_2x4x16 = np.add.reduceat(spectra_2x4x16, np.arange(0, len(spectra_2x4x16), 10))

    dataTable='Det_3x3_12062403_data'
    detDB.getColumn(dataTable,'Spectrum__IntArray')
    spectra_3x3=np.array([np.fromstring(x[0],sep=',',dtype='int') for x in detDB.data])
    detDB.getColumn(dataTable,'Time')
    Time_3x3=np.array(detDB.data)[:,0]
    dateTime_3x3 = np.array([dt.datetime.fromtimestamp(ts) for ts in Time_3x3])
    detDB.getColumn(dataTable,'Live_Time')
    Live_Time_3x3=np.array(detDB.data)[:,0]
    
    cpms_3x3 = np.array([np.sum(x) for x in spectra_3x3])
    cps_3x3 = np.add.reduceat(cpms_3x3, np.arange(0, len(cpms_3x3), 10))
    dates_cps_3x3=np.array([dt.datetime.fromtimestamp(ts) for ts in Time_3x3[::10]])
    spectra_seconds_3x3 = np.add.reduceat(spectra_3x3, np.arange(0, len(spectra_3x3), 10))
    '''
    Energy Calibration Section
    '''
    print Time_3x3[0]
    print Time_3x3[12048]
    print Time_3x3[12391]
    print Time_3x3[18984]
    print Time_3x3[19945]
    print Time_3x3[-1]

    # Imported from Lab Notes
    time0 = [1326052843,1326054082,1326054838] #Estimated Posix source introduction
    time1 = [1326054048,1326054742,1326055940] #Estimated Posix source departure
    source_names = ['K-40','Th-232','Eu-152']
    
    calibration_energies_pairs_2x4x16 = {} # {energy_keV:channel,...}
    calibration_energies_pairs_3x3 = {}
        
    for i in range(len(time0)):
        
        start = time0[i]
        end = time1[i]
        
        index_2x4x16 = np.where( (Time_2x4x16>start) & (Time_2x4x16<end) )[0]
        counts = spectra_2x4x16[index_2x4x16]
        final_spectra_2x4x16 = np.sum(counts,axis=0)
        
        index_3x3 = np.where( (Time_3x3>start) & (Time_3x3<end) )[0]
        counts = spectra_3x3[index_3x3]
        final_spectra_3x3 = np.sum(counts,axis=0)
    
        peaks,peak_windows_lower,peak_windows_upper = getIsotopeEnergies_large(source_names[i])
        print '\t Fitting isotope:\t%s'%(source_names[i])
        for z in range(len(peaks)):

            j,k,l = (peaks[z],peak_windows_lower[z],peak_windows_upper[z])
            indexs = np.where( (energy > k) & (energy < l) )
                        
            # Detector 2x4x16 -- 217
            spec = final_spectra_2x4x16[indexs[0][0]:indexs[0][len(indexs[0])-1]]
            e = energy[indexs[0][0]:indexs[0][len(indexs[0])-1]]
            a=spec[-1]
            b=np.max(spec)-spec[0]
            c=-1.0*np.abs(spec[0]-spec[-1])/(energy[-1]-energy[0]) 
            mu=e[np.argmax(spec)]
            sigma=0.9
            p0=[a,b,c,mu,sigma]
            popt,pcov=curve_fit(peakFunc,e,spec,p0=p0)
            popt[4]=np.abs(popt[4])
            channel = (np.abs(energy-e[np.argmax(peakFunc(e,*popt))])).argmin()
            if j not in calibration_energies_pairs_2x4x16.keys():
                calibration_energies_pairs_2x4x16[j] = [channel]
            else:
                calibration_energies_pairs_2x4x16[j].append(channel)
            
            # Detector 3x3 -- 3x3
            spec = final_spectra_3x3[indexs[0][0]/2:indexs[0][len(indexs[0])-1]/2]
            channel_spec = np.arange(indexs[0][0]/2,indexs[0][len(indexs[0])-1]/2)
            a=spec[-1]
            b=np.max(spec)-spec[0]
            c=-1.0*np.abs(spec[0]-spec[-1])/(channel_spec[-1]-channel_spec[0]) 
            mu=channel_spec[np.argmax(spec)]
            sigma=1.5
            p0=[a,b,c,mu,sigma]
            popt,pcov=curve_fit(peakFunc,channel_spec,spec,p0=p0)
            popt[4]=np.abs(popt[4])
            channel = (np.abs(np.arange(0,len(final_spectra_3x3))-channel_spec[np.argmax(peakFunc(channel_spec,*popt))])).argmin()
            if j not in calibration_energies_pairs_3x3.keys():
                calibration_energies_pairs_3x3[j] = [channel]
            else:
                calibration_energies_pairs_3x3[j].append(channel)
          
            plotIt = True
            
            if plotIt:
                '''
                Plots fitted data on top of spectra
                '''
                string = str(j)+' gauss fit'
                plt.figure(i)
                plt.plot(e,peakFunc(e,*popt),label=string,zorder=10)
                
        if plotIt:
            plt.figure(i);plt.plot(energy,final_spectra_2x4x16,label=source_names[i],zorder=0)
            plt.yscale('log')
            plt.legend()
            plt.grid(True,which='both')
            
        plotDiagnostics = False   
        
        if plotDiagnostics:
            '''
            Used for Diagnositcs
            Plots two spectra: 
                (1) using guess-energy calibration pairs
                (2) using channel number
            '''
            fig,(ax1,ax2) = plt.subplots(nrows=2)
            ax1.plot(energy,final_spectra_2x4x16,label=source_names[i])
            ax1.set_yscale('log')
            ax1.legend(fancybox=True,shadow=True)
            ax1.grid(True,which='both')
            ax1.set_xlabel('Energy Guess (keV)')
            ax1.set_ylabel('Intensity')
            ax2.plot(final_spectra_2x4x16,label=source_names[j])
            ax2.set_yscale('log')
            ax2.legend(fancybox=True,shadow=True)
            ax2.grid(True,which='both')
            ax2.set_xlabel('Channel')
            ax2.set_ylabel('Intensity')
            fig.tight_layout()
                
    '''
    Creating energy pairs
    '''   
    energy_values_2x4x16 = [] ; energy_channel_2x4x16 = []
    for i in calibration_energies_pairs_2x4x16:
        if len(calibration_energies_pairs_2x4x16[i])==1:
            energy_values_2x4x16.append(i)
            energy_channel_2x4x16.append(int(calibration_energies_pairs_2x4x16[i][0]))
        else:
            energy_values_2x4x16.append(i)
            energy_channel_2x4x16.append(int(np.mean(calibration_energies_pairs_2x4x16[i])))
    energy_pairs_2x4x16=generate_calibration_axis(energy_channel_2x4x16,energy_values_2x4x16)
    
    energy_values_3x3 = [] ; energy_channel_3x3 = []
    for i in calibration_energies_pairs_3x3:
        if len(calibration_energies_pairs_3x3[i])==1:
            energy_values_3x3.append(i)
            energy_channel_3x3.append(int(calibration_energies_pairs_3x3[i][0]))
        else:
            energy_values_3x3.append(i)
            energy_channel_3x3.append(int(np.mean(calibration_energies_pairs_3x3[i])))
    energy_pairs_3x3=generate_calibration_axis(energy_channel_3x3,energy_values_3x3)
    
    '''
    us = UnivariateSpline(channels[0:-1],energies[0:-1],k=2,s=0.0,ext=0)
    interp = us(np.arange(0,channels[-2]))
    ext = np.polyfit(np.array([channels[-2],channels[-1]]),np.array([energies[-2],energies[-1]]),1)
    extrap = np.polyval(ext, np.arange(channels[-2], 2048))
    cs = np.concatenate([interp,extrap])
    x_range = np.arange(len(cs))
    fig,axes = plt.subplots(2,1)
    axes[0].plot(x_range, np.array(cs), label='Cubic Spline')
    axes[0].plot(channels,np.array(energies),'r*',label='Fitted Energies')
    axes[0].grid(alpha=0.5)
    axes[0].set_xlabel('Channel')
    axes[0].set_ylabel('Energy')
    axes[1].plot(x_range, np.array((x_range))/np.array(cs), label='Cubic Spline')
    axes[1].plot(channels,np.array(channels)/np.array(energies),'r*',label='Fitted Energies')
    axes[1].grid(alpha=0.5)
    axes[1].set_xlabel('Channel')
    axes[1].set_ylabel('Channel/Energy')
    fig.tight_layout()  
    '''
    
if (plotIt or plotDiagnostics):
    plt.show()