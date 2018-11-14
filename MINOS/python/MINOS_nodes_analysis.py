import numpy as np
import os
import h5py
import matplotlib

import matplotlib.pyplot as plt
import MINOS_analysis
from MINOS_analysis import node_timeseries,ksigma_RIO_alg,linFunc,node_colors,fixed_threshold_alg,sprt_alg
import matplotlib.dates as md
import datetime as dt
import time
from itertools import groupby
from operator import itemgetter
from scipy.optimize import curve_fit

directory_path = '/Volumes/Ian External HD/Node Data/'
csv_path = '/Volumes/Ian External HD/Node Data/'
nodes = [x for x in os.listdir(directory_path) if 'MUSE' in x]
outPath = directory_path+'figures/'

# Y-min/Y-max for CPS plots
counts_range = {'MUSE01':[1300,2100],'MUSE04':[1500,2900],\
    'MUSE06':[500,1800],'MUSE10':[1300,2400],\
    'MUSE11':[1500,2900],'MUSE12':[1300,2500]}

nodes = ['MUSE01','MUSE04','MUSE06','MUSE11','MUSE12']
Node_counts = {}
Node_times = {}
Node_spectra = {}
Node_datetimes = {}
Node_livetimes = {}
for i in nodes:
    print 'Reading %s'%i
    inPath = directory_path+str(i)+'/'#+'/hdf5Files/'
    dbFiles = [x for x in os.listdir(inPath) if '_ALG' not in x]
    #dbFiles = [x for x in dbFiles if (('09-26T' in x)or('09-27T' in x)or('09-28T' in x)or('09-29T' in x)or('09-30T' in x) or ('09-31T' in x)or('10-01T' in x))]
    dbFiles = [x for x in dbFiles if (('08-08T' in x)or('08-09T' in x)or('08-10T' in x)or('08-11T' in x)or('08-12T' in x) or ('08-13T' in x)or('08-14T' in x))]
    #dbFiles = [x for x in dbFiles if (('08-15T' in x)or('08-16T' in x)or('08-17T' in x)or('08-18T' in x)or('08-19T' in x) or ('08-20T' in x)or('08-21T' in x))]

    save_string = '08_08-08_14'

    first_time_times = True; first_time_spectra = True; first_time_livetimes = True

    for dbFile in dbFiles:

        hdf5 = h5py.File(inPath+dbFile,'r')
        print dbFile
        for j in hdf5.keys():

            if '2x4x16Spectra' in j:
                spectra = np.array(hdf5.get(str(j)))
                counts = np.sum(np.array(hdf5.get(str(j))),axis=1)
                if not first_time_spectra:
                    Node_counts[i] = np.concatenate([Node_counts[i],counts])
                    Node_spectra[i] = np.concatenate([Node_spectra[i],spectra])
                else:
                    Node_counts[i] = counts
                    Node_spectra[i] = spectra
                    first_time_spectra = False

            if '2x4x16Times' in j:
                if not first_time_times:
                    Node_times[i] = np.concatenate([Node_times[i],np.array(hdf5.get(str(j)))])
                else:
                    Node_times[i] = np.array(hdf5.get(str(j)))
                    first_time_times = False

            if '2x4x16LiveTimes' in j:
                if not first_time_livetimes:
                    Node_livetimes[i] = np.concatenate([Node_livetimes[i],np.array(hdf5.get(str(j)))])
                else:
                    Node_livetimes[i] = np.array(hdf5.get(str(j)))
                    first_time_livetimes = False

        Node_datetimes[i] = np.array([dt.datetime.fromtimestamp(x) for x in Node_times[i]])

    plotIt = False
    k40_peak = []
    k40_r2 = []
    th232_peak = []
    th232_r2 = []
    calc_time = []
    try:
        print lkjsadflkj
        for i in Node_spectra:
            step = 600
            step0=0

            for k in range(step,len(Node_spectra[i]),step):
                time = Node_times[i][k]
                spec = list(Node_spectra[i][step0:k])
                spectra = np.sum(spec,axis=0)

                # Fitting K40
                eMin = 420
                eMax = 540
                sigma = 3
                indexs=np.arange(eMin,eMax)
                xs=indexs
                ys=spectra[indexs]
                a=np.min(ys)
                b=np.max(ys)/(sigma*np.sqrt(2.0*np.pi))
                c=(ys[-1]-ys[0])/(eMax-eMin)
                mu=500
                popt,pcov=curve_fit(MINOS_analysis.gausswLine,xs,ys,p0=[a,b,c,mu,sigma])
                k40_peak.append(popt[3])
                calc_time.append(time)
                # Calculating R2 value for fit
                calc_ys = MINOS_analysis.gausswLine(xs,*popt)
                ss_res = np.sum((ys - calc_ys) ** 2)
                ss_tot = np.sum((ys - np.mean(ys)) ** 2)
                r2 = 1 - (ss_res / ss_tot)
                k40_r2.append(r2)

                # Fitting Th-232
                eMin = 800
                eMax = 960
                sigma = 4
                indexs=np.arange(eMin,eMax)
                xs_th=indexs
                ys_th=spectra[indexs]
                max = np.argmax(ys_th)
                eMin = max+810-30
                eMax = max+810+30
                a=np.min(ys)
                b=np.max(ys)/(sigma*np.sqrt(2.0*np.pi))
                c=(ys[-1]-ys[0])/(eMax-eMin)
                mu=870
                popt_th,pcov_th=curve_fit(MINOS_analysis.gausswLine,xs_th,ys_th,p0=[a,b,c,mu,sigma])
                th232_peak.append(popt_th[3])
                # Calculating R2 value for fit
                calc_ys = MINOS_analysis.gausswLine(xs_th,*popt)
                ss_res = np.sum((ys_th- calc_ys) ** 2)
                ss_tot = np.sum((ys_th - np.mean(ys_th)) ** 2)
                r2 = 1 - (ss_res / ss_tot)
                th232_r2.append(r2)

                if False:
                    fig,ax = plt.subplots()
                    ax.plot(spectra)
                    ax.set_yscale('log')
                    ax.plot(xs,MINOS_analysis.gausswLine(xs,*popt),lw=2,marker='',ls='--',color='k',label=str(popt[3]))
                    ax.plot(xs_th,MINOS_analysis.gausswLine(xs_th,*popt_th),lw=2,marker='',ls='--',color='k',label=str(popt_th[3]))
                    plt.grid(True,which='both',alpha=0.5)
                    title_string = i+'-'+str(time)
                    plt.title(title_string)
                    plt.savefig(outPath+i+'/'+i+'-'+str(time)+'.png',format='png',dpi=200)
                    plt.close('all')

                step0=k
            fig, axes = plt.subplots(nrows=2)
            axes[0].plot(k40_peak, 'o', markersize=0.75, label='K40 Peak Channel');
            axes[0].grid(alpha=0.5);
            axes[0].set_ylabel('Channel')
            axes[1].plot(th232_peak, 'ro', markersize=0.75, label='Th232 Peak Channel');
            axes[1].grid(alpha=0.5);
            axes[1].set_ylabel('Channel')
            axes[0].legend(fancybox=True);
            axes[1].legend(fancybox=True)
            axes[0].set_title('Peak locations every 5-minutes')
            plt.savefig(outPath+i+'/' + i + '_peak_tracking_'+save_string+'.png', format='png', dpi=200)
    except:
        pass
    runAlarms = False

    if runAlarms:
        back_windows = 20
        fore_windows = 8
        ksigma_threshold = 10
        sprt_threshold = 1100
        roi_counts,ksigma_alarms,ratio = ksigma_RIO_alg(Node_spectra[i],fore_windows,back_windows,ksigma_threshold,'spectra')
        sprt_values,sprt_alarms = sprt_alg(Node_counts[i],Node_livetimes[i],fore_windows,back_windows,sprt_threshold)

        sprt_indexs = []
        for k, g in groupby(enumerate(sprt_alarms), lambda (i, x): i-x):
            sprt_indexs.append(map(itemgetter(1), g))

        if ( (len(ksigma_alarms)>0)or(len(sprt_alarms)>0) ):
            fig,ax = plt.subplots(figsize=(8,8), nrows=3,sharex=True)
            ax[0].plot(Node_datetimes[i],roi_counts,c=node_colors(i),linewidth=0.75,label=str(i))
            ax[0].set_ylim(counts_range[i][0],counts_range[i][1])
            ax[0].set_title(str(i))
            ax[0].grid(alpha=0.5,linewidth=0.5)
            ax[0].set_xlim(Node_datetimes[i][0],Node_datetimes[i][-1])
            ax[0].set_ylabel('CPS')

            ax[1].plot(Node_datetimes[i][(len(Node_datetimes[i])-len(ratio))::],ratio,label='k-sigma',zorder=10)
            ax[1].set_ylim(0,ksigma_threshold+5)
            datetime_alarms = [Node_datetimes[i][x] for x in ksigma_alarms]
            [ax[1].axvline(x,c='r',alpha=0.35,label='Alarm',zorder=1) for x in datetime_alarms]
            ax[1].axhline(ksigma_threshold,linewidth=2,color='k',linestyle='--')
            ax[1].set_xlim(Node_datetimes[i][0],Node_datetimes[i][-1])
            ax[1].grid(alpha=0.5,linewidth=0.5)
            ax[1].set_ylabel('k-sigma')
            ax[1] = plt.gca()
            ax[1].xaxis.set_major_formatter(md.DateFormatter('%m/%d %H:%M'))
            ax[1].set_xlabel('Time (month/day hour:minute)')
            
            ax[2].plot(Node_datetimes[i][(len(Node_datetimes[i])-len(sprt_values))::],np.array(sprt_values),label='SPRT',zorder=10)
            if len(sprt_indexs)>0:
                for index in sprt_indexs:
                    ax[2].axvspan(Node_datetimes[i][index[0]],Node_datetimes[i][index[-1]],color='red',alpha=0.3)
            ax[2].axhline(sprt_threshold,linewidth=2,color='k',linestyle='--')
            ax[2].set_xlim(Node_datetimes[i][0],Node_datetimes[i][-1])
            ax[2].grid(alpha=0.5,linewidth=0.5)
            ax[2].set_ylabel('SPRT')
            ax[2] = plt.gca()
            ax[2].xaxis.set_major_formatter(md.DateFormatter('%m/%d %H:%M'))
            ax[2].set_xlabel('Time (month/day hour:minute)')
            for ax in fig.axes:
                matplotlib.pyplot.sca(ax)
                plt.xticks(rotation=25)
            plt.subplots_adjust(hspace=0.05)
        plt.savefig(outPath+i+'/' + i +'_Alarms_'+save_string+'.png',format='png',dpi=200)
fig,ax=plt.subplots(6,sharex=True)
count=0
for i in Node_counts:
    print i
    ax[count].plot(Node_datetimes[i],Node_counts[i],label=str(i),c=node_colors(i))
    ax[count].set_xlim(Node_datetimes[i][0], Node_datetimes[i][-1])
    ax[count].set_ylim(counts_range[i][0], counts_range[i][1])
    ax[count].grid(alpha=0.5, linewidth=0.5)
    ax[count].legend(loc='upper right',fancybox=False,shadow=True)
    count+=1
ax[5] = plt.gca()
ax[5].xaxis.set_major_formatter(md.DateFormatter('%m/%d'))
ax[5].set_xlabel('Date (month/day)')
plt.show()
#node_timeseries(Node_times,Node_counts,nodes)
if False:
    plt.savefig(outPath+'test.pdf',format='pdf',dpi=200)
if False:
    plt.show()

