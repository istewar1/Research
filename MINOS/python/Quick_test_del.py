import numpy as np
import os
import h5py
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import MINOS_analysis
from MINOS_analysis import node_timeseries, ksigma_RIO_alg, linFunc, node_colors, fixed_threshold_alg, sprt_alg
import matplotlib.dates as md
import datetime as dt
import time
from itertools import groupby
from operator import itemgetter
from scipy.optimize import curve_fit

directory_path = '/Volumes/Ian External HD/Node Data/Alarm_Analysis/'
csv_path = directory_path
nodes = [x for x in os.listdir(directory_path) if 'MUSE' in x]
outPath = directory_path + 'Figures/'

df =  pd.read_excel('/Users/i6o/MINOS/Node_Correlation_Events.xlsx')
alarm_events = np.array(df['Posix'])[24:]
alarm_events = alarm_events[np.logical_not(np.isnan(alarm_events))]

# Y-min/Y-max for CPS plots
counts_range = {'MUSE01': [1100, 2100], 'MUSE04': [1200, 2900], \
                'MUSE06': [500, 1800], 'MUSE10': [1100, 2400], \
                'MUSE11': [1500, 2900], 'MUSE12': [1100, 2500]}

nodes = ['MUSE01', 'MUSE04', 'MUSE06', 'MUSE10', 'MUSE11', 'MUSE12']

Node_counts = {}
Node_times = {}
Node_spectra = {}
Node_datetimes = {}
Node_livetimes = {}

for event in alarm_events:
    date = dt.datetime.fromtimestamp(event)
    print date
    month, day, hour = date.month, date.day, date.hour
    start, end = event-100,event+100
    Node_counts = {}
    Node_times = {}
    Node_spectra = {}
    Node_datetimes = {}
    Node_livetimes = {}
    for i in nodes:
        print 'Reading %s' % i
        inPath = directory_path + str(i) + '/' + '/hdf5Files/'
        dbFiles = [x for x in os.listdir(inPath) if '_ALG' not in x]
        # dbFiles = [x for x in dbFiles if (('09-26T' in x)or('09-27T' in x)or('09-28T' in x)or('09-29T' in x)or('09-30T' in x) or ('09-31T' in x)or('10-01T' in x))]
        # dbFiles = [x for x in dbFiles if (('2018-10-01T' in x)or('2018-10-02T' in x)or('2018-10-03T' in x)or('2018-10-04T' in x)or('2018-10-05T' in x)or('2018-10-06T' in x)or('2018-10-07T' in x))]
        if day < 10:
            if month < 10:
                string = '2018-0%i-0%iT'%(month,day)
            elif month >= 10:
                string = '2018-%i-0%iT' % (month, day)
        else:
            if month<10:
                string = '2018-0%i-%iT' % (month, day)
            elif month>=10:
                string = '2018-%i-%iT' % (month, day)


        dbFiles = [x for x in dbFiles if (string in x)]

        save_string = string

        first_time_times = True;
        first_time_spectra = True;
        first_time_livetimes = True

        for dbFile in dbFiles:

            hdf5 = h5py.File(inPath + dbFile, 'r')
            print dbFile
            for j in hdf5.keys():

                if '2x4x16Spectra' in j:
                    spectra = np.array(hdf5.get(str(j)))
                    counts = np.sum(np.array(hdf5.get(str(j))), axis=1)
                    if not first_time_spectra:
                        Node_counts[i] = np.concatenate([Node_counts[i], counts])
                        Node_spectra[i] = np.concatenate([Node_spectra[i], spectra])
                    else:
                        Node_counts[i] = counts
                        Node_spectra[i] = spectra
                        first_time_spectra = False

                if '2x4x16Times' in j:
                    if not first_time_times:
                        Node_times[i] = np.concatenate([Node_times[i], np.array(hdf5.get(str(j)))])
                    else:
                        Node_times[i] = np.array(hdf5.get(str(j)))
                        first_time_times = False

                if '2x4x16LiveTimes' in j:
                    if not first_time_livetimes:
                        Node_livetimes[i] = np.concatenate([Node_livetimes[i], np.array(hdf5.get(str(j)))])
                    else:
                        Node_livetimes[i] = np.array(hdf5.get(str(j)))
                        first_time_livetimes = False

            Node_datetimes[i] = np.array([dt.datetime.fromtimestamp(x) for x in Node_times[i]])

    indexes = np.where(np.logical_and(Node_times[i]>=start, Node_times[i]<=end))

    date_start = dt.datetime.fromtimestamp(start)
    month_start,day_start,hour_start,minute_start,second_start = date_start.month,date_start.day,date_start.hour,date_start.minute,date_start.second
    date_end = dt.datetime.fromtimestamp(end)
    month_end,day_end,hour_end,minute_end,second_end = date_end.month,date_end.day,date_end.hour,date_end.minute,date_end.second

    save_string = 'Correlation_Event_DATE-%i-%i-%i_TIME-%i-%i-%i_%i-%i-%i'%(date.year,date.month,date.day,hour_start,minute_start,second_start,hour_end,minute_end,second_end)

    fig, ax = plt.subplots(6, sharex=True, figsize=(10, 10))
    count = 0
    try:
        for i in nodes:
            print i
            indexes = np.where(np.logical_and(Node_times[i] >= start, Node_times[i] <= end))[0]
            ax[count].plot(Node_datetimes[i][indexes], Node_counts[i][indexes], label=str(i), c=node_colors(i))
            ax[count].set_xlim(Node_datetimes[i][indexes][0], Node_datetimes[i][indexes][-1])
            #ax[count].set_ylim(counts_range[i][0], counts_range[i][1])
            ax[count].grid(alpha=0.5, linewidth=0.5)
            ax[count].legend(loc='upper right', fancybox=False, shadow=True)
            count += 1
        ax[5] = plt.gca()
        ax[5].xaxis.set_major_formatter(md.DateFormatter('%m/%d %H:%M:%S'))
        ax[5].set_xlabel('Date (Hour:Minute:Second)')
        ax[0].set_title(save_string)
        for ax in fig.axes:
            matplotlib.pyplot.sca(ax)
            plt.xticks(rotation=25)
        # node_timeseries(Node_times,Node_counts,nodes)
        if True:
            plt.savefig(outPath + 'CPS_' + save_string + '.png', format='png', dpi=200)
    except:
        print date
    if False:
        plt.show()
    plt.close('all')
