'''
Creates individual HDF5s for cross correlation events from Figures
'''

import numpy as np
import matplotlib.pyplot as plt
import os
import time as t
import h5py
import datetime as dt

inPath   = '/Volumes/Ian External HD/Node Data/Alarm_Analysis/Figures/correlationCPS/'
outPath  = '/Volumes/Ian External HD/Node Data/Alarm_Analysis/Figures/'
nodePath = '/Volumes/Ian External HD/Node Data/Alarm_Analysis/'
dbFiles  = [x for x in os.listdir(inPath) if '.png' in x]
nodes    = ['MUSE01','MUSE04','MUSE06','MUSE10','MUSE11','MUSE12']
starts = []; ends = []
for dbFile in dbFiles:
    i0    = dbFile.split('_')
    year  = int(i0[3].split('-')[1])
    month = int(i0[3].split('-')[2])
    day   = int(i0[3].split('-')[3])
    hour_start   = int(i0[4].split('-')[1])
    minute_start = int(i0[4].split('-')[2])
    second_start = int(i0[4].split('-')[3])
    hour_end     = int(i0[5].split('-')[0])
    minute_end   = int(i0[5].split('-')[1])
    second_end   = int(i0[5].split('-')[2].split('.')[0])
    #
    date_start = dt.datetime(year,month,day,hour_start,minute_start,second_start)
    date_end   = dt.datetime(year,month,day,hour_end,minute_end,second_end)
    #
    posix_start = t.mktime(date_start.timetuple())-200
    posix_end = t.mktime(date_end.timetuple())+200
    #
    starts.append(posix_start)
    ends.append(posix_end)
    #
    hdf5 = h5py.File(outPath+dbFile.replace('png','hdf5'),'w')
    #
    for node in nodes:
        node_files = [x for x in os.listdir(nodePath+node+'/hdf5Files/') if 'ALG' not in x]
        if ((int(month)<10) or (int(day)<10)):
            if int(day)<10:
                if int(month)<10:
                    node_files = [x for x in node_files if ('%s-0%s-0%s' % (str(year), str(month), str(day)) in x)]
                else:
                    node_files = [x for x in node_files if ('%s-%s-0%s' % (str(year), str(month), str(day)) in x)]
            else:
                    node_files = [x for x in node_files if ('%s-0%s-%s' % (str(year), str(month), str(day)) in x)]
        else:
            node_files = [x for x in node_files if ('%s-%s-%s'%(str(year),str(month),str(day)) in x)]

        check = True
        for node_file in node_files:
            f = h5py.File(nodePath+node+'/hdf5Files/'+node_file,'r')

            time = np.array(f['2x4x16Times'])

            if posix_end < time[-1] and check:
                spectra = np.array(f['2x4x16Spectra'])
                livetime = np.array(f['2x4x16LiveTimes'])
                event_index = np.where( (time>posix_start)&(time<posix_end) )[0]
                # writing to output file
                dset = hdf5.create_dataset(node+'_2x4x16Times', (len(time[event_index]), ), dtype='f', compression="gzip")
                dset[...] = time[event_index]
                dset = hdf5.create_dataset(node+'_2x4x16LiveTimes', (len(livetime[event_index]),), dtype='f', compression="gzip")
                dset[...] = livetime[event_index]
                dset = hdf5.create_dataset(node+'_2x4x16Spectra', spectra[event_index].shape, dtype='f', compression="gzip")
                dset[...] = spectra[event_index]
                check = False
            else:
                pass
import pandas as pd
df = pd.DataFrame({"start": starts, "end": ends})
df.to_csv(outPath+"cross_correlation_times.csv", index=False)
