import numpy as np
import matplotlib.pyplot as plt
import os
import datetime as dt

inPath = '/Volumes/Ian External HD/Node Data/Alarm_Analysis/Figures/correlationCPS/'
dbFiles = [x for x in os.listdir(inPath) if '.png' in x]

for dbFile in dbFiles:
    i0    = dbFile.split('_')
    year  = i0[3].split('-')[1]
    month = i0[3].split('-')[2]
    day   = i0[3].split('-')[3]
    hour_start   = i0[4].split('-')[1]
    minute_start = i0[4].split('-')[2]
    second_start = i0[4].split('-')[3]
    hour_end     = i0[5].split('-')[0]
    minute_end   = i0[5].split('-')[1]
    second_end   = i0[5].split('-')[2].split('.')[0]
