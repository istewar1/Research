import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as md
import os
import pandas as pd
import numpy as np

def gainstab_csv(inPath,csv_file):
    '''
    Description: Returns matplotlib figure from .csv
    :param inPath:          /path/to/file
    :param csv_file:        name of .csv file

    :return: matplotlib figure of K40 & Th232 channel
    '''
    df = pd.read_csv(inPath+csv_file,header=0)
    fig, ax = plt.subplots(nrows=3,sharex=True)
    dates = [dt.datetime.fromtimestamp(ts) for ts in df.time]
    datetimes = np.array([md.date2num(v) for v in dates])
    ax[0].plot(datetimes,df['total_counts'],label='Avg.10s CPS')
    ax[0].set_ylabel('CPS')
    ax[0].grid(alpha=0.5)
    xfmt = md.DateFormatter('%m/%d %H:%M')
    ax[0].xaxis.set_major_formatter(xfmt)
    ax[0].legend(loc='upper right')
    ax[1].plot(datetimes,df['k40_peak'],label='K-40')
    ax[1].set_ylabel('Peak Channel No.')
    ax[1].xaxis.set_major_formatter(xfmt)
    ax[1].grid(alpha=0.5)
    ax[1].legend(loc='upper right')
    ax[2].plot(datetimes,df['th232_peak'],label='Th-232')
    ax[2].set_ylabel('Peak Channel No.')
    ax[2].xaxis.set_major_formatter(xfmt)
    ax[2].grid(alpha=0.5);
    ax[2].legend(loc='upper right')
    ax[2].set_xlabel('Date (Month/Day  Hour:Minute)')
    plt.setp( ax[1].xaxis.get_majorticklabels(), rotation=25 )
    plt.subplots_adjust(hspace=0.)
    return fig

if __name__ == '__main__':
    inPath = '/Users/i6o/Downloads/'
    dbFiles = [x for x in os.listdir(inPath) if '.csv' in x]
    fig = gainstab_csv(inPath,dbFiles[0])
    plt.show()