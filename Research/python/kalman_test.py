import numpy as np
import matplotlib.pyplot as plt
from itertools import izip, count
import matplotlib.pyplot as plt
import pandas as pd
import collections
from numba import jit
import h5py

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
    m = y1.size
    n = x2.size - 1
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

    # first(0th) X1 bin
    ibin = 0
    i_place_lower = i_place[ibin]
    i_place_upper = i_place[ibin + 1]
    if i_place_lower == 0 and i_place_upper > 0:
        if i_place_upper == 1:
            y2_index = i_place_upper - 1
            y2[y2_index] += y1[ibin]

        if i_place_upper == 2:
            print 'first bin 2 overlaps'

    # redistributing X1 bins with [1,m-1] indeces
    for ibin in range(1, m - 1):
        if y1[ibin] == 0:
            continue
        i_place_lower = i_place[ibin]
        i_place_upper = i_place[ibin + 1]
        x1min = x1[ibin]
        x1max = x1[ibin + 1]
        # Stop at the end
        if i_place_lower >= n - 2:
            return y2

        # X1 bin fully inside X2 bin:
        if i_place_lower == i_place_upper:
            y2_index = i_place_upper - 1
            y2[y2_index] += y1[ibin]

        # X1 bin overlaps w/ two X2 bins:
        if i_place_lower == i_place_upper - 1:
            # X2 value that "splits" the X1 bin
            x2_ov1 = x2[i_place_lower]
            # Roll the dice y1[ibin]-times
            for irand in range(0, int(y1[ibin])):
                probValue = np.random.uniform(x1min, x1max)
                if probValue < x2_ov1:
                    # send photon to lower bin :))
                    y2_index = i_place_upper - 2
                    y2[y2_index] += 1
                else:
                    y2_index = i_place_upper - 1
                    y2[y2_index] += 1
        # X1 bin overplaps w/ three X2 bins
        if i_place_lower == i_place_upper - 2:
            print '2 overlap bins'

    return y2

class LP_Filter_Detector(object):
    def __init__(self,data,window_size=10,sigma_value=5):
        self.data = data
        self.window_size = window_size
        self.sigma_value = sigma_value

    def moving_average(self):
        window = np.ones(int(self.window_size))/float(self.window_size)
        return np.convolve(self.data, window, 'same')

    @jit
    def explain_anomalies(self):
        avg = self.moving_average().tolist()
        residual = self.data - avg
        # Calculate the variation in the distribution of the residual
        std = np.std(residual)
        return {'standard_deviation': round(std, 3),
                'anomalies_dict': collections.OrderedDict([(index, y_i) for
                                                           index, y_i, avg_i in izip(count(), self.data, avg)
                  if (y_i > avg_i + (self.sigma_value*std)) | (y_i < avg_i - (self.sigma_value*std))])}

    def explain_anomalies_rolling_std(self):
        avg = self.moving_average()
        avg_list = avg.tolist()
        residual = self.data - avg
        # Calculate the variation in the distribution of the residual
        testing_std = pd.rolling_std(residual, self.window_size)
        testing_std_as_df = pd.DataFrame(testing_std)
        rolling_std = testing_std_as_df.replace(np.nan,
                                      testing_std_as_df.ix[self.window_size - 1]).round(3).iloc[:,0].tolist()
        std = np.std(residual)
        return {'stationary standard_deviation': round(std, 3),
                'anomalies_dict': collections.OrderedDict([(index, y_i)\
                  for index, y_i, avg_i, rs_i in izip(count(),self.data, avg_list, rolling_std)\
                  if (y_i > avg_i + (self.sigma_value * rs_i)) | (y_i < avg_i - (self.sigma_value * rs_i))])}

    def alarm_indexes(self):
        events = self.explain_anomalies()
        x_anomaly = np.fromiter(events['anomalies_dict'].iterkeys(), dtype=int, count=len(events['anomalies_dict']))
        #
        # Finding sequential alarm events and grouping together if less than <time_concatenate> seconds apart
        #
        time_concatenate = 12
        diff = np.where(np.array([x_anomaly[i]-x_anomaly[i-1] for i in range(1,len(x_anomaly))])>time_concatenate)[0]
        loc = 0; grouped_events = []
        for i in diff:
            if loc == 0:
                grouped_events.append(np.array(x_anomaly[0:i + 1]))
            elif len(x_anomaly[loc:i]) == 0:
                grouped_events.append(np.array(x_anomaly[loc]))
            elif i == diff[-1]:
                grouped_events.append(np.array(x_anomaly[i - 1::]))
            else:
                grouped_events.append(x_anomaly[loc:i + 1])
            loc = i + 1
        #
        # Adding a few seconds on either end of alarm
        #
        final_events = []
        for event in grouped_events:
            if np.array(event).shape == ():
                new_event = np.insert(np.array(event), 0, event - 2)
                new_event = np.insert(new_event, len(new_event), new_event[-1] + 2)
                final_events.append(new_event)
            elif (event[0] > 10):
                new_event = np.insert(np.array(event), 0, event[0] - 2)
                new_event = np.insert(new_event, len(new_event), new_event[-1] + 4)
                final_events.append(new_event)
        #
        return final_events

    @jit(nopython=False)
    def plot_results(self,x, y,text_xlabel="X Axis", text_ylabel="Y Axis", applying_rolling_std=False,show_plot=False):
        fig,ax = plt.subplots(figsize=(13,6))
        ax.plot(x, y, "k.")
        y_av = self.moving_average()
        ax.plot(x, y_av, color='green')
        ax.set_xlabel(text_xlabel)
        ax.set_ylabel(text_ylabel)

        # Query for the anomalies and plot the same
        events = {}
        if applying_rolling_std:
            events = self.explain_anomalies_rolling_std()
        else:
            events = self.explain_anomalies()

        x_anomaly = np.fromiter(events['anomalies_dict'].iterkeys(), dtype=int, count=len(events['anomalies_dict']))
        y_anomaly = np.fromiter(events['anomalies_dict'].itervalues(), dtype=float,
                                                count=len(events['anomalies_dict']))
        ax.plot(x_anomaly, y_anomaly, "r*", markersize=12)
        events = self.alarm_indexes()
        [plt.axvspan(x[0], x[-1], alpha=0.35, color='red') for x in events]

        # add grid and lines and enable the plot
        ax.grid(True)
        ax.set_axisbelow(True)

        if show_plot:
            plt.show(True)

        return

    def make_alarm_plots(self,spectra,time,livetime,events,path,energy_pairs=np.arange(1,1025)):
        import datetime as dt
        import matplotlib.dates as md
        #
        events = self.alarm_indexes()
        #
        for i in events:
            #
            time_added = 60
            CPS = np.sum(spectra[i[0]-time_added:i[-1]+time_added],axis=1)
            CPS_times = np.array([dt.datetime.fromtimestamp(x) for x in time[i[0]-time_added:i[-1]+time_added]])
            #
            alarm_range = np.arange(i[0]-time_added,i[-1]+time_added)
            start = CPS_times[np.where(alarm_range >= i[0])[0][0]]       # Alarm start for shading
            end_index = np.where(alarm_range  >= i[-1])[0][0]
            if end_index >= len(CPS_times):
                end_index = len(CPS_times)
            end   = CPS_times[end_index]        # Alarm end for shading
            #
            spectrum_alarm       = np.sum(spectra[i[0]:i[-1]],axis=0)/sum(livetime[i[0]:i[-1]])
            background_timeframe = np.arange(i[0]-60,i[0]-30)
            background_spectrum  = np.sum(spectra[background_timeframe],axis=0)/sum(livetime[background_timeframe])
            spectrum_residual    = np.abs(spectrum_alarm - background_spectrum)
            #
            #
            fig,ax = plt.subplots(2,1,figsize=(12,8))
            ax[0].plot(CPS_times,CPS,label='CPS')
            ax[0].axvspan(start,end, alpha=0.45, color='red',label='Alarm Region')
            ax[0].legend(loc='upper right',fancybox=False,shadow=True)
            ax[0].grid(True,alpha=0.5); ax[0].set_axisbelow(True)
            ax[0].set_ylim(550,2500)
            ax[0].set_xlabel('Time (Hour:Minute:Second)'); ax[0].set_ylabel('Counts per second')
            title_string = 'Alarm Event on %i-%i-%i at %i:%i:%i'%(start.month,start.day,start.year,start.hour,start.minute,start.second)
            ax[0].set_title(title_string)
            date_formatter = md.DateFormatter('%H:%M:%S')
            ax[0].xaxis.set_major_formatter(date_formatter)
            #
            ax[1].plot(energy_pairs,background_spectrum,'b',label='Background',linewidth=1.2)
            ax[1].plot(energy_pairs,spectrum_alarm,'r',label='Alarm',linewidth=1)
            ax[1].plot(energy_pairs,spectrum_residual,'g',label='Abs. Residual',linewidth=0.8)
            ax[1].legend(loc='upper right', fancybox=False, shadow=True)
            ax[1].grid(True, alpha=0.5)
            ax[1].set_axisbelow(True)
            ax[1].set_yscale('log')
            ax[1].set_ylabel('Time-normalized\nFrequency')
            ax[1].set_xlabel('Energy (keV)')
            #
            date0 = dt.datetime.fromtimestamp(time[i[0]])
            day,month,year,hour,minute,second=date0.day,date0.month,date0.year,date0.hour,date0.minute,date0.second
            title_string = 'Alarm_%i-%i-%i_at_%i-%i-%i.png'%(day,month,year,hour,minute,second)
            plt.show()
            plt.savefig(path+title_string,dpi=200)
            plt.close('all')
            #

if __name__ == '__main__':
    import os
    filename = '/Volumes/Ian External HD/RSL NOVArray/data/August_2018/Aug1-Aug8_2018/Full Sensor Reports.2018-08-01 04_00_00+0000 - 2018-08-08 04_00_00+0000.4 rsl Sensors (4).hdf5'
    spectra = np.array(h5py.File(filename)['digiBASE-RH 17178703/RadiationReading/spectrum/'])
    time = np.array(h5py.File(filename)['digiBASE-RH 17178703/RadiationReading/time/'])
    livetime = np.array(h5py.File(filename)['digiBASE-RH 17178703/RadiationReading/real_time/'])
    CPS = np.sum(spectra,axis=1)
    times = np.arange(1,len(CPS)+1)
    #
    # LP_Filter_Detection(DATA,WINDOW_SIZE,SIGMA_VALUE)
    #
    LPF = LP_Filter_Detector(CPS,10,3)
    #
    # Grabbing alarm indexes using previous window_size and sigma_value
    indexes = LPF.alarm_indexes()
    #
    # LoaEneryding energy pair data
    energy_path = '/Volumes/Ian External HD/RSL NOVArray/energy_pairs/'
    energy_files = [x for x in os.listdir(energy_path) if '.csv' in x]
    energy_pairs = pd.read_csv(energy_path+energy_files['17178703' in energy_files],sep=',',header=None)
    save_path = '/Volumes/Ian External HD/RSL NOVArray/data/August_2018/figures/'
    LPF.make_alarm_plots(spectra,time,livetime,indexes,save_path,energy_pairs)
    #
    