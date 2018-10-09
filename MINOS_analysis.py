import numpy as np
import sqlite3

def node_colors(name):
    """
    Description: Provides colors for Nodes for consistency across analysis
    Parameters
    -------
    :para name: {string}    : Name of node (e.g. 'MUSE01')

    Returns
    -------
    :return:    {string}    : String color for matplotlib or seaborn colormap
    """
    colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']
    node_names = ['MUSE01','MUSE04','MUSE06','MUSE10','MUSE11','MUSE12']
    return colors[node_names.index(name)]


def node_timeseries(data_x,data_y,data_labels):
    """
    Description: Creates 3-d plot for nodes (Counts,Time,Node name)
    Parameters
    -------
    :param  data_x:      {dictionary}    : data along x-axis to plot (e.g. CPS)
    :param  data_y:      {dictionary}    : data along y-axis to plot (e.g. Datetime)
    :param  data_labels: {list/array}    : labels for each component in data
    :return:
    """
    import datetime as dt
    import matplotlib.dates as md
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    n = len(data_labels)
    fig = plt.figure()
    ax = plt.subplot(projection='3d')
    counter = n+1
    for i in range(n):
        values = data_x[data_labels[i]] # Dates
        values = [dt.datetime.fromtimestamp(ts) for ts in values]
        x = np.array([md.date2num(v) for v in values])
        z = data_y[data_labels[i]]      # Counts
        y = np.ones(x.size)*(i+1)       # Node Label (name)
        ax.plot(x, y, z, zorder = counter)
        if False:
            ax.set_xlabel('Time\n(Hour:Minute:Second)',fontweight='bold')
        else:
            ax.set_xlabel('Time\n(Month/Day)',fontweight='bold')
        ax.set_ylabel('Node Name',fontweight='bold')
        ax.set_zlabel('CPS',fontweight='bold')
        ax.xaxis.labelpad = 12
        ax.yaxis.labelpad = 10
        ax.zaxis.labelpad = 5
        ax1 = plt.gca()
        if True:
            xfmt = md.DateFormatter('%H:%M:%S')
        else:
            xfmt = md.DateFormatter('%m/%d')
        ax1.xaxis.set_major_formatter(xfmt)
        counter=counter-1

    ax.axes.set_yticks(range(1,n+1))
    ax.axes.set_yticklabels(data_labels)
    return fig

def ksigma_RIO_alg(data,fore,back,threshold,isotope):
    '''
    Description: Performs alarm algorithm
    :param data: {array} : array of spectral data
    :param fore: {int}   : foreground length
    :param back: {int}   : background length
    :param threshold: {int} : ratio threshold to alarm

    :return roi_counts:  {list}  : counts in ROI window
    :return alarm_index: {array} : index of alarm
    '''
    if isotope == ('Ar-41' or 'Ar41' or '41Ar'):
        window = tuple((440,470))
        background_lower = tuple((window[0]-60,window[0]-10))
        background_upper = tuple((window[1]+80,window[1]+120))
        x_coords = (background_lower[0],background_lower[1],background_upper[0],background_upper[1])
    if isotope ==  ('rain'):
        window = []
    if isotope == ('spectra' or 'full spectrum' or 'spectrum' or 'full'):
        isotope = 'full'
    alarm_index = []
    s = []
    roi_counts = []
    for j in data:

        if isotope not in 'full':
            peak_values = j[window[0]:window[1]]
            roi_counts.append(sum(peak_values))
            y_coords = [j[x] for x in x_coords]
            linear_coef = np.polyfit(x_coords,y_coords,1)
            fit = np.poly1d(linear_coef)
            fitted_data = fit(np.linspace(window[0],window[1]))
            value = sum([x-y for x,y in zip(peak_values,fitted_data)])
        else:
            value = sum(j)
            roi_counts.append(value)

        if value > 0:
            s.append(value)
        else:
            s.append(0)

    dead_zone = 5 #separation between for and back

    alarm_ratio = []
    for i in range(len(s)):
        if i > (fore+back+5):
            foreground = np.average(s[(i-fore):i])
            background = np.average(s[(i-(back+fore+dead_zone)):(i-(fore+dead_zone))])
            if 'background_std' not in vars():
                background_std = np.std(s[(i-(back+fore+dead_zone)):(i-(fore+dead_zone))])
            ratio = (float(foreground)-float(background))/background_std
            alarm_ratio.append(ratio)

            if ratio>threshold:
                alarm_index.append(i)
            else:
                background_std = np.std(s[(i-(back+fore+dead_zone)):(i-(fore+dead_zone))])

    return roi_counts,alarm_index,alarm_ratio

def fixed_threshold_alg(data,time,threshold):
    """
    Description: Performs alarm algorithm

    Parameters
    ----------
    :param data: {array} : array of spectral data
    :param threshold: {int} : ratio threshold to alarm
    :param time:      {int} : initial time for calc fixed threshold

    Returns
    -------
    :return std:    {int}   : CPS thresholds
    """

    s = []
    for j in data[0:time]:
        s.append(sum(j))
    std = np.std(s)
    alarm_value = threshold * std
    alarm_index = [x for x in range(len(data)) if (alarm_value <= data[x])]
    return alarm_index

def sprt_alg(counts,livetime,fore,back,threshold):
    """

    Parameters
    ----------
    :param counts:      {array} : CPS data
    :param livetime:    {array} : Live time data
    :param fore:        {int}   : foreground window size
    :param back:        {int}   : background window size
    :param threshold:   {int}   : threshold parameter

    Returns
    -------
    :return sprt:   {list}  : SPRT values from LLR
    :return indexs: {list}  : Data indexs for SPRT>threshold
    """
    dead_zone = 3; indexs = []; sprt_values = []; sprt_old = 0
    for i in range(len(counts)):
        if i > (fore+back+dead_zone):
            frgC = np.sum(counts[(i-fore):i])
            frgLT = np.sum(livetime[(i-fore):i])
            bkgrdC = np.sum(counts[(i-(back+fore+dead_zone)):(i-(fore+dead_zone))])
            bkgrdLT = np.sum(livetime[(i-(back+fore+dead_zone)):(i-(fore+dead_zone))])
            foreground = float(frgC/frgLT)
            background = float(bkgrdC/bkgrdLT)
            sprt=(foreground-background)-foreground*np.log10(foreground/background)+sprt_old
            if sprt<0:
                sprt = 0.0
            sprt_values.append(sprt)
            sprt_old = sprt
            if sprt>=threshold:
                indexs.append(i)
    #print 'len(sprt): %i\type(sprt): %s'%(len(sprt_values),type(sprt_values))
    return sprt_values,indexs

def gausswLine(x,a,b,c,mu,sigma): # linear + gaussian
    """

    Parameters
    ----------
    x
    a
    b
    c
    mu
    sigma

    Returns
    -------

    """
    return (
        #line
        a+c*(x-mu)+
        #photopeak
        (b/(sigma*np.sqrt(2.0*np.pi)))*np.exp(-(x-mu)**2/(2.0*sigma**2))
        )

def linFunc(x,b):
    '''

    Parameters
    ----------
    x
    b

    Returns
    -------

    '''
    return b*x

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

def rebin_Archer(x1, y1, x2):
    """
    Rebin histogram values y1 from old bin edges x1 to new edges x2.

     :param x1: : m+1 array of old bin edges.
     :param y1: m array of old histogram values. This is the total number in each bin, not an average.
     :param x2: : n+1 array of new bin edges.

    :return y2: n array of rebinned histogram values.

    The rebinning algorithm assumes that the counts in each old bin are
    uniformly distributed in that bin.

    Bins in x2 that are entirely outside the range of x1 are assigned 0.
    """

    x1 = np.asarray(x1)
    y1 = np.asarray(y1)
    x2 = np.asarray(x2)


    # allocating y2 vector
    m  = y1.size
    n  = x2.size - 1
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

    #print len(x1),len(y1),len(y2)
    #print len(i_place)
    #print i_place[980:len(i_place)]

    # first(0th) X1 bin
    ibin=0
    i_place_lower=i_place[ibin]
    i_place_upper=i_place[ibin+1]
    if i_place_lower==0 and i_place_upper>0:
        if i_place_upper==1:
            y2_index = i_place_upper - 1
            y2[y2_index] += y1[ibin]
            #print 'First bin %d counts go to y2 bin-%d' % (y1[ibin],y2_index)
            #print y2
        if i_place_upper==2:
            print 'first bin 2 overlaps'

    # redistributing X1 bins with [1,m-1] indeces
    for ibin in range(1,m-1):
        if y1[ibin]==0:
            continue
        i_place_lower=i_place[ibin]
        i_place_upper=i_place[ibin+1]
        x1min = x1[ibin]
        x1max = x1[ibin+1]
        #Stop at the end
        if i_place_lower >= n-2:
            return y2

        # X1 bin fully inside X2 bin:
        if i_place_lower == i_place_upper:
            y2_index = i_place_upper - 1
            #print 'y2_index %d -- i_place_lower=%d -- ibin=%d' % (y2_index,i_place_lower, ibin)
            y2[y2_index] += y1[ibin]
            #print '%d-th bin %d counts go to y2 bin-%d: X1 overlaps with X2' % (ibin,y1[ibin],y2_index)
            #print y2
        # X1 bin overlaps w/ two X2 bins:
        if i_place_lower == i_place_upper - 1:
            # X2 value that "splits" the X1 bin
            x2_ov1 = x2[i_place_lower]
            # Roll the dice y1[ibin]-times
            for irand in range(0, int(y1[ibin])):
                probValue = np.random.uniform(x1min,x1max)
                #print 'rand-%d probV=%f x2_ov1=%f' %(irand,probValue,x2_ov1)
                if probValue < x2_ov1:
                    # send photon to lower bin :))
                    y2_index = i_place_upper - 2
                    y2[y2_index] += 1
                else:
                    y2_index = i_place_upper - 1
                    y2[y2_index] += 1
            #print '%d-th bin %d counts go to y2 bin-%d: X1 is completely in X2' % (ibin,y1[ibin],y2_index)
            #print y2
        # X1 bin overplaps w/ three X2 bins
        if i_place_lower == i_place_upper - 2:
            print '2 overlap bins'

    return y2

if __name__ == '__main__':
    print '\n<<< No script to run for python __main__ >>>'
    pass