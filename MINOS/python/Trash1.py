import numpy as np
import h5py
import datetime as dt

class create_n42(object):
    def convert_hdf5(self,full_path_to_input,outPath,start_time=0.0,end_time=1E16):
        '''
        Description: Converts .hdf5 from Dr. Nuclear (Nicholson) to .n42

        Parameters
        ----------
        full_path_to_input:  {string} full path to .hdf5 file
        outPath:             {string} full path to output .n42 file
        start_time:          {int}    POSIX start time for .n42 conversion
        end_time:            {int}    POSIX end time for .n42 conversion

        Returns .n42 file at outPath
        -------
        '''
        hdf5 = h5py.File(full_path_to_input,'r')
        #
        '''
        self.spectra  = hdf5['2x4x16Spectra']
        self.time     = hdf5['2x4x16Times']
        self.livetime = hdf5['2x4x16LiveTimes']
        '''
        det_name = 2
        self.spectra  = hdf5['WIND'+str(det_name)+'_Spectra']
        self.time     = hdf5['WIND'+str(det_name)+'_Times']
        self.livetime = hdf5['WIND'+str(det_name)+'_LiveTimes']
        #
        if start_time!=0.0 or end_time!=1E16:
            start,end = np.where(np.array(self.time)>=start_time)[0][0],np.where(np.array(self.time)<=end_time)[0][-1]
            self.spectra   = np.array(self.spectra[start:end])
            self.time      = np.array(self.time[start:end])
            self.livetime  = np.array(self.livetime[start:end])
        #
        self.filename = full_path_to_input.split('/')[-1].replace('.hdf5','.n42')
        self.node_name = str(full_path_to_input.split('/')[-1].split('-')[0])
        #
        #   Initializing .n42 header prior to creating list mode spectra
        #
        f = open(outPath+'WIND'+str(det_name)+'_'+self.filename, "w+")
        #
        f.write('<RadInstrumentData n42DocUUID="b53e471b-28c4-497e-9810-ae7299e230fc" xmlns="http://physics.nist.gov/N42/2011/N42">\n')
        f.write('<RadInstrumentInformation id="sensor9b7996618a05e5a953de0b7345c415f">	<RadInstrumentManufacturerName>Canberra</RadInstrumentManufacturerName>\n')
        f.write('\t<RadInstrumentIdentifier>%s</RadInstrumentIdentifier>\n'%self.node_name)
        f.write('</RadInstrumentInformation><RadDetectorInformation id="gamma9b7996618a05e5a953de0b7345c415f">		<RadDetectorCategoryCode>Gamma</RadDetectorCategoryCode>\n')
        f.write('\t<RadDetectorKindCode>NaI</RadDetectorKindCode>\n')
        f.write('\t<RadDetectorLengthValue>0.0000001</RadDetectorLengthValue>\n')
        f.write('\t<RadDetectorDiameterValue>0.0000001</RadDetectorDiameterValue>\n')
        f.write('</RadDetectorInformation>\n')
        f.write('<EnergyCalibration id="EnCal9b7996618a05e5a953de0b7345c415f">	<CoefficientValues>3.89694510e+00 1.40542060e+00 5.41037604e-05</CoefficientValues>\n')
        f.write('</EnergyCalibration>\n')
        # Beginning measurements (1-second intervals)
        self.time = [self.time[-1]]
        self.spectra = [np.sum(self.spectra,axis=0)]
        self.livetime = [sum(self.livetime)]
        measurement = self.RadMeasurement(self.time,self.spectra,self.livetime)
        #
        [f.write(line) for line in measurement]
        #
        f.write('</RadInstrumentData>')
        f.close()
        return

    def RadMeasurement(self,time_i,spectra_i,livetime_i):
        '''
        Description: Writes each sample (1-second) spectra in RadMeasurement .xml format
        Parameters
        ----------
        time
        spectra
        livetime

        Returns
        -------

        '''
        string = []
        for m in range(len(time_i)):
            t,s,l = time_i[m],spectra_i[m],livetime_i[m]
            measurement_name = str(self.node_name)+'_'+str(t)
            string.append('<RadMeasurement id="%s">\n'%(measurement_name))
            string.append('\t<MeasurementClassCode>Foreground</MeasurementClassCode>\n')
            date_m = dt.datetime.fromtimestamp(t)
            year,month,day,hour,minute,second=str(date_m.year),str(date_m.month),str(date_m.day),str(date_m.hour),str(date_m.minute),str(date_m.second)
            #
            if month<10: month='0'+month
            if day<10: day='0'+day
            if hour<10: hour='0'+hour
            if minute<10: minute='0'+minute
            if second<10: second='0'+second
            #
            string.append('\t<StartDateTime>%s-%s-%sT%s:%s:%s.000Z</StartDateTime>\n'%(year,month,day,hour,minute,second))
            string.append('\t<RealTimeDuration>PT1.0S</RealTimeDuration>\n')
            string.append('\t<Spectrum id="%.1f" radDetectorInformationReference="gamma9b7996618a05e5a953de0b7345c415f" energyCalibrationReference="EnCal9b7996618a05e5a953de0b7345c415f">\n'%float(t))
            string.append('\t\t<LiveTimeDuration>PT1.0S</LiveTimeDuration>\n')
            string.append('\t\t<ChannelData> ')
            # Inserting spectrum for interval
            for bin in s:
                string.append(' %i'%bin)
            string.append('</ChannelData>\n\t</Spectrum>\n')
            #
            string.append('</RadMeasurement>\n')
        return string

if __name__ == '__main__':
    #
    import os
    full_input_path  = '/Volumes/IAN USB/WIND/GADRAS dat Creation/2018-11-21/'
    full_output_path = full_input_path+'WIND2/'
    start_posix_time = 1533098110.0
    end_posix_time   = 1533099010.0
    # Call class
    n42 = create_n42()
    # Call class function
    #n42.convert_hdf5(full_input_path,full_output_path,start_posix_time,end_posix_time)
    filename = [x for x in os.listdir(full_input_path) if '.hdf5' in x]
    #
    for f in filename:
        n42.convert_hdf5(full_input_path+f,full_output_path)