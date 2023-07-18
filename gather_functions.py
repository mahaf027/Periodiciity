'''
This document contains the functions that acquire data from sources and create LightCurve instances
(or lists of LightCurve instances) from them.
'''

#Imports
import os
import csv
from datetime import datetime
from datetime import timedelta
from datetime import timezone
import numpy as np
import cdflib
from cdasws import CdasWs

from lightcurve_class import LightCurve
from utility import read_l2_hres, lc_sort_func

os.environ['CDF_LIB'] = './lib'

#Constants
SPEED_OF_LIGHT = 3e8 #m/s
METERS_PER_AU = 149597870700

#Gather functions
def gather_psp(start_time, end_time, power_variable, auto='averages', shift_1_au=False):
    '''
    Scrapes PSP data from cdaweb and generates a list of LightCurve objects in the given time
    constraints.
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The offset-aware UTC start of the time interval of interest.
    end_time : datetime.Datetime
        The offset-aware UTC end of the time interval of interest.
    power_variable : String
        Specifies which variable to collect. 'V1V2', 'V3V4', or 'combination'.
    auto : String (optional)
        Specifies whether to collect auto_averages or auto_peaks data.
    shift_1_au : Boolean (optional)
        Toggles whether start and end times are shifted earlier to account for PSP being closer to
        the sun than Earth is.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        The list of LightCurve objects created from the PSP data.
    '''
    cdas = CdasWs()
    #Calculating the PSP --> 1 au signal time offset if shift1au specified
    offset = timedelta(seconds=0)
    #Position is only reported daily, so we get the position on the day of start_time
    rad_start = datetime(start_time.year, start_time.month, start_time.day, 0, 0,
                            tzinfo=timezone.utc)
    rad_end = rad_start + timedelta(days=2)
    if shift_1_au:
        status, data = cdas.get_data('PSP_HELIO1DAY_POSITION', ['RAD_AU'], rad_start, rad_end)
        if status['http']['status_code'] == 200:
            rad = data['RAD_AU'][0]
            dist_m = (1 - rad) * METERS_PER_AU
            dist_light_sec = dist_m / SPEED_OF_LIGHT
            offset = timedelta(seconds=dist_light_sec)
        else:
            print('Connection to cdaweb unsuccessful. No offset subtracted.')
    #Adjusting start_time and end_time accordingly
    start = start_time - offset
    end = end_time - offset
    freq_ranges = ['hfr', 'lfr']
    lc_list = []
    v12 = f'_auto_{auto}_ch0_V1V2'
    v34 = f'_auto_{auto}_ch1_V3V4'
    seed = 'psp_fld_l2_rfs_'
    for freq in freq_ranges:
        status, data = cdas.get_data(f'{seed}{freq}'.upper(), [f'{seed}{freq}{v12}',
                                        f'{seed}{freq}{v34}'], start, end)
        if status['http']['status_code'] == 200:
            ts12 = [dt.replace(tzinfo=timezone.utc).timestamp()
                    for dt in data[f'epoch_{freq}{v12}']]
            ts34 = [dt.replace(tzinfo=timezone.utc).timestamp()
                    for dt in data[f'epoch_{freq}{v34}']]
            power12 = np.transpose([p for p in data[f'{seed}{freq}{v12}']])
            power34 = np.transpose([p for p in data[f'{seed}{freq}{v34}']])
            #Frequencies are given in Hz, we want them in MHz
            freqs12 = [f * 1e-6 for f in data[f'frequency_{freq}{v12}'][0]]
            freqs34 = [f * 1e-6 for f in data[f'frequency_{freq}{v34}'][0]]
            if power_variable == 'V1V2':
                for i, freq in enumerate(freqs12):
                    lc_list.append(LightCurve(ts12, power12[i], f'PSPV12{auto[0]}',
                                                input_frequency=freq))
            elif power_variable == 'V3V4':
                for i, freq in enumerate(freqs34):
                    lc_list.append(LightCurve(ts34, power34[i], f'PSPV34{auto[0]}',
                                                input_frequency=freq))
            elif power_variable == 'combination':
                #The frequency lists and timestamps should always be the same for both channels
                #We assume this and arbitrarily pick the V1V2 versions
                for i, freq in enumerate(freqs12):
                    powers = [power12[i], power34[i]]
                    powers = np.transpose(powers)
                    power_comb = [np.sqrt((p[0] ** 2) + (p[1] ** 2)) for p in powers]
                    lc_list.append(LightCurve(ts12, power_comb, f'PSPComb{auto[0]}',
                                                input_frequency=freq))
            else:
                print("Please enter a valid power variable. 'V1V2', 'V3V4', or 'combination'")
        else:
            print('Connection to cdaweb failed')
    lc_list.sort(key=lc_sort_func)
    return lc_list

def gather_wind_rad1(start_time, end_time):
    '''
    Reads WIND RAD1 data from files and generates a list of LightCurve objects in the given time
    constraints.
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The offset-aware UTC start of the time interval of interest.
    end_time : datetime.Datetime
        The offset-aware UTC end of the time interval of interest.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        The list of LightCurve objects created from the RAD1 data.
    '''
    file = f"./DATA/WIND/wi_wa_rad1_l2_{start_time.strftime('%Y%m%d')}_v01.dat"
    header, data = read_l2_hres(file)
    #In RAD1 data, some frequencies are measured 4 times per sweep. We want data from only those
    #frequencies in order to have the highest resolution possible
    all_freqs = data[0]['FREQ']
    freqs = {}
    for idx, freq in enumerate(all_freqs):
        if all_freqs.count(freq) == 4:
            freqs[idx] = freq #kHz
    #Within each sweep, individual frequencies are measured in bursts of 8 measurements
    #This is important later
    #Storing the timestamps of the beginning of each sweep and all of the power measurements, as
    #well as the time within the sweep of each measurement
    timestamps = []
    power = []
    offsets = []
    for idx, _ in enumerate(header):
        year = header[idx]['YEAR']
        month = header[idx]['MONTH']
        day = header[idx]['DAY']
        hour = header[idx]['HOUR']
        minute = header[idx]['MINUTE']
        second = header[idx]['SECOND']
        msmt_time = datetime(year, month, day, hour, minute, second, tzinfo=timezone.utc)
        if msmt_time >= start_time - timedelta(minutes=3) and msmt_time <= end_time:
            timestamps.append(msmt_time.timestamp())
            power.append(data[idx]['VSPAL'][::2])
            offsets.append(data[idx]['TSPAL'][::2])
    #Creating LightCurves for each frequency measured 4 times per sweep
    lc_list_v0 = []
    freq_indices = freqs.keys()
    for idx in freq_indices:
        temp_ts = []
        temp_power = []
        for i, timestamp in enumerate(timestamps):
            if timestamp + offsets[i][idx*8] >= start_time.timestamp() and \
            timestamp + offsets[i][idx*8] <= end_time.timestamp():
                for j in range(8):
                    #offset = round(np.mean(offsets[i][idx*8:8 + idx*8]))
                    temp_ts.append(timestamp + offsets[i][j + idx*8])
                    temp_power.append(power[i][j + idx*8])
                    #temp_power.append(np.mean(power[i][idx*8:8 + idx*8]))
        lc_list_v0.append(LightCurve(temp_ts, temp_power, 'WIND',
                                     input_frequency=freqs[idx] * 1e-3))
    lc_list_v0.sort(key = lc_sort_func)
    #Now we have 4 LightCurves at each frequency. We want to consolidate these into one each
    lc_list = []
    for i in range(int(len(lc_list_v0) / 4)):
        lc_comb = lc_list_v0[4*i] + lc_list_v0[1 + 4*i] + lc_list_v0[2 + 4*i] + lc_list_v0[3 + 4*i]
        lc_list.append(lc_comb)
    return lc_list

def gather_wind_rad2(start_time, end_time):
    '''
    Reads WIND RAD2 data from files and generates a list of LightCurve objects in the given time
    constraints.
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The offset-aware UTC start of the time interval of interest.
    end_time : datetime.Datetime
        The offset-aware UTC end of the time interval of interest.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        The list of LightCurve objects created from the rad2 data.
    '''
    lc_list = []
    file = f"./DATA/WIND/wi_wa_rad2_l2_{start_time.strftime('%Y%m%d')}_v01.dat"
    header, data = read_l2_hres(file)
    #Setting up timestamps for LightCurve objects
    freqs = []
    power = []
    timestamps = []
    for idx, _ in enumerate(header):
        year = header[idx]['YEAR']
        month = header[idx]['MONTH']
        day = header[idx]['DAY']
        hour = header[idx]['HOUR']
        minute = header[idx]['MINUTE']
        second = header[idx]['SECOND']
        msmt_time = datetime(year,month,day,hour,minute,second,tzinfo=timezone.utc)
        if msmt_time >= start_time and msmt_time <= end_time:
            timestamps.append(msmt_time.timestamp())
            power.append(data[idx]['VSPAL'][::2])
    #Creating LightCurves
    power = np.array(power)
    power = power.transpose()
    freqs = data[0]['FREQ'] #kHz
    for idx, freq in enumerate(freqs):
        lc_list.append(LightCurve(timestamps, power[idx], 'WIND', input_frequency=freq * 1e-3))
    # Sorting LightCurves by frequency
    lc_list.sort(key = lc_sort_func)
    return lc_list

def gather_wind(start_time, end_time):
    '''
    Calls the rad1 and rad2 gather functions to generate a full list of WIND Light Curve objects.
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The offset-aware UTC start of the time interval of interest.
    end_time : datetime.Datetime
        The offset-awrae UTC end of the time interval of interest.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        The list of LightCurve objects created from WIND data.
    '''
    lc_list = gather_wind_rad1(start_time, end_time) + gather_wind_rad2(start_time, end_time)
    return lc_list

def gather_stereo_a(start_time, end_time):
    '''
    Reads STEREO A data from cdaweb and generates a list of LightCurve objects in the given time
    constraints.
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The offset-aware UTC start of the time interval of interest.
    end_time : datetime.Datetime
        The offset-aware UTC end of the time interval of interest.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        The list of LightCurve objects created from STEREO A data.
    '''
    cdas = CdasWs()
    lc_list = []
    #Getting the data from cdaweb
    status, data = cdas.get_data('STEREO_LEVEL2_SWAVES', ['avg_intens_ahead'], start_time,
                                    end_time)
    #Unpacking the data only if the connection was successful
    if status['http']['status_code'] == 200:
        timestamps = [dt.replace(tzinfo=timezone.utc).timestamp() for dt in data['Epoch']]
        freqs = [f * 1e-3 for f in data['frequency']] #MHz
        power_by_freq = np.transpose(data['avg_intens_ahead'])
        for i, freq in enumerate(freqs):
            lc_list.append(LightCurve(timestamps, power_by_freq[i], 'STEREO_A',
                                        input_frequency=freq))
    else:
        print('Connection to cdaweb failed.')
    lc_list.sort(key = lc_sort_func)
    return lc_list

def gather_stereo_b(start_time, end_time):
    '''
    Reads STEREO B data from cdaweb and generates a list of LightCurve objects in the given time
    constraints.
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The offset-aware UTC start of the time interval of interest.
    end_time : datetime.Datetime
        The offset-aware UTC end of the time interval of interest.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        The list of LightCurve objects created from STEREO B data.
    '''
    cdas = CdasWs()
    lc_list = []
    #Getting the data from cdaweb
    status, data = cdas.get_data('STEREO_LEVEL2_SWAVES', ['avg_intens_behind'], start_time,
                                    end_time)
    #Unpacking the data only if the connection was successful
    if status['http']['status_code'] == 200:
        timestamps = [dt.replace(tzinfo=timezone.utc).timestamp() for dt in data['Epoch']]
        freqs = [f * 1e-3 for f in data['frequency']] #MHz
        power_by_freq = np.transpose(data['avg_intens_ahead'])
        for i, freq in enumerate(freqs):
            lc_list.append(LightCurve(timestamps, power_by_freq[i], 'STEREO_B',
                                        input_frequency=freq))
    else:
        print('Connection to cdaweb failed.')
    lc_list.sort(key = lc_sort_func)
    return lc_list

def gather_aia(start_time, end_time, channel, reg):
    '''
    Creates a LightCurve object from csv files created by Reed Masek. The csv files consist of a
    Unix timestamp column and a "value" column (used as power).
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The offset-aware UTC start of the time interval of interest.
    end_time : datetime.Datetime
        The offset-aware UTC end of the time interval of interest.
    channel : int
        The wavelength of the AIA channel being analyzed (eg. 171 or 211) (Angstrom)
    reg : int
        The arbitrarily assigned number of the region being analyzed
    --------
    Returns:
    --------
    lc_aia : LightCurve
        The LightCurve object created from the AIA data.
    '''
    #Generating the filename given the parameters entered
    fname = f"lc_{start_time.strftime('%Y%m%dT%H%M')}00.0-{end_time.strftime('%Y%m%dT%H%M')}00.0" \
            + f"_{channel}_N21_reg{reg}.csv"
    #Stores the data in a 2D numpy array. First column is time, 2nd column is value
    with open("./DATA/AIA/" + fname, 'r', newline='', encoding='utf8') as file:
        reader = csv.reader(file)
        data = list(reader)
    data = np.transpose(data)
    timestamps = [int(time.split('.')[0]) for time in data[0]]
    power = [int(val.split('.')[0]) for val in data[1]]
    #Creating a LightCurve object from the timestamps and detrended power values
    lc_aia = LightCurve(timestamps, power, "AIA", input_wavelength=channel, input_region=str(reg))
    lc_aia.detrend(21)
    return lc_aia

def gather_aia_fname(fname):
    '''
    Creates a LightCurve object from an AIA csv file created by Reed Masek. This version of the
    method uses the filename directly instead of creating it from input parameters.
    -----------
    Parameters:
    -----------
    fname : str
        The name of the data file from which we wish to make a LightCurve.
    --------
    Returns:
    --------
    lc_aia : LightCurve
        The LightCurve object created from the AIA data.
    '''
    #Stores the data in a 2D numpy array. First column is time, 2nd column is value
    with open("./DATA/AIA/" + fname, 'r', newline='', encoding='utf8') as file:
        reader = csv.reader(file)
        data = list(reader)
    data = np.transpose(data)
    timestamps = [int(time.split('.')[0]) for time in data[0]]
    power = [int(val.split('.')[0]) for val in data[1]]
    #Getting the wavelength (in A) from the file name
    wavelength = int(fname[39:42])
    #Creating a LightCurve object from the timestamps and detrended power values
    lc_aia = LightCurve(timestamps, power, "AIA", input_wavelength=wavelength,
                            input_region=fname[-5])
    lc_aia.detrend(51)
    return lc_aia

def gather_psp_local_files(start_time, end_time, v_str, power_variable):
    '''
    Generates a list of LightCurve objects in the given time constraints from a file already
    downloaded to ./DATA/PSP.
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The offset-aware UTC start of the time interval of interest.
    end_time : datetime.Datetime
        The offset-aware UTC end of the time interval of interest.
    v_str : String
        Specifies whether the file to be read is v01, v02, v03, etc.
    power_variable : String
        Specifies whether to ingest ch0_V1V2, ch1_V3V4, or a combination of the 2.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        The list of LightCurve objects.
    '''
    #Defining lists for data storage
    freq_ranges = ['hfr','lfr']
    lc_list = []
    #Toggling through the HFR and LFR ranges to get all of the available frequencies
    for freq_range in freq_ranges:
        power = []
        timestamps = []
        freqs = []
        file_name = f"psp_fld_l2_rfs_{freq_range}_{start_time.strftime('%Y%m%d')}_v{v_str}.cdf"
        file_location = './DATA/PSP/'
        #Reading through the CDF file to get data
        try:
            cdf_file = cdflib.CDF(file_location + file_name)
        except OSError:
            print("File not available.")
            raise
        #Gathering relevant variables from the CDF file
        v12_data = cdf_file.varget(f'psp_fld_l2_rfs_{freq_range}_auto_averages_ch0_V1V2',
                                    expand=False)
        v34_data = cdf_file.varget(f'psp_fld_l2_rfs_{freq_range}_auto_averages_ch1_V3V4',
                                    expand=False)
        if power_variable == "combination":
            if len(v12_data) != len(v34_data):
                print("ch0_V1V2 and ch1_V3V4 are not compatible. Please choose one or the other")
            else:
                for i, _ in enumerate(v12_data):
                    vector_sum_magnitude = np.sqrt((v12_data[i] ** 2) + (v34_data[i] ** 2))
                    power.append(vector_sum_magnitude)
        elif power_variable == "V1V2":
            for value in v12_data:
                power.append(value)
        elif power_variable == "V3V4":
            for value in v34_data:
                power.append(value)
        else:
            print("Please choose a valid power variable: 'V1V2', 'V3V4', or 'combination'")
            break
        epoch_times = cdf_file.varget('epoch_' + freq_range,expand=False)
        #Timestamps are needed in creating light curves. The PSP data gives times in nanoseconds
        #starting from 2000/01/01 close to noon, whereas the format needed is Unix time
        offset = datetime(2000,1,1,11,58,55,816000, tzinfo=timezone.utc).timestamp() #in seconds
        for epoch in epoch_times:
            timestamps.append((epoch/1e9) + offset)
        if power_variable == "V3V4":
            freqs = cdf_file.varget(f'frequency_{freq_range}_auto_averages_ch1_V3V4')[0]
        elif power_variable == "V1V2":
            freqs = cdf_file.varget(f'frequency_{freq_range}_auto_averages_ch0_V1V2')[0]
        else:
            #Defaults to the V3V4 frequency range for a combination
            freqs = cdf_file.varget(f'frequency_{freq_range}_auto_averages_ch1_V3V4')[0]
        #Converting to datetimes for constriction of data to specified time range
        dts = [datetime.fromtimestamp(int(timestamp), tz=timezone.utc) for timestamp in timestamps]
        constricted_dts = [dt for dt in dts if (dt > start_time and dt < end_time)]
        constricted_timestamps = [dt.timestamp() for dt in constricted_dts]
        #The power must be accordingly restricted to match up with the timestamps correctly.
        constricted_power = np.array([power[idx] for idx, dt in enumerate(dts) \
                                        if (dt > start_time and dt < end_time)])
        for idx, freq in enumerate(freqs):
            if power_variable == "V1V2":
                lc_list.append(LightCurve(constricted_timestamps,constricted_power[:,idx],
                            'PSPV12', input_frequency=freq * 1e-6))
            elif power_variable == "V3V4":
                lc_list.append(LightCurve(constricted_timestamps,constricted_power[:,idx],
                            'PSPV34', input_frequency=freq * 1e-6))
            elif power_variable == "combination":
                lc_list.append(LightCurve(constricted_timestamps,constricted_power[:,idx],
                            'PSPComb', input_frequency=freq * 1e-6))
    # Sorting LightCurves by frequency
    lc_list.sort(key = lc_sort_func)
    return lc_list

def gather_nustar(start_time, end_time, fpm, region, spacing=12):
    '''
    Uses Reed's NuSTAR data csvs to create LightCurve objects from the count rate calculation.
    -----------
    Parameters:
    -----------
    start_time : datetime.Datetime
        The start of the interval of the file for which we wish to create a lightcurve (UTC).
    end_time : datetime.Datetime
        The end of the interval of the file for which we wish to create a lightcurve (UTC).
    fpm : Str ('A' or 'B')
        Specifies the detector on NuSTAR whose data we wish to access.
    region : int
        Specifies which AIA region the NuSTAR data coincides with, if we are investigating both
        concurrently.
    spacing : int (optional)
        Specifies which time bin width (in seconds) of the data in the file
    --------
    Returns:
    --------
    lc_nustar : LightCurve
        The lightcurve object created from the data.
    '''
    start_str = start_time.strftime('%H%M%S')
    end_str = end_time.strftime('%H%M%S')
    #NuSTAR data are saved to a subdirectory specifiying the date
    date = start_time.date().strftime('%Y.%m.%d')
    fname = f'./DATA/NuSTAR/{date}/region{region}_lightcurve_{start_str}_{end_str}_fpm{fpm}_' + \
            f'frames{spacing}s_2.5-12keV.csv'
    #Stores the data in a 2D numpy array. 1st column is time, 4th column is LightCurve value
    with open(fname, 'r', newline='', encoding='utf8') as file:
        reader = csv.reader(file)
        data = list(reader)
    data = np.transpose(data)
    timestamps = []
    power = []
    for i, timestamp in enumerate(data[0][1:-1]):
        #The timestamps in Reed's files include decimals; we truncate these
        ts_int = int(timestamp.split('.')[0])
        #We want the middle of the time bins, so we add half the sample spacing to the front edge
        timestamps.append(ts_int + int(spacing / 2))
        #The first row of the data file is text, so we increment i to get the right power value
        power.append(float(data[3][i+1]))
    #Getting the energy of the lightcurve from the filename convention
    energy = fname.split('_')[-1].split('c')[0][:-1]
    lc_nustar = LightCurve(timestamps, power, f'NuSTAR_fpm{fpm}', input_wavelength=energy,
                            input_region=str(region))
    return lc_nustar
