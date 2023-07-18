'''
This documents defines the methods used the control the flow of information given to the main
function, based on the source we wish to analyze.
'''

#Imports
import gather_functions as gf
from utility import make_dts_from_input, get_lc_at_freq, str_to_dt
from lightcurve_class import split_lc_list

def prompt_time():
    '''
    Asks questions about the start and end times of intervals for creating LightCurves. Capable of
    creating multiple intervals on the same date and handling roll-over from one day to the next.
    Also, if only one interval is give, the user is prompted with the option to split that interval
    into 2 subintervals about a time in the middle.
    --------
    Returns:
    --------
    intervals : list of [datetime.Datetime, datetime.Datetime] lists
        The datetime objects describing the start and end times of each interval given.
    split_info : list
        Contains a Boolean toggling whether or not splitting is to be performed, and the time about
        which the splitting is to happen (if necessary).
    '''
    date = input('Enter date as YYYY.mm.dd\n')
    intervals = []
    move_on = ''
    while move_on != 'done':
        start = input('Enter start time (UTC) as HH:MM or HH:MM:SS\n')
        end = input('Enter end time (UTC) as HH:MM or HH:MM:SS\n')
        intervals.append(make_dts_from_input(date, start, end))
        move_on = input("Press Enter to give another time interval, or enter 'done' to move on\n")
    #Offering the option to split an interval into subintervals, if only one is given
    if len(intervals) == 1:
        bool_split = bool(int(input('Would you like to split this interval into 2 ' + \
                                        'subintervals? Enter 0 for no, 1 for yes\n')))
        split_info = [bool_split]
        if bool_split:
            split_str = input('Enter the time about which to split the interval as HH:MM or ' + \
                                'HH:MM:SS\n')
            split_dt = str_to_dt(f'{date} {split_str}')
            if split_dt <= intervals[0][0] or split_dt >= intervals[0][1]:
                print('Enter a time that is between the inteval start and end times. Exitting.')
                return
            split_info.append(split_dt)
    else:
        split_info = [False]
    return intervals, split_info

def prompt_psp():
    '''
    Asks the flow of questions necessary to run routines on PSP data, and acquires the requested
    LightCurves.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        A list containing the LightCurves specified by the user's answers.
    '''
    all_lcs = []
    intervals, split_info = prompt_time()
    variable = input("Which PSP variable would you like? Enter 'V1V2', 'V3V4', or 'combination'\n")
    method = int(input('Enter 0 to gather data from cdaweb, or 1 to gather from a local file\n'))
    #The following questions depend on which method is to be performed
    if method == 0:
        auto = input("Would you like averages or peaks data? Enter 'averages' or 'peaks'\n")
        bool_shift = bool(int(input('Would you like to adjust the times to account for the ' + \
                                        'Earth-PSP distance? Enter 0 for no, 1 for yes\n')))
        for interval in intervals:
            all_lcs.append(gf.gather_psp(interval[0], interval[1], variable, auto, bool_shift))
    elif method == 1:
        v_str = input("What number appears after v in the file name? Enter as '01', '02', etc.\n")
        for interval in intervals:
            try:
                all_lcs.append(gf.gather_psp_local_files(interval[0], interval[1], v_str, variable))
            except FileNotFoundError:
                print(f"No v{v_str} PSP data file found for the " + \
                                f"interval beginning at {interval[0].strftime('%Y.%m.%d %H:%M')}")
                continue
    else:
        print('Enter a 0 or 1. Exitting.')
        return
    #Picking out LightCurves at the specified frequencies
    freqs = []
    input_freq = ''
    while input_freq != 'done':
        input_freq = input("Enter frequencies (in MHz) one at a time, or 'done' to move on\n")
        if input_freq != 'done':
            freqs.append(float(input_freq))
    lc_list = []
    #LightCurves at the given frequencies are returned for all time intervals specified earlier
    for lc_group in all_lcs:
        for freq in freqs:
            lc_list.append(get_lc_at_freq(lc_group, freq))
    #Splitting the LightCurves into subintervals if requested
    if split_info[0]:
        lc_list_split = split_lc_list(lc_list, split_info[1])
        lc_list = lc_list_split
    return lc_list

def prompt_aia():
    '''
    Asks the flow of questions necessary to run routines on AIA data, and acquires the requested
    LightCurves.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        A list containing the LightCurves specified by the user's answers.
    '''
    bool_fname = bool(int(input('Would you like to input filenames to skip questions? Enter 0 ' + \
                                'for no, 1 for yes\n')))
    if bool_fname:
        fname_list = []
        fname = ''
        while fname != 'done':
            fname = input("Enter a filename (path not necessary), or 'done' to move on\n")
            if fname != 'done':
                fname_list.append(fname)
        lc_list = []
        for file in fname_list:
            lc_list.append(gf.gather_aia_fname(file))
    else:
        intervals, split_info = prompt_time()
        #Getting the channel(s) to be analyzed
        channels = []
        channel = ''
        while channel != 'done':
            channel = input("Enter wavelengths in A one at a time, or 'done' to move on\n")
            if channel != 'done':
                channels.append(int(channel))
        #Getting the regions(s) to be analyzed
        regions = []
        region = ''
        while region != 'done':
            region = input("Enter region numbers one at a time, or 'done' to move on\n")
            if region != 'done':
                regions.append(int(region))
        #Making LightCurves for all regions at all channels for each interval
        lc_list = []
        for interval in intervals:
            for channel in channels:
                for region in regions:
                    try:
                        lc_list.append(gf.gather_aia(interval[0], interval[1], channel, region))
                    except FileNotFoundError:
                        print(f"No AIA data file found for the {channel} A, region {region} " + \
                                f"interval beginning at {interval[0].strftime('%Y.%m.%d %H:%M')}")
                        continue
        #Splitting the LightCurves into subintervals if requested
        if split_info[0]:
            lc_list_split = split_lc_list(lc_list, split_info[1])
            lc_list = lc_list_split
        return lc_list

def prompt_nustar():
    '''
    Asks the flow of questions necessary to run routines on NuSTAR data, and acquires the requested
    LightCurves.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        A list containing the LightCurves specified by the user's answers.
    '''
    intervals, split_info = prompt_time()
    fpm = input("Which detector? Enter 'A', 'B', or 'both'\n")
    if fpm == 'both':
        fpms = ['A', 'B']
    elif fpm == 'A':
        fpms = ['A']
    elif fpm == 'B':
        fpms = ['B']
    else:
        print("Please answer 'A', 'B', or 'both'. Exitting")
        return
    #Getting the regions(s) to be analyzed
    regions = []
    region = ''
    while region != 'done':
        region = input("Enter region number, or 'done' to move on\n")
        if region != 'done':
            regions.append(int(region))
    lc_list = []
    for interval in intervals:
        for fpm_val in fpms:
            for region in regions:
                try:
                    new_lc = gf.gather_nustar(interval[0], interval[1], fpm_val, region)
                    new_lc.detrend(51)
                    lc_list.append(new_lc)
                except FileNotFoundError:
                    print(f"No NuSTAR data file found for the FPM {fpm_val}, region {region} " + \
                                f"interval beginning at {interval[0].strftime('%Y.%m.%d %H:%M')}")
                    continue
    #Splitting the LightCurves into subintervals if requested
    if split_info[0]:
        lc_list_split = split_lc_list(lc_list, split_info[1])
        lc_list = lc_list_split
    return lc_list

def prompt_wind():
    '''
    Asks the flow of questions necessary to run routines on WIND data, and acquires the requested
    LightCurves.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        A list containing the LightCurves specified by the user's answers.
    '''
    intervals, split_info = prompt_time()
    #Picking out LightCurves at the desired frequencies
    freqs = []
    input_freq = ''
    while input_freq != 'done':
        input_freq = input("Enter frequencies (in MHz) one at a time, or 'done' to move on\n")
        if input_freq != 'done':
            freqs.append(float(input_freq))
    lc_list = []
    for interval in intervals:
        lc_group = gf.gather_wind(interval[0], interval[1])
        for freq in freqs:
            lc_list.append(get_lc_at_freq(lc_group, freq))
    #Splitting the LightCurves into subintervals if requested
    if split_info[0]:
        lc_list_split = split_lc_list(lc_list, split_info[1])
        lc_list = lc_list_split
    return lc_list

def prompt_stereo():
    '''
    Asks the flow of questions necessary to run routines on STEREO data, and acquires the requested
    LightCurves.
    --------
    Returns:
    --------
    lc_list : list of LightCurve
        A list containing the LightCurves specified by the user's answers.
    '''
    intervals, split_info = prompt_time()
    all_lcs = []
    #Checking which STEREO satellite to obtain data from
    ab_sel = input("Which STEREO satellite would you like data from? Enter 'A' or 'B'\n")
    if ab_sel == 'A':
        for interval in intervals:
            all_lcs.append(gf.gather_stereo_a(interval[0], interval[1]))
    elif ab_sel == 'B':
        for interval in intervals:
            all_lcs.append(gf.gather_stereo_b(interval[0], interval[1]))
    else:
        print("Enter 'A' or 'B'. Exitting.")
        return
    #Picking out LightCurves at the desired frequencies
    freqs = []
    input_freq = ''
    while input_freq != 'done':
        input_freq = input("Enter frequencies (in MHz) one at a time, or 'done' to move on\n")
        if input_freq != 'done':
            freqs.append(float(input_freq))
    lc_list = []
    for lc_group in all_lcs:
        for freq in freqs:
            lc_list.append(get_lc_at_freq(lc_group, freq))
    #Splitting the LightCurves into subintervals if requested
    if split_info[0]:
        lc_list_split = split_lc_list(lc_list, split_info[1])
        lc_list = lc_list_split
    return lc_list

def prompt_methods():
    '''
    Asks the questions necessary to determine which methods to run, and how to run them.
    --------
    Returns:
    --------
    method_params : 2D array
        Array containing the instructions to run or not run a method, as well as the parameters
        for the run if necessary. 1st entry is peak finder, 2nd is wavelet, 3rd is FFT.
    bool_show : Boolean
        Toggles whether or not to show results.
    bool_save : Boolean
        Toggles whether or not to save results.
    '''
    method_params = [[], [], []]
    #Peak finder
    bool_pf = bool(int(input('Do you want to run the peak finder? Enter 0 for no, 1 for yes\n')))
    method_params[0].append(bool_pf)
    #Wavelet
    bool_wavelet = bool(int(input('Do you want to run the wavelet? Enter 0 for no, 1 for yes\n')))
    method_params[1].append(bool_wavelet)
    #Asking follow up questions and saving the answers if necessary
    if bool_wavelet:
        bool_flatten = bool(int(input('Do you want to perform flattening as a preprocessing ' + \
                                        'step? Enter 0 for no, 1 for yes\n')))
        method_params[1].append(bool_flatten)
    #FFT
    bool_fft = bool(int(input('Do you want to run the FFT? Enter 0 for no, 1 for yes\n')))
    method_params[2].append(bool_fft)
    #Asking about showing and saving results (goes for all methods)
    bool_show = bool(int(input('Do you want to display results? Enter 0 for no, 1 for yes\n')))
    bool_save = bool(int(input('Do you want to save results? Enter 0 for no, 1 for yes\n')))
    return method_params, bool_show, bool_save
