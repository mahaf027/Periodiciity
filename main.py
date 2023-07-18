'''
This file defines and runs the function that performs periodicity analysis routines following the
parameters input via prompts given to the user.
'''

#Imports
import prompting

def main():
    '''
    Prompts the user for the necessary information, and runs periodicity analysis routines
    accordingly.
    '''
    done_str = ''
    while done_str != 'done':
        #Getting the necessary information and creating LightCurves
        source = input('Enter the instrument whose data you wish to analyze\n')
        if source.upper() == 'PSP':
            lc_list = prompting.prompt_psp()
        elif source.upper() == 'AIA':
            lc_list = prompting.prompt_aia()
        elif source.upper() == 'NUSTAR':
            lc_list = prompting.prompt_nustar()
        elif source.upper() == 'WIND':
            lc_list = prompting.prompt_wind()
        elif source.upper() == 'STEREO':
            lc_list = prompting.prompt_wind()
        else:
            print("Enter a valid instrument name. Accepted values: 'psp', 'aia', 'nustar', " + \
                    "'wind', and 'stereo'")
            return
        #Asking which periodicity methods to run and the necessary follow-up questions
        method_params, bool_show, bool_save = prompting.prompt_methods()
        #Running the methods as specified on the LightCurves in lc_list
        for lightcurve in lc_list:
            #This is a parameter for the FFT that is changed if the wavelet is run
            period_ax_range = 0
            #Peak finder
            if method_params[0][0]:
                _, _, _, _ = lightcurve.peak_finder(180, show_results=bool_show,
                                                    save_results=bool_save)
                lightcurve.plot_time_series(show_peaks=True, input_dist=180, show_plot=bool_show,
                                            save_plot=bool_save)
            #FFT
            if method_params[1][0]:
                period_ax_range = lightcurve.wavelet(bool_show, bool_save, 'Morlet',
                                                        method_params[1][1])
            #Wavelet
            if method_params[2][0]:
                lightcurve.fft_periodicity(bool_show, bool_save, period_ax_range)
        #Giving the option to do the whole process again
        done_str = input("Press Enter to repeat this process, or enter 'done' to exit\n")

if __name__ == '__main__':
    main()
