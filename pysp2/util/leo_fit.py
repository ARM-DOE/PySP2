import numpy as np
from scipy.optimize import curve_fit
import xarray as xr
from .peak_fit import _gaus, _do_fit_records
#from pysp2.util.peak_fit import _do_fit_records, _gaus

def beam_shape(my_binary, beam_position_from='split point', Globals=None):
    """
    
    Calculates the beam shape needed to determine the laser intensity profile
    used for leading-edge-only fits. The beam position is first determined
    from the position of the peak. The beam profile is then calculated by
    positioning all the peaks to that position. The beam shape is then the mean
    of all the normalized profiles arround that position.
    
    Parameters
    ----------
    my_binary : xarray Dataset
            Dataset with gaussian fits. The input data set is retuned by 
            pysp2.util.gaussian_fit().

    Globals: DMTGlobals structure or None
        DMTGlobals structure containing calibration coefficients and detector 
        signal limits.
        
    beam_position_from : str
           'peak maximum' = construct the beam profile arround the maximum peak
           poistion. The maximum peak position is determied from the peak-height
           weighted average peak position.
           'split point' = construct the beam profile around the split position.
           The split position is taken from the split detector. Use this if split
           detector is installed.
           
    Returns
    -------
    my_binary : xarray Dataset
             Dataset with additional statistics and information about the 
             laser beam profile, splitpoint positions relative to the beam 
             profile etc. These are needed for the actual leo_fit() function.
             All variables that are added to the xarray Dataset begin with 
             "leo_". These leo_ variables are available for all particles, hence
             making the leo fit possible for incandesence particles as well.
             
    """
        
    num_base_pts_2_avg = 10 #take from globals
    moving_average_window = 5 #make an argument out of this
    max_amplitude_fraction = 0.033 #make and argument out of this
    #bins = my_binary['columns']
    
    # median_peak_width = my_binary['PkFWHM_ch0'].median().values / 2.35482
    # median_peak_pos = my_binary['FtPos_ch0'].median().values
    # peak_pos_ok = np.logical_and(my_binary['FtPos_ch0'].values >= median_peak_pos - median_peak_width, 
    #                              my_binary['FtPos_ch0'].values <= median_peak_pos + median_peak_width)
    
    
    #change this so that it uses scat_reject_key instead.
    #scatterin signal ok = True
    scatter_high_gain_accepted = np.logical_and.reduce((
        my_binary['PkFWHM_ch0'].values > Globals.ScatMinWidth,
        my_binary['PkFWHM_ch0'].values < Globals.ScatMaxWidth,
        my_binary['PkHt_ch0'].values > Globals.ScatMinPeakHt1,
        my_binary['PkHt_ch0'].values < Globals.ScatMaxPeakHt1,
        my_binary['FtPos_ch0'].values < Globals.ScatMaxPeakPos,
        my_binary['FtPos_ch0'].values > Globals.ScatMinPeakPos))

    scatter_low_gain_accepted = np.logical_and.reduce((
        my_binary['PkFWHM_ch4'].values > Globals.ScatMinWidth,
        my_binary['PkFWHM_ch4'].values < Globals.ScatMaxWidth,
        my_binary['PkHt_ch4'].values > Globals.ScatMinPeakHt2,
        my_binary['PkHt_ch4'].values < Globals.ScatMaxPeakHt2,
        my_binary['FtPos_ch4'].values < Globals.ScatMaxPeakPos,
        my_binary['FtPos_ch4'].values > Globals.ScatMinPeakPos))

        
    #no incandesence signal = True
    no_incand_trigged = np.logical_and(
        my_binary['PkHt_ch1'].values < Globals.IncanMinPeakHt1,
        my_binary['PkHt_ch5'].values < Globals.IncanMinPeakHt2)
    
    #Particles that only scatter light
    only_scattering_high_gain = np.logical_and.reduce((scatter_high_gain_accepted,
                                                       no_incand_trigged))
    only_scattering_low_gain = np.logical_and.reduce((scatter_low_gain_accepted,
                                                       no_incand_trigged))
    
    #iloc = "event_index" and "index" with the scattering only particle events
    iloc_high_gain = np.argwhere(only_scattering_high_gain).flatten()
    iloc_low_gain = np.argwhere(only_scattering_low_gain).flatten()
    
    print('High gain scattering particles for beam analysis :: ',
          np.sum(only_scattering_high_gain))
    print('Low gain scattering particles for beam analysis :: ',
          np.sum(only_scattering_low_gain))
    
    #make an xarray of the purely scattering particles
    my_high_gain_scatterers = my_binary.sel(index = only_scattering_high_gain,
                                            event_index = only_scattering_high_gain)
    my_low_gain_scatterers = my_binary.sel(index = only_scattering_low_gain,
                                            event_index = only_scattering_low_gain)

    
    #numpy array for the normalized beam profiels
    my_high_gain_profiles = np.zeros((my_high_gain_scatterers.sizes['index'],
                                    my_high_gain_scatterers.sizes['columns'])) \
                                    * np.nan
    my_low_gain_profiles = np.zeros((my_low_gain_scatterers.sizes['index'],
                                    my_low_gain_scatterers.sizes['columns'])) \
                                    * np.nan
    
    #my_high_gain_profiles_ = np.zeros_like(my_high_gain_profiles)*np.nan
    
    mean_high_gain_max_peak_pos = int(np.nanmean(my_high_gain_scatterers['PkPos_ch0'].values))
    mean_high_gain_split_pos_float = np.nanmean(my_high_gain_scatterers['PkSplitPos_ch3'].values)
    mean_high_gain_split_pos = np.round(mean_high_gain_split_pos_float).astype(np.int32)
    #mean_high_gain_split_pos = int(np.nanmean(my_high_gain_scatterers['PkSplitPos_ch3'].values))
    mean_low_gain_max_peak_pos = int(np.nanmean(my_low_gain_scatterers['PkPos_ch4'].values))
    mean_low_gain_split_pos_float = np.nanmean(my_low_gain_scatterers['PkSplitPos_ch7'].values)
    mean_low_gain_split_pos = np.round(mean_low_gain_split_pos_float).astype(np.int32)
    #mean_low_gain_split_pos = int(np.nanmean(my_low_gain_scatterers['PkSplitPos_ch7'].values))

    #cross to center
    high_gain_c2c = my_high_gain_scatterers['FtPos_ch0'].values - my_high_gain_scatterers['PkSplitPos_ch3'].values
    low_gain_c2c = my_low_gain_scatterers['FtPos_ch4'].values - my_low_gain_scatterers['PkSplitPos_ch7'].values
    
    #high_gain_split_position = my_high_gain_scatterers['PkSplitPos_ch3'].values
    #low_gain_split_position = my_low_gain_scatterers['PkSplitPos_ch7'].values
    
    #loop through all particle events (THIS CAN BE MADE SMARTER WITH MUCH OF 
    #THE CALCULATIONS BEFORE THE LOOP) --> TBD later
    for i in my_high_gain_scatterers['event_index']:
        data_ch0 = my_high_gain_scatterers['Data_ch0'].sel(event_index=i).values
        #base level
        base_ch0 = np.mean(data_ch0[0:num_base_pts_2_avg])
        #peak height
        peak_height_ch0 = data_ch0.max() - base_ch0
        #max peak position
        peak_pos_ch0 = data_ch0.argmax()
        #split position
        split_pos_ch3 = my_high_gain_scatterers['PkSplitPos_ch3'].sel(event_index=i).values
        if split_pos_ch3 > peak_pos_ch0:
            continue
        #normalize the profile to range [~0,1]
        profile_ch0 = (data_ch0 - base_ch0) / peak_height_ch0
        #insert the profile as it was recorded (no shifting due to PEAK POSITION or PSD POSITION)
        #my_high_gain_profiles_[i,:] = profile
        #distance to the mean beam peak position
        if beam_position_from == 'peak maximum':
            peak_difference_ch0 = mean_high_gain_max_peak_pos - peak_pos_ch0
        elif beam_position_from == 'split point':
            peak_difference_ch0 = mean_high_gain_split_pos - split_pos_ch3
        #insert so that the peak is at the right position (accounts for 
        #particles travelling at different speeds)
        if peak_difference_ch0 > 0:
            my_high_gain_profiles[i, peak_difference_ch0:] = profile_ch0[:len(data_ch0) - 
                                                                peak_difference_ch0]
        elif peak_difference_ch0 < 0:
            my_high_gain_profiles[i, :len(data_ch0)+peak_difference_ch0] = profile_ch0[-peak_difference_ch0:]
        else:
            my_high_gain_profiles[i, :] = profile_ch0

    for i in my_low_gain_scatterers['event_index']:
        data_ch4 = my_low_gain_scatterers['Data_ch4'].sel(event_index=i).values
        #base level
        base_ch4 = np.mean(data_ch4[0:num_base_pts_2_avg])
        #peak height
        peak_height_ch4 = data_ch4.max() - base_ch4
        #max peak position
        peak_pos_ch4 = data_ch4.argmax()
        #split position
        split_pos_ch7 = my_low_gain_scatterers['PkSplitPos_ch7'].sel(event_index=i).values
        if split_pos_ch7 > peak_pos_ch4:
            continue
        #normalize the profile to range [~0,1]
        profile_ch4 = (data_ch4 - base_ch4) / peak_height_ch4
        #insert the profile as it was recorded (no shifting due to PEAK POSITION or PSD POSITION)
        #my_high_gain_profiles_[i,:] = profile
        #distance to the mean beam peak position
        if beam_position_from == 'peak maximum':
            peak_difference_ch4 = mean_low_gain_max_peak_pos - peak_pos_ch4
        elif beam_position_from == 'split point':
            peak_difference_ch4 = mean_low_gain_split_pos - split_pos_ch7
        #insert so that the peak is at the right position (accounts for 
        #particles travelling at different speeds)
        if peak_difference_ch4 > 0:
            my_low_gain_profiles[i, peak_difference_ch4:] = profile_ch4[:len(data_ch4) - 
                                                                peak_difference_ch4]
        elif peak_difference_ch4 < 0:
            my_low_gain_profiles[i, :len(data_ch4)+peak_difference_ch4] = profile_ch4[-peak_difference_ch4:]
        else:
            my_low_gain_profiles[i, :] = profile_ch4


    #MOVING AVERAGE OF THE BEAM PROFILE TO FIND THE DISTANCE BETWEEN THE SPLIT POINT AND THE POINT IN THE 
    #LASER BEAM WHERE PARTICLES CAN START TO EVAPORATE
    
    #moving average of the beam shape with a window of moving_average_window
    moving_high_gain_profile_window = np.lib.stride_tricks.sliding_window_view(my_high_gain_profiles, 
                                                                     moving_average_window, axis=0)
    moving_avg_high_gain_profiles_ = np.nanmedian(moving_high_gain_profile_window,axis=2)
    moving_low_gain_profile_window = np.lib.stride_tricks.sliding_window_view(my_low_gain_profiles, 
                                                                     moving_average_window, axis=0)
    moving_avg_low_gain_profiles_ = np.nanmedian(moving_low_gain_profile_window,axis=2)
    
    moving_avg_high_gain_profiles = np.zeros_like(moving_avg_high_gain_profiles_) * np.nan
    moving_avg_high_gain_split_to_leo_pos = np.zeros(moving_avg_high_gain_profiles_.shape[0]) * np.nan
    moving_avg_high_gain_max_leo_pos = np.zeros(moving_avg_high_gain_profiles_.shape[0]) * np.nan

    moving_avg_low_gain_profiles = np.zeros_like(moving_avg_low_gain_profiles_) * np.nan
    moving_avg_low_gain_split_to_leo_pos = np.zeros(moving_avg_low_gain_profiles_.shape[0]) * np.nan
    moving_avg_low_gain_max_leo_pos = np.zeros(moving_avg_low_gain_profiles_.shape[0]) * np.nan

    moving_avg_high_gain_max_leo_amplitude_factor = np.zeros(moving_avg_high_gain_profiles_.shape[0]) * np.nan
    moving_avg_low_gain_max_leo_amplitude_factor = np.zeros(moving_avg_low_gain_profiles_.shape[0]) * np.nan

    for i in range(moving_avg_high_gain_profiles_.shape[0]):
        i_profile = moving_avg_high_gain_profiles_[i,:]
        i_max = np.nanargmax(i_profile)
        i_range = i_profile[i_max] - np.nanmin(i_profile[:i_max])
        moving_avg_high_gain_profiles[i,:] =  (i_profile - np.nanmin(i_profile[:i_max])) / i_range
        #interpolate here to get the exact position in fraction (not integer) :: which posiiton (float) is the 0.03 cross in
        #and skip if it is the where the baseline is calculated
        above_max_leo_pos = np.ndarray.flatten(np.argwhere(moving_avg_high_gain_profiles[i,:] >= max_amplitude_fraction))
        moving_avg_high_gain_max_leo_pos[i] = above_max_leo_pos[above_max_leo_pos>num_base_pts_2_avg].min()-1
        
        moving_avg_high_gain_split_to_leo_pos[i] = moving_avg_high_gain_max_leo_pos[i] - mean_high_gain_split_pos
        moving_avg_high_gain_max_leo_amplitude_factor[i] = 1./ moving_avg_high_gain_profiles[i, np.round(moving_avg_high_gain_max_leo_pos[i]).astype(int)]

    for i in range(moving_avg_low_gain_profiles_.shape[0]):
        i_profile = moving_avg_low_gain_profiles_[i,:]
        i_max = np.nanargmax(i_profile)
        i_range = i_profile[i_max] - np.nanmin(i_profile[:i_max])
        moving_avg_low_gain_profiles[i,:] =  (i_profile - np.nanmin(i_profile[:i_max])) / i_range
        #interpolate here to get the exact position in fraction (not integer) :: which posiiton (float) is the 0.03 cross in
        #and skip if it is the where the baseline is calculated
        above_max_leo_pos = np.ndarray.flatten(np.argwhere(moving_avg_low_gain_profiles[i,:] >= max_amplitude_fraction))
        moving_avg_low_gain_max_leo_pos[i] = above_max_leo_pos[above_max_leo_pos>num_base_pts_2_avg].min()-1
        
        moving_avg_low_gain_split_to_leo_pos[i] = moving_avg_low_gain_max_leo_pos[i] - mean_low_gain_split_pos
        moving_avg_low_gain_max_leo_amplitude_factor[i] = 1./ moving_avg_low_gain_profiles[i, np.round(moving_avg_low_gain_max_leo_pos[i]).astype(int)]

    #cleaning up
    moving_avg_high_gain_max_leo_pos = np.where(moving_avg_high_gain_max_leo_pos < num_base_pts_2_avg, 
                                                np.nan, moving_avg_high_gain_max_leo_pos)    
    moving_avg_high_gain_split_to_leo_pos = np.where(moving_avg_high_gain_split_to_leo_pos < -30. ,
                                                   np.nan, moving_avg_high_gain_split_to_leo_pos)
    moving_avg_high_gain_max_leo_amplitude_factor = np.where(moving_avg_high_gain_max_leo_pos < num_base_pts_2_avg,
                                                         np.nan, moving_avg_high_gain_max_leo_amplitude_factor)

    moving_avg_low_gain_max_leo_pos = np.where(moving_avg_low_gain_max_leo_pos < num_base_pts_2_avg, 
                                                np.nan, moving_avg_low_gain_max_leo_pos)    
    moving_avg_low_gain_split_to_leo_pos = np.where(moving_avg_low_gain_split_to_leo_pos < -30. ,
                                                   np.nan, moving_avg_low_gain_split_to_leo_pos)
    moving_avg_low_gain_max_leo_amplitude_factor = np.where(moving_avg_low_gain_max_leo_pos < num_base_pts_2_avg,
                                                         np.nan, moving_avg_low_gain_max_leo_amplitude_factor)

    
    #moving average of beam width
    moving_high_gain_beam_width = np.lib.stride_tricks.sliding_window_view(my_high_gain_scatterers['PkFWHM_ch0'].values, 
                                                                     moving_average_window, axis=0)
    moving_low_gain_beam_width = np.lib.stride_tricks.sliding_window_view(my_low_gain_scatterers['PkFWHM_ch4'].values, 
                                                                     moving_average_window, axis=0)
    moving_median_high_gain_beam_width = np.nanmedian(moving_high_gain_beam_width,axis=1)
    moving_median_low_gain_beam_width = np.nanmedian(moving_low_gain_beam_width,axis=1)
    
    #Moving leo_Base
    #moving_high_gain_base = np.lib.stride_tricks.sliding_window_view(my_high_gain_scatterers['Base_ch0'].values, 
    #                                                                 moving_average_window, axis=0)
    #moving_low_gain_base = np.lib.stride_tricks.sliding_window_view(my_low_gain_scatterers['Base_ch4'].values, 
    #                                                                 moving_average_window, axis=0)
    #moving_median_high_gain_base = np.nanpercentile(moving_high_gain_base, 10,axis=1)
    #moving_median_low_gain_base = np.nanpercentile(moving_low_gain_base, 10,axis=1)
    
    #Moving cross to centre (c2c)
    moving_high_gain_c2c = np.lib.stride_tricks.sliding_window_view(high_gain_c2c, 
                                                                     moving_average_window, axis=0)
    moving_low_gain_c2c = np.lib.stride_tricks.sliding_window_view(low_gain_c2c, 
                                                                     moving_average_window, axis=0)

    moving_median_high_gain_c2c = np.nanmedian(moving_high_gain_c2c,axis=1)
    moving_median_low_gain_c2c = np.nanmedian(moving_low_gain_c2c,axis=1)
    
    #JB LATER
    leo_AmpFactor_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan 
    leo_AmpFactor_ch0[iloc_high_gain[:-moving_average_window+1]] = moving_avg_high_gain_max_leo_amplitude_factor
    leo_AmpFactor_ch4 = np.zeros(scatter_low_gain_accepted.shape)*np.nan 
    leo_AmpFactor_ch4[iloc_low_gain[:-moving_average_window+1]] = moving_avg_low_gain_max_leo_amplitude_factor

    
    leo_PkFWHM_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_PkFWHM_ch0[iloc_high_gain[:-moving_average_window+1]] = moving_median_high_gain_beam_width
    leo_PkFWHM_ch4 = np.zeros(scatter_low_gain_accepted.shape)*np.nan
    leo_PkFWHM_ch4[iloc_low_gain[:-moving_average_window+1]] = moving_median_low_gain_beam_width

    #TAKE FROM ACTUAL PARTICLE TRACE WITH INCANDESENCE (NOT NEEDED IN ACTUAL LEO FIT?)
    #leo_Base_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    #leo_Base_ch0[iloc_high_gain[:-moving_average_window+1]] = moving_median_high_gain_base
    #leo_Base_ch4 = np.zeros(scatter_low_gain_accepted.shape)*np.nan
    #leo_Base_ch4[iloc_low_gain[:-moving_average_window+1]] = moving_median_low_gain_base
    
    #JB OK
    leo_c2c_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_c2c_ch0[iloc_high_gain[:-moving_average_window+1]] = moving_median_high_gain_c2c
    leo_c2c_ch4 = np.zeros(scatter_low_gain_accepted.shape)*np.nan
    leo_c2c_ch4[iloc_low_gain[:-moving_average_window+1]] = moving_median_low_gain_c2c
    
    #JB OK
    leo_split_to_leo_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_split_to_leo_ch0[iloc_high_gain[:-moving_average_window+1]] = moving_avg_high_gain_split_to_leo_pos
    leo_split_to_leo_ch4 = np.zeros(scatter_low_gain_accepted.shape)*np.nan
    leo_split_to_leo_ch4[iloc_low_gain[:-moving_average_window+1]] = moving_avg_low_gain_split_to_leo_pos
    
    """
    WHAT IS NEEDED:
    * GAUSSIAN WIDTH FOR ALL PARTICLES - FROM PROFILES (DONE)
    * PEAK POSITION FOR ALL PARTICLES - FROM C2C (DONE)
    * LAST BIN TO USE FOR LEO FIT (from split detector)  (DONE)
    * do also for LG.
    * TAKE BASELINE SCAT FROM ACTUAL TRACERS IN LEO_FIT
    """
    
    output_ds = my_binary.copy()
    
    output_ds['leo_AmpFactor_ch0'] = (('event_index'), leo_AmpFactor_ch0)
    output_ds['leo_AmpFactor_ch0'] = output_ds['leo_AmpFactor_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_AmpFactor_ch4'] = (('event_index'), leo_AmpFactor_ch4)
    output_ds['leo_AmpFactor_ch4'] = output_ds['leo_AmpFactor_ch4'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkFWHM_ch0'] = (('event_index'), leo_PkFWHM_ch0)
    output_ds['leo_PkFWHM_ch4'] = (('event_index'), leo_PkFWHM_ch4)

    #distance from cross-to-centre (split point to laser maximum intensity). 
    #This comes from scattering only partilces
    output_ds['leo_PkPos_ch0'] = (('event_index'), leo_c2c_ch0)
    output_ds['leo_PkPos_ch4'] = (('event_index'), leo_c2c_ch4)
    #interpolate to all particles (including, and most importantly, rBC containing particles)
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkPos_ch4'] = output_ds['leo_PkPos_ch4'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    #add split position to cross-to-centre distance to get location of the (estimated) 
    #location of the peak maximum. This is needed for the LEO-fit.
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'] + my_binary['PkSplitPos_ch3']
    output_ds['leo_PkPos_ch4'] = output_ds['leo_PkPos_ch4'] + my_binary['PkSplitPos_ch7']
    
    #First add t_alpha_to_split data
    output_ds['leo_EndPos_ch0'] = (('event_index'), leo_split_to_leo_ch0)
    output_ds['leo_EndPos_ch4'] = (('event_index'), leo_split_to_leo_ch4)
    #interpolate to all particles
    output_ds['leo_EndPos_ch0'] = output_ds['leo_EndPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_EndPos_ch4'] = output_ds['leo_EndPos_ch4'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")

    #calculate the position at which the leo fits should end based on the split position 
    #(of all particles)
    output_ds['leo_EndPos_ch0'] = output_ds['leo_EndPos_ch0'] + output_ds['PkSplitPos_ch3']
    output_ds['leo_EndPos_ch4'] = output_ds['leo_EndPos_ch4'] + output_ds['PkSplitPos_ch7']
    #cleaning up
    output_ds['leo_EndPos_ch0'].values = np.where(output_ds['leo_EndPos_ch0'].values < num_base_pts_2_avg,
                                                  np.nan, output_ds['leo_EndPos_ch0'].values)
    output_ds['leo_EndPos_ch4'].values = np.where(output_ds['leo_EndPos_ch4'].values < num_base_pts_2_avg,
                                                  np.nan, output_ds['leo_EndPos_ch4'].values)
        
    output_ds['leo_PkFWHM_ch0'] = output_ds['leo_PkFWHM_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkFWHM_ch4'] = output_ds['leo_PkFWHM_ch4'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkPos_ch4'] = output_ds['leo_PkPos_ch4'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
        
    return output_ds

def leo_fit(my_binary,Globals=None):
    bins = my_binary['columns'].astype('float').values
    #number of points at the beginning to use for base line average
    num_base_pts_2_avg = 10
    #max_amplitude_fraction = 0.03
    
    #get the information about the gaussian fits
    pos_ch0 = my_binary['leo_PkPos_ch0'].values
    pos_ch4 = my_binary['leo_PkPos_ch4'].values

    #amplitude = my_binary['PkHt_ch0'].values
    width_ch0 = my_binary['leo_PkFWHM_ch0'].values / 2.35482 #2*sqrt(2*log(2))
    width_ch4 = my_binary['leo_PkFWHM_ch4'].values / 2.35482 #2*sqrt(2*log(2))
    data_ch0 = my_binary['Data_ch0'].values
    data_ch4 = my_binary['Data_ch4'].values
    
    #mean of the first num_base_pts_2_avg points
    #leo_base_ch0 = np.mean(data_ch0[:, 0:num_base_pts_2_avg], axis=1)
    #leo_base_ch0 = my_binary['leo_Base_ch0'].values
    data_ch0_sorted = np.sort(data_ch0[:, 0:num_base_pts_2_avg], axis=1)
    data_ch4_sorted = np.sort(data_ch4[:, 0:num_base_pts_2_avg], axis=1)
    leo_base_ch0 = np.min(data_ch0_sorted[:, 0:int(num_base_pts_2_avg)], axis=1)
    leo_base_ch4 = np.min(data_ch4_sorted[:, 0:int(num_base_pts_2_avg)], axis=1)
    
    leo_fit_max_pos_ch0 = my_binary['leo_EndPos_ch0'].astype(int).values
    leo_fit_max_pos_ch4 = my_binary['leo_EndPos_ch4'].astype(int).values
    leo_AmpFactor_ch0 = my_binary['leo_AmpFactor_ch0'].values
    leo_AmpFactor_ch4 = my_binary['leo_AmpFactor_ch4'].values
    leo_PkHt_ch0 = np.zeros_like(my_binary['PkHt_ch0'].values)*np.nan
    leo_PkHt_ch4 = np.zeros_like(my_binary['PkHt_ch4'].values)*np.nan

    #High gain
    for i in range(my_binary.sizes['event_index']):
        #max_value = data_ch0[i,:].max() - data_ch0[i,:].min()
        #bins_ = bins[num_base_pts_2_avg:leo_fit_max_pos[i]]
        bins_ = bins[leo_fit_max_pos_ch0[i]-3:leo_fit_max_pos_ch0[i]]
        if len(bins_) < 2:
            leo_PkHt_ch0[i] = np.nan
            continue
        #signals
        data_ch0_ = data_ch0[i, num_base_pts_2_avg:leo_fit_max_pos_ch0[i]]
        data_ch0_ = data_ch0[i, leo_fit_max_pos_ch0[i]-3:leo_fit_max_pos_ch0[i]]
        leo_coeff, var_matrix = curve_fit(
            lambda x, a: _gaus(x, a, pos_ch0[i], width_ch0[i], leo_base_ch0[i]), 
            bins_[:], data_ch0_[:], p0=[data_ch0[i,:].max()], maxfev=100, 
            ftol=1e-5, method='lm' ) #, bounds=(0, 1e6)) #, method='lm'
        leo_PkHt_ch0[i] = leo_coeff[0]
        leo_PkHt_ch0[i] = (data_ch0[i, leo_fit_max_pos_ch0[i]] - leo_base_ch0[i]) * leo_AmpFactor_ch0[i]
    
    #Low gain
    for i in range(my_binary.sizes['event_index']):
        #max_value = data_ch0[i,:].max() - data_ch0[i,:].min()
        #bins_ = bins[num_base_pts_2_avg:leo_fit_max_pos[i]]
        bins_ = bins[leo_fit_max_pos_ch4[i]-3:leo_fit_max_pos_ch4[i]]
        if len(bins_) < 2:
            leo_PkHt_ch4[i] = np.nan
            continue
        #signals
        data_ch4_ = data_ch4[i, num_base_pts_2_avg:leo_fit_max_pos_ch4[i]]
        data_ch4_ = data_ch4[i, leo_fit_max_pos_ch4[i]-3:leo_fit_max_pos_ch4[i]]
        leo_coeff, var_matrix = curve_fit(
            lambda x, a: _gaus(x, a, pos_ch4[i], width_ch4[i], leo_base_ch4[i]), 
            bins_[:], data_ch4_[:], p0=[data_ch4[i,:].max()], maxfev=100, 
            ftol=1e-5, method='lm' ) #, bounds=(0, 1e6)) #, method='lm'
        leo_PkHt_ch4[i] = leo_coeff[0]
        leo_PkHt_ch4[i] = (data_ch4[i, leo_fit_max_pos_ch4[i]] - leo_base_ch4[i]) * leo_AmpFactor_ch4[i]
    
    #Only positive values, negative values are set to nan.
    leo_PkHt_ch0 = np.where(leo_PkHt_ch0>0, leo_PkHt_ch0, np.nan)
    leo_PkHt_ch4 = np.where(leo_PkHt_ch4>0, leo_PkHt_ch4, np.nan)

    output_ds = my_binary.copy()
    output_ds['leo_FtAmp_ch0'] = (('index'), leo_PkHt_ch0)
    output_ds['leo_FtAmp_ch4'] = (('index'), leo_PkHt_ch4)
    output_ds['leo_Base_ch0'] = (('index'), leo_base_ch0)
    output_ds['leo_Base_ch4'] = (('index'), leo_base_ch4)

    return output_ds

