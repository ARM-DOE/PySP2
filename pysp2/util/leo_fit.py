import numpy as np
from scipy.optimize import curve_fit
import xarray as xr
#from .peak_fit import _gaus, _do_fit_records
from pysp2.util.peak_fit import _do_fit_records

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
           The split position is taken from the split detector. Not working yet.
           
    Returns
    -------
    coeff : numpy array with gaussian fit [amplitude, peakpos, width, base] 
            of the beam shape profile.
    beam_profile : numpy array with the beam profile calculated from the mean
                   of all the profiles. The array has as many data points as 
                   the input data.

    """
        
    num_base_pts_2_avg = 10 #take from globals
    moving_average_window = 5 #make an argument out of this
    max_amplitude_fraction = 0.033 #make and argument out of this
    bins = my_binary['columns']
    
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
        
    #no incandesence signal = True
    no_incand_trigged = np.logical_and(
        my_binary['PkHt_ch1'].values < Globals.IncanMinPeakHt1,
        my_binary['PkHt_ch5'].values < Globals.IncanMinPeakHt2)
    
    #Particles that only scatter light
    only_scattering_high_gain = np.logical_and.reduce((scatter_high_gain_accepted,
                                                       no_incand_trigged))
    #iloc = "event_index" and "index" with the scattering only particle events
    iloc = np.argwhere(only_scattering_high_gain).flatten()
    
    print('High gain scattering particles for beam analysis :: ',
          np.sum(only_scattering_high_gain))
    
    #make an xarray of the purely scattering particles
    my_high_gain_scatterers = my_binary.sel(index = only_scattering_high_gain,
                                            event_index = only_scattering_high_gain)
    
    #numpy array for the normalized beam profiels
    my_high_gain_profiles = np.zeros((my_high_gain_scatterers.sizes['index'],
                                    my_high_gain_scatterers.sizes['columns'])) \
                                    * np.nan
    
    my_high_gain_profiles_ = np.zeros_like(my_high_gain_profiles)*np.nan
    
    mean_high_gain_max_peak_pos = int(np.nanmean(my_high_gain_scatterers['PkPos_ch0'].values))
    mean_high_gain_split_pos_float = np.nanmean(my_high_gain_scatterers['PkSplitPos_ch3'].values)
    mean_high_gain_split_pos = np.round(mean_high_gain_split_pos_float).astype(np.int32)
    #mean_high_gain_split_pos = int(np.nanmean(my_high_gain_scatterers['PkSplitPos_ch3'].values))
        
    #cross to center
    high_gain_c2c = my_high_gain_scatterers['FtPos_ch0'].values - my_high_gain_scatterers['PkSplitPos_ch3'].values
        
    high_gain_split_position = my_high_gain_scatterers['PkSplitPos_ch3'].values
    
    #loop through all particle events (THIS CAN BE MADE SMARTER WITH MUCH OF 
    #THE CALCULATIONS BEFORE THE LOOP) --> TBD later
    for i in my_high_gain_scatterers['event_index']:
        data = my_high_gain_scatterers['Data_ch0'].sel(event_index=i).values
        #base level
        base = np.mean(data[0:num_base_pts_2_avg])
        #peak height
        peak_height = data.max() - base
        #max peak position
        peak_pos = data.argmax()
        #split position
        split_pos = my_high_gain_scatterers['PkSplitPos_ch3'].sel(event_index=i).values
        #normalize the profile to range [~0,1]
        profile = (data - base) / peak_height
        #insert the profile as it was recorded (no shifting due to PEAK POSITION or PSD POSITION)
        #my_high_gain_profiles_[i,:] = profile
        #distance to the mean beam peak position
        if beam_position_from == 'peak maximum':
            peak_difference = mean_high_gain_max_peak_pos - peak_pos
        elif beam_position_from == 'split point':
            peak_difference = mean_high_gain_split_pos - split_pos
        #insert so that the peak is at the right position (accounts for 
        #particles travelling at different speeds)
        if peak_difference > 0:
            my_high_gain_profiles[i, peak_difference:] = profile[:len(data) - 
                                                                peak_difference]
        elif peak_difference < 0:
            my_high_gain_profiles[i, :len(data)+peak_difference] = profile[-peak_difference:]
        else:
            my_high_gain_profiles[i, :] = profile

    #MOVING AVERAGE OF THE BEAM PROFILE TO FIND THE DISTANCE BETWEEN THE SPLIT POINT AND THE POINT IN THE 
    #LASER BEAM WHERE PARTICLES CAN START TO EVAPORATE
    
    #moving average of the beam shape with a window of moving_average_window
    moving_high_gain_profile_window = np.lib.stride_tricks.sliding_window_view(my_high_gain_profiles, 
                                                                     moving_average_window, axis=0)
#    moving_avg_high_gain_profiles_ = np.nanmean(moving_high_gain_profile_window,axis=2)
    moving_avg_high_gain_profiles_ = np.nanmedian(moving_high_gain_profile_window,axis=2)
    
    # #
    # moving_high_gain_FtPos = np.lib.stride_tricks.sliding_window_view(my_high_gain_scatterers['FtPos_ch0'].values, 
    #                                                                  moving_average_window, axis=0)
    # moving_high_gain_FtPos_ = np.nanmedian(moving_high_gain_FtPos,axis=1)

    
    moving_avg_high_gain_profiles = np.zeros_like(moving_avg_high_gain_profiles_) * np.nan
    moving_avg_high_gain_split_to_leo_pos = np.zeros(moving_avg_high_gain_profiles_.shape[0]) * np.nan
    moving_avg_high_gain_max_leo_pos = np.zeros(moving_avg_high_gain_profiles_.shape[0]) * np.nan

    moving_avg_high_gain_max_leo_amplitude_factor = np.zeros(moving_avg_high_gain_profiles_.shape[0]) * np.nan

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

    #cleaning up
    moving_avg_high_gain_max_leo_pos = np.where(moving_avg_high_gain_max_leo_pos < num_base_pts_2_avg, 
                                                np.nan, moving_avg_high_gain_max_leo_pos)    
    moving_avg_high_gain_split_to_leo_pos = np.where(moving_avg_high_gain_split_to_leo_pos < -30. ,
                                                   np.nan, moving_avg_high_gain_split_to_leo_pos)
    moving_avg_high_gain_max_leo_amplitude_factor = np.where(moving_avg_high_gain_max_leo_pos < num_base_pts_2_avg,
                                                         np.nan, moving_avg_high_gain_max_leo_amplitude_factor)

    # leo_pos_unique = np.unique(moving_avg_high_gain_split_to_leo_pos)
    # mean_leo_pos_unique = np.zeros_like(leo_pos_unique)
    # mean_leo_pos_unique_num = np.zeros_like(leo_pos_unique)
    # for i,pos in enumerate(leo_pos_unique):
    #     bl = moving_avg_high_gain_split_to_leo_pos == pos
    #     mean_leo_pos_unique[i] = np.nanmedian(moving_avg_high_gain_max_leo_amplitude_factor[bl])
    #     mean_leo_pos_unique_num[i] = sum(bl)

    #moving average of beam width
    moving_high_gain_beam_width = np.lib.stride_tricks.sliding_window_view(my_high_gain_scatterers['PkFWHM_ch0'].values, 
                                                                     moving_average_window, axis=0)
    #moving_avg_high_gain_beam_width = np.nanmedian(moving_high_gain_beam_width,axis=1)
    moving_median_high_gain_beam_width = np.nanmedian(moving_high_gain_beam_width,axis=1)
    
    #Moving leo_Base_ch0
    moving_high_gain_base = np.lib.stride_tricks.sliding_window_view(my_high_gain_scatterers['Base_ch0'].values, 
                                                                     moving_average_window, axis=0)
    moving_median_high_gain_base = np.nanpercentile(moving_high_gain_base, 10,axis=1)
    
    #JB OK
    #Moving c2c high gain
    moving_high_gain_c2c = np.lib.stride_tricks.sliding_window_view(high_gain_c2c, 
                                                                     moving_average_window, axis=0)
    moving_median_high_gain_c2c = np.nanmedian(moving_high_gain_c2c,axis=1)
    
    #JB LATER
    #leo_FtMaxPosAmpFactor_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan 
    #leo_FtMaxPosAmpFactor_ch0[iloc[:-moving_average_window+1]] = moving_avg_high_gain_max_leo_amplitude_factor
    
    #JB OK
    leo_PkFWHM_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_PkFWHM_ch0[iloc[:-moving_average_window+1]] = moving_median_high_gain_beam_width
    
    #TAKE FROM ACTUAL PARTICLE TRACE WITH INCANDESENCE (NOT NEEDED IN ACTUAL LEO FIT?)
    leo_Base_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_Base_ch0[iloc[:-moving_average_window+1]] = moving_median_high_gain_base
    
    #JB OK
    leo_c2c_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_c2c_ch0[iloc[:-moving_average_window+1]] = moving_median_high_gain_c2c
    
    #JB OK
    leo_split_to_leo_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_split_to_leo_ch0[iloc[:-moving_average_window+1]] = moving_avg_high_gain_split_to_leo_pos
    
    """
    WHAT IS NEEDED:
    * GAUSSIAN WIDTH FOR ALL PARTICLES - FROM PROFILES (DONE)
    * PEAK POSITION FOR ALL PARTICLES - FROM C2C (DONE)
    * LAST BIN TO USE FOR LEO FIT (from split detector)  (DONE)
    * do also for LG.
    * TAKE BASELINE SCAT FROM ACTUAL TRACERS IN LEO_FIT
    """
    
    output_ds = my_binary.copy()
    
    output_ds['leo_AmpFactor_ch0'] = (('event_index'), np.zeros_like(leo_PkFWHM_ch0) + 
                                      np.nanmean(moving_avg_high_gain_max_leo_amplitude_factor))
    output_ds['leo_PkFWHM_ch0'] = (('event_index'), leo_PkFWHM_ch0)
    
    output_ds['leo_Base_ch0'] = (('event_index'), leo_Base_ch0)
    output_ds['leo_Base_ch0'] = output_ds['leo_Base_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    
    #distance from cross-to-centre (split point to laser maximum intensity). 
    #This comes from scattering only partilces
    output_ds['leo_PkPos_ch0'] = (('event_index'), leo_c2c_ch0)
    #interpolate to all particles (including, and most importantly, rBC containing particles)
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    #add split position to cross-to-centre distance to get location of the (estimated) 
    #location of the peak maximum. This is needed for the LEO-fit.
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'] + my_binary['PkSplitPos_ch3']
    
    #First add t_alpha_to_split data
    output_ds['leo_EndPos_ch0'] = (('event_index'), leo_split_to_leo_ch0)
    #interpolate to all particles
    output_ds['leo_EndPos_ch0'] = output_ds['leo_EndPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    #calculate the position at which the leo fits should end based on the split position 
    #(of all particles)
    output_ds['leo_EndPos_ch0'] = output_ds['leo_EndPos_ch0'] + output_ds['PkSplitPos_ch3']
    #cleaning up
    output_ds['leo_EndPos_ch0'].values = np.where(output_ds['leo_EndPos_ch0'].values < num_base_pts_2_avg,
                                                  np.nan, output_ds['leo_EndPos_ch0'].values)
    
    
#same as leo_PkPos_ch0
#    output_ds['leo_FtMaxPos_ch0'] = output_ds['leo_FtMaxPos_ch0'].interpolate_na(dim="event_index", 
#                                                        method="nearest", fill_value="extrapolate")
#    output_ds['leo_FtMaxPosAmpFactor_ch0'] = output_ds['leo_FtMaxPosAmpFactor_ch0'].interpolate_na(dim="event_index", 
#                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkFWHM_ch0'] = output_ds['leo_PkFWHM_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    
    
    #ToDo:
    #done: c2c time series (data from scattering particles but extrapolted to all paricles) 
    #done: leo_PkPos_ch0 (calculated from c2c) --> position of the max amplitude calculated from the split and extrapolated c2c
    #done: leo_FtEndPos_ch0 (from split)
    
    return output_ds

def leo_fit(my_binary,Globals=None):
    bins = my_binary['columns'].astype('float').values
    #number of points at the beginning to use for base line average
    num_base_pts_2_avg = 10
    #max_amplitude_fraction = 0.03
    
    # bl_scattering_ok = my_binary['ScatRejectKey'].values == 0
    # bl_only_scattering_particles = np.logical_and(bl_scattering_ok, 
    #                                               my_binary['PkHt_ch1'].values < Globals.IncanMinPeakHt1)

    # #Particles that only scatter light and ch0 not saturated 
    # bl_only_scattering_particles_ch0 = np.logical_and(my_binary['PkHt_ch0'].values < Globals.ScatMaxPeakHt1,
    #                                                   bl_only_scattering_particles)
    
    #split to peak height difference (in bins) for scattering only particles
    #split_to_peak_high_gain = my_binary['PkPos_ch0'].values - my_binary['PkSplitPos_ch3'].values
    #For particles with inandesence signal, set to NaN since the peak needn't 
    #be where the laser intensity is the highest, so se to NaN
    #split_to_peak_high_gain[~bl_only_scattering_particles_ch0] = np.nan
    
    #get the information about the gaussian fits
    pos = my_binary['leo_PkPos_ch0'].values
    #amplitude = my_binary['PkHt_ch0'].values
    width = my_binary['leo_PkFWHM_ch0'].values / 2.35482 #2*sqrt(2*log(2))
    data_ch0 = my_binary['Data_ch0'].values
    
    #mean of the first num_base_pts_2_avg points
    #leo_base_ch0 = np.mean(data_ch0[:, 0:num_base_pts_2_avg], axis=1)
    leo_base_ch0 = my_binary['leo_Base_ch0'].values
    data_ch0_sorted = np.sort(data_ch0[:, 0:num_base_pts_2_avg], axis=1)
    leo_base_ch0_ = np.min(data_ch0_sorted[:, 0:int(num_base_pts_2_avg)], axis=1)
        
    leo_fit_max_pos = my_binary['leo_EndPos_ch0'].astype(int).values
    leo_AmpFactor_ch0 = my_binary['leo_AmpFactor_ch0'].values
    leo_PkHt_ch0 = np.zeros_like(my_binary['PkHt_ch0'].values)*np.nan
    leo_PkHt_ch0_ = np.zeros_like(my_binary['PkHt_ch0'].values)*np.nan

    for i in range(my_binary.sizes['event_index']):
        max_value = data_ch0[i,:].max() - data_ch0[i,:].min()
        #bins_ = bins[num_base_pts_2_avg:leo_fit_max_pos[i]]
        bins_ = bins[leo_fit_max_pos[i]-3:leo_fit_max_pos[i]]
        if len(bins_) < 2:
            leo_PkHt_ch0[i] = np.nan
            continue
        #signals
        data_ch0_ = data_ch0[i, num_base_pts_2_avg:leo_fit_max_pos[i]]
        data_ch0_ = data_ch0[i, leo_fit_max_pos[i]-3:leo_fit_max_pos[i]]
        leo_coeff, var_matrix = curve_fit(
            lambda x, a: _gaus(x, a, pos[i], width[i], leo_base_ch0[i]), 
            bins_[:], data_ch0_[:], p0=[data_ch0[i,:].max()], maxfev=100, 
            ftol=1e-5, method='lm' ) #, bounds=(0, 1e6)) #, method='lm'
        leo_PkHt_ch0[i] = leo_coeff[0]
        leo_PkHt_ch0_[i] = (data_ch0[i, leo_fit_max_pos[i]] - leo_base_ch0_[i]) * leo_AmpFactor_ch0[i]
        
    output_ds = my_binary.copy()
    output_ds['leo_FtAmp_ch0'] = (('index'), leo_PkHt_ch0_)
    output_ds['leo_FtAmp_ch0_'] = (('index'), leo_PkHt_ch0_)
    output_ds['leo_Base_ch0'] = (('index'), leo_base_ch0)
    output_ds['leo_Base_ch0_'] = (('index'), leo_base_ch0_)

    
    return output_ds

#my_binary = pysp2.util.beam_shape(my_binary, beam_position_from='peak maximum', Globals=global_settings)
