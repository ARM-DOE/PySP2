import numpy as np
from scipy.optimize import curve_fit
import xarray as xr
from .peak_fit import _gaus, _do_fit_records
#from pysp2.util.peak_fit import _do_fit_records

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
        
    num_base_pts_2_avg = 10
    moving_average_window = 5
    max_amplitude_fraction = 0.03
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
    my_high_gain_profiles = np.zeros((my_high_gain_scatterers.dims['index'],
                                    my_high_gain_scatterers.dims['columns'])) \
                                    * np.nan
    
    my_high_gain_profiles_ = np.zeros_like(my_high_gain_profiles)*np.nan
    
    mean_high_gain_max_peak_pos = int(np.nanmean(my_high_gain_scatterers['PkPos_ch0'].values))
    mean_high_gain_split_pos = np.round(np.nanmean(my_high_gain_scatterers['PkSplitPos_ch3'].values)).astype(np.int32)
    #mean_high_gain_split_pos = int(np.nanmean(my_high_gain_scatterers['PkSplitPos_ch3'].values))
        
    #cross to center
    high_gain_c2c = my_high_gain_scatterers['FtPos_ch0'].values - my_high_gain_scatterers['PkSplitPos_ch3'].values
    high_gain_mean_c2c = np.mean(high_gain_c2c) #mean cross to centre
    
    high_gain_split_position = my_high_gain_scatterers['PkSplitPos_ch3'].values
    
    #loop through all particle events (THIS CAN BE MADE SMARTER WITH MUCH OF THE CALCULATIONS BEFORE THE LOOP)
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
        my_high_gain_profiles_[i,:] = profile
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

    # MOVING AVERAGE OF THE BEAM PROFILE TO FIND THE DISTANCE BETWEEN THE SPLIT POINT AND THE POINT IN THE 
    #LASER BEAM WHERE PARTICLES CAN START TO EVAPORATE
    
    #moving average of the beam shape with a window of moving_average_window
    moving_high_gain_profile_window = np.lib.stride_tricks.sliding_window_view(my_high_gain_profiles, 
                                                                     moving_average_window, axis=0)
    moving_avg_high_gain_profiles_ = np.nanmean(moving_high_gain_profile_window,axis=2)
    
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
        moving_avg_high_gain_max_leo_pos[i] = np.argwhere(moving_avg_high_gain_profiles[i,:] >= max_amplitude_fraction).min()-1
        moving_avg_high_gain_split_to_leo_pos[i] = moving_avg_high_gain_max_leo_pos[i] - high_gain_split_position[i]
        moving_avg_high_gain_max_leo_amplitude_factor[i] = 1./ moving_avg_high_gain_profiles[i, np.round(moving_avg_high_gain_max_leo_pos[i]).astype(int)]

    #cleaning up
    moving_avg_high_gain_max_leo_amplitude_factor = np.where(moving_avg_high_gain_max_leo_pos < num_base_pts_2_avg,
                                                             np.nan, moving_avg_high_gain_max_leo_amplitude_factor)
    
    moving_avg_high_gain_max_leo_pos = np.where(moving_avg_high_gain_max_leo_pos < num_base_pts_2_avg, 
                                                np.nan, moving_avg_high_gain_max_leo_pos)
    

    #moving average of beam width
    moving_high_gain_beam_width = np.lib.stride_tricks.sliding_window_view(my_high_gain_scatterers['PkFWHM_ch0'].values, 
                                                                     moving_average_window, axis=0)
    #moving_avg_high_gain_beam_width = np.nanmedian(moving_high_gain_beam_width,axis=1)
    moving_median_high_gain_beam_width = np.nanmedian(moving_high_gain_beam_width,axis=1)
    
    moving_high_gain_c2c = np.lib.stride_tricks.sliding_window_view(high_gain_c2c, 
                                                                     moving_average_window, axis=0)
    moving_median_high_gain_c2c = np.nanmedian(moving_high_gain_c2c,axis=1)

    
    leo_FtMaxPosAmpFactor_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan #will later be converted to int
    leo_FtMaxPosAmpFactor_ch0[iloc[:-moving_average_window+1]] = moving_avg_high_gain_max_leo_amplitude_factor
    
    leo_PkFWHM_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_PkFWHM_ch0[iloc[:-moving_average_window+1]] = moving_median_high_gain_beam_width
    
    leo_c2c_ch0 = np.zeros(scatter_high_gain_accepted.shape)*np.nan
    leo_c2c_ch0[iloc[:-moving_average_window+1]] = moving_median_high_gain_c2c
    
    #leo_PkPos_ch0 = np.zeros(scatter_high_gain_accepted.shape) * np.nan
    #leo_PkPos_ch0[iloc[:-moving_average_window+1]] = 
    
    output_ds = my_binary.copy()
    #output_ds['leo_FtMaxPos_ch0'] = (('event_index'), leo_FtMaxPos_ch0)
    output_ds['leo_FtMaxPosAmpFactor_ch0'] = (('event_index'), leo_FtMaxPosAmpFactor_ch0)
    output_ds['leo_PkFWHM_ch0'] = (('event_index'), leo_PkFWHM_ch0)
    #output_ds['leo_PkPos_ch0'] = (('event_index'), leo_PkPos_ch0)
    
    # ADD HERE: CALCULATE leo_PkPos_ch0 from split_point + moving_median_high_gain_c2c for all particles
    output_ds['leo_PkPos_ch0'] = (('event_index'), leo_c2c_ch0)
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'] + my_binary['PkSplitPos_ch3']
    
    output_ds['leo_FtMaxPos_ch0'] = output_ds['leo_FtMaxPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_FtMaxPosAmpFactor_ch0'] = output_ds['leo_FtMaxPosAmpFactor_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkFWHM_ch0'] = output_ds['leo_PkFWHM_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    output_ds['leo_PkPos_ch0'] = output_ds['leo_PkPos_ch0'].interpolate_na(dim="event_index", 
                                                        method="nearest", fill_value="extrapolate")
    
    #c2c time series (data from scattering particles but extrapolted to all paricles) 
    #leo_PkPos_ch0 (calculated from c2c) --> position of the max amplitude calculated from the split and extrapolated c2c
    #leo_FtMaxPos_ch0 (from split)
    
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
    split_to_peak_high_gain = my_binary['PkPos_ch0'].values - my_binary['PkSplitPos_ch3'].values
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
    #mean of the lowest 3 points
    leo_base_ch0 = np.mean(data_ch0[:, 0:num_base_pts_2_avg], axis=1)
    
    leo_fit_max_pos = my_binary['leo_FtMaxPos_ch0'].astype(int).values
    leo_FtMaxPosAmpFactor_ch0 = my_binary['leo_FtMaxPosAmpFactor_ch0'].values
    leo_PkHt_ch0 = np.zeros_like(my_binary['PkHt_ch0'].values)*np.nan
    leo_PkHt_ch0_ = np.zeros_like(my_binary['PkHt_ch0'].values)*np.nan

    for i in range(my_binary.dims['event_index']):
        #NAGOYA STYLE
        fractional_peak_height_ch0 = data_ch0[i, leo_fit_max_pos[i]] - leo_base_ch0[i]
        estimated_peak_height_ch0 = fractional_peak_height_ch0 * leo_FtMaxPosAmpFactor_ch0[i]
        leo_PkHt_ch0_[i] = estimated_peak_height_ch0
        max_value = data_ch0[i,:].max() - data_ch0[i,:].min()
        bins_ = bins[:leo_fit_max_pos[i]]
        #signals
        data_ch0_ = data_ch0[i, :leo_fit_max_pos[i]]
        leo_coeff, var_matrix = curve_fit(
            lambda x, a: _gaus(x, a, pos[i], width[i], leo_base_ch0[i]), 
            bins_[:], data_ch0_[:], p0=[max_value], maxfev=40, 
            ftol=1e-3, method='lm' ) #, bounds=(0, 1e6)) #, method='lm'
        leo_PkHt_ch0[i] = leo_coeff[0]
    
    output_ds = my_binary.copy()
    output_ds['leo_FtAmp_ch0'] = (('index'), leo_PkHt_ch0)
    output_ds['leo_FtAmp_ch0_'] = (('index'), leo_PkHt_ch0_)
    output_ds['leo_FtMaxPos_ch0'] = (('index'), leo_fit_max_pos)
    output_ds['leo_Base_ch0'] = (('index'), leo_base_ch0)
    
    #my_high_gain_scatterers['FtAmp_ch0'].plot(marker='.',linestyle='')
    #my_high_gain_scatterers['leo_FtAmp_ch0_'].plot(marker='+',linestyle='')
    #my_high_gain_scatterers['leo_FtAmp_ch0'].plot(marker='x',linestyle='')
    
    return output_ds

#my_binary = pysp2.util.beam_shape(my_binary, beam_position_from='peak maximum', Globals=global_settings)
