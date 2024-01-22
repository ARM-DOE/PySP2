import numpy as np
from .peak_fit import _gaus
from scipy.optimize import curve_fit


def beam_shape(my_binary, beam_position_from='peak maximum', Globals=None):
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
           'split detector' = construct the beam profile around the split position.
           The split position is taken from the split detector. Not working yet.
           
    Returns
    -------
    coeff : numpy array with gaussian fit [amplitude, peakpos, width, base] 
            of the beam shape profile.
    beam_profile : numpy array with the beam profile calculated from the mean
                   of all the profiles. The array has as many data points as 
                   the input data.

    """
        
    num_base_pts_2_avg = 5
    bins=my_binary['columns']
    
    #boolean array for ok particles in the high gain channel
    scatter_high_gain_accepted = np.logical_and.reduce((
        my_binary['PkFWHM_ch0'].values > Globals.ScatMinWidth,
        my_binary['PkFWHM_ch0'].values < Globals.ScatMaxWidth,
        my_binary['PkHt_ch0'].values > Globals.ScatMinPeakHt1,
        my_binary['PkHt_ch0'].values < Globals.ScatMaxPeakHt1,
        my_binary['FtPos_ch0'].values < Globals.ScatMaxPeakPos,
        my_binary['FtPos_ch0'].values > Globals.ScatMinPeakPos))
    
    #boolean array that is True for particles that have not been triggered by 
    #the high gain incandesence channel
    no_incand_trigged = my_binary['PkHt_ch1'].values < Globals.IncanMinPeakHt1
    
    #find good events that only scatter light
    only_scattering_high_gain = np.logical_and.reduce((scatter_high_gain_accepted,
                                                       no_incand_trigged))
        
    print('High gain scattering particles for beam analysis :: ',
          np.sum(only_scattering_high_gain))
    
    #make an xarray of the purely scattering particles
    my_high_gain_scatterers = my_binary.sel(index = only_scattering_high_gain,
                                            event_index = only_scattering_high_gain)
    
    #numpy array for the normalized beam profiels
    my_high_gain_profiles = np.zeros((my_high_gain_scatterers.dims['index'],
                                    my_high_gain_scatterers.dims['columns'])) \
                                    * np.nan
    
    #weighted mean of beam peak position. Weight is scattering amplitude.
    high_gain_peak_pos = int(
        np.sum(np.multiply(my_high_gain_scatterers['PkPos_ch0'].values,
        my_high_gain_scatterers['PkHt_ch0'].values))/ \
                            np.sum(my_high_gain_scatterers['PkHt_ch0'].values))
    
    #loop through all particle events
    for i in my_high_gain_scatterers['event_index']:
        data = my_high_gain_scatterers['Data_ch0'].sel(event_index=i).values
        #base level
        base = np.mean(data[0:num_base_pts_2_avg])
        #peak height
        peak_height = data.max()-base
        #peak position
        peak_pos = data.argmax()
        #normalize the profile to range [0,1]
        profile = (data - base) / peak_height
        #distance to the mean beam peak position
        peak_difference = high_gain_peak_pos - peak_pos
        #insert so that the peak is at the right position (accounts for 
        #particles travelling at different speeds)
        if peak_difference > 0:
            my_high_gain_profiles[i, peak_difference:] = profile[:len(data) - 
                                                                peak_difference]
        elif peak_difference < 0:
            my_high_gain_profiles[i, :len(data)+peak_difference] = profile[-peak_difference:]
        else:
            my_high_gain_profiles[i, :] = profile

    #get the beam profile
    beam_profile = np.nanmean(my_high_gain_profiles, axis=0)
    #find values that are lower than 5% of the max value.
    low_values = np.argwhere(beam_profile < 0.05)
    #fit the gaussian curve to the beginning of the profile only. The tail 
    #can deviate from zero substantially and is not of interest.
    fit_to = low_values[low_values > high_gain_peak_pos].min()
    
    #initial guess
    p0 = np.array([beam_profile.max() - beam_profile.min(), 
                   np.argmax(beam_profile), 20., 
                   np.nanmin(beam_profile)]).astype(float)
    #fit gaussian curve
    #coeff[amplitude, peakpos, width , baseline]
    coeff, var_matrix = curve_fit(_gaus, bins[:fit_to], beam_profile[:fit_to], 
                                  p0=p0, method='lm', maxfev=40, ftol=1e-3)
    
    return coeff, beam_profile