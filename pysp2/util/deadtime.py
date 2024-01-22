import numpy as np
import xarray as xr

def deadtime(my_binary, my_ini):
    """
    Function to calculate Deadtime Bias in the SP2 data as published in
    Joshua. P. Schwarz, J. M. Katich, S. L. Lee, D. S. Thomson & L. A. Watts 
    (2022) “Invisible bias” in the single particle soot photometer due to 
    trigger deadtime, Aerosol Science and Technology, 56:7, 623-635, 
    DOI: 10.1080/02786826.2022.2064265 
    

    Parameters
    ----------
    my_binary : xarrray Dataset
        xarray Dataset that is the output of pysp2.util.gaussian_fit(my_sp2b, my_ini)
        Dataset with gaussian fits, peak heights etc.
        
    my_ini : dict
        The .ini file loaded as a dict. my_ini=pysp2.io.read_config(ini_file)

    Returns
    -------
    my_binary : xarrray Dataset
        Returns the input xarray dataset with one additional variable:
        DeadtimeRelativeBias. This can be used to correct for Deadtime Bias 
        in the SP2 at high concentrations. 
    """
    #numpy array to store the relative bias in
    bias = np.zeros(my_binary.dims['event_index'], )
    #scattering triggered
    scatter_triggered = np.any(np.isfinite(np.vstack((my_binary['PkHt_ch0'].values,
                                                    my_binary['PkHt_ch4'].values))), axis=0)
    #incandesence triggered
    incandesence_triggered = np.any(np.isfinite(np.vstack((my_binary['PkHt_ch1'].values,
                                                    my_binary['PkHt_ch5'].values))), axis=0)
    digitization_rate = int(my_ini['Acquisition']['Samples/Sec'])
    buffer_length = int(my_ini['Acquisition']['Scan Length']) #how many data points in one buffer
    buffer_length_seconds = buffer_length/digitization_rate #in seconds
    event_length_points = int(my_ini['Acquisition']['Points per Event']) 
    pretrigger_length_points = int(my_ini['Acquisition']['Pre-Trig Points']) 
    #intentialanny skipped particles, notation from Eq. 2 in doi:10.1080/02786826.2022.2064265
    S_S = int(my_ini['Acquisition']['1 of every']) 
    #find locations where the buffer changes (i.e. where the TimeWave values changes)
    #get unique time stamps
    unique_buffer_time_stamps = np.unique(my_binary['TimeWave'].values) 
    #find locations where those unique time stamps should be inserted
    ind = np.searchsorted(my_binary['TimeWave'].values, unique_buffer_time_stamps) 
    #add the end to the list of indexes if needed
    if ind[-1] < my_binary.dims['event_index']:
        ind = np.append(ind, my_binary.dims['event_index'])
    for i in range(len(ind)-1):
        #from index
        fr = ind[i]
        #to index
        to = ind[i+1]
        #number of scattering particle windows saved in the buffer
        N_S = scatter_triggered[fr:to].sum()
        #number of incandesence particle windows saved in the buffer
        N_I = incandesence_triggered[fr:to].sum()
        #length of time, in seconds, for a single measurement of the digitizer
        t_b = 1. / digitization_rate
        #"PW is the total number of data points (for each channel of data) in the window"
        P_W = event_length_points
        #Buffer length in seconds
        T_B = buffer_length_seconds     
        #Fraction triggered (Eq. 2) in doi:10.1080/02786826.2022.2064265
        F_T = (N_S * S_S+N_I) * P_W * t_b / T_B
        #pretrigger length data points
        P_PT = pretrigger_length_points
        #T_P is the total points in each window
        T_P = event_length_points
        #relative bias as in Eq. 3 in doi:10.1080/02786826.2022.2064265
        B_rel = -(P_PT / T_P) * F_T
        #save the bias value to the whole buffer
        bias[fr:to] = B_rel
    #add the bias containing array to the xarray
    my_binary['DeadtimeRelativeBias'] = xr.DataArray(bias, dims=['event_index'])
    return my_binary