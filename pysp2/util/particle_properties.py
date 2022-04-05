import numpy as np
import xarray as xr
import datetime
import pandas as pd

from .DMTGlobals import DMTGlobals

def calc_diams_masses(input_ds, debug=True, factor=1.0, Globals=None):
    """
    Calculates the scattering and incandescence diameters/BC masses for each particle.

    Parameters
    ----------
    input_ds: ACT Dataset
        The ACT Dataset containing the processed SP2 data.
    debug: boolean
        If true, print out particle rejection statistics
    factor: float
        Multiply soot masses by this factor for AquaDag calibation. Use
        1.3 for NSA.
    Globals: DMTGlobals structure or None
        DMTGlobals structure containing calibration coefficients. Set to
        None to use default values for MOSAiC.

    Returns
    -------
    output_ds: ACT Dataset
        The ACT Dataset containing the scattering/incadescence diameters.
    """
    rejectMinScatTotal = 0
    rejectWidthTotal = 0
    rejectFatPeakTotal = 0
    rejectFtPosTotal = 0
    if Globals is None:
        Globals = DMTGlobals()

    PkHt_ch0 = np.nanmax(np.stack([input_ds['PkHt_ch0'].values, input_ds['FtAmp_ch0'].values]), axis=0)
    PkHt_ch4 = np.nanmax(np.stack([input_ds['PkHt_ch4'].values, input_ds['FtAmp_ch4'].values]), axis=0)
    accepted = np.logical_and.reduce((
        input_ds['PkFWHM_ch0'].values > Globals.ScatMinWidth,
        input_ds['PkFWHM_ch0'].values < Globals.ScatMaxWidth,
        input_ds['PkHt_ch0'].values > Globals.ScatMinPeakHt1,
        input_ds['FtPos_ch0'].values < Globals.ScatMaxPeakPos,
        input_ds['FtPos_ch0'].values > Globals.ScatMinPeakPos))
    numScatFlag = np.sum(accepted)

    rejectMinScatTotal += np.sum(PkHt_ch0 < Globals.ScatMinPeakHt1)
    rejectWidthTotal += np.sum(np.logical_or(
        input_ds['PkFWHM_ch0'].values < Globals.ScatMinWidth, input_ds['PkFWHM_ch0'].values > Globals.ScatMaxWidth))
    rejectFatPeakTotal += np.sum(np.less_equal(input_ds['FtAmp_ch0'].values, input_ds['PkFWHM_ch0'].values))
    rejectFtPosTotal += np.sum(np.logical_or(
        input_ds['FtPos_ch0'].values > Globals.ScatMaxPeakPos, input_ds['FtPos_ch0'].values < Globals.ScatMinPeakPos))

    if debug:
        print("Number of scattering particles accepted = %d" % numScatFlag)
        print("Number of scattering particles rejected for min. peak height = %d" % rejectMinScatTotal)
        print("Number of scattering particles rejected for peak width = %d" % rejectWidthTotal)
        print("Number of scattering particles rejected for fat peak = %d" % rejectFatPeakTotal)
        print("Number of scattering particles rejected for peak pos. = %d" % rejectFtPosTotal)

    PkHt_ch1 = input_ds['PkHt_ch1'].values
    PkHt_ch5 = input_ds['PkHt_ch5'].values
    width = input_ds['PkEnd_ch1'].values - input_ds['PkStart_ch1'].values
    accepted_incand = width >= Globals.IncanMinWidth
    accepted_incand = np.logical_and(accepted_incand,
                                     input_ds['PkHt_ch2'].values >= Globals.IncanMinPeakHt1)
    accepted_incand = np.logical_and(accepted_incand,
                                     input_ds['PkHt_ch1'].values >= Globals.IncanMinPeakHt1)
    sat_incand = np.logical_and(
        PkHt_ch1 >= Globals.IncanMaxPeakHt1, accepted_incand)
    unsat_incand = np.logical_and(accepted_incand, PkHt_ch1 < Globals.IncanMaxPeakHt1)
    unsat_incand = np.logical_and.reduce((unsat_incand,
                                          input_ds['IncanRatioch1ch2'].values > Globals.IncanMinPeakRatio,
                                          input_ds['IncanRatioch1ch2'].values < Globals.IncanMaxPeakRatio,
                                          np.abs(input_ds['IncanPkOffsetch1ch2'].values) < Globals.IncanMaxPeakOffset,
                                          input_ds['PkPos_ch1'].values > Globals.IncanMinPeakPos,
                                          input_ds['PkPos_ch1'].values < Globals.IncanMaxPeakPos))
    sat_incand = np.logical_and.reduce((sat_incand, PkHt_ch5 <= Globals.IncanMaxPeakHt1,
                                        input_ds['IncanRatioch5ch6'].values > Globals.IncanMinPeakRatio,
                                        input_ds['IncanRatioch5ch6'].values < Globals.IncanMaxPeakRatio,
                                        np.abs(input_ds['IncanPkOffsetch5ch6'].values) < Globals.IncanMaxPeakOffset,
                                        input_ds['PkPos_ch5'].values > Globals.IncanMinPeakPos,
                                        input_ds['PkPos_ch5'].values < Globals.IncanMaxPeakPos))
    accepted_incand = np.logical_or(unsat_incand, sat_incand)

    Scat_not_sat = 1e-18 * (Globals.c0Scat1 + Globals.c1Scat1*PkHt_ch0 + Globals.c2Scat1*PkHt_ch0**2)
    Scat_sat = 1e-18 * (Globals.c0Scat2 + Globals.c1Scat2*PkHt_ch4 + Globals.c2Scat2*PkHt_ch4**2)
    Scatter = np.where(PkHt_ch0 < Globals.ScatMaxPeakHt1, Scat_not_sat, Scat_sat)
    Scatter = np.where(accepted, Scatter, np.nan)

    output_ds = input_ds.copy()
    output_ds['Scatter'] = (('index'), Scatter)
    output_ds['logScatter'] = (('index'), np.log10(Scatter))
    output_ds['ScatDiaSO4'] = xr.DataArray(1000 * (-0.015256 + 16.835 * Scatter ** 0.15502),
                                           dims=['index'])
    output_ds['ScatMassSO4'] = xr.DataArray(0.5236e-9 * Globals.densitySO4 * output_ds['ScatDiaSO4'].data ** 3,
                                            dims=['index'])
    output_ds['ScatDiaBC50'] = xr.DataArray(1000*(0.013416 + 25.066*(Scatter**0.18057)),
                                            dims=['index'])

    sootMass_not_sat = factor * 1e-3 * (
        Globals.c0Mass1 + Globals.c1Mass1*PkHt_ch1 + Globals.c2Mass1*PkHt_ch1**2)
    sootDiam_not_sat = (sootMass_not_sat/(0.5236e-9*Globals.densityBC))**(1./3.)
    sootMass_sat = factor * 1e-3 * (
        Globals.c0Mass2 + Globals.c1Mass2*PkHt_ch5 + Globals.c2Mass2*PkHt_ch5**2)
    sootDiam_sat = (sootMass_sat/(0.5236e-9*Globals.densityBC))**(1./3.)
    sootMass_not_sat = np.where(accepted_incand, sootMass_not_sat, np.nan)
    sootMass_sat = np.where(accepted_incand, sootMass_sat, np.nan)
    sootMass = np.where(sat_incand, sootMass_sat, np.nan)
    sootMass[unsat_incand] = sootMass_not_sat[unsat_incand]
    sootDiam = np.where(sat_incand, sootDiam_sat, np.nan)
    sootDiam[unsat_incand] = sootDiam_not_sat[unsat_incand]

    output_ds['sootMass'] = (('index'), sootMass)
    output_ds['sootDiam'] = (('index'), sootDiam)
    output_ds['ScatDiaSO4'] = output_ds['ScatDiaSO4'].where(accepted)
    output_ds['ScatMassSO4'] = output_ds['ScatMassSO4'].where(accepted)

    return output_ds


def process_psds(particle_ds, hk_ds, config, deltaSize=0.005, num_bins=199, avg_interval=10):
    """
    Processes the Scattering and BC mass size distributions:

    Parameters
    ----------
    particle_ds: xarrray Dataset
        The xarray Dataset containing the particle statistics (.dat file information).
    hk_ds: xarray Dataset
        The xarray Dataset containing the housekeeping variables
    config: dict
        The .ini file loaded as a dict.
    deltaSize: float
        The size distribution bin width in microns.
    num_bins: int
        The number of size bins
    avg_interval: int
        The time in seconds to average the concentrations into.

    Returns
    -------
    psd_ds: xarray Dataset
        The xarray Dataset containing the time-averaged particle statistics.
    """
    DMTGlobal = DMTGlobals()
    time_bins = np.arange(hk_ds['Timestamp'].values[0],
                          hk_ds['Timestamp'].values[-1],
                          avg_interval)
    time_wave = particle_ds['DateTimeWaveUTC'].values
    if len(time_bins) == 0:
        time_bins = np.array([
            particle_ds['DateTimeWaveUTC'].values[0],
            particle_ds['DateTimeWaveUTC'].values[0] + avg_interval])

    flow = hk_ds['Sample Flow LFE'].values
    dcycle = hk_ds['Duty Cycle'].values
    for i in range(1, len(flow)):
        if flow[i] == 0:
            flow[i] = flow[i - 1]

    flow = flow / 60. * (dcycle / 100.)
    flow = np.where(flow == 0, 2, flow)
    ChmPress = hk_ds['Chamber Pressure'].values
    ChmTemp = hk_ds['Chamber Temp'].values
    ScatDiaBC50 = particle_ds['ScatDiaBC50'].values / 1000.
    ScatMassSO4 = particle_ds['ScatMassSO4'].values
    ScatDiaSO4 = particle_ds['ScatDiaSO4'].values / 1000.
    sootMass = particle_ds['sootMass'].values
    SizeIncandOnly = particle_ds['sootDiam'].values / 1000.
    SpecSizeBins = 0.01 + np.arange(0, num_bins, 1) * deltaSize
    ScatNumEnsembleBC = np.zeros((len(time_bins[:-1]), num_bins))
    ScatMassEnsembleBC = np.zeros_like(ScatNumEnsembleBC)
    IncanNumEnsemble = np.zeros((len(time_bins[:-1]), num_bins))
    IncanMassEnsemble = np.zeros_like(ScatNumEnsembleBC)
    ScatNumEnsemble = np.zeros((len(time_bins[:-1]), num_bins))
    ScatMassEnsemble = np.zeros_like(ScatNumEnsembleBC)
    scatter_accept = ~np.isnan(particle_ds['Scatter'].values)
    incand_accept = ~np.isnan(sootMass)
    try:
        OneOfEvery = int(config['Acquisition']['1 of Every'])
    except:
        OneOfEvery = 1
    
    IncandPos = particle_ds['PkPos_ch1'].values
    PkHt_ch1 = particle_ds['PkHt_ch1'].values
    incan_sat = PkHt_ch1 > DMTGlobal.IncanMaxPeakHt1
    IncandPos[incan_sat] = particle_ds['PkPos_ch5'].values[incan_sat]

    ScatPos = particle_ds['PkPos_ch0'].values
    PkHt_ch0 = particle_ds['PkHt_ch0'].values
    scat_sat = PkHt_ch0 > DMTGlobal.ScatMaxPeakHt1
    ScatPos[scat_sat] = particle_ds['PkPos_ch4'].values[scat_sat]

    NumFracBC = np.zeros_like(time_bins[:-1])
    NumFracBCSat = np.zeros_like(time_bins[:-1])
    NumConcScat1 = np.zeros_like(time_bins[:-1])
    NumConcIncan2 = np.zeros_like(time_bins[:-1])
    NumConcIncanScat = np.zeros_like(time_bins[:-1])
    NumConcIncanSat = np.zeros_like(time_bins[:-1])
    NumConcScatSat = np.zeros_like(time_bins[:-1])
    NumConcTotal = np.zeros_like(time_bins[:-1])
    MassScat2 = np.zeros_like(time_bins[:-1])
    MassScat2total = np.zeros_like(time_bins[:-1])
    MassIncand2 = np.zeros_like(time_bins[:-1])
    MassIncand2Sat = np.zeros_like(time_bins[:-1])
    for t in range(len(time_bins) - 1):
        parts_time = np.logical_and(
            time_wave >= time_bins[t], time_wave <= time_bins[t + 1])
        if np.sum(parts_time) == 0:
            continue
        times_hk = np.argwhere(np.logical_and(
                       hk_ds['Timestamp'].values >= time_bins[t],
                       hk_ds['Timestamp'].values <= time_bins[t + 1]))
        if len(times_hk) == 0:
            continue

        for i in range(num_bins):
            the_particles = np.logical_and.reduce((parts_time, scatter_accept,
                            ScatDiaBC50 >= SpecSizeBins[i] - deltaSize / 2,
                            ScatDiaBC50 < SpecSizeBins[i] + deltaSize / 2))
            ScatNumEnsembleBC[t, i] = OneOfEvery * np.sum(the_particles)
            the_particles = np.logical_and.reduce((parts_time, scatter_accept,
                            ScatDiaSO4 >= SpecSizeBins[i] - deltaSize / 2,
                            ScatDiaSO4 < SpecSizeBins[i] + deltaSize / 2))
            ScatNumEnsemble[t, i] = OneOfEvery * np.sum(the_particles)
            ScatMassEnsemble[t, i] = OneOfEvery * np.sum(ScatMassSO4[the_particles])
            the_particles = np.logical_and.reduce((parts_time, incand_accept,
                            SizeIncandOnly >= SpecSizeBins[i] - deltaSize / 2,
                            SizeIncandOnly < SpecSizeBins[i] + deltaSize / 2))
            IncanNumEnsemble[t, i] = OneOfEvery * np.sum(the_particles)
            IncanMassEnsemble[t, i] = OneOfEvery * np.sum(sootMass[the_particles])

        scat_parts = np.logical_and(scatter_accept, parts_time)
        incan_parts = np.logical_and(incand_accept, parts_time)
        ConcIncanCycle = OneOfEvery * np.sum(incan_parts)
        ConcTotalCycle = OneOfEvery * (np.sum(incan_parts) + np.sum(scat_parts))
        ConcScatCycle = OneOfEvery * np.sum(scat_parts)
        ConcScatSatCycle = OneOfEvery * np.sum(np.logical_and(scat_parts, incan_sat))
        ConcIncanScatCycle = OneOfEvery * np.sum(np.logical_and(scat_parts, incan_parts))
        ConcIncanSatCycle = OneOfEvery * np.sum(np.logical_and(incan_parts, incan_sat))
        massAvgScatCycle = OneOfEvery * np.sum(
            ScatMassSO4[np.argwhere(np.logical_and(scat_parts, ~scat_sat))])
        massAvgScatSatCycle = OneOfEvery * np.sum(
            ScatMassSO4[np.argwhere(np.logical_and(scat_parts, scat_sat))])
        massAvgIncandCycle = OneOfEvery * np.sum(
            sootMass[np.argwhere(np.logical_and(~incan_sat, incan_parts))])
        massAvgIncandSatCycle = OneOfEvery * np.sum(
            sootMass[np.argwhere(np.logical_and(incan_sat, incan_parts))])
        fracBCnum = ConcIncanCycle / ConcTotalCycle
        if ConcTotalCycle < 5:
            fracBCnum = 0
        fracBCnumSat = ConcIncanSatCycle / ConcTotalCycle
        if ConcIncanCycle < 5:
            fracBCnumSat = 0
        HalfPeriod = np.ones_like(times_hk).astype(float)
        if len(HalfPeriod) > 1:
            HalfPeriod[0] = 0.5        
            HalfPeriod[-1] = 0.5
        
        sum_flow = np.sum(HalfPeriod * flow[times_hk] * ChmPress[times_hk] / (273.15 + ChmTemp[times_hk]))
        FlowCycle = sum_flow * (DMTGlobal.TempSTP / DMTGlobal.PressSTP)
        
        if FlowCycle > 0:
            NumConcIncan2[t] = ConcIncanCycle / FlowCycle
            NumConcIncanSat[t] = ConcIncanSatCycle / FlowCycle
            NumConcIncanScat[t] = ConcIncanScatCycle / FlowCycle
            NumConcScatSat[t] = ConcScatSatCycle / FlowCycle
            NumConcScat1[t] = ConcScatCycle / FlowCycle
            NumConcTotal[t] = ConcTotalCycle / FlowCycle
            MassIncand2[t] = 1000 * massAvgIncandCycle / FlowCycle
            MassIncand2Sat[t] = 1000 * massAvgIncandSatCycle / FlowCycle
            MassScat2[t] = 1000 * massAvgScatCycle / FlowCycle
            MassScat2total[t] = 1000 * (massAvgScatSatCycle + massAvgScatCycle) / FlowCycle
            NumFracBC[t] = fracBCnum / FlowCycle
            NumFracBCSat[t] = fracBCnumSat / FlowCycle
            IncanNumEnsemble[t, :] = IncanNumEnsemble[t, :] / FlowCycle
            IncanMassEnsemble[t, :] = IncanMassEnsemble[t, :] / FlowCycle
            ScatNumEnsembleBC[t, :] = ScatNumEnsembleBC[t, :] / FlowCycle
            ScatMassEnsembleBC[t, :] = ScatMassEnsembleBC[t, :] / FlowCycle
            ScatNumEnsemble[t, :] = ScatNumEnsemble[t, :] / FlowCycle
            ScatMassEnsemble[t, :] = ScatMassEnsemble[t, :] / FlowCycle

    MassIncand2total = xr.DataArray(MassIncand2Sat + MassIncand2, dims=('time'))
    MassIncand2total.attrs["long_name"] = "Incandescence mass concentration (total)"
    MassIncand2total.attrs["standard_name"] = "mass_concentration"
    MassIncand2total.attrs["units"] = "ng m-3"
    base_time = pd.Timestamp('1904-01-01')
    time = np.array([(base_time + datetime.timedelta(seconds=x)).to_datetime64() for x in time_bins[:-1]])
    time = xr.DataArray(time, dims=('time'))
    time_wave = xr.DataArray(time_bins[:-1], dims=('time'))
    time_wave.attrs["long_name"] = "Time"
    time_wave.attrs["units"] = "seconds since 01-01-1904 00:00:00 UTC"

    NumConcIncan = xr.DataArray(NumConcIncan2, dims=('time'))
    NumConcIncan.attrs["long_name"] = "Total number concentration (incandescence)"
    NumConcIncan.attrs["standard_name"] = "incandescence_number_concentration"
    NumConcIncan.attrs["units"] = "cm-3"

    NumConcIncanScat = xr.DataArray(NumConcIncanScat, dims=('time'))
    NumConcIncanScat.attrs["long_name"] = "Number concentration detected by both incadesence and scatering detectors"
    NumConcIncanScat.attrs["standard_name"] = "number_concentration_detected_by_both"
    NumConcIncanScat.attrs["units"] = "cm-3"

    NumConcScat = xr.DataArray(NumConcScat1, dims=('time'))
    NumConcScat.attrs["long_name"] = "Total number concentration (scattering)"
    NumConcScat.attrs["standard_name"] = "scatter_number_concentration"
    NumConcScat.attrs["units"] = "cm-3"

    NumConcScatSat = xr.DataArray(NumConcScatSat, dims=('time'))
    NumConcScatSat.attrs["long_name"] = "Total number concentration (scattering, saturated)"
    NumConcScatSat.attrs["standard_name"] = "scatter_number_concentration"
    NumConcScatSat.attrs["units"] = "cm-3"

    NumConcTotal = xr.DataArray(NumConcTotal, dims=('time'))
    NumConcTotal.attrs["long_name"] = "Total number concentration"
    NumConcTotal.attrs["standard_name"] = "number_concentration"
    NumConcTotal.attrs["units"] = "cm-3"

    MassIncand = xr.DataArray(MassIncand2total, dims=('time'))
    MassIncand.attrs["long_name"] = "Incandescence mass concentration"
    MassIncand.attrs["standard_name"] = "mass_concentration"
    MassIncand.attrs["units"] = "ng m-3"

    MassIncand2Sat = xr.DataArray(MassIncand2Sat, dims=('time'))
    MassIncand2Sat.attrs["long_name"] = "Incandescence mass concentration (saturated)"
    MassIncand2Sat.attrs["standard_name"] = "mass_concentration"
    MassIncand2Sat.attrs["units"] = "ng m-3"

    MassScat2 = xr.DataArray(MassScat2, dims=('time'))
    MassScat2.attrs["long_name"] = "Scattering mass concentration (cm-3)"
    MassScat2.attrs["standard_name"] = "scatter_mass_concentration"
    MassScat2.attrs["units"] = "ng cm-3"

    NumFracBC = xr.DataArray(NumFracBC, dims=('time'))
    NumFracBC.attrs["long_name"] = "Number fraction of black carbon"
    NumFracBC.attrs["standard_name"] = "number_fraction_black_carbon"
    NumFracBC.attrs["units"] = "1"

    ScatNumEnsemble = xr.DataArray(ScatNumEnsemble, dims=('time', 'num_bins'))
    ScatNumEnsemble.attrs["long_name"] = "Scattering number distribution"
    ScatNumEnsemble.attrs["standard_name"] = "scattering_number_distribution"
    ScatNumEnsemble.attrs["units"] = "cm-3 per bin"
   
    IncanNumEnsemble = xr.DataArray(IncanNumEnsemble, dims=('time', 'num_bins'))
    IncanNumEnsemble.attrs["long_name"] = "Incandescence number distribution"
    IncanNumEnsemble.attrs["standard_name"] = "incandescence_number_distribution"
    IncanNumEnsemble.attrs["units"] = "cm-3 per bin"

    ScatMassEnsemble = xr.DataArray(ScatMassEnsemble, dims=('time', 'num_bins'))
    ScatMassEnsemble.attrs["long_name"] = "Scattering mass distribution"
    ScatMassEnsemble.attrs["standard_name"] = "scattering_mass_distribution"
    ScatMassEnsemble.attrs["units"] = "ng m-3 per bin"

    IncanMassEnsemble = xr.DataArray(IncanMassEnsemble, dims=('time', 'num_bins'))
    IncanMassEnsemble.attrs["long_name"] = "Incandesence mass distribution"
    IncanMassEnsemble.attrs["standard_name"] = "incadesence_mass_distribution"
    IncanMassEnsemble.attrs["units"] = "ng m-3 per bin"

    ScatNumEnsembleBC = xr.DataArray(ScatNumEnsembleBC, dims=('time', 'num_bins'))
    ScatNumEnsembleBC.attrs["long_name"] = "Scattering number distribution (black carbon)"
    ScatNumEnsembleBC.attrs["standard_name"] = "scattering_number_distribution (black carbon)"
    ScatNumEnsembleBC.attrs["units"] = "cm-3 per bin"

    SpecSizeBins = xr.DataArray(SpecSizeBins, dims=('num_bins'))
    SpecSizeBins.attrs["long_name"] = "Spectra size bin centers"
    SpecSizeBins.attrs["standard_name"] = "particle_diameter"
    SpecSizeBins.attrs["units"] = "nm"
    psd_ds = xr.Dataset({'time': time,
                         'TimeWave': time_wave,
                         'NumConcIncan': NumConcIncan,
                         'NumConcIncanScat': NumConcIncanScat,
                         'NumConcTotal': NumConcTotal,
                         'NumConcScatSat': NumConcScatSat,
                         'NumConcScat': NumConcScat,
                         'MassIncand2': MassIncand2,
                         'MassIncand2Sat': MassIncand2Sat,
                         'MassIncand2total': MassIncand2total,
                         'MassScat2': MassScat2,
                         'ScatNumEnsemble': ScatNumEnsemble,
                         'ScatMassEnsemble': ScatMassEnsemble,
                         'IncanNumEnsemble': IncanNumEnsemble,
                         'IncanMassEnsemble': IncanMassEnsemble,
                         'ScatNumEnsembleBC': ScatNumEnsembleBC,
                         'NumFracBC': NumFracBC})

    return psd_ds
