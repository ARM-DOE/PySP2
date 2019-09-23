import act
from .DMTGlobals import DMTGlobals

def calc_diams_masses(input_ds):
    """
    Calculates the scattering and incadescence diameters/BC masses for each particle.

    Parameters
    ----------
    input_ds: ACT Dataset
        The ACT Dataset containing the processed SP2 data.

    Returns
    -------
    output_ds: ACT Dataset
        The ACT Dataset containing the scattering/incadescence diameters.
    """
    rejectMinScatTotal = 0
    rejectWidthTotal = 0
    rejectFatPeakTotal = 0
    rejectFtPosTotal = 0
    Globals = DMTGlobals()
    max_pk0 = np.max(np.stack([input_ds['PkHt_ch0'].values, input_ds['FtAmp_ch0'].values)]), axis=0)
    max_pk4 = np.max(np.stack([input_ds['PkHt_ch4'].values, input_ds['FtAmp_ch4'].values)]), axis=0)
    accepted = np.logical_and.reduce((input_ds['PkHt_ch0'] > Globals.ScatMinPeakHt,
                                      input_ds['PkFWHM_ch0'] > Globals.ScatMinWidth,
                                      input_ds['PkFWHM_ch0'] < Globals.ScatMaxWidth,
                                      input_ds['FtPos_ch0'] < Globals.ScatMaxPeakPos,
                                      input_ds['FtPos_ch0'] >= Globals.ScatMinPeakPos,
                                      np.greater(input_ds['FtAmp_ch0'], input_ds['PkFWHM_ch0'])))
    numScatFlag = np.sum(accepted)
    rejectMinScatTotal += np.sum(input_ds['PkHt_Ch0'] < Globals.ScatMinPeakHt)
    rejectWidthTotal += np.sum(np.logical_or(
        input_ds['PkFWHM_ch0'] < Globals.ScatMinWidth, input_ds['PkFWHM_ch0'] > Globals.ScatMaxWidth))
    rejectFatPeakTotal += np.sum(np.less_equal(input_ds['FtAmp_ch0'], input_ds['PkFWHM_ch0']))
    rejectFtPosTotal += np.sum(np.logical_or(
        input_ds['FtPos_ch0'] > Globals.ScatMaxPeakPos, input_ds['FtPos_ch0'] < Globals.ScatMinPeakPos))
    accepted_incand = np.logical_and(input_ds['PkHt_ch1'] > Globals.IncanMinPeakHt,
                                     input_ds['PkHt_ch2'] > Globals.IncanMinPeakHt)
    numMinCh2reject = np.sum(input_ds['PkHt_ch2'] < Globals.IncanMinPeakHt)
    saturated_incand = np.logical_and(
        accepted_incand, np.logical_or(input_ds['PkHt_ch1'] > Globals.IncanMaxPeakHt,
                                       input_ds['PkHt_ch5'] > Globals.IncanMaxPeakHt)
    IncandPos = np.where(input_ds['PkHt_ch1'] < Globals.IncanMaxPeakHt, input_ds['PkPos_ch1'], input_ds['PkPos_ch5'])
    IncandPos[~accepted_incand] = np.nan
    saturated_scattering = np.logical_or(input_ds['PkHt_ch0'] > Globals.ScatMaxPeakHt,
                                         input_ds['PkHt_ch4'] > Globals.ScatMaxPeakHt)
    use_ch0 = np.logical_and(input_ds['PkHt_ch0'] < Globals.ScatMaxPeakHt, input_ds['PkHt_ch0'] > 5000.)
    ScatPos = np.where(use_ch0, input_ds['FtPos_ch0'], input_ds['FtPos_ch4'])
    ScatPos[~accepted] = np.nan
    ScatPos[np.logical_and(~use_ch0, input_ds['PkHt_ch4'] < 2500.)] = np.nan
    accepted_incand = np.logical_and(accepted_incand, input_ds['PkFWHM_ch1'] >= Globals.IncanMinWidth)
    Ch1BaseWidthRejects = np.sum(input_ds['PkFWHM_ch1'] < Globals.IncanMinWidth)
    accepted_incand = np.logical_and.reduce(
        (accepted_incand, input_ds['IncanRatioCh1Ch2'] >= Globals.IncanMinPeakRatio,
         input_ds['IncanRatioCh1Ch2'] <= Globals.IncanMaxPeakRatio))
    TempRatioRejects = np.sum(np.logical_or(input_ds['IncanRatioCh1Ch2'] < Globals.IncanMinPeakRatio,
         input_ds['IncanRatioCh1Ch2'] > Globals.IncanMaxPeakRatio))
    accepted_incand = np.logical_and.reduce(
        (accepted_incand, input_ds['IncanRatioCh5Ch6'] >= Globals.IncanMinPeakRatio,
         input_ds['IncanRatioCh5Ch6'] <= Globals.IncanMaxPeakRatio))
    TempRatioRejects += np.sum(np.logical_or(input_ds['IncanRatioCh5Ch6'] < Globals.IncanMinPeakRatio,
         input_ds['IncanRatioCh5Ch6'] > Globals.IncanMaxPeakRatio))
    accepted_incand = np.logical_and.reduce(
        (accepted_incand, ~np.logical_and(input_ds['PkHt_ch1'] < Globals.IncanMaxPeakHt,
                                          np.abs(input_ds['IncanPeakOffsetCh1Ch2']) < Globals.IncanMaxPeakOffset)))
    IncanPeakDiffReject = np.sum(np.logical_and(input_ds['PkHt_ch1'] < Globals.IncanMaxPeakHt,
                                          np.abs(input_ds['IncanPeakOffsetCh1Ch2']) > Globals.IncanMaxPeakOffset))
    accepted_incand = np.logical_and.reduce(
        (accepted_incand, ~np.logical_and(input_ds['PkHt_ch1'] > Globals.IncanMaxPeakHt,
                                          np.abs(input_ds['IncanPeakOffsetCh5Ch6']) < Globals.IncanMaxPeakOffset)))
    IncanPeakDiffReject += np.sum(np.logical_and(input_ds['PkHt_ch1'] > Globals.IncanMaxPeakHt,
                                          np.abs(input_ds['IncanPeakOffsetCh5Ch6']) > Globals.IncanMaxPeakOffset))
    peak_loc1_reject = np.logical_and(accepted_incand, np.logical_or(input_ds['PkPos_ch1'] < Globals.IncanMinPeakPos,
                                      input_ds['PkPos_ch1'] > Globals.IncanMaxPeakPos))
    peak_loc5_reject = np.logical_and.reduce(
        (accepted_incand, input_ds['PkHt_ch1'] > Globals.IncanMaxPeakHt, input_ds['PkHt_ch5'] < Globals.IncanMaxPeakHt,
         np.logical_or(input_ds['PkPos_ch5'] < Globals.IncanMinPeakPos, input_ds['PkPos_ch5'] > Globals.IncanMaxPeakPos)))
    PeakLoc1Rejects = np.sum(peak_loc1_reject) + np.sum(peak_loc5_reject)
    accepted_incand = np.logical_and.reduce((accepted_incand, ~peak_loc1_reject, ~peak_loc5_reject))
    Scat_not_sat = 1e-18*(Globals.c0Scat1 + c1Scat1*input_ds['PkHt_ch0'] + c2Scat1*input_ds['PkHt_ch0']**2)
    Scat_sat = 1e-18*(Globals.c0Scat2 + c1Scat2*input_ds['PkHt_ch4'] + c2Scat2*input_ds['PkHt_ch4']**2)
    Scatter = np.where(input_ds['PkHt_ch0'] < Globals.ScatMaxPeakHt, Scat_not_sat, Scat_sat)
    Scatter[~accepted] = 0
    input_ds['ScatDiaSO4'] = 1000*(-0.015256 + 16.835*Scatter**0.15502)
    input_ds['ScatMassSO4'] = 0.5236e-9*Globals.densitySO4*input_ds['ScatDiaSO4']**3
    input_ds['ScatDiaBC50'] = 1000*(0.013416 + 25.066*(Scatter**0.18057))
    sootMass_not_sat = 1e-3*(Globals.c0Mass2 + Globals.c1Mass2*input_ds['PkHt_ch5'] + Globals.c)





