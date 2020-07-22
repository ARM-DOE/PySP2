import act
import numpy as np

from .DMTGlobals import DMTGlobals

def calc_diams_masses(input_ds, debug=True):
    """
    Calculates the scattering and incadescence diameters/BC masses for each particle.

    Parameters
    ----------
    input_ds: ACT Dataset
        The ACT Dataset containing the processed SP2 data.
    debug: boolean
        If true, print out particle rejection statistics

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
    PkHt_ch0 = np.nanmax(np.stack([input_ds['PkHt_ch0'].values, input_ds['FtAmp_ch0'].values]), axis=0)
    PkHt_ch4 = np.nanmax(np.stack([input_ds['PkHt_ch4'].values, input_ds['FtAmp_ch4'].values]), axis=0)
    accepted = np.logical_and.reduce((PkHt_ch0 > Globals.ScatMinPeakHt1,
                                      input_ds['PkFWHM_ch0'].values > Globals.ScatMinWidth,
                                      input_ds['PkFWHM_ch0'].values < Globals.ScatMaxWidth,
                                      input_ds['FtPos_ch0'].values < Globals.ScatMaxPeakPos,
                                      input_ds['FtPos_ch0'].values >= Globals.ScatMinPeakPos,
                                      np.greater(input_ds['FtAmp_ch0'].values, input_ds['PkFWHM_ch0'].values)))
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

    accepted_incand = np.logical_and(input_ds['PkHt_ch1'].values >= Globals.IncanMinPeakHt1,
                                     input_ds['PkHt_ch2'].values >= Globals.IncanMinPeakHt2)
    numMinCh2reject = np.sum(~accepted_incand)
    already_rejected = numMinCh2reject
    peak_width = input_ds['PkEnd_ch1'].values - input_ds['PkStart_ch1'].values
    peak_width[peak_width < 0] = 0.
    accepted_incand = np.logical_and(accepted_incand, peak_width >= Globals.IncanMinWidth)
    Ch1BaseWidthRejects = np.sum(~accepted_incand) - numMinCh2reject
    already_rejected += Ch1BaseWidthRejects
    accepted_incand = np.logical_and.reduce(
        (accepted_incand, ~np.logical_and(input_ds['PkHt_ch1'].values > Globals.IncanMaxPeakHt1,
         np.logical_or(input_ds['IncanRatioch1ch2'].values < Globals.IncanMinPeakRatio,
                       input_ds['IncanRatioch1ch2'].values > Globals.IncanMaxPeakRatio))))
    TempRatioRejects = np.sum(~accepted_incand) - already_rejected
    accepted_incand = np.logical_and.reduce(
        (accepted_incand, ~np.logical_and.reduce((
            input_ds['PkHt_ch1'].values > Globals.IncanMaxPeakHt1, input_ds['PkHt_ch5'].values < Globals.IncanMaxPeakHt1,
            input_ds['IncanRatioch5ch6'].values >= Globals.IncanMinPeakRatio,
            input_ds['IncanRatioch5ch6'].values <= Globals.IncanMaxPeakRatio))))

    TempRatioRejects += np.sum(~accepted_incand) - already_rejected
    already_rejected += TempRatioRejects
    accepted_incand = np.logical_and.reduce(
        (accepted_incand, ~np.logical_and(input_ds['PkHt_ch1'].values < Globals.IncanMaxPeakHt1,
                                          np.abs(input_ds['IncanPkOffsetch1ch2'].values) > Globals.IncanMaxPeakOffset)))
    IncanPeakDiffReject = np.sum(~accepted_incand) - already_rejected
    accepted_incand = np.logical_and.reduce(
        (accepted_incand, ~np.logical_and.reduce((
            input_ds['PkHt_ch1'].values > Globals.IncanMaxPeakHt1, input_ds['PkHt_ch5'].values < Globals.IncanMaxPeakHt1,
            np.abs(input_ds['IncanPkOffsetch5ch6'].values) > Globals.IncanMaxPeakOffset))))
    IncanPeakDiffReject += np.sum(~accepted_incand) - already_rejected
    already_rejected += IncanPeakDiffReject
    peak_loc1_reject = np.logical_and.reduce((accepted_incand, input_ds['PkHt_ch1'].values < Globals.IncanMaxPeakHt1,
                                      np.logical_or(input_ds['PkPos_ch1'].values < Globals.IncanMinPeakPos,
                                      input_ds['PkPos_ch1'].values > Globals.IncanMaxPeakPos)))
    peak_loc5_reject = np.logical_and.reduce(
        (accepted_incand, input_ds['PkHt_ch1'].values > Globals.IncanMaxPeakHt1,
         input_ds['PkHt_ch5'].values < Globals.IncanMaxPeakHt1,
         np.logical_or(input_ds['PkPos_ch5'].values < Globals.IncanMinPeakPos,
                       input_ds['PkPos_ch5'].values > Globals.IncanMaxPeakPos)))

    accepted_incand = np.logical_and.reduce((accepted_incand, ~peak_loc1_reject, ~peak_loc5_reject))
    PeakLoc1Rejects = np.sum(~accepted_incand) - already_rejected
    if debug:
        print("Number of incandescent particles accepted = %d" % np.sum(accepted_incand))
        print('Number rejected due to min peak height = %d' % numMinCh2reject)
        print("Number rejected due to ch1. base width = %d" % Ch1BaseWidthRejects)
        print("Number rejected due to peak ratio = %d" % TempRatioRejects)
        print("Number rejected due to peak difference = %d" % IncanPeakDiffReject)
        print("Number rejected due to peak location = %d" % PeakLoc1Rejects)

    Scat_not_sat = 1e-18*(Globals.c0Scat1 + Globals.c1Scat1*PkHt_ch0 + Globals.c2Scat1*PkHt_ch0**2)
    Scat_sat = 1e-18*(Globals.c0Scat2 + Globals.c1Scat2*PkHt_ch4 + Globals.c2Scat2*PkHt_ch4**2)
    Scatter = np.where(PkHt_ch0 < Globals.ScatMaxPeakHt1, Scat_not_sat, Scat_sat)
    Scatter = Scat_not_sat
    Scatter = np.where(accepted, Scatter, np.nan)

    output_ds = input_ds.copy()
    output_ds['Scatter'] = (('event_index'), Scatter)
    output_ds['logScatter'] = (('event_index'), np.log10(Scatter))
    output_ds['ScatDiaSO4'] = (('event_index'), 1000*(-0.015256 + 16.835*Scatter**0.15502))
    output_ds['ScatMassSO4'] = (('event_index'), 0.5236e-9*Globals.densitySO4*output_ds['ScatDiaSO4']**3)
    output_ds['ScatDiaBC50'] = (('event_index'), 1000*(0.013416 + 25.066*(Scatter**0.18057)))
    sootMass_not_sat = 1e-3*(Globals.c0Mass1 + Globals.c1Mass1*input_ds['PkHt_ch1'] + Globals.c2Mass1*input_ds['PkHt_ch1']**2)
    sootDiam_not_sat = (sootMass_not_sat/(0.5236e-9*Globals.densityBC))**(1./3.)
    sootMass_sat = 1e-3*(Globals.c0Mass2 + Globals.c1Mass2*input_ds['PkHt_ch5'] + Globals.c2Mass2*input_ds['PkHt_ch5']**2)
    sootDiam_sat = (sootMass_sat/(0.5236e-9*Globals.densityBC))**(1./3.)

    output_ds['sootMass'] = (('event_index'),
                             np.where(input_ds['PkHt_ch1'] > Globals.IncanMaxPeakHt1, sootMass_sat, sootMass_not_sat))
    output_ds['sootDiam'] = (('event_index'),
                             np.where(input_ds['PkHt_ch1'] > Globals.IncanMaxPeakHt1, sootDiam_sat, sootDiam_not_sat))
    output_ds['ScatDiaSO4'] = (('event_index'), output_ds['ScatDiaSO4'].where(accepted))
    output_ds['ScatMassSO4'] = (('event_index'), output_ds['ScatMassSO4'].where(accepted))
    output_ds['sootMass'] = output_ds['sootMass'].where(accepted_incand)
    output_ds['sootDiam'] = output_ds['sootDiam'].where(accepted_incand)

    return output_ds