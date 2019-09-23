"""
Housekeeping data visualization routines
"""
import matplotlib.pyplot as plt

def plot_hk_timeseries(hk_df, var_name, **kwargs):
    """
    Plots a housekeeping data timeseries.


    Parameters
    ----------
    hk_df: ACT dataset
        The ACT dataset storing the housekeeping data.
    var_name: str or list
        The name of the variable(s) to plot

    Additional keyword arguments are passed into
    :func:`matplotlib.pyplot.plot`.

    Returns
    -------
    Display:
        The matplotlib axis handle of the plot.

    """



    return ax