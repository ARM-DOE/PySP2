"""
Housekeeping data visualization routines
"""
import matplotlib.pyplot as plt

def plot_hk_timeseries(hk_df, var_name, **kwargs):
    """
    Plots a housekeeping data timeseries.


    Parameters
    ----------
    hk_df: pandas DataFrame
        The pandas DataFrame storing the housekeeping data.
    var_name: str or list
        The name of the variable(s) to plot

    Additional keyword arguments are passed into
    :func:`matplotlib.pyplot.plot`.

    Returns
    -------
    ax: matplotlib axis handle
        The matplotlib axis handle of the plot.

    """

    if isinstance(var_name, str):
        ax = hk_df[var_name].plot(**kwargs)
        ax.set_ylabel(var_name)
    elif isinstance(var_name, list):
        ax = hk_df[var_name[0]].plot(**kwargs, label=var_name[0])
        for v in var_name[1:]:
            hk_df[v].plot(**kwargs, label=v, ax=ax)
        ax.legend()
    else:
        raise TypeError("var_name must be str or list!")

    return ax