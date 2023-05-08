from matplotlib import dates
from matplotlib.backend_bases import MouseButton
import numpy as np


class DataEditor(object):

    """
    This class lets you remove data in a visual way. You can zoom and pan using
    the normal matplotlib window buttons.
    To delete data, press [d] to start deleting. Then, select data (left mouse button)
    to remove data between the clicks. After selecting two points,
    the data about to be deleted will be marked red. Press [y]es to accept
    removing that data or [n]o if you don't want to delete it.

    Use the editor like this:
    #make the plot
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(ncols=1)
    #plot the data
    line,=my_psds['NumConcIncan'].plot(ax=axs)
    #edit the data
    browser=DataEditor(fig,axs,line)

    Parameters
    ----------
    fig: matplotlib.pyplot.figure
    ax: matplotlib.pyplot.figure.axes
    line: matplotlib.lines.Line2D

    Returns
    -------
    browser.data_edits['x_range']: list of numpy.arrays with boundaries in numpy.datetime64
        format that should be deleted. E.g.
        [[date1,date2],
         [date3,date4],
         ...]
    indicates that data should be deleted from date1 to date2 and from date3 to date4.
    The dates are in numpy.datetime64[us] format
    """
    def __init__(self, fig, ax, line):
        self.click_buffer = []
        self.x_range_highlighted = False
        self.data_edits = {'x_range': []}
        self.fig = fig
        self.ax = ax
        self.line = line
        self.fig.canvas.mpl_connect('key_press_event', self.on_press)
        self.ax.set_title('press [d] to delete x-range')
        self.fig.canvas.draw()

    def on_press(self, event):
        """
        This is the callback event handler for the DataEditor interface

        Parameters
        ----------
        event: matplotlib event
            The event to handle.

        """
        if event.key not in ('d', 'y', 'n'):
            print('unknown key pressed')
            return
        if event.key == 'd':
            self.ax.set_title('pick x range to delete (left click)')
            self.fig.canvas.draw()
            self.binding_id = self.fig.canvas.mpl_connect(
                'button_press_event', self.pick_x_range)
        elif event.key == 'y':
            if self.x_range_highlighted is True:
                x1, y1 = self.click_buffer[0]
                x2, y2 = self.click_buffer[1]
                self.click_buffer = []
                self.remove_x_range_points[0].remove()
                self.ax.set_title('')
                self.fig.canvas.draw()
                self.x_range_highlighted = False
                x_range = dates.num2date([min([x1, x2]), max([x1, x2])])
                # remove tzinfo and convert to numpy.datetime64 since psds
                # xarray index is in numpy.datetime64
                x_range = [np.datetime64(
                    d.replace(tzinfo=None)) for d in x_range]
                self.data_edits['x_range'].append(x_range)
                lnx = dates.date2num(self.line.get_xdata())
                lny = self.line.get_ydata()
                bl = np.logical_and(lnx > (min(x1, x2)), lnx < (max(x1, x2)))
                lny[bl] = np.nan
                # this will propagate through the plot to the xarray variable
                # that made the plot in the first place.
                self.line.set_ydata(lny)
                self.fig.canvas.draw()
            else:
                print('pick your points first by pressing [d]')
        elif event.key == 'n':
            print('discarding')
            self.remove_x_range_points[0].remove()
            self.ax.set_title('press [d] to delete x-range')
            self.fig.canvas.draw()
            self.x_range_highlighted = False
            self.click_buffer = []

    def pick_x_range(self, event):
        if event.button is MouseButton.LEFT:
            x = event.xdata
            y = event.ydata
            if len(self.click_buffer) == 1:
                self.click_buffer.append((x, y))
                self.fig.canvas.mpl_disconnect(self.binding_id)
                x1, y1 = self.click_buffer[0]
                x2, y2 = self.click_buffer[1]         
                lnx = dates.date2num(self.line.get_xdata())
                lny = self.line.get_ydata()
                bl = np.logical_and(lnx > (min(x1, x2)), lnx < (max(x1, x2)))
                self.remove_x_range_points = self.ax.plot(
                    lnx[bl], lny[bl], '.r', linestyle='')
                print(x1, x2, np.sum(bl))
                self.ax.set_title('remove %i points [y/n]' % (np.sum(bl)))
                self.fig.canvas.draw()
                self.x_range_highlighted = True
            else:
                self.click_buffer = [(x, y)]
        if event.button is MouseButton.RIGHT:
            self.click_buffer = []
            print('buffer reset. pick two new points.')
