import matplotlib.pyplot as plt
from matplotlib import dates
import numpy as np
from matplotlib.backend_bases import MouseButton

class DataEditor(object):

    """
    This class lets you remove data in a visual way. You can zoom and pan using
    the normal matplotlib functions.
    To delete data, press [d] and select (using the left button on the mouse) 
    the edges of the data that you want to remove. After two point the data to 
    be deleted will be marked. Press [y]es to accept removing that data.
    
    use the class like this:
    #make the plot
    fig, axs = plt.subplots(ncols=1)
    #plot the data
    line,=my_psds['NumConcScat'].plot(ax=axs)
    #edit the data
    browser=DataEditor(fig,axs,line)
    
    
    Parameters
    ----------
    fig: matplotlib.pyplot.figure
    axs: matplotlib.pyplot.figure.axes
    axs: matplotlib.pyplot.figure.axes
    line: matplotlib.lines.Line2D

    Returns
    -------
    DataEditor.edits['x_range']: list of numpy.arrays with boundaries in numpy.datetime64 
        format that should be deleted. E.g. 
        [[date1,date2],
         [date3,date4],
         ...]
    indicates that data should be deleted from date1 to date2 and from date3 to date4.
    The dates are in numpy.datetime64[us] format
        
        
    """
    def __init__(self,fig,axs,line):
            self.click_buffer=[]
            self.x_range_highlighted=False
            self.data_edits={'x_range':[],'y_range':[]}
            self.fig=fig
            self.axs=axs
            self.line=line
            self.fig.canvas.mpl_connect('key_press_event', self.on_press)
            self.axs.set_title('press [d] to delete x-range')
            self.fig.canvas.draw()

    def on_press(self, event):
        print('on_press invoked')
        if event.key not in ('d','y','n'):
            print('?')
            return
        if event.key == 'd':
            self.axs.set_title('pick x range to delete (left click)')
            self.fig.canvas.draw()
            self.binding_id=self.fig.canvas.mpl_connect('button_press_event', self.pick_x_range)
            #print('pick two points (left click)')
        elif event.key == 'y':
            if self.x_range_highlighted==True:
                x1,y1=self.click_buffer[0]
                x2,y2=self.click_buffer[1]
                self.click_buffer=[]
                self.remove_x_range_points[0].remove()
                self.axs.set_title('')
                self.fig.canvas.draw()
                self.x_range_highlighted=False
                #ToDo: remove tzinfo here!
                x_range=dates.num2date([min([x1,x2]),max([x1,x2])])
                #remove tzinfo and convert to numpy.datetime64 since psds 
                #xarray index is in numpy.datetime64
                x_range=[np.datetime64(d.replace(tzinfo=None)) for d in x_range]
                self.data_edits['x_range'].append(x_range)
                print(self.data_edits['x_range'])
                #update the original line with removed points set to nan
                lnx=dates.date2num(self.line.get_xdata())
                lny=self.line.get_ydata()
                bl=np.logical_and(lnx>(min(x1,x2)),lnx<(max(x1,x2)))
                lny[bl]=np.nan
                #this will propagate through the plot to the xarray variable
                #that made the plot in the first place.
                self.line.set_ydata(lny)
                self.fig.canvas.draw()
            else:
                print('pick your points first by pressing [d]')
        elif event.key == 'n':
            print('not saving')
            self.remove_x_range_points[0].remove()
            self.axs.set_title('')
            self.fig.canvas.draw()
            self.x_range_highlighted=False
            self.click_buffer=[]

    def pick_x_range(self, event):
        if event.button is MouseButton.LEFT:
            x = event.xdata
            y = event.ydata
            if len(self.click_buffer)==1:
                self.click_buffer.append((x,y))
                #print('Got 2nd point')
                self.fig.canvas.mpl_disconnect(self.binding_id)
                x1,y1=self.click_buffer[0]
                x2,y2=self.click_buffer[1]
                
                lnx=dates.date2num(self.line.get_xdata())
                lny=self.line.get_ydata()
                bl=np.logical_and(lnx>(min(x1,x2)),lnx<(max(x1,x2)))
                self.remove_x_range_points=self.axs.plot(lnx[bl],lny[bl],'.r',linestyle='')
                print(x1,x2,np.sum(bl))
                self.axs.set_title('remove %i points [y/n]'%(np.sum(bl)))
                self.fig.canvas.draw()
                self.x_range_highlighted=True
                #self.click_buffer=[]
            else:
                self.click_buffer=[(x,y)]
        if event.button is MouseButton.RIGHT:
            self.click_buffer=[]
            print('buffer reset')

