import matplotlib.pyplot as plt
from matplotlib import dates
import numpy as np
from matplotlib.backend_bases import MouseButton

class DataEditor:
    """
    Visually edit bad data from the time series data.
    Start editing a variable by providing a particle size distribution series like this:
    edits=pysp2.vis.edit_variable(my_psds,'NumConcIncan')
    After the plot has opened, press d to delete a range using the mouse (left button). 
    When the selected range is highlighted, press [y]es to delet or [n]o to discard.
    
    """
    def __init__(self,my_psds,key,fig,axs,line):
            self.click_buffer=[]
            self.x_range_highlighted=False
            self.data_edits={'x_range':[],'y_range':[]}
            self.key=key
            self.my_psds=my_psds
            self.fig=fig
            self.axs=axs
            self.line=line
            #print('__init__',key)

    def on_press(self, event):
        if event.key not in ('d','y','n'):
            print('?')
            return
        if event.key == 'd':
            print('pick x range to delete (left click)')
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

def edit_variable(my_psds,variable):
    '''
    

    Parameters
    ----------
    my_psds : TYPE
        DESCRIPTION.
    variable : TYPE
        DESCRIPTION.

    Returns
    -------
    edits : TYPE
        DESCRIPTION.

    '''
    fig, axs = plt.subplots(ncols=1)
    line,=my_psds[variable].plot(ax=axs)
    axs.set_yscale('log')
    browser = DataEditor(my_psds,variable,fig,axs,line)
    fig.canvas.mpl_connect('key_press_event', browser.on_press)
    edits=browser.data_edits
    return edits