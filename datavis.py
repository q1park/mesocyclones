import numpy as np

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.transforms as transforms
import matplotlib.axis as maxis
import matplotlib.spines as mspines
import matplotlib.path as mpath
from matplotlib.projections import register_projection

from matplotlib.widgets import Slider  # import the Slider widget
from ipywidgets import *

from wxcalc import *

class animate_wxblocks(object):
    def __init__(self, wxblocks, var):
        istep = 0
        #title = (idatetime + timedelta(hours = istep)).strftime('%Y%m%d-%H')
        #self.wxblocks = wxblocks
        self.var = var 
        self.steps = np.array(wxblocks.columns.get_level_values(level=0).unique())
        self.data = np.array([wxblocks[x, self.var].values for x in self.steps])
        
        self.lon = np.array(wxblocks.columns.get_level_values(level=2).unique())
        self.lat = np.array(wxblocks.index.get_level_values(level=0).unique())

        self.fig, self.ax = plt.subplots(figsize = (6,5) )
        title = self.steps[istep]
        self.ax.set_title(title)

        self.sliderax = self.fig.add_axes([0.1, 0.02, 0.65, 0.04])
        self.slider = Slider(self.sliderax, 'Slice', 0, len(self.steps)-1, valinit=istep, valstep=1)
        self.slider.on_changed(self.update)

        self.im = self.ax.pcolormesh(self.lon, self.lat, self.data[istep], 
                                     vmin=self.data.min(), vmax=self.data.max())
        self.fig.colorbar(self.im, ax = self.ax)
        
    def update(self, i):
        newtitle = self.steps[int(i)]
        self.ax.set_title(newtitle)

        self.im.set_array(self.data[int(i),:-1,:-1].ravel())
        self.fig.canvas.draw()

    def show(self):
        plt.show()
        
        
        
def thetas(theta, presvals):
    return ((theta + thermo.ZEROCNK) / (np.power((1000. / presvals),thermo.ROCP))) - thermo.ZEROCNK


# indices = {'SBCAPE': [int(sfcpcl.bplus), 'J/kg'],\
#            'SBCIN': [int(sfcpcl.bminus), 'J/kg'],\
#            'SBLCL': [int(sfcpcl.lclhght), 'm AGL'],\
#            'SBLFC': [int(sfcpcl.lfchght), 'm AGL'],\
#            'SBEL': [int(sfcpcl.elhght), 'm AGL'],\
#            'SBLI': [int(sfcpcl.li5), 'C'],\
#            'MLCAPE': [int(mlpcl.bplus), 'J/kg'],\
#            'MLCIN': [int(mlpcl.bminus), 'J/kg'],\
#            'MLLCL': [int(mlpcl.lclhght), 'm AGL'],\
#            'MLLFC': [int(mlpcl.lfchght), 'm AGL'],\
#            'MLEL': [int(mlpcl.elhght), 'm AGL'],\
#            'MLLI': [int(mlpcl.li5), 'C'],\
#            'MUCAPE': [int(mupcl.bplus), 'J/kg'],\
#            'MUCIN': [int(mupcl.bminus), 'J/kg'],\
#            'MULCL': [int(mupcl.lclhght), 'm AGL'],\
#            'MULFC': [int(mupcl.lfchght), 'm AGL'],\
#            'MUEL': [int(mupcl.elhght), 'm AGL'],\
#            'MULI': [int(mupcl.li5), 'C'],\
#            '0-1 km SRH': [int(srh1km[0]), 'm2/s2'],\
#            '0-1 km Shear': [int(utils.comp2vec(sfc_1km_shear[0], sfc_1km_shear[1])[1]), 'kts'],\
#            '0-3 km SRH': [int(srh3km[0]), 'm2/s2'],\
#            'Eff. SRH': [int(effective_srh[0]), 'm2/s2'],\
#            'EBWD': [int(ebwspd), 'kts'],\
#            'PWV': [round(params.precip_water(prof), 2), 'inch'],\
#            'K-index': [int(params.k_index(prof)), ''],\
#            'STP(fix)': [round(stp_fixed, 1), ''],\
#            'SHIP': [round(ship, 1), ''],\
#            'SCP': [round(scp, 1), ''],\
#            'STP(cin)': [round(stp_cin, 1), '']}




def plotsounding(prof, parcel, title, fig, ax, extra = ''):
    name, ylabel, xlabel = '', 'P (mb)', 'T (C)'
    pcltype = 0
    
    if parcel == 'SB':
        pcltype = 1
        name = 'Surface Based'
    elif parcel == 'MU':
        pcltype = 3
        name = 'Most Unstable'
    elif parcel == 'ML':
        pcltype = 4
        name = 'Mixed Layer'
    else:
        print("Invalid parcel type")
        
    if len(extra) > 0:
        name = name + '\n' + extra
        
    pad = 5 # in points
      
    ax.annotate(title, xy=(0.5, 0.), xytext=(0, pad),
                xycoords=ax.title, textcoords='offset points',
                size='large', ha='center', va='baseline');
    ax.annotate(name, xy=(0.8, 0.875), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline');
    ax.annotate(ylabel, xy=(0., 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size='large', ha='center', va='center', rotation=90);
    ax.annotate(xlabel, xy=(0.5, -0.2), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline');
        
    pcl = params.parcelx( prof, flag=pcltype ) # Surface Parcel
    
    ax.grid(True)

    pmax = 1000
    pmin = 10
    dp = -10
    presvals = np.arange(int(pmax), int(pmin)+dp, dp)

    # plot the moist-adiabats
    for t in np.arange(-10,45,5):
        tw = []
        for p in presvals:
            tw.append(thermo.wetlift(1000., t, p))
        ax.semilogy(tw, presvals, 'k-', alpha=.2)



    # plot the dry adiabats
    for t in np.arange(-50,110,10):
        ax.semilogy(thetas(t, presvals), presvals, 'r-', alpha=.2)
    
    #plt.title(' OAX 140616/1900 (Observed)', fontsize=14, loc='left')
    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dicatated by the typical meteorological plot
    ax.semilogy(prof.tmpc, prof.pres, 'r', lw=2)
    ax.semilogy(prof.dwpc, prof.pres, 'g', lw=2)
    ax.semilogy(pcl.ttrace, pcl.ptrace, 'k-.', lw=2)

    # An example of a slanted line at constant X
    ax.axvline(0, color='b', linestyle='--')
    ax.axvline(-20, color='b', linestyle='--')

    # Disables the log-formatting that comes with semilogy
    ax.yaxis.set_major_formatter(plt.ScalarFormatter());
    ax.set_yticks(np.linspace(100,1000,10));
    ax.set_ylim(1050,100);

    ax.xaxis.set_major_locator(plt.MultipleLocator(10));
    ax.set_xlim(-50,50);
    
    
    
# The sole purpose of this class is to look at the upper, lower, or total
# interval as appropriate and see what parts of the tick to draw, if any.
class SkewXTick(maxis.XTick):
    def draw(self, renderer):
        if not self.get_visible(): return
        renderer.open_group(self.__name__)

        lower_interval = self.axes.xaxis.lower_interval
        upper_interval = self.axes.xaxis.upper_interval

        if self.gridOn and transforms.interval_contains(
                self.axes.xaxis.get_view_interval(), self.get_loc()):
            self.gridline.draw(renderer)

        if transforms.interval_contains(lower_interval, self.get_loc()):
            if self.tick1On:
                self.tick1line.draw(renderer)
            if self.label1On:
                self.label1.draw(renderer)

        if transforms.interval_contains(upper_interval, self.get_loc()):
            if self.tick2On:
                self.tick2line.draw(renderer)
            if self.label2On:
                self.label2.draw(renderer)

        renderer.close_group(self.__name__)


# This class exists to provide two separate sets of intervals to the tick,
# as well as create instances of the custom tick
class SkewXAxis(maxis.XAxis):
    def __init__(self, *args, **kwargs):
        maxis.XAxis.__init__(self, *args, **kwargs)
        self.upper_interval = 0.0, 1.0

    def _get_tick(self, major):
        return SkewXTick(self.axes, 0, '', major=major)

    @property
    def lower_interval(self):
        return self.axes.viewLim.intervalx

    def get_view_interval(self):
        return self.upper_interval[0], self.axes.viewLim.intervalx[1]


# This class exists to calculate the separate data range of the
# upper X-axis and draw the spine there. It also provides this range
# to the X-axis artist for ticking and gridlines
class SkewSpine(mspines.Spine):
    def _adjust_location(self):
        trans = self.axes.transDataToAxes.inverted()
        if self.spine_type == 'top':
            yloc = 1.0
        else:
            yloc = 0.0
        left = trans.transform_point((0.0, yloc))[0]
        right = trans.transform_point((1.0, yloc))[0]

        pts  = self._path.vertices
        pts[0, 0] = left
        pts[1, 0] = right
        self.axis.upper_interval = (left, right)


# This class handles registration of the skew-xaxes as a projection as well
# as setting up the appropriate transformations. It also overrides standard
# spines and axes instances as appropriate.
class SkewXAxes(Axes):
    # The projection must specify a name.  This will be used be the
    # user to select the projection, i.e. ``subplot(111,
    # projection='skewx')``.
    name = 'skewx'

    def _init_axis(self):
        #Taken from Axes and modified to use our modified X-axis
        self.xaxis = SkewXAxis(self)
        self.spines['top'].register_axis(self.xaxis)
        self.spines['bottom'].register_axis(self.xaxis)
        self.yaxis = maxis.YAxis(self)
        self.spines['left'].register_axis(self.yaxis)
        self.spines['right'].register_axis(self.yaxis)

    def _gen_axes_spines(self):
        spines = {'top':SkewSpine.linear_spine(self, 'top'),
                  'bottom':mspines.Spine.linear_spine(self, 'bottom'),
                  'left':mspines.Spine.linear_spine(self, 'left'),
                  'right':mspines.Spine.linear_spine(self, 'right')}
        return spines

    def _set_lim_and_transforms(self):
        """
        This is called once when the plot is created to set up all the
        transforms for the data, text and grids.
        """
        rot = 30

        #Get the standard transform setup from the Axes base class
        Axes._set_lim_and_transforms(self)

        # Need to put the skew in the middle, after the scale and limits,
        # but before the transAxes. This way, the skew is done in Axes
        # coordinates thus performing the transform around the proper origin
        # We keep the pre-transAxes transform around for other users, like the
        # spines for finding bounds
        self.transDataToAxes = self.transScale + (self.transLimits +
                transforms.Affine2D().skew_deg(rot, 0))

        # Create the full transform from Data to Pixels
        self.transData = self.transDataToAxes + self.transAxes

        # Blended transforms like this need to have the skewing applied using
        # both axes, in axes coords like before.
        self._xaxis_transform = (transforms.blended_transform_factory(
                    self.transScale + self.transLimits,
                    transforms.IdentityTransform()) +
                transforms.Affine2D().skew_deg(rot, 0)) + self.transAxes

# Now register the projection with matplotlib so the user can select
# it.
register_projection(SkewXAxes)