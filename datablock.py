import os
import re
import sys
import h5py
import pygrib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dateutil import parser
from scipy.interpolate import LinearNDInterpolator as lin_interp

soundvars = ['TMP', 'DWPT', 'HGT', 'RH', 'WDIR', 'WSPD']
preslist = np.arange(1000, 75, -25)
#preslist = np.arange(1000, 0, -200)
degCtoK = 273.15

import sharppy
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
def parseSPC(spc_file):
    ## read in the file
    data = np.array([l.strip() for l in spc_file.split('\n')])

    ## necessary index points
    title_idx = np.where( data == '%TITLE%')[0][0]
    start_idx = np.where( data == '%RAW%' )[0] + 1
    finish_idx = np.where( data == '%END%')[0]

    ## create the plot title
    data_header = data[title_idx + 1].split()
    location = data_header[0]
    time = data_header[1][:11]

    ## put it all together for StringIO
    full_data = '\n'.join(data[start_idx[0] : finish_idx[0]][:])
    sound_data = StringIO( full_data )

    ## read the data into arrays
    p, h, T, Td, wdir, wspd = np.genfromtxt( sound_data, delimiter=',', comments="%", unpack=True )
    rh = thermo.relh(p, T, Td)
    return [T, Td, h, rh, wdir, wspd, p]

def grbname(name, res, dt, fcst, ext):
    filedt = "%s%s%s_%s" % (dt.year, str(dt.month).zfill(2), 
                             str(dt.day).zfill(2), str(dt.hour*100).zfill(4) )
    label = "%s_%s_%s_%s.%s" % (name, str(res), filedt, str(fcst).zfill(3), ext)
    return label

def calc_dwpt(tempc, rh):
    parvap_pres = 6.112*rh*np.exp((17.67*tempc)/(243.12+tempc))
    
    logterm = np.log(parvap_pres / 611.2)
    dwpt_num = (17.67 - logterm)*degCtoK + 243.5*logterm
    dwpt_den = 17.67 - logterm
    return np.divide(dwpt_num, dwpt_den) - degCtoK

def calc_wspd(u, v):
    return np.sqrt(u**2 + v**2)

def calc_wdir(u, v):
    return (np.arctan2(u,v) + np.pi)*(180./np.pi)

def std_grid(gridrange, degres):
    lonmin, lonmax = gridrange[0,0], gridrange[0,1]
    latmin, latmax = gridrange[1,0], gridrange[1,1]
    longrid = np.linspace(lonmin, lonmax, num = 1 + int((lonmax - lonmin)/degres))
    latgrid = np.linspace(latmin, latmax, num = 1 + int((latmax - latmin)/degres))
    
    return longrid, latgrid

def interp_grid(grblon, grblat, grbval, stdlon, stdlat):
    flatlon=grblon.flatten()
    flatlat=grblat.flatten()
    flatval=grbval.flatten()
    func = lin_interp(list(zip(flatlon, flatlat ) ), flatval, fill_value = min(flatval))

    return func(*np.meshgrid(stdlon, stdlat));

class wxblocks(object):
    def __init__(self, gridrange, degres):
        """Class to create, modify, and store weather data"""
        self.gridlon, self.gridlat = std_grid(gridrange, degres)
        self.res = degres
        self.lonmin, self.lonmax = min(self.gridlon), max(self.gridlon)
        self.latmin, self.latmax = min(self.gridlat), max(self.gridlat)
        self.gridmesh = np.meshgrid(self.gridlon, self.gridlat)
        
        self.blocks = pd.DataFrame()

    def info(self):
        """Prints information about wxpanel"""
        if self.blocks.shape == (0, 0):
            print("Empty wxblock!")
        
        else:
            names = list(self.blocks.columns.names)[:-1]
            steps = list(self.blocks.columns.get_level_values(0).unique())
            params = list(self.blocks.columns.get_level_values(1).unique())
            
            datatype = ''
            if names[0]=='pres':
                datatype = 'sounding'
            elif names[0]=='datetime':
                datatype = 'time series'
                
            print("datatype: ", datatype )
            print("steps: ", steps )
            print("params: ", params )
            
    def grb2time(self, menu, datadir, varlist):
        dtlist = []
        for i, istep in menu.iterrows():
            dtlist.append(istep['datetime'])

        levels = pd.MultiIndex.from_product([dtlist, varlist, list(self.gridlon)],
                                            names = ['datetime', 'var', 'coord'])
        self.blocks = pd.DataFrame(index=list(self.gridlat), 
                                   columns=levels).sort_index(axis=1)
        
        for i, istep in menu.iterrows():
            file = grbname(*list(menu.iloc[i][:5]))
            gribpy = pygrib.open(datadir + file)
            dt = istep['datetime']
            
            for j, jvar in enumerate(varlist):
                short = re.split('(\d+)', jvar)[0]
                lev = int(re.split('(\d+)', jvar)[1] )
                griddata =  gribpy.select(shortName=short, level=lev )[0]
                
                gridvals = griddata.values
                griblat, griblon = griddata.latlons()
                gridded = interp_grid(griblon, griblat, 
                                      gridvals, 
                                      self.gridlon, self.gridlat)
                self.blocks[dt, jvar] = gridded
            #print(i, i%2)
            if i%2==0:
                print("finished ", str(i), ", ", str(len(dtlist)-i), " to go" )
        pass
            
    def grb2sound(self, gribpy):
        """
        Extracts an interpolated grid of values for the variables in 
        soundvars at the pressures in preslist from the data in gribfile
        """
        rawlist = ['gh', 't', 'r', 'u', 'v']
        gridlist = []
        for i, ivar in enumerate(rawlist):
            ilevel = []
            #gridlist.append(gribpy.select(shortName=ivar, level=tuple(preslist) ) )
            for j, jpres in enumerate(preslist):
                ilevel.append(gribpy.select(shortName=ivar, level=jpres )[0])
                
            gridlist.append(ilevel)
        levels = pd.MultiIndex.from_product([preslist, soundvars, list(self.gridlon)],
                                            names = ['pres', 'var', 'coord'])
        self.blocks = pd.DataFrame(index=list(self.gridlat), 
                                   columns=levels).sort_index(axis=1)
        
        griblat, griblon = gridlist[0][0].latlons()
        for i, ival in enumerate(preslist):
            for j, jval in enumerate(soundvars):
                gridvals = np.array([])
                
                if jval=='HGT':
                    gridvals = gridlist[0][i].values
                elif jval=='TMP':
                    gridvals = gridlist[1][i].values - degCtoK
                elif jval=='DWPT':
                    gridvals = calc_dwpt(gridlist[1][i].values - degCtoK,
                                         gridlist[2][i].values)
                elif jval=='WDIR':
                    gridvals = calc_wdir(gridlist[3][i].values,
                                         gridlist[4][i].values)
                elif jval=='WSPD':
                    gridvals = calc_wspd(gridlist[3][i].values,
                                         gridlist[4][i].values)
                    gridvals = gridvals*1.94384
                elif jval=='RH':
                    gridvals = gridlist[2][i].values
                    
                gridded = interp_grid(griblon, griblat, 
                                      gridvals, 
                                      self.gridlon, self.gridlat)
                self.blocks[ival, jval] = gridded
                #print("variable/pressure: ", jval, "/", ival, "mb")
            print("pressure: ", ival, " mb")
        pass
    
    def soundpoint(self, inlon, inlat):
        lonpt = self.gridlon[np.abs(self.gridlon-inlon).argmin()]
        latpt = self.gridlat[np.abs(self.gridlat-inlat).argmin()]
        
        if np.abs(lonpt-inlon)>self.res or np.abs(latpt-inlat)>self.res:
            return print("Error: Point not found on grid")
        else:
            idx = pd.IndexSlice
            soundpt = self.blocks.loc[latpt, idx[:,:,lonpt]]
            
            soundlist = []
            for i, val in enumerate(soundvars):
                soundlist.append(np.flip(np.array(soundpt.loc[idx[:,val]]), axis=0))
            
            soundlist.append(np.flip(np.array(soundpt.index.get_level_values(0).unique()), axis=0))
            return soundlist
    
    def save_h5(self, h5name):
        if os.path.exists(h5name):
            sys.exit('File already exists!')
    
        with h5py.File(h5name, 'w') as f:
            names = list(self.blocks.columns.names)[:-1]
            steps = list(self.blocks.columns.get_level_values(0).unique())
            params = list(self.blocks.columns.get_level_values(1).unique())
            
            if names[0]=='pres':
                f.attrs['type'] = 'sounding'
            elif names[0]=='datetime':
                f.attrs['type'] = 'time series'
                
            f.create_dataset('longitudes', data = self.blocks.columns.get_level_values(2).unique())
            f.create_dataset('latitudes', data = self.blocks.index)
            g = f.create_group('parameters')
                
            for i, param in enumerate(params):
                var = g.create_group(param)
                    
                for j, step in enumerate(steps):
                    var.create_dataset(str(step), data = self.blocks[step, param])
    
    def load_h5(self, h5name, params=[], start='', end=''):
        with h5py.File(h5name, 'r') as f:
        
            if len(params)==0:
                params = list(f['parameters'].keys())

            steplist = [f['parameters/' + x].keys() for x in params]
            steplist = [val for sublist in steplist for val in sublist]
            steplist = list(dict.fromkeys(steplist))
        
            gridlon = f['longitudes'][:]
            gridlat = f['latitudes'][:]
        
            blocktype = ''
            if f.attrs['type'] == 'sounding':
                steplist = list(map(int, steplist))
                blocktype = 'pres'
            
            if f.attrs['type'] == 'time series':
                steplist = [parser.parse(str(x)) for x in steplist]
                blocktype = 'datetime'
            
            if start=='':
                start=min(steplist)
            if end=='':
                end=max(steplist)
            
            dflevels = pd.MultiIndex.from_product([steplist, params, gridlon],
                                                  names = [blocktype, 'var', 'coord'])
            self.blocks = pd.DataFrame(index=gridlat, columns=dflevels).sort_index(axis=1)

            for i, param in enumerate(params):
                for j, step in enumerate(steplist):
                    data = f['parameters/' + param + '/' + str(step)]
                    self.blocks[step, param] = data[:]
    
    def plot2fig(self, step, var, fig, ax):
        plotmesh = ax.pcolormesh(self.gridmesh[0], self.gridmesh[1], self.blocks[step, var])
        plotbar = fig.colorbar(plotmesh, ax = ax)
        return plotmesh, plotbar
    
from matplotlib.widgets import Slider  # import the Slider widget
from ipywidgets import *

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

        self.fig, self.ax = plt.subplots(figsize = (10,8) )
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
