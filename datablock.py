import os
import re
import sys
import h5py
import pygrib
import numpy as np
import pandas as pd

from dateutil import parser
from scipy.interpolate import LinearNDInterpolator as lin_interp

from wxcalc import *
from dataform import *

soundvars = ['TMP', 'DWPT', 'HGT', 'RH', 'WDIR', 'WSPD', 'PRES']
soundunits = [r'$C^o$', r'$km$', r'$gpm$', '%', r'deg', r'$kts$', r'$mb$']
preslist = np.arange(1000, 75, -25)
#preslist = np.arange(1000, 0, -100)

def std_grid(gridrange, gridres):
    lonmin, lonmax = gridrange[0,0], gridrange[0,1]
    latmin, latmax = gridrange[1,0], gridrange[1,1]
    longrid = np.linspace(lonmin, lonmax, num = 1 + int((lonmax - lonmin)/gridres))
    latgrid = np.linspace(latmin, latmax, num = 1 + int((latmax - latmin)/gridres))
    
    return longrid, latgrid

def interp_grid(grblon, grblat, grbval, stdlon, stdlat):
    flatlon=grblon.flatten()
    flatlat=grblat.flatten()
    flatval=grbval.flatten()
    func = lin_interp(list(zip(flatlon, flatlat ) ), flatval, fill_value = min(flatval))

    return func(*np.meshgrid(stdlon, stdlat));

class wxblocks(object):
    def __init__(self, gridrange, gridres):
        """Class to create, modify, and store weather data"""
        self.gridrange = gridrange
        self.gridres = gridres
        
        self.gridlon, self.gridlat = std_grid(self.gridrange, self.gridres)
        self.blocks = pd.DataFrame()
        self.raw = [pd.DataFrame(), datetime(9999,1,1,0)]

    def info(self):
        """Prints information about wxpanel"""
        print("###########################")
        print("### grid info: ")
        print("range / res: ")
        print(self.gridrange, self.gridres )
        print("###########################\n")
        if self.blocks.shape == (0, 0):
            print("###########################")
            print("### no time series data ###")
            print("###########################\n")
        else:
            names = list(self.blocks.columns.names)[:-1]
            datetimes = list(self.blocks.columns.get_level_values(0).unique())
            variables = list(self.blocks.columns.get_level_values(1).unique())
            print("###########################")
            print("### block info: ")
            print("datetimes: ", datetimes )
            print("variables: ", variables )
            print("###########################\n")
        
        if self.raw[0].shape == (0, 0):
            print("############################")
            print("### no raw variable data ###")
            print("############################\n")
        
        else:
            names = list(self.raw[0].columns.names)[:-1]
            pressures = list(self.raw[0].columns.get_level_values(0).unique())
            variables = list(self.raw[0].columns.get_level_values(1).unique())
            timestamp = str(self.raw[1])
            print("############################")
            print("### raw info: ")
            print("timestamp: ", timestamp)
            print("pressures: ", pressures )
            print("params: ", variables )
            print("############################")
    
    def dfstruct(self, level1, level2):
        levelnames = []
        if type(level1[0]) == datetime:
            levelnames = ['datetime', 'variable', 'coordinate']
        else:
            levelnames = ['pressure', 'variable', 'coordinate']
            
            
        levels = pd.MultiIndex.from_product([level1, level2, list(self.gridlon)],
                                            names = levelnames)
        return pd.DataFrame(index=list(self.gridlat), 
                            columns=levels).sort_index(axis=1)
            
    
    def grb2block(self, menu, datadir, varlist):
        self.blocks = self.dfstruct(list(menu['datetime']), varlist)
        
        for i, istep in menu.iterrows():
            file = grbname(*list(menu.iloc[i][:5]))
            gribpy = pygrib.open(datadir + file)
            dt = istep['datetime']
            
            for j, jvar in enumerate(varlist):
                name, strlevel = re.split('(\d+)', jvar)[:-1]
                gribitem =  gribpy.select(shortName=short, level=int(strlevel) )[0]

                gridded = interp_grid(gribitem.latlons()[1], gribitem.latlons()[0], 
                                      gribitem.values, self.gridlon, self.gridlat)
                self.blocks[dt, jvar] = gridded
            if i%2==0:
                print("finished ", str(i), ", ", str(len(dtlist)-i), " to go" )
        pass

    def grb2raw(self, gribfile, datadir):
        """
        Extracts an interpolated grid of values for the variables in 
        soundvars at the pressures in preslist from the data in gribfile
        """
        rawlist = ['gh', 't', 'r', 'u', 'v', 'w']
        self.raw[0] = self.dfstruct(preslist, rawlist)
        self.raw[1] = parsefile(gribfile)[2]
        gribpy = pygrib.open(datadir + gribfile )
        
        for j, jpres in enumerate(preslist):
            for i, ivar in enumerate(rawlist):
                gribitem = gribpy.select(shortName=ivar, level=jpres )[0]
                
                gridded = interp_grid(gribitem.latlons()[1], gribitem.latlons()[0], 
                                      gribitem.values, self.gridlon, self.gridlat)
                self.raw[0][jpres, ivar] = gridded
            print("added level: ", jpres, " mb")
        pass    
    
    def raw2sound(self, inlon, inlat):
        lonpt = self.gridlon[np.abs(self.gridlon-inlon).argmin()]
        latpt = self.gridlat[np.abs(self.gridlat-inlat).argmin()]
        
        if np.abs(lonpt-inlon)>self.gridres or np.abs(latpt-inlat)>self.gridres:
            return print("Error: Point not found on grid")
        else:
            idx = pd.IndexSlice
            rawpoint = self.raw[0].loc[latpt, idx[:,:,lonpt]]
            
            hgtarr = np.flip(np.array(rawpoint.loc[idx[:,'gh']]), axis=0)###/9.81
            tmparr = np.flip(np.array(rawpoint.loc[idx[:,'t']]), axis=0)-degCtoK
            rharr = np.flip(np.array(rawpoint.loc[idx[:,'r']]), axis=0)
            uarr = np.flip(np.array(rawpoint.loc[idx[:,'u']]), axis=0)*mstokts
            varr = np.flip(np.array(rawpoint.loc[idx[:,'v']]), axis=0)*mstokts
            warr = np.flip(np.array(rawpoint.loc[idx[:,'w']]), axis=0)
            presarr = np.flip(np.array(rawpoint.index.get_level_values(0).unique()), axis=0)
            
            soundlist = []
            for i, ivar in enumerate(soundvars):
                if ivar=='HGT':
                    soundlist.append(hgtarr)
                elif ivar=='TMP':
                    soundlist.append(tmparr)
                elif ivar=='DWPT':
                    soundlist.append(calc_dwpt(tmparr, rharr) )
                elif ivar=='WDIR':
                    soundlist.append(calc_wdir(uarr, varr) )
                elif ivar=='WSPD':
                    soundlist.append(calc_wspd(uarr, varr) )
                elif ivar=='RH':
                    soundlist.append(rharr )
                elif ivar=='PRES':
                    soundlist.append(presarr )
            return soundlist
        
    def raw2block(self, *args):
        sharpargs = [x for x in args if x in features]
        customargs = [x for x in args if x in customfeats]
        junk = [x for x in args if x not in features + customfeats]
        
        if len(junk)>0:
            print("Ignoring unknown variables: ", junk)
                
        if self.blocks.shape == (0, 0):
            self.blocks = self.dfstruct([self.raw[1]], [(sharpargs + customargs)[0]] )
            
        if len(sharpargs)>0:
            sharparr = np.zeros((len(self.raw[0].index), len(sharpargs) ) )
        
            for i, ilon in enumerate(self.raw[0].columns.get_level_values(2).unique() ):
                for j, jlat in enumerate(self.raw[0].index):
                    sounddict = dict(zip(soundvars, self.raw2sound(ilon, jlat)) )
                
                    prof = sharppyprof(sounddict)
                    sharpvals = sharpcalc(sharpargs, prof)
                    sharparr[j, :] = np.array(sharpvals)
        
                for j, jvar in enumerate(sharpargs):
                    self.blocks[self.raw[1], jvar, ilon] = sharparr[:, j]
                if int(i%2) == 0:
                    print(ilon, " finished")
                
        if 'UPHEL' in customargs:
            print("calculating updraft helicity...")
            
            levels = list(self.raw[0].columns.levels[0])
            nlev, nlon, nlat = len(levels), len(self.gridlon), len(self.gridlat)
            idx = pd.IndexSlice
            
            gharr = self.raw[0].loc[:, idx[:,'gh',:]].values.reshape(nlat, nlev, nlon).swapaxes(0,1)
            uarr = self.raw[0].loc[:, idx[:,'u',:]].values.reshape(nlat, nlev, nlon).swapaxes(0,1)
            varr = self.raw[0].loc[:, idx[:,'v',:]].values.reshape(nlat, nlev, nlon).swapaxes(0,1)
            warr = self.raw[0].loc[:, idx[:,'w',:]].values.reshape(nlat, nlev, nlon).swapaxes(0,1)
            
            uphelarr = udhelicity(gharr, uarr, varr, warr, self.gridres) 
            
            for i, ilon in enumerate(self.raw[0].columns.get_level_values(2).unique() ):
                self.blocks[self.raw[1], 'UPHEL', ilon] = uphelarr[:, i]
                
                if int(i%2) == 0:
                    print(ilon, " finished")
            
        self.blocks = self.blocks.sort_index(axis=1)
        ###INEFFICIENT!!!
        pass
    
    def savehdf(self, h5name, datype, varlist = [], overwrite = False):
        if datype not in ['3D', '2D']:
            sys.exit('Invalid data type! Must be 3D or 2D')

        names = list(self.blocks.columns.names)[:-1] if datype == '2D' \
        else list(self.raw[0].columns.names)[:-1]

        dts = list(self.blocks.columns.get_level_values(0).unique()) if datype == '2D' \
        else [self.raw[1]]

        levels = [] if datype == '2D' \
        else list(self.raw[0].columns.get_level_values(0).unique()) 

        allvars = list(self.blocks.columns.get_level_values(1).unique()) if datype == '2D' \
        else list(self.raw[0].columns.get_level_values(1).unique())

        params = [x for x in varlist if x in allvars] if len(varlist) > 0 else allvars
        junk = [x for x in varlist if x not in allvars]

        if len(junk) > 0:
            print ("Ignoring unknown variables: ", junk)

        with h5py.File(h5name, 'a') as f:

            if os.path.exists(h5name) and len(list(f.keys() ) ) > 0:
                print('File exists.. adding to existing file')

                boolrange = (f.attrs['gridrange'] == self.gridrange).all()
                boolres = f.attrs['gridres'] == self.gridres

                if not boolrange or not boolres:
                    print(self.gridrange, self.gridres)
                    print("Should be: ")
                    print(f.attrs['gridrange'], f.attrs['gridres'])
                    sys.exit("Inconsistent grid sizes")

            else:
                f.attrs['gridrange'] = self.gridrange
                f.attrs['gridres'] = self.gridres

            for i, idt in enumerate(dts):
                if str(idt) not in f.keys():
                    g = f.create_group(str(idt))
                else:
                    g = f[str(idt)]

                for j, jparam in enumerate(params):
                    if jparam not in f[str(idt)].keys():
                        gg = g.create_group(jparam)
                    else:
                        gg = g[jparam]

                    if len(levels) == 0:
                        if jparam in gg.keys():
                            if overwrite:
                                gg[jparam][:] = self.blocks[idt, jparam]
                            else:
                                print('Will not overwrite data: ', str(idt), jparam)
                        else:
                            gg.create_dataset(jparam, data = self.blocks[idt, jparam])

                    else:
                        if 'levels' not in gg.keys():
                            ggg = gg.create_group('levels')
                        else:
                            ggg = gg['levels']

                        for k, klevel in enumerate(levels):
                            if str(klevel) in ggg.keys():
                                if overwrite:
                                    ggg[str(klevel)][:] = self.raw[0][klevel, jparam]
                                else:
                                    print('Will not overwrite data: ', str(idt), str(klevel))

                            else:
                                ggg.create_dataset(str(klevel), data = \
                                                   self.raw[0][klevel, jparam])
                                
    def loadhdf(self, h5name, datype, varlist = []):
        if datype not in ['3D', '2D']:
            sys.exit('Invalid data type! Must be 3D or 2D')

        with h5py.File(h5name, 'r') as f:
            if len(f.keys()) == 0:
                sys.exit("Empty file")

            boolrange = (f.attrs['gridrange'] == self.gridrange).all()
            boolres = f.attrs['gridres'] == self.gridres

            keydts = list(f.keys())
            keyparams2D = []
            keyparams3D = []
            keylevels3D = []

            for i, idt in enumerate(keydts ):
                g = f[idt]
                params2D = []
                params3D = []

                for j, jparam in enumerate(list(g.keys() ) ):
                    gg = g[jparam]
                    if jparam in gg.keys():
                        params2D.append(jparam)
                    if 'levels' in list(gg.keys()):
                        ggg = gg['levels']
                        params3D.append(jparam)

                        if len(keylevels3D) == 0:
                            keylevels3D = list(ggg.keys())
                        else:
                            if len(ggg.keys()) != len(keylevels3D):
                                print("Warning: Different variables contain different \
                                number of levels")

                keyparams2D.append(params2D)
                keyparams3D.append(params3D)

                print("Contents for datetime: ", idt)
                print("2D parameters: ", params2D)
                print("3D parameters: ", params3D)

            if not boolrange or not boolres:
                if self.blocks.shape == (0, 0):
                    print("Safe to adjust grid size (empty block):")
                    print(f.attrs['gridrange'], f.attrs['gridres'])

                    self.gridrange = f.attrs['gridrange']
                    self.gridres = f.attrs['gridres']
                    self.gridlon, self.gridlat = std_grid(self.gridrange, 
                                                                    self.gridres)
                else:
                    print("Incompatible grid dimensions:")
                    print(f.attrs['gridrange'], f.attrs['gridres'])
                    print("Grid range / res must be:")
                    print(self.gridrange, self.gridres)
                    sys.exit()

            if datype == '3D':
                if len(keydts) > 0:
                    print("Multiple datetimes: ", keydts)
                    print("Loading datetime: ", keydts[0])   
                params = keyparams3D[0]
                levels = [int(x) for x in keylevels3D]
                g = f[keydts[0]]

                self.raw[0] = self.dfstruct(levels, params)
                self.raw[1] = parser.parse(keydts[0])

                for j, jparam in enumerate(params):
                    gg = g[jparam]

                    for k, klevel in enumerate(levels):
                        data = gg['levels/' + str(klevel)]
                        self.raw[0][klevel, jparam] = data[:]

            if datype == '2D':
                params = keyparams2D[0] 

                if len(varlist) > 0:
                    params = [x for x in keyparams2D[0] if x in varlist]
                    junk = [x for x in keyparams2D[0] if x not in varlist]
                    if len(junk) > 0:
                        print("Ignoring unknown variables:", junk)

                dts = [parser.parse(x) for x in keydts]

                print(list(self.blocks.columns.get_level_values(0).unique() ) )
                
                

                for i, idt in enumerate(dts):
                    print("Loading parameters: ", params, idt)
                    for j, jparam in enumerate(params):
                        data = f[str(idt) + '/' + jparam + '/' + jparam]
                        
                        if self.blocks.shape == (0,0):
                            self.blocks = self.dfstruct(dts, params)
                            self.blocks[idt, jparam] = data[:]
                        else:
                            lons = self.blocks.columns.get_level_values(2).unique()
                            
                            for k, klon in enumerate(lons):
                                self.blocks[idt, jparam, klon] = data[:, k]
    
    def plot2fig(self, step, var, fig, ax, coord = [], title = ''):
        rawlist = ['gh', 't', 'r', 'u', 'v', 'w']
        gridmesh = np.meshgrid(self.gridlon, self.gridlat)
        
        plotmesh = ax.pcolormesh(gridmesh[0], gridmesh[1], self.raw[0][step, var]) if var in rawlist \
        else ax.pcolormesh(gridmesh[0], gridmesh[1], self.blocks[step, var])
        
        plotbar = fig.colorbar(plotmesh, ax = ax)

        if len(coord)>0:
            for ilabel, jll in coord:
                ax.plot([jll[0]], [jll[1]], marker='o', markersize=6, color="red")
                ax.text(jll[0] + 0.3, jll[1] + 0.3, ilabel, fontsize=10, color="white")
                
        pad = 5 # in points
        
        if len(title) > 0:
            ax.annotate(title, xy=(0.5, 0.), xytext=(0, pad),
                xycoords=ax.title, textcoords='offset points',
                size='large', ha='center', va='baseline');
        
        return plotmesh, plotbar
    

