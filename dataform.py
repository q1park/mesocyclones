import re
import os
import glob
import datetime
import subprocess
import pandas as pd

from datetime import datetime
from datetime import timedelta

import requests
import urllib3
urllib3.disable_warnings()
####CAREFUL!####

grb1to2 = '/home/q1park/grib2/grb1to2/grb1to2.pl'

def checkurl(url):
    request = requests.get(url)
    return request.status_code == 200

def str2dt(strdt):
    if len(strdt)==12:
        year, month, day, hour = int(strdt[:4]), int(strdt[4:6]), int(strdt[6:8]), int(strdt[8:10])
        return datetime(year, month, day, hour)
    else:
        print("Incorrect input format: 1999010123 and not: ", strdt)

def parsefile(filename):
    filetags = filename.split('.', 1)[0]
    fileexts = filename.split('.', 1)[-1]
    filetags = re.split('\_', filetags)
    
    h5type=[]
    if (len(filetags)==6):
        h5type = filetags[:1]
        filetags = filetags[1:]

    name, res = filetags[0], filetags[1]
    dt = str2dt(filetags[2] + filetags[3])
    fchr = int(filetags[4])
    return h5type + [name, res, dt, fchr, fileexts]

def grbname(name, res, dt, fcst, ext):
    filedt = "%s%s%s_%s" % (dt.year, str(dt.month).zfill(2), 
                             str(dt.day).zfill(2), str(dt.hour*100).zfill(4) )
    label = "%s_%s_%s_%s.%s" % (name, str(res), filedt, str(fcst).zfill(3), ext)
    return label

def grbconvert(datadir, files = 'all'):
    filetags = []
    fileexts = []
    
    if files=='all':
        filetags, fileexts = listdir(datadir, 'grb')
    elif isinstance(files,list):
        filetags = files
        fileexts = ['grb']*len(files)
        
    for i, val in enumerate(filetags):
        filename = "%s.%s" % (filetags[i], fileexts[i])
        
        if os.path.exists(datadir + '/' + filename + '.grb2'):
            print(val + " seems converted already... converting again")
            
        converter_cmd = "%s %s/%s" % (grb1to2, datadir, filename)
        subprocess.call(converter_cmd, shell=True)
        print("Converted file: " + val)

class gribmenu():
    def __init__(self):
        """Class to get grib files from NOAA site"""
        self.menu = pd.DataFrame(columns=['name', 'res', 'datetime', 'fchr', 'ext', 'url'])

    def createmenu(self, model, dtrange):
        """Prints list of available files for date specified"""
        urlbase = 'https://nomads.ncdc.noaa.gov/data/'
        modelurl = '%s%s/' % (urlbase, model)

        if not checkurl(modelurl):
            print('Invalid model? No website for: ' + modelurl) 
        
        idt, fdt = datetime.date(dtrange[0]), datetime.date(dtrange[1])  
        dtlist = [idt + timedelta(days=x) \
                  for x in range(0, (fdt-idt+timedelta(days=1)).days)]

       
        for dt in dtlist:
            strdt = dt.strftime('%Y%m%d')

            remd5 = ''
            if dt < datetime.date(datetime(2006, 1, 1) ):
                remd5 = '.*\/(.*)\n'
            else:
                remd5 = '\s(.*)\n'

            urlmenu = '%s%s/%s/%s.%s' % (modelurl, strdt[:-2], strdt, 'md5sum', strdt)
            if not checkurl(urlmenu):
                urlmenu = '%s%s/%s/%s.%s' % (modelurl, strdt[:-2], strdt, strdt, 'md5sum')
            if not checkurl(urlmenu):
                print('Invalid request or empty listing at: ' + urlmenu) 
            else:
                http = urllib3.PoolManager()
                getmenu = http.request('GET', urlmenu)
                menu = http.request('GET', urlmenu).data.decode('utf-8')
                menu = re.findall(remd5, menu)
                menu = [x.strip(' ') for x in menu]
                menu = [x for x in menu if "grb" in x]
                for i, val in enumerate(menu):
                    urlfile = '%s%s/%s/%s' % (modelurl, strdt[:-2], strdt, val)
                    menuitem = parsefile(val)
                    menuitem.append(str(urlfile))
                    self.menu.loc[len(self.menu)] = menuitem
    
    def dirmenu(self, datdir, ext):
        filelist = glob.glob(datdir+'/*.'+ext)
        filelist = [re.split('\\/', x)[::-1][0] for x in filelist]

        for i, val in enumerate(filelist):
            menuitem = parsefile(val)
            menuitem.append('N/A')
            self.menu.loc[len(self.menu)] = menuitem

    def filtermenu(self, **kwargs):
        if 'name' in kwargs:
            self.menu = self.menu.drop(self.menu[self.menu.name!=kwargs.get('name')].index)
        if 'res' in kwargs:
            self.menu = self.menu.drop(self.menu[self.menu.res!=kwargs.get('res')].index)
        if 'fchr' in kwargs:
            self.menu = self.menu.drop(self.menu[self.menu.fchr!=kwargs.get('fchr')].index)
        if 'start' in kwargs:
            self.menu = self.menu.drop(self.menu[self.menu.datetime<kwargs.get('start')].index)
        if 'end' in kwargs:
            self.menu = self.menu.drop(self.menu[self.menu.datetime>kwargs.get('end')].index)
        if 'ext' in kwargs:
            self.menu = self.menu.drop(self.menu[self.menu.ext>kwargs.get('ext')].index)                 
        self.menu = self.menu.reset_index(drop=True)
        pass
    
    def menuselect(self, **kwargs):
        selectdf = self.menu
        if 'name' in kwargs:
            selectdf = selectdf[(selectdf['name'] == kwargs.get('name'))]
        if 'res' in kwargs:
            selectdf = selectdf[(selectdf['res'] == kwargs.get('res'))]
        if 'fchr' in kwargs:
            selectdf = selectdf[(selectdf['fchr'] == kwargs.get('fchr'))]
        if 'datetime' in kwargs:
            selectdf = selectdf[(selectdf['datetime'] == kwargs.get('datetime'))]
        if 'ext' in kwargs:
            selectdf = selectdf[(selectdf['ext'] == kwargs.get('ext'))]
            
        selectdf = selectdf.reset_index(drop=True)
        if len(selectdf)==1:
            return grbname(*list(selectdf.loc[0][:5]))

        else:
            print("Number of matching files: ", len(selectdf))
            
    def checkonline(self):
        intcheck = 0
        for row, cols in self.menu.iterrows():
            argstag = [cols['name'], cols['res'], cols['datetime'], cols['fchr'], cols['ext']]
            filetag = grbname(*argstag)
            fileurl = cols['url']
            
            if not checkurl(fileurl):
                intcheck = intcheck + 1
                print("Missing from website: ", fileurl)
            if row%5==0:
                print("Good up to row: ", row)
        print(intcheck, " files missing from website")
        pass
    
    def checklocal(self, checkdir = './'):
        intcheck = 0
        for row, cols in self.menu.iterrows():
            argstag = [cols['name'], cols['res'], cols['datetime'], cols['fchr'], cols['ext']]
            filetag = grbname(*argstag)
            
            if not os.path.exists(checkdir + filetag):
                intcheck = intcheck + 1
                print("In ", checkdir, " missing file: ", filetag)
        print("Missing ", intcheck, " files total")
        pass
    
    def loadmenu(self, csvfile):
        self.menu = pd.read_csv(csvfile, index_col=False)
        self.menu['datetime'] = pd.to_datetime(self.menu['datetime'])
        self.menu['fchr'] = self.menu['fchr'].astype(int)
    
    def savemenu(self, menuname):
        self.menu.to_csv(menuname, index=False)
        