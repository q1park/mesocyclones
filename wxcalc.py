import numpy as np

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
    
degCtoK = 273.15
mstokts = 1.94384

features = ['SBCAPE', 'SBCIN', 'SBLCL', 'SBLFC', 'SBEL', 'SBLI',
            'MLCAPE', 'MLCIN', 'MLLCL', 'MLLFC', 'MLEL', 'MLLI', 
            'MUCAPE', 'MUCIN', 'MULCL', 'MULFC', 'MUEL', 'MULI', 
            '0-1 km SRH', '0-1 km Shear', '0-3 km SRH', 'Eff. SRH', 
            'EBWD', 'PWV', 'K-index', 
            'STP(fix)', 'SHIP', 'SCP', 'STP(cin)']

customfeats = ['UPHEL']

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

def sharppyprof(dictsound):
    profsharp = profile.create_profile(profile='default', 
                                       pres=dictsound.get('PRES'),
                                       hght=dictsound.get('HGT'),
                                       tmpc=dictsound.get('TMP'),
                                       dwpc=dictsound.get('DWPT'),
                                       wspd=dictsound.get('WSPD'),
                                       wdir=dictsound.get('WDIR'),
                                       missing=-9999.00, strictQC=False)
    return profsharp
    
def sharpcalc(feats, profsharp):
    sfcpcl = params.parcelx( profsharp, flag=1 ) # Surface Parcel
    mupcl = params.parcelx( profsharp, flag=3 ) # Most-Unstable Parcel
    mlpcl = params.parcelx( profsharp, flag=4 ) # 100 mb Mean Layer Parcel
    
    sfc = profsharp.pres[profsharp.sfc]
    p3km = interp.pres(profsharp, interp.to_msl(profsharp, 3000.))
    p6km = interp.pres(profsharp, interp.to_msl(profsharp, 6000.))
    p1km = interp.pres(profsharp, interp.to_msl(profsharp, 1000.))
    mean_3km = winds.mean_wind(profsharp, pbot=sfc, ptop=p3km)
    sfc_6km_shear = winds.wind_shear(profsharp, pbot=sfc, ptop=p6km)
    sfc_3km_shear = winds.wind_shear(profsharp, pbot=sfc, ptop=p3km)
    sfc_1km_shear = winds.wind_shear(profsharp, pbot=sfc, ptop=p1km)
    
    srwind = params.bunkers_storm_motion(profsharp)
    srh3km = winds.helicity(profsharp, 0, 3000., stu = srwind[0], stv = srwind[1])
    srh1km = winds.helicity(profsharp, 0, 1000., stu = srwind[0], stv = srwind[1])

    ship = params.ship(profsharp)
    eff_inflow = params.effective_inflow_layer(profsharp)
    ebot_hght = interp.to_agl(profsharp, interp.hght(profsharp, eff_inflow[0]))
    etop_hght = interp.to_agl(profsharp, interp.hght(profsharp, eff_inflow[1]))
    effective_srh = winds.helicity(profsharp, ebot_hght, etop_hght, stu = srwind[0], stv = srwind[1])
    ebwd = winds.wind_shear(profsharp, pbot=eff_inflow[0], ptop=eff_inflow[1])
    ebwspd = utils.mag( ebwd[0], ebwd[1] )
    
    scp = params.scp(mupcl.bplus, effective_srh[0], ebwspd)
    stp_cin = params.stp_cin(mlpcl.bplus, effective_srh[0], ebwspd, mlpcl.lclhght, mlpcl.bminus)
    stp_fixed = params.stp_fixed(sfcpcl.bplus, sfcpcl.lclhght, srh1km[0], 
                             utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1])

    scp = 0.0 if np.isnan(float(scp) ) else scp
    stp_cin = 0.0 if np.isnan(float(stp_cin) ) else stp_cin
    stp_fixed = 0.0 if np.isnan(float(stp_fixed) ) else stp_fixed
    
    dictfeats = {'SBCAPE': [sfcpcl.bplus, 'J/kg'],\
                 'SBCIN': [sfcpcl.bminus, 'J/kg'],\
                 'SBLCL': [sfcpcl.lclhght, 'm AGL'],\
                 'SBLFC': [sfcpcl.lfchght, 'm AGL'],\
                 'SBEL': [sfcpcl.elhght, 'm AGL'],\
                 'SBLI': [sfcpcl.li5, 'C'],\
                 'MLCAPE': [mlpcl.bplus, 'J/kg'],\
                 'MLCIN': [mlpcl.bminus, 'J/kg'],\
                 'MLLCL': [mlpcl.lclhght, 'm AGL'],\
                 'MLLFC': [mlpcl.lfchght, 'm AGL'],\
                 'MLEL': [mlpcl.elhght, 'm AGL'],\
                 'MLLI': [mlpcl.li5, 'C'],\
                 'MUCAPE': [mupcl.bplus, 'J/kg'],\
                 'MUCIN': [mupcl.bminus, 'J/kg'],\
                 'MULCL': [mupcl.lclhght, 'm AGL'],\
                 'MULFC': [mupcl.lfchght, 'm AGL'],\
                 'MUEL': [mupcl.elhght, 'm AGL'],\
                 'MULI': [mupcl.li5, 'C'],\
                 '0-1 km SRH': [srh1km[0], 'm2/s2'],\
                 '0-1 km Shear': [utils.comp2vec(sfc_1km_shear[0], sfc_1km_shear[1])[1], 'kts'],\
                 '0-3 km SRH': [srh3km[0], 'm2/s2'],\
                 'Eff. SRH': [effective_srh[0], 'm2/s2'],\
                 'EBWD': [ebwspd, 'kts'],\
                 'PWV': [params.precip_water(profsharp), 'inch'],\
                 'K-index': [params.k_index(profsharp), ''],\
                 'STP(fix)': [stp_fixed, ''],\
                 'SHIP': [ship, ''],\
                 'SCP': [scp, ''],\
                 'STP(cin)': [stp_cin, '']}
    featvals = []
    for ifeat in feats:
        featvals.append(dictfeats.get(ifeat)[0])
    
    return featvals

def calc_vertvort(u, v, dx):
    """
    Parameters 
    u (ms^-1)
    v (ms^-1)
    dx (m)
    
    Returns
    Vertical vorticity (s^-1) 
    """
    dvdx = np.gradient(v,dx,dx)[1] # calc dvdx
    dudy = np.gradient(u,dx,dx)[0] # calc dudy
    vertvort = dvdx - dudy
    return vertvort

from scipy.interpolate import interp1d

def interp_generic(level, coords, data):
    """
    Parameters
    level - int/float of level to interpolate too
    coords - 3d array of co-ordinates (height/pres, lat, lon)
    data - 3d array of data on coords
    
    Returns
    Value at interpolated level in 2d array
    """
    out = np.zeros((np.shape(data)[1],np.shape(data)[2])) # create array to hold interpolated data
    for i in range(np.shape(data)[1]):
        for j in range(np.shape(data)[2]):
            f = interp1d(coords[:,i,j], data[:,i,j], kind='linear', fill_value=np.nan, bounds_error=False)
            out[i,j] = f(level) # do interpolation
    return out

def udhelicity(heights, U, V, W, dll):
    nhgt, dhgt = 6, 1000
    dz = dhgt / (nhgt - 1)
    
    u2km = np.zeros((nhgt, U.shape[1], U.shape[2]))
    v2km = np.zeros((nhgt, V.shape[1], V.shape[2]))
    u3km = np.zeros((nhgt, U.shape[1], U.shape[2]))
    v3km = np.zeros((nhgt, V.shape[1], V.shape[2]))
    u4km = np.zeros((nhgt, U.shape[1], U.shape[2]))
    v4km = np.zeros((nhgt, V.shape[1], V.shape[2]))
    #u5km = np.zeros((nhgt, U.shape[1], U.shape[2]))
    #v5km = np.zeros((nhgt, V.shape[1], V.shape[2]))
    w2km = np.zeros((nhgt, W.shape[1], W.shape[2]))
    w3km = np.zeros((nhgt, W.shape[1], W.shape[2]))
    w4km = np.zeros((nhgt, W.shape[1], W.shape[2]))
    #w5km = np.zeros((nhgt, W.shape[1], W.shape[2]))
    zeta2km = np.zeros((nhgt, U.shape[1], U.shape[2]))
    zeta3km = np.zeros((nhgt, U.shape[1], U.shape[2]))
    zeta4km = np.zeros((nhgt, U.shape[1], U.shape[2]))
    
    for i in range(nhgt): # loop through to interpolate to levels and store in array
        print("Interpolating...doing loop ", i, "of ", (nhgt-1))
        increment = i*dz
        u2km[i] = interp_generic(2000+increment, heights, U)
        v2km[i] = interp_generic(2000+increment, heights, V)
        u3km[i] = interp_generic(3000+increment, heights, U)
        v3km[i] = interp_generic(3000+increment, heights, V)
        u4km[i] = interp_generic(4000+increment, heights, U)
        v4km[i] = interp_generic(4000+increment, heights, V)
        #u5km[i] = interp_generic(5000+increment, heights, U)
        #v5km[i] = interp_generic(5000+increment, heights, V)
        w2km[i] = interp_generic(2000+increment, heights, W)
        w3km[i] = interp_generic(3000+increment, heights, W)
        w4km[i] = interp_generic(4000+increment, heights, W)
        #w5km[i] = interp_generic(2000+increment, heights, W)
        zeta2km[i] = calc_vertvort(u2km[i], v2km[i], dll)
        zeta3km[i] = calc_vertvort(u3km[i], v3km[i], dll)
        zeta4km[i] = calc_vertvort(u4km[i], v4km[i], dll)
        #zeta5km[i] = calc_vertvort(u5km[i], v5km[i], dll)
        
    # calc the layer mean
    w2to3 = np.mean(w2km, axis=0)
    w3to4 = np.mean(w3km, axis=0)
    w4to5 = np.mean(w4km, axis=0)
    zeta2to3 = np.mean(zeta2km, axis=0)
    zeta3to4 = np.mean(zeta3km, axis=0)
    zeta4to5 = np.mean(zeta4km, axis=0)
    # calc the 2-5km UH
    return ( w2to3*zeta2to3 + w3to4*zeta3to4 + w4to5*zeta4to5 ) * 1000