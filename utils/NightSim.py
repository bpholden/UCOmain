import sys
import re
sys.path.append("../Main")
#from ExposureCalc import *

import numpy as np
import astroplan
import astropy.units
import astropy.time
import ExposureCalculations
import Generate_Errors 

def sun_times(dt_obj):
    '''
    Compute the sunset and sunrise times for a given date and observer location.

    Parameters:
    dt_obj : datetime.datetime
        The datetime object representing the date for which to compute sunset and sunrise times.

    Returns:
    tuple
        A tuple containing the sunset time, sunrise time, and the observer object.
    '''
    lat  = '37:20:33.1'
    lon = '-121:38:17.7'

    # Minimum observation to observe things at
    horizon = -9.0*astropy.units.deg
    cur_date = astropy.time.Time(dt_obj, format='datetime', scale='utc')

    elevation = 1274 * astropy.units.m
    pressure = 0.870 * astropy.units.bar # typical value

    apf_obs = astroplan.Observer(latitude=lat, longitude=lon, elevation=elevation,
                                 pressure=pressure, name='APF', timezone='US/Pacific')

    sunset_time = apf_obs.sun_set_time(cur_date, horizon=horizon, which='next')
    sunrise_time = apf_obs.sun_rise_time(cur_date, horizon=horizon, which='next')

    return sunset_time, sunrise_time, apf_obs

def make_obs_sample(fn):
    '''
    Load observational sample data from a file.

    Parameters:
    fn : str
        The filename containing the observational sample data.

    Returns:
    tuple
        A tuple containing the slowdown and FWHM arrays.
    '''
    slow,fwhm = np.loadtxt(fn,unpack=True)
    return slow, fwhm

def gen_seeing(nsize=200,val=-1):
    '''
    Generate a sequence of seeing deviation values based on atmospheric conditions.

    Parameters:
    nsize : int, optional
        Number of values to generate (default is 200).
    val : float, optional
        Control parameter for seeing conditions (default is -1, which triggers random selection).

    Returns:
    numpy.ndarray
        Array of generated seeing deviation values.
    '''
    if val < 0:
        val = np.random.uniform(size=1)
    alpha = 0.52

    if val < 0.9:
        mean = np.random.normal(loc=8.,scale=1.0,size=1)
        rms  = np.random.normal(loc=1.9,scale=1.0,size=1)
    else:
        mean = np.random.normal(loc=19.0,scale=1.0,size=1)
        rms  = np.random.normal(loc=4.5,scale=1.0,size=1)
    real_rms = np.sqrt((1-alpha**2) * rms**2)
    deviates = np.random.normal(loc=0,scale=real_rms,size=nsize)
    for i in np.arange(1,nsize):
        deviates[i] += alpha*deviates[i-1]
    deviates += mean

#    deviates = np.random.normal(loc=mean,scale=rms,size=nsize)
    return deviates

def gen_seeing_el(deviate,el):
    '''
    Adjust the seeing deviation based on the elevation angle.

    Parameters:
    deviate : float
        The initial seeing deviation.
    el : astropy.units.Quantity
        Elevation angle in degrees.

    Returns:
    float
        The adjusted seeing deviation.
    '''
    zd = 90 - el
    deviate += (0.0903544076597*zd +  -0.00172591889888*zd*zd + 3.3157238117e-05*zd*zd*zd)
    return deviate

def gen_clouds(nsize=200,val=-1):
    '''
    Generate a sequence of cloud slowdown values based on atmospheric conditions.

    Parameters:
    nsize : int, optional
        Number of values to generate (default is 200).
    val : float, optional
        Control parameter for cloud conditions (default is -1, which triggers random selection).

    Returns:
    numpy.ndarray
        Array of generated cloud slowdown values.
    '''
    if val < 0:
        val = np.random.uniform(size=1)
    alpha = 0.353
    if val < 0.7:
        mean = np.random.normal(loc=1.,scale=0.1,size=1)
        rms  = np.random.normal(loc=0.05,scale=0.01,size=1)
    elif val > 0.9:
        mean = np.random.normal(loc=2.0,scale=0.4,size=1)
        rms  = np.random.normal(loc=0.5,scale=0.1,size=1)
    else:
        mean = np.random.normal(loc=1.7,scale=0.1,size=1)
        rms  = np.random.normal(loc=0.3,scale=0.05,size=1)

    real_rms = np.sqrt((1-alpha**2) * rms**2)
    deviates = np.random.normal(loc=0,scale=real_rms,size=nsize)
    for i in np.arange(1,nsize):
        deviates[i] += alpha*deviates[i-1]
    deviates += mean
    deviates[deviates < 0.3] = 0.3
    return deviates

def rand_obs_sample(slows,fwhms):
    '''
    Randomly sample a slowdown and FWHM value from the provided arrays.

    Parameters:
    slows : array-like
        Array of slowdown values.
    fwhms : array-like
        Array of FWHM values.

    Returns:
    tuple
        A tuple containing a randomly selected slowdown and FWHM value.
    '''
    ls = len(slows) -1
    lf = len(fwhms) -1
    sindx = np.random.randint(0,ls)
    findx = np.random.randint(0,lf)
    return slows[sindx], fwhms[findx]

def compute_el(curtime,star,apf_obs):
    '''
    Compute the altitude and azimuth of a star at a given time and observer location.

    Parameters:
    curtime : float
        Current time in Julian date format.
    star : astropy.coordinates.FixedBody
        The star object.
    apf_obs : astropy.coordinates.EarthLocation
        The observer object.

    Returns:
    actel : float
        Altitude of the star in degrees.
    actaz : float
        Azimuth of the star in degrees.
    '''
    star_pos = apf_obs.altaz(curtime, star)
    return star_pos.alt.value, star_pos.az.value


def check_date(datestr):
    '''
    Check if a date string is in the format YYYY/MM/DD or YYYY-MM-DD and represents a valid date.

    Parameters:
    datestr : str
        The date string to check.

    Returns:
    bool
        True if the date string is valid, False otherwise.
    '''
    match = re.match(r"(\d{4})(\/|\-)(\d{1,2})(\/|\-)(\d{1,2})",datestr)

    if not match:
        return False
    if int(match.group(3)) < 1 or int(match.group(3)) > 12:
        return False
    if int(match.group(5)) < 1 or int(match.group(5)) > 31:
        return False

    return True


def compute_simulation(result, curtime, star, apf_obs, slowdowns, fwhms, owner):
    '''
    Compute the simulation for a given observation result.

    Parameters:
    result : dict
        Dictionary containing observation result information.
    curtime : float
        Current time in Julian date format.
    star : astropy.coordinates.FixedBody
        The star object.
    apf_obs : astroplan.Observer
        The observer object.
    slowdowns : array-like
        Array of slowdown values.
    fwhms : array-like
        Array of FWHM values.
    owner : str
        Owner of the observation.

    Returns:
    tuple
        A tuple containing the updated current time, last FWHM, last slowdown, and output string.
    '''
    actel, actaz = compute_el(curtime, star, apf_obs)
    actslow, actfwhm = rand_obs_sample(slowdowns, fwhms)
    actfwhm = gen_seeing_el(actfwhm, actel)
    lastfwhm = actfwhm
    lastslow = actslow
    meterrate = ExposureCalculations.getEXPMeter_Rate(result['VMAG'],
                                                      result['BV'],actel,actfwhm,result['DECKER'])
    meterrate *= 1 + 0.11*np.random.randn(1)[0]
    meterrate /= actslow
    specrate = ExposureCalculations.getSpec_Rate(result['VMAG'],
                                                 result['BV'],actel,actfwhm,result['DECKER'])
    specrate *= 1 + 0.11*np.random.randn(1)[0]
    specrate /= actslow
    metertime = result['COUNTS'] / meterrate
    exp_time = result['EXP_TIME']
    barycentertime = curtime
    if metertime < exp_time:
        fexptime = float(metertime)
    else:
        fexptime = float(exp_time)

    curtime += astropy.time.TimeDelta(fexptime+40., format='sec')
    barycentertime += astropy.time.TimeDelta(fexptime/2., format='sec')
    totcounts = fexptime * specrate

    #precision, true_error = Generate_Errors.compute_real_uncertainty(totcounts,result['BV'])
    if actaz < 180:
        actel *= -1.
    outstr = "%s %s %.5f %.1f %.1f %.2f %.2f %.2f %.2f %s" % \
    (result['NAME'] , curtime.isot, barycentertime.jd,\
      fexptime, totcounts, actel, actaz, actfwhm, actslow, owner)
    print (outstr)

    return curtime, lastfwhm, lastslow, outstr



def init_sim_vals():
    '''
    Initialize simulation values.

    Returns:
    tuple
        A tuple containing a dictionary of simulation values and a list of simulation names.
    '''
    vals = dict()

    sim_names = ('name','date','time','jd','etime','I2','El','Az','fwhm','slow','owner')
    for n in sim_names:
        vals[n] = []
    return vals, sim_names


def read_sim_lines(lines,sim_names,sim_vals):
    '''
    Read simulation lines and populate simulation values.

    Parameters:
    lines : list of str
        List of lines from the simulation output file.
    sim_names : list of str
        List of simulation parameter names.
    sim_vals : dict
        Dictionary to store simulation values.

    Returns:
    None
    '''
    for l in lines:
        if re.search(r"\A\#",l):
            continue
        d = l.split()
        if len(d) == 11:
            for i, sim_name in enumerate(sim_names):
                sim_vals[sim_name].append(d[i])

    return

def sum_owner_times(vals):
    '''
    Sum the observation times, elevations, and number of exposures for each owner.

    Parameters:
    vals : dict
        Dictionary containing simulation values.

    Returns:
    tuple
        A tuple containing dictionaries for total times, average elevations, and number of exposures per owner.
    '''

    owner_tots = dict()
    owner_els = dict()
    owner_nexps = dict()
    for o in vals['owner']:
        if o not in owner_tots:
            owner_tots[o] = 0.
            owner_els[o] = 0.
            owner_nexps[o] = 0
    owner_tots['total'] = 0.
    owner_els['total'] = 0.
    owner_nexps['total'] = 0

    for i in range(0,len(vals['owner'])):
        owner_tots[vals['owner'][i]] += float(vals['etime'][i])
        owner_els[vals['owner'][i]] += float(vals['El'][i])
        owner_nexps[vals['owner'][i]] += 1
        owner_tots['total'] += float(vals['etime'][i])
        owner_els['total'] += float(vals['El'][i])
        owner_nexps['total'] += 1

    for o in owner_els:
        if owner_nexps[o] > 0:
            owner_els[o] /= owner_nexps[o]

    return owner_tots, owner_els, owner_nexps
