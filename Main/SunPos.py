import datetime

import astroplan
import numpy as np
import astropy.units
import astropy.time

import SchedulerConsts

def make_APF_obs():
    '''
    apf_obs = make_APF_obs(dt, horizon=str(TARGET_ELEVATION_MIN))
    dt - datetime object
    horizon - string of horizon in degrees
    apf_obs - returns astroplan.Observer object for the time dt with the horizon set to horizon
    '''
    # Generate an astroplan observer for the APF
    lat  = '37:20:33.1'
    lon = '-121:38:17.7'
    elevation = 1274 * astropy.units.m
    pressure = 0.870 * astropy.units.bar # typical value

    apf_obs = astroplan.Observer(latitude=lat, longitude=lon, elevation=elevation,
                                 pressure=pressure, name='APF', timezone='US/Pacific')

    return apf_obs

def compute_sunset_rise(dt, horizon='0'):
    '''
    sunset, sunrise = compute_sunset_rise(dt, horizon='0')
    dt - datetime object
    horizon - string of horizon in degrees
    computes time in seconds before sunset and next sunrise from dt
    '''
    apf_obs = make_APF_obs()
    curr_time = astropy.time.Time(dt, format='datetime')
    sunset_time = apf_obs.sun_set_time(curr_time, horizon=float(horizon)*astropy.units.deg, which='next')
    sunrise_time = apf_obs.sun_rise_time(curr_time, horizon=float(horizon)*astropy.units.deg, which='next')
    sunset = float(sunset_time.unix - curr_time.unix)
    sunrise = float(sunrise_time.unix - curr_time.unix)
    return sunset, sunrise

def compute_sunset(dt, horizon='0'):
    '''
    sunset = compute_sunset(dt, horizon='0')
    dt - datetime object
    horizon - string of horizon in degrees
    helper to compute just sunset, calls compute_sunset_rise
    '''
    sunset, _ = compute_sunset_rise(dt, horizon=horizon)
    return sunset

def compute_sunrise(dt, horizon='0'):
    '''
    sunrise = compute_sunrise(dt, horizon='0')
    dt - datetime object
    horizon - string of horizon in degrees
    helper to compute just sunrise, calls compute_sunset_rise
    '''
    _, sunrise = compute_sunset_rise(dt, horizon=horizon)
    return sunrise

def sun_el_check(star_table, apf_obs, horizon='-18'):
    '''
    sun_el_check = sun_el_check(star_table, stars, idx, apf_obs, dt, horizon='0')
    star_table - astropy table of targets
    stars - list of astroplan.FixedTarget objects
    idx - index of target in star_table
    apf_obs - astroplan.Observer object
    dt - datetime object
    horizon - string of horizon in degrees
    sun_el_check - boolean
    '''
    bright_enough = np.ones(len(star_table['Vmag']), dtype=bool)

    sun = apf_obs.sun_altaz(datetime.datetime.utcnow())
    sun_el = float(sun.alt.value)

    faint = star_table['Vmag'] > SchedulerConsts.SLOWDOWN_VMAG_LIM
    faint &= star_table['too'] is False

    if sun_el > float(horizon):
        bright_enough[faint] = False

    return bright_enough
