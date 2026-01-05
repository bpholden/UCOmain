import ephem
import numpy as np

import SchedulerConsts

def make_APF_obs(dt, horizon=str(SchedulerConsts.TARGET_ELEVATION_MIN)):
    '''
    apf_obs = make_APF_obs(dt, horizon=str(TARGET_ELEVATION_MIN))
    dt - datetime object
    horizon - string of horizon in degrees
    apf_obs - returns ephem.Observer object for the time dt with the horizon set to horizon
    '''
    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = horizon
    apf_obs.date = dt

    return apf_obs

def compute_sunset_rise(dt, horizon='0'):
    '''
    sunset, sunrise = compute_sunset_rise(dt, horizon='0')
    dt - datetime object
    horizon - string of horizon in degrees
    computes time in seconds before sunset and next sunrise from dt
    '''
    apf_obs = make_APF_obs(dt, horizon=horizon)
    sunset = apf_obs.next_setting(ephem.Sun())
    sunset -= ephem.Date(dt)
    sunset *= 86400.0 # convert to seconds

    sunrise = apf_obs.next_rising(ephem.Sun())
    sunrise -= ephem.Date(dt)
    sunrise *= 86400.0 # convert to seconds
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
    stars - list of ephem.FixedBody objects
    idx - index of target in star_table
    apf_obs - ephem.Observer object
    dt - datetime object
    horizon - string of horizon in degrees
    sun_el_check - boolean
    '''
    bright_enough = np.ones(len(star_table['Vmag']), dtype=bool)

    sun = ephem.Sun()
    sun.compute(apf_obs)
    sun_el = np.degrees(sun.alt)

    faint = star_table['Vmag'] > SchedulerConsts.SLOWDOWN_VMAG_LIM

    if sun_el > float(horizon):
        bright_enough[faint] = False

    return bright_enough
