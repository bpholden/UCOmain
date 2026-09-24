# Observability.py
# pure array/astro filters on the star table, split out of UCOScheduler.py
# nothing here holds scheduler state
from __future__ import print_function

import datetime

import numpy as np
import ephem

import SchedulerConsts
import ScriptobsLine
import SunPos
import Visible

try:
    from apflog import apflog
except:
    from fake_apflog import apflog

def compute_datetime(ctime):
    '''
    dt = compute_datetime(ctime)
    ctime - can be a float, datetime, or ephem.Date, else UT now is used
    dt - datetime object appropriate for ctime.
    '''
    if isinstance(ctime, float):
        dt = datetime.datetime.utcfromtimestamp(int(ctime))
    elif isinstance(ctime, datetime.datetime):
        dt = ctime
    elif isinstance(ctime, ephem.Date):
        dt = ctime.datetime()
    else:
        #punt and use current UT
        dt = datetime.datetime.utcnow()
    return dt

def tot_exp_times(star_table, targ_num):
    '''
    totexptimes = tot_exp_times(star_table, targ_num)
    star_table - astropy table of targets
    targ_num - number of targets

    totexptimes - numpy array of total exposure times
    '''
    totexptimes = np.zeros(targ_num, dtype=float)

    nobs = np.ones(targ_num)
    multiples = (star_table['night_cad'] > 0)  & (star_table['night_obs'] == 0)
    nobs[multiples] = star_table['night_nexp'][multiples]

    totexptimes = nobs*(star_table['texp'] * star_table['nexp'] + 40 * (star_table['nexp']-1))
    totexptimes += (nobs-1)*star_table['night_cad']*86400

    return totexptimes

def time_check(star_table, totexptimes, dt, start_time=None):
    """ time_check = time_check(star_table, totexptimes, dt, hour_table)
    star_table - astropy table of targets
    totexptimes - numpy array of total exposure times
    dt - datetime object
    time_check - numpy array of booleans
    values are determined by whether or not the target can be observed in the time left
    """
    maxexptime = SunPos.compute_sunrise(dt,horizon='-9')
    maxfaintexptime = SunPos.compute_sunrise(dt,horizon='-18')
    if maxfaintexptime > maxexptime:
        maxfaintexptime = 0

    if start_time is not None:
        # dt is a UT datetime object, start_time is a UT time stamp
        # however strftime assumes that the dt is in local time
        # JFC, this is a mess
        utc_offset = datetime.datetime.utcnow() - datetime.datetime.now()
        curr_time = float(dt.strftime('%s')) - utc_offset.total_seconds()
        if curr_time < start_time:
            maxexptime = start_time - curr_time
            maxfaintexptime = start_time - curr_time

    if maxexptime < SchedulerConsts.TARGET_EXPOSURE_TIME_MIN:
        maxexptime = SchedulerConsts.TARGET_EXPOSURE_TIME_MIN
        # this will try a target in case we get lucky
        # bright stars often have longer than
        # necessary exposure times, relying on the
        # exposure meter, we will not make this modification for faint stars
        # we need to handle cases where night_cad > 0

    started_multiples = (star_table['night_cad'] > 0) & (star_table['night_obs'] == 1)
    if np.any(started_multiples):
        cadence_check = ephem.julian_date(dt) - star_table['lastobs']
        waiting = cadence_check < (star_table['night_cad'] - SchedulerConsts.BUFFER )
        if np.any(waiting):
            maxexptimes = (star_table['night_cad'] - cadence_check) * 86400
            maxfaintexptimes = (star_table['night_cad'] - cadence_check) * 86400
            try:
                maxexptime = np.min(maxexptimes[waiting & started_multiples]) + SchedulerConsts.BUFFERSEC - 180
            except ValueError:
                # this means we have double observations we are waiting for
                # but they are in the selection window
                # they will be selected in the priorities method
                # so we should use the usual maximum exposure time
                pass
            try:
                maxfaintexptime = np.min(maxfaintexptimes[waiting & started_multiples])
                maxfaintexptime += SchedulerConsts.BUFFERSEC - 180
            except ValueError:
                pass

    faint = star_table['Vmag'] > SchedulerConsts.SLOWDOWN_VMAG_LIM
    faint &= star_table['too'] is False
    time_good = totexptimes <= maxexptime
    time_good_faint = totexptimes <= maxfaintexptime

    time_good[faint] = time_good_faint[faint]

    return time_good

def condition_cuts(moon, seeing, slowdown, star_table):
    """ available = condition_cuts(moon, seeing, slowdown, star_table)

    Checks if columns are in the star_table, then cuts on those, returns a boolean numpy array

    available - Boolean numpy array of available targets

    moon - phase value from pyephem, ranges from 0 to 100 (a percentage)
    seeing - size in pixels
    transparency - magnitudes of extinction

    """

    available = np.ones(len(star_table['ra']), dtype=bool)

    if 'seeing' in star_table.colnames:
        available = (star_table['seeing']/0.109 > seeing) & available

    if 'moon' in star_table.colnames and float(moon.alt) > 0:
        available = (star_table['moon'] > moon.moon_phase) & available

    if 'transparency' in star_table.colnames:
        ext = 2.5 * np.log10(slowdown)
        available = (star_table['transparency'] > ext) & available


    return available

def behind_moon(moon,ras,decs):
    '''
    moon_check = behind_moon(moon,ras,decs)
    moon - pyephem moon object
    ras - numpy array of right ascensions in radians
    decs - numpy array of declinations in radians
    moon_check - numpy array of booleans, True if the target is too close to the moon
    '''
    md = SchedulerConsts.TARGET_MOON_DIST_MAX - SchedulerConsts.TARGET_MOON_DIST_MIN
    min_moon_dist = ((moon.phase / 100.) * md) + SchedulerConsts.TARGET_MOON_DIST_MIN
    moon_dist = np.arccos(np.cos(moon.dec) * np.cos(decs) * np.cos(moon.ra - ras)
                          + np.sin(moon.dec) * np.sin(decs)) # values in radians

    moon_check = np.degrees(moon_dist) > min_moon_dist

    return moon_check

def template_conditions(moon, seeing, slowdown):
    """ istrue = template_condition(moon, seeing, slowdown)

    Checks to see if moon, seeing and slowdown factor are within template conditions

    istrue - a simple Boolean

    moon - phase value from pyephem, ranges from 0 to 100 (a percentage)
    seeing - size in pixels
    slowdown - relative to clear

    """

    if seeing < SchedulerConsts.SEEING_TEMP and slowdown < SchedulerConsts.SLOWDOWN_TEMP:
        apflog("moon.phase=%.2f moon.alt=%.2f" % (moon.phase,moon.alt),echo=True,level='info')
        if moon.phase < 50 and float(moon.alt) < 0:
            return True
        if moon.phase < 25 and float(moon.alt) < 0.7:
            return True

    return False

def find_closest(ras, decs, ra, dec):
    '''
    find_closest(ras, decs, ra, dec)

    ras - numpy array of right ascensions in radians
    decs - numpy array of declinations in radians
    ra - right ascension in radians
    dec - declination in radians

    min_ind - index of the closest target

    searches for the closest target in ras, decs to ra, dec
    '''

    distances = np.sqrt((ra - ras)**2 + (dec - decs)**2)

    min_ind = distances.argmin()

    return min_ind

def find_Bstars(star_table,idx, bstars):
    '''
    find_Bstars(star_table,idx,bstars)

    star_table - astropy table of targets
    idx - index of target in star_table
    bstars - numpy array of booleans

    near_idx - index of the closest B star to template start time
    end_idx - index of the closest B star to template end time
    '''

    near_idx = find_closest(star_table['ra'][bstars], star_table['dec'][bstars],\
                            star_table['ra'][idx], star_table['dec'][idx])

    end_idx = find_closest(star_table['ra'][bstars], star_table['dec'][bstars],\
                            (star_table['ra'][idx]+15*np.pi/180.), star_table['dec'][idx])


    return near_idx,end_idx

def enough_time_templates(star_table, stars, idx, apf_obs, dt):
    '''
    enough_time_templates(star_table, stars, idx, apf_obs, dt)
    star_table - astropy table of targets
    stars - list of ephem.FixedBody objects
    idx - index of target in star_table
    apf_obs - ephem.Observer object
    dt - datetime object

    enough_time_templates - boolean

    Computes the time needed for a template observation
    and checks if there is enough time left before sunrise.
    '''

    count = ScriptobsLine.num_template_exp(star_table['Vmag'][idx])

    tot_time = count * 1200

    tot_time += 210 + (2*40 + 40*(star_table['nexp'][idx]-1)) + 2400 
    # two B star exposures + three 70 second acquisitions and the actual observation readout times
    vis, _, _ = Visible.visible(apf_obs, [stars[idx]], [tot_time])
    time_left_before_sunrise = SunPos.compute_sunrise(dt, horizon='-18')

    try:
        apflog("enough_time_templates(): time for obs= %.1f  time until sunrise= %.1f " % (tot_time, time_left_before_sunrise),echo=True)
    except:
        apflog("enough_time_templates(): cannot log times!?!",echo=True)

    if tot_time < time_left_before_sunrise  and vis and time_left_before_sunrise < 14*3600.:
        return True
    else:
        return False
