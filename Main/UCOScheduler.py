# UCOScheduler_V1.py
from __future__ import print_function
import os

import time
import datetime
import subprocess

import numpy as np
import astropy
import astropy.io
import astropy.table
import ephem
import ParseUCOSched
import SchedulerConsts

try:
    from apflog import *
    import ktl
except:
    from fake_apflog import *
import Visible

# a global
last_objs_attempted = []

# some constants
ACQUIRE = 'A'
BLANK = 'B'
FIRST = '1'
LAST = 'L'

BUFFERSEC = 600
BUFFER = BUFFERSEC / (24.*60*60)

def zero_last_objs_attempted():
    """
    zero_last_objs_attempted()

    Sets the global last_objs_attempted to an empty list.
    """
    global last_objs_attempted
    last_objs_attempted = []
    return

def need_cal_star(star_table, observed, priorities):
    """
    need_cal_star(star_table, priorities)

    Returns True if there is a calibration star in the star_table.
    """

    # this is kind of clunky
    if observed is None or observed.sheetns is None:
        return priorities

    cal_stars = np.zeros_like(priorities, dtype=bool)
    for sheetn in observed.sheetns:
        if np.any(star_table['need_cal'][star_table['sheetn'] == sheetn] == "Y"):
            # need to check if we need a cal, ie. the program that needs cals had targets
            # observed
            notdone = True
            cal_star_inds = (star_table['cal_star'] == 'Y') & (star_table['sheetn'] == sheetn)
            cal_star_names = star_table['name'][cal_star_inds]
            for cal_star_name in cal_star_names:
                if cal_star_name in observed.names:
                    notdone = False
            if notdone:
                cal_stars = cal_stars | cal_star_inds

    if np.any(cal_stars):
        priorities[cal_stars] = np.max(priorities) + 1

    return priorities

def compute_priorities(star_table, cur_dt, observed=None, hour_table=None, rank_table=None):
    """
    new_pri = compute_priorities(star_table, cur_dt,
                                    hour_table=None, rank_table=None)

    Computes the priorities for the targets in star_table.
    This is a function of the current time, the last time the target was observed,
    the cadence of the target, the current hour table and the rank table.
    """
    # make this a function, have it return the current priorities, than change
    # references to the star_table below into references to the current priority list
    new_pri = np.zeros_like(star_table['pri'])

    # new priorities will be
    new_pri[star_table['pri'] == 1] += 0
    new_pri[star_table['pri'] == 2] -= 20
    new_pri[star_table['pri'] == 3] -= 40

    cadence_check = ephem.julian_date(cur_dt) - star_table['lastobs']
    good_cadence = cadence_check > star_table['cad']
    bad_cadence = np.logical_not(good_cadence)

    started_doubles = star_table['night_cad'] > 0
    started_doubles = started_doubles & (star_table['night_obs'] > 0)
    started_doubles = started_doubles & (star_table['night_obs'] < star_table['night_nexp'])
    if np.any(started_doubles):
        redo = started_doubles & (cadence_check > (star_table['night_cad'] - BUFFER))
        redo = redo & (cadence_check < (star_table['night_cad'] + BUFFER))
    else:
        redo = np.zeros(1,dtype=bool)

    if hour_table is not None:
        too_much = hour_table['cur']  > hour_table['tot']
        done_sheets = hour_table['sheetn'][too_much]
    else:
        done_sheets = False

    if done_sheets is not False:
        done_sheets_str = " ".join(list(done_sheets))
        apflog("The following sheets are finished for the night: %s " %
               (done_sheets_str), echo=True)

    cadence_check /= star_table['cad']
    bad_pri = np.floor(cadence_check * 100)
    bad_pri = np.int_(bad_pri)

    done_all = star_table['nobs'] >= star_table['totobs']
    new_pri[done_all] = 0

    if rank_table is not None:
        for sheetn in rank_table['sheetn']:
            if sheetn not in done_sheets:
                cur = star_table['sheetn'] == sheetn
                new_pri[cur & good_cadence] += rank_table['rank'][rank_table['sheetn'] == sheetn]
                new_pri[cur & bad_cadence] += bad_pri[cur & bad_cadence]
            else:
                cur = star_table['sheetn'] == sheetn
                new_pri[cur & good_cadence] += 100
                new_pri[cur & bad_cadence] += bad_pri[cur & bad_cadence]

    if np.any(redo):
        new_pri[redo] = np.max(rank_table['rank'])

    new_pri = need_cal_star(star_table, observed, new_pri)

    return new_pri

def update_hour_table(hour_table, observed, dt, outfn='hour_table', outdir=None):
    '''
    update_hour_table(hour_table, observed, dt, outfn='hour_table', outdir=None)

    Updates hour_table with history of observations.
    observed is the observed log
    dt is the current datetime
    outfn is the output filename, defaults to hour_table
    outdir is the output directory, defaults to current working directory

    '''

    if not outdir :
        outdir = os.getcwd()

    outfn = os.path.join(outdir, outfn)

    hours = dict()

    # observed objects have lists as attributes
    # reverse time order, so most recent target observed is first.

    observed.reverse()

    nobj = len(observed.names)
    for i in range(0,nobj):
        own = observed.owners[i]
        if own not in list(hours):
            hours[own] = 0.0

    cur = dt
    for i in range(0,nobj):
        hr, mn = observed.times[i]
        prev = datetime.datetime(dt.year, dt.month, dt.day, hr, mn)
        diff = cur - prev
        hourdiff = diff.days * 24 + diff.seconds / 3600.
        if hourdiff > 0:
            hours[observed.owners[i]] += hourdiff
            cur = prev

    for ky in list(hours.keys()):
        if ky == 'public':
            hour_table['cur'][hour_table['sheetn'] == 'RECUR_A100'] = hours[ky]
        else:
            hour_table['cur'][hour_table['sheetn'] == ky] = hours[ky]

    try:
        hour_table.write(outfn,format='ascii',overwrite=True)
    except Exception as e:
        apflog("Cannot write table %s: %s %s" % (outfn, type(e), e), level='error', echo=True)

    observed.reverse()

    return hour_table


def make_hour_table(rank_table, dt, outfn='hour_table', outdir=None, hour_constraints=None):
    """

    hour_table = make_hour_table(rank_table, dt, outfn='hour_table', outdir=None, hour_constraints=None)

    Makes an hour table from the rank table and the current datetime.
    Writes it to outfn in outdir.
    The dt is a datetime object used to compute the length of the night.
    If hour_constraints is not None, it is a dictionary with keys 'runname' and 'left'
    which is checked against the default values in the hour table, and the final values
    are the lesser of the two.
    """

    if not outdir :
        outdir = os.getcwd()

    outfn = os.path.join(outdir,outfn)

    if os.path.exists(outfn):
        hour_table =  astropy.table.Table.read(outfn,format='ascii')
        return hour_table

    # file does not exist to make it from scratch using the fracs

    hour_table = astropy.table.Table([rank_table['sheetn'], \
                                      rank_table['frac']], names=['sheetn','frac'])

    sunset,sunrise = compute_sunset_rise(dt,horizon='-9')
    if sunrise < sunset:
        sunrise += 86400
    tot = sunrise - sunset
    tot /= 3600.

    hour_table['tot'] =np.abs(tot*hour_table['frac'])
    hour_table['cur'] =0.0*hour_table['frac']

    if hour_constraints is not None:
        if 'runname' in list(hour_constraints.keys()) and 'left' in list(hour_constraints.keys()):
            for runname in hour_constraints['runname']:
                if hour_constraints['left'][hour_constraints['runname']==runname] < hour_table['tot'][hour_table['sheetn']==runname]:
                    hour_table['tot'][hour_table['sheetn']==runname] = hour_constraints['left'][hour_constraints['runname']==runname]
                elif hour_constraints['left'][hour_constraints['runname']==runname] < 0:
                    hour_table['tot'][hour_table['sheetn']==runname] = -1.0

    try:
        hour_table.write(outfn,format='ascii')
    except Exception as e:
        apflog("Cannot write table %s: %s %s" % (outfn, type(e), e), level='error', echo=True)
    return hour_table

def find_time_left():
    """
    time_left = find_time_left()

    Uses the timereport/time_left command to find the time left each program has.
    Writes the output to a table and returns it.

    This is slow, so it should only be called once per night.

    """

    cmd = "/usr/local/lick/bin/timereport/time_left"
    if os.path.exists(cmd):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while p.poll() is None:
            time.sleep(1)
        out, err = p.communicate()
        if len(err):
            return None

        sheetns = []
        left = []
        alloc = []
        used = []
        lines = out.split('\n')
        if len(lines) <= 1:
            return None
        for ln in lines[1:]:
            d = ln.split(",")
            if len(d) >= 2:
                sheetns.append(d[0].strip())
                left.append(d[1].strip())
                alloc.append(d[2].strip())
                used.append(d[3].strip())

        rv = astropy.table.Table([sheetns,left,alloc,used], names=["runname","left","alloc","used"])

        return rv

    return None


def make_rank_table(sheet_table_name, outfn='rank_table', outdir=None, hour_constraints=None):
    """
    make_rank_table(sheet_table_name, outfn='rank_table', outdir=None, hour_constraints=None)

    Makes a rank table. The sheet_table_name is the name of the sheet which contains the
    current semester's rank table.
    The outfn is the output filename, defaults to rank_table, and the outdir is the output
    directory, defaults to the current working directory.

    If hour_constraints is not None, it is a dictionary with keys 'runname' and 'left'
    which is checked against the default values in the hour table, and the final values
    are the lesser of the two.
    If hour_constraints is None, calls find_time_left() to get the time left for each program.

    """
    if not outdir :
        outdir = os.getcwd()

    outfn = os.path.join(outdir,outfn)
    if os.path.exists(outfn):
        rank_table = astropy.table.Table.read(outfn,format='ascii')
        # the Booleans are now strings, so we have to convert back
        bs = [ True if sb == 'True' else False for sb in rank_table['too'] ]
        rank_table['too'] = bs
    else:
        sheetns, ranks, fracs, asciitoos = ParseUCOSched.parse_rank_table(sheet_table_name=sheet_table_name)
        if sheetns is None or len(sheetns) == 0:
            return None
            # this should result in this function being called again but with the
            # backup table being used
        toos = [ True if str(a) == 'y' else False for a in asciitoos ]

        rank_table= astropy.table.Table([sheetns,ranks,fracs,toos], names=['sheetn','rank','frac','too'])

        if hour_constraints:
            time_left = hour_constraints
        else:
            time_left = find_time_left()

        if time_left is not None:
            if 'runname' in list(time_left.keys()) and 'left' in list(time_left.keys()):
                for runname in time_left['runname']:
                    if float(time_left['left'][time_left['runname']==runname]) < 0:
                        rank_table['rank'][rank_table['sheetn']==runname] = -1000

        try:
            rank_table.write(outfn,format='ascii')
        except Exception as e:
            apflog("Cannot write table %s: %s %s" % (outfn, type(e), e), level='error', echo=True)

    return rank_table


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
    maxexptime = compute_sunrise(dt,horizon='-9')
    maxfaintexptime = compute_sunrise(dt,horizon='-18')
    if maxfaintexptime > maxexptime:
        maxfaintexptime = 0

    if start_time is not None:
        maxexptime = start_time - float(dt.strftime('%s'))
        maxfaintexptime = start_time - float(dt.strftime('%s'))

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
        waiting = cadence_check < (star_table['night_cad'] - BUFFER )
        if np.any(waiting):
            maxexptimes = (star_table['night_cad'] - cadence_check) * 86400
            maxfaintexptimes = (star_table['night_cad'] - cadence_check) * 86400
            try:
                maxexptime = np.min(maxexptimes[waiting & started_multiples]) + BUFFERSEC - 180
            except ValueError:
                # this means we have double observations we are waiting for
                # but they are in the selection window
                # they will be selected in the priorities method
                # so we should use the usual maximum exposure time
                pass
            try:
                maxfaintexptime = np.min(maxfaintexptimes[waiting & started_multiples])
                maxfaintexptime += BUFFERSEC - 180
            except ValueError:
                pass

    faint = star_table['Vmag'] > SchedulerConsts.SLOWDOWN_VMAG_LIM
    time_good = totexptimes <= maxexptime
    time_good_faint = totexptimes <= maxfaintexptime

    time_good[faint] = time_good_faint[faint]

    return time_good

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


def make_scriptobs_line(star_table_row, t, decker="W", I2="Y", owner='public', focval=0, coverid='', temp=False):
    """ given a name, a row in a star table and a do_flag, will generate
    a scriptobs line as a string:
    line = make_scriptobs_line(star_table_row, t, decker="W",I2="Y")

    star_table_row -contains all of the data needed for the line except
    t - a datetime object, this is used to fill in the uth and utm fields
    decker - one character field for the decker, defaults to "W"
    I2 - one character field for whether or not the Iodine cell is in, must be "Y" or "N"
    temp - a boolean for whether or not this is a template observation
    """

    # Add the RA as three elements, HR, MIN, SEC
    rastr = "%s %s %s " % (star_table_row['RA hr'],
                           star_table_row['RA min'],
                           star_table_row['RA sec'])

    # Add the DEC as three elements, DEG, MIN, SEC
    decstr = "%s %s %s " % (star_table_row['Dec deg'],
                            star_table_row['Dec min'],
                            star_table_row['Dec sec'])
    # Start with the target name
    ret = "%s %s %s 2000 " % (str(star_table_row['name']),rastr, decstr)

    # Proper motion RA and DEC
    ret += 'pmra=%.4f ' % (star_table_row['pmRA'])
    ret += 'pmdec=%.4f ' % (star_table_row['pmDEC'])
    # V Mag
    ret += 'vmag=%.2f ' % (star_table_row['Vmag'])

    # T Exp
    if temp:
        ret += 'texp=1200 '
    else:
        ret += 'texp=%d ' % int(star_table_row['texp'])

    # I2
    ret += 'I2=%s ' % (I2)
    # lamp
    ret += 'lamp=none '
    # start time
    ret += 'uth=%02d utm=%02d ' % (int(t.hour),int(t.minute))

    # Exp Count
    if star_table_row['expcount'] > SchedulerConsts.EXP_LIM:
        ret += 'expcount=%.3g ' % (SchedulerConsts.EXP_LIM)
    elif temp:
        ret += 'expcount=%.3g ' % (1e9)
    else:
        ret += 'expcount=%.3g ' % (star_table_row['expcount'])
    # Decker
    ret += 'decker=%s ' % (decker)
    # do flag
    if star_table_row['do']:
        ret += 'do=Y '
    else:
        ret += 'do= '
    # Count
    if temp:
        count = num_template_exp(star_table_row['Vmag'])
    else:
        count = int(star_table_row['nexp'])

    ret += 'count=%d ' % (count)

    ret += 'foc=%d ' % (int(focval))

    if owner != '':
        if owner == 'RECUR_A100':
            owner = 'public'
        ret += 'owner=%s ' % str(owner)

    if coverid != '':
        ret += 'coverid=%s ' % str(coverid)

    ret += 'binning=%s ' % str(star_table_row['binning'])

#    if star_table_row['mode'] != None:
#        if star_table_row['mode'] == BLANK:
#            ret += ' blank=Y'
#        elif star_table_row['mode'] == ACQUIRE:
#            ret += ' guide=Y'
#    else:
#        ret += ''

#    raoff  = star_table_row['raoff']
#    decoff = star_table_row['decoff']
#    if raoff == 'None':
#        raoff = ''
#    if decoff == 'None':
#        decoff = ''
#    if raoff is not '' and decoff is not '':
#        ret += ' raoff=' + str(raoff) + ' decoff=' + str(decoff)

    return str(ret)

def calc_elevations(stars, observer):
    '''
    els = calc_elevations(stars, observer)
    stars - list of ephem.FixedBody objects
    observer - ephem.Observer object
    els - numpy array of elevations in degrees of stars for observer
    '''
    els = []
    for s in stars:
        observer.date = ephem.Date(observer.date)
        s.compute(observer)
        cur_el = np.degrees(s.alt)
        els.append(cur_el)
    return np.array(els)

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


def make_APF_obs(dt, horizon=str(SchedulerConsts.TARGET_ELEVATION_MIN)):
    '''
    apf_obs = make_APF_obs(dt, horizon=str(TARGET_ELEVATION_MIN))
    dt - datetime object
    horizon - string of horizon in degrees
    apf_obs - ephem.Observer object for the time dt with the horizon set to horizon
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
    computes time in seconds before sunset
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
    helper to compute just sunset, calls compute_sunset_rise
    '''
    sunset, _ = compute_sunset_rise(dt, horizon=horizon)
    return sunset

def compute_sunrise(dt, horizon='0'):
    '''
    sunrise = compute_sunrise(dt, horizon='0')
    helper to compute just sunrise, calls compute_sunset_rise
    '''
    _, sunrise = compute_sunset_rise(dt, horizon=horizon)
    return sunrise


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

    if 'moon' in star_table.colnames:
        available = (star_table['moon'] > moon.moon_phase) & available

    if 'transparency' in star_table.colnames:
        ext = 2.5 * np.log10(slowdown)
        available = (star_table['transparency'] > ext) & available


    return available


def template_conditions(moon, seeing, slowdown):
    """ istrue = template_condition(moon, seeing, slowdown)

    Checks to see if moon, seeing and slowdown factor are within template conditions

    istrue - a simple Boolean

    moon - phase value from pyephem, ranges from 0 to 100 (a percentage)
    seeing - size in pixels
    slowdown - relative to clear

    """

    if seeing < SchedulerConsts.SEEING_TEMP and slowdown < SchedulerConsts.SLOWDOWN_TEMP:
        apflog("moon.phase=%.2f moon.alt=%.2f" % (moon.phase,moon.alt),echo=True,level='debug')
        if moon.phase < 50 and float(moon.alt) < 0:
            return True
        elif moon.phase < 25 and float(moon.alt) < 0.7:
            return True
        else:
            return False
    else:
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

def num_template_exp(vmag):
    '''
    num_template_exp(vmag)

    vmag - V magnitude of target
    count - number of exposures for a template observation

    '''
    count = 7

    if vmag > 10:
        count = 9

    elif vmag  < 8:
        count = 5


    return count

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

    count = num_template_exp(star_table['Vmag'][idx])

    tot_time = count * 1200

    tot_time += 210 + (2*40 + 40*(star_table['nexp'][idx]-1)) + 2400 # two B star exposures + three 70 second acquisitions and the actual observation readout times
    vis, star_elevations, scaled_els = Visible.visible(apf_obs, [stars[idx]], [tot_time])
    time_left_before_sunrise = compute_sunrise(dt, horizon='-18')

    try:
        apflog("enough_time_templates(): time for obs= %.1f  time until sunrise= %.1f " % (tot_time, time_left_before_sunrise),echo=True)
    except:
        apflog("enough_time_templates(): cannot log times!?!",echo=True)

    if tot_time < time_left_before_sunrise  and vis and time_left_before_sunrise < 14*3600.:
        return True
    else:
        return False


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


def make_obs_block(star_table, idx, dt, focval):
    '''

    make_obs_block(star_table, idx, dt, focval)

    star_table - astropy table of targets
    idx - index of target in star_table
    dt - datetime object
    focval - focus value

    rv - list of scriptobs lines for an obsblock
    '''

    rv = []

    cur_obsblock = star_table['obsblock'][idx]

    allinblock = star_table['obsblock'] == cur_obsblock
    allinblock = allinblock & (star_table['sheetn'] == star_table['sheetn'][idx])

    if np.any(star_table['mode'][allinblock] == FIRST):
        first = star_table['mode'][allinblock] == FIRST
    elif np.any(star_table['mode'][allinblock] == ACQUIRE):
        first = star_table['mode'][allinblock] == ACQUIRE
    else:
        first = None

    if np.any(star_table['mode'][allinblock] == LAST):
        last = star_table['mode'][allinblock] == LAST
    else:
        last = None

    rest = star_table['mode'][allinblock] != FIRST
    rest = rest & (star_table['mode'][allinblock] != ACQUIRE)
    rest = rest & (star_table['mode'][allinblock] != LAST)
    rest_idxs, = np.where(rest)


    if np.any(first):
        first_idxs, = np.where(first)
        for idx in first_idxs:
            scriptobs_line = make_scriptobs_line(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                                owner=star_table['sheetn'][allinblock][idx], \
                                                I2=star_table['I2'][allinblock][idx], focval=focval)
            rv.append(scriptobs_line)

    for idx in rest_idxs:
        scriptobs_line = make_scriptobs_line(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                               owner=star_table['sheetn'][allinblock][idx], \
                                               I2=star_table['I2'][allinblock][idx], focval=focval)
        rv.append(scriptobs_line)

    if np.any(last):
        last_idxs, = np.where(last)
        for idx in last_idxs:
            scriptobs_line = make_scriptobs_line(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                                owner=star_table['sheetn'][allinblock][idx], \
                                                I2=star_table['I2'][allinblock][idx], focval=focval)
            rv.append(scriptobs_line)


    rv.reverse()
    rv[0] += ' # obsblock=%s end' % (cur_obsblock)
    return rv

def make_result(stars, star_table, totexptimes, final_priorities, dt, idx, focval=0, bstar=False, mode=''):
    '''

    make_result(stars, star_table, totexptimes, final_priorities, dt, idx, focval=0, bstar=False, mode='')

    stars - list of ephem.FixedBody objects
    star_table - astropy table of targets
    totexptimes - numpy array of total exposure times
    final_priorities - numpy array of final priorities
    dt - datetime object
    idx - index of target in star_table
    focval - focus value
    bstar - boolean, True if target is a B star
    mode - string, mode of observation

    res - dictionary of target information
    '''
    res = dict()

    res['RA'] = stars[idx].a_ra
    res['DEC'] = stars[idx].a_dec
    res['PM_RA'] = star_table['pmRA'][idx]
    res['PM_DEC'] = star_table['pmDEC'][idx]
    res['VMAG'] = star_table['Vmag'][idx]
    res['BV'] = star_table['B-V'][idx]
    res['COUNTS'] = star_table['expcount'][idx]
    res['EXP_TIME'] = star_table['texp'][idx]
    res['NEXP'] = star_table['nexp'][idx]
    res['TOTEXP_TIME'] = totexptimes[idx]
    res['NAME'] = star_table['name'][idx]
    res['PRI'] = final_priorities[idx]
    res['DECKER'] = star_table['decker'][idx]
    res['I2'] = star_table['I2'][idx]
    res['BINNING'] = star_table['binning'][idx]
    res['isTemp'] = False
    res['isBstar'] = bstar
    res['mode'] = ''
    res['owner'] = star_table['sheetn'][idx]

#    if np.ma.is_masked(star_table[idx]['obsblock']):
#        res['obsblock'] = ''

    res['SCRIPTOBS'] = []
    if bstar:
        scriptobs_line = make_scriptobs_line(star_table[idx], dt, decker=res['DECKER'], \
                                         owner=res['owner'], I2='N', \
                                            focval=0)
        scriptobs_line = scriptobs_line + " # end"
        res['SCRIPTOBS'].append(scriptobs_line)

    scriptobs_line = make_scriptobs_line(star_table[idx], dt, decker=res['DECKER'], \
                                         owner=res['owner'], I2=star_table['I2'][idx], \
                                            focval=focval)

    scriptobs_line = scriptobs_line + " # end"
    res['SCRIPTOBS'].append(scriptobs_line)
#    else:
#        res['obsblock'] = star_table['obsblock'][idx]
#        res['SCRIPTOBS'] = make_obs_block(star_table, idx, dt, focval)

    return res

def last_attempted():
    """

    last_attempted()

    failed_obs - string of the last object attempted

    searches for the last object attempted in the apftask ktl variables
    SCRIPTOBS_LINE and SCRIPTOBS_LINE_RESULT

    If the last object was not observed successfully,
    returns the name of the object

    Returns the last object attempted to be observed
    if the observation failed.
    If it cannot read the keyword, returns None.

    """
    failed_obs = None

    try:
        last_line = ktl.read("apftask", "SCRIPTOBS_LINE")
        last_obj = last_line.split()[0]
    except:
        return None


    try:
        last_result = ktl.read("apftask", "SCRIPTOBS_LINE_RESULT", binary=True)
    except:
        return None

    apflog( "last_attempted(): Last objects attempted %s" % (last_obj), echo=True)
    # 3 is success
    if last_result != 3:
        failed_obs = last_obj
        apflog( "last_attempted(): Failed to observe %s" % (last_obj), echo=True)

    return failed_obs


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
    moon_dist = np.degrees(np.sqrt((moon.ra - ras)**2 + (moon.dec - decs)**2))

    moon_check = moon_dist > min_moon_dist

    return moon_check

def config_defaults(owner):
    '''
    config_defaults(owner)
    owner - string, owner of the targets

    config - dictionary of default values for the config
    '''

    config = dict()
    config['I2'] = 'Y'
    config['decker'] = 'W'
    config['mode'] = ''
    config['obsblock'] = ''
    config['Bstar'] = 'N'
    config['owner'] = owner
    config['inst'] = 'levy'
    config['raoff'] = ''
    config['decoff'] = ''

    return config

def get_next(ctime, seeing, slowdown, bstar=False, template=False, \
                sheetns=["RECUR_A100",], owner='public', \
                outfn="googledex.dat", toofn="too.dat", \
                outdir=None, focval=0, inst='', \
                rank_sheetn='rank_table', start_time=None):
    """ Determine the best target to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """

    global last_objs_attempted

    if not outdir:
        outdir = os.getcwd()

    dt = compute_datetime(ctime)

    config = config_defaults(owner)

    apflog( "get_next(): Finding target for time %s" % (dt), echo=True)

    if slowdown > SchedulerConsts.SLOWDOWN_MAX:
        log_str = "get_next(): Slowndown value of %f " % (slowdown)
        log_str += "exceeds maximum of %f at time %s" % (SchedulerConsts.SLOWDOWN_MAX, dt)
        apflog(log_str , echo=True)
        return None

    try:
        apfguide = ktl.Service('apfguide')
        stamp = apfguide['midptfin'].read(binary=True)
        ptime = datetime.datetime.utcfromtimestamp(stamp)
    except:
        if type(dt) == datetime.datetime:
            ptime = dt
        else:
            ptime = datetime.datetime.utcfromtimestamp(int(time.time()))

    apflog("get_next(): Updating star list with previous observations", echo=True)
    observed, star_table = ParseUCOSched.update_local_starlist(ptime,\
                                                               outfn=outfn, toofn=toofn, \
                                                                observed_file="observed_targets")

    rank_table = make_rank_table(rank_sheetn)
    hour_table = make_hour_table(rank_table, ptime)

    if hour_table is not None:
        hour_table = update_hour_table(hour_table, observed, ptime)

    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians

    if star_table is None:
        apflog("get_next(): Parsing the star list", echo=True)
        star_table, stars = ParseUCOSched.parse_UCOSched(sheetns=sheetns, \
                                                         outfn=outfn, outdir=outdir, \
                                                            config=config)
    else:
        stars = ParseUCOSched.gen_stars(star_table)
    targ_num = len(stars)

    # List of targets already observed

    last_failure = last_attempted()
    if last_failure is not None:
        last_objs_attempted.append(last_failure)

    ###
    # Need to update the googledex with the lastObserved date for observed targets
    # Scriptobs line uth utm can be used for this
    # Need to convert a uth and utm to a JD quickly.
    # timedelta = now - uth,utm : minus current JD?
    ###

    apf_obs = make_APF_obs(dt)

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    do_templates = template and template_conditions(moon, seeing, slowdown)

    apflog("get_next(): Will attempt templates = %s" % str(do_templates) ,echo=True)
    # Note which of these are B-Stars for later.
    bstars = (star_table['Bstar'] == 'Y')|(star_table['Bstar'] == 'y')

    apflog("get_next(): Computing exposure times", echo=True)
    totexptimes = tot_exp_times(star_table, targ_num)

    available = np.ones(targ_num, dtype=bool)
    cur_elevations = np.zeros(targ_num, dtype=float)
    scaled_elevations = np.zeros(targ_num, dtype=float)

    # Is the target behind the moon?

    moon_check = behind_moon(moon, star_table['ra'], star_table['dec'])
    available = available & moon_check
    log_str = "get_next(): Moon visibility check - stars rejected = "
    log_str += "%s" % ( np.asarray(star_table['name'][np.logical_not(moon_check)]))
    apflog(log_str, echo=True)

    sun_el_good = sun_el_check(star_table, apf_obs, horizon='-18')
    available = available & sun_el_good

    # other condition cuts (seeing, transparency, moon phase)
    cuts = condition_cuts(moon, seeing, slowdown, star_table)
    available = available & cuts

    if len(last_objs_attempted)>0:
        for n in last_objs_attempted:
            attempted = star_table['name'] == n
            available = available & np.logical_not(attempted) # Available and not observed

    if bstar:
        # We just need a B star
        apflog("get_next(): Selecting B stars", echo=True)
        available = available & bstars
        shiftwest = False
    else:
        apflog("get_next(): Culling B stars", echo=True)
        available = available & np.logical_not(bstars)
        shiftwest = True

    # Is the exposure time too long?
    apflog("get_next(): Removing really long exposures", echo=True)
    time_good = time_check(star_table, totexptimes, dt, start_time=start_time)

    available = available & time_good
    if np.any(available) is False:
        apflog( "get_next(): Not enough time left to observe any targets", level="error", echo=True)
        return None


    apflog("get_next(): Computing star elevations",echo=True)
    fstars = [s for s,_ in zip(stars,available) if _ ]
    vis, star_elevations, scaled_els = Visible.visible(apf_obs, fstars, \
                                                       totexptimes[available], shiftwest=shiftwest)
    currently_available = available
    if len(star_elevations) > 0:
        currently_available[available] = currently_available[available] & vis
    else:
        apflog( "get_next(): Couldn't find any suitable targets!", level="error", echo=True)
        return None

    cur_elevations[available] += star_elevations[vis]
    scaled_elevations[available] += scaled_els[vis]

    if slowdown > SchedulerConsts.SLOWDOWN_THRESH or seeing > SchedulerConsts.SEEING_THRESH:
        bright_enough = star_table['Vmag'] < SchedulerConsts.SLOWDOWN_VMAG_LIM
        available = available & bright_enough

    if not do_templates:
        available = available & (star_table['only_template'] == 'N')

    # Now just sort by priority, then cadence. Return top target
    if len(star_table['name'][available]) < 1:
        apflog( "get_next(): Couldn't find any suitable targets!", level="error", echo=True)
        return None

    final_priorities = compute_priorities(star_table,dt,
                                             rank_table=rank_table,
                                             hour_table=hour_table)

    try:
        pri = max(final_priorities[available])
        sort_i = (final_priorities == pri) & available
    except:
        apflog( "get_next(): Couldn't find any suitable targets!", level="error", echo=True)
        return None

    if bstar:
        sort_j = cur_elevations[sort_i].argsort()[::-1]
        focval=2
    else:
        sort_j = scaled_elevations[sort_i].argsort()[::-1]

    allidx, = np.where(sort_i)
    idx = allidx[sort_j][0]

    t_n = star_table['name'][idx]
    o_n = star_table['sheetn'][idx]
    p_n = final_priorities[idx]

    apflog("get_next(): selected target %s for program %s at priority %.0f" % (t_n, o_n, p_n) )
    nmstr= "get_next(): star names %s" % (np.asarray(star_table['name'][sort_i][sort_j]))
    pristr= "get_next(): star priorities %s" % (np.asarray(final_priorities[sort_i][sort_j]))
    mxpristr= "get_next(): max priority %d" % (pri)
    shstr= "get_next(): star sheet names %s" % (np.asarray(star_table['sheetn'][sort_i][sort_j]))
    if bstar:
        elstr= "get_next(): Bstar current elevations %s" % (cur_elevations[sort_i][sort_j])
    else:
        elstr= "get_next(): star scaled elevations %s" % (scaled_elevations[sort_i][sort_j])
    apflog(nmstr, echo=True)
    apflog(shstr, echo=True)
    apflog(pristr, echo=True)
    apflog(mxpristr, echo=True)
    apflog(elstr, echo=True)

    stars[idx].compute(apf_obs)

    take_template = do_templates and star_table['Template'][idx] == 'N' \
        and star_table['I2'][idx] == 'Y'
    if star_table['only_template'][idx] == 'Y' and do_templates:
        take_template = True

    res =  make_result(stars, star_table, totexptimes, final_priorities, dt, \
                       idx, focval=focval, bstar=bstar, mode=config['mode'])
    if take_template and bstar is False:
        bidx, bfinidx = find_Bstars(star_table, idx, bstars)

        if enough_time_templates(star_table,stars,idx,apf_obs,dt):
            bline = make_scriptobs_line(star_table[bstars][bidx], dt, \
                                        decker="N", I2="Y", owner=res['owner'], focval=2)
            line  = make_scriptobs_line(star_table[idx], \
                                        dt, decker="N", I2="N", owner=res['owner'], temp=True)
            bfinline = make_scriptobs_line(star_table[bstars][bfinidx], dt,\
                                            decker="N", I2="Y", owner=res['owner'], focval=0)
            res['SCRIPTOBS'] = []
            res['SCRIPTOBS'].append(bfinline + " # temp=Y end")
            res['SCRIPTOBS'].append(line + " # temp=Y")
            res['SCRIPTOBS'].append(bline + " # temp=Y")
            res['isTemp'] = True
            apflog("Attempting template observation of %s" % (star_table['name'][idx]), echo=True)

    return res

if __name__ == '__main__':

    t_dt = datetime.datetime.now()

    cfn = os.path.join('.','time_left.csv')
    if os.path.exists(cfn):
        hour_constraints = astropy.io.ascii.read(cfn)
    else:
        hour_constraints = None

    RANK_TABLEN='2024A_ranks'
    trank_table = make_rank_table(RANK_TABLEN, hour_constraints=hour_constraints)

    thour_table = make_hour_table(trank_table, t_dt, hour_constraints=hour_constraints)

    tsheet_list = list(trank_table['sheetn'][trank_table['rank'] > 0])


    try:
        ktl.write('apftask', 'SCRIPTOBS_LINE_RESULT', 3, binary=True)
    except:
        pass

    # For some test input what would the best target be?
    OTFN = "observed_targets"
    ot = open(OTFN, "w")
    starttime = time.time()
    result = get_next(starttime, 7.99, 0.4, bstar=True, \
                      sheetns=tsheet_list, rank_sheetn=RANK_TABLEN)
    while len(result['SCRIPTOBS']) > 0:
        ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
    ot.close()

    for i in range(5):

        result = get_next(starttime, 7.99, 0.4, bstar=False, sheetns=tsheet_list, \
                         template=True, rank_sheetn=RANK_TABLEN)
        #result = smartList("tst_targets", time.time(), 13.5, 2.4)

        if result is None:
            print("Get None target")

        while len(result["SCRIPTOBS"]) > 0:
            ot = open(OTFN, "a+")
            while len(result['SCRIPTOBS']) > 0:
                ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
            ot.close()
            starttime += result["TOTEXP_TIME"]

    print("Done")
    ot.close()

    print("Testing a failure")
    try:
        ktl.write('apftask', 'SCRIPTOBS_LINE_RESULT', 2, binary=True)
    except:
        pass
    result = get_next(starttime, 7.99, 0.4, bstar=False, sheetns=tsheet_list, \
                     template=True, rank_sheetn=RANK_TABLEN)


    print("Testing templates")

    tstar_table, tstars = ParseUCOSched.parse_UCOSched(sheetns=tsheet_list, \
                                                     outfn='googledex.dat', outdir=".", \
                                                        config=config_defaults('public'))
    tidx, = np.asarray(tstar_table['name'] == '185144').nonzero()
    tidx = tidx[0]
    tbstars = (tstar_table['Bstar'] == 'Y')|(tstar_table['Bstar'] == 'y')
    tbidx, tbfinidx = find_Bstars(tstar_table, tidx, tbstars)
    tbline = make_scriptobs_line(tstar_table[tbstars][tbidx], t_dt, \
                                decker="N", I2="Y", owner='public', focval=2)
    tline  = make_scriptobs_line(tstar_table[tidx], t_dt, \
                                decker="N", I2="N", owner='public', temp=True)
    tbfinline = make_scriptobs_line(tstar_table[tbstars][tbfinidx], t_dt, \
                                   decker="N", I2="Y", owner='public', focval=0)
    temp_res= []
    temp_res.append(tbfinline + " # temp=Y end")
    temp_res.append(tline + " # temp=Y")
    temp_res.append(tbline + " # temp=Y")
    out_r = [print(r) for r in temp_res]
