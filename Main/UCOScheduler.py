# UCOScheduler_V1.py
from __future__ import print_function
import os
import sys
import time
import re
from datetime import datetime, timedelta

import numpy as np
import astropy
import astropy.io
import ephem
from ExposureCalculations import getI2_M, getI2_K, getEXPMeter, getEXPMeter_Rate, getEXPTime
import ParseUCOSched
import Coords
from SchedulerConsts import * # I know

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

def computePriorities(star_table,available,cur_dt,observed=None,hour_table=None,rank_table=None):
    # make this a function, have it return the current priorities, than change references to the star_table below into references to the current priority list
    new_pri = np.zeros_like(star_table['APFpri'])
    new_pri += star_table['APFpri']

    if hour_table is not None:
        too_much = hour_table['cur']  > hour_table['tot']
        done_sheets = hour_table['sheetn'][too_much]
    else:
        done_sheets = []

    if rank_table is not None:
        for sheetn in rank_table['sheetn']:
            if sheetn not in done_sheets:
                cur = star_table['sheetn'] == sheetn
                new_pri[cur] += rank_table['rank'][rank_table['sheetn'] == sheetn]
            else:
                apflog("Sheet %s has exceeded it's allocation for the night" % (sheetn),echo=True)

    return new_pri

def updateHourTable(hour_table,observed,dt,outfn='hour_table',outdir=None):
    '''
    updateHourTableobserved_logs,outfn='hour_table')

    Updates hour_table with history of observations.

    '''

    if not outdir :
        outdir = os.getcwd()

    outfn = os.path.join(outdir,outfn)

    hours = dict()

    # observed objects have lists as attributes
    # in reverse time order, so most recent target observed is first.
    nobj = len(observed.names)
    for i in range(0,nobj):
            own = observed.owners[i]
            if own not in hours.keys():
                    hours[own] = 0.0

    cur = dt
    for i in range(0,nobj):
            hr, mn = observed.times[i]
            prev = datetime(dt.year,dt.month,dt.day,hr,mn)
            diff = cur - prev
            hourdiff = (diff.days * 24 + diff.seconds / 3600.)
            if hourdiff > 0:
                hours[observed.owners[i]] += hourdiff
            cur = prev

    for ky in hours.keys():
        if ky == 'public':
            hour_table['cur'][hour_table['sheetn'] == 'RECUR_A100'] = hours[ky]
        else:
            hour_table['cur'][hour_table['sheetn'] == ky] = hours[ky]

    try:
        hour_table.write(outfn,format='ascii',overwrite=True)
    except Exception as e:
        apflog("Cannot write table %s: %s" % (outfn,e),level='error',echo=True)

    return hour_table


def makeHourTable(sheet_table_name,dt,outfn='hour_table',outdir=None,frac_fn='frac_table',hour_constraints=None):

    if not outdir :
        outdir = os.getcwd()

    outfn = os.path.join(outdir,outfn)

    if os.path.exists(outfn):
        hour_table =  astropy.table.Table.read(outfn,format='ascii')
        return hour_table

    frac_fn = os.path.join(outdir,frac_fn)
    if os.path.exists(frac_fn):
        frac_table = astropy.table.Table.read(frac_fn,format='ascii')
    else:
        sheetns, fracs = ParseUCOSched.parseFracTable(sheet_table_name=sheet_table_name,outfn=frac_fn)
        frac_table = []
        for i in range(0,len(fracs)):
            frow = []
            frow.append(sheetns[i])
            frow.append(fracs[i])
            frac_table.append(frow)
        frac_table = astropy.table.Table(rows=frac_table,names=['sheetn','frac'])
        try:
            frac_table.write(frac_fn,format='ascii')
        except Exception as e:
            apflog("Cannot write table %s: %s" % (frac_fn,e),level='error',echo=True)


    hour_table= astropy.table.Table(frac_table,names=['sheetn','frac'])

    sunset,sunrise = computeSunsetRise(dt,horizon='-9')
    if sunrise < sunset:
        sunrise += 86400
    tot = sunrise - sunset
    tot /= 3600.

    hour_table['tot'] =np.abs(tot*hour_table['frac'])
    hour_table['cur'] =0.0*hour_table['frac']

    hour_table['tot'][hour_table['tot'] < 1.0] = 1.0
    
    if hour_constraints is not None:
        if 'runname' in hour_constraints.keys() and 'left' in hour_constraints.keys():
            for runname in hour_constraints['runname']:
                if hour_constraints['left'][hour_constraints['runname']==runname] < hour_table['tot'][hour_table['sheetn']==runname]:
                    hour_table['tot'][hour_table['sheetn']==runname] = hour_constraints['left'][hour_constraints['runname']==runname]
                elif hour_constraints['left'][hour_constraints['runname']==runname] < 0:
                    hour_table['tot'][hour_table['sheetn']==runname] = -1.0

    try:
        hour_table.write(outfn,format='ascii')
    except Exception as e:
        apflog("Cannot write table %s: %s" % (outfn,e),level='error',echo=True)
    return hour_table


def timeCheck(star_table,totexptimes,dt,hour_table):

    maxexptime = TARGET_EXPOSURE_TIME_MAX
    time_left_before_sunrise = computeSunrise(dt,horizon='-9')
    if maxexptime > time_left_before_sunrise:
        maxexptime = time_left_before_sunrise
    if maxexptime < TARGET_EXPOSURE_TIME_MIN:
        maxexptime = TARGET_EXPOSURE_TIME_MIN
        # this will try a target in case we get lucky
        
    time_check = totexptimes <= maxexptime

    hour_table['left'] = hour_table['tot'] - hour_table['cur']
    hour_table['left'] *= 3600.

    
    program_times = np.zeros_like(totexptimes)
    for sheetn in hour_table['sheetn']:
        program_times[star_table['sheetn'] == sheetn] += hour_table['left'][hour_table['sheetn'] == sheetn]

    time_check = time_check & (totexptimes < program_times)

    return time_check
        
def makeRankTable(sheet_table_name,outfn='rank_table',outdir=None):

    if not outdir :
        outdir = os.getcwd()

    outfn = os.path.join(outdir,outfn)
    if os.path.exists(outfn):
        rank_table = astropy.table.Table.read(outfn,format='ascii')
    else:
        sheetns, ranks = ParseUCOSched.parseRankTable(sheet_table_name=sheet_table_name)

        rank_table= astropy.table.Table([sheetns,ranks],names=['sheetn','rank'])
        try:
            rank_table.write(outfn,format='ascii')
        except Exception as e:
            apflog("Cannot write table %s: %s" % (outfn,e),level='error',echo=True)

    return rank_table


def makeScriptobsLine(star_table_row, t, decker="W", I2="Y", owner='public', focval=0, coverid='',temp=False):
    """ given a name, a row in a star table and a do_flag, will generate a scriptobs line as a string
    line = makeScriptobsLine(star_table_row, t, decker="W",I2="Y")
    star_table_row -contains all of the data needed for the line except

    t - a datetime object, this is used to fill in the uth and utm fields,
    decker - one character field for the decker, defaults to "W"
    I2 - one character field for whether or not the Iodine cell is in, must be "Y" or "N"
    temp - a boolean for whether or not this is a template observation
    """

    """Takes a line from the star table and generates the appropriate line to pass to scriptobs. """
    # Start with the target name

    # Add the RA as three elements, HR, MIN, SEC
    rastr = "%s %s %s " % (star_table_row['RA hr'],star_table_row['RA min'],star_table_row['RA sec'])

    # Add the DEC as three elements, DEG, MIN, SEC
    decstr = "%s %s %s " % (star_table_row['Dec deg'],star_table_row['Dec min'],star_table_row['Dec sec'])

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
    if star_table_row['expcount'] > EXP_LIM:
        ret += 'expcount=%.3g ' % (EXP_LIM)
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
        if star_table_row['Vmag'] > 10:
            count = 9
        elif star_table_row['Vmag'] < 8:
            count = 5
        else:
            count = 7
    else:
        count = int(star_table_row['APFnshots'])

    ret += 'count=%d ' % (count)

    ret += 'foc=%d ' % (int(focval))

    if owner != '':
        if owner == 'RECUR_A100':
            owner = 'public'
        ret += 'owner=%s ' % str(owner)

    if coverid != '':
        ret += 'coverid=%s ' % str(coverid)

    if star_table_row['mode'] != None:
        if star_table_row['mode'] == BLANK:
            ret += ' blank=Y'
        elif star_table_row['mode'] == ACQUIRE:
            ret += ' guide=Y'
    else:
        ret += ''

    raoff  = star_table_row['raoff']
    decoff = star_table_row['decoff']
    if raoff == 'None':
        raoff = ''
    if decoff == 'None':
        decoff = ''
#    if raoff is not '' and decoff is not '':
#        ret += ' raoff=' + str(raoff) + ' decoff=' + str(decoff)

    return str(ret)

def calc_elevations(stars, observer):
    els = []
    for s in stars:
        observer.date = ephem.Date(observer.date)
        s.compute(observer)
        cur_el = np.degrees(s.alt)
        els.append(cur_el)
    return np.array(els)

def computeDatetime(ctime):
    if type(ctime) == float:
        dt = datetime.utcfromtimestamp(int(ctime))
    elif type(ctime) == datetime:
        dt = ctime
    elif type(ctime) == ephem.Date:
        dt = ctime.datetime()
    else:
        #punt and use current UT
        dt = datetime.utcfromtimestamp(int(time.time()))
    return dt


def makeAPFObs(dt,horizon=str(TARGET_ELEVATION_MIN)):
    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = horizon
    apf_obs.date = dt

    return apf_obs

def computeSunsetRise(dt,horizon='0'):
    # computes time in seconds before sunset
    apf_obs = makeAPFObs(dt,horizon=horizon)
    sunset = apf_obs.next_setting(ephem.Sun())
    sunset -= ephem.Date(dt)
    sunset *= 86400.0 # convert to seconds

    sunrise = apf_obs.next_rising(ephem.Sun())
    sunrise -= ephem.Date(dt)
    sunrise *= 86400.0 # convert to seconds
    return sunset, sunrise

def computeSunset(dt,horizon='0'):

    sunset, sunrise = computeSunsetRise(dt,horizon=horizon)
    return sunset

def computeSunrise(dt,horizon='0'):
    sunset, sunrise = computeSunsetRise(dt,horizon=horizon)
    return sunrise


def conditionCuts(moon,seeing,slowdown,star_table):
    """ available = conditionCuts(moon, seeing, slowdown, star_table)

    Checks if columns are in the star_table, then cuts on those, returns a boolean numpy array

    available - Boolean numpy array of available targets

    moon - phase value from pyephem, ranges from 0 to 100 (a percentage)
    seeing - size in pixels
    transparency - magnitudes of extinction

    """

    if 'seeing' in star_table.colnames:
        available = star_table['seeing']/0.109 > seeing
    else:
        available = np.ones(len(star_table['ra'], dtype=bool))

    if 'moon' in star_table.colnames:
        available = (star_table['moon'] > moon) & available

    if 'transparency' in star_table.colnames:
        ext = 2.5 * np.log10(slowdown)
        available = (star_table['transparency'] > ext) & available


    return available


def templateConditions(moon, seeing, slowdown):
    """ istrue = conditionCuts(moon, seeing, slowdown)

    Checks to see if moon, seeing and slowdown factor are within template conditions

    istrue - a simple Boolean

    moon - phase value from pyephem, ranges from 0 to 100 (a percentage)
    seeing - size in pixels
    slowdown - relative to clear

    """

    if seeing < SEEING_TEMP and slowdown < SLOWDOWN_TEMP:
        apflog("moon.phase=%.2f moon.alt=%.2f" % (moon.phase,moon.alt),echo=True,level='debug')
        if moon.phase < 50 and float(moon.alt) < 0:
            return True
        elif moon.phase < 25 and float(moon.alt) < 0.7:
            return True
        else:
            return False
    else:
        return False

def findClosest(ras,decs,ra,dec):

    distances = np.sqrt((ra - ras)**2 + (dec - decs)**2)

    min_ind = distances.argmin()

    return min_ind


def enoughTimeTemplates(star_table,stars,idx,apf_obs,dt):
    tot_time = star_table['APFnshots'][idx]*star_table['texp'][idx]
    tot_time += 210 + (2*40 + 40*(star_table['APFnshots'][idx]-1)) + 2400 # two B star exposures + three 70 second acquisitions and the actual observation readout times
    vis, star_elevations, fin_els, scaled_els = Visible.visible(apf_obs,[stars[idx]],[tot_time])
    time_left_before_sunrise = computeSunrise(dt,horizon='-9')

    try:
        apflog("enoughTimeTemplates(): time for obs= %.1f  time until sunrise= %.1f " % (tot_time,  time_left_before_sunrise),echo=True)
    except:
        apflog("enoughTimeTemplates(): cannot log times!?!",echo=True)

    if tot_time < time_left_before_sunrise  and vis and time_left_before_sunrise < 14*3600.:
        return True
    else:
        return False


def findBstars(star_table,idx, bstars):

    near_idx = findClosest(star_table['ra'][bstars],star_table['dec'][bstars],star_table['ra'][idx],star_table['dec'][idx])

    end_idx = findClosest(star_table['ra'][bstars],star_table['dec'][bstars],(star_table['ra'][idx]+15*np.pi/180.),star_table['dec'][idx])


    return near_idx,end_idx


def makeObsBlock(star_table, idx, dt, focval):

    rv = []

    cur_obsblock = star_table['obsblock'][idx]

    allinblock = (star_table['obsblock'] == cur_obsblock)
    allinblock = allinblock & (star_table['sheetn'] == star_table['sheetn'][idx])

    if np.any(star_table['mode'][allinblock] == FIRST):
        first = (star_table['mode'][allinblock] == FIRST)
    elif np.any(star_table['mode'][allinblock] == ACQUIRE):
        first = (star_table['mode'][allinblock] == ACQUIRE)
    else:
        first = None

    if np.any(star_table['mode'][allinblock] == LAST):
        last = (star_table['mode'][allinblock] == LAST)
    else:
        last = None

    rest = (star_table['mode'][allinblock] != FIRST)
    rest = rest & (star_table['mode'][allinblock] != ACQUIRE)
    rest = rest & (star_table['mode'][allinblock] != LAST)
    rest_idxs, = np.where(rest)


    if np.any(first):
        first_idxs, = np.where(first)
        for idx in first_idxs:
            scriptobs_line = makeScriptobsLine(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                                owner=star_table['sheetn'][allinblock][idx], \
                                                I2=star_table['I2'][allinblock][idx], focval=focval)
            rv.append(scriptobs_line)

    for idx in rest_idxs:
        scriptobs_line = makeScriptobsLine(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                               owner=star_table['sheetn'][allinblock][idx], \
                                               I2=star_table['I2'][allinblock][idx], focval=focval)
        rv.append(scriptobs_line)

    if np.any(last):
        last_idxs, = np.where(last)
        for idx in last_idxs:
            scriptobs_line = makeScriptobsLine(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                                owner=star_table['sheetn'][allinblock][idx], \
                                                I2=star_table['I2'][allinblock][idx], focval=focval)
            rv.append(scriptobs_line)


    rv.reverse()
    rv[0] += ' # obsblock=%s end' % (cur_obsblock)
    return(rv)

def makeResult(stars,star_table,totexptimes,dt,idx,focval=0,bstar=False,mode=''):
    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table['pmRA'][idx]
    res['PM_DEC'] = star_table['pmDEC'][idx]
    res['VMAG']   = star_table['Vmag'][idx]
    res['BV']     = star_table['B-V'][idx]
    res['COUNTS'] = star_table['expcount'][idx]
    res['EXP_TIME'] = star_table['texp'][idx]
    res['NEXP'] = star_table['APFnshots'][idx]
    res['TOTEXP_TIME'] = totexptimes[idx]
    res['NAME']   = star_table['name'][idx]
    res['PRI']    = star_table['APFpri'][idx]
    res['DECKER'] = star_table['decker'][idx]
    res['I2']     = star_table['I2'][idx]
    res['isTemp'] =    False
    res['isBstar'] =    bstar
    res['mode']   =   ''
    res['owner'] =    star_table['sheetn'][idx]

    if np.ma.is_masked(star_table[idx]['obsblock']):
        res['obsblock'] = ''
        
        res['SCRIPTOBS'] = []
        scriptobs_line = makeScriptobsLine(star_table[idx], dt, decker=res['DECKER'], owner=res['owner'], I2=star_table['I2'][idx], focval=focval)
        scriptobs_line = scriptobs_line + " # end"
        res['SCRIPTOBS'].append(scriptobs_line)
    else:
        res['obsblock'] = star_table['obsblock'][idx]
        res['SCRIPTOBS'] = makeObsBlock(star_table, idx, dt, focval)

    return res

def lastAttempted(observed):

    failed_obs = None

    try:
        lastresult = ktl.read("apftask","SCRIPTOBS_LINE")
        lastobj = lastresult.split()[0]
    except:
        return None

    if lastobj not in observed.names:
        apflog( "lastAttempted(): Last objects attempted %s" % (lastobj),echo=True)
        failed_obs = lastobj

    return failed_obs


def behindMoon(moon,ras,decs):
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN
    moonDist = np.degrees(np.sqrt((moon.ra - ras)**2 + (moon.dec - decs)**2))

    apflog("behindMoon(): Culling stars behind the moon",echo=True)
    moon_check = moonDist > minMoonDist

    return moon_check

def configDefaults(owner):
    config = dict()
    config['I2'] = 'Y'
    config['decker']='W'
    config['mode']=''
    config['obsblock']=''
    config['Bstar']='N'
    config['owner']=owner
    config['inst']='levy'
    config['raoff'] = ''
    config['decoff'] = ''

    return config

def getNext(ctime, seeing, slowdown, bstar=False,template=False,sheetns=["RECUR_A100"],owner='public',outfn="googledex.dat",toofn="too.dat",outdir=None,focval=0,inst='',rank_sheetn='rank_table',frac_sheet=None):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """

    global last_objs_attempted

    if not outdir:
        outdir = os.getcwd()

    dt = computeDatetime(ctime)

    config = configDefaults(owner)

    apflog( "getNext(): Finding target for time %s" % (dt),echo=True)

    if slowdown > SLOWDOWN_MAX:
        apflog( "getNext(): Slowndown value of %f exceeds maximum of %f at time %s" % (slowdown,SLOWDOWN_MAX,dt),echo=True)
        return None

    try:
        apfguide = ktl.Service('apfguide')
        stamp = apfguide['midptfin'].read(binary=True)
        ptime = datetime.utcfromtimestamp(stamp)
    except:
        if type(dt) == datetime:
            ptime = dt
        else:
            ptime = datetime.utcfromtimestamp(int(time.time()))

    apflog("getNext(): Updating star list with previous observations",echo=True)
    observed, star_table = ParseUCOSched.updateLocalStarlist(ptime,outfn=outfn,toofn=toofn,observed_file="observed_targets")

    hour_table = None
    if frac_sheet is not None:
        hour_table = makeHourTable(frac_sheet,dt)
        hour_table = updateHourTable(hour_table,observed,dt)

    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians

    if star_table is None:
        apflog("getNext(): Parsing the star list",echo=True)
        star_table, stars = ParseUCOSched.parseUCOSched(sheetns=sheetns,outfn=outfn,outdir=outdir,config=config)
    else:
        stars = ParseUCOSched.genStars(star_table)
    targNum = len(stars)

    # List of targets already observed

    lastfailure = lastAttempted(observed)
    if lastfailure is not None:
        last_objs_attempted.append(lastfailure)
    if len(last_objs_attempted) == 5:
        apflog( "getNext(): 5 failed acquisition attempts",level="warn",echo=True)


    ###
    # Need to update the googledex with the lastObserved date for observed targets
    # Scriptobs line uth utm can be used for this
    # Need to convert a uth and utm to a JD quickly.
    # timedelta = now - uth,utm : minus current JD?
    ###

    apf_obs = makeAPFObs(dt)
    # APF latitude in radians
    apf_lat = (37 + 20/60. + 33.1/3600.) * np.pi/180.

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    do_templates = template and templateConditions(moon, seeing, slowdown)


    apflog("getNext(): Parsed the Googledex...",echo=True)


    apflog("getNext(): Finding B stars",echo=True)
    # Note which of these are B-Stars for later.
    bstars = (star_table['Bstar'] == 'Y')|(star_table['Bstar'] == 'y')

    # Distance to stay away from the moon


    totexptimes = np.zeros(targNum, dtype=float)
    totexptimes = star_table['texp'] * star_table['APFnshots'] + 40 * (star_table['APFnshots']-1)

    available = np.ones(targNum, dtype=bool)
    cur_elevations = np.zeros(targNum, dtype=float)
    scaled_elevations = np.zeros(targNum, dtype=float)

    # Is the target behind the moon?
    moon_check = behindMoon(moon,star_table['ra'],star_table['dec'])
    available = available & moon_check
    if len(last_objs_attempted)>0:
        for n in last_objs_attempted:
            attempted = (star_table['name'] == n)
            available = available & np.logical_not(attempted) # Available and not observed

    cadence_check = (ephem.julian_date(dt) - star_table['lastobs'])
    good_cadence = cadence_check >  star_table['APFcad']
    available = available & good_cadence

    if bstar:
        # We just need a B star
        apflog("getNext(): Selecting B stars",echo=True)
        available = available & bstars

    else:

        apflog("getNext(): Culling B stars",echo=True)
        available = available & np.logical_not(bstars)

    # Calculate the exposure time for the target
    # Want to pass the entire list of targets to this function

    apflog("getNext(): Computing exposure times",echo=True)
    exp_counts = star_table['expcount']

    # Is the exposure time too long?
    apflog("getNext(): Removing really long exposures",echo=True)
    time_check = timeCheck(star_table,totexptimes,dt,hour_table)

    available = available & time_check
    if np.any(available) == False:
        apflog( "getNext(): Not enough time left to observe any targets",level="error",echo=True)
        return None


    apflog("getNext(): Computing star elevations",echo=True)
    fstars = [s for s,_ in zip(stars,available) if _ ]
    vis,star_elevations,fin_star_elevations, scaled_els = Visible.visible(apf_obs, fstars, totexptimes[available],shiftwest=True)
    currently_available = available
    if len(star_elevations) > 0:
        currently_available[available] = currently_available[available] & vis
    else:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None
    
    cur_elevations[available] += star_elevations[vis]
    scaled_elevations[available] += scaled_els[vis]

    if slowdown > SLOWDOWN_THRESH or seeing > SEEING_THRESH:
        bright_enough = star_table['Vmag'] < SLOWDOWN_VMAG_LIM
        available = available & bright_enough

    # Now just sort by priority, then cadence. Return top target
    if len(star_table['name'][available]) < 1:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None


    final_priorities = computePriorities(star_table,available,dt,rank_table=makeRankTable(rank_sheetn),hour_table=hour_table,observed=observed)

    try:
        pri = max(final_priorities[available])
        sort_i = (final_priorities == pri) & available
    except:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
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

    apflog("getNext(): selected target %s for program %s at priority %.0f" % (t_n,o_n,p_n) )
    nmstr= "getNext(): star names %s" % (np.asarray(star_table['name'][sort_i][sort_j]))
    shstr= "getNext(): star sheet names %s" % (np.asarray(star_table['sheetn'][sort_i][sort_j]))
    if bstar:
        elstr= "getNext(): Bstar current elevations %s" % (cur_elevations[sort_i][sort_j])
    else:
        elstr= "getNext(): star scaled elevations %s" % (scaled_elevations[sort_i][sort_j])
    apflog(nmstr,echo=True)
    apflog(shstr,echo=True)
    apflog(elstr,echo=True)

    stars[idx].compute(apf_obs)

    res =  makeResult(stars,star_table,totexptimes,dt,idx,focval=focval,bstar=bstar,mode=config['mode'])
    if do_templates and star_table['Template'][idx] == 'N' and star_table['I2'][idx] == 'Y':
        bidx,bfinidx = findBstars(star_table,idx,bstars)

        if enoughTimeTemplates(star_table,stars,idx,apf_obs,dt):
            bline = makeScriptobsLine(star_table[bidx],dt,decker="N",I2="Y", owner=res['owner'],focval=2)
            line  = makeScriptobsLine(star_table[idx],dt,decker="N",I2="N", owner=res['owner'],temp=True)
            bfinline = makeScriptobsLine(star_table[bfinidx],dt,decker="N",I2="Y",owner=res['owner'],focval=0)
            res['SCRIPTOBS'] = []
            res['SCRIPTOBS'].append(bfinline + " # temp=Y end")
            res['SCRIPTOBS'].append(line + " # temp=Y")
            res['SCRIPTOBS'].append(bline + " # temp=Y")
            res['isTemp'] = True
            apflog("Attempting template observation of %s" % (star_table['name'][idx]),echo=True)

    return res

if __name__ == '__main__':

    dt = datetime.now()

    cfn = os.path.join('.','time_left.csv')
    if os.path.exists(cfn):
        hour_constraints = astropy.io.ascii.read(cfn)
    else:
        hour_constraints = None

    frac_tablen='2020B_frac'
    hour_table = makeHourTable(frac_tablen,dt,hour_constraints=hour_constraints)

    rank_tablen='2020B_ranks'
    rank_table = makeRankTable(rank_tablen)

#    sheetn=["2018B"]
    sheetn="RECUR_A100,2020B_A000,2020B_A008,2020B_A009,2020B_A010"

    # For some test input what would the best target be?
    otfn = "observed_targets"
    ot = open(otfn,"w")
    starttime = time.time()
    result = getNext(starttime, 7.99, 0.4, bstar=True,sheetns=sheetn.split(","),rank_sheetn=rank_tablen,frac_sheet=frac_tablen)
    ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
    ot.close()
    starttime += 400
    for i in range(5):

        result = getNext(starttime, 7.99, 0.4, bstar=False,sheetns=sheetn,template=True,rank_sheetn=rank_tablen,frac_sheet=frac_tablen)
        #result = smartList("tst_targets", time.time(), 13.5, 2.4)

        if result is None:
            print("Get None target")
        else:
            for k in result:
                print(k, result[k])
        while len(result["SCRIPTOBS"]) > 0:
            ot = open(otfn,"a")
            ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
            ot.close()
            starttime += result["EXP_TIME"]

    print("Done")
