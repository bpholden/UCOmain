# UCSCScheduler_V2.py
from __future__ import print_function
import os
import sys
import time
import re
from datetime import datetime, timedelta

import numpy as np
import ephem
from ExposureCalculations import getI2_M, getI2_K, getEXPMeter, getEXPMeter_Rate, getEXPTime
import ParseGoogledex
import ObservedLog
import Coords
from SchedulerConsts import * # I know

try:
    from apflog import *
    import ktl
except:
    from fake_apflog import *
import Visible


###
# arGGGGGG!! a global
###

last_objs_attempted = []

def computeMaxTimes(exp_times,maxtimes):
    fintimes = np.zeros_like(exp_times)
    user_selected_max = (exp_times > maxtimes)&(maxtimes>0)
    fintimes[user_selected_max] = maxtimes[user_selected_max]
    fintimes[exp_times < maxtimes] = exp_times[exp_times < maxtimes]
    fintimes[maxtimes<=0] = exp_times[maxtimes<=0]

    return fintimes


def computePriorities(star_table,available,cur_dt,flags,frac_table=None):
    # make this a function, have it return the current priorities, than change references to the star_table below into references to the current priority list
    if any(star_table[available, DS_DUR] > 0):
        new_pri = np.zeros_like(star_table[:, DS_APFPRI])
        new_pri[available] += star_table[available,DS_APFPRI]
        delta_pri = np.zeros_like(new_pri[available])
        timedependent, = np.where(star_table[available, DS_DUR] > 0)
        for tdinx in timedependent:
            sdt = datetime(cur_dt.year,cur_dt.month,cur_dt.day,int(star_table[:, DS_UTH][available][tdinx]),int(star_table[:, DS_UTM][available][tdinx]),0)
            durdelt = timedelta(0,star_table[:, DS_DUR][available][tdinx],0)
            if (cur_dt - sdt < durdelt) and (cur_dt - sdt > timedelta(0,0,0) ):
                delta_pri[tdinx] += PRI_DELTA
        new_pri[available] += delta_pri
    elif frac_table is not None:
        new_pri = np.zeros_like(star_table[:, DS_APFPRI])
        new_pri[available] += star_table[available,DS_APFPRI]
        too_much = frac_table[:,DS_FT_CUR]  > frac_table[:,DS_FT_TOT]
        done_sheets = frac_table[too_much,DS_FT_NAMES]
        for sn in done_sheets:
            bad = star_table[:, DS_SHEETN] == sn
            new_pri[bad] = 0
        
    else:
        new_pri = star_table[:, DS_APFPRI]
    return new_pri


def parseStarlist(starlist):
    """ Parse a scriptobs-compatible starlist for the scheduler.

    names, star_table, lines, stars = parseStarlist(starlist)
    starlist - a filename

    names - a list of stars in the starlist
    star_table - a numpy array
    lines - a list of strings that can be used for scriptobs input
    stars - a list of pyEphem objects
    """
    names = []
    lines = []
    stars = []
    star_table = []
    try:
        f = open(starlist,'r')
    except IOError:
        apflog("Warning: Could not open %s. No target can be selected." % starlist,echo=True)
        return None
    else:
        for line in f:
            if not re.search("\A\#", line):
                ls = line.split()
                names.append(ls[0])
                row = []
                # RA value in radians
                row.append(getRARad(ls[1], ls[2], ls[3]))
                # Dec value in radians
                row.append(getDECRad(ls[4], ls[5], ls[6]))
                # PM RA
                row.append(float(ls[8].split('=')[-1]))
                # PM Dec
                row.append(float(ls[9].split('=')[-1]))
                # V mag
                row.append(float(ls[10].split('=')[-1]))
                # Exposure time
                row.append(float(ls[11].split('=')[-1]))
                # Desired Counts
                row.append(float(ls[16].split('=')[-1]))
                # Filler not used here
                row.append(0.)
                row.append(0.)
                # Number of exposures
                row.append(int(ls[19].split('=')[-1]))

                star_table.append(row)

                # Save the scriptobs line for later
                lines.append(line)

                # Generate a pyEphem object for this target
                star = ephem.FixedBody()
                star.name = ls[0]
                star._ra = ephem.hours(":".join([ls[1], ls[2], ls[3]]))
                star._dec = ephem.degrees(":".join([ls[4], ls[5], ls[6]]))
                stars.append(star)

    return names, np.array(star_table), lines, stars



def makeScriptobsLine(name, row, do_flag, t, decker="W",I2="Y",owner='Vogt',focval=0):
    """ given a name, a row in a star table and a do_flag, will generate a scriptobs line as a string
    line = makeScriptobsLine(name, row, do_flag, t, decker="W",I2="Y")
    name - name of star, first column in line
    row - star_table row for star that begins with name, cotains all of the data needed for the line except
    do_flag - a string for whether or not scriptob needs to do a pointing check before slewing to the target
    t - a datetime object, this is used to fill in the uth and utm fields,
    decker - one character field for the decker, defaults to "W"
    I2 - one character field for whether or not the Iodine cell is in, must be "Y" or "N"
    """

    """Takes a line from the star table and generates the appropriate line to pass to scriptobs. """
    # Start with the target name
    ret = name + ' '
    # Add the RA as three elements, HR, MIN, SEC
    rastr = Coords.getCoordStr(np.degrees(row[DS_RA]), isRA=True)
    ret += rastr + ' '
    # Add the DEC as three elements, DEG, MIN, SEC
    decstr = Coords.getCoordStr(np.degrees(row[DS_DEC]))
    ret += decstr + ' '
    # Epoch
    ret += '2000 '
    # Proper motion RA and DEC
    ret += 'pmra=' + str(row[DS_PMRA]) + ' '
    ret += 'pmdec=' + str(row[DS_PMDEC]) + ' '
    # V Mag
    ret += 'vmag=' + str(row[DS_VMAG]) + ' '
    # T Exp
    if row[DS_EXPT] > MAX_EXPTIME:
        ret += 'texp=%d ' % (int(MAX_EXPTIME))
    elif row[DS_EXPT] <= MIN_EXPTIME:
        ret += 'texp=%d ' % (int(MIN_EXPTIME))
    else:
        ret += 'texp=' + str(int(row[DS_EXPT])) + ' '
    # I2
    ret += 'I2=%s ' % (I2)
    # lamp
    ret += 'lamp=none '
    # start time
    ret += 'uth=' + str(t.hour) + ' '
    ret += 'utm=' + str(t.minute) + ' '
    # Exp Count
    if row[DS_COUNTS] > 3e9:
        ret += 'expcount=%.3g' % (3e9) + ' '
    else:
        ret += 'expcount=%.3g' % (row[DS_COUNTS]) + ' '
    # Decker
    ret += 'decker=%s ' % (decker)
    # do flag
    if do_flag:
        ret += 'do=Y '
    else:
        ret += 'do= '
    # Count
    ret += 'count=' + str(int(row[DS_NSHOTS]))

    ret += ' foc=' + str(int(focval))

    if owner != '':
        ret += ' owner=' + str(owner)

    return ret


def calculateUCSCExposureTime(vmag, i2counts, elevation, seeing, bmv, deckers):
    """ calculateUCSCExposureTime uses the recipe from Burt et al. (2015) to compute the exposure time for a target.

    exp_time, exp_counts, i2counts = calculateUCSCExposureTime(vmag, precision, elevation, seeing, bmv, decker="W")
    vmag - numpy array of V magnitudes (Johnson filter, Vega mags)
    i2counts - the required number of median Iodine cell counts, this is calculated from the precision and color of the star, this, in effect, sets the exposure time.
    elevation - elevation of the star above the horizon at the start of the exposure
    seeing - FWHM of the seeing in pixels on the guider
    bmv - (B - V) for the star (both Johnson filters, Vega zeropoint)
    decker - apeture

    exp_time - a numpy array of times in seconds, are integer values
    exp_counts - values for the exposure meter, this can be a floating point value



    """
    vmag = np.array(vmag)
    i2counts = np.array(i2counts)
    bmv = np.array(bmv)


	# Now lets calculate the exposure times

    # minimum I2 counts so exposures are not rejected by P. Butler's DRP
    mini2_idx = np.where(i2counts < MIN_I2)
    if len(mini2_idx) > 0:
        i2counts[mini2_idx] = MIN_I2

	# Exposure Meter counts to reach desired I2 counts
    exp_counts = getEXPMeter(i2counts, bmv)
    #	exp_counts = 1e9
	# Exposure time to reach desired I2 counts
    exp_time = getEXPTime(i2counts, vmag, bmv, elevation, seeing, deckers)

    return exp_time, exp_counts, i2counts


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

def compute_sunset_n_rise(dt,horizon='0'):
    # computes time in seconds before sunset
    apf_obs = makeAPFObs(dt,horizon=horizon)
    sunset = apf_obs.next_setting(ephem.Sun())
    sunset -= ephem.Date(dt)
    sunset *= 86400.0 # convert to seconds
    
    sunrise = apf_obs.next_rising(ephem.Sun())
    sunrise -= ephem.Date(dt)
    sunrise *= 86400.0 # convert to seconds
    return sunset, sunrise

def compute_sunset(dt,horizon='0'):

    sunset, sunrise = compute_sunset_n_rise(dt,horizon=horizon)
    return sunset
    
def compute_sunrise(dt,horizon='0'):
    sunset, sunrise = compute_sunset_n_rise(dt,horizon=horizon)
    return sunrise

def smartList(starlist, ctime, seeing, slowdown,outdir = None):
    """ Determine the best target to observe from the provided scriptobs-compatible starlist.
        Here the best target is defined as an unobserved target (ie not in observed targets )
        that is visible above 30 degrees elevation. Higher elevation targets are prefered,
        but those that rise above 85 degrees will be regected to avoid slewing through the zenith. """
    # Convert the unix timestamp into a python datetime

    dt = computeDatetime(ctime)

    if not outdir:
        outdir = os.getcwd()
    observed, times, temps = ObservedLog.getObserved(os.path.join(outdir, "observed_targets"))

    apf_obs = makeAPFObs(dt)
    # APF latitude in radians
    apf_lat = apf_obs.lat

    # Parse the starlist
    try:
        sn, star_table, lines, stars = parseStarlist(starlist)
    except ValueError:
        # This will be raised when the starlist could not be parsed successfully.
        apflog( "No target could be selected because starlist could not be parsed.", level="warn", echo=True)
        return None
    targNum = len(sn)

    # Minimum Brightness based on conditions
    VMAX = 14

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    # Distance to stay away from the moon [Between 15 and 25 degrees]

    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN

    available = np.ones(targNum, dtype=bool)
    moon_check = np.where(moonDist > minMoonDist, True, False)
    available = available & moon_check

    # If seeing is bad, only observe bright targets ( Large VMAG is dim star )
    brightenough = np.where(star_table[:, DS_VMAG] < VMAX, True, False)
    available = available & brightenough

    obs_length = star_table[:,DS_EXPT] * star_table[:,DS_NSHOTS] + 45 * (star_table[:,DS_NSHOTS]-1)
    vis, star_elevations, fin_els = Visible.is_visible(apf_obs, stars, obs_length)
    available = available & vis

    done = [ True if n in observed else False for n in sn ]
    availableandnotdone = available & np.logical_not(done)

    if not any(availableandnotdone):
        apflog( "All visible targets have been observed", level="warn", echo=True)
        (good,) = np.where(available)
    else:
        (good,) = np.where(availableandnotdone)

    delta = fin_els[good] - star_elevations[good]
    neg   = np.where(delta < 0)
    pos   = np.where(delta >= 0)
    inv_els = fin_els[good]
    inv_els[neg]  = fin_els[good[neg]] - 90
    inv_els[pos] = 90 - fin_els[good[pos]]
    sort_fin_els_idx = np.argsort(inv_els)
    idx = good[sort_fin_els_idx[0]]

    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table[idx, DS_PMRA]
    res['PM_DEC'] = star_table[idx, DS_PMDEC]
    res['VMAG']   = star_table[idx, DS_VMAG]
    res['BV'] = 0.6
    res['PRI'] = 10.
    res['SCORE'] = 1.0
    res['COUNTS'] = star_table[idx, DS_COUNTS]
    res['EXP_TIME'] = star_table[idx, DS_EXPT]
    res['NEXP'] = star_table[idx, DS_NSHOTS]
    res['NAME']   = sn[idx]
    res['SCRIPTOBS'] = [lines[idx]]
    return res

def format_expmeter(exp_counts, nexp, exptime):

    nexp[nexp <= 0] = 1
    nexp[nexp > MAX_NSHOTS] = MAX_NSHOTS
    exps = np.zeros_like(exp_counts)
    exp_counts *= 1.1
    long_idx = exp_counts > MAX_EXPMETER
    toofew_idx = exp_counts/nexp > MAX_EXPMETER
    nexp[toofew_idx] = np.ceil((exp_counts[toofew_idx]/MAX_EXPMETER) + 1)
    exp_counts[long_idx] = MAX_EXPMETER
    exps[exps < nexp] = nexp[exps < nexp]
    exps[exps > MAX_NSHOTS] = MAX_NSHOTS
    exp_counts /= exps
    # this is basically so that really bright stars do not get really short exposures
    # certain targets we want to hit the exposure meter threshold per exposure not for the total
    # this violates the spirit of the scheduler but leads to consistent exposures for bright standard
    # stars like 185144 or tau Ceti

    really_short = exptime < MIN_TOTOBS
    exp_counts[really_short] = MAX_EXPMETER/2.
    return exp_counts, exps

def format_time(total, i2counts, nexp, mintime, maxtime, hitthemall=False):
    total = np.array(total)
    times = np.zeros(len(total))
    exps  = np.ones(len(total))

    middle_idx = (total > mintime ) &(total < maxtime)
    times[middle_idx] = maxtime[middle_idx] # pad out to make it more likely exposure meter threshold sets actual limit
    exps[middle_idx] = 1

    max_idx = total > maxtime
    if hitthemall:
        exps[max_idx] = 1
    else:
        exps[max_idx] = np.ceil(total[max_idx]/maxtime[max_idx])
    times[max_idx] = maxtime[max_idx]

    short_idx = total < mintime
    times[short_idx] = mintime[short_idx]  # pad out to make it more likely exposure meter threshold sets actual limit
    exps[short_idx] = np.ceil(mintime[short_idx]/(total[short_idx]+READOUT)) + 1

    # I am not sure why these are falling through
    times[times <= 0] = mintime[times <= 0]

    exps[exps < nexp] = nexp[exps < nexp]

    return times, exps

def templateConditions(moon, seeing, slowdown):

    if seeing < 20 and slowdown < 0.7:
        apflog("moon.phase=%.2f moon.alt=%.2f" % (moon.phase,moon.alt),echo=True,level='debug')
        if moon.phase < 50 and float(moon.alt) < 0.7:
            return True
        else:
            return False
    else:
        return False

def findClosest(ras,decs,ra,dec):


    distances = np.sqrt((ra - ras)**2 + (dec - decs)**2)

    min_ind = distances.argmin()

    return min_ind

def makeTempRow(star_table,ind,bstar=False):

    row = []

    row.append(star_table[ind, DS_RA])
    row.append( star_table[ind, DS_DEC])
    row.append(star_table[ind, DS_PMRA])
    row.append(star_table[ind, DS_PMDEC])
    row.append(star_table[ind, DS_VMAG])
    row.append(1200)
    row.append(1e9)
    row.append( star_table[ind, DS_APFPRI])
    row.append(0)
    if bstar:
        row.append(2)
    else:
        if star_table[ind, DS_VMAG] > 10:
            row.append(7)
        elif star_table[ind, DS_VMAG] < 8:
            row.append(3)
        else:
            row.append(5)
    return row

def enoughTime(star_table,stars,idx,row,apf_obs,dt):
    tot_time = row[DS_NSHOTS]*row[DS_EXPT]
    tot_time += 210 + (2*40 + 40*(row[DS_NSHOTS]-1)) # two B star exposures + three 70 second acquisitions and the actual observation readout times
    vis, star_elevations, fin_els = Visible.is_visible(apf_obs,[stars[idx]],[tot_time])
    time_left_before_sunrise = compute_sunrise(dt,horizon='-9')

    try:
        apflog( "enoughTime(): time for obs= %.1f  time until sunrise= %.1f " % (tot_time,  time_left_before_sunrise),echo=True)
    except:
        apflog("enoughTime(): cannot log times!?!",echo=True)
        
    if tot_time < time_left_before_sunrise  and vis and time_left_before_sunrise < 14*3600.:
        return True
    else:
        return False
        
    
def findBstars(snames,star_table,idx, bstars):

    near_idx = findClosest(star_table[:,DS_RA][bstars],star_table[:,DS_DEC][bstars],star_table[idx,DS_RA],star_table[idx,DS_DEC])
    row = makeTempRow(star_table[bstars],near_idx,bstar=True)
    
    end_idx = findClosest(star_table[:,DS_RA][bstars],star_table[:,DS_DEC][bstars],(star_table[idx,DS_RA]+15*np.pi/180.),star_table[idx,DS_DEC])
    finrow = makeTempRow(star_table[bstars],end_idx,bstar=True)
    
    return snames[near_idx],row,snames[end_idx],finrow


def makeResult(stars,star_table,flags,totexptimes,i2cnts,sn,dt,idx,focval=0):
    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table[idx, DS_PMRA]
    res['PM_DEC'] = star_table[idx, DS_PMDEC]
    res['VMAG']   = star_table[idx, DS_VMAG]
    res['BV']     = star_table[idx, DS_BV]
    res['COUNTS'] = star_table[idx, DS_COUNTS]
    res['EXP_TIME'] = star_table[idx, DS_EXPT]
    res['NEXP'] = star_table[idx, DS_NSHOTS]
    res['TOTEXP_TIME'] = totexptimes[idx]
    res['I2CNTS'] = i2cnts[idx]
    res['NAME']   = sn[idx]
    res['SCORE']  = star_table[idx,DS_NSHOTS]
    res['PRI']    = star_table[idx, DS_APFPRI]
    res['DECKER'] = flags['decker'][idx]
    res['I2'] =    flags['I2'][idx]
    res['isTemp'] =    False
    res['owner'] =    flags['owner'][idx]
    res['SCRIPTOBS'] = []
    scriptobs_line = makeScriptobsLine(sn[idx], star_table[idx,:], flags['do'][idx], dt, decker=flags['decker'][idx], I2=flags['I2'][idx], owner=flags['owner'][idx],focval=focval) + " # end"
    res['SCRIPTOBS'].append(scriptobs_line)
    return res


def getNext(ctime, seeing, slowdown, bstar=False,template=False,sheetns=["Bstars"],owner='public',outfn="googledex.dat",toofn="too.dat",outdir=None,focval=0):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """

    if not outdir:
        outdir = os.getcwd()

    dt = computeDatetime(ctime)

    confg = dict()
    confg['I2'] = 'Y'
    confg['decker']='W'


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

    observed = ParseGoogledex.updateLocalGoogledex(ptime,googledex_file=os.path.join(outdir,"googledex.dat"), observed_file=os.path.join(outdir,"observed_targets"))

    # List of targets already observed

    global last_objs_attempted
    try:
        lastline = ktl.read("apftask","SCRIPTOBS_LINE")
        if not bstar:             # otherwise from previous night
            lastobj = lastline.split()[0]
        else:
            lastobj = None

    except:
        lastobj = None

    if lastobj:
        if lastobj not in observed and lastobj not in last_objs_attempted:
            last_objs_attempted.append(lastobj)
            
            apflog( "getNext(): Last objects attempted %s" % (last_objs_attempted),echo=True)

            if len(last_objs_attempted) > 5:
                apflog( "getNext(): 5 failed acquisition attempts",echo=True)
                last_objs_attempted = []
                return None
        else:
            last_objs_attempted = []
            # we had a succes so we are zeroing this out

            
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

    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians

    apflog("getNext(): Parsing the Googledex...",echo=True)
    config={'I2': 'Y', 'decker': 'W', 'owner' : owner}
    sn, star_table, flags, stars = ParseGoogledex.parseGoogledex(sheetns=sheetns,outfn=outfn,outdir=outdir,config=config)
    sn = np.array(sn)
    deckers = np.array(flags['decker'])
    targNum = len(sn)
    
    apflog("getNext(): Parsed the Googledex...",echo=True)

    # Note which of these are B-Stars for later.
    bstars = np.array([ True if 'HR' in n else False for n in sn ], dtype=bool)

    apflog("getNext(): Finding B stars",echo=True)


    # Distance to stay away from the moon
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN

    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))

    available = np.ones(targNum, dtype=bool)
    totexptimes = np.zeros(targNum, dtype=float)
    cur_elevations = np.zeros(targNum, dtype=float)
    scaled_elevations = np.zeros(targNum, dtype=float)
    i2cnts = np.zeros(targNum, dtype=float)

    # Is the target behind the moon?

    apflog("getNext(): Culling stars behind the moon",echo=True)
    moon_check = moonDist > minMoonDist
    available = available & moon_check

    #    totobs_check = (star_table[:,DS_NOB] < star_table[:,DS_TOT]) | (star_table[:,DS_TOT] <= 0)
    #    available = available & totobs_check

    # We just need a B star, so restrict our math to those
    if bstar:
        
        apflog("getNext(): Selecting B stars",echo=True)
        available = available & bstars

        f = available
        
        apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations = Visible.is_visible(apf_obs, fstars, [400]*len(bstars[f]))

        available[f] = available[f] & vis
        cur_elevations[np.where(f)] += star_elevations[np.where(vis)]

        star_table[available, DS_COUNTS] = 2e9
        star_table[available, DS_EXPT] = 900
        star_table[available, DS_NSHOTS] = 2
        totexptimes[available] = 400

    # Just need a normal star for observing
    else:
        # Available and not a BStar

        apflog("getNext(): Culling B stars",echo=True)
        available = np.logical_and(available, np.logical_not(bstars))

        # has the star been observed - commented out as redundant with cadence
        if len(last_objs_attempted)>0:
            for n in last_objs_attempted:
                attempted = (sn == n)
                available = available & np.logical_not(attempted) # Available and not observed

        # Calculate the exposure time for the target
        # Want to pass the entire list of targets to this function
        f = available

        apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations, scaled_els = Visible.is_visible_se(apf_obs, fstars, [0]*len(fstars))
#        vis,star_elevations,fin_star_elevations = Visible.is_visible( apf_obs,fstars,[0]*len(fstars))
        available[f] = available[f] & vis
        f = available
        fstars = [s for s,_ in zip(stars,f) if _ ]

        apflog("getNext(): Computing exposure times",echo=True)
        exp_times, exp_counts, i2counts = calculateUCSCExposureTime( star_table[f,DS_VMAG], \
                                            star_table[f,DS_I2CNTS], star_elevations[np.array(vis)], seeing, \
                                            star_table[f,DS_BV], deckers[f])

        exp_times = exp_times * slowdown
        totexptimes[f] += computeMaxTimes(exp_times,star_table[f, DS_MAX])
        i2cnts[f] += i2counts

        apflog("getNext(): Formating exposure times",echo=True)
        mxtime = np.zeros_like(star_table[f,DS_MAX])
        mxtime += MAX_EXPTIME
        shorter = (star_table[f,DS_MAX] < MAX_EXPTIME)&(star_table[f,DS_MAX] >0)
        mxtime[shorter] = star_table[f,DS_MAX][shorter]
        star_table[f, DS_EXPT], exps = format_time(totexptimes[f],i2counts,star_table[f, DS_NSHOTS],star_table[f, DS_MIN],mxtime)


        apflog("getNext(): Formating exposure meter",echo=True)
        star_table[f, DS_COUNTS], star_table[f, DS_NSHOTS] = format_expmeter(exp_counts,exps, totexptimes[f])

        # Is the exposure time too long?

        apflog("getNext(): Removing really long exposures",echo=True)
        time_check = np.where( exp_times < TARGET_EXPOSURE_TIME_MAX, True, False)

        available[f] = available[f] & time_check
        f = available

        # Is the star currently visible?
        apflog("getNext(): Computing stars visibility",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations, scaled_els = Visible.is_visible_se(apf_obs, fstars, exp_times)
        if vis != []:
            available[f] = available[f] & vis
        cur_elevations[np.where(f)] += star_elevations[np.where(vis)]
        scaled_elevations[np.where(f)] += scaled_els[np.where(vis)]


    # Now just sort by priority, then cadence. Return top target
    if len(sn[available]) < 1:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None


    final_priorities = computePriorities(star_table,available,dt,flags)

    cadence_check = (ephem.julian_date(dt) - star_table[:, DS_LAST]) / star_table[:, DS_CAD]
    good_cadence = np.where(cadence_check >  1.0, True, False)
    good_cadence_available = available & good_cadence

    if any(good_cadence_available):
        try:
            pri = max(final_priorities[good_cadence_available])
            sort_i = (final_priorities == pri) & good_cadence_available
        except:
            pri = max(final_priorities[available])
            sort_i = (final_priorities == pri) & available
    elif any(available):
        apflog( "getNext(): No new stars available, going back to the previously observed list.",level="warn",echo=True)
        pri = max(final_priorities[available])
        sort_i = (final_priorities == pri) & available
    else:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None

    starstr = "getNext(): star table available: %s" % (sn[sort_i])
    apflog(starstr,echo=True)

    starstr = "getNext(): star table available priorities: %s" % (final_priorities[sort_i])
    apflog(starstr,echo=True)

    if bstar:
        sort_j = cur_elevations[sort_i].argsort()[::-1]
        focval=2
    else:
        sort_j = scaled_elevations[sort_i].argsort()[::-1]
        cstr= "getNext(): cadence check: %s" % (cadence_check[sort_i][sort_j][0])
        apflog(cstr,echo=True)

    t_n = sn[sort_i][sort_j][0]

    elstr= "getNext(): star elevations %s" % (cur_elevations[sort_i][sort_j])
    apflog(elstr,echo=True)

    t_n = sn[sort_i][sort_j][0]

    apflog("getNext(): selected target %s" % (t_n) )

    idx, = np.where(sn == t_n)
    idx = idx[0]

    stars[idx].compute(apf_obs)
    cstr= "getNext(): cadence check: %f (%f %f %f)" % (((ephem.julian_date(dt) - star_table[idx, DS_LAST]) / star_table[idx, DS_CAD]), ephem.julian_date(dt), star_table[idx, DS_LAST], star_table[idx, DS_CAD])
    apflog(cstr,echo=True)

    res =  makeResult(stars,star_table,flags,totexptimes,i2cnts,sn,dt,idx,focval=focval)
    if do_templates and flags['template'][idx] == 'N' and flags['I2'][idx] == 'Y':
        bname,brow,bnamefin,browfin = findBstars(sn,star_table,idx,bstars)
        row = makeTempRow(star_table,idx)
        if enoughTime(star_table,stars,idx,row,apf_obs,dt):
            bline = makeScriptobsLine(bname,brow,'',dt,decker="N",I2="Y", owner='public',focval=2)
            line  = makeScriptobsLine(sn[idx],row,flags['do'][idx],dt,decker="N",I2="N", owner=flags['owner'][idx])
            bfinline = makeScriptobsLine(bnamefin,browfin,'',dt,decker="N",I2="Y", owner=flags['owner'][idx],focval=2)
            res['SCRIPTOBS'] = []
            res['SCRIPTOBS'].append(bfinline + " # temp=Y end")
            res['SCRIPTOBS'].append(line + " # temp=Y")
            res['SCRIPTOBS'].append(bline + " # temp=Y")
            res['isTemp'] = True
            apflog("Attempting template observation of %s" % (sn[idx]),echo=True)

    return res

if __name__ == '__main__':

#    sheetn=["2018B"]
    sheetn="Bstars,A003_PRobertson_2019B,A006_PDalba_2019B,A007_HIsaacson_2019B,A009_MKosiarek_2019B,A011_SKane_2019B,A012_SKane_2019B,A015_AHoward_2019B,A013_ASiemion_2019B,A000_BWelsh_2019B,A001_ICzekala_2019B,A002_ICzekala_2019B,A004_PRobertson_2019B,A007_HIsaacson_2019B,A008_BHolden_2019B,A014_SVogt_2019B,A015_TBrandt_2019B"

    # For some test input what would the best target be?
    otfn = "observed_targets"
    ot = open(otfn,"w")
    starttime = time.time()
    result = getNext(starttime, 7.99, 0.4, bstar=True,sheetns=sheetn.split(","))
    ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
    ot.close()
    starttime += 400
    for i in range(5):

        result = getNext(starttime, 7.99, 0.4, bstar=False,sheetns=sheetn,template=True)
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
