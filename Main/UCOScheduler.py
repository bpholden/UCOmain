# UCOScheduler_V1.py
from __future__ import print_function
import os

import time
import datetime

import numpy as np
import ephem

import ParseUCOSched
import SchedulerConsts
import Observability
import ScriptobsLine
import SunPos
import UCOTargets
import Visible

try:
    from apflog import apflog
    import ktl
except:
    from fake_apflog import apflog

# a global
last_objs_attempted = []

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

def compute_priorities(star_table, cur_dt, observed=None, hour_table=None, rank_table=None, do_templates=False):
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
    really_bad_candence = cadence_check < .7

    started_doubles = star_table['night_cad'] > 0
    started_doubles = started_doubles & (star_table['night_obs'] > 0)
    started_doubles = started_doubles & (star_table['night_obs'] < star_table['night_nexp'])
    if np.any(started_doubles):
        redo = started_doubles & (cadence_check > (star_table['night_cad'] - SchedulerConsts.BUFFER))
        redo = redo & (cadence_check < (star_table['night_cad'] + SchedulerConsts.BUFFER))
    else:
        redo = np.zeros(1,dtype=bool)

    done_sheets = False
    if hour_table is not None:
        too_much = hour_table['cur']  > hour_table['tot']
        done_sheets = hour_table['sheetn'][too_much]

    if done_sheets is not False:
        done_sheets_str = " ".join(list(done_sheets))
        apflog("The following sheets are finished for the night: %s " %
               (done_sheets_str), echo=True)

    bad_pri = np.ones_like(star_table['pri'])

    if rank_table is not None:
        for sheetn in rank_table['sheetn']:
            if sheetn not in done_sheets:
                cur = star_table['sheetn'] == sheetn
                new_pri[cur & good_cadence] += rank_table['rank'][rank_table['sheetn'] == sheetn]
                new_pri[cur & bad_cadence] += 2*bad_pri[(cur & bad_cadence)]
                new_pri[cur & really_bad_candence] += bad_pri[(cur & really_bad_candence)]
            else:
                cur = star_table['sheetn'] == sheetn
                new_pri[cur & good_cadence] += 100
                new_pri[cur & bad_cadence] += 2*bad_pri[(cur & bad_cadence)]
                new_pri[cur & really_bad_candence] += bad_pri[(cur & really_bad_candence)]

    done_all = (star_table['nobs'] >= star_table['totobs']) & (star_table['totobs'] > 0)
    new_pri[done_all] = 0

    if np.any(redo):
        new_pri[redo] = np.max(rank_table['rank'])

    new_pri[star_table['cal_star'] == 'Y'] = 0

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
    res['isTOO'] = star_table['too'][idx]
    res['mode'] = ''
    res['owner'] = star_table['sheetn'][idx]

#    if np.ma.is_masked(star_table[idx]['obsblock']):
#        res['obsblock'] = ''

    res['SCRIPTOBS'] = []
    if bstar:
        scriptobs_line = ScriptobsLine.make_scriptobs_line(star_table[idx], dt, decker=res['DECKER'], \
                                         owner=res['owner'], I2='N', \
                                            focval=0)
        # we hard code the focval to skip it because 
        # we will focus on the observation created 
        # in the line below
        scriptobs_line = scriptobs_line + " # end"
        res['SCRIPTOBS'].append(scriptobs_line)

    scriptobs_line = ScriptobsLine.make_scriptobs_line(star_table[idx], dt, decker=res['DECKER'], \
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


def get_next(ctime, seeing, slowdown, ucotargets, \
                bstar=False, do_templates=False, \
                do_too=False, owner='public', \
                outfn="googledex.dat", toofn="too.dat", \
                outdir=None, focval=0, inst='', \
                start_time=None):
    """ Determine the best target to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """

    global last_objs_attempted

    if not outdir:
        outdir = os.getcwd()

    dt = Observability.compute_datetime(ctime)

    config = ScriptobsLine.config_defaults(owner)

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
    observed, ucotargets.star_table = ParseUCOSched.update_local_starlist(ptime,\
                                                               outfn=outfn, toofn=toofn, \
                                                                observed_file="observed_targets")

    ucotargets.make_hour_table()

    if ucotargets.hour_table is not None:
        ucotargets.hour_table = update_hour_table(ucotargets.hour_table, observed, ptime)
    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians

    if ucotargets.star_table is None:
        apflog("get_next(): Parsing the star list", echo=True)
        ucotargets.make_star_table()
    ucotargets.append_too_column()

    stars = ParseUCOSched.gen_stars(ucotargets.star_table)
    targ_num = len(stars)

    last_failure = last_attempted()
    if last_failure is not None:
        last_objs_attempted.append(last_failure)

    ###
    # Need to update the googledex with the lastObserved date for observed targets
    # Scriptobs line uth utm can be used for this
    # Need to convert a uth and utm to a JD quickly.
    # timedelta = now - uth,utm : minus current JD?
    ###

    apf_obs = SunPos.make_APF_obs(dt)

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    template_conditions_met = Observability.template_conditions(moon, seeing, slowdown)
    do_templates = do_templates and template_conditions_met

    apflog("get_next(): Will attempt templates = %s" % str(do_templates) ,echo=True)
    # Note which of these are B-Stars for later.
    bstars = (ucotargets.star_table['Bstar'] == 'Y')|(ucotargets.star_table['Bstar'] == 'y')

    if bstar and np.any(bstars) is False:
        apflog("get_next(): No B stars listed in target sheets!", label='Error', echo=True)
        return None

    apflog("get_next(): Computing exposure times", echo=True)
    totexptimes = Observability.tot_exp_times(ucotargets.star_table, targ_num)

    available = np.ones(targ_num, dtype=bool)
    cur_elevations = np.zeros(targ_num, dtype=float)
    scaled_elevations = np.zeros(targ_num, dtype=float)

    # Is the target behind the moon?

    moon_check = Observability.behind_moon(moon, ucotargets.star_table['ra'], ucotargets.star_table['dec'])
    available = available & moon_check
    log_str = "get_next(): Moon visibility check - stars rejected = "
    log_str += "%s" % ( np.asarray(ucotargets.star_table['name'][np.logical_not(moon_check)]))
    apflog(log_str, echo=True)

    sun_el_good = SunPos.sun_el_check(ucotargets.star_table, apf_obs, horizon='-18')
    available = available & sun_el_good

    # other condition cuts (seeing, transparency, moon phase)
    cuts = Observability.condition_cuts(moon, seeing, slowdown, ucotargets.star_table)
    available = available & cuts

    if len(last_objs_attempted)>0:
        for n in last_objs_attempted:
            attempted = ucotargets.star_table['name'] == n
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

    if do_too is False:
        apflog("get_next(): Selecting TOO targets", echo=True)
        not_too = ucotargets.star_table['too'] == False
        available = available & not_too

    # Is the exposure time too long?
    apflog("get_next(): Removing really long exposures", echo=True)
    time_good = Observability.time_check(ucotargets.star_table, totexptimes, dt, start_time=start_time)

    available = available & time_good
    if np.any(available) is False:
        apflog( "get_next(): Not enough time left to observe any targets", level="error", echo=True)
        return None

    # Compute the elevations of the stars

    apflog("get_next(): Computing star elevations",echo=True)
    fstars = [s for s,_ in zip(stars,available) if _ ]
    vis, star_elevations, scaled_els = Visible.visible(apf_obs, fstars, \
                                                       totexptimes[available],
                                                       shiftwest=shiftwest
    )

    currently_available = available
    if len(star_elevations) > 0:
        currently_available[available] = currently_available[available] & vis
    else:
        apflog( "get_next(): Couldn't find any suitable targets!", level="error", echo=True)
        return None

    cur_elevations[available] += star_elevations[vis]
    scaled_elevations[available] += scaled_els[vis]

    if slowdown > SchedulerConsts.SLOWDOWN_THRESH or seeing > SchedulerConsts.SEEING_THRESH:
        bright_enough = ucotargets.star_table['Vmag'] < SchedulerConsts.SLOWDOWN_VMAG_LIM
        available = available & bright_enough

    if not do_templates:
        available = available & (ucotargets.star_table['only_template'] == 'N')
    # Now just sort by priority, then cadence. Return top target
    if len(ucotargets.star_table['name'][available]) < 1:
        apflog( "get_next(): Couldn't find any suitable targets!", level="error", echo=True)
        return None

    final_priorities = compute_priorities(ucotargets.star_table, dt,
                                             rank_table=ucotargets.rank_table,
                                             hour_table=ucotargets.hour_table,
                                             do_templates=do_templates,
                                             observed=observed)

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

    t_n = ucotargets.star_table['name'][idx]
    o_n = ucotargets.star_table['sheetn'][idx]
    p_n = final_priorities[idx]

    apflog("get_next(): selected target %s for program %s at priority %.0f" % (t_n, o_n, p_n) )
    nmstr= "get_next(): star names %s" % (np.asarray(ucotargets.star_table['name'][sort_i][sort_j]))
    pristr= "get_next(): star priorities %s" % (np.asarray(final_priorities[sort_i][sort_j]))
    mxpristr= "get_next(): max priority %d" % (pri)
    shstr= "get_next(): star sheet names %s" % (np.asarray(ucotargets.star_table['sheetn'][sort_i][sort_j]))
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

    take_template = do_templates and ucotargets.star_table['Template'][idx] == 'N' \
        and ucotargets.star_table['I2'][idx] == 'Y'
    if ucotargets.star_table['only_template'][idx] == 'Y' and do_templates:
        take_template = True

    res =  make_result(stars, ucotargets.star_table, totexptimes, final_priorities, dt, \
                       idx, focval=focval, bstar=bstar, mode=config['mode'])
    if take_template and bstar is False:
        bidx, bfinidx = Observability.find_Bstars(ucotargets.star_table, idx, bstars)

        if Observability.enough_time_templates(ucotargets.star_table,stars,idx,apf_obs,dt):
            decker= "N"
            line  = ScriptobsLine.make_scriptobs_line(ucotargets.star_table[idx], \
                                        dt, decker=decker, I2="N", owner=res['owner'], temp=True)
            if "decker=W" in line:
                decker = "W"
            bline = ScriptobsLine.make_scriptobs_line(ucotargets.star_table[bstars][bidx], dt, \
                                        decker=decker, I2="Y", owner=res['owner'], focval=2)
            bfinline = ScriptobsLine.make_scriptobs_line(ucotargets.star_table[bstars][bfinidx], dt,\
                                            decker=decker, I2="Y", owner=res['owner'], focval=0)
            res['SCRIPTOBS'] = []
            res['SCRIPTOBS'].append(bfinline + " # temp=Y end")
            res['SCRIPTOBS'].append(line + " # temp=Y")
            res['SCRIPTOBS'].append(bline + " # temp=Y")
            res['isTemp'] = True
            res['DECKER'] = decker
            apflog("Attempting template observation of %s" % (ucotargets.star_table['name'][idx]), echo=True)

    res['template_conditions_met'] = template_conditions_met
    return res

def test_basic_ops(ucotargets):
    """
    test_basic_ops()
    """

    # Test the basic operations of the scheduler
    # This is a test function to see if the basic operations work
    # It will not be run in production

    try:
        ktl.write('apftask', 'SCRIPTOBS_LINE_RESULT', 3, binary=True)
    except:
        pass

    # For some test input what would the best target be?
    OTFN = "observed_targets"
    ot = open(OTFN, "w")
    starttime = time.time()
    result = get_next(starttime, 7.99, 0.4, ucotargets, bstar=True, \
                      do_templates=False)
    while len(result['SCRIPTOBS']) > 0:
        ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
    ot.close()

    for i in range(5):

        result = get_next(starttime, 7.99, 0.4, ucotargets, bstar=False, \
                         do_templates=False)
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

    return starttime

def test_failure(starttime, ucotargets):
    '''
    test_failure(starttime, ucotargets)
    starttime - time to start the test
    ucotargets - UCOTargets object
    '''
    print("Testing a failure")
    try:
        ktl.write('apftask', 'SCRIPTOBS_LINE_RESULT', 2, binary=True)
    except:
        pass
    result = get_next(starttime, 7.99, 0.4, ucotargets, bstar=False, \
                     do_templates=True, )
    print(result)
    print("Nonsensical start time")
    result = get_next(starttime, 7.99, 0.4, ucotargets, bstar=True, \
                     do_templates=True, start_time=1)
    print(result)
    return

def test_templates(ucotargets):
    """
    test_templates(ucotargets)
    ucotargets - UCOTargets object
    """
    print("Testing templates")
    t_dt = datetime.datetime.now()
    tstar_table, _ = ParseUCOSched.parse_UCOSched(ucotargets.rank_table, \
                                                     outfn='googledex.dat', outdir=".", \
                                                        config=ScriptobsLine.config_defaults('public'))
    tidx, = np.asarray(tstar_table['name'] == '185144').nonzero()
    tidx = tidx[0]
    tbstars = (tstar_table['Bstar'] == 'Y')|(tstar_table['Bstar'] == 'y')
    tbidx, tbfinidx = Observability.find_Bstars(tstar_table, tidx, tbstars)
    decker = "N"
    tline  = ScriptobsLine.make_scriptobs_line(tstar_table[tidx], t_dt, \
                                decker=decker, I2="N", owner='public', temp=True)
    if "decker=W" in tline:
        decker = "W"
    tbline = ScriptobsLine.make_scriptobs_line(tstar_table[tbstars][tbidx], t_dt, \
                                decker=decker, I2="Y", owner='public', focval=2)

    tbfinline = ScriptobsLine.make_scriptobs_line(tstar_table[tbstars][tbfinidx], t_dt, \
                                   decker=decker, I2="Y", owner='public', focval=0)
    temp_res= []
    temp_res.append(tbfinline + " # temp=Y end")
    temp_res.append(tline + " # temp=Y")
    temp_res.append(tbline + " # temp=Y")
    out_r = [print(r) for r in temp_res]
    print(out_r)
    print("Done")

def test_main():
    """
    test_main()
    """

    # Test the basic operations of the scheduler
    # This is a test function to see if the basic operations work
    # It will not be run in production

    RANK_TABLEN='2025B_ranks_operational'

    class Opt:
        def __init__(self):
            self.test = True
            self.time_left = "/home/holden/time_left.csv"
            self.rank_table = RANK_TABLEN

    uco_targets = UCOTargets.UCOTargets(Opt())

    # this calls make_rank_table
    uco_targets.make_hour_constraints()
    # this calls make_hour_table 
    uco_targets.make_hour_table()

    starttime = test_basic_ops(uco_targets)

    test_failure(starttime, uco_targets)
    test_templates(uco_targets)

if __name__ == '__main__':

    test_main()
