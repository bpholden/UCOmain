# UCOScheduler_V1.py
from __future__ import print_function
import os

import time
import datetime

import numpy as np
import ephem

import SchedulerConsts
import Observability
import ScriptobsLine
import SunPos
import Target
import Visible

try:
    from apflog import apflog
    import ktl
except:
    from fake_apflog import apflog

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

    res - Target
    '''
    res = Target.Target.from_star_table(star_table, idx, stars[idx], totexptimes[idx],
                                        final_priorities[idx], bstar=bstar)

#    if np.ma.is_masked(star_table[idx]['obsblock']):
#        res['obsblock'] = ''

    if bstar:
        scriptobs_line = ScriptobsLine.make_scriptobs_line(star_table[idx], dt, decker=res.decker, \
                                         owner=res.owner, I2='N', \
                                            focval=0)
        # we hard code the focval to skip it because 
        # we will focus on the observation created 
        # in the line below
        scriptobs_line = scriptobs_line + " # end"
        res.scriptobs.append(scriptobs_line)

    scriptobs_line = ScriptobsLine.make_scriptobs_line(star_table[idx], dt, decker=res.decker, \
                                         owner=res.owner, I2=star_table['I2'][idx], \
                                            focval=focval)

    scriptobs_line = scriptobs_line + " # end"
    res.scriptobs.append(scriptobs_line)
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


class UCOScheduler(object):
    '''
    Selects the next target to observe from a UCOTargetTables object.
    Holds the run-state that must survive between calls, such as the
    list of objects that recently failed to be observed.

    track_failures - if True, get_next records the last object attempted
    when it failed and skips it on later calls, until
    zero_last_objs_attempted is called.

    tot_temps - the most template observations get_next will return over
    the scheduler's lifetime; None means no limit.

    do_too - whether ToO targets may be selected. get_next turns it off
    after returning a ToO, so one is taken per time it is turned on.

    '''
    def __init__(self, targets, owner='public', outdir=None,
                 do_templates=False, do_too=False, start_time=None,
                 outfn='googledex.dat', toofn='too.dat', track_failures=False,
                 tot_temps=None):
        self.targets = targets
        self.owner = owner
        self.outdir = outdir or os.getcwd()
        self.outfn = outfn
        self.toofn = toofn
        self.do_templates = do_templates
        self.do_too = do_too
        self.start_time = start_time
        self.track_failures = track_failures
        self.tot_temps = tot_temps

        # run-state, was a module global
        self.last_objs_attempted = []
        # templates returned so far, counted against tot_temps
        self.n_temps = 0

        # per-call scratch, kept for logging and for Observe to inspect
        self.observed = None
        self.apf_obs = None
        self.moon = None
        self.stars = None
        self.result = None
        self.template_conditions_met = False

    def __repr__(self):
        return "<UCOScheduler targets=%s owner=%s>" % (self.targets, self.owner)

    def zero_last_objs_attempted(self):
        self.last_objs_attempted = []

    def record_last_attempt(self):
        last_failure = last_attempted()
        if last_failure is not None:
            self.last_objs_attempted.append(last_failure)

    def get_next(self, ctime, seeing, slowdown, bstar=False, focval=0,
                 do_templates=None, do_too=None):
        """ Determine the best target to observe for the given input.
            Takes the time, seeing, and slowdown factor.
            Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
        """
        if do_templates is None:
            do_templates = self.do_templates
        if self.tot_temps is not None and self.n_temps >= self.tot_temps:
            do_templates = False
        if do_too is None:
            do_too = self.do_too

        dt = Observability.compute_datetime(ctime)

        config = ScriptobsLine.config_defaults(self.owner)

        apflog( "get_next(): Finding target for time %s" % (dt), echo=True)

        if slowdown > SchedulerConsts.SLOWDOWN_MAX:
            log_str = "get_next(): Slowndown value of %f " % (slowdown)
            log_str += "exceeds maximum of %f at time %s" % (SchedulerConsts.SLOWDOWN_MAX, dt)
            apflog(log_str , echo=True)
            return None

        ptime = self._previous_obs_time(dt)
        self._refresh_tables(ptime)
        star_table = self.targets.star_table
        stars = self.stars
        targ_num = len(stars)

        if self.track_failures:
            self.record_last_attempt()

        self._sky_state(dt)

        self.template_conditions_met = Observability.template_conditions(self.moon, seeing, slowdown)
        do_templates = do_templates and self.template_conditions_met

        apflog("get_next(): Will attempt templates = %s" % str(do_templates) ,echo=True)
        # Note which of these are B-Stars for later.
        bstars = (star_table['Bstar'] == 'Y')|(star_table['Bstar'] == 'y')

        if bstar and not np.any(bstars):
            apflog("get_next(): No B stars listed in target sheets!", level='error', echo=True)
            return None

        apflog("get_next(): Computing exposure times", echo=True)
        totexptimes = Observability.tot_exp_times(star_table, targ_num)

        found = self._available(dt, seeing, slowdown, bstar, bstars, totexptimes,
                                do_too, do_templates)
        if found is None:
            return None
        available, cur_elevations, scaled_elevations = found

        final_priorities = compute_priorities(star_table, dt,
                                                 rank_table=self.targets.rank_table,
                                                 hour_table=self.targets.hour_table,
                                                 do_templates=do_templates,
                                                 observed=self.observed)

        idx = self._select(available, final_priorities, bstar, cur_elevations, scaled_elevations)
        if idx is None:
            return None
        if bstar:
            focval = 2

        stars[idx].compute(self.apf_obs)

        take_template = do_templates and star_table['Template'][idx] == 'N' \
            and star_table['I2'][idx] == 'Y'
        if star_table['only_template'][idx] == 'Y' and do_templates:
            take_template = True

        res =  make_result(stars, star_table, totexptimes, final_priorities, dt, \
                           idx, focval=focval, bstar=bstar, mode=config['mode'])
        if take_template and bstar is False:
            self._add_template(res, idx, dt, bstars)
        if res.is_temp:
            self.n_temps += 1
        if res.is_too:
            self.do_too = False

        res.template_conditions_met = self.template_conditions_met
        self.result = res
        return res

    def _previous_obs_time(self, dt):
        '''
        Time of the previous observation, from the guider, falling back to dt.

        '''
        try:
            apfguide = ktl.Service('apfguide')
            stamp = apfguide['midptfin'].read(binary=True)
            ptime = datetime.datetime.utcfromtimestamp(stamp)
        except:
            if type(dt) == datetime.datetime:
                ptime = dt
            else:
                ptime = datetime.datetime.utcfromtimestamp(int(time.time()))
        return ptime

    def _refresh_tables(self, ptime):
        '''
        Fold the observed log into the tables and regenerate the ephem objects.

        '''
        apflog("get_next(): Updating star list with previous observations", echo=True)
        self.observed = self.targets.update_from_observed(ptime, outfn=self.outfn, toofn=self.toofn)

        self.targets.make_hour_table()

        self.targets.update_hour_table(self.observed, ptime)
        # Parse the Googledex
        # Note -- RA and Dec are returned in Radians

        if self.targets.star_table is None:
            apflog("get_next(): Parsing the star list", echo=True)
            self.targets.make_star_table()
        self.targets.append_too_column()

        self.stars = self.targets.gen_stars()

    def _sky_state(self, dt):
        '''
        Set the observatory and the moon for dt.

        '''
        self.apf_obs = SunPos.make_APF_obs(dt)

        # Calculate the moon's location
        self.moon = ephem.Moon()
        self.moon.compute(self.apf_obs)

    def _available(self, dt, seeing, slowdown, bstar, bstars, totexptimes, do_too, do_templates):
        '''
        Apply the visibility and condition cuts.
        Returns (available, cur_elevations, scaled_elevations), or None
        if no target survives.

        '''
        star_table = self.targets.star_table
        targ_num = len(self.stars)

        available = np.ones(targ_num, dtype=bool)
        cur_elevations = np.zeros(targ_num, dtype=float)
        scaled_elevations = np.zeros(targ_num, dtype=float)

        # Is the target behind the moon?

        moon_check = Observability.behind_moon(self.moon, star_table['ra'], star_table['dec'])
        available = available & moon_check
        log_str = "get_next(): Moon visibility check - stars rejected = "
        log_str += "%s" % ( np.asarray(star_table['name'][np.logical_not(moon_check)]))
        apflog(log_str, echo=True)

        sun_el_good = SunPos.sun_el_check(star_table, self.apf_obs, horizon='-18')
        available = available & sun_el_good

        # other condition cuts (seeing, transparency, moon phase)
        cuts = Observability.condition_cuts(self.moon, seeing, slowdown, star_table)
        available = available & cuts

        if len(self.last_objs_attempted)>0:
            for n in self.last_objs_attempted:
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

        if do_too is False:
            apflog("get_next(): Selecting TOO targets", echo=True)
            not_too = star_table['too'] == False
            available = available & not_too

        # Is the exposure time too long?
        apflog("get_next(): Removing really long exposures", echo=True)
        time_good = Observability.time_check(star_table, totexptimes, dt, start_time=self.start_time)

        available = available & time_good
        if not np.any(available):
            apflog( "get_next(): Not enough time left to observe any targets", level="error", echo=True)
            return None

        # Compute the elevations of the stars

        apflog("get_next(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(self.stars,available) if _ ]
        vis, star_elevations, scaled_els = Visible.visible(self.apf_obs, fstars, \
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
            bright_enough = star_table['Vmag'] < SchedulerConsts.SLOWDOWN_VMAG_LIM
            available = available & bright_enough

        if not do_templates:
            available = available & (star_table['only_template'] == 'N')
        # Now just sort by priority, then cadence. Return top target
        if len(star_table['name'][available]) < 1:
            apflog( "get_next(): Couldn't find any suitable targets!", level="error", echo=True)
            return None

        return available, cur_elevations, scaled_elevations

    def _select(self, available, final_priorities, bstar, cur_elevations, scaled_elevations):
        '''
        Pick the highest priority available target, breaking ties on elevation.
        Returns the index into the star table, or None.

        '''
        star_table = self.targets.star_table

        try:
            pri = max(final_priorities[available])
            sort_i = (final_priorities == pri) & available
        except:
            apflog( "get_next(): Couldn't find any suitable targets!", level="error", echo=True)
            return None

        if bstar:
            sort_j = cur_elevations[sort_i].argsort()[::-1]
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

        return idx

    def _add_template(self, res, idx, dt, bstars):
        '''
        Replace the scriptobs lines in res with a template sequence,
        if there is enough time for one.

        '''
        star_table = self.targets.star_table
        bidx, bfinidx = Observability.find_Bstars(star_table, idx, bstars)

        if Observability.enough_time_templates(star_table,self.stars,idx,self.apf_obs,dt):
            decker= "N"
            line  = ScriptobsLine.make_scriptobs_line(star_table[idx], \
                                        dt, decker=decker, I2="N", owner=res.owner, temp=True)
            if "decker=W" in line:
                decker = "W"
            bline = ScriptobsLine.make_scriptobs_line(star_table[bstars][bidx], dt, \
                                        decker=decker, I2="Y", owner=res.owner, focval=2)
            bfinline = ScriptobsLine.make_scriptobs_line(star_table[bstars][bfinidx], dt,\
                                            decker=decker, I2="Y", owner=res.owner, focval=0)
            res.set_template_lines([bfinline + " # temp=Y end",
                                    line + " # temp=Y",
                                    bline + " # temp=Y"], decker)
            apflog("Attempting template observation of %s" % (star_table['name'][idx]), echo=True)

