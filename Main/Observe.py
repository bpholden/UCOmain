from __future__ import print_function
import datetime
import math
import os
import os.path
import shutil
import re
import sys
import threading
import time

try:
    import ktl
    import APF as APFLib
    import APFTask
except:
    pass

import APFControl
from apflog import apflog
import UCOScheduler as ds
import ExposureCalculations
import SchedulerConsts

DMLIM = 1140

class Observe(threading.Thread):
    """ Observe(apf, opt, totTemps=4, task='master')
        The Observe class is a thread
        that runs the observing process.
    """
    def __init__(self, apf, opt, totTemps=4, task='master'):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.apf = apf
        self.task = task
        if opt.name:
            self.user = opt.name
        else:
            self.user = 'apf'
        if opt.owner:
            self.owner = opt.owner
        else:
            self.owner = 'public'

        if opt.windshield:
            self.windshield_mode = opt.windshield
        else:
            self.windshield_mode = 'auto'

        self.name = 'Observe'
        self.signal = True
        self.scriptobs = None

        self.BV = None
        self.VMAG = None
        self.blank = False
        self.decker = "W"

        self.obs_B_star = not opt.no_Bstar
        self.last_obs_success = True
        self.last_obs_finished = True
        self.star_failures = 0

        if opt.fixed:
            self.fixedList = opt.fixed
        else:
            self.fixedList = None
        if opt.sheet:
            self.sheetn = opt.sheet
        else:
            self.sheetn = None
        if opt.rank_table:
            self.rank_tablen = opt.rank_table
        else:
            self.rank_tablen = None
        if opt.start:
            try:
                self.starttime = float(opt.start)
            except ValueError as e:
                apflog("ValueError: %s" % (e), echo=True, level='error')
                self.starttime = None
        else:
            self.starttime = None
        if opt.raster:
            self.raster = opt.raster
        else:
            self.raster = False
        if opt.test:
            self.debug = opt.test
        else:
            self.debug = False
        self.doTemp = True
        self.nTemps = 0
        self.focval = 0
        self.totTemps = totTemps

        self.target = None
        self.fixedtarget = None
        self.bad_weather = False

        self.apftask = ktl.Service('apftask')
        self.lineresult = self.apftask['SCRIPTOBS_LINE_RESULT']
        self.lineresult.monitor()
        self.observed = self.apftask['SCRIPTOBS_OBSERVED']
        self.observed.monitor()
        self.phase = self.apftask['SCRIPTOBS_PHASE']
        self.phase.monitor()
        self.selected = None
        self.notify_focus_failure = True

        if opt.sheet is False:
            sheetlist = self.apftask['MASTER_SHEETLIST'].read().split(",")
            if len(sheetlist) > 0:
                self.sheetn = sheetlist
            else:
                self.sheetn = ["RECUR_A100",]

        self.canOpen = True

    def append_selected(self, curstr):
        try:
            self.selected = open("selected_targets", "a+")
        except:
            self.selected = None
        else:
            out_line = "%s %s\n" % (str(datetime.datetime.utcnow()), curstr)
            self.selected.write(out_line)
            self.selected.close()

        return

    def check_scriptobs_messages(self):
        message = self.apf.message.read()
        mtch = re.search("ERR/UCAM", message)
        if mtch:
            # uh oh
            apflog("scriptobs has failed post UCAM recovery", level="error", echo=True)
            # reboot warsaw
            rv = self.apf.ucam_restart()
            if rv:
                self.apf.message.write("")
                return True
            else:
                return False

        mtch = re.search("ERR/WIND", message)
        if mtch and self.apf.ok2open is True:
            # uh oh
            apflog("scriptobs has failed - checking servos", level="error", echo=True)
            rv = self.check_servos()
            if rv is False:
                return False
            self.apf.message.write("")
        return True

    def check_obs_success(self):
        """ Observe.check_obs_success()
            checks the value of SCRIPTOBS_LINE_RESULT to see if the last observation suceeded.
        """
        retval = False

        if self.lineresult.binary == 3:
            retval = True
        else:
            if "ERR/WIND"  in self.apf.robot["MASTER_MESSAGE"].read():
                apflog("Windshield error, check for faults", echo=True, level='error')
                return retval
            if self.lineresult.binary == 2:
                apflog("Observation failed at phase %s" % (self.phase.ascii), \
                       echo=True, level='warn')
        return retval

    def check_obs_finished(self):
        """ Observe.check_obs_finished()
            checks the value of SCRIPTOBS_LINE to see if we are on the last line of the block
            checks SCRIPTOBS_LINE_RESULT and SCRIPTOBS_OBSERVED to see if the last line is done
        """
        retval = False

        mtch = re.search("end\Z", self.apf.line.read())
        if self.apf.ldone.read(binary=True) == 0 or mtch:
            retval = True
        return retval


    def check_star(self, haveobserved):
        """ Observe.obsBstar(haveobserved)
            if observing has begun, and the last observation was a success,
            set Observe.obsBstar to false, writes master_obsbstar to
            the current value of obsBstar
            The variable OBSBSTAR still overrides
        """
        self.obs_B_star = ktl.read('apftask', 'MASTER_OBSBSTAR', binary=True)

        if haveobserved and self.last_obs_success:
            self.obs_B_star = False
            try:
                ktl.write('apftask', 'MASTER_OBSBSTAR', self.obs_B_star, binary=True)
            except Exception as e:
                apflog("Error: Cannot communicate with apftask: %s" % (e), level="error")
            self.star_failures = 0
        else:
            self.star_failures += 1
            if self.star_failures%3 == 0:
                log_str = "%d failures of observing a star in a row " % (self.star_failures)
                log_str += "- suggesting homing telescope or closing for the night"
                apflog(log_str, echo=True, level='timed_alert')

    def check_servos(self):
        """ Observe.check_servos()
            checks for servo faults and power cycles the telescope if necessary
        """
        _, running = self.apf.find_robot()
        if running:
            self.apf.kill_robot(now=True)

        chk_done = "$checkapf.MOVE_PERM == true"
        result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=600)
        if result:

            rv = self.apf.servo_failure()
            if rv:

                rv = self.apf.power_down_telescope()
                if rv:
                    apflog("Power cycled telescope", echo=True)
                else:
                    apflog("Failure power cycling telescope", echo=True, level="alert")

                return rv

            apflog("No current servo faults", echo=True)
            return True

        elif result is False and "DomeShutter" in self.apf.is_open()[1]:
            apflog("Error: After 10 min move permission did not return, and the dome is still open.", level='error', echo=True)
            self.apf.close(force=True)
            return False

    def check_files(self, outfn='googledex.dat'):
        """ Observe.check_files(outfn='googledex.dat')
            checks for the existence of a file, and if it exists, makes a backup
        """
        outdir = os.getcwd()
        fullpath = os.path.join(outdir, outfn)
        if os.path.isfile(fullpath):
            return

        # make it so
        backup = fullpath + ".1"
        try:
            shutil.copyfile(backup, fullpath)
        except Exception as e:
            err_str = "Cannot copy %s to %s: %s %s" % (backup, fullpath, type(e), e)
            apflog(err_str, echo=True, level='error')

        return

    def should_start_list(self):
        """ Observe.should_start_list()
            should we start a fixed observing list or not? true if start 
            time is None or if w/in + 1 hour - 0.5 hours of start time
        """
        if self.starttime is None:
            return True
        ct = time.time()
        if ct > self.starttime and ct - self.starttime < 3600:
            return True
        if ct < self.starttime and self.starttime - ct < 600:
            return True
        return False



    ####
    # run is the main event loop, for historical reasons
    # it has its own functions that are local in scope
    ####


    def run(self):
        """ Observe.run() - runs the observing
        """

        def calc_slowdown():

            if self.blank:
                return self.apf.robot["MASTER_SLOWDOWN"].read()

            if self.BV is None:
                apflog("Warning!: Ended up in get_target() with no B Magnitude value, slowdown can't be computed.", echo=True)
                self.BV = 0.6 # use a default average

            if self.VMAG is None:
                apflog("Warning!: Ended up in get_target() with no V Magnitude value, slowdown can't be computed.", echo=True)
                return 5

            if self.apf.avg_fwhm < 1.0:
                apflog("Warning!: AVG_FWHM = %4.2f. By Odin's beard that seems low."\
                        % self.apf.avg_fwhm, echo=True)
                return SchedulerConsts.SLOWDOWN_MAX

            slowdown = 1
            apflog("Calculating expected counts")
            apflog("self.VMAG [%4.2f] - self.BV [%4.2f] - self.apf.ael [%4.2f]"\
                    % (self.VMAG, self.BV, self.apf.ael))
            exp_cnts_sec = ExposureCalculations.getEXPMeter_Rate(self.VMAG, \
                                                                 self.BV, self.apf.ael, \
                                                                    self.apf.avg_fwhm, \
                                                                        self.decker)
            try:
                if self.apf.countrate <= 0:
                    try:
                        self.apf.countrate = self.apf.ccountrate
                    except:
                        self.apf.countrate = -1.0
                if self.apf.countrate*10 < self.apf.ccountrate:
                    self.apf.countrate = self.apf.ccountrate
                slowdown = exp_cnts_sec / self.apf.countrate
                if slowdown < 0:
                    slowdown = 1
                    apflog("Countrate non-sensical %g" % (self.apf.countrate), \
                           echo=True, level='warn')
                    self.apf.kcountrate.monitor(start=False)
                    self.apf.kcountrate.monitor(start=True)
                    self.apf.kcountrate.callback(self.apf.countrate_mon)
                    # yes this happened.
                if slowdown < SchedulerConsts.SLOWDOWN_MIN:
                    slowdown = SchedulerConsts.SLOWDOWN_MIN
                    apflog("slowdown too low, countrate= %g" % (self.apf.countrate), \
                           echo=True, level='debug')
                    # yes this happened.
                if slowdown > SchedulerConsts.SLOWDOWN_MAX:
                    slowdown = SchedulerConsts.SLOWDOWN_MAX
                    apflog("slowdown too high, countrate= %g" % (self.apf.countrate), \
                           echo=True, level='debug')
            except ZeroDivisionError:
                apflog("Current countrate was 0. Slowdown will be set to 1.", echo=True)
                slowdown = 1

            apflog("countrate = %.2f, ccountrate = %.2f" % (self.apf.countrate, \
                                                            self.apf.ccountrate))
            apflog("slowdown factor = %4.2f" % slowdown, echo=True)
            APFLib.write(self.apf.robot["MASTER_SLOWDOWN"], slowdown)
            return slowdown

        def pop_next():
            '''
            pop_next() - pops the next target from the target queue if any
            '''

            curstr = None

            if self.target is not None and 'SCRIPTOBS' in list(self.target.keys()):
                tlist = self.target["SCRIPTOBS"]
                if len(tlist) > 0:
                    apflog("get_target(): Going through remaining target queue.", echo=True)
                    curstr = tlist.pop()
                    return curstr

            if self.fixedtarget is not None and 'SCRIPTOBS' in list(self.fixedtarget.keys()):
                tlist = self.fixedtarget["SCRIPTOBS"]
                if len(tlist) > 0:
                    apflog("get_target(): Going through fixed starlist.", echo=True)
                    curstr = tlist.pop()
                else:
                    apflog("get_target(): Finished fixed starlist.", echo=True)
                    self.fixedtarget = None

            return curstr

        def empty_queue():
            '''
            empty_queue() - empties the target queue
            '''
            if self.target is not None and 'SCRIPTOBS' in list(self.target.keys()):
                while len(self.target["SCRIPTOBS"]) > 0:
                    self.target["SCRIPTOBS"].pop()

            return

        # This is called when an observation finishes, and selects the next target
        def get_target():

            if self.check_obs_finished():
                apflog("get_target(): Scriptobs phase is input, determining next target.", echo=True)
            else:
                apflog("get_target(): Not at end of block but out of targets.", echo=True)

            self.obs_B_star = ktl.read("apftask", "MASTER_OBSBSTAR", binary=True)
            apflog("get_target(): Setting obsBstar to %s" % (str(self.obs_B_star)), echo=True)

            if self.scriptobs is None:
                apflog("Called get_target, but there is not instance of scriptobs associated with %s. This is an error condition." % (self.name), level='error', echo=True)
                ripd, running = self.apf.find_robot()
                if running:
                    apflog("Attempting to kill the existing robot, %d" % (ripd), level='error', echo=True)
                    self.apf.kill_robot()
                return

            # Calculate the slowdown factor.
            slowdown = calc_slowdown()

            # Check for a valid seeing measurment. If there isn't one, use a default
            if self.apf.avg_fwhm == 0.:
                apflog("get_target(): Warning AVG_FWHM is 0. A default value of 15 will be used in its place.", echo=True)
                seeing = 15
            else:
                seeing = float(self.apf.avg_fwhm)
                apflog("get_target(): Current AVG_FWHM = %5.2f" % seeing)

            if self.apf.hatch_correct() == False:
                apflog("get_target(): Error setting hatch position.", level='Alert')
                return

            if self.apf.init_guide_cam() == False:
                apflog("get_target(): Error initializing guide camera.", echo=True, level='warn')
                if not self.apf.gcam_power.binary:
                    return
            self.apf.update_windshield(self.windshield_mode)
            self.focval = self.apf.set_autofoc_val()

            # setup a B star observation if needed
            # if not B star observation, look at current stack of
            # observations and see if anything is left
            if self.obs_B_star:
                self.apf.autofoc.write("robot_autofocus_enable")
            else:
                curstr = pop_next()
                # if there is a target in the queue, use it
                if curstr:
                    self.append_selected("%s avgfwhm=%05.2f slowdown=%04.2f" % \
                                         (curstr, seeing, slowdown))
                    self.scriptobs.stdin.write(curstr + '\n')
                    return

            self.check_files()

            self.target = ds.get_next(time.time(), seeing, slowdown, bstar=self.obs_B_star, \
                                         sheetns=self.sheetn, owner=self.owner,  \
                                         template=self.doTemp, focval=self.focval, \
                                         rank_sheetn=self.rank_tablen,\
                                         start_time=self.starttime)

            if self.target is None:
                log_str = "No acceptable target was found. "
                log_str += "Since there does not seem to be anything to observe"
                log_str += "%s will now shut down." % (self.name)
                apflog(log_str, echo=True)
                # Send scriptobs EOF to finish execution -
                # wouldn't want to leave a zombie scriptobs running
                self.scriptobs.stdin.close()
                self.apf.close()
                if self.fixedList is None:
                    APFLib.write(self.apf.ldone, 0)
                self.apf.countrate = -1.0
                # sleep for a half hour to see if the clouds blow by
                APFTask.waitfor(self.task, True, timeout=60*30)
                return

            apflog("Observing target: %s" % self.target['NAME'], echo=True)
            APFTask.set(self.task, suffix="MESSAGE", value="Observing target", wait=False)
            cur_line = self.target["SCRIPTOBS"].pop()
            cur_line = cur_line.strip()
            out_line = "%s avgfwhm=%05.2f slowdown=%04.2f" % (cur_line, seeing, slowdown )
            self.append_selected(out_line)

            apflog("Binning = %s" % self.apf.ucam['BINNING'].read(),echo=True)

            try:
                self.scriptobs.stdin.write(cur_line + '\n')
            except IOError as e:
                apflog("Cannot observe target %s: IOError: %s"\
                        % (self.target['NAME'], e), echo=True, level='error')
                return


            # Set the Vmag and B-V mag of the latest target
            self.VMAG = self.target["VMAG"]
            self.BV = self.target["BV"]
            self.decker = self.target["DECKER"]
            istemp = str(self.target['isTemp'])
            if self.target["mode"] == 'B' or self.target["mode"] == 'A':
                self.blank = True
            else:
                self.blank = False

            apflog("get_target(): V=%.2f  B-V=%.2f Pri=%.2f "\
                    % (self.VMAG, self.BV, self.target["PRI"]))
            apflog("get_target(): FWHM=%.2f  Slowdown=%.2f  Countrate=%.2f"\
                    % (self.apf.avg_fwhm, slowdown, self.apf.countrate))

            apflog("get_target(): Target= %s Temp=%s" % (self.target["NAME"], istemp))
            apflog("get_target(): Counts=%.2f  EXPTime=%.2f  Nexp=%d"\
                    % (self.target["COUNTS"], self.target["EXP_TIME"], self.target["NEXP"]))
            if self.target['isTemp']:
                self.nTemps += 1
                if self.nTemps >= self.totTemps:
                    self.doTemp = False

            return

        # opens the dome & telescope, if sunset is True calls open at sunset, else open at night
        def opening(sunel, sunset=False):
            if self.canOpen == False:
                apflog("We cannot open, so not trying", level='Error', echo=True)
                return False
            if self.apf.robot["SLEW_ALLOWED"].read(binary=True) == False:
                apflog("Opening: Slewing not allowed, so not opening", echo=True)
                self.canOpen = False
                return False
            when = "night"
            if sunset:
                when = "sunset"
            mstr = "Open at %s" % (when)
            APFTask.set(self.task, suffix="MESSAGE", value=mstr, wait=False)

            result = self.apf.ucam_status()
            if result is False:
                apflog("Failure in UCAM status and restart!", level='timed_alert', echo=True)
            else:
                apflog("UCAM OK", echo=True)

            apflog("Running open at %s as sunel = %4.2f" % (when, float(sunel)), echo=True)
            apfopen, _ = self.apf.is_open()
            if apfopen:
                self.apf.dm_reset()
            else:
                self.apf.dm_zero()

            result = self.apf.openat(sunset=sunset)
            apflog("opening completed with result %s" % (result), echo=True)
            if result is False:
                apflog("opening hasn't successfully opened. Current sunel = %4.2f" % (float(sunel)), level='warn', echo=True)
                if float(sunel) < SchedulerConsts.SUNEL_ENDLIM:
                    result = self.apf.openat(sunset=sunset)
                    if not result and self.apf.openOK['binary']:
                        apflog("Error: opening has failed twice, likely needs intervention.", level='Alert', echo=True)
                        self.apf.close()
                        self.canOpen = False

            self.apf.check_FCUs()
            self.apf.dm_reset()
            empty_queue()

            return result


        # closing
        def closing(force=False):
            if running:
                self.apf.kill_robot(now=True)

            self.append_selected("closing")

            APFTask.set(self.task, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())

            rv = self.apf.disable_inst()
            rv = self.apf.close(force=force)
            if rv:
                return
            rv = self.apf.servo_failure()
            if rv:
                ostr = "Servo failure detected, power cycling telescope"
                apflog(ostr, level="timed_alert", echo=True)
                rv = self.apf.power_down_telescope()
                if rv:
                    apflog("Power cycled telescope", echo=True)
                else:
                    apflog("Failure power cycling telescope", echo=True, level="alert")

            self.apf.check_FCUs()
            ds.zero_last_objs_attempted()
            self.star_failures = 0
            return


        def check_tel_state():
            slewing = '$eostele.AZSSTATE == Slewing  or  $eostele.ELSSTATE == Slewing'
            tracking = '$eostele.AZSSTATE == Tracking and $eostele.ELSSTATE == Tracking'
            istracking = ktl.Expression(tracking)
            isslewing = ktl.Expression(slewing)

            if istracking.evaluate() or isslewing.evaluate():
                rv = True
            else:
                rv = False
            return rv

        def start_telescope():
            '''This starts up the telescope if the Az drive is disabled or the E-Stop State is True
            If the telescope is just disabled, the start up procedure for a new version 
            of scriptobs should clear that state.
            '''
            rv = False

            eosdome = ktl.Service('eosdome')
            isenabled = eosdome['AZDRVENA'].read(binary=True)
            isstopped = eosdome['ESTOPST'].read(binary=True)
            fullstop = eosdome['SWESTOP'].read(binary=True)
            if fullstop:
                rv = False
                # cannot start the telescope
            else:
                # we can!
                if isstopped:
                    eosdome['ESTOPCMD'].write('ResetEStop')
                if isenabled is False:
                    eosdome['AZENABLE'].write('Enable')
                isenabled = eosdome['AZDRVENA'].read(binary=True)
                isstopped = eosdome['ESTOPST'].read(binary=True)
                if isenabled and isstopped is False:
                    rv = True
                else:
                    rv = False

            return rv

        def read_starlist_file():
            tot = 0
            if self.fixedList is None:
                return 0
            self.fixedtarget = dict()
            self.fixedtarget["SCRIPTOBS"] = []
            apflog("Reading star list fixedlist %s" % (self.fixedList), echo=True)
            with open(self.fixedList, 'r') as f:
                for line in f:
                    sline = line.strip()
                    if sline == '':
                        continue
                    elif sline[0] == '#':
                        continue
                    else:
                        tot += 1
                        self.fixedtarget["SCRIPTOBS"].append(sline)
            self.fixedtarget["SCRIPTOBS"].reverse()


            if tot == 0:
                apflog("Error: starlist %s is empty" % (self.fixedList), level="error")
                self.fixedList = None
                self.starttime = None
                self.target = None
            else:
                apflog("%d total starlist lines and %d lines done." % (tot, self.apf.ldone))

            return tot


        # starts an instance of scriptobs
        def start_scriptobs():
            # Update the last obs file and hitlist if needed

            APFTask.set(self.task, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())
            self.apf.validate_UCAM_outputs()
            self.apf.update_windshield(self.windshield_mode)
            self.apf.apftask_mon(self.apf.metxfersta)
            self.apf.apftask_mon(self.apf.apfteqsta)
            self.apf.status_clear()

            _, running = self.apf.find_robot()
            if running:
                apflog("Scriptobs is already running yet start_scriptobs was called", level="warn", echo=True)
                return
            rv = self.check_scriptobs_messages()
            if rv is False:
                return

            expr = "$checkapf.MOVE_PERM = True and $checkapf.INSTR_PERM = True"
            perms = APFTask.waitFor(self.task, True, expression=expr,timeout=1200)
            if perms is False:
                apflog("Cannot start an instance of scriptobs because do not have permission", echo=True, level='error')
                return

            if self.fixedList is not None and self.should_start_list():
                # We wish to observe a fixed target list, in it's original order
                if not os.path.exists(self.fixedList):
                    apflog("Error: starlist %s does not exist" % (self.fixedList), level="error")
                    self.fixedList = None
                    self.starttime = None
                    APFLib.write(self.apf.robot["MASTER_STARLIST"], "")
                    APFLib.write(self.apf.robot["MASTER_UTSTARTLIST"], "")

                # this reads in the list and appends it to self.target

                tot = read_starlist_file()

                if self.apf.ldone == tot:
                    APFLib.write(self.apf.robot["MASTER_STARLIST"], "")
                    APFLib.write(self.apf.robot["MASTER_UTSTARTLIST"], "")
                    self.fixedList = None
                    self.starttime = None
                    self.target = None
                    if not self.apf.test:
                        APFTask.set(self.task, suffix="STARLIST", value="")
                    apflog("Finished fixed list on line %d, will start dynamic scheduler" % int(self.apf.ldone), echo=True)
                else:
                    apflog("Found Fixed list %s" % self.fixedList, echo=True)
                    apflog("Starting fixed list on line %d" % int(self.apf.ldone), echo=True)
                    self.fixedList = None

            else:
                if self.BV is None:
                    apflog("No B-V value at the moment", echo=True)
                    #self.BV = 0.028
                if self.VMAG is None:
                    apflog("No VMag value at the moment", echo=True)
                    #self.VMAG = None
                # We wish to observe with the dynamic scheduler
            _, running = self.apf.find_robot()
            if running is False:
                apflog("Starting an instance of scriptobs for observing.", echo=True)
                self.scriptobs = self.apf.start_robot()
                # Don't let the watcher run over the robot starting up
                APFTask.waitFor(self.task, True, timeout=10)




            return




        ###############################

        # Actual Watching loop
        apflog("Beginning observing process....", echo=True)
        self.apf.dm_zero()
        haveobserved = False
        failstart = 0
        do_msg = 0

        while self.signal:
            # Check on everything
            if self.apf.sun_rising():
                rising = True
                sunel_lim = SchedulerConsts.SUNEL_ENDLIM
            else:
                rising = False
                sunel_lim = SchedulerConsts.SUNEL_STARTLIM

            _, running = self.apf.find_robot()
            cursunel = self.apf.sunel
            current_msg = APFTask.get("master", ["MESSAGE"])
            focusing = self.apf.focussta['binary'] < 3
            self.apf.status_clear()
            # Check and close for weather

            self.bad_weather = self.apf.dew_too_close or not self.apf.openOK \
                or not self.apf.gcam_power.binary

            if self.apf.is_open()[0] and self.bad_weather:
                closetime = datetime.datetime.now()
                APFTask.set(self.task, suffix="MESSAGE", \
                            value="Closing for weather or instrument issues", wait=False)
                apflog("No longer ok to open: %s." % (closetime), echo=True)
                apflog("OPREASON: " + self.apf.checkapf["OPREASON"].read(), echo=True)
                apflog("WEATHER: " + self.apf.checkapf['WEATHER'].read(), echo=True)
                apflog("CLOSE TO DEW POINT: %s" % (str(self.apf.dew_too_close)), echo=True)
                apflog("Guider camera power: %s" % ("ON" if self.apf.gcam_power.binary else "OFF"), \
                       echo=True)
                closing()

            self.apf.userkind.read(timeout=1)
            if self.apf.userkind.binary != 3:
                if do_msg == 0:
                    msg = "checkapf.USERKIND is no longer Robotic, instead %s" % (self.apf.userkind.ascii)
                    apflog(msg, echo=True, level='error')
                    APFTask.set(self.task, suffix="MESSAGE", value=msg, wait=False)
                    do_msg += 1
                APFTask.waitFor(self.task, True, timeout=60)
                continue
            elif self.apf.userkind.binary == 3:
                do_msg = 0


            # Check the slowdown factor to close for clouds
            if self.VMAG is not None and self.BV is not None and False:
                slow = calc_slowdown()
                APFTask.set(self.task, suffix="MESSAGE", value="FWHM = %.2f and slowdown %.2f" % (self.apf.avg_fwhm, slow), wait=False)
                if slow > 16:
                    # The slowdown is too high, we should close up and wait.
                    APFTask.set(self.task, suffix="MESSAGE", value="Closing for clouds", wait=False)
                    apflog("Slowdown factor of %.2f is too high. Waiting 30 min to check again." % slow, echo=True)
                    closing()
                    APFTask.waitfor(self.task, True, timeout=60*30)
                    self.VMAG = None
                    self.BV = None
                    self.apf.countrate = 0


            # If scriptobs is running and waiting for input, give it a target
            if running and (float(cursunel) < sunel_lim) and (self.apf.sop.read().strip() == "Input"):
                apflog("Entering target section", echo=True)
                if self.fixedList is None or not self.should_start_list():
                    self.last_obs_success = self.check_obs_success()
                    self.check_star(haveobserved)

                    APFTask.set(self.task, suffix="MESSAGE", value="Calling get_target", wait=False)
                    apflog("Scriptobs phase is input ( dynamic scheduler ), calling get_target.")
                    get_target()
                    APFTask.waitfor(self.task, True, timeout=15)

                    haveobserved = True
                elif self.starttime is not None and self.should_start_list():
                    apflog("Observing a fixed list called %s" % (self.fixedList), echo=True)
                    tot = read_starlist_file()
                    if tot == 0:
                        apflog("Error: starlist %s is empty" % (self.fixedList), level="error")
                        self.fixedList = None
                        self.starttime = None
                        self.target = None
                    else:
                        apflog("%d total starlist lines and %d lines done." % (tot, self.apf.ldone))
                        while len(self.target["SCRIPTOBS"]) > 0:
                            self.scriptobs.stdin.write(self.target["SCRIPTOBS"].pop() + '\n')
                        self.fixedList = None


            # If the sun is rising and we are finishing an observation
            # Send scriptobs EOF. This will shut it down after the observation
            if float(cursunel) >= sunel_lim and running:
                if self.scriptobs is None:
                    apflog("Robot claims to be running, but no scriptobs instance can be found. Instead calling kill_robot().", echo=True)
                    self.apf.kill_robot()
                else:
                    self.scriptobs.stdin.close()
                    self.apf.kill_robot()
                omsg = "Stopping scriptobs"
                if current_msg['MESSAGE'] != omsg:
                    APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)

            # If the sun is rising and scriptobs has stopped, run closeup
            if float(cursunel) > sunel_lim and not running and rising:
                outstr = "Closing due to sun elevation. Sunel = % 4.2f" % float(cursunel)
                apflog(outstr, echo=True)
                APFTask.set(self.task, suffix="MESSAGE", value=outstr, wait=False)
                if self.apf.is_open()[0]:
                    msg = "APF is open, closing due to sun elevation = %4.2f" % float(cursunel)
                    closing()
                else:
                    msg = "Telescope was already closed when sun got to %4.2f" % float(cursunel)

                if self.apf.is_open()[0]:
                    apflog("Error: Closeup did not succeed", level='error', echo=True)

                self.exitMessage = msg
                self.stop()

            # Open
            if self.apf.openOK and self.canOpen and not self.bad_weather:
                if not self.apf.is_ready_observing()[0] and float(cursunel) < SchedulerConsts.SUNEL_HOR:
                    if float(cursunel) > sunel_lim and not rising:
                        APFTask.set(self.task, suffix="MESSAGE", value="Open at sunset", wait=False)
                        success = opening(cursunel, sunset=True)
                        if success is False:
                            if self.apf.openOK:
                                apflog("Error: Cannot open the dome", level="timed_alert", echo=True)
                            else:
                                # lost permision during opening, happens more often than you think
                                apflog("Error: No longer have opening permission", level="error", echo=True)

                        else:
                            rv = self.apf.evening_star()
                            if not rv:
                                apflog("evening star targeting and telescope focus did not work", level='warn', echo=True)

                            chk_done = "$eostele.SUNEL < %f" % (SchedulerConsts.SUNEL_STARTLIM*math.pi/180.0)
                            while float(cursunel) > SchedulerConsts.SUNEL_STARTLIM and not rising:
                                outstr = "Sun is setting and sun at elevation of %.3f" % (float(cursunel))
                                apflog(outstr, level='info', echo=True)
                                result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=60)
                                self.apf.dm_reset()
                                if self.apf.openOK['binary'] is False:
                                    closetime = datetime.datetime.now()
                                    APFTask.set(self.task, suffix="MESSAGE", value="Closing for weather", wait=False)
                                    apflog("No longer ok to open.", echo=True)
                                    apflog("OPREASON: " + self.apf.checkapf["OPREASON"].read(), echo=True)
                                    apflog("WEATHER: " + self.apf.checkapf['WEATHER'].read(), echo=True)
                                    closing()
                                    break


                    elif not rising or (rising and float(cursunel) < (sunel_lim - 5)) and self.canOpen and not self.bad_weather:
                        success = opening(cursunel)
                        omsg = "Opening at %s" % (cursunel)
                        APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)

                    else:
                        success = True
                    if success is False:
                        self.apf.close()
                        if self.apf.openOK:
                            omsg = "Error: Cannot open the dome"
                            APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)
                            apflog(omsg, echo=True, level='error')
                        else:
                            apflog("Error: Lost permission during opening", echo=True)

                # If we can open, try to set stuff up so the vent doors can be controlled by apfteq
                if not rising and not self.apf.is_open()[0] and float(cursunel) > SchedulerConsts.SUNEL_HOR:
                    APFTask.set(self.task, suffix="MESSAGE", value="Powering up for APFTeq", wait=False)
                    if self.apf.clear_estop():
                        try:
                            APFLib.write(self.apf.dome['AZENABLE'], 'enable', timeout=10)
                        except:
                            apflog("Error: Cannot enable AZ drive", level="error")

                        self.apf.set_apfteq_mode('Evening')
                        vent_open = "$eosdome.VD4STATE = VENT_OPENED"
                        result = APFTask.waitfor(self.task, True, expression=vent_open, timeout=180)
                        if result:
                            try:
                                APFLib.write(self.apf.dome['AZENABLE'], 'disable', timeout=10)
                            except:
                                apflog("Error: Cannot disable AZ drive", level="warn", echo=True)

                        else:
                            if self.apf.openOK:
                                ostr="Error: Vent doors did not open, is apfteq and eosdome running correctly?"
                                apflog(ostr, level='info', echo=True)
                            else:
                                apflog("Error: Lost permission during attempt at opening",\
                                        level='info', echo=True)
                    else:
                        apflog("Error: Cannot clear emergency stop, sleeping for 600 seconds",\
                                level="error")
                        APFTask.waitFor(self.task, True, timeout=600)

            else:
                pass

            # Check for servo errors
            if not self.apf.slew_allowed.read(binary=True) and self.apf.is_ready_observing()[0]:
                apflog("Likely amplifier failure, may power cycle telescope",\
                        echo=True, level='error')
                rv = self.check_servos()

            # If we are open and scriptobs isn't running, start it up
            if self.apf.is_ready_observing()[0] and not running \
                and float(cursunel) <= sunel_lim and self.apf.openOK and not focusing:
                rv = APFTask.waitFor(self.task, False,\
                                      expression="$apftask.CALIBRATE_STATUS == 'Running'",\
                                          timeout=1)
                if rv:
                    try:
                        APFTask.abort("CALIBRATE")
                    except ktl.ktlError:
                        apflog("Warning: CALIBRATE still running after abort, this happens",\
                                echo=True, level='warn')
                    except Exception as e:
                        apflog("Error: Cannot abort CALIBRATE: %s" % (e), echo=True, level='error')
                    rv = APFTask.waitFor(self.task, False,\
                                          expression="$apftask.CALIBRATE_STATUS != 'Running'",\
                                              timeout=300)
                    if rv is False:
                        apflog("Error: CALIBRATE did not stop", echo=True, level='warn')
                    else:
                        apflog("CALIBRATE has stopped", echo=True)

                APFTask.set(self.task, suffix="MESSAGE", value="Starting scriptobs", wait=False)

                result = self.apf.enable_obs_inst()
                if result is False:
                    apflog("Cannot enable instrument", level='warn', echo=True)
                    result = self.apf.enable_obs_inst()
                    if not result:
                        apflog("Error: cannot enable instrument twice.", level='timed_alert', echo=True)
                        return result
                else:
                    apflog("Instrument OK", echo=True)

                rv = check_tel_state()
                if rv is False:
                    # this means that the telescope is not slewing and is not tracking
                    rv = start_telescope()
                    # if needed, will power up the Az drive and clear the estop state
                    if rv is False:
                        ostr = "Error: Telescope is not tracking or slewing, cannot start up telescope."
                        apflog(ostr, level='timed_alert', echo=True)
                        closing(force=True)

                start_scriptobs()
                expr = "$apftask.SCRIPTOBS_STATUS == 'Running'"
                if not APFTask.waitFor(self.task, True, expression=expr, timeout=10):
                    failstart += 1
                    if failstart % 11 == 0 and failstart > 0:
                        lvl = "timed_alert"
                    else:
                        lvl = "warn"
                    ostr = "Scriptobs is not running just after being started!"
                    apflog(ostr, level=lvl, echo=True)
                    APFTask.set(self.task, suffix="MESSAGE", value=ostr, wait=False)
                omsg = "Starting scriptobs"
                if current_msg['MESSAGE'] != omsg:
                    APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)


            # Keep an eye on the deadman timer if we are open
            if self.apf.is_open()[0] and self.apf.dmtime <= DMLIM:
                self.apf.dm_reset()

            if not self.apf.is_open()[0] and not rising:
                omsg = "Waiting for sunset"
                if current_msg['MESSAGE'] != omsg:
                    APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)
                APFTask.waitFor(self.task, True, timeout=5)

            if  self.apf.is_open()[0] and float(cursunel) > sunel_lim and not rising:
                omsg = "Waiting for sunset"
                if current_msg['MESSAGE'] != omsg:
                    APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)
                APFTask.waitFor(self.task, True, timeout=5)


    def stop(self):
        self.signal = False
        self.apf.kill_robot()

if __name__ == "__main__":

    class Test:
        def __init__(self):
            self.owner = 'public'
            self.name = 'apf'
            self.windshield = 'auto'
            self.fixed = None
            self.sheet = 'RECUR_A100'
            self.rank_table = '2024B_ranks'
            self.start = None
            self.test = True
            self.raster = False
            self.no_Bstar = False


    parent = 'example'
    APFTask.establish(parent, os.getpid())

    t_opt = Test()

    t_apf = APFControl.APF(task=parent, test=True)
    APFTask.waitFor(parent, True, timeout=2)
    print(str(t_apf))

    observe = Observe(t_apf, t_opt, task=parent)
    APFTask.waitFor(parent, True, timeout=2)
    observe.start()
    while observe.signal:
        try:
            dt = datetime.datetime.now()
            print(dt)
            APFTask.wait(parent, True, timeout=100)
        except KeyboardInterrupt:
            apflog("%s has been killed by user." % (observe.name), echo=True)
            sys.exit()
    print("Done")
