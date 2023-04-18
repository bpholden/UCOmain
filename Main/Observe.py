from __future__ import print_function
from datetime import datetime, timedelta
import os
import os.path
import shutil
import signal
from select import select
import re
import subprocess
import sys
import threading
import time


import numpy as np

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

        self.obsBstar = True
        self.lastObsSuccess = True
        self.lastObsFinished = True
        self.starFailures = 0

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

        self.apftask = ktl.Service('apftask')
        self.lineresult = self.apftask['SCRIPTOBS_LINE_RESULT']
        self.lineresult.monitor()
        self.observed = self.apftask['SCRIPTOBS_OBSERVED']
        self.observed.monitor()
        self.selected = None

        if opt.sheet is False:
            sheetlist = self.apftask['MASTER_SHEETLIST'].read().split(",")
            if len(sheetlist) > 0:
                self.sheetn = sheetlist
            else:
                self.sheetn = ["RECUR_A100",]
                
        self.canOpen = True
        self.badWeather = False

    def append_selected(self, curstr):
        try:
            self.selected = open("selected_targets", "a+")
        except:
            self.selected = None
        else:
            out_line = "%s %s\n" % (str(datetime.utcnow()), curstr)
            self.selected.write(out_line)
            self.selected.close()

        return

    def checkScriptobsMessages(self):
        message = self.apf.message.read()
        mtch = re.search("ERR/UCAM", message)
        if mtch:
            # uh oh
            apflog("scriptobs has failed post UCAM recovery", level="error", echo=True)
            # reboot warsaw
            rv = self.apf.ucamRestart()
            if rv:
                self.apf.message.write("")
                return True
            else:
                return False

        mtch = re.search("ERR/WIND", message)
        if mtch:
            # uh oh
            apflog("scriptobs has failed - checking servos", level="error", echo=True)
            rv = self.checkServos()
            if rv is False:
                return False
        return True

    def checkObsSuccess(self):
        """ Observe.checkObsSuccess()
            checks the value of SCRIPTOBS_LINE_RESULT to see if the last observation suceeded.
        """
        retval = False

        if self.observed.read(binary=True) == True:
            retval = True
        return retval

    def checkObsFinished(self):
        """ Observe.checkObsFinished()
            checks the value of SCRIPTOBS_LINE to see if we are on the last line of the block
            checks SCRIPTOBS_LINE_RESULT and SCRIPTOBS_OBSERVED to see if the last line is done
        """
        retval = False

        mtch = re.search("end\Z", self.apf.line.read())
        if self.apf.ldone.read(binary=True) == 0 or mtch:
            retval = True
        return retval


    def checkStar(self, haveobserved):
        """ Observe.obsBstar(haveobserved)
            if observing has begun, and the last observation was a success, set Observe.obsBstar to false, writes master_obsbstar to
            the current value of obsBstar
            The variable OBSBSTAR still overrides
        """
        self.obsBstar = ktl.read('apftask', 'MASTER_OBSBSTAR', binary=True)

        if haveobserved and self.lastObsSuccess:
            self.obsBstar = False
            try:
                ktl.write('apftask', 'MASTER_OBSBSTAR', self.obsBstar, binary=True)
            except Exception as e:
                apflog("Error: Cannot communicate with apftask: %s" % (e), level="error")
            self.starFailures = 0
        else:
            self.starFailures += 1
            if self.starFailures%3 == 0:
                apflog("%d failures of observing a star in a row - suggesting homing telescope or closing for the night" % (self.starFailures), echo=True, level='Alert')

    def checkServos(self):

        _, running = self.apf.findRobot()
        if running:
            self.apf.killRobot(now=True)

        chk_done = "$checkapf.MOVE_PERM == true"
        result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=600)
        if result:

            rv = self.apf.servoFailure()
            if rv:

                rv = self.apf.powerDownTelescope()
                if rv:
                    apflog("Power cycled telescope", echo=True)
                else:
                    apflog("Failure power cycling telescope", echo=True, level="alert")

                return rv
            
            apflog("No current servo faults", echo=True)
            
        elif result is False and "DomeShutter" in self.apf.isOpen()[1]:
            apflog("Error: After 10 min move permission did not return, and the dome is still open.", level='error', echo=True)
            self.apf.close(force=True)
            return False

    def checkFiles(self, outfn='googledex.dat'):
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
        
    def shouldStartList(self):
        """ Observe.shouldStartList()
            should we start a fixed observing list or not? true if start time is None or if w/in + 1 hour - 0.5 hours of start time
        """
        if self.starttime == None:
            return True
        ct = time.time()
        if ct > self.starttime and ct - self.starttime < 3600:
            return True
        if ct < self.starttime and self.starttime - ct < 1800:
            return True
        return False



    ####
    # run is the main event loop, for historical reasons it has its own functions that are local in scope
    ####


    def run(self):
        """ Observe.run() - runs the observing
        """

        def calcSlowdown():

            if self.blank:
                return self.apf.robot["MASTER_SLOWDOWN"].read()

            if self.BV is None:
                apflog("Warning!: Ended up in getTarget() with no B Magnitude value, slowdown can't be computed.", echo=True)
                self.BV = 0.6 # use a default average

            if self.VMAG is None:
                apflog("Warning!: Ended up in getTarget() with no V Magnitude value, slowdown can't be computed.", echo=True)
                return 5

            if self.apf.avg_fwhm < 1.0:
                apflog("Warning!: AVG_FWHM = %4.2f. By Odin's beard that seems low." % self.apf.avg_fwhm, echo=True)
                return SchedulerConsts.SLOWDOWN_MAX
            
            slowdown = 1
            apflog("Calculating expected counts")
            apflog("self.VMAG [%4.2f] - self.BV [%4.2f] - self.apf.ael [%4.2f]" % (self.VMAG, self.BV, self.apf.ael))
            exp_cnts_sec = ExposureCalculations.getEXPMeter_Rate(self.VMAG, self.BV, self.apf.ael, self.apf.avg_fwhm, self.decker)
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
                    apflog("Countrate non-sensical %g" % (self.apf.countrate), echo=True, level='warn')
                    self.apf.counts.monitor(start=False)
                    self.apf.counts.monitor(start=True)
                    self.apf.counts.callback(self.apf.countmon)
                    # yes this happened.
                if slowdown < SchedulerConsts.SLOWDOWN_MIN:
                    slowdown = SchedulerConsts.SLOWDOWN_MIN
                    apflog("slowdown too low, countrate= %g" % (self.apf.countrate), echo=True, level='debug')
                    # yes this happened.
                if slowdown > SchedulerConsts.SLOWDOWN_MAX:
                    slowdown = SchedulerConsts.SLOWDOWN_MAX
                    apflog("slowdown too high, countrate= %g" % (self.apf.countrate), echo=True, level='debug')
            except ZeroDivisionError:
                apflog("Current countrate was 0. Slowdown will be set to 1.", echo=True)
                slowdown = 1

            apflog("countrate = %.2f, ccountrate = %.2f" % (self.apf.countrate, self.apf.ccountrate))
            apflog("slowdown factor = %4.2f" % slowdown, echo=True)
            APFLib.write(self.apf.robot["MASTER_SLOWDOWN"], slowdown)
            return slowdown

        def popNext():

            curstr = None

            if self.target is not None and 'SCRIPTOBS' in list(self.target.keys()):
                tlist = self.target["SCRIPTOBS"]
                if len(tlist) > 0:
                    apflog("getTarget(): Going through remaining target queue.", echo=True)
                    curstr = tlist.pop() 
                    return curstr

            if self.fixedtarget is not None and 'SCRIPTOBS' in list(self.fixedtarget.keys()):
                tlist = self.fixedtarget["SCRIPTOBS"]
                if len(tlist) > 0:
                    apflog("getTarget(): Going through fixed starlist.", echo=True)
                    curstr = tlist.pop() 
                else:
                    apflog("getTarget(): Finished fixed starlist.", echo=True)
                    self.fixedtarget = None

            return curstr

        # This is called when an observation finishes, and selects the next target
        def getTarget():
            APFLib.write(self.apf.ucam["RECORD"], "Yes") # safe / sorry

            if self.apf.nerase != 2:
                self.apf.nerase.write(2, binary=True)

            if self.checkObsFinished():
                apflog("getTarget(): Scriptobs phase is input, determining next target.", echo=True)
            else:
                apflog("getTarget(): Not at end of block but out of targets.", echo=True)


            self.obsBstar = ktl.read("apftask", "MASTER_OBSBSTAR", binary=True)
            apflog("getTarget(): Setting obsBstar to %s" % (str(self.obsBstar)), echo=True)

            if self.scriptobs is None:
                apflog("Called getTarget, but there is not instance of scriptobs associated with %s. This is an error condition." % (self.name), level='error', echo=True)
                ripd, running = self.apf.findRobot()
                if running:
                    apflog("Attempting to kill the existing robot, %d" % (ripd), level='error', echo=True)
                    self.apf.killRobot()
                return

            # Calculate the slowdown factor.
            slowdown = calcSlowdown()

            # Check for a valid seeing measurment. If there isn't one, use a default
            if self.apf.avg_fwhm == 0.:
                apflog("getTarget(): Warning AVG_FWHM is 0. A default value of 15 will be used in its place.", echo=True)
                seeing = 15
            else:
                seeing = float(self.apf.avg_fwhm)
                apflog("getTarget(): Current AVG_FWHM = %5.2f" % seeing)

            if self.apf.hatchCorrect() == False:
                apflog("getTarget(): Error setting hatch position.", level='Alert')
                return
            
            self.apf.initGuideCam()
            self.apf.updateWindshield(self.windshield_mode)


            # setup a B star observation if needed
            # if not B star observation, look at current stack of observations and see if anything is left
            if self.obsBstar:
                self.apf.autofoc.write("robot_autofocus_enable")
            else:
                curstr = popNext()
                if curstr:
                    self.append_selected("%s avgfwhm=%05.2f slowdown=%04.2f" % (curstr, seeing, slowdown))
                    self.scriptobs.stdin.write(curstr + '\n')
                    return

            self.focval = self.apf.setAutofocVal()

            self.checkFiles()
            delta_t = time.time() - self.apf.lastopen.binary
            self.target = ds.getNext(time.time(), seeing, slowdown, bstar=self.obsBstar, \
                                         sheetns=self.sheetn, owner=self.owner, template=self.doTemp, \
                                         focval=self.focval, rank_sheetn=self.rank_tablen, delta_t=delta_t)

            if self.target is None:
                apflog("No acceptable target was found. Since there does not seem to be anything to observe, %s will now shut down." % (self.name), echo=True)
                # Send scriptobs EOF to finish execution - wouldn't want to leave a zombie scriptobs running
                self.scriptobs.stdin.close()
                self.apf.close()
                if self.fixedList is None:
                    APFLib.write(self.apf.ldone, 0)
                self.apf.countrate = -1.0
                # sleep for a half hour to see if the clouds blow by
                APFTask.waitfor(self.task, True, timeout=60*30)
                return

            apflog("Observing target: %s" % self.target['NAME'], echo=True)
            APFTask.set(self.task, suffix="MESSAGE", value="Observing target: %s"  % self.target['NAME'], wait=False)
            cur_line = self.target["SCRIPTOBS"].pop()
            cur_line = cur_line.strip()
            out_line = "%s avgfwhm=%05.2f slowdown=%04.2f" % (cur_line, seeing, slowdown )
            self.append_selected(out_line)

            self.apf.ucam['BINNING'].write(self.target['BINNING'], timeout=0.1)
            apflog("Binning = %s" % self.apf.ucam['BINNING'].read(),echo=True)
            
            try:
                self.scriptobs.stdin.write(cur_line + '\n')
            except IOError as e:
                apflog("Cannot observe target %s: IOError: %s" % (self.target['NAME'], e), echo=True, level='error')
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

            apflog("getTarget(): V=%.2f  B-V=%.2f Pri=%.2f " % (self.VMAG, self.BV, self.target["PRI"]))
            apflog("getTarget(): FWHM=%.2f  Slowdown=%.2f  Countrate=%.2f" % (self.apf.avg_fwhm, slowdown, self.apf.countrate))

            apflog("getTarget(): Target= %s Temp=%s" % (self.target["NAME"], istemp))
            apflog("getTarget(): Counts=%.2f  EXPTime=%.2f  Nexp=%d" % (self.target["COUNTS"], self.target["EXP_TIME"], self.target["NEXP"]))
            if self.target['isTemp']:
                self.nTemps += 1
                if self.nTemps >= self.totTemps:
                    self.doTemp = False

            return

        # opens the dome & telescope, if sunset is True calls open at sunset, else open at night
        def opening(sunel, sunset=False):
            if self.canOpen is False:
                apflog("We cannot open, so not trying", level='Error', echo=True)
                return False
            when = "night"
            if sunset:
                when = "sunset"
            mstr = "Open at %s" % (when)
            APFTask.set(self.task, suffix="MESSAGE", value=mstr, wait=False)

            result = self.apf.ucamStatus()
            if result is False:
                apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
            else:
                apflog("UCAM OK", echo=True)

            apflog("Running open at %s as sunel = %4.2f" % (when, float(sunel)), echo=True)
            apfopen, what = self.apf.isOpen()
            if apfopen:
                self.apf.DMReset()
            else:
                self.apf.DMZero()

            result = self.apf.openat(sunset=sunset)
            apflog("opening completed with result %s" % (result), echo=True)
            if result == False:
                apflog("opening hasn't successfully opened. Current sunel = %4.2f" % (float(sunel)), level='warn', echo=True)
                if float(sunel) < SchedulerConsts.SUNEL_ENDLIM:
                    result = self.apf.openat(sunset=sunset)
                    if not result and self.apf.openOK['binary']:
                        apflog("Error: opening has failed twice, likely needs intervention.", level='Alert', echo=True)
                        self.apf.close()
                        self.canOpen = False

            self.apf.checkFCs()

            if datetime.now().strftime("%p") == 'PM':
                setting = True
            else:
                setting = False
            self.apf.DMReset()

            return result


        # closing
        def closing(force=False):
            if running:
                self.apf.killRobot(now=True)

            self.append_selected("closing")

            APFTask.set(self.task, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())

            rv = self.apf.disableInst()
            rv = self.apf.close(force=force)
            if rv:
                return
            rv = self.apf.servoFailure()
            if rv:
                apflog("Servo Failure, cannot close and power off telescope ", level="alert", echo=True)
                rv = self.apf.powerDownTelescope()
                if rv:
                    apflog("Power cycled telescope", echo=True)
                else:
                    apflog("Failure power cycling telescope", echo=True, level="alert")

            self.apf.checkFCs()

            return


        def checkTelState():
            slewing = '$eostele.AZSSTATE == Slewing  or  $eostele.ELSSTATE == Slewing'
            tracking = '$eostele.AZSSTATE == Tracking and $eostele.ELSSTATE == Tracking'
            istracking = ktl.Expression(tracking)
            isslewing = ktl.Expression(slewing)

            if istracking.evaluate() or isslewing.evaluate():
                rv = True
            else:
                rv = False
            return rv

        def startTelescope():
            '''This starts up the telescope if the Az drive is disabled or the E-Stop State is True
            If the telescope is just disabled, the start up procedure for a new version of scriptobs should clear that state.
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

        def readStarlistFile():
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
        def startScriptobs():
            # Update the last obs file and hitlist if needed

            APFTask.set(self.task, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())

            self.apf.updateWindshield(self.windshield_mode)
            ripd, running = self.apf.findRobot()
            if running:
                apflog("Scriptobs is already running yet startScriptobs was called", level="warn", echo=True)
                return
            rv = self.checkScriptobsMessages()
            if rv is False:
                return

            expr = "$checkapf.MOVE_PERM = True and $checkapf.INSTR_PERM = True"
            perms = APFTask.waitFor(self.task, True, expression=expr,timeout=1200)
            if perms is False:
                apflog("Cannot start an instance of scriptobs because do not have permission", echo=True, level='error')
                return

            apflog("Starting an instance of scriptobs", echo=True)
            if self.fixedList is not None and self.shouldStartList():
                # We wish to observe a fixed target list, in it's original order
                if not os.path.exists(self.fixedList):
                    apflog("Error: starlist %s does not exist" % (self.fixedList), level="error")
                    self.fixedList = None
                    self.starttime = None
                    APFLib.write(self.apf.robot["MASTER_STARLIST"], "")
                    APFLib.write(self.apf.robot["MASTER_UTSTARTLIST"], "")

                # this reads in the list and appends it to self.target

                tot = readStarlistFile()

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
            _, running = self.apf.findRobot()
            if running is False:
                apflog("Starting an instance of scriptobs for dynamic observing.", echo=True)
                self.scriptobs = self.apf.startRobot()
                # Don't let the watcher run over the robot starting up
                APFTask.waitFor(self.task, True, timeout=10)




            return




        ###############################

        # Actual Watching loop
        apflog("Beginning observing process....", echo=True)
        self.apf.DMZero()
        haveobserved = False
        failstart = 0
        do_msg = 0

        self.apf.validateUCAMoutputs()
        
        while self.signal:
            # Check on everything
            if self.apf.sunRising():
                rising = True
                sunel_lim = SchedulerConsts.SUNEL_ENDLIM
            else:
                rising = False
                sunel_lim = SchedulerConsts.SUNEL_STARTLIM

            ripd, running = self.apf.findRobot()
            cursunel = self.apf.sunel
            current_msg = APFTask.get("master", ["MESSAGE"])
            focusing = (self.apf.focussta['binary'] < 3)
            calibrating = (self.apf.calsta['binary'] < 3)

            # Check and close for weather
            self.badweather = self.apf.dewTooClose or not self.apf.openOK

            if self.apf.isOpen()[0] and self.badweather:
                closetime = datetime.now()
                APFTask.set(self.task, suffix="MESSAGE", value="Closing for weather", wait=False)
                apflog("No longer ok to open.", echo=True)
                apflog("OPREASON: " + self.apf.checkapf["OPREASON"].read(), echo=True)
                apflog("WEATHER: " + self.apf.checkapf['WEATHER'].read(), echo=True)
                apflog("CLOSE TO DEW POINT: %s" % (str(self.apf.dewTooClose)), echo=True)
                closing()

            self.apf.userkind.read(timeout=1)
            if self.apf.userkind.binary != 3:
                if do_msg == 0:
                    msg = "checkapf.USERKIND is no longer Robotic, instead %s" % (self.apf.userkind.ascii)
                    apflog(msg, echo=True, level='error')
                    APFTASK.set(self.task, suffix="MESSAGE", value=msg, wait=False)
                    do_msg += 1
                APFTask.waitFor(self.task, True, timeout=60)
                continue
            elif self.apf.userkind.binary == 3:
                do_msg = 0
                
                    
            # Check the slowdown factor to close for clouds
            if self.VMAG is not None and self.BV is not None and False:
                slow = calcSlowdown()
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
                if self.fixedList is None or not self.shouldStartList():
                    self.lastObsSuccess = self.checkObsSuccess()
                    self.checkStar(haveobserved)

                    APFTask.set(self.task, suffix="MESSAGE", value="Calling getTarget", wait=False)
                    apflog("Scriptobs phase is input ( dynamic scheduler ), calling getTarget.")
                    getTarget()
                    APFTask.waitfor(self.task, True, timeout=15)

                    haveobserved = True
                elif self.starttime is not None and self.shouldStartList():
                    apflog("Observing a fixed list called %s" % (self.fixedList), echo=True)
                    tot = readStarlistFile()
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
                APFTask.set(self.task, suffix="MESSAGE", value="Last call", wait=False)
                if self.scriptobs is None:
                    apflog("Robot claims to be running, but no self.scriptobs instance can be found. Instead calling killRobot().", echo=True)
                    self.apf.killRobot()
                else:
                    self.scriptobs.stdin.close()
                    self.apf.killRobot()
                omsg = "Stopping scriptobs"
                if current_msg['MESSAGE'] != omsg:
                    APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)

            # If the sun is rising and scriptobs has stopped, run closeup
            if float(cursunel) > sunel_lim and not running and rising:
                apflog("Closing due to sun elevation. Sunel = % 4.2f" % float(cursunel), echo=True)
                APFTask.set(self.task, suffix="MESSAGE", value="Closing, sun is rising", wait=False)
                if self.apf.isOpen()[0]:
                    msg = "APF is open, closing due to sun elevation = %4.2f" % float(cursunel)
                    closing()
                else:
                    msg = "Telescope was already closed when sun got to %4.2f" % float(cursunel)

                if self.apf.isOpen()[0]:
                    apflog("Error: Closeup did not succeed", level='error', echo=True)

                self.exitMessage = msg
                self.stop()

            # Open
            if self.apf.openOK and self.canOpen and not self.badweather:
                APFTask.phase(self.task, "Observing")
                if not self.apf.isReadyForObserving()[0] and float(cursunel) < SchedulerConsts.SUNEL_HOR:
                    if float(cursunel) > sunel_lim and not rising:
                        APFTask.phase(self.task, "Watching")
                        APFTask.set(self.task, suffix="MESSAGE", value="Open at sunset", wait=False)
                        success = opening(cursunel, sunset=True)
                        if success is False:
                            if self.apf.openOK:
                                apflog("Error: Cannot open the dome", level="alert", echo=True)
                            else:
                                # lost permision during opening, happens more often than you think
                                apflog("Error: No longer have opening permission", level="error", echo=True)

                        else:
                            rv = self.apf.eveningStar()
                            if not rv:
                                apflog("evening star targeting and telescope focus did not work", level='warn', echo=True)

                            chk_done = "$eostele.SUNEL < %f" % (SchedulerConsts.SUNEL_STARTLIM*np.pi/180.0)
                            while float(cursunel) > SchedulerConsts.SUNEL_STARTLIM and not rising:
                                outstr = "Sun is setting and sun at elevation of %.3f" % (float(cursunel))
                                apflog(outstr, level='info', echo=True)
                                result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=60)
                                self.apf.DMReset()
                                if self.apf.openOK['binary'] is False:
                                    closetime = datetime.now()
                                    APFTask.set(self.task, suffix="MESSAGE", value="Closing for weather", wait=False)
                                    apflog("No longer ok to open.", echo=True)
                                    apflog("OPREASON: " + self.apf.checkapf["OPREASON"].read(), echo=True)
                                    apflog("WEATHER: " + self.apf.checkapf['WEATHER'].read(), echo=True)
                                    closing()
                                    break


                    elif not rising or (rising and float(cursunel) < (sunel_lim - 5)) and self.canOpen and not self.badweather:
                        success = opening(cursunel)
                        omsg = "Opening at %s" % (cursunel)
                        APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)

                    else:
                        success = True
                    if success == False:
                        self.apf.close()
                        if self.apf.openOK:
                            apflog("Error: Cannot open the dome", echo=True, level='error')
                        else:
                            apflog("Error: Lost permission during opening", echo=True)

                # If we can open, try to set stuff up so the vent doors can be controlled by apfteq
                if not rising and not self.apf.isOpen()[0] and float(cursunel) > SchedulerConsts.SUNEL_HOR:
                    APFTask.set(self.task, suffix="MESSAGE", value="Powering up for APFTeq", wait=False)
                    APFTask.phase(self.task, "Watching")
                    
                    if self.apf.clearestop():
                        try:
                            APFLib.write(self.apf.dome['AZENABLE'], 'enable', timeout=10)
                        except:
                            apflog("Error: Cannot enable AZ drive", level="error")

                        self.apf.setTeqMode('Evening')
                        vent_open = "$eosdome.VD4STATE = VENT_OPENED"
                        result = APFTask.waitfor(self.task, True, expression=vent_open, timeout=180)
                        if result:
                            try:
                                Apflib.write(self.apf.dome['AZENABLE'], 'disable', timeout=10)
                            except:
                                apflog("Error: Cannot disable AZ drive", level="warn", echo=True)

                        else:
                            if self.apf.openOK:
                                apflog("Error: Vent doors did not open, is apfteq and eosdome running correctly?", level='info', echo=True)
                            else:
                                apflog("Error: Lost permission during attempt at opening", level='info', echo=True)
                    else:
                        apflog("Error: Cannot clear emergency stop, sleeping for 600 seconds", level="error")
                        APFTask.waitFor(self.task, True, timeout=600)

            else:
                pass

            # Check for servo errors
            if not self.apf.slew_allowed.read(binary=True) and self.apf.isReadyForObserving()[0]:
                apflog("Likely amplifier failure, may power cycle telescope", echo=True, level='alert')
                rv = self.checkServos()

            # If we are open and scriptobs isn't running, start it up
            if self.apf.isReadyForObserving()[0] and not running and float(cursunel) <= sunel_lim and self.apf.openOK and not focusing:
                calstat = APFTask.get('CALIBRATE', ['STATUS'])
                if calstat['STATUS'] in ['Running', 'Pausing', 'Paused']:
                    APFTask.abort("CALIBRATE")

                APFTask.set(self.task, suffix="MESSAGE", value="Starting scriptobs", wait=False)
                
                result = self.apf.enableObsInst()
                if result == False:
                    apflog("Cannot enable instrument", level='warn', echo=True)
                    result = self.apf.enableObsInst()
                    if not result:
                        apflog("Error: cannot enable instrument twice.", level='alert', echo=True)
                        return result
                else:
                    apflog("Instrument OK", echo=True)
                    
                rv = checkTelState()
                if rv is False:
                    # this means that the telescope is not slewing and is not tracking
                    rv = startTelescope()
                    # if needed, will power up the Az drive and clear the estop state
                    if rv == False:
                        apflog("Telescope stopped and cannot be restarted", level='Alert', echo=True)
                        closing(force=True)

                startScriptobs()
                if not APFTask.waitFor(self.task, True, expression="$apftask.SCRIPTOBS_STATUS == 'Running'", timeout=10):
                    failstart += 1
                    if failstart % 11 == 0 and failstart > 0:
                        lvl = "Alert"
                    else:
                        lvl = "warn"
                    apflog("scriptobs is not running just after being started!", level=lvl, echo=True)
                    APFTask.set(self.task, suffix="MESSAGE", value="scriptobs is not running just after being started!", wait=False)
                omsg = "Starting scriptobs"
                if current_msg['MESSAGE'] != omsg:
                    APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)


            # Keep an eye on the deadman timer if we are open
            if self.apf.isOpen()[0] and self.apf.dmtime <= DMLIM:
                APFTask.set(self.task, suffix="MESSAGE", value="Reseting DM timer", wait=False)
                self.apf.DMReset()
#                apflog("The APF is open, the DM timer is clicking down, and scriptobs is %s." % ( str(running)), level="debug")

            if not self.apf.isOpen()[0] and not rising:
                omsg = "Waiting for sunset"
                if current_msg['MESSAGE'] != omsg:
                    APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)
                APFTask.waitFor(self.task, True, timeout=5)
                
            if  self.apf.isOpen()[0] and float(cursunel) > sunel_lim and not rising:
                omsg = "Waiting for sunset"
                if current_msg['MESSAGE'] != omsg:
                    APFTask.set(self.task, suffix="MESSAGE", value=omsg, wait=False)
                APFTask.waitFor(self.task, True, timeout=5)


    def stop(self):
        self.signal = False
        threading.Thread._Thread__stop(self)


if __name__ == "__main__":

    class Test:
        pass


    parent = 'example'
    APFTask.establish(parent, os.getpid())

    opt = Test()
    opt.owner = 'public'
    opt.name = 'apf'
    opt.windshield = 'auto'
    opt.fixed = None
    opt.sheet = 'Bstars'
    opt.rank_table = '2022A_ranks'
    opt.start = None
    opt.test = True
    opt.raster = False

    apf = APFControl.APF(task=parent, test=True)
    APFTask.waitFor(parent, True, timeout=2)
    print(str(apf))

    observe = Observe(apf, opt, task=parent)
    APFTask.waitFor(parent, True, timeout=2)
    observe.start()
    while observe.signal:
        try:
            dt = datetime.now()
            print(dt)
            APFTask.wait(parent, True, timeout=100)
        except KeyboardInterrupt:
            apflog("%s has been killed by user." % (observe.name), echo=True)
            sys.exit()
        except:
            apflog("%s killed by unknown." % (observe.name), echo=True)
            sys.exit()
