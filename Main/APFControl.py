from __future__ import print_function
# Class definition for an APF object which tracks the state of the telescope.

import subprocess
import time
import os
import os.path
import math
import datetime

import numpy as np

try:
    from apflog import apflog
    import ktl
    import APF as APFLib
    import APFTask
    import Exposure
except:
    from fake_apflog import *


slowlim = 100
wxtimeout = datetime.timedelta(seconds=1800)
DEWARMAX = 8600
DEWARMIN = 8300


if "LROOT" in os.environ:
    LROOT = os.environ["LROOT"]
else:
    LROOT = '/usr/local/lick'
SCRIPTDIR = os.path.join(LROOT,'bin/robot/')


def apftask_do(cmd, debug=True, cwd='./'):
    newcmd = "apftask do %s" % (cmd)
    rv, retcode = cmd_exec(newcmd, debug=debug, cwd=cwd)
    return rv, retcode


def cmd_exec(cmd, debug=False, cwd='./'):
    apflog("Executing Command: %s" % repr(cmd), echo=True)

    args = cmd.split()

    try:
        p = subprocess.Popen(args, stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=cwd)
    except OSError as e:
        apflog("command %s does not exist: %s" % (cmd,e),echo=True)
        return False, -1
    except Exception as e:
        apflog("command %s failed: %s" % (cmd,e),echo=True)
        return False, -1

    while p.poll() is None:
        if debug:
            l = p.stdout.readline().rstrip('\n')
            apflog(l, echo=True)

    out, err = p.communicate()
    if debug: apflog(out, echo=True)
    if len(err): apflog(err, echo=True)
    ret_code = p.returncode
    if ret_code == 0:
        return True, ret_code
    else:
        return False, ret_code





class APF:
    """ Class which creates a monitored state object to track the condition of the APF telescope. """

    # Initial seeing conditions


    # KTL Services and Keywords
    robot        = ktl.Service('apftask')
    vmag         = robot['SCRIPTOBS_VMAG']
    ldone        = robot['SCRIPTOBS_LINES_DONE']
    line         = robot['SCRIPTOBS_LINE']
    sop          = robot['SCRIPTOBS_PHASE']
    message      = robot['SCRIPTOBS_MESSAGE']
    autofoc      = robot["SCRIPTOBS_AUTOFOC"]
    slew_allowed = robot['SLEW_ALLOWED']
    observed     = robot['SCRIPTOBS_OBSERVED']
    line_result  = robot['SCRIPTOBS_LINE_RESULT']

    apfteqsta    = robot['APFTEQ_STATUS']
    metxfersta   = robot['METSXFER_STATUS']
    calsta       = robot['CALIBRATE_STATUS']
    focussta     = robot['FOCUSINSTR_STATUS']
    slewsta      = robot['SLEW_STATUS']
    opensta      = robot['OPENUP_STATUS']
    closesta     = robot['CLOSEUP_STATUS']
    shuttersta   = robot['SHUTTERS_STATUS']
    focustelsta  = robot['FOCUSTEL_STATUS']

    ucamcmd      = robot['UCAMLAUNCHER_UCAM_COMMAND']
    lastopen     = robot['OPENUP_LAST_SUCCESS']

    ucam       = ktl.Service('apfucam')
    outfile    = ucam['OUTFILE']
    elapsed    = ucam['ELAPSED']
    obsnum     = ucam['OBSNUM']
    event      = ucam['EVENT']
    combo_ps   = ucam['COMBO_PS']
    ctalk      = ucam['CTALKTO']
    nerase     = ucam['NERASE']
    disp0sta   = ucam['DISP0STA']

    apfschedule= ktl.Service('apfschedule')

    guide      = ktl.Service('apfguide')
    counts     = guide['COUNTS']
    kcountrate     = guide['COUNTRATE']
    avg_fwhm   = guide['AVG_FWHM']

    motor      = ktl.Service('apfmot')
    decker     = motor['DECKERNAM']
#    deckerord  = motor['DECKERORD']
    dewarfoc   = motor["DEWARFOCRAW"]
    hatchpos   = motor["HATCHPOS"]
    ucampower  = motor['UCAMPOWER']
    gcam_power = motor['GCAMPOWER']

    eosgcam    = ktl.Service('eosgcam')
    fits3pre   = eosgcam('FITS3PRE')
    fits3dir   = eosgcam('FITS3DIR')
    save3d     = eosgcam('SAVE3D')
    fits2pre   = eosgcam('FITS2PRE')
    fits2dir   = eosgcam('FITS2DIR')
    save2d     = eosgcam('SAVE2D')
    gexptime   = eosgcam('GEXPTIME')
    sumframe   = eosgcam('SUMFRAME')

    apfmon     = ktl.Service('apfmon')
    ucamd0sta  = apfmon['UCAMDSTA0STA']

    apfminimon = ktl.Service('apfminimon')

    def __init__(self, task="example", test=False):
        """ Initilize the current state of APF. Setup the callbacks and monitors necessary for automated telescope operation."""
        # Set up the calling task that set up the monitor and if this is a test instance
        self.test = test
        self.task = task
        self.desired_outfile = None
        self.old_size = 0

        self.cwd        = os.getcwd()
        self.slowdown   = 0.0
        self.ncountrate = 0
        self.countrate = 0.0
        self.ccountrate = 0.0

        try:
            self.eosgcam['GENABLE'].write(True,binary=True)
        except:
            apflog("Cannot write True to eosgcam.GENABLE, issue with guider and/or dresden",level='error',echo=True)


        self.apfstas = []
        for n in range(1,8):
            kwnm = 'apfmon%dsta'  % (n)
            kw = self.apfminimon[kwnm]
            try:
                kw.monitor()
                kw.callback(self.mini_mon_mon)
                self.apfstas.append(kw)
            except Exception as e:
                apflog("Cannot monitor keyword %s: %s" % (kwnm,e),echo=True, level='warn')

        # Set the callbacks and monitors

        self.kcountrate.monitor()
        self.kcountrate.callback(self.countrate_mon)

        self.elapsed.monitor()

        self.obsnum.monitor()
        self.obsnum.callback(self.update_last_obs)

        self.event.monitor()
        self.event.callback(self.event_mon)

        self.hatchpos.monitor()
        self.ucampower.monitor()
        self.nerase.monitor()
        self.gcam_power.monitor()

        self.calsta.monitor() 
        self.focussta.monitor() 

        self.counts.monitor()
        self.vmag.monitor()
        self.ldone.monitor()
        self.counts.monitor()
        self.decker.monitor()
        self.avg_fwhm.monitor()
        self.dewarfoc.monitor()
        self.apfteqsta.monitor()
        self.metxfersta.monitor()

    def __str__(self):
        # Determine if the sun rising / setting check is working
        now = datetime.datetime.now()
        s = ''
        s += "At %s state of system is:\n" % str(now)
        s += "Slowdown = %5.2f x\n" % self.slowdown
        s += "countrate = %5.2g cts/s\n" % self.countrate
        s += "kcountrate = %5.2g cts/s\n" % self.kcountrate
        s += "ncountrate = %d frames \n" % self.ncountrate
        s += "elapsed = %5.2f sec \n" % self.elapsed
        s += "Guider camera power is %s\n" % ("ON" if self.gcam_power.binary else "OFF")

        rpid, rr = self.find_robot()

        if rr:
            s += "Robot is running as %s\n" % (rpid)
        else:
            s += "Robot is not running\n"

        stasum = ''
        for kw in self.apfstas:
            stasum += '%s %s ' % (kw['name'],kw['ascii'])

        s += 'APFmon is %s' % (stasum)

        return s

    ## callbacks - these are the callbacks used to measure quantities or flip states
    ## these are always running

    def ucam_dispatch_mon(self):
        if self.ucamd0sta['populated'] == False:
            return
        try:
            apfmon_stat = self.ucamd0sta['binary']
            if apfmon_stat == 4:
                # modify -s apfucam DISP0DWIM="ksetMacval DISP0STA READY"
                if self.disp0sta.read(binary=True,timeout=2) == 0:
                    self.ucam['DISP0DWIM'].write("ksetMacval DISP0STA READY")
        except:
            return

        return

    def count_mon(self, counts):
        if counts['populated'] == False:
            return
        try:
            cnts = counts['binary']
        except:
            return
        try:
            time = float(self.elapsed.read(binary=True,timeout=2))
        except:
            return

        try:
            self.ccountrate = cnts/time
        except:
            return


    def countrate_mon(self, kcountrate):
        if kcountrate['populated'] == False:
            return

        try:
            ctr = float(kcountrate['binary'])
        except:
            apflog("Cannot read apfguide.COUNTRATE",level='warn',echo=True)
            return
        self.countrate *=  (1.0*self.ncountrate)/(self.ncountrate+1)
        self.countrate += ctr/(self.ncountrate+1)
        self.ncountrate += 1
        return

    def event_mon(self, event):
        if event['populated'] == False:
            return

        try:
            eventval = event['binary']
        except:
            return

        if eventval == 0 or eventval == 7 :
            self.ncountrate = 0

        #if eventval == 4:
        #    self.set_tel_foc()

        try:
            cnts = float(self.counts.read(binary=True,timeout=2))
            time = float(self.elapsed.read(binary=True,timeout=2))
        except:
            return
        try:
            self.ccountrate = cnts/time
        except:
            return




    def apftask_mon(self, status):
        if status['populated'] == False:
            return
        try:
            status_val = status['binary']
        except:
            return

        hosts = dict()
        hosts['METSXFER'] = 'frankfurt.ucolick.org'
        hosts['APFTEQ'] = 'hamburg.ucolick.org'


        if status_val > 0:
            # apftask status values are $(TASKNAME)_status
            taskname_status = status['name'].lower().strip()
            taskname_split =  taskname_status.split("_")
            taskname = taskname_split[0]
            apflog("%s has status %s" % (taskname,status['ascii']),level='error',echo=True)

            # now we need the runhost

            runhost = hosts[taskname.upper()]

            master_runhost_keyword = 'MASTER_RUNHOST'
            current_host = ktl.read('apftask',master_runhost_keyword)

            apfcmd = os.path.join(LROOT,"bin/apf")
            restart = '%s restart %s' % (apfcmd,taskname)
            if runhost == current_host:
                cmdlist = restart.split()
            else:
                cmdlist = ["ssh", "-f", runhost, restart]

            try:
                p = subprocess.check_output(cmdlist,stderr=subprocess.STDOUT)
            except Exception as e:
                apflog("Cannot restart %s on %s: %s" % (taskname,runhost,e),level='error',echo=True)
                return
            apflog("%s should be restarted" % (taskname),echo=True)
        return


    def restart(self,name,host):
        apfcmd = os.path.join(LROOT,"bin/apf")
        restart = '%s restart %s' % (apfcmd,name)
        cmdlist = ["ssh", "-f", host, restart]
        try:
            p = subprocess.check_output(cmdlist,stderr=subprocess.STDOUT)
        except Exception as e:
            apflog("Cannot restart %s on %s: %s" % (name,host,e),level='error',echo=True)
        return


    def mini_mon_mon(self,sta,host="hamburg"):
        if sta['populated'] == False:
            return
        try:
            sta_val = sta['binary']
        except:
            return

        if sta_val > 3:
            # warning or higher
            nmsta = sta['name'].lower()
            name = nmsta[0:7] # this relies on the fact that all of the STA 
            # variables are serviceSTA and service is
            self.restart(name,host)
        return

    def apfmon_mon(self,sta,host="shadow"):
        if sta['populated'] == False:
            return
        try:
            sta_val = sta['binary']
        except:
            return

        if sta_val > 3:
            # warning or higher
            nmsta = sta['name'].lower()
            name = "apf" + nmsta[0:7]
            # this relies on the fact that all of the STA
            # variables are serviceSTA and service is
            self.restart(name,host)
        return

    def apftask_status_mon(self,sta):
        if sta['populated'] is False:
            return
        try:
            sta_val = sta['binary']
        except:
            return

        if sta_val >= 3:
            # exited
            nmsta = sta['name'].lower()
            name, _ = nmsta.split("_")
            pid = ktl.read('apftask',name + "_PID", binary=True)
            if pid > 0:
                # we have a contradiction
                # the process is not running but the pid is stale
                try:
                    ktl.write('apftask',name + "_PS_STATE", '')
                except:
                    apflog("Cannot write to apftask.%s_PS_STATE" % (name),level='error',echo=True)
                    pass
        return
    ## end of callbacks for monitoring stuff


    ## these are various methods, there are a LOT of them
    ##



    def init_keyword(self, keyword, value, timeout=None):
        """

        init_keyword(keyword, value, timeout=None)

        Writes the keyword to the value given a timeout.

        """

        success = False
        trials = 0
        while not success and trials < 5:
            # this is aimed at the guider
            # the keyword for the exposure time
            # cannot be set during an exposure
            # so we try a few times
            try:
                if timeout:
                    keyword.write(value,wait=True,timeout=timeout)
                else:
                    keyword.write(value,wait=True)
            except Exception as e:
                trials += 1
            else:
                success = True
        return success

    def init_gexptime(self):
        """
        init_gexptime()

        Sets guider exposure time to 1 second.

        """

        timeout = None
        if self.gexptime.binary > 1:
            timeout = self.gexptime.binary

        ret_val = self.init_keyword(self.gexptime, 1, timeout=timeout)
        return ret_val

    def init_sumframe(self):
        """
        init_sumframe()

        Sets sumframe to 1.

        """
        ret_val = self.init_keyword(self.sumframe, 1)
        return ret_val


    def init_guide_cam(self):
        """
        init_guide_cam()

        Sets the guider camera to sensible defaults.

        """

        if self.gcam_power.binary is False:
            return False

        ret_val = True
        self.save3d.write(False,binary=True)
        self.fits3pre.write('')
        if self.gexptime.read(binary=True) >= 1:
            ret_val = self.init_sumframe()
            if ret_val:
                ret_val = self.init_gexptime()

        else:
            ret_val = self.init_gexptime()
            time.sleep(0.2)
            # the guider is not very smart
            if ret_val:
                ret_val = self.init_sumframe()

        return ret_val

    def validate_UCAM_outputs(self):
        """
        validate_UCAM_outputs()

        Checks the output file from the UCAM to make sure it matches
        the values it should have.

        """
        if self.outfile.read() != self.desired_outfile:
            apflog("Output filename is %s and not the current date %s" % \
                   (self.desired_outfile, self.outfile),level='error',echo=True)
            self.outfile.write(self.outfile)

        if self.obsnum < self.robot["MASTER_LAST_OBS_UCSC"]:
            apflog("Output file number is %s which is less than the last logged value %s"\
                    % (self.obsnum, self.robot["MASTER_LAST_OBS_UCSC"]),level='error',echo=True)
            self.obsnum.write(self.robot["MASTER_LAST_OBS_UCSC"])

        return


    def set_observer_info(self, num=10000, name='Robot', owner='public'):
        """
        set_observer_info(num=10000, name='Robot', owner='public')

        This sets the apfucam.OBSERVER keyword,
        the apfucam.OBSNUM keyword,
        and the checkapf.OWNRHINT keyword.
        Also sets up other items just to be safe, such as
        the output directory, binning, and the file prefix.

        """
        if self.test: return
        apflog("Setting science camera parameters.")
        self.ucam('OBSERVER').write(name)
        self.apfschedule('OWNRHINT').write(owner)
        self.outfile.write(name)
        self.desired_outfile = name
        self.ucam('OUTDIR').write('/data/apf/')
        if num:
            self.obsnum.write(str(num))
        self.robot['UCAMLAUNCHER_UCAM_PCC'].write(0)

        bstr = "%d,%d" % (1,1)
        self.ucam['BINNING'].write(bstr)

        apflog("Updated science camera parameters:")
        apflog("Observer = %s" % self.ucam('OBSERVER').read(),echo=True)
        apflog("Ownrhint = %s" % self.apfschedule('OWNRHINT').read(),echo=True)
        apflog("Output directory = %s" % self.ucam('OUTDIR').read(),echo=True)
        apflog("Binning = %s" % self.ucam('BINNING').read(),echo=True)
        apflog("File prefix = %s" % self.outfile.read(), echo=True)
        apflog("Observation number = %s" % self.obsnum.read(), echo=True)

        return

    def instr_permit(self):
        """
        instr_permit()

        Do we have instrument permission and is the user kind correct?
        """
        try:
            while not self.instr_perm.read(binary=True,timeout=2) or \
                self.userkind.read(binary=True,timeout=2) != 3:
                apflog("Waiting for instrument permission to be true and userkind to be robotic")
                APFTask.waitfor(self.task, True, expression="$checkapf.INSTR_PERM = true", \
                                timeout=600)
                APFTask.waitfor(self.task, True, expression="$checkapf.USERKIND = robotic", \
                                timeout=600)
        except Exception as e:
            apflog("Cannot communicate with checkapf: %s" % (e),echo=True,level='alert')
            return False
        else:
            return True

    def turn_off_lamps(self):
        """
        turn_off_lamps()

        Does what it says on the tine. All lamps, even the ones not used.
        """
        for lamp in ("HALOGEN2","HALOGEN1","THORIUM1","THORIUM2"):
            try:
                rv = ktl.write("apfmot",lamp,"Off",wait=False)
            except Exception as e:
                apflog("Exception: %s" % (e),echo=True,level="alert")
                rv = False
            if rv is False:
                apflog("Cannot turn off lamp %s" % (lamp),echo=True,level="alert")
        return rv


    def write_stages(self,stagelist,component,state):
        """
        write_stages(stagelist, component, state)

        For every stage in stagelist, write the
        component to the value state.
        The timeout is 10 seconds.
        Returns True on success, False otherwise.

        """
        rv = True
        for stage in stagelist:
            curkwd = stage + component
            try:
                crv = ktl.write("apfmot",curkwd,state,wait=True,timeout=10)
                if crv is False:
                    rv = False
            except:
                rv = False
        return rv

    def enable_obs_inst(self):
        '''
        enable_obs_inst(self)
        This function enables the instrument for observing by setting certain
        stages to be on at the end of moves, while leaving others off at
        the end of moves.
        '''
        stagelist = ['CALMIRROR','CALSOURCE','IODINE','GUIDEFOC']
        rv1 = self.write_stages(stagelist,'MOE','Off')
        rv2 = self.write_stages(stagelist,'MOD','Pos')

        stagelist = ['ADC','DECKER','DEWARFOC']
        rv3 = self.write_stages(stagelist,'MOE','On')
        rv4 = self.write_stages(stagelist,'MOD','Pos')

        retval = rv1 and rv2 and rv3 and rv4
        return retval

    def enable_cal_inst(self):
        '''
        enable_cal_inst(self)
        This function enables the instrument for calibrations by setting certain
        stages to be on at the end of moves, while leaving others off at
        the end of moves. The difference between enable_obs_inst and this
        function is that the ADC is off. This is because the ADC is not
        use for calibrations.
        '''
        retval = True

        stagelist = ['ADC','CALMIRROR','CALSOURCE','IODINE','GUIDEFOC']
        rv1 = self.write_stages(stagelist,'MOE','Off')
        rv2 = self.write_stages(stagelist,'MOD','Pos')

        stagelist = ['DECKER','DEWARFOC']
        rv3 = self.write_stages(stagelist,'MOE','On')
        rv4 = self.write_stages(stagelist,'MOD','Pos')

        retval = rv1 and rv2 and rv3 and rv4
        return retval


    def disable_inst(self):
        '''
        disable_inst(self)
        This disables the instrument, turning off all stages except the
        dewar focus stage. That stage is under a gravity load and
        needs to be on to prevent the stage from sliding.
        '''
        stagelist = ['ADC','GUIDEFOC','CALMIRROR','CALSOURCE','DECKER','IODINE']
        rv1 = self.write_stages(stagelist,'MOE','Off')
        rv2 = self.write_stages(stagelist,'MOO','Off')
        rv3 = self.write_stages(['DEWARFOC'],'MOE','On')
        return rv1 and rv2 and rv3

    def turn_off_inst(self):
        '''
        turn_off_inst(self)
        Turns off all motor stages.
        '''
        stagelist = ['ADC','GUIDEFOC','CALMIRROR','CALSOURCE','IODINE','DECKER','DEWARFOC']
        rv1 = self.write_stages(stagelist,'MOE','Off')
        rv2 = self.write_stages(stagelist,'MOO','Off')

        return rv1 and rv2


    def hatch_correct(self):
        '''
        hatch_correct(self)

        The hatch occasionally fails to open, leaving it in an unknown state.
        This function tries to fix that by closing the hatch and then
        opening the hatch again.
        '''
        if self.hatchpos['populated'] == False:
            return

        try:
            curval = self.hatchpos.binary
        except Exception as e:
            apflog("Exception in hatch_correct: %s %s" % (type(e), e), level='error')
            return

        if curval == 0:
            self.hatchpos.write(2,wait=False)
            if self.hatchpos.waitFor("==2",timeout=10):
                self.hatchpos.write(1,wait=False)
                self.hatchpos.waitFor("==1",timeout=10)
            else:
                return False

        return True

    def status_clear(self):
        """
        status_clear()

        Clears the PS status of any APFTask where
        the PS status does not match the actual status.
        """
        for kw in (self.slewsta, self.calsta, self.focussta, \
                   self.shuttersta, self.opensta, self.closesta,\
                    self.focustelsta):
            #self.apftask_status_mon(kw)
            pass


    def focusinstr(self, log_error_level='Alert'):
        """
        focusinstr()

        A wrapper for the instrument focus method.
        This ensures that the permissions are available,
        turns on various stages, makes sure that the UCAM is configured,
        and checks the output of the focus.

        """
        self.instr_permit()
        rv = self.enable_cal_inst()
        if rv is False:
            try:
                ip = self.checkapf['INSTR_PERM'].read(timeout=2)
            except:
                ip = 'Unknown'
            apflog("Cannot enable instrument to move stages but instr_perm is %s" % (ip), level='alert',echo=True)
            return rv
        try:
            owner = self.apfschedule('OWNRHINT').read(timeout=10)
        except Exception as e:
            apflog("Cannot communicate with apfschedule %s" % (e), level='alert',echo=True)
        else:
            self.apfschedule('OWNRHINT').write('public')

        self.validate_UCAM_outputs()

        lastfocus_dict = APFTask.get("focusinstr", ["lastfocus","nominal","useref"])
        if float(lastfocus_dict["lastfocus"]) > DEWARMAX or float(lastfocus_dict["lastfocus"]) < DEWARMIN:
            lastfocus_dict["lastfocus"] =  lastfocus_dict["nominal"]
        if bool(lastfocus_dict["useref"]) is False:
            self.robot['FOCUSINSTR_USEREF'].write(True,binary=True)

        result, msg = self.run_focusinstr()
        if result:
            apflog(msg,echo=True)
        else:
            if not self.instr_perm.read(binary=True):
                msg += "Focusinstr has failed because of loss of instrument permission, waiting for return"
                wait = True
            else:
                wait = False
            # this complication is just to ensure the logging happens before waiting for permissions
            apflog(msg,echo=True,level=log_error_level)
            if wait:
                self.instr_permit()

        dewarfocraw = self.dewarfoc.read(binary=True)

        if  (dewarfocraw > DEWARMAX or dewarfocraw < DEWARMIN):
            apflog("Focusinstr has failed, result = %s, current focus is value = %d, and last value was %s." % ( str(result),dewarfocraw,lastfocus_dict["lastfocus"]), level='error', echo=True)
            APFLib.write("apfmot.DEWARFOCRAW", lastfocus_dict["lastfocus"])
            return False

        try:
            self.apfschedule('OWNRHINT').write(owner,timeout=10)
        except Exception as e:
            apflog("Cannot communicate with apfschedule %s" % (e), level='alert',echo=True)
        try:
            self.decker.write("W (1.00:3.0)",wait=False)
            self.decker.waitFor(" == 'W (1.00:3.0)'",timeout=120)
        except Exception as e:
            apflog("Cannot communicate with apfmot.DECKERNAM %s" % (e), level='alert',echo=True)

        return result

    def calibrate(self, script, time):
        """
        calibrate(script, time)

        Runs the calibrate shell script with the specified script as
        the option.
        The time variable must be "pre" or "post" which is also passed
        to the calibrate script.

        """
        self.validate_UCAM_outputs()

        s_calibrate = os.path.join(SCRIPTDIR,"calibrate")
        if self.test:
            apflog("Test Mode: calibrate %s %s." % (script, time))
            APFTask.waitFor(self.task, True, timeout=10)
            return True

        if time != 'pre' and time != 'post':
            apflog("Couldn't understand argument %s, nothing was done." % time)
            return False

        rv = self.enable_cal_inst()
        if rv is False:
            try:
                ip = self.instr_perm.read()
            except:
                ip = 'Unknown'
            apflog("Cannot enable instrument to move stages but instr_perm is %s" % (ip), level='alert',echo=True)
            return rv

        try:
            APFLib.write("apfmot.DEWARFOCRAW",ktl.read("apftask","FOCUSINSTR_LASTFOCUS",binary=True))
        except:
            apflog("Cannot read the last best fitting focus value or write the dewar focus value", level='error')
        if self.dewarfoc > DEWARMAX or self.dewarfoc < DEWARMIN:
            apflog("Warning: The dewar focus is currently %d. This is outside the typical range of acceptable values." % (self.dewarfoc), level = "error", echo=True)
            return False
        apflog("Running calibrate %s %s" % (script, time), level = 'info')

        try:
            owner = self.apfschedule('OWNRHINT').read(timeout=10)
        except Exception as e:
            apflog("Cannot communicate with apfschedule %s" % (e), level='alert',echo=True)
        else:
            self.apfschedule('OWNRHINT').write('public')

        cmd = '%s %s %s' % (s_calibrate,script, time)
        result, code = apftask_do(cmd,debug=True,cwd=os.getcwd())
        if not result:
            apflog("%s %s failed with return code %d" % (s_calibrate, script, code),echo=True)
        expression="($apftask.CALIBRATE_STATUS != 0) and ($apftask.CALIBRATE_STATUS != 1) "
        if not APFTask.waitFor(self.task,True,expression=expression,timeout=30):
            apflog("%s %s failed to exit" % (s_calibrate,script),echo=True)

        try:
            self.apfschedule('OWNRHINT').write(owner,timeout=10)
        except Exception as e:
            apflog("Cannot communicate with apfschedule %s" % (e), level='alert',echo=True)

        return result

    def run_focusinstr(self,flags="-b"):
        """
        run_focusinstr(flags="-b")

        Runs the instrument focus routine with flags option.
        Flags for the instrument focus script as passed to the
        script, and are not checked. The default starts the focus
        from the beginning, ignoring if the script failed in the
        previous attempt.

        This is called by the method focusinstr()
        """

        msg = ""
        if self.test:
            APFTask.waitFor(self.task, True, timeout=10)
            msg = "Test Mode: Would be running focusinstr."
            return True, msg

        supplies = ('PS1_48V_ENA', 'PS2_48V_ENA')
        for keyword in supplies:
            try:
                value = self.motor[keyword].read(binary=True,timeout=2)
                if value != 1:
                    self.motor[keyword].write('Enabled', wait=False)
            except Exception as e:
                apflog("Cannot read status of PS's:  %s"  % e,level='alert', echo=True)
                return False, "Cannot read status of PS's:  %s" % (e)

        apflog("Running focusinstr routine.",echo=True)

        execstr = " ".join(['focusinstr',flags])
        cmd = os.path.join(SCRIPTDIR,execstr)
        result, code = apftask_do(cmd,cwd=os.getcwd())

        expression="($apftask.FOCUSINSTR_STATUS != 3)"
        if ktl.waitFor(expression=expression, timeout=.1):
            try:
                resultd = APFTask.get('FOCUSINSTR',('LASTFOCUS','PHASE'))
                if resultd['PHASE'] == 'Cleanup':
                    msg += 'focusinstr failed in or after cleanup, proceeding with value %s ' % (str(resultd['LASTFOCUS']))
                    result = True
                else:
                    msg += 'focusinstr failed in %s, focusinstr needs to be rerun ' % (resultd['PHASE'])
            except Exception as e:
                msg += "focusinstr failed, exited with %s: %s %s" % (self.robot['focusinstr_status'].read(),type(e),e)
                result = False

        expression="($apftask.FOCUSINSTR_MEASURED == 1)"
        if not ktl.waitFor(expression=expression, timeout=1):
            msg += "focusinstr fit to the data either failed or did not occur "
            result = False

        return result, msg



    def update_last_obs(self,obsnum):
        """ If the last observation was a success,
        this function updates the file storing the last
        observation number and the hit_list which is
        required by the dynamic scheduler."""

        if obsnum['populated']:
            if obsnum >= 10000:
                APFLib.write(self.robot["MASTER_LAST_OBS_UCSC"], obsnum)

        return


    def start_robot(self,observation=None,skip=False,raster=False):
        """Start an instance of scriptobs. Returns the result from subprocess.Popen()."""
        # For running in test mode
        if self.test:
            apflog("Would start robot",echo=True)
            if observation is not None:
                apflog("Would be taking observation in starlist %s" % observation,echo=True)
            APFTask.waitFor(self.task, True, timeout=10)
            return None

        # Make sure the telescope autofocus is enabled
        APFLib.write(self.autofoc, "robot_autofocus_enable")
        chk_foc = '$apftask.SCRIPTOBS_AUTOFOC == robot_autofocus_enable'
        result = APFTask.waitFor(self.task, False, chk_foc, timeout=60)
        if not result:
            apflog("Error setting scriptobs_autofoc", level='error',echo=True)
            return None

        # Make sure APFTEQ is in night mode for observations
        if self.teqmode.read() != 'Night':
            self.set_apfteq_mode('Night')

            # Check the instrument focus for a reasonable value
        if self.dewarfoc > DEWARMAX or self.dewarfoc < DEWARMIN:
            lastfit_dewarfoc = ktl.read("apftask","FOCUSINSTR_LASTFOCUS",binary=True)
            log_str = "Warning: The dewar focus is currently %d. " % (self.dewarfoc)
            log_str += "This is outside the typical range of acceptable values."
            log_str += "Resetting to last derived value %d" % (lastfit_dewarfoc)
            apflog(log_str, level = "error", echo=True)
            APFLib.write("apfmot.DEWARFOCRAW",lastfit_dewarfoc)

        # check on weirdness for UCAM host post-reboot
        self.ucam_dispatch_mon()

        telstate = self.tel['TELSTATE'].read()
        if telstate == 'Disabled':
            rv, _ = apftask_do(os.path.join(SCRIPTDIR,"slew --hold"))
            if not rv:
                return None
        # Start scriptobs

        outfile = open("robot.log", 'a')
        if raster:
            args = ['/home/holden/src/scriptobs_offset']
        else:
            if skip:
                args = [os.path.join(SCRIPTDIR,'scriptobs'), '-dir', os.getcwd(),'-skip']
            else:
                args = [os.path.join(SCRIPTDIR,'scriptobs'), '-dir', os.getcwd()]


        p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=outfile, stderr=outfile)

        return p

    def find_robot(self):
        """Trys to find a running instance of scriptobs.
            Returns the PID along with a boolean representing
            if the robot was succesfully found."""
        rpid = self.robot['SCRIPTOBS_PID'].read(binary=True)
        if rpid == '' or rpid == -1:
            return rpid, False

        return rpid, True

    def kill_robot(self, now=False):
        """ In case during an exposure there is a need to stop the robot and close up."""
        apflog("Terminating scriptobs")
        if now:
            apflog("Abort exposure, terminating scriptobs now.")
        else:
            if not self.ucam['EVENT_STR'].read() == "ControllerReady":
                apflog("Waiting for current exposure to finish.")
                time_out_val = self.ucam['REMAINING'].read(binary=True)+2
                self.ucam['EVENT_STR'].waitfor(" = ReadoutBegin", timeout=time_out_val)
                self.sop.waitfor(" = Input", timeout=2)

        ripd, running = self.find_robot()
        if running:
            apflog("Killing scriptobs %s" % (str(ripd)))
            try:
                APFLib.write(self.robot['SCRIPTOBS_CONTROL'], "abort")
            except Exception as e:
                errstr = "Cannot abort scriptobs: %s" % (e)
                apflog(errstr,level="Warn",echo=True)

    def turn_on_ucam(self):

        self.ucampower.write('On',wait=False)
        rv = self.ucampower.waitFor('== On',timeout=30)
        APFTask.wait(self.task, True, timeout=5)
        if rv is False:
            apflog('Cannot power on UCam',level='Alert',echo=True)
            return False

        ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","stop")
        APFTask.wait(self.task, True, timeout=1)
        ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","run")
        APFTask.wait(self.task, True, timeout=1)
        rv = ktl.read("apftask","UCAMLAUNCHER_UCAM_STATUS",binary=True)
        if rv == 0:
            apflog("UCAM software did not start, trying again",echo=True,level='Warn')
            ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","stop")
            APFTask.wait(self.task, True, timeout=1)
            ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","run")
            APFTask.wait(self.task, True, timeout=1)
            rv = ktl.read("apftask","UCAMLAUNCHER_UCAM_STATUS",binary=True)
            if rv == 0:
                apflog("UCAM software did not start",echo=True,level='Alert')
                return False
        return True

    def test_bias(self):
        """
        test_bias()

        Takes a single test bias exposure.
        The function will turn on the UCAM and start the UCAM software if needed.
        The output will be test_#.fits, where # is an integer that increments to
        not overwrite the previous bias if it exists.
        """

        if self.ucampower is False:
            rv = self.turn_on_ucam()
            if rv is False:
                apflog("Cannot power on UCAM",level='alert')
                return False

        apfschedule = ktl.Service('apfschedule')
        # check if the focusinstr or calibrate tasks are already running
        if ktl.read('apftask','FOCUSINSTR_PID',binary=True) > 0:
            return False
        if ktl.read('apftask','CALIBRATE_PID',binary=True) > 0:
            return False
        if ktl.read('apftask','SCRIPTOBS_PID',binary=True) > 0:
            return False
        # create exposure object
        exp = Exposure.Exposure(0,"bias",count=1,record="yes",parent=self.task,dark=True)

        combval = exp.comb.binary
        # Is the UCAM ok?
        if combval > 0:
            return False

        # read original values and save them
        try:
            outdir = exp.apfucam['OUTDIR'].read()
            orignam = exp.outfile.read()
            orignum = exp.obsnum.read()
        except:
            return False


        ofn  = 'test_'
        obsn = 1
        ffn = "%s%d.fits" % (ofn,obsn)
        fpath = os.path.join(outdir,ffn)
        while os.path.exists(fpath):
            obsn += 1
            ffn = "%s%d.fits" % (ofn,obsn)
            fpath = os.path.join(outdir,ffn)

        # write test values
        try:
            exp.outfile.write(ofn)
            exp.obsnum.write(str(obsn))
            apfschedule['OWNRHINT'].write('unknown')
        except:
            return False

        # final file name
        ffn = exp.outfile.read() + exp.obsnum.read() + '.fits'
        fpath = os.path.join(outdir,ffn)

        apflog("Taking a test bias image called %s" % (ffn),echo=True)
        # take two pictures
        try:
            c = exp.expose(waitlast=True)
        except:
            rv = False

        apfschedule['OWNRHINT'].write('public')

        # the number of exposures (c) should be 1
        if os.path.exists(fpath) and c == 1:
            rv = True
        else:
            rv = False
        apflog("File located at %s is %s" % (fpath,str(os.path.exists(fpath))),echo=True)

        # restore original values
        try:
            exp.outfile.write(orignam)
            exp.obsnum.write(orignum)
        except:
            rv = False
        return rv



    def ucam_power_cycle(self, fake=False):
        """
        Power cycles the UCAM power supply.
        This forks a script execution.
        """
        if fake:
            apflog("would have executed %s" % (os.path.join(SCRIPTDIR,"robot_power_cycle_ucam")))
            return True

        val = subprocess.call(os.path.join(SCRIPTDIR,"robot_power_cycle_ucam"))
        if val > 0:
            apflog("power cycle of UCAM failed",level='alert')
            return False

        return True

    def ucam_reboot(self,fake=False):
        """
        ucam_reboot

        Reboots the UCAM host using the UCAMLAUNCHER

        """
        if fake:
            apflog("Would have rebooted UCAM host ",echo=True)
            return True

        apftask = ktl.Service('apftask')
        command = apftask['UCAMLAUNCHER_UCAM_COMMAND']
        ucamstat = apftask['UCAMLAUNCHER_UCAM_STATUS']

        try:
            command.write("Stop")
            apflog("Stopping UCAM software",echo=True)

            self.combo_ps.waitFor(" == MissingProcesses",timeout=30)
            command.write("Reboot")
            apflog("Rebooting UCAM host",echo=True)

        except:
            apflog("UCAM status bad, cannot restart",level='alert')
            return False

        # we cannot wait on the UCAM or UCAMLAUNCHER keywords,
        # the dispatcher is down
        APFTask.wait(self.task,False,timeout=240)

        command = apftask['UCAMLAUNCHER_UCAM_COMMAND']
        ucamstat = apftask['UCAMLAUNCHER_UCAM_STATUS']

        try:
            command.write("Run")
            ucamstat.waitFor(" == running",timeout=300)
            apflog("UCAM software running",echo=True)

        except:
            apflog("UCAM status bad, cannot restart",level='alert')
            return False

        nv = self.combo_ps.waitFor(" == Ok",timeout=30)
        apflog("UCAM software combo_ps keyword OK",echo=True)
        if nv:
            return nv

        apflog("UCAM host reboot failure, combo_ps still not ok" , level="alert", echo=True)

        self.obsnum.monitor()
        self.obsnum.callback(self.update_last_obs)

        self.event.monitor()
        self.event.callback(self.event_mon)


    def ucam_restart(self, fake=False):
        """
        Restarts the UCAM software using the ucamlauncher.
        If it is fails, reboots UCAM host.
        """
        # modify -s apftask UCAMLAUNCHER_UCAM_COMMAND=Stop
        if fake:
            # would have restarted software
            apflog("Would have restarted UCAM software ")
            return True

        try:
            apflog("Stop and restarting UCAM software",echo=True)
            ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","stop")
            if self.combo_ps.waitFor(" == MissingProcesses",timeout=30):
                ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","run")
                nv = self.combo_ps.waitFor(" == Ok",timeout=30)
                if nv:
                    return True    
                apflog("UCAM  restart failure, combo_ps still not ok" , level="error", echo=True)
        except:
            apflog("UCAM status bad, cannot restart",level='alert')
            return False

        rv = self.ucam_reboot()

        return rv



    def ucam_status(self, fake=False):
        """
        Checks the status of the UCAM software.

        If needed, it will power cycle the UCAM controller or
        reboot the UCAM host.

        """
        if self.ctalk.read(binary=True) > 0:
            rv = self.ucam_power_cycle(fake=fake)
            return rv

        if self.combo_ps.read(binary=True) > 0:
            # brains!
            rv = self.ucam_restart(fake=fake)
            return rv

        try:
            ucamsta0 = self.ucam['DISP0STA'].read(binary=True)
            ucamsta1 = self.ucam['DISP1STA'].read(binary=True)
        except Exception as e:
            apflog('apfucam.DISPSTA failure, apfucam likely not running: %s' % (e),echo=True,level='Alert')
            rv = self.ucam_restart(fake=fake)
            return rv
        else:
            if ucamsta1 > 0 and ucamsta0 > 0:
                # Things are still starting up
                if ucamsta0 > 2:
                    # failure to connect
                    rv = self.ucam_restart(fake=fake)
                    return rv

                rv = APFTask.waitfor(self.task, True, expression="$apfucam.DISP0STA = 0 & $apfucam.DISP1STA = 0", timeout=600)
                return rv

        ucamlaunch_sta = self.robot['UCAMLAUNCHER_UCAM_STATUS'].read(binary=True)
        if ucamlaunch_sta == 0:
            try:
                self.robot['UCAMLAUNCHER_UCAM_COMMAND'].write('Run')
            except Exception as e:
                apflog('Failed when writing apftask.UCAMLAUNCHER_UCAM_COMMAND to Run: %s' % (e),echo=True,level='Crit')
                return False

        return True
    
    def ucam_watchdir(self):
        '''
        
        ucam_watchdir()

        This simply watches the data directory to make sure that it increases in size
        when exposures are taken. 
        '''

        outdir = self.ucam['OUTDIR'].read()
        cmd_str = ['du','-s',outdir]
        try:
            dir_size = subprocess.check_output(cmd_str)
        except Exception as e:
            apflog("Cannot compute size of %s: %s %s" % (outdir, type(e), e))
            return False

        if dir_size > self.old_size:
            self.old_size = dir_size
            return True
        return False

if __name__ == '__main__':
    print("Testing telescope monitors, grabbing and printing out current state.")

    task = 'example'

    APFTask.establish(task, os.getpid())
    apf = APF(task=task,test=False)

    # Give the monitors some time to start up
    APFTask.waitFor(task, True,timeout=2)

    print(str(apf))

    while True:
        APFTask.wait(task,True,timeout=10)
        print(str(apf))
        apf.robot['FOCUSTEL_STARTFOCUS'].write(apf.pred_tel_focus())
