from __future__ import print_function
# Class definition for an APF object which tracks the state of the telescope.

import subprocess
import time
import os
import os.path
import math
from datetime import datetime, timedelta

import numpy as np

try:
    from apflog import apflog
    import ktl
    import APF as APFLib
    import APFTask
    import Exposure
except:
    from fake_apflog import *


windlim = 40.0
slowlim = 100
WINDSHIELD_LIMIT = 10. # mph at the APF
FOCUSTIME = 3600. # minimum time before checking telescope focus
TEMP_LIMIT = 35. # deg F at the APF
wxtimeout = timedelta(seconds=1800)
SUNEL_HOR = -3.2
DEWARMAX = 8650
DEWARMIN = 8350
TELFOCUSMIN = -0.00096
TELFOCUSMAX = -0.00060
#TELFOCUSTYP = -0.83529
TELFOCUSTYP = -0.76529
TELFOCUSMAXOFF = 0.00002

if "LROOT" in os.environ:
    LROOT = os.environ["LROOT"]
else:
    LROOT = '/usr/local/lick'
SCRIPTDIR = os.path.join(LROOT,'bin/robot/')


def apftaskDo(cmd, debug=True, cwd='./'):
    newcmd = "apftask do %s" % (cmd)
    rv, retcode = cmdexec(newcmd, debug=debug, cwd=cwd)
    return rv, retcode


def cmdexec(cmd, debug=False, cwd='./'):
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

    cwd        = os.getcwd()
    slowdown   = 0.0
    ncountrate = 0
    countrate = 0.0
    ccountrate = 0.0

    # Initial Wind conditions
    wslist = []

    # Initial temps


    dewlist = []
    dewTooClose = False

    # KTL Services and Keywords
    tel        = ktl.Service('eostele')
    sunel      = tel('SUNEL')
    ael        = tel('AEL')
    aaz        = tel('AAZ')
    aafocus    = tel('AAFOCUS')
    focus      = tel('FOCUS')
    faenable   = tel('FAENABLE')

    dome       = ktl.Service('eosdome')
    rspos      = dome('RSCURPOS')
    fspos      = dome('FSCURPOS')
    shclosed   = dome('SHCLOSED')

    eostdio    = ktl.Service('eostdio')
    mcopen     = eostdio('MCOPEN')

    checkapf   = ktl.Service('checkapf')
    ok2open    = checkapf('OPEN_OK')
    userkind   = checkapf('USERKIND')
    dmtimer    = checkapf('DMTIME')
    whatsopn   = checkapf('WHATSOPN')
    mv_perm    = checkapf('MOVE_PERM')
    instr_perm = checkapf('INSTR_PERM')
    chk_close  = checkapf('CHK_CLOSE')

    apfmet     = ktl.Service('met3apf')
    wx         = apfmet('M5WIND')
    airtemp    = apfmet('M5OUTEMP')
    down       = apfmet('M5DOWN')
    altwx      = apfmet('M3WIND')

    eosti8k    = ktl.Service('eosti8k')
    m2tempkw   = eosti8k('TM2CSUR')
    m2airkw    = eosti8k('TM2CAIR')
    m1tempkw   = eosti8k('TM1S210')
    taveragekw = eosti8k('TAVERAGE')
    t045kw     = eosti8k('TTRUS045')
    t135kw     = eosti8k('TTRUS135')
    t225kw     = eosti8k('TTRUS225')
    t315kw     = eosti8k('TTRUS315')


    eoscool    = ktl.Service('eoscool')
    dewpt      = eoscool('DEWPAVG3')
    temp3now   = eoscool('TEMPNOW3')
    temp4now   = eoscool('TEMPNOW4')

    robot        = ktl.Service('apftask')
    vmag         = robot['SCRIPTOBS_VMAG']
    ldone        = robot['SCRIPTOBS_LINES_DONE']
    line         = robot['SCRIPTOBS_LINE']
    sop          = robot['SCRIPTOBS_PHASE']
    message      = robot['SCRIPTOBS_MESSAGE']
    autofoc      = robot["SCRIPTOBS_AUTOFOC"]
    slew_allowed = robot['SLEW_ALLOWED']
    observed     = robot['SCRIPTOBS_OBSERVED']
    apfteqsta    = robot['APFTEQ_STATUS']
    metxfersta   = robot['METSXFER_STATUS']
    calsta       = robot['CALIBRATE_STATUS']
    focussta     = robot['FOCUSINSTR_STATUS']
    ucamcmd      = robot['UCAMLAUNCHER_UCAM_COMMAND']
    lastopen     = robot['OPENUP_LAST_SUCCESS']
    msg = ""

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

    apfteq     = ktl.Service('apfteq')
    teqmode    = apfteq['MODE']

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
                kw.callback(self.miniMonMon)
                self.apfstas.append(kw)
            except Exception as e:
                apflog("Cannot monitor keyword %s: %s" % (kwnm,e),echo=True, level='warn')

        # Set the callbacks and monitors
        self.ok2open.monitor()
        self.ok2open.callback(self.okmon)

        self.dmtimer.monitor()
        self.dmtimer.callback(self.dmtimemon)

        self.kcountrate.monitor()
        self.kcountrate.callback(self.countratemon)

        self.elapsed.monitor()

        self.obsnum.monitor()
        self.obsnum.callback(self.updateLastObs)

        self.event.monitor()
        self.event.callback(self.eventmon)

        self.hatchpos.monitor()
        self.ucampower.monitor()
        self.nerase.monitor()
        self.gcam_power.monitor()

        self.down.monitor()
        self.whatsopn.monitor()

        self.mon_lists = dict()
        self.avg_lists = dict()

        for kw in (self.m1tempkw,self.m2tempkw,self.m2airkw,self.taveragekw,\
                   self.t045kw,self.t135kw,self.t225kw,self.t315kw,self.temp3now,\
                    self.temp4now,self.wx,self.altwx,self.airtemp):
            self.mon_lists[kw['name']] = []
            self.avg_lists[kw['name']] = None
            kw.monitor()
            kw.callback(self.listMon)
            kw.read()

        self.dewpt.monitor()
        self.dewpt.callback(self.dewPtMon)

        self.counts.monitor()
        self.teqmode.monitor()
        self.vmag.monitor()
        self.ldone.monitor()
        self.counts.monitor()
        self.decker.monitor()
        self.avg_fwhm.monitor()
        self.dewarfoc.monitor()
        self.mv_perm.monitor()
        self.chk_close.monitor()
        self.slew_allowed.monitor()

        self.sunel.monitor()
        self.aaz.monitor()
        self.ael.monitor()
        self.fspos.monitor()
        self.rspos.monitor()
        self.focus.monitor()
        self.aafocus.monitor()
        self.faenable.monitor()

        self.apfteqsta.monitor()
        self.apfteqsta.callback(self.apftaskMon)

        self.metxfersta.monitor()
        self.metxfersta.callback(self.apftaskMon)

        self.calsta.monitor()
        self.focussta.monitor()

        self.lastopen.monitor()

        # Grab some initial values for the state of the telescope

        self.wx.read()
        self.altwx.read()
        self.dewpt.read()
        self.counts.read()
        self.ok2open.read()
        self.avgtemps = np.asarray([self.avg_lists[nm] for nm in \
                                     ('TM1S210','TM2CSUR','TAVERAGE',\
                                      'TM2CAIR','TEMPNOW3','TEMPNOW4')])



    def __str__(self):
        # Determine if the sun rising / setting check is working
        now = datetime.now()
        self.avgTelTemps()
        s = ''
        s += "At %s state of telescope is:\n" % str(now)
        s += "Sun elevation = %4.2f %s\n" % (self.sunel, "Rising" if self.sunRising() else "Setting")
        s += "Telescope -- AZ=%4.2f  EL=%4.2f \n" % (self.aaz, self.ael)
        s += "Front/Rear Shutter=%4.2f / %4.2f\n"%(self.fspos, self.rspos)
        s += "Wind = %3.1f mph (APF) %3.1f mph (Shane) \n" % (np.average(self.mon_lists['M5WIND']),np.average(self.mon_lists['M3WIND']))
        s += "Slowdown = %5.2f x\n" % self.slowdown
        s += "Last open time = %.2f sec\n" % (self.lastopen.binary)
        s += "Time since opening = %6.2f sec\n" % (time.time() - self.lastopen.binary)
        s += "Msg = %s\n" % self.msg
        s += "countrate = %5.2g cts/s\n" % self.countrate
        s += "kcountrate = %5.2g cts/s\n" % self.kcountrate
        s += "ncountrate = %d frames \n" % self.ncountrate
        s += "elapsed = %5.2f sec \n" % self.elapsed
        s += "M1 = %5.2f deg C M2 = %5.2f deg C Tel Avg = %5.2f deg C M2 Air = %5.2f deg C FCU3 = %5.2f deg C FCU4 = %5.2f deg C\n" % tuple(self.avgtemps)
        s += "Dewpt = %5.2f deg C Teq Mode - %s\n" % (np.average(self.dewlist),self.teqmode)
        s += "Too close to the dewpoint? = %s\n" % self.dewTooClose
        s += "Guider camera power is %s\n" % ("ON" if self.is_gcam_power else "OFF")
        s += "M2 Focus Value = % 4.3f\n" % (float(self.aafocus['binary'])*1000.0)
        s += "M2 Focus Value = % 4.3f (focus kwd)\n" % (float(self.focus['binary'])*1000.0)
        s += "Preferred M2 Focus Value =  % 4.3f\n" % (float(self.predTelFocus())*1000.0)
        s += "Okay to open = %s -- %s\n" % (self.openOK['ascii'], self.checkapf['OPREASON'].read() )
        s += "Current Weather = %s\n" % self.checkapf['WEATHER'].read()

        isopen, what = self.isOpen()
        if isopen:
            s += "Currently open: %s\n" % what
        else:
            s += "Not currently open\n"
        ripd, rr = self.findRobot()
        if rr:
            s += "Robot is running as %s\n" % (ripd)
        else:
            s += "Robot is not running\n"
        focval = self.setAutofocVal()
        s += "Focus value for scriptobs = %d\n" % focval
        windshield = self.updateWindshield("auto")
        s += "Windshield state = %s\n" % windshield

        stasum = ''
        for kw in self.apfstas:
            stasum += '%s %s ' % (kw['name'],kw['ascii'])

        s += 'APFmon is %s' % (stasum)

        return s

    ## callbacks - these are the callbacks used to measure quantities or flip states
    ## these are always running

    def ucamdispatchmon(self):
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

    def countmon(self, counts):
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


    def countratemon(self, kcountrate):
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

    def eventmon(self, event):
        if event['populated'] == False:
            return

        try:
            eventval = event['binary']
        except:
            return

        if eventval == 0 or eventval == 7 :
            self.ncountrate = 0

        #if eventval == 4:
        #    self.setTelFoc()

        try:
            cnts = float(self.counts.read(binary=True,timeout=2))
            time = float(self.elapsed.read(binary=True,timeout=2))
        except:
            return
        try:
            self.ccountrate = cnts/time
        except:
            return

    # Callback for ok2open permission
    # -- Check that if we fall down a logic hole we don't error out
    def okmon(self,ok2open):
        if ok2open['populated'] == False:
            return
        try:
            ok = ok2open # historical
        except Exception as e:
            apflog("Exception in okmon for checkapf.OPEN_OK: %s" % (e), level='error')
            return
        try:
            if self.mv_perm.read(binary=False) == False:
                ok = False
        except Exception as e:
            apflog("Exception in okmon for checkapf.MOVE_PERM: %s" % (e), level='error')
            return
        try:
            if not self.userkind.read(binary=True) == 3:
                ok = False
        except Exception as e:
            apflog("Exception in okmon checkapf.USERKIND: %s" % (e), level='error')
            return
        self.openOK = ok
        return

    def listMon(self,keyword):
        if keyword['populated'] == False:
            return

        name = keyword['name']

        try:
            curval = float(keyword['binary'])
        except Exception as e:
            apflog("Exception in listMon: %s" % (e), level='error')
            return

        if self.mon_lists[name] == []:
            self.mon_lists[name] = [curval]*300
            self.avg_lists[name] = curval
        else:
            self.mon_lists[name].append(curval)
            self.mon_lists[name] = self.mon_lists[name][-300:]
            self.avg_lists[name] = np.average(self.mon_lists[name])

        return

    # Callback for Deadman timer
    def dmtimemon(self,dmtime):
        if dmtime['populated'] == False:
            return
        try:
            self.dmtime = dmtime
        except Exception as e:
            apflog("Exception in dmtimemon: %s %s" % (type(e), e), level='error')

        return

    def dewPtMon(self,dew):
        if dew['populated'] == False:
            return
        try:
            dewpt = float(dew['binary'])
        except:
            return

        if self.dewlist == []:
            self.dewlist = [dewpt]
        else:
            self.dewlist.append(dewpt)
            if len(self.dewlist) >= 300:
                self.dewlist = self.dewlist[-300:]

                dewlist = np.asarray(self.dewlist)

                curdew = np.average(dewlist)
                curm2 = np.average(np.asarray(self.mon_lists['TM2CSUR']))
                curm2air = np.average(np.asarray(self.mon_lists['TM2CAIR']))

                if self.dewTooClose:
                    if curm2air - curdew > 3 or curm2 - curdew > 5:
                        self.dewTooClose = False
                else:
                    if curm2air - curdew < 2 or curm2 - curdew < 4:
                        self.dewTooClose = True
                        apflog("M2 temperatures too close to dew point",echo=True,level='warning')

        return


    def apftaskMon(self,status):
        if status['populated'] == False:
            return
        try:
            status_val = status['binary']
        except:
            return

        hosts = dict()
        hosts['METSXFER'] = 'frankfurt.ucolick.org'
        hosts['APFTEQ'] = 'hamburg.ucolick.org'


        if status_val > 2:
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


    def miniMonMon(self,sta,host="hamburg"):
        if sta['populated'] == False:
            return
        try:
            sta_val = sta['binary']
        except:
            return

        if sta_val > 3:
            # warning or higher
            nmsta = sta['name'].lower()
            name = nmsta[0:7] # this relies on the fact that all of the STA variables are serviceSTA and service is
            self.restart(name,host)
        return

    def apfMonMon(self,sta,host="shadow"):
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

    ## end of callbacks for monitoring stuff


    ## these are various methods, there are a LOT of them
    ##

    def sunRising(self):
        # the sun also rises
        now = datetime.now()
        if now.strftime("%p") == 'AM':
            return True
        else:
            return False


    def init_keyword(self, keyword, value, timeout=None):
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
                apflog("Cannot write %s: %s %s" % (keyword.full_name, type(e), e),\
                       level='warn',echo=True)
                trials += 1
            else:
                success = True
        return success

    def init_gexptime(self):

        timeout = None
        if self.gexptime.binary > 1:
            timeout = self.gexptime.binary

        ret_val = self.init_keyword(self.gexptime, 1, timeout=timeout)
        return ret_val

    def init_sumframe(self):
        ret_val = self.init_keyword(self.sumframe, 1)
        return ret_val


    def initGuideCam(self):

        if self.is_gcam_power is False:
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

    def validateUCAMoutputs(self):

        if self.outfile.read() != self.desired_outfile:
            apflog("Output filename is %s and not the current date %s" % (self.outfile, self.desired_outfile),level='error',echo=True)
            self.outfile.write(self.desired_outfile)

        if self.obsnum < self.robot["MASTER_LAST_OBS_UCSC"]:
            apflog("Output file number is %s which is less than the last logged value %s" % (self.obsnum, self.robot["MASTER_LAST_OBS_UCSC"]),level='error',echo=True)
            self.obsnum.write(self.robot["MASTER_LAST_OBS_UCSC"])

        return

    def avgTelTemps(self):
        """
        self.avgTelTemp()
        sets the APFControl.avgtemps attribute.
        This will be a numpy array with the values in the order of:

        * M1 temp (eostele.TM1S210)
        * M2 temp (eostele.TM2CSUR)
        * Telescope average truss temp (eostele.TAVERAGE)
        * M2 air temp (eostele.TM2CAIR)
        * Air temp at FCU2 (eostele.TEMPNOW3)
        * Air temp at FCU3 (eostele.TEMPNOW4)

        The values for the trusses themselves are also included.
        This checks the values and if one is out of nominal range,
        after taking into the offset between the truss temperature
        and the usual average, it computes a new value for the
        truss average using only the trusses with nominal temperatures.

        """

        #self.avgtemps = [self.avg_lists[nm] for nm in ('TM1S210','TM2CSUR','TAVERAGE','TM2CAIR','TEMPNOW3','TEMPNOW4')]
        self.avgtemps = []
        self.avgtemps.append(self.m1tempkw.read(binary=True))
        self.avgtemps.append(self.m2tempkw.read(binary=True))
        self.avgtemps.append(self.taveragekw.read(binary=True))
        self.avgtemps.append(self.m2airkw.read(binary=True))
        self.avgtemps.append(self.temp3now.read(binary=True))
        self.avgtemps.append(self.temp4now.read(binary=True))

        # TAVERAGE is the average of the four trusses
        # temp_names = ["TTRUS135","TTRUS225","TTRUS045","TTRUS315"]
        # offsets = dict()
        # offsets['TTRUS135'] = -0.18
        # offsets['TTRUS225'] =  0.33
        # offsets['TTRUS045'] = -0.32
        # offsets['TTRUS315'] =  0.15

        # replaceavg = 0.0
        # n_good = 0
        # for tnm in temp_names:
        #     curtemp = np.average(self.mon_lists[tnm])
        #     if np.abs(curtemp - self.avgtemps[2] - offsets[tnm]) < 1:
        #         #bad
        #         replaceavg += curtemp
        #         n_good += 1
        # if n_good < 4 and n_good > 0:
        #     replaceavg /= n_good

        self.avgtemps = np.asarray(self.avgtemps)
        return

    def predTelFocus(self):

        self.avgTelTemps()
        # m1 m2 tavg m2air tf3 tf4

        # values as of July 3 2022
        # Final_focus_temp_fits.ipynb which (along with the data) should get checked in
        slopes = np.asarray([-0.007930, 0.021897, 0.011436, -0.008821, 0.001830, -0.017760])
        zeropoint_temps = np.asarray([15.549, 14.195, 14.610, 13.285, 15.716, 15.797])

        predfoc = np.sum(slopes*(self.avgtemps-zeropoint_temps)) + TELFOCUSTYP
        predfoc /= 1000.0 # convert to meters
        return predfoc


    # Function for checking what is currently open on the telescope
    def isOpen(self):
        """Returns the state of checkapf.WHATSOPN as a tuple (bool, str)."""
        try:
            whatstr = str(self.whatsopn)
            what = whatstr.split()
        except:
            apflog("checkapf.WHATSOPN returned a value that str.split cannot split",level='warn',echo=True)
            return False, ''
        if hasattr(what,'__iter__'):
            if "DomeShutter" in what or "MirrorCover" in what or "Vents" in what:
                return True, what
            else:
                return False, ''
        else:
            return False, ''

    # Fucntion for checking what is currently open on the telescope

    def isReadyForObservingDirect(self):
        what = ''
        rv = False
        try:
            ismcopen = self.mcopen.read(binary=True,timeout=2)
        except:
            return False, ''
        try:
            isshutterclosed = self.shclosed.read(binary=True,timeout=2)
        except:
            return False, ''

        if ismcopen:
            what = what + "MirrorCover"
            rv = True
        if isshutterclosed == False:
            if len(what) > 0 :
                what = what + " DomeShutter"
            else:
                what = "DomeShutter"
            if rv:
                rv = True
        else:
            rv = False
        return rv, what

    def isReadyForObserving(self):
        """Returns the state of checkapf.WHATSOPN as a tuple (bool, str)."""
        try:
            whatstr = str(self.whatsopn)
            what = whatstr.split()
        except:
            apflog("checkapf.WHATSOPN returned a value that str.split cannot split",level='warn',echo=True)
            return self.isReadyForObservingDirect()


        if hasattr(what,'__iter__'):
            if "DomeShutter" in what and "MirrorCover" in what:
                return True, what
            else:
                return False, ''
        else:
            return self.isReadyForObservingDirect()

    def setObserverInfo(self, num=10000, name='Robot', owner='public'):
        if self.test: return
        apflog("Setting science camera parameters.")
        self.ucam('OBSERVER').write(name)
        self.apfschedule('OWNRHINT').write(owner)
        self.outfile.write(name)
        self.desired_outfile = name
        self.ucam('OUTDIR').write('/data/apf/')
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

    def instrPermit(self):
        try:
            while not self.instr_perm.read(binary=True,timeout=2) or self.userkind.read(binary=True,timeout=2) != 3:
                apflog("Waiting for instrument permission to be true and userkind to be robotic")
                APFTask.waitfor(self.task, True, expression="$checkapf.INSTR_PERM = true", timeout=600)
                APFTask.waitfor(self.task, True, expression="$checkapf.USERKIND = robotic", timeout=600)
        except Exception as e:
            apflog("Cannot communicate with checkapf: %s" % (e),echo=True,level='alert')
            return False
        else:
            return True

    def turnOffLamps(self):
        for lamp in ("HALOGEN2","HALOGEN1","THORIUM1","THORIUM2"):
            try:
                rv = ktl.write("apfmot",lamp,"Off",wait=False)
            except Exception as e:
                apflog("Exception: %s" % (e),echo=True,level="alert")
                rv = False
            if rv is False:
                apflog("Cannot turn off lamp %s" % (lamp),echo=True,level="alert")
        return rv


    def writeStages(self,stagelist,component,state):
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
        rv1 = self.writeStages(stagelist,'MOE','Off')
        rv2 = self.writeStages(stagelist,'MOD','Pos')

        stagelist = ['ADC','DECKER','DEWARFOC']
        rv3 = self.writeStages(stagelist,'MOE','On')
        rv4 = self.writeStages(stagelist,'MOD','Pos')

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
        rv1 = self.writeStages(stagelist,'MOE','Off')
        rv2 = self.writeStages(stagelist,'MOD','Pos')

        stagelist = ['DECKER','DEWARFOC']
        rv3 = self.writeStages(stagelist,'MOE','On')
        rv4 = self.writeStages(stagelist,'MOD','Pos')

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
        rv1 = self.writeStages(stagelist,'MOE','Off')
        rv2 = self.writeStages(stagelist,'MOO','Off')
        rv3 = self.writeStages(['DEWARFOC'],'MOE','On')
        return rv1 and rv2 and rv3

    def turn_off_inst(self):
        '''
        turn_off_inst(self)
        Turns off all motor stages.
        '''
        stagelist = ['ADC','GUIDEFOC','CALMIRROR','CALSOURCE','IODINE','DECKER','DEWARFOC']
        rv1 = self.writeStages(stagelist,'MOE','Off')
        rv2 = self.writeStages(stagelist,'MOO','Off')

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

    def focusinstr(self, obsnum=None):
        self.instrPermit()
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

        self.validateUCAMoutputs()

        lastfocus_dict = APFTask.get("focusinstr", ["lastfocus","nominal"])
        if float(lastfocus_dict["lastfocus"]) > DEWARMAX or float(lastfocus_dict["lastfocus"]) < DEWARMIN:
            lastfocus_dict["lastfocus"] =  lastfocus_dict["nominal"]

        result, msg = self.runFocusinstr()
        if result:
            apflog(msg,echo=True)
        else:
            if not self.instr_perm.read(binary=True):
                msg += "Focusinstr has failed because of loss of instrument permission, waiting for return"
                wait = True
            else:
                wait = False
            # this complication is just to ensure the logging happens before waiting for permissions
            apflog(msg,echo=True,level='Alert')
            if wait:
                self.instrPermit()

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

        self.validateUCAMoutputs()

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
        result, code = apftaskDo(cmd,debug=True,cwd=os.getcwd())
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

    def runFocusinstr(self,flags="-b"):
        """Runs the focus routine appropriate for the user."""

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
        result, code = apftaskDo(cmd,cwd=os.getcwd())

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


    def findStar(self):
        ra = self.tel['RA'].read()
        dec = self.tel['DEC'].read()
        rah,ram,ras = ra.split(":")
        decd,decm,decs = dec.split(":")

        cmd = os.path.join(SCRIPTDIR,"closest")
        cmdargs =  [cmd, rah,ram, ras, decd,decm,decs, "5","1","8"]
        sfncat = os.path.join(LROOT,"data/apf/StarCatalog.dat")
        try:
            starcat = open(sfncat)
        except:
            apflog("Cannot open file %s" % (sfncat), level="warn",echo=True)
            return False
        #$line[3] $line[4] $line[5] $line[6] $line[7] 5 1 "8"] < /usr/local/lick/data/apf/StarCatalog.dat"
        p = subprocess.Popen(cmdargs, stdin=starcat, stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=os.path.curdir)
        out, err = p.communicate()
        ret_code = p.returncode
        if ret_code != 0:
            apflog(out,echo=True)
            apflog(err, level="warn",echo=True)
            return False

        return out.split()

    def slew(self,star):

        cmd = os.path.join(SCRIPTDIR,'slewlock')
        try:
            ra = float(star[1])
            ra *= 3.819718
            dec = float(star[2])
            dec *= 57.295779
        except:
            return False
        cmd +=  ' %s %s %f %f %s %s %d ' % ("reference",star[0],ra, dec,star[4],star[5],210)
        if self.test:
            apflog("Would slew by executing %s" %(cmd), echo=True)
            return True

        apflog("Slewing by executing %s" %(cmd), echo=True)
        result, code = apftaskDo(cmd,cwd=os.path.curdir,debug=True)
        if not result:
            apflog("Failed at slewlock - check logs could be a slew failure or a failure to acquire: %s" %(code), level="error", echo=True)

        return result

    def saveMovie(self):
        now = datetime.now()
        self.fits3pre.write('%d%02d%02d_%s_' % (now.year,now.month,now.day, self.tel['TARGNAME'].read()))
        self.fits3dir.write('/data/apfguide')
        self.save3d.write(True)
        return

    def stopMovie(self):
        self.save3d.write(False)
        self.fits3dir.write('/tmp/')
        return

    def runFocusTel(self):
        """Runs the telescope focus routine."""
        el = self.tel['EL'].read(binary=True)
        cfspos = self.fspos.read(binary=True)
        crspos = self.rspos.read(binary=True)

        if abs(el - cfspos) < 2.5 or abs(el - crspos) < 2.5:
            apflog("Cannot focus, telescope too close to shutter", level="warn", echo=True)
            return False

        if self.test:
            APFTask.waitFor(self.task, True, timeout=10)
            apflog("Test Mode: Would be running focus_telescope.",echo=True)
            return True

        self.saveMovie()

        apflog("Running focus_telescope routine.",echo=True)
        cmd = os.path.join(SCRIPTDIR,'focus_telescope -c %.3f' % (float(self.predTelFocus())*1000.0))
        result, code = apftaskDo(cmd,cwd=os.path.curdir)

        self.stopMovie()

        try:
            self.guide['MODE'].write('Guide')
        except:
            apflog('Cannot modify apfguide.MODE to Guide.',level='error',echo=True)

        if not result:
            apflog("focustel failed with code %d" % code, echo=True)
            expression="($apftask.FOCUSINSTR_STATUS != 0) and ($apftask.FOCUSINSTR_STATUS != 1) "
            if not APFTask.waitFor(self.task,True,expression=expression,timeout=30):
                apflog("focus_telescope failed to exit" ,echo=True)
            return result

        return True

    def runAutoexposure(self,ind=5):

        cmd = os.path.join(SCRIPTDIR,'autoexposure')
        istr = "%d" % (ind)
        cmdargs = cmd

        result, code = apftaskDo(cmdargs,cwd=os.path.curdir)

        if not result:
            apflog("autoexposure failed with code %d" % code, echo=True)
        return result

    def runCenterup(self):

        cmd = os.path.join(SCRIPTDIR,'centerup')
        result, code = apftaskDo(cmd,cwd=os.path.curdir)
        if not result:
            apflog("centerup failed with code %d" % code, echo=True)
        return result

    def findStarfocusTel(self):
        """ This finds a star in the catalog of pointing reference stars close
        to the current position of the telescope.
        It then runs slewlock to slew to the star, set up the guider, and center
        the star on the slit.
        Finally, once the star is acquired, it runs the telescope focus routine.
        """
        star = self.findStar()
        if not star:
            apflog("Cannot find star near current position!?",level='error',echo=True)
            return False
        apflog("Targeting telescope on %s" % star[0], echo=True)
        try:
            self.vmag.write(star[6])
        except Exception as e:
            apflog("Cannot write SCRIPTOBS_VMAG: %s" % (e), level='error',echo=True)
        try:
            sline = "%s %f %f pmra=%s pmdec=%s vmag=%s # end" % (star[0],float(star[1])*3.819718,float(star[2])*57.295779,star[4],star[5],star[6])
            self.line.write(sline)
        except Exception as e:
            apflog("Cannot write SCRIPTOBS_LINE: %s" % (e), level='error',echo=True)

        try:
            self.robot['SCRIPTOBS_LINE_RESULT'].write(0)
            self.robot['SCRIPTOBS_OBSERVED'].write(False)
        except Exception as e:
            apflog("Cannot write 0 to SCRIPTOBS_LINE_RESULT or False to SCRIPTOBS_OBSERVED: %s" % (e), level='warn', echo=True)


        predfocus  = self.predTelFocus()
        self.robot['FOCUSTEL_STARTFOCUS'].write(predfocus)
        self.robot['FOCUSTEL_LASTFOCUS'].write(predfocus)
        self.focus.write(self.predTelFocus(), binary=True, wait=False)
        if self.slew(star):
            return self.runFocusTel()
        return False


    def setTeqMode(self, mode):
        apflog("Setting TEQMode to %s" % mode)
        if self.test:
            apflog("Would be setting TEQMode to %s" % mode)
            return
        self.teqmode.write(mode,wait=False)
        result = self.teqmode.waitfor('== %s' % mode, timeout=60)
        if not result:
            apflog("Error setting the TEQMODE.")
            raise RuntimeError("Couldn't set TEQ mode")


    def clearestop(self):
        if self.test: return True
        if self.mv_perm.binary == False:
            apflog("Waiting for permission to move...", echo=True)
            chk_move = "$checkapf.MOVE_PERM == true"
            result = APFTask.waitFor(self.task, False, chk_move, timeout=600)
            if not result:
                apflog("Can't open. No move permission.",echo=True)
                return False

        cmd = os.path.join(SCRIPTDIR,'clear_estop')
        result, code = apftaskDo(cmd,debug=True,cwd=os.getcwd())
        if result:
            try:
                estopstate = self.dome.read('ESTOPST',binary=True)
                if estopstate:
                    return False
                else:
                    return True
            except:
                return False
        else:
            return False

    def statesSet(self):
        # there are three states - but we do not care about ESTOPST, that is will be cleared in openatsunset/openatnight
        if self.dome['ECLOSEST']:
            return True
        if self.dome['ESECURST']:
            return True
        return False


    def homeTelescope(self):
        cmd = os.path.join(SCRIPTDIR,"slew") + " --home"
        rv, rc = apftaskDo(cmd)
        try:
            homed = self.apfmon('ELHOMERIGHTSTA').read(binary=True,timeout=2)
        except Exception as e:
            apflog("cannot read apfmon keyword ELHOMERIGHTSTA: %s" % (e),level='Alert',echo=True)
            return False
        else:
            if rc == 0 and homed == 2:
                return True
            else:
                apflog("cannot home telescope",level='Alert',echo=True)
                return False

    def checkHome(self,home=True):
        try:
            homed = self.apfmon('ELHOMERIGHTSTA').read(binary=True,timeout=2)
        except Exception as e:
            apflog("apfmon.ELHOMERIGHTSTA cannot be read: %s" % (e),level='Alert',echo=True)
            return False
        if homed == 2:
            return True
        else:
            if homed == 5 or homed == 6:
                if home:
                    self.homeTelescope()
                else:
                    apflog("Telescope needs to be homed",level='Alert',echo=True)
                    return False
            else:
                apflog("apfmon.ELHOMERIGHTSTA value is %d" % (homed),level='Alert',echo=True)
                return False

        return False

    def openat(self, sunset=False):
        """Function to ready the APF for observing. Calls either openatsunset or openatnight.
           This function will attempt to open successfully twice. If both attempts
           fail, then it will return false, allowing the master to register the error
           and behave accodingly. Otherwise it will return True. """
        # If this is a test run, just return True
        if self.test: return True

        if not self.ok2open:
            # This should really never happen. In case of a temporary condition, we give
            # a short waitfor rather than immediatly exiting.
            chk_open = "$checkapf.OPEN_OK == true"
            result = APFTask.waitFor(self.task, False, chk_open, timeout=30)
            if not result:
                apflog("Tried calling openat with OPEN_OK = False. Can't open.", echo=True)
                apflog(self.checkapf["OPREASON"].read(), echo=True)
                return False

        if float(self.sunel) > SUNEL_HOR:
            apflog("Sun is still up. Current sunel = %4.2f. Can't open." % self.sunel, echo=True)
            return False

        if self.mv_perm.binary == False:
            apflog("Waiting for permission to move...", echo=True)
            chk_move = "$checkapf.MOVE_PERM == true"
            result = APFTask.waitFor(self.task, False, chk_move, timeout=600)
            if not result:
                apflog("Can't open. No move permission.",echo=True)
                return False

        if self.statesSet():
            apflog("An unusal emergency state is set.", level="error",echo=True)
            return False

        # Everything seems acceptable, so lets try opening
        if sunset:
            cmd = os.path.join(SCRIPTDIR,'openatsunset')
        else:
            cmd = os.path.join(SCRIPTDIR,'openatnight')

        # Make two tries at opening. If they both fail return False so the caller can act
        # accordingly.
        result, code = apftaskDo(cmd)
        if not result:
            apflog("First openup attempt has failed. Exit code = %d. After a pause, will make one more attempt." % code,echo=True)
            APFTask.waitFor(self.task, True, timeout=10)
            result, code = apftaskDo(cmd)
            if not result:
                apflog("Second openup attempt also failed. Exit code %d. Giving up." % code,echo=True)
                return False
        rv = self.checkHome()
        if rv == False:
            return False
        try:
            APFLib.write("eostele.FOCUS",ktl.read("apftask","FOCUSTEL_LASTFOCUS",binary=True))
        except:
            apflog("Cannot move secondary focus.",level="error")
            return False
        return True

    def powerDownTelescope(self):
        """Checks that we have the proper permission and dome is closed, then resets telescope power."""
        if self.test: return True
        cmd = os.path.join(SCRIPTDIR,"power_down_telescope")
        self.DMReset()
        if self.mv_perm.binary == False:
                apflog("Waiting for permission to move")
        chk_mv = '$checkapf.MOVE_PERM == true'
        result = APFTask.waitFor(self.task, False, chk_mv, timeout=300)
        if not result:
            apflog("Didn't have move permission after 5 minutes.", echo=True)
            return False
        # one last check

        apflog("Running power_down_telescope script")
        result, code = apftaskDo(cmd)
        if result:
            return True
        else:
            apflog("power_down_telescope has failed. Human intervention likely required.", level='alert', echo=True)
            return False


    def servoFailure(self):
        """checks for amplifier faults
        Checks all seven moving components.
        """
        servo_failed = False
        prefixs = ["AZ","EL","FA","FB","FC","TR" ]
        for pr in prefixs:
            nm = pr + "AMPFLT"
            val = self.tel[nm].read(binary=True)
            if val:
                servo_failed = True
                apflog("Error: Servo Amplifier Fault: %s %s" % (nm,val), level="alert", echo=True)

        for pr in prefixs:
            nm = pr + "FERROR"
            val = self.tel[nm].read(binary=True)
            if val:
                servo_failed = True
                apflog("Error: Fatal Following Error: %s %s" % (nm,val), level="alert", echo=True)

        return servo_failed

    def close(self, force=False):
        """Checks that we have the proper permission, then runs the closeup script.
        On failures retries, and also sends alerts. Good for waking people up.
        """

        if self.test: return True
        cmd = os.path.join(SCRIPTDIR,"closeup")
        if force:
            apflog("Calling a single instance of closeup. Will return regardless of result.", echo=True)
            result, code = apftaskDo(cmd)
            return result
        if self.mv_perm.binary == False:
            if self.chk_close.binary == True:
                apflog("Waiting for checkapf to close up")
            else:
                apflog("Waiting for permission to move")
        chk_mv = '$checkapf.MOVE_PERM == true'
        result = APFTask.waitFor(self.task, False, chk_mv, timeout=1200)
        if not result:
            apflog("Didn't have move permission after 20 minutes. Going ahead with closeup.", echo=True)
            return False
        try:
            if self.apfmon['FRONT_SHUTTER_CLOSEUPSTA'].read(binary=True,timeout=2) != 2:
                # this is a check to see if the front shutter got caught running
                # away, if so do not send any more shutter commands
                apflog("Dome Shutters maybe running away!", level='error', echo=True)
                return False
        except:
            apflog("Cannot communicate with apfmon1, proceeding anyway (fingers crossed)", level='error', echo=True)


        apflog("Running closeup script")
        attempts = 0
        close_start = datetime.now()
        while (datetime.now() - close_start).seconds < 1800:
            result = APFTask.waitFor(self.task, False, chk_mv, timeout=300)
            if not result:
                apflog("Didn't have move permission after 5 minutes. ", echo=True)
                break
            attempts += 1
            result, code = apftaskDo(cmd)
            if not result:
                apflog("Closeup failed with exit code %d" % code, echo=True)
                if self.servoFailure() or self.slew_allowed.read(binary=True) is False:
                    apflog("Likely Servo amplifier failure, may power cycle telescope",echo=True,level='alert')
                    rv = self.powerDownTelescope()
                    if rv:
                        apflog("Power cycled telescope",echo=True,level="error")
                    else:
                        apflog("Failure power cycling telescope",echo=True,level="alert")
                if attempts == 7:
                    lstr = "Closeup has failed %d times consecutively. Human intervention likely required." % (attempts)
                    areopen, whatsopen = self.isOpen()
                    if areopen == True:
                        # truly dire, the telescope is open
                        apflog(lstr, level='Alert', echo=True)
                    else:
                        # telescope powered on, and possibly in the wrong place, but not open
                        apflog(lstr, level='error', echo=True)
                APFTask.waitFor(self.task, True, timeout=30)
            else:
                break
        if result:
            try:
                APFTask.set(self.task, suffix='LAST_CLOSE', value=time.time())
            except:
                apflog("cannot write apftask.MASTER_LAST_CLOSE",level='warn',echo=True)
            return True
        else:
            apflog("Closeup could not successfully complete.")
            return False
        return False

    def updateLastObs(self,obsnum):
        """ If the last observation was a success,
        this function updates the file storing the last
        observation number and the hit_list which is
        required by the dynamic scheduler."""

        if obsnum['populated']:
            if obsnum >= 10000:
                APFLib.write(self.robot["MASTER_LAST_OBS_UCSC"], obsnum)

        return


    def setTelFoc(self):
        """
        Sets the telescope focus to the predicted value returned by
        predTelFocus()
        """

        predfocus  = self.predTelFocus()
        self.robot['FOCUSTEL_STARTFOCUS'].write(predfocus)

        if self.mv_perm and self.faenable['binary'] == 1:
            try:
                self.focus.write(predfocus,binary=True,wait=False)
                self.robot['MASTER_MESSAGE'].write("Wrote %f to eostele.Focus" % (predfocus*1000.) )
            except Exception as e:
                apflog("Cannot write eostele.FOCUS: %s" % (e), level="error", echo=True)

    def setAutofocVal(self):
        """ APFControl.setAutofocVal()
            tests when the last time the telescope was focused, if more than FOCUSTIME enable focus check
        """

        # check last telescope focus
        lastfoc = self.robot['FOCUSTEL_LAST_SUCCESS'].read(binary=True)

        if self.sunRising() and (self.sunel.read(binary=True) > -20):
            self.autofoc.write("robot_autofocus_disable")
            return 0
        
        if time.time() - lastfoc < 3600:
            self.autofoc.write("robot_autofocus_disable")
            return 0

        predfocus  = self.predTelFocus()
        self.robot['FOCUSTEL_STARTFOCUS'].write(predfocus)
        focus_diff = math.fabs(predfocus - self.focus['binary'])

        #self.focus.write(predfocus,binary=True,wait=False)
        self.autofoc.write("robot_autofocus_enable")
        focval = 1
        APFTask.set(self.task, suffix="MESSAGE", value="Current telescope focus more than %6.3f microns from predicted, setting to %.3f." % (focus_diff*1e6,predfocus*1000), wait=False)

        return focval

    def updateWindshield(self, state):
        """Checks the current windshielding mode.
        If the input state is auto, makes sure the mode is set properly based on wind speed and temperature.
        Otherwise, the input state defines the mode.
        """
        currMode = self.robot["SCRIPTOBS_WINDSHIELD"].read().strip().lower()
        rv = currMode
        if state == 'on':
            if currMode != 'enable':
                apflog("Setting scriptobs_windshield to Enable")
                rv = "Enable"
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], rv)
        elif state == 'off':
            if currMode != 'disable':
                apflog("Setting scriptobs_windshield to Disable")
                rv = "Disable"
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], rv)

        else:
            # State must be auto, so check wind and temperature.
            # This state enables or disables windshielding based on the wind speed and the outside temperature
            if self.down > 0:
                wvel = self.avg_lists['M3WIND']
            else:
                wvel = self.avg_lists['M5WIND']

            apflog("Current median wind speed is %.2f with the limit %.2f" % (wvel,WINDSHIELD_LIMIT), level='debug')
            if currMode == 'enable' and wvel <= WINDSHIELD_LIMIT and float(self.avg_lists['M5OUTEMP']) > TEMP_LIMIT:
                apflog("Setting scriptobs_windshield to Disable")
                rv = "Disable"
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], rv)

            if currMode == 'disable' and (wvel > WINDSHIELD_LIMIT or float(self.avg_lists['M5OUTEMP']) < TEMP_LIMIT):
                apflog("Setting scriptobs_windshield to Enable")
                rv = "Enable"
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], rv)

        return rv


    def checkFCs(self):
        """
        If the dome is open the FCs should be off. Sometimes they are commanded off but do not turn off.
        This will attempt to fix that by checking if the dome is open, checking if the FCs are on, and
        then turning them on, then off. Which is the only way to turn them off.
        """

        for fc in ('FC2','FC3'):
            if 'Vents' in self.whatsopn['ascii'] or 'DomeShutter' in self.whatsopn['ascii']:
                # fcs should be off
                if self.dome[fc].read(binary=True):
                    try:
                        self.dome[fc + 'CMD'].write(True)
                        time.sleep(.1)
                        self.dome[fc + 'CMD'].write(False)
                    except:
                        pass
            else:
                # fcs should be on
                if self.dome[fc].read(binary=False):
                    try:
                        self.dome[fc + 'CMD'].write(False)
                        time.sleep(.1)
                        self.dome[fc + 'CMD'].write(True)
                    except:
                        pass

        return

    def eveningStar(self):
        """Aim the APF at the desired target. This calls prep-obs, slewlock, and focus-telescope."""
        if self.isOpen()[0] == False:
            apflog("APF is not open. Can't target a star while closed.",level='error',echo=True)
            return
        self.DMReset()

        if self.calsta['binary'] < 3 or self.focussta['binary'] < 3:
            apflog('Focusinstr and/or Calibrate are running, will skip evening star observation. focusinstr=%s calibrate=%s' % (self.calsta,self.focussta),echo=True)
            return

        # check on weirdness for UCAM host post-reboot
        self.ucamdispatchmon()

        # Call prep-obs
        apflog("Calling prep-obs.",echo=True)
        prepobs = os.path.join(SCRIPTDIR,'prep-obs') + ' --evening'
        result, ret_code = apftaskDo(prepobs)
        if result == False:
            # try again
            self.DMReset()
            result, ret_code = apftaskDo(prepobs)
            if result is False:
                apflog("Prep-obs returned error code %d. Targeting object has failed." % (ret_code),level='error',echo=True)
                return

        self.decker.write('W',wait=False)

        self.DMReset()
        apflog("Slewing to lower el",echo=True)
        result, ret_code = apftaskDo('slew -e 75')
        if result == False:
            apflog("Slew returned error code %d. Targeting object has failed." % (ret_code),level='error',echo=True)
            return
        # Slew to the specified RA and DEC, set guide camera settings, and centerup( Slewlock )
        # Focus the telescope - all of this, including finding the star, is done in focusTel
        self.DMReset()
        if self.findStarfocusTel():
            try:
                self.robot['SCRIPTOBS_LINE_RESULT'].write(3)
                self.robot['SCRIPTOBS_OBSERVED'].write(True)
                self.guide['CLEARSUMS'].write('now')
                self.guide['CLEARSUMS'].write('gstate')
            except Exception as e:
                apflog("Cannot write 3 to SCRIPTOBS_LINE_RESULT or True to SCRIPTOBS_OBSERVED: %s" % (e), level='warn', echo=True)
            return True
        else:
            try:
                self.robot['SCRIPTOBS_LINE_RESULT'].write(2)
            except Exception as e:
                apflog("Cannot write 2 to SCRIPTOBS_LINE_RESULT: %s" % (e), level='warn', echo=True)
            return False

    def DMReset(self):
        try:
            APFLib.write(self.checkapf['ROBOSTATE'], "master operating",timeout=10)
        except Exception as e:
            try:
                ukind = self.userkind.read()
            except:
                ukind = "Unknown"
            ostr = "Error: Cannot write to ROBOSTATE, USERKIND = %s, reason: %s" % (ukind,e)
            apflog(ostr,level='error',echo=True)

    def DMZero(self):
        try:
            if self.checkapf['DMTIME'].read(binary=True) < 1:
                APFLib.write(self.dmtimer, -1,timeout=10)
        except Exception as e:
            ostr = "Error: cannot touch DM Timer: %s " %( e)
            apflog(ostr,level='error',echo=True)

    def findRobot(self):
        """Trys to find a running instance of scriptobs.
            Returns the PID along with a boolean representing
            if the robot was succesfully found."""
        rpid = self.robot['SCRIPTOBS_PID'].read(binary=True)
        if rpid == '' or rpid == -1:
            return rpid, False
        else:
            return rpid, True

    def startRobot(self,observation=None,skip=False,raster=False):
        """Start an instance of scriptobs. Returns the result from subprocess.Popen()."""
        # For running in test mode
        if self.test:
            apflog("Would start robot",echo=True)
            if observation is not None:
                apflog("Would be taking observation in starlist %s" % observation,echo=True)
            APFTask.waitFor(self.task, True, timeout=10)
            return

        # Make sure the telescope autofocus is enabled
        APFLib.write(self.autofoc, "robot_autofocus_enable")
        chk_foc = '$apftask.SCRIPTOBS_AUTOFOC == robot_autofocus_enable'
        result = APFTask.waitFor(self.task, False, chk_foc, timeout=60)
        if not result:
            apflog("Error setting scriptobs_autofoc", level='error',echo=True)
            return

        # Make sure APFTEQ is in night mode for observations
        if self.teqmode.read() != 'Night':
            self.setTeqMode('Night')

            # Check the instrument focus for a reasonable value
        if self.dewarfoc > DEWARMAX or self.dewarfoc < DEWARMIN:
            lastfit_dewarfoc = ktl.read("apftask","FOCUSINSTR_LASTFOCUS",binary=True)
            apflog("Warning: The dewar focus is currently %d. This is outside the typical range of acceptable values. Resetting to last derived value %d" % (self.dewarfoc,lastfit_dewarfoc), level = "error", echo=True)
            APFLib.write("apfmot.DEWARFOCRAW",lastfit_dewarfoc)

        # check on weirdness for UCAM host post-reboot
        self.ucamdispatchmon()



        telstate = self.tel['TELSTATE'].read()
        if telstate == 'Disabled':
            rv, retc = apftaskDo(os.path.join(SCRIPTDIR,"slew --hold"))
            if not rv:
                return rv
        rv, retc = apftaskDo(os.path.join(SCRIPTDIR,"prep-obs"))
        if not rv:
            return rv
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


    def killRobot(self, now=False):
        """ In case during an exposure there is a need to stop the robot and close up."""
        apflog("Terminating scriptobs")
        if now:
            apflog("Abort exposure, terminating scriptobs now.")
        else:
            if not self.ucam['EVENT_STR'].read() == "ControllerReady":
                apflog("Waiting for current exposure to finish.")
                self.ucam['EVENT_STR'].waitfor(" = ReadoutBegin", timeout=1200)

        ripd, running = self.findRobot()
        if running:
            apflog("Killing scriptobs %s" % (str(ripd)))
            try:
                APFLib.write(self.robot['SCRIPTOBS_CONTROL'], "abort")
            except Exception as e:
                errstr = "Cannot abort scriptobs: %s" % (e)
                apflog(errstr,level="Warn",echo=True)


    def testBias(self):

        apftask = ktl.Service('apftask')

        if self.ucampower == False:
            self.ucampower.write('On',wait=False)
            rv = self.ucampower.waitFor('== On',timeout=30)
            APFTask.wait(self.task, True, timeout=5)
            if rv == False:
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

        apfschedule = ktl.Service('apfschedule')
        # check if the focusinstr or calibrate tasks are already running
        if ktl.read('apftask','FOCUSINSTR_PID',binary=True) > 0 or ktl.read('apftask','CALIBRATE_PID',binary=True) > 0 or ktl.read('apftask','SCRIPTOBS_PID',binary=True) > 0 :
            return None

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



    def ucamPowercycle(self, fake=False):

        if fake:
            apflog("would have executed %s" % (os.path.join(SCRIPTDIR,"robot_power_cycle_ucam")))
            return True
        else:
            val = subprocess.call(os.path.join(SCRIPTDIR,"robot_power_cycle_ucam"))
            if val > 0:
                apflog("power cycle of UCAM failed",level='alert')
                return False

            return True

        return True



    def ucamReboot(self,fake=False):
        if fake:
            apflog("Would have rebooted UCAM host ",echo=True)
            return True

        apftask = ktl.Service('apftask')
        command = apftask['UCAMLAUNCHER_UCAM_COMMAND']
        ucamstat = apftask['UCAMLAUNCHER_UCAM_STATUS']
        status = apftask['UCAMLAUNCHER_STATUS']

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
        status = apftask['UCAMLAUNCHER_STATUS']

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
        else:
            apflog("UCAM host reboot failure, combo_ps still not ok" , level="alert", echo=True)

        self.obsnum.monitor()
        self.obsnum.callback(self.updateLastObs)

        self.event.monitor()
        self.event.callback(self.eventmon)


    def ucamRestart(self, fake=False):

        # modify -s apftask UCAMLAUNCHER_UCAM_COMMAND=Stop
        if fake:
            # would have restarted software
            apflog("Would have restarted UCAM software ")
            return True
        else:
            try:
                apflog("Stop and restarting UCAM software",echo=True)
                ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","stop")
                if self.combo_ps.waitFor(" == MissingProcesses",timeout=30):
                    ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","run")
                    nv = self.combo_ps.waitFor(" == Ok",timeout=30)
                    if nv:
                        return nv
                    else:
                        apflog("UCAM  restart failure, combo_ps still not ok" , level="error", echo=True)
            except:
                apflog("UCAM status bad, cannot restart",level='alert')
                return False

            self.ucamReboot()


        return False



    def ucamStatus(self, fake=False):

        if self.ctalk.read(binary=True) > 0:
            rv = self.ucamPowercycle(fake=fake)
            return rv

        if self.combo_ps.read(binary=True) > 0:
            # brains!
            rv = self.ucamRestart(fake=fake)
            return rv

        try:
            ucamsta0 = self.ucam['DISP0STA'].read(binary=True)
            ucamsta1 = self.ucam['DISP1STA'].read(binary=True)
        except Exception as e:
            apflog('apfucam.DISPSTA failure, apfucam likely not running: %s' % (e),echo=True,level='Alert')
            rv = self.ucamRestart(fake=fake)
            return rv
        else:
            if ucamsta1 > 0 and ucamsta0 > 0:
                # Things are still starting up
                if ucamsta0 > 2:
                    # failure to connect
                    rv = self.ucamRestart(fake=fake)
                    return rv
                else:
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
        apf.robot['FOCUSTEL_STARTFOCUS'].write(apf.predTelFocus())
