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
except:
    from fake_apflog import *

from APFControl import apftask_do, restart, APF

WINDLIM = 40.0
SLOWLIM = 100
WINDSHIELD_LIMIT = 4.5 # mps at the APF
FOCUSTIME = 3600. # minimum time before checking telescope focus
TEMP_LIMIT = 35. # deg F at the APF
wxtimeout = datetime.timedelta(seconds=1800)
SUNEL_HOR = -3.2
TELFOCUSMIN = -0.00096
TELFOCUSMAX = -0.00060
# this value comes an average over many measurements of the telescope focus
#TELFOCUSTYP = -0.83529
#TELFOCUSTYP = -0.76529
TELFOCUSTYP = -0.69029
TELFOCUSMAXOFF = 0.00002

if "LROOT" in os.environ:
    LROOT = os.environ["LROOT"]
else:
    LROOT = '/usr/local/lick'
SCRIPTDIR = os.path.join(LROOT,'bin/robot/')

class TelescopeControl:
    '''
    TelescopeControl
    This class is used to control the telescope and its associated
    equipment.
    '''
    def __init__(self, apf=None, test=False):
        self.apf = apf
        self.test = test
        self.task = None
        if self.apf is not None:
            self.task = self.apf.task

        self.lastfocuscheck = datetime.datetime.now() - datetime.timedelta(days=1)
        self.lastfocusval = None

        self.dew_too_close = False
        self.openOK = False

        self.tel        = ktl.Service('eostele')
        self.sunel      = self.tel('SUNEL')
        self.ael        = self.tel('AEL')
        self.aaz        = self.tel('AAZ')
        self.aafocus    = self.tel('AAFOCUS')
        self.focus      = self.tel('FOCUS')
        self.faenable   = self.tel('FAENABLE')
        self.telstate   = self.tel('TELSTATE')

        self.dome       = ktl.Service('eosdome')
        self.rspos      = self.dome('RSCURPOS')
        self.fspos      = self.dome('FSCURPOS')
        self.shclosed   = self.dome('SHCLOSED')

        self.eostdio    = ktl.Service('eostdio')
        self.mcopen     = self.eostdio('MCOPEN')

        self.checkapf   = ktl.Service('checkapf')
        self.ok2open    = self.checkapf('OPEN_OK')
        self.userkind   = self.checkapf('USERKIND')
        self.dmtimer    = self.checkapf('DMTIME')
        self.whatsopn   = self.checkapf('WHATSOPN')
        self.mv_perm    = self.checkapf('MOVE_PERM')
        self.chk_close  = self.checkapf('CHK_CLOSE')

        self.apfmet     = ktl.Service('apftempest')
        self.wx         = self.apfmet('WINDAV')
        self.airtemp    = self.apfmet('TEMP')
        self.down       = self.apfmet('STATUS')

        self.eosti8k    = ktl.Service('eosti8k')
        self.m2tempkw   = self.eosti8k('TM2CSUR')
        self.m2airkw    = self.eosti8k('TM2CAIR')
        self.m1tempkw   = self.eosti8k('TM1S210')
        self.taveragekw = self.eosti8k('TAVERAGE')
        self.t045kw     = self.eosti8k('TTRUS045')
        self.t135kw     = self.eosti8k('TTRUS135')
        self.t225kw     = self.eosti8k('TTRUS225')
        self.t315kw     = self.eosti8k('TTRUS315')

        self.eosgcam    = ktl.Service('eosgcam')
        self.fits3pre   = self.eosgcam('FITS3PRE')
        self.fits3dir   = self.eosgcam('FITS3DIR')
        self.save3d     = self.eosgcam('SAVE3D')

        self.guide      = ktl.Service('apfguide')

        self.eoscool    = ktl.Service('eoscool')
        self.dewpt      = self.eoscool('DEWPAVG3')
        self.temp3now   = self.eoscool('TEMPNOW3')
        self.temp4now   = self.eoscool('TEMPNOW4')

        self.robot        = ktl.Service('apftask')
        self.slew_allowed = self.robot['SLEW_ALLOWED']
        self.autofoc      = self.robot["SCRIPTOBS_AUTOFOC"]
        self.apfteqsta    = self.robot['APFTEQ_STATUS']
        self.metxfersta   = self.robot['METSXFER_STATUS']
        self.slewsta      = self.robot['SLEW_STATUS']
        self.opensta      = self.robot['OPENUP_STATUS']
        self.closesta     = self.robot['CLOSEUP_STATUS']
        self.shuttersta   = self.robot['SHUTTERS_STATUS']
        self.focustelsta  = self.robot['FOCUSTEL_STATUS']
        self.lastopen     = self.robot['OPENUP_LAST_SUCCESS']

        self.apfteq     = ktl.Service('apfteq')
        self.teqmode    = self.apfteq['MODE']

        self.apfmon     = ktl.Service('apfmon')

        self.ok2open.monitor()
        self.ok2open.callback(self.ok_mon)

        self.dmtimer.monitor()

        self.teqmode.monitor()
        self.down.monitor()
        self.whatsopn.monitor()
        self.mv_perm.monitor()
        self.chk_close.monitor()
        self.slew_allowed.monitor()
        self.apfteqsta.monitor()
        self.metxfersta.monitor()

        self.sunel.monitor()
        self.aaz.monitor()
        self.ael.monitor()
        self.fspos.monitor()
        self.rspos.monitor()
        self.focus.monitor()
        self.aafocus.monitor()
        self.faenable.monitor()
        self.lastopen.monitor()

       # Initial Wind conditions
        self.wslist = []
        self.mon_lists = dict()
        self.avg_lists = dict()
        self.dewlist = []

        for kw in (self.m1tempkw,self.m2tempkw,self.m2airkw,self.taveragekw,\
                   self.t045kw,self.t135kw,self.t225kw,self.t315kw,self.temp3now,\
                    self.temp4now,self.wx,self.airtemp):
            self.mon_lists[kw['name']] = []
            self.avg_lists[kw['name']] = None
            kw.monitor()
            kw.callback(self.list_mon)
            kw.read()

        self.dewpt.monitor()
        self.dewpt.callback(self.dew_pt_mon)

        for kw in (self.slewsta,\
                   self.shuttersta, self.opensta, self.closesta,\
                    self.focustelsta):
            kw.monitor()

        # Grab some initial values for the state of the telescope

        self.wx.read()
        #self.altwx.read()
        self.dewpt.read()
        self.ok2open.read()
        self.avgtemps = np.asarray([self.avg_lists[nm] for nm in \
                                     ('TM1S210','TM2CSUR','TAVERAGE',\
                                      'TM2CAIR','TEMPNOW3','TEMPNOW4')])

    def __str__(self):
        # Determine if the sun rising / setting check is working
        now = datetime.datetime.now()
        self.avg_tel_temps()
        s = ''
        s += "At %s state of telescope is:\n" % str(now)
        s += "Sun elevation = %4.2f %s\n" % (self.sunel, "Rising" if self.sun_rising() else "Setting")
        s += "Telescope -- AZ=%4.2f  EL=%4.2f \n" % (self.aaz, self.ael)
        s += "Front/Rear Shutter=%4.2f / %4.2f\n"%(self.fspos, self.rspos)
        s += "Wind = %3.1f mps (APF)  \n" % (np.average(self.mon_lists['WINDAV']))
        s += "Last open time = %.2f sec\n" % (self.lastopen.binary)
        s += "Time since opening = %6.2f sec\n" % (time.time() - self.lastopen.binary)
        s += "M1 = %5.2f deg C M2 = %5.2f deg C Tel Avg = %5.2f deg C M2 Air = %5.2f deg C FCU3 = %5.2f deg C FCU4 = %5.2f deg C\n" % tuple(self.avgtemps)
        s += "Dewpt = %5.2f deg C Teq Mode - %s\n" % (np.average(self.dewlist),self.teqmode)
        s += "Too close to the dewpoint? = %s\n" % self.dew_too_close
        s += "M2 Focus Value = % 4.3f\n" % (float(self.aafocus['binary'])*1000.0)
        s += "M2 Focus Value = % 4.3f (focus kwd)\n" % (float(self.focus['binary'])*1000.0)
        s += "Preferred M2 Focus Value =  % 4.3f\n" % (float(self.pred_tel_focus())*1000.0)
        s += "Okay to open = %s -- %s\n" % (self.openOK, self.checkapf['OPREASON'].read() )
        s += "Current Weather = %s\n" % self.checkapf['WEATHER'].read()
        s += "Slew Allowed = %s\n" % self.slew_allowed['binary']
        s += "Is the system sane? = %s\n" % self.check_sanity()

        isopen, what = self.is_open()
        if isopen:
            s += "Currently open: %s\n" % what
        else:
            s += "Not currently open\n"

        focval = self.set_autofoc_val()
        s += "Focus value for scriptobs = %d\n" % focval
        windshield = self.update_windshield("auto")
        s += "Windshield state = %s\n" % windshield

        return s

    def list_mon(self,keyword):
        """
        list_mon(keyword)
        This is a callback function for the monitor keyword
        keyword. It is called when the keyword changes.
        It updates the list of values for the keyword and
        computes the average value for the keyword.
        """
        if keyword['populated'] is False:
            return

        name = keyword['name']

        try:
            curval = float(keyword['binary'])
        except Exception as e:
            apflog("Exception in list_mon: %s" % (e), level='error')
            return

        if self.mon_lists[name] == []:
            self.mon_lists[name] = [curval]*300
            self.avg_lists[name] = curval
        else:
            self.mon_lists[name].append(curval)
            self.mon_lists[name] = self.mon_lists[name][-300:]
            self.avg_lists[name] = np.average(self.mon_lists[name])

        return

    def dew_pt_mon(self, dew):
        '''
        dew_pt_mon(dew)
        This is a callback function for the dew point.
        If the mirror temp is too close, logs that and
        sets the dew_too_close flag.
        '''
        if dew['populated'] is False:
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

                if self.dew_too_close:
                    if curm2air - curdew > 3 or curm2 - curdew > 5:
                        self.dew_too_close = False
                else:
                    if curm2air - curdew < 2 or curm2 - curdew < 4:
                        self.dew_too_close = True
                        apflog("M2 temperatures too close to dew point",echo=True,level='warning')

        return

    # Callback for ok2open permission
    # -- Check that if we fall down a logic hole we don't error out
    def ok_mon(self, ok2open):
        """
        ok_mon(ok2open)
        This is a callback function for the ok2open keyword.
        It also monitors other keywords as a unified check
        on whether or not it is ok to open the telescope.
        """
        if ok2open['populated'] is False:
            return
        try:
            ok = ok2open['binary'] # historical
        except Exception as e:
            apflog("Exception in ok_mon for checkapf.OPEN_OK: %s" % (e), level='error')
            return
        try:
            if self.mv_perm.read(binary=False) is False:
                ok = False
        except Exception as e:
            apflog("Exception in ok_mon for checkapf.MOVE_PERM: %s" % (e), level='error')
            return
        try:
            if not self.userkind.read(binary=True) == 3:
                ok = False
        except Exception as e:
            apflog("Exception in ok_mon checkapf.USERKIND: %s" % (e), level='error')
            return
        self.openOK = ok
        return

    def check_sanity(self):
        """
        check_sanity()

        Checks the sanity of the APF system.  This includes checking
        the weather, the dew point, and the wind speed.

        """
        ret_val = True
        # check the sanity of the EOS dispatchers
        pc_keywords = ['PC_EOSCOOLSTA',
                        'PC_EOSCTRLSTA',
                        'PC_EOSDOMESTA',
                        'PC_EOSGCAMSTA',
                        'PC_EOSMETSSTA',
                        'PC_EOSTDIOSTA',
                        'PC_EOSTELESTA',
                        'PC_EOSTI8KSTA',]

        pc_servers = {
            'eosgcam': 'dresden',
            'eostele': 'hamburg',
            'eosdome': 'hamburg',
            'eoscool': 'hamburg',
            'eosctrl' :'hamburg',
            'eosmets': 'hamburg',
            'eostdio': 'hamburg',
            'eosti8k': 'hamburg',
        }

        for kw in pc_keywords:
            try:
                pc_kw = self.apfmon[kw]
                kw_val = pc_kw.read(binary=True)
                if kw_val > 4:
                    # this is a warning or error
                    apflog("PC keyword %s has value %s, recommend restarting" %\
                            (kw,pc_kw['ascii']),level='Crit',echo=True)
                    srv_name = kw[3:-3].lower()
                    restart(srv_name,host=pc_servers[srv_name.lower()])
                    ret_val = False
            except Exception as e:
                apflog("Cannot monitor keyword %s: %s" % (kw,e),echo=True, level='warn')
                ret_val = False

        return ret_val

    def apftask_mon(self, status):
        if status['populated'] is False:
            return
        try:
            status_val = status['binary']
        except:
            return

        # if the task is not running, the runhost
        # keyword is empty, so we need to hardcode 
        # the runhosts for the various tasks
        hosts = dict()
        hosts['METSXFER'] = 'frankfurt.ucolick.org'
        hosts['APFTEQ'] = 'bremen.ucolick.org'


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


    def sun_rising(self):
        """
        sun_rising()

        Is the sun rising, ie., is it after midnight.

        """
        # the sun also rises
        now = datetime.datetime.now()
        if now.strftime("%p") == 'AM':
            return True
        else:
            return False

    def avg_tel_temps(self):
        """
        self.avg_tel_temps()
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

        #self.avgtemps = [self.avg_lists[nm] for nm in:
        # ('TM1S210','TM2CSUR','TAVERAGE','TM2CAIR','TEMPNOW3','TEMPNOW4')]
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

    def pred_tel_focus(self):
        """
        pred_tel_focus()

        Predicts the telescope focus based on the values of
        various telescope temperatures.
        """
        self.avg_tel_temps()
        # m1 m2 tavg m2air tf3 tf4

        # values as of July 3 2022
        # Final_focus_temp_fits.ipynb which (along with the data) should get checked in
        slopes = np.asarray([-0.007930, 0.021897, 0.011436, -0.008821, 0.001830, -0.017760])
        zeropoint_temps = np.asarray([15.549, 14.195, 14.610, 13.285, 15.716, 15.797])

        predfoc = np.sum(slopes*(self.avgtemps-zeropoint_temps)) + TELFOCUSTYP
        predfoc /= 1000.0 # convert to meters
        return predfoc


    # Function for checking what is currently open on the telescope
    def is_open(self):
        """
        is_open()

        Returns the state of checkapf.WHATSOPN as a tuple (bool, str).
        """
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

    def is_ready_observing_direct(self):
        '''
        
        is_ready_observing_direct()
        Checks the state of the mirror cover and dome shutter
        to determine if the telescope is ready for observing.
        Queries the keywords directly rather than going
        through checkapf.WHATSOPN.'''
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
        if isshutterclosed is False:
            if len(what) > 0 :
                what = what + " DomeShutter"
            else:
                what = "DomeShutter"
            if rv:
                rv = True
        else:
            rv = False
        return rv, what

    def is_ready_observing(self):
        """
        is_ready_observing()

        Returns the state of checkapf.WHATSOPN as a tuple (bool, str).
        """
        try:
            whatstr = str(self.whatsopn)
            what = whatstr.split()
        except:
            apflog("checkapf.WHATSOPN returned a value that str.split cannot split",level='warn',echo=True)
            return self.is_ready_observing_direct()


        if hasattr(what,'__iter__'):
            if "DomeShutter" in what and "MirrorCover" in what:
                return True, what
            else:
                return False, ''
        else:
            return self.is_ready_observing_direct()

    def find_star(self):
        """
        find_star()

        Given the telescopes current position, find the closest bright star.
        """
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
        # example line
        #$line[3] $line[4] $line[5] $line[6] $line[7] 5 1 "8"] < /usr/local/lick/data/apf/StarCatalog.dat"
        p = subprocess.Popen(cmdargs, stdin=starcat, stdout=subprocess.PIPE,\
                             stderr=subprocess.PIPE,cwd=os.path.curdir)
        out, err = p.communicate()
        ret_code = p.returncode
        if ret_code != 0:
            apflog(out,echo=True)
            apflog(err, level="warn",echo=True)
            return False

        return out.split()

    def slew(self,star):
        """
        slew(star)

        Executes the slewlock script to a star (not the slew).
        A star object is list and these are assumed to be in the format
        of the bright star list used for pointing model acquisition.
        So the name is the first element, the RA in radians is second,
        the Dec in radians is third, and so forth.

        """
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
        result, code = apftask_do(cmd,cwd=os.path.curdir,debug=True)
        if not result:
            errstr = "Failed at slewlock - check logs could be a slew failure " 
            errstr += "or a failure to acquire: %s" % (code)
            apflog(errstr, level="error", echo=True)

        return result


    def save_movie(self):
        """
        save_movie()

        Starts the recording of a guider move.

        """
        now = datetime.datetime.now()
        self.fits3pre.write('%d%02d%02d_%s_' % (now.year,now.month,now.day,\
                                                 self.tel['TARGNAME'].read()))
        self.fits3dir.write('/data/apfguide')
        self.save3d.write(True)
        return

    def stop_movie(self):
        """
        stop_movie()

        Stops recording a guider movie.

        """
        self.save3d.write(False)
        self.fits3dir.write('/tmp/')
        return


    def run_focustel(self):
        """
        Runs the telescope focus routine.

        """
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

        self.save_movie()

        apflog("Running focus_telescope routine.",echo=True)
        focus_tel_str = "focus_telescope -c %.3f" % (float(self.pred_tel_focus())*1000.0)
        cmd = os.path.join(SCRIPTDIR,focus_tel_str)
        result, code = apftask_do(cmd,cwd=os.path.curdir)

        self.stop_movie()

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

    def run_autoexposure(self, ind=5):
        """
        runs the autoexposure script to configure the guider
        exposure time for the target.
        """
        cmd = os.path.join(SCRIPTDIR,'autoexposure')
        istr = "%d" % (ind)
        cmdargs = cmd

        result, code = apftask_do(cmdargs,cwd=os.path.curdir)

        if not result:
            apflog("autoexposure failed with code %d" % code, echo=True)
        return result

    def run_centerup(self):
        """
        run_centerup()

        Runs the centerup script to center up the telescope on the
        target in the guider.
        """
        cmd = os.path.join(SCRIPTDIR,'centerup')
        result, code = apftask_do(cmd,cwd=os.path.curdir)
        if not result:
            apflog("centerup failed with code %d" % code, echo=True)
        return result

    def find_star_focustel(self):
        """ This finds a star in the catalog of pointing reference stars close
        to the current position of the telescope.
        It then runs slewlock to slew to the star, set up the guider, and center
        the star on the slit.
        Finally, once the star is acquired, it runs the telescope focus routine.
        """
        star = self.find_star()
        if not star:
            apflog("Cannot find star near current position!?",level='error',echo=True)
            return False
        apflog("Targeting telescope on %s" % star[0], echo=True)
        try:
            self.robot['SCRIPTOBS_VMAG'].write(star[6])
        except Exception as e:
            apflog("Cannot write SCRIPTOBS_VMAG: %s" % (e), level='error',echo=True)
        try:
            sline = "%s %f %f pmra=%s pmdec=%s vmag=%s # end" %\
              (star[0],float(star[1])*3.819718,float(star[2])*57.295779,star[4],star[5],star[6])
            self.robot['SCRIPTOBS_LINE'].write(sline)
        except Exception as e:
            apflog("Cannot write SCRIPTOBS_LINE: %s" % (e), level='error',echo=True)

        try:
            self.robot['SCRIPTOBS_LINE_RESULT'].write(0)
            self.robot['SCRIPTOBS_OBSERVED'].write(False)
        except Exception as e:
            errstr = "Cannot reset SCRIPTOBS_LINE_RESULT or SCRIPTOBS_OBSERVED: %s" % (e)
            apflog(errstr, level='error', echo=True)

        apftask_do(os.path.join(SCRIPTDIR,"enable_dome_azdrive"))

        predfocus  = self.pred_tel_focus()
        self.robot['FOCUSTEL_STARTFOCUS'].write(predfocus)
        self.robot['FOCUSTEL_LASTFOCUS'].write(predfocus)
        self.focus.write(self.pred_tel_focus(), binary=True, wait=False)
        if self.slew(star):
            return self.run_focustel()
        return False


    def set_apfteq_mode(self, mode):
        """
        set_apfteq_mode()

        Sets the mode for APFTEQ to the commanded value.
        Does not check if value is allowed, but does if the switch was
        successful.

        """

        apflog("Setting TEQMode to %s" % mode)
        if self.test:
            apflog("Would be setting TEQMode to %s" % mode)
            return
        self.teqmode.write(mode,wait=False)
        result = self.teqmode.waitfor('== %s' % mode, timeout=60)
        if not result:
            apflog("Error setting the TEQMODE.")
            raise RuntimeError("Couldn't set TEQ mode")


    def clear_estop(self):
        """
        clear_estop()

        Clears the Estop state by running the clear_estop script.

        """

        if self.test: return True
        if self.mv_perm.binary is False:
            apflog("Waiting for permission to move...", echo=True)
            chk_move = "$checkapf.MOVE_PERM == true"
            result = APFTask.waitFor(self.task, False, chk_move, timeout=600)
            if not result:
                apflog("Can't open. No move permission.",echo=True)
                return False

        cmd = os.path.join(SCRIPTDIR,'clear_estop')
        result, code = apftask_do(cmd,debug=True,cwd=os.getcwd())
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

    def state_set(self):
        """
        state_set()

        This checks if certain emergency stop states are set.

        """

        # there are three states - but we do not care about ESTOPST,
        # that is will be cleared in openatsunset/openatnight
        if self.dome['ECLOSEST']:
            return True
        if self.dome['ESECURST']:
            return True
        return False


    def home_telescope(self):
        """
        home_telescope()

        Homes the telescope using the slew script.
        Then verifies homing occurred.

        """
        cmd = os.path.join(SCRIPTDIR,"slew") + " --home"
        rv, rc = apftask_do(cmd)
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

    def check_home(self,home=True):
        """
        check_home(home=True)

        Checks if the telescope has been homed.
        If not, homes if the variable home=True.

        """
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
                    self.home_telescope()
                else:
                    apflog("Telescope needs to be homed",level='Alert',echo=True)
                    return False
            else:
                apflog("apfmon.ELHOMERIGHTSTA value is %d" % (homed),level='Alert',echo=True)
                return False

        return False

    def openat(self, sunset=False):
        """
        Function to ready the APF for observing. Calls either openatsunset or openatnight.
           This function will attempt to open successfully twice. If both attempts
           fail, then it will return false, allowing the master to register the error
           and behave accodingly. Otherwise it will return True.
        """
        # If this is a test run, just return True
        if self.test:
            return True

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

        if self.mv_perm.binary is False:
            apflog("Waiting for permission to move...", echo=True)
            chk_move = "$checkapf.MOVE_PERM == true"
            result = APFTask.waitFor(self.task, False, chk_move, timeout=600)
            if not result:
                apflog("Can't open. No move permission.",echo=True)
                return False

        if self.state_set():
            apflog("An unusual emergency state is set.", level="alert",echo=True)
            return False

        # Everything seems acceptable, so lets try opening
        if sunset:
            cmd = os.path.join(SCRIPTDIR,'openatsunset')
        else:
            cmd = os.path.join(SCRIPTDIR,'openatnight')

        # Make two tries at opening. If they both fail return False so the caller can act
        # accordingly.
        result, code = apftask_do(cmd)
        if not result:
            errstr = "First openat attempt has failed. Exit code = %d." % code
            errstr += " After a pause, will make one more attempt."
            apflog(errstr,echo=True)
            APFTask.waitFor(self.task, True, timeout=10)
            result, code = apftask_do(cmd)
            if not result:
                errstr = "Second openup attempt also failed. Exit code = %d. Giving up." % code
                apflog(errstr, echo=True)
                return False
        rv = self.check_home()
        if rv is False:
            return False
        try:
            APFLib.write("eostele.FOCUS",ktl.read("apftask","FOCUSTEL_LASTFOCUS",binary=True))
        except:
            apflog("Cannot move secondary focus.",level="error")
            return False
        return True

    def power_down_telescope(self):
        """
        Checks that we have the proper permission and dome is closed, then resets telescope power.
        This core is a separate script.
        """
        if self.test: return True
        cmd = os.path.join(SCRIPTDIR,"power_down_telescope")
        self.dm_reset()
        if self.mv_perm.binary is False:
            apflog("Waiting for permission to move")
        chk_mv = '$checkapf.MOVE_PERM == true'
        result = APFTask.waitFor(self.task, False, chk_mv, timeout=300)
        if not result:
            apflog("Didn't have move permission after 5 minutes.", echo=True)
            return False
        # one last check

        apflog("Running power_down_telescope script")
        result, _ = apftask_do(cmd)
        if result:
            return True

        apflog("power_down_telescope has failed. Human intervention likely required.",
                level='alert', echo=True)
        return False


    def servo_failure(self):
        """checks for amplifier faults
        Checks all seven moving components.
        """
        servo_failed = False

        estopstate = self.dome.read('ESTOPST',binary=True)
        if estopstate:
            return False

        msg = ""
        prefixs = ["AZ","EL","FA","FB","FC","TR" ]
        for pr in prefixs:
            nm = pr + "AMPFLT"
            val = self.tel[nm].read(binary=True)
            if val:
                servo_failed = True
                msg += "Error: Servo Amplifier Fault: " + str(nm) + " " + str(val) + "\n"

        for pr in prefixs:
            nm = pr + "FERROR"
            val = self.tel[nm].read(binary=True)
            if val:
                servo_failed = True
                msg = "Error: Following Error Fault: " + str(nm) + " " + str(val) + "\n"
        if msg != "":
            apflog(msg, level='error', echo=True)
            self.robot['MASTER_MESSAGE'].write(msg)

        return servo_failed

    def close(self, force=False):
        """Checks that we have the proper permission, then runs the closeup script.
        On failures retries, and also sends alerts. Good for waking people up.
        """

        if self.test: 
            return True
        cmd = os.path.join(SCRIPTDIR,"closeup")
        if force:
            apflog("Calling a single instance of closeup. Will return regardless of result.",\
                    echo=True)
            result, code = apftask_do(cmd)
            return result
        if self.mv_perm.binary is False:
            if self.chk_close.binary is True:
                apflog("Waiting for checkapf to close up")
            else:
                apflog("Waiting for permission to move")
        chk_mv = '$checkapf.MOVE_PERM == true'
        result = APFTask.waitFor(self.task, False, chk_mv, timeout=1200)
        if not result:
            apflog("Didn't have move permission after 20 minutes. Going ahead with closeup.",\
                    echo=True)
            return False
        try:
            if self.apfmon['FRONT_SHUTTER_CLOSEUPSTA'].read(binary=True,timeout=2) != 2:
                # this is a check to see if the front shutter got caught running
                # away, if so do not send any more shutter commands
                apflog("Dome Shutters maybe running away!", level='error', echo=True)
                return False
        except:
            apflog("Cannot communicate with apfmon1, proceeding anyway (fingers crossed)", \
                 level='error', echo=True)


        apflog("Running closeup script")
        attempts = 0
        close_start = datetime.datetime.now()
        while (datetime.datetime.now() - close_start).seconds < 1800 and attempts < 8:
            result = APFTask.waitFor(self.task, False, chk_mv, timeout=300)
            if not result:
                apflog("Didn't have move permission after 5 minutes. ", echo=True)
                break
            attempts += 1
            result, code = apftask_do(cmd)
            if not result:
                apflog("Closeup failed with exit code %d" % code, echo=True)
                if self.servo_failure() or self.slew_allowed.read(binary=True) is False:
                    apflog("Likely Servo amplifier failure, may power cycle telescope",\
                           echo=True,level='error')
                    rv = self.power_down_telescope()
                    if rv:
                        apflog("Power cycled telescope",echo=True,level="error")
                    else:
                        apflog("Failure power cycling telescope",echo=True,level="error")
                if attempts == 7:
                    lstr = "Closeup has failed %d times consecutively. " % (attempts)
                    lstr += "Human intervention likely required."
                    areopen, _ = self.is_open()
                    if areopen is True:
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

        apflog("Closeup could not successfully complete.")
        return False


    def set_tel_foc(self):
        """
        Sets the telescope focus to the predicted value returned by
        pred_tel_focus()
        """

        predfocus  = self.pred_tel_focus()
        self.robot['FOCUSTEL_STARTFOCUS'].write(predfocus)

        if self.mv_perm and self.faenable['binary'] == 1:
            try:
                self.focus.write(predfocus,binary=True,wait=False)
                msg = "Wrote %f to eostele.Focus" % (predfocus*1000.) # in mm
                self.robot['MASTER_MESSAGE'].write(msg)
            except Exception as e:
                apflog("Cannot write eostele.FOCUS: %s" % (e), level="error", echo=True)

    def set_autofoc_val(self):
        """ APFControl.set_autofoc_val()
            tests when the last time the telescope was focused,
            if more than FOCUSTIME enable focus check
        """

        # check last telescope focus
        lastfoc = self.robot['FOCUSTEL_LAST_SUCCESS'].read(binary=True)

        if self.sun_rising() and (self.sunel.read(binary=True) < -20):
            self.autofoc.write("robot_autofocus_disable")
            return 0

        lastopen = self.robot['OPENUP_LAST_SUCCESS'].read(binary=True)

        if lastfoc < lastopen:
            self.autofoc.write("robot_autofucs_enable")
            return 2

        if time.time() - lastfoc < 3600:
            self.autofoc.write("robot_autofocus_disable")
            return 0

        predfocus  = self.pred_tel_focus()
        self.robot['FOCUSTEL_STARTFOCUS'].write(predfocus)
        focus_diff = math.fabs(predfocus - self.focus['binary'])
        focus_diff *= 1e3

        #self.focus.write(predfocus,binary=True,wait=False)
        predfocus *= 1e3
        self.autofoc.write("robot_autofocus_enable")
        focval = 1

        ostr = "Current telescope focus more than %5.3f mm" % focus_diff
        ostr += "from predicted, setting to %5.3f." % predfocus
        APFTask.set(self.task, suffix="MESSAGE", value=ostr, wait=False)

        return focval

    def update_windshield(self, state):
        """Checks the current windshielding mode.
        If the input state is auto, makes sure the mode is set properly 
        based on wind speed and temperature.
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
            # This state enables or disables windshielding based on the 
            # wind speed and the outside temperature
            #if self.down > 0:
            #    wvel = self.avg_lists['M3WIND']
            #else:
            wvel = self.avg_lists['WINDAV']

            apflog("Current median wind speed is %.2f with the limit %.2f" % \
                   (wvel,WINDSHIELD_LIMIT), level='debug')
            if currMode == 'enable' and wvel <= WINDSHIELD_LIMIT and \
                float(self.avg_lists['TEMP']) > TEMP_LIMIT:
                apflog("Setting scriptobs_windshield to Disable")
                rv = "Disable"
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], rv)

            if currMode == 'disable' and (wvel > WINDSHIELD_LIMIT or \
                                          float(self.avg_lists['TEMP']) < TEMP_LIMIT):
                apflog("Setting scriptobs_windshield to Enable")
                rv = "Enable"
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], rv)

        return rv


    def check_FCUs(self, check_apfmon=False):
        """
        If the dome is open the FCs should be off. Sometimes they are commanded 
        off but do not turn off.
        This will attempt to fix that by checking if the dome is open, checking 
        if the FCs are on, and then turning them on, then off. Which is the only
        way to turn them off.
        """

        if check_apfmon:
            status_value = self.apfmon['FC_STATUSSTA'].read(binary=True)
            if status_value < 4:
                return

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

    def run_prepobs(self, evening=False):
        """
        run_prepobs()
        Runs the prep-obs script to prepare the telescope for observing.

        evening: If True, run prep-obs with the --evening flag.
        """
        apflog("Calling prep-obs.",echo=True)
        prepobs = os.path.join(SCRIPTDIR,'prep-obs')
        if evening:
            prepobs += " --evening"
        result, ret_code = apftask_do(prepobs)
        if result is False:
            # try again
            self.dm_reset()
            result, ret_code = apftask_do(prepobs)
            if result is False:
                log_str = "Prep-obs returned error code %d. " % (ret_code)
                log_str += "Targeting object has failed."
                apflog(log_str,level='error',echo=True)
                return
        return


    def evening_star(self):
        """Aim the APF at the desired target. This calls prep-obs, slewlock, and focus-telescope."""
        if self.is_open()[0] is False:
            apflog("APF is not open. Can't target a star while closed.",level='error',echo=True)
            return
        self.dm_reset()

        if self.apf.calsta['binary'] < 3 or self.apf.focussta['binary'] < 3:
            log_str = 'Focusinstr and/or Calibrate are running, will skip evening star observation.'
            log_str += ' focusinstr=%s calibrate=%s' % (self.apf.calsta, self.apf.focussta)
            apflog(log_str,echo=True)
            return

        # check on weirdness for UCAM host post-reboot
        self.apf.ucam_dispatch_mon()
        self.run_prepobs(evening=True)

        self.dm_reset()
        apflog("Slewing to lower el",echo=True)
        result, ret_code = apftask_do('slew -e 75')
        if result is False:
            apflog("Slew returned error code %d. Targeting object has failed."\
                    % (ret_code),level='error',echo=True)
            return
        # Slew to the specified RA and DEC, set guide camera settings, and centerup( Slewlock )
        # Focus the telescope - all of this, including finding the star, is done in focusTel
        self.dm_reset()
        if self.find_star_focustel():
            try:
                self.robot['SCRIPTOBS_LINE_RESULT'].write(3)
                self.robot['SCRIPTOBS_OBSERVED'].write(True)
                self.guide['CLEARSUMS'].write('now')
                self.guide['CLEARSUMS'].write('gstate')
            except Exception as e:
                log_str = "Cannot write 3 to SCRIPTOBS_LINE_RESULT or True"
                log_str += " to SCRIPTOBS_OBSERVED: %s" % (e)
                apflog(log_str, level='warn', echo=True)
            return True
        else:
            try:
                self.robot['SCRIPTOBS_LINE_RESULT'].write(2)
            except Exception as e:
                apflog("Cannot write 2 to SCRIPTOBS_LINE_RESULT: %s" % (e), level='warn', echo=True)
            return False

    def dm_reset(self):
        """
        dm_reset()

        Writes to ROBOSTATE keyword, effecitively tapping the deadman switch.

        """
        try:
            APFLib.write(self.checkapf['ROBOSTATE'], "master operating",timeout=10)
        except Exception as e:
            try:
                ukind = self.userkind.read()
            except:
                ukind = "Unknown"
            ostr = "Error: Cannot write to ROBOSTATE, USERKIND = %s, reason: %s" % (ukind,e)
            apflog(ostr,level='error',echo=True)

    def dm_zero(self):
        """
        dm_zero()

        Sets deadman timer to 0 (-1).
        If the telescope or dome is open, checkapf will force it closed.
        """
        try:
            if self.checkapf['DMTIME'].read(binary=True) < 1:
                APFLib.write(self.dmtimer, -1,timeout=10)
        except Exception as e:
            ostr = "Error: cannot touch DM Timer: %s " %( e)
            apflog(ostr,level='error',echo=True)

def main():
    """
    Main function to test the telescope monitors.
    """
    print("Testing telescope monitors, grabbing and printing out current state.")

    task = 'example'

    APFTask.establish(task, os.getpid())
    apf = APF(task=task,test=False)
    tel = TelescopeControl(apf, test=False)

    # Give the monitors some time to start up
    APFTask.waitFor(task, True,timeout=2)

    print(str(tel))

    while True:
        APFTask.wait(task,True,timeout=10)
        print(str(tel))


if __name__ == '__main__':
    main()
