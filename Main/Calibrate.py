from __future__ import print_function
from datetime import datetime, timedelta
import os
import os.path
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
from apflog import *

AVERAGE_INSTRFOC = 8522

class Calibrate(threading.Thread):
    def __init__(self, apf, name, wait_time, calfile, phase_index=0, task='master', test=False, possible_phases=['Init','Focus','Cal-Pre','Watching','Cal-Post','Focus-Post']):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.apf = apf
        self.task = task
        self.wait_time = wait_time
        self.user = name
        self.owner = 'public'
        self.test = test
        self.calfile = calfile
        self.possible_phases = possible_phases
        self.phase_index = phase_index

        self.name = 'Calibrate'
        self.signal = True
        self.start()


    def testBias(self):

        if self.test:
            apflog("Would have taken a single bias frame using APFControl.testBias()",echo=True)
        else:
            result = self.apf.testBias()
            if result == None:
                apflog("Focusinstr or calibrate or scriptobs are running?!", level='Error', echo=True)
            if result == False:
                # this is a UCAM problem
                rv = self.apf.ucamRestart()
                if rv == False:
                    apflog("Failure in UCAM status and restart!", level='Alert', echo=True)


        result = self.apf.ucamStatus()
        if result is False:
            apflog("Failure in UCAM status and restart!", level='Alert', echo=True)

        return result

    def focusInstr(self):

        if self.test:
            apflog("Would have run APFControl.focusinstr",echo=True)
            result = True
        else:

            result = self.apf.focusinstr()
            apflog("Focus has finished.",echo=True)

        if not self.test:
            APFTask.set(self.task, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())

        return result

    def calibrate(self,phase):

        eostele = ktl.Service('eostele')
        sunel = eostele["sunel"].read()
        sunel = float(sunel)
        if sunel < 3:
            apflog("Not Starting calibrate %s script, sun too low." % (phase), level='Info', echo=True)
            return True

        apflog("Starting calibrate %s script." % (phase), level='Info', echo=True)
        if self.test:
            apflog("Would have waited for permission (APFControl.instrPermit()) for phase %s" % (phase),echo=True)
        else:
            self.apf.instrPermit()

        if self.test:
            apflog("Would have run APFControl.ucamStatus() for phase %s" % (phase),echo=True)
            result = True
        else:
            result = self.apf.ucamStatus()
            if result is False:
                apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
                return False

        time = phase[4:].lower()

        if self.test:
            apflog("Would have run APFControl.calibrate for time %s" % (time),echo=True)
            result = True
        else:
            result = self.apf.calibrate(script = self.calfile, time = time)
            if not self.test:
                APFTask.set(self.task, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())

        if result == False:
            apflog("Calibrate Pre has failed. Trying again",level='warn',echo=True)
            self.apf.instrPermit()
            result = self.apf.calibrate(script = opt.calibrate, time = 'pre')
            if not result:
                apflog("Error: Calibrate Pre has failed twice. Observer is exiting.",level='error',echo=True)
                self.apf.turnOffLamps()

        return result


    def run(self):

        apflog("Will start with phase %d %s after %.1f seconds" % (self.phase_index,self.possible_phases[self.phase_index],self.wait_time), echo=True)

        if self.wait_time > 0:
            APFTask.wait(self.task, True, timeout=self.wait_time)

        start = self.phase_index
        end = self.possible_phases.index('Watching')

        for pi in range(start,end):
            self.phase_index = pi
            cur_phase = self.possible_phases[pi]
            APFTask.phase(self.task, self.possible_phases[pi])
            apflog("Phase now %s %d" % (cur_phase,pi),echo=True)
            
            if pi == 0:
                result = self.testBias()
                if result is False:
                    self.stop()
                    return
            elif pi == 1:
                result = self.focusInstr()
            elif pi == 2:
                result = self.calibrate(cur_phase)

            if result:
                apflog("Phase %s is complete" % cur_phase,echo=True)
            else:
                apflog("Phase %s failed" % cur_phase,echo=True)
                return

        self.stop()
        return


    def stop(self):
        self.signal = False
        threading.Thread._Thread__stop(self)


if __name__ == "__main__":

    task = 'example'
    APFTask.establish(task, os.getpid())
    apf = APFControl.APF(task=task,test=True)
    # Give the monitors some time to start up
    APFTask.waitFor(task, True,timeout=2)
    print(str(apf))

    calibrate = Calibrate(apf,'public',5,task=task,test=True)
    while calibrate.signal:
        try:
            APFTask.wait(task,True,timeout=1)
        except KeyboardInterrupt:
            apflog("%s has been killed by user." % (calibrate.name), echo=True)
            calibrate.stop()
            sys.exit()
        except:
            apflog("%s killed by unknown." % (calbrate.name), echo=True)
            calibrate.stop()
            sys.exit()
