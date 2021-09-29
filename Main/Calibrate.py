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

class Observe(threading.Thread):
    def __init__(self, apf, name, stime, task='master'):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.apf = apf
        self.task = task
        self.stime = stime
        self.user = name
        self.owner = 'public'

        self.name = 'Calibrate'
        self.signal = True
        self.start()


    def testBias(self):
        
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

        if result:
            result = self.apf.focusinstr()
            self.apflog("Focus has finished. Setting phase to Cal-Pre")
            
        if not debug:
            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())

        return result

    def calibratePre(self):
                apflog("Starting calibrate pre script.", level='Info', echo=True)
        apf.instrPermit()

        result = apf.ucamStatus()
        if result is False:
            apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
            shutdown()
        
        result = apf.calibrate(script = opt.calibrate, time = 'pre')
        if not debug:
            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())

        if result == False:
            apflog("Calibrate Pre has failed. Trying again",level='warn',echo=True)
            apf.instrPermit()
            result = apf.calibrate(script = opt.calibrate, time = 'pre')
            if not result:
                apflog("Error: Calibrate Pre has failed twice. Observer is exiting.",level='error',echo=True)
                apf.turnOffLamps()
                sys.exit(2)

        phase_index += 1

        APFTask.phase(parent, possible_phases[phase_index])
        apflog("Calibrate Pre has finished. Setting phase to %s." % phase)

        
    def run(self):

        now = time.time()
        timeout = int(self.stime-now)
        if now < self.stime:
            APFTask.wait(parent, True, timeout=timeout)
        
        bias_result = self.testBias()
        if bias_result:
            focus_result = self.focusInstr()

            if focus_result:
                phase_index += 1
        
                APFTask.phase(parent, possible_phases[phase_index])
                apflog("Phase now %s" % phase)

                calibrate_result = self.calibratePre()

        return
        

