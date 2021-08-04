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
        
    def run(self):

        now = time.time()
        if now < self.stime:
            APFTask.wait(parent, True, timeout=self.stime-now)
        
