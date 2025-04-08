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

windlim = 40.0
slowlim = 100
WINDSHIELD_LIMIT = 10. # mph at the APF
FOCUSTIME = 3600. # minimum time before checking telescope focus
TEMP_LIMIT = 35. # deg F at the APF
wxtimeout = datetime.timedelta(seconds=1800)
SUNEL_HOR = -3.2
TELFOCUSMIN = -0.00096
TELFOCUSMAX = -0.00060
# this value comes an average over many measurements of the telescope focus
#TELFOCUSTYP = -0.83529
TELFOCUSTYP = -0.76529
TELFOCUSMAXOFF = 0.00002

if "LROOT" in os.environ:
    LROOT = os.environ["LROOT"]
else:
    LROOT = '/usr/local/lick'
SCRIPTDIR = os.path.join(LROOT,'bin/robot/')

class TelescopeControl:
    def __init__(self, apf=None):
        self.apf = apf
        
        self.lastfocuscheck = datetime.datetime.now() - datetime.timedelta(days=1)
        self.lastfocusval = None

        self.tel        = ktl.Service('eostele')
        self.sunel      = self.tel('SUNEL')
        self.ael        = self.tel('AEL')
        self.aaz        = self.tel('AAZ')
        self.aafocus    = self.tel('AAFOCUS')
        self.focus      = self.tel('FOCUS')
        self.faenable   = self.tel('FAENABLE')

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
        self.instr_perm = self.checkapf('INSTR_PERM')
        self.chk_close  = self.checkapf('CHK_CLOSE')

        self.apfmet     = ktl.Service('met3apf')
        self.wx         = self.apfmet('M5WIND')
        self.airtemp    = self.apfmet('M5OUTEMP')
        self.down       = self.apfmet('M5DOWN')
        self.altwx      = self.apfmet('M3WIND')

        self.eosti8k    = ktl.Service('eosti8k')
        self.m2tempkw   = self.eosti8k('TM2CSUR')
        self.m2airkw    = self.eosti8k('TM2CAIR')
        self.m1tempkw   = self.eosti8k('TM1S210')
        self.taveragekw = self.eosti8k('TAVERAGE')
        self.t045kw     = self.eosti8k('TTRUS045')
        self.t135kw     = self.eosti8k('TTRUS135')
        self.t225kw     = self.eosti8k('TTRUS225')
        self.t315kw     = self.eosti8k('TTRUS315')


        self.eoscool    = ktl.Service('eoscool')
        self.dewpt      = self.eoscool('DEWPAVG3')
        self.temp3now   = self.eoscool('TEMPNOW3')
        self.temp4now   = self.eoscool('TEMPNOW4')
