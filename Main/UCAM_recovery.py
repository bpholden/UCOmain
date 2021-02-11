from __future__ import print_function
import os
import sys
import time
import re
from datetime import datetime, timedelta

try:
    from apflog import *
    import ktl
    import APFTask
except:
    from fake_apflog import *



class UCAM_recovery():

    def __init__(self, task='example', test=False):
        self.test = test
        self.task = task

        self.apfucam = ktl.Service('apfucam')
        self.apftask = ktl.Service('apftask')
        self.apfmot = ktl.Service('apfmot')
        
        self.combo_ps = self.apfucam['combo_ps']
        self.combo_ps.monitor()
        self.ctalkto = self.apfucam['ctalkto']
        self.ctalkto.monitor()
        self.ucam_command = self.apftask['UCAMLAUNCHER_UCAM_COMMAND']
        self.ucam_command.monitor()
        self.ucam_status = self.apftask['UCAMLAUNCHER_UCAM_STATUS']
        self.ucam_status.monitor()


    def __repr__(self):
        rstr = "task=%s combops=%s ctalkto=%s command=%s status=%s " % (self.task,self.combo_ps,self.ctalkto,self.ucam_command,self.ucam_status)
        return rstr
                
    def reboot_warsaw(self):

        if self.ucam_status.read(binary=True) == 1:
            rv = self.stop_ucam_software()
            if rv is False:
                apflog("UCAM software did not stop running, rebooting anyway",level='Error',echo=True)
            
        self.ucam_command.write(2)
        rv = APFTask.waitfor(self.task, True, expression="$apftask.UCAMLAUNCHER_STATUS == Running", timeout=120)
        if rv:
                # yay!
            self.ucam_command.write(1)
            rv = self.power_cycle_fousb()
            return rv
        else:
            # this is bad
            apflog("UCAM host not re-booted",level='Alert',echo=True)

        return False

    def power_cycle_fousb(self):

        self.apfmot['FOUSB_POWER'].write(0)
        rv = APFTask.waitfor(self.task, True, expression="$apfmot.FOUSB_POWER == Off", timeout=10)
        if rv:
            APFTask.wait(self.task,True,timeout=1)
        else:
            apflog("Cannot power cycle FOUSB",level='error',echo=True)
            return False
        self.apfmot['FOUSB_POWER'].write(1)
        rv = APFTask.waitfor(self.task, True, expression="$apfmot.FOUSB_POWER == On",timeout=10)
        if rv:
            APFTask.wait(self.task,True,timeout=1)
        else:
            apflog("Cannot power cycle FOUSB",level='error',echo=True)
            return False
        return rv

    def stop_ucam_software(self):
        
        if self.ucam_status.read(binary=True) != 0:
            self.ucam_command.write(0)
            rv = APFTask.waitfor(self.task, True, expression="$apftask.UCAMLAUNCHER_UCAM_STATUS == 'Not running'", timeout=3)
            return rv
        else:
            return False

    def start_ucam_software(self):
        
        if self.ucam_status.read(binary=True) == 0:
            self.ucam_command.write(1)
            rv = APFTask.waitfor(self.task, True, expression="$apftask.UCAMLAUNCHER_UCAM_STATUS == 'Running'", timeout=3)
            return rv
        else:
            return False
        
if __name__ == "__main__":
    uc = UCAM_recovery()
    print(uc)
    
