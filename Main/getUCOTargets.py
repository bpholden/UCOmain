from __future__ import print_function
import time
import os
import threading
import sys

import numpy as np
import ktl
import APFTask

from apflog import apflog
import ParseUCOSched
import UCOTargets


class getUCOTargets(threading.Thread):

    def __init__(self, uco_targets, task='master', wait_time=0):
        threading.Thread.__init__(self)

        self.uco_targets = uco_targets
        self.debug = uco_targets.debug
        self.wait_time = wait_time
        self.task = task
        self.apftask = ktl.Service('apftask')
        self.timeout = 1200
        self.reading = False
        self.signal = True
        self.too = None
        self.start()

    def run(self):

        expression = '$eostele.SUNEL < -0.0'
        if self.debug:
            expression = '$eostele.SUNEL < 100.0'
        APFTask.waitFor(self.task, True, expression=expression, timeout=self.wait_time)

        if self.signal is False:
            return

        if self.uco_targets.rank_table_name is None:
            return

        self.uco_targets.make_hour_table()

        if self.signal is False:
            return

        if self.uco_targets.sheets is None and self.uco_targets.rank_table:
            self.uco_targets.sheets = list(self.uco_targets.rank_table['sheetn'][self.uco_targets.rank_table['rank'] > 0])

        sheetlist_name = 'MASTER_SHEETLIST'
        if self.debug:
            sheetlist_name = 'EXAMPLE_VAR_1'

        try:
            self.apftask.write(sheetlist_name,",".join(self.uco_targets.sheets),timeout=2)
        except Exception as e:
            apflog("Cannot write apftask.%s: %s" % (sheetlist_name, e), level='warn',echo=True)

        if self.uco_targets.rank_table and 'too' in self.uco_targets.rank_table.columns:
            if np.any(self.uco_targets.rank_table['too']):
                self.too = list(self.uco_targets.rank_table['sheetn'][self.uco_targets.rank_table['too']])

        if self.signal is False:
            return

        self.uco_targets.make_star_table()

        while self.signal and self.too is not None and not self.debug and False:

            if APFTask.waitfor(self.task, False, expression='apftask.SCRIPTOBS_PHASE==Observing',
                               timeout=self.timeout):

                self.reading = True
                try:
                    ParseUCOSched.parse_TOO(too_sheetns=self.too, outfn=self.star_tab,
                                            outdir=os.getcwd(), prilim=self.prilim,
                                            certificate=self.certificate)
                except Exception as e:
                    apflog("Error: Cannot download %s: %s" % (self.too,e),level="error")
                self.reading = False
            APFTask.waitFor(self.task, True, expression='apftask.SCRIPTOBS_PHASE==Input',
                            timeout=self.timeout)

        self.signal = False
        
def main():

    class Opt:
        def __init__(self):
            self.test = True
            self.time_left = "/home/holden/time_left.csv"
            self.rank_table = '2025B_ranks_operational'


    task = 'example'
    APFTask.establish(task, os.getpid())
    opt = Opt()
    uco_targets = UCOTargets.UCOTargets(opt)
    print(uco_targets)
    gt = getUCOTargets(uco_targets, task=task)
    while gt.signal:
        try:
            APFTask.wait(task, True, timeout=100)
        except KeyboardInterrupt:
            apflog("%s has been killed by user." % (gt.uco_targets), echo=True)
            sys.exit()
    print("Done")


if __name__ == "__main__":
    main()
