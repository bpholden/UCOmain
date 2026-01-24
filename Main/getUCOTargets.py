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
        self.start()

    def run(self):
        apflog("Starting getUCOTargets for %s" % (self.uco_targets.rank_table_name), echo=True)
        expression = '$eostele.SUNEL < -0.0'
        if self.debug:
            expression = '$eostele.SUNEL < 100.0'
        apflog("Waiting for %s to be true" % (expression), echo=True)
        APFTask.waitFor(self.task, True, expression=expression, timeout=self.wait_time)

        if self.signal is False:
            return
        apflog("getUCOTargets making rank and hour tables", echo=True)
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

        if self.signal is False:
            return
        apflog("getUCOTargets making star table", echo=True)
        self.uco_targets.make_star_table()

        too_sheetlist_name = 'MASTER_TOOSHEETLIST'
        if self.debug:
            too_sheetlist_name = 'EXAMPLE_VAR_2'
        try:
            if self.uco_targets.too_sheets is not None:
                self.apftask.write(too_sheetlist_name, ",".join(self.uco_targets.too_sheets), timeout=2)
            else:
                self.apftask.write(too_sheetlist_name, "", timeout=2)
        except Exception as e:
            apflog("Cannot write apftask.%s: %s" % (too_sheetlist_name, e), level='warn',echo=True)

        while self.signal and self.uco_targets.too_sheets is not None and not self.debug and False:

            if APFTask.waitfor(self.task, False, expression='apftask.SCRIPTOBS_PHASE==Observing',
                               timeout=self.timeout):

                self.reading = True
                try:
                    ParseUCOSched.parse_TOO(too_sheetns=self.uco_targets.too_sheets, outfn=self.uco_targets.star_table_name,
                                            outdir=os.getcwd(), prilim=self.uco_targets.prilim,
                                            certificate=self.uco_targets.certificate)
                except Exception as e:
                    apflog("Error: Cannot download %s: %s" % (self.uco_targets.too_sheets,e),level="error")
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
