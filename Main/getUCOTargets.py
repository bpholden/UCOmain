from __future__ import print_function
from datetime import datetime
import os
import os.path
import threading
import shutil

import astropy.io.ascii
import numpy as np

import ktl
import APFTask

from apflog import apflog
import ParseUCOSched
import UCOScheduler as ds


class getUCOTargets(threading.Thread):

    def __init__(self, opt, task='master',prilim=0.5,certificate='UCSC_Dynamic_Scheduler-4f4f8d64827e.json',wait_time=0):
        threading.Thread.__init__(self)

        self.task = task
        self.apftask = ktl.Service('apftask')
        self.rank_table = opt.rank_table
        self.time_left = opt.time_left

        self.too = None
        self.sheets = None
        self.signal = True
        self.timeout = 1200
        self.proceed = False
        self.proceed_timeout = 120
        self.reading = False
        self.wait_time = wait_time

        self.prilim = prilim
        self.certificate = certificate

        if opt.test:
            self.debug = opt.test
        else:
            self.debug = False

        if self.rank_table is None and self.sheets is None:
            apflog("Error: no rank table and no sheet list provided", level='error')
            return

        self.start()

    def stop(self):
        self.signal = False
        threading.Thread._Thread__stop(self)

    def run(self):

        # these tests are against binary values
        # and binary values are in radians
        # the test is for -4 degrees
        expression = '$eostele.SUNEL < -0.0'
        APFTask.waitFor(self.task, True, expression=expression, timeout=self.wait_time)

        if self.signal is False:
            return

        if self.debug:
            print("Would have downloaded %s" % (self.rank_table))
        else:
            if self.time_left is None :
                hour_constraints=None
            else:
                if os.path.exists(self.time_left):
                    try:
                        hour_constraints = astropy.io.ascii.read(self.time_left)
                    except Exception as e:
                        hour_constraints = None
                        apflog("Error: Cannot read file of time left %s : %s" % (opt.time_left,e))
                else:
                    hour_constraints = None

            try:
                rank_table = ds.make_rank_table(self.rank_table,outdir=os.getcwd(),hour_constraints=hour_constraints)
            except Exception as e:
                apflog("Error: Cannot download rank_table?! %s" % (e),level="error")
                # goto backup
                if os.path.exists("rank_table.1"):
                    shutil.copyfile("rank_table.1","rank_table")
                    try:
                        rank_table = ds.make_rank_table(sheet_table_name=opt.rank_table,outdir=os.getcwd(),hour_constraints=hour_constraints)
                    except Exception as e:
                        apflog("Error: Cannot reuse rank_table?! %s" % (e),level="error")
                        rank_table = None

            if self.signal is False:
                return

            if self.sheets is None:
                self.sheets = list(rank_table['sheetn'][rank_table['rank'] > 0])

            try:
                self.apftask.write('MASTER_SHEETLIST',",".join(self.sheets),timeout=2)
            except Exception as e:
                apflog("Cannot write apftask.MASTER_SHEETLIST: %s" % (e), level='warn',echo=True)

            if 'too' in rank_table.columns:
                if np.any(rank_table['too']):
                    self.too = list(rank_table['sheetn'][rank_table['too']])

            if self.debug:
                print("Would have downloaded %s" % (" ".join(self.sheets)))
                print("TOO sheets found are: %s" % (" ".join(self.too)))
            else:
                if self.signal is False:
                    return
                try:
                    star_table,stars = ParseUCOSched.parse_UCOSched(sheetns=self.sheets,outfn='googledex.dat',
                                                                   outdir=os.getcwd(),prilim=self.prilim,certificate=self.certificate)
                except Exception as e:
                    apflog("Error: Cannot download googledex?! %s %s" % (type(e), e),level="error")
                    # goto backup
                    if os.path.exists("googledex.dat.1"):
                        shutil.copyfile("googledex.dat.1","googledex.dat")

            if self.debug:
                print("Would have made hour table")
            else:
                try:
                    hour_table = ds.make_hour_table(rank_table,datetime.now(),hour_constraints=hour_constraints)
                except Exception as e:
                    hour_table = None
                    apflog("Error: Cannot make hour_table?! %s" % (e),level="error")

        while self.signal and self.too is not None:

            if APFTask.waitfor(self.task,False,expression='apftask.SCRIPTOBS_PHASE==Observing',timeout=self.timeout):

                self.reading = True
                try:
                    ParseUCOSched.parse_TOO(too_sheetns=self.too,outfn='googledex.dat',outdir=os.getcwd(),prilim=self.prilim,certificate=self.certificate)
                except Exception as e:
                    apflog("Error: Cannot download %s: %s" % (self.too,e),level="error")
                self.reading = False

        self.stop()

        return


if __name__ == "__main__":

    class Opt:
        pass


    task = 'example'
    APFTask.establish(task,os.getpid())
    opt = Opt()
    opt.test = True
    opt.time_left = "/home/holden/time_left.csv"
    opt.rank_table = '2023A_ranks'

    get_targs = getUCOTargets(opt, task=task)
