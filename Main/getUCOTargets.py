from __future__ import print_function
from datetime import datetime, timedelta
import os
import os.path
import signal
from select import select
import threading
import time

import astropy.io.ascii

import ktl
import APFTask

from apflog import apflog
import ParseUCOSched
import UCOScheduler as ds


class getUCOTargets(threading.Thread):

    def __init__(self, opt, task='master',prilim=0.5,certificate='UCSC_Dynamic_Scheduler-4f4f8d64827e.json',wait_time=0):
        threading.Thread.__init__(self)

        self.task = task

        if opt.too:
            self.too = opt.too
        else:
            self.too = None

        if opt.rank_table:
            self.rank_table = opt.rank_table
        else:
            self.rank_table = None

        if opt.time_left:
            self.time_left = opt.time_left
        else:
            self.time_left = None

        if opt.sheet:
            self.sheets = opt.sheet
        else:
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
        expression = '$eostele.SUNEL < -0.0698'
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
                rank_table = ds.makeRankTable(self.rank_table,outdir=os.getcwd(),hour_constraints=hour_constraints)
            except Exception as e:
                apflog("Error: Cannot download rank_table?! %s" % (e),level="error")
                # goto backup
                if os.path.exists("rank_table.1"):
                    shutil.copyfile("rank_table.1","rank_table")
                    try:
                        rank_table = ds.makeRankTable(sheet_table_name=opt.rank_table,outdir=os.getcwd(),hour_constraints=hour_constraints)
                    except Exception as e:
                        apflog("Error: Cannot reuse rank_table?! %s" % (e),level="error")
                        rank_table = None

            if self.signal is False:
                return

            if self.sheets is None:
                self.sheets = list(rank_table['sheetn'][rank_table['rank'] > 0])

            if self.debug:
                print("Would have downloaded %s" % (self.sheets))
            else:
                if self.signal is False:
                    return
                try:
                    star_table,stars = ParseUCOSched.parseUCOSched(sheetns=self.sheets,outfn='googledex.dat',
                                                                   outdir=os.getcwd(),prilim=self.prilim,certificate=self.certificate)
                except Exception as e:
                    apflog("Error: Cannot download googledex?! %s" % (e),level="error")
                    # goto backup
                    if os.path.exists("googledex.dat.1"):
                        shutil.copyfile("googledex.dat.1","googledex.dat")

            if self.debug:
                print("Would have made hour table")
            else:
                try:
                    hour_table = ds.makeHourTable(rank_table,datetime.now(),hour_constraints=hour_constraints)
                except Exception as e:
                    hour_table = None
                    apflog("Error: Cannot make hour_table?! %s" % (e),level="error")

        while self.signal and self.too is not None:

            if APFTask.waitfor(self.task,False,expression='apftask.SCRIPTOBS_PHASE==Observing',timeout=self.timeout):

                self.reading = True
                try:
                    ParseUCOSched.parseTOO(too_sheetns=self.too,outfn='googledex.dat',outdir=os.getcwd(),prilim=self.prilim,certificate=self.certificate)
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
    opt.test = False
    opt.too = None
    opt.time_left = "/home/holden/time_left.csv"
    opt.rank_table = '2022A_ranks'

    get_targs = getUCOTargets(opt, task=task)
