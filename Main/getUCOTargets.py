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
        if opt.sheet:
            self.sheets = opt.sheet
        else:
            self.sheets = ['RECUR_A100']
        if opt.too:
            self.too = opt.too
        else:
            self.too = None

        if opt.rank_table:
            self.rank_table = opt.rank_table
        else:
            self.rank_table = None

        if opt.frac_table:
            self.frac_table = opt.frac_table
        else:
            self.frac_table = None

        if opt.time_left:
            self.time_left = opt.time_left
        else:
            self.time_left = None

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

        if self.signal:
            if self.debug:
                print("Would have downloaded %s" % (self.sheets))
            else:
                try:
                    star_table,stars = ParseUCOSched.parseUCOSched(sheetns=self.sheets,outfn='googledex.dat',
                                                                   outdir=os.getcwd(),prilim=self.prilim,certificate=self.certificate)
                except Exception as e:
                    apflog("Error: Cannot download googledex?! %s" % (e),level="error")
                    # goto backup
                    if os.path.exists("googledex.dat.1"):
                        shutil.copyfile("googledex.dat.1","googledex.dat")
            if self.debug:
                print("Would have downloaded %s" % (self.rank_table))
            else:
                try:
                    rank_table = ds.makeRankTable(self.rank_table,outdir=os.getcwd())
                except Exception as e:
                    apflog("Error: Cannot download rank_table?! %s" % (e),level="error")
                    # goto backup
                    if os.path.exists("rank_table.1"):
                        shutil.copyfile("rank_table.1","rank_table")

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

            if self.debug:
                print("Would have downloaded %s" % (self.frac_table))
            else:
                try:
                    hour_table = ds.makeHourTable(self.frac_table,datetime.now(),outdir=os.getcwd(),hour_constraints=hour_constraints)
                except Exception as e:
                    hour_table = None
                    apflog("Error: Cannot download frac_table?! %s" % (e),level="error")
                if hour_table is None and os.path.exists("frac_table.1"):
                    shutil.copyfile("frac_table.1","frac_table")
                    hour_table = ds.makeHourTable(self.frac_table,datetime.now(),outdir=os.getcwd(),frac_table='frac_table',hour_constraints=hour_constraints)

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
    opt.rank_table = '2021B_ranks'
    opt.frac_table = '2021B_frac'
    opt.sheet = "RECUR_A100,2021B_A000,2021B_A001,2021B_A002,2021B_A003,2021B_A005,2021B_A006,2021B_A007,2021B_A008,2021B_A010,2021B_A011".split(",")


    get_targs = getUCOTargets(opt, task=task)
