from __future__ import print_function
from datetime import datetime, timedelta
import os
import os.path
import signal
from select import select
import threading
import time

import apflog
import ParseUCOSched

class getUCOTargets(threading.Thread):

    def __init__(self, opt, task='master',prilim=0.5,certificate='UCSC_Dynamic_Scheduler-4f4f8d64827e.json'):
        threading.Thread.__init__(self)
        
        self.task = task
        if opt.sheet:
            self.sheets = opt.sheet
        else:
            self.sheets = ['RECUR_A100']
        if opt.too:
            self.too = opt.too
        else:
            opt.too = None

        if opt.rank_table:
            self.rank_tablen = opt.rank_table
        else:
            self.rank_tablen = None
            
        if opt.frac_table:
            self.frac_tablen = opt.frac_table
        else:
            self.frac_tablen = None

        if opt.time_left:
            self.time_left = opt.time_left
        else:
            self.time_left = None

        self.signal = True
        self.timeout = 1200
        self.proceed = False
        self.proceed_timeout = 120
        self.reading = False
        
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

        
        if self.signal:
            try:
                star_table,stars = ParseUCOSched.parseUCOSched(sheetns=self.sheets,outfn='googledex.dat',
                                                                   outdir=os.getcwd(),prilim=self.prilim,certificate=self.certificate)
            except Exception as e:
                apflog("Error: Cannot download googledex?! %s" % (e),level="error")
                # goto backup
                if os.path.exists("googledex.dat.1"):
                    shutil.copyfile("googledex.dat.1","googledex.dat")

            try:
                rank_table = ds.makeRankTable(sheet_table_name=self.rank_table,outdir=os.getcwd(),certificate=self.certificate)
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
                    
            try:
                hour_table = ds.makeHourTable(self.frac_table,datetime.now(),outdir=os.getcwd(),hour_constraints=hour_constraints,
                                                  certificate=self.certificate)
            except Exception as e:
                apflog("Error: Cannot download frac_table?! %s" % (e),level="error")

        while self.signal and self.too is not None:
            time.sleep(self.timeout)
            # there is a horrible race condition here
            # we need to download the TOO when getNext is not running
            # download TOO
            #
            while self.proceed is False:
                time.sleep(self.proceed_timeout)

            self.reading = True
            try:
                ParseUCOSched.parseTOO(too_sheetns=self.too,outfn='googledex.dat',outdir=os.getcwd(),prilim=self.prilim,certificate=self.certificate)
            except Exception as e:
                apflog("Error: Cannot download %s: %s" % (self.too,e),level="error")
            self.reading = False
            
        return
