from __future__ import print_function
import datetime
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
import SchedulerConsts

class getUCOTargets(threading.Thread):

    def __init__(self, opt, task='master', prilim=0.5, wait_time=0):
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
        self.certificate = SchedulerConsts.DEFAULT_CERT
        self.star_tab = 'googledex.dat' # historical

        if opt.test:
            self.debug = opt.test
        else:
            self.debug = False

        if self.rank_table is None and self.sheets is None:
            apflog("Error: no rank table and no sheet list provided", level='error')
            return

        self.start()

    def copy_backup(self, file_name):
        old_name = file_name + ".1"
        if os.path.exists(old_name):
            shutil.copyfile(old_name, file_name)
            return True
        return False

    def append_too_column(self, tab, rank_table):
        if tab is None or rank_table is None:
            return
        if 'too' not in rank_table.columns:
            return
        too_sheets =  rank_table['sheetn'][rank_table['too']]

        tab['too'] = np.zeros(len(tab), dtype=bool)
        for sn in too_sheets:
            idxs = tab['sheetn'] == sn
            tab['too'][idxs] = True

    def run(self):

        expression = '$eostele.SUNEL < -0.0'
        if self.debug:
            expression = '$eostele.SUNEL < 100.0'
        APFTask.waitFor(self.task, True, expression=expression, timeout=self.wait_time)

        if self.signal is False:
            return
        
        if self.rank_table is None and self.sheets is None:
            return

        hour_constraints = None
        if os.path.exists(self.time_left):
            try:
                hour_constraints = astropy.io.ascii.read(self.time_left)
            except Exception as e:
                hour_constraints = None
                apflog("Error: Cannot read file of time left %s : %s" % (self.time_left,e))

        rank_table = None
        try:
            rank_table = ds.make_rank_table(self.rank_table, outdir=os.getcwd(),\
                                            hour_constraints=hour_constraints)
        except Exception as e:
            apflog("Error: Cannot download rank_table?! %s" % (e),level="error")

        if rank_table is None:
        # goto backup
            if self.copy_backup(self.rank_table):
                try:
                    rank_table = ds.make_rank_table(self.rank_table,
                                                    outdir=os.getcwd(),
                                                    hour_constraints=hour_constraints)
                except Exception as e:
                    apflog("Error: Cannot reuse rank_table?! %s" % (e),level="error")
                    rank_table = None


        if self.signal is False:
            return

        if self.sheets is None and rank_table:
            self.sheets = list(rank_table['sheetn'][rank_table['rank'] > 0])

        sheetlist_name = 'MASTER_SHEETLIST'
        if self.debug:
            sheetlist_name = 'EXAMPLE_VAR_1'

        try:
            self.apftask.write(sheetlist_name,",".join(self.sheets),timeout=2)
        except Exception as e:
            apflog("Cannot write apftask.%s: %s" % (sheetlist_name, e), level='warn',echo=True)

        if rank_table and 'too' in rank_table.columns:
            if np.any(rank_table['too']):
                self.too = list(rank_table['sheetn'][rank_table['too']])

        if self.signal is False:
            return

        tab = None
        try:
            tab, _ = ParseUCOSched.parse_UCOSched(rank_table, outfn=self.star_tab,
                                                        outdir=os.getcwd(),
                                                        prilim=self.prilim,
                                                        certificate=self.certificate)
        except Exception as e:
            apflog("Error: Cannot download googledex?! %s %s" % (type(e), e),level="error")

        if tab is None:
            self.copy_backup(self.star_tab)

        try:
            _ = ds.make_hour_table(rank_table,datetime.datetime.now(),
                                            hour_constraints=hour_constraints)
        except Exception as e:
            apflog("Error: Cannot make hour_table?! %s" % (e),level="error")

        while self.signal and self.too is not None and not self.debug:

            if APFTask.waitfor(self.task,False,expression='apftask.SCRIPTOBS_PHASE==Observing',
                               timeout=self.timeout):

                self.reading = True
                try:
                    ParseUCOSched.parse_TOO(too_sheetns=self.too, outfn=self.star_tab,
                                            outdir=os.getcwd(), prilim=self.prilim,
                                            certificate=self.certificate)
                except Exception as e:
                    apflog("Error: Cannot download %s: %s" % (self.too,e),level="error")
                self.reading = False

        return


if __name__ == "__main__":

    class Opt:
        pass


    task = 'example'
    APFTask.establish(task, os.getpid())
    opt = Opt()
    opt.test = True
    opt.time_left = "/home/holden/time_left.csv"
    opt.rank_table = '2025B_ranks_operational'

    get_targs = getUCOTargets(opt, task=task)
