import os
import shutil
import datetime

import astropy.io.ascii
import numpy as np

import ParseUCOSched
import SchedulerConsts
import UCOScheduler as ds

try: 
    from apflog import apflog
except ImportError:
    from fake_apflog import *

class UCOTargets(object):

    def __init__(self, opt, prilim=0.5):

        self.rank_table_name = opt.rank_table
        self.time_left_name = opt.time_left
        self.debug = opt.test
        self.star_tab_name = 'googledex.dat' # historical
        self.star_tab = None
        self.rank_table = None
        self.hour_table = None
        self.too = None
        self.sheets = None
        self.hour_constraints = None
        
        self.prilim = prilim
        self.certificate = SchedulerConsts.DEFAULT_CERT
        

        if opt.test:
            self.debug = opt.test
        else:
            self.debug = False

        if self.rank_table_name is None:
            apflog("Error: no rank table provided", level='error')
            return

    def __repr__(self):
        return "<UCOTargets rank_table=%s star_tab=%s>" % (self.rank_table_name, self.star_tab_name)


    def copy_backup(self, file_name):
        '''
        Copy a backup file if it exists.
        
        '''

        old_name = file_name + ".1"
        if os.path.exists(old_name):
            shutil.copyfile(old_name, file_name)
            return True
        return False

    def append_too_column(self):
        '''
        Append a 'too' column to star table based on rank_table info.
        
        '''
        if self.star_tab is None or self.rank_table is None:
            return
        if 'too' not in self.rank_table.columns:
            return
        too_sheets =  self.rank_table['sheetn'][self.rank_table['too']]

        self.star_tab['too'] = np.zeros(len(self.star_tab), dtype=bool)
        for sn in too_sheets:
            idxs = self.star_tab['sheetn'] == sn
            self.star_tab['too'][idxs] = True

    def make_hour_constraints(self):
        '''
        Read hour constraints from file if available.
        
        '''
        if self.rank_table_name is None:
            return

        if os.path.exists(self.time_left_name):
            try:
                self.hour_constraints = astropy.io.ascii.read(self.time_left_name)
            except Exception as e:
                apflog("Error: Cannot read file of time left %s : %s" % (self.time_left_name, e))


    def make_hour_table(self):
        '''
        Make hour table from rank table and hour constraints.
        
        '''
        if self.rank_table is None:
            return

        try:
            self.hour_table = ds.make_hour_table(self.rank_table,datetime.datetime.now(),
                                            hour_constraints=self.hour_constraints)
        except Exception as e:
            apflog("Error: Cannot make hour_table?! %s" % (e),level="error")


    def make_rank_table(self):
        '''
        Get the rank table google sheet if available.
        if not, try to use a backup copy.

        '''

        try:
            self.rank_table = ds.make_rank_table(self.rank_table_name, outdir=os.getcwd(),\
                               hour_constraints=self.hour_constraints)
        except Exception as e:
            apflog("Error: Cannot download rank_table?! %s" % (e),level="error")

        if self.rank_table is None:
            # goto backup
            if self.copy_backup(self.rank_table_name):
                try:
                    self.rank_table = ds.make_rank_table(self.rank_table_name,
                                                    outdir=os.getcwd(),
                                                    hour_constraints=self.hour_constraints)
                except Exception as e:
                    apflog("Error: Cannot reuse rank_table?! %s" % (e),level="error")


    def make_star_table(self):
        '''
        Get the star table google sheet if available.
        if not, try to use a backup copy.

        '''

        try:
            self.star_tab, _ = ParseUCOSched.parse_UCOSched(self.rank_table, outfn=self.star_tab_name,
                                                        outdir=os.getcwd(),
                                                        prilim=self.prilim,
                                                        certificate=self.certificate)
        except Exception as e:
            apflog("Error: Cannot download googledex?! %s" % (e),level="error")

        if self.star_tab is None:
            # goto backup
            if self.copy_backup(self.star_tab_name):
                try:
                    self.star_tab, _ = ParseUCOSched.parse_UCOSched(self.rank_table, outfn=self.star_tab_name,
                                                        outdir=os.getcwd(),
                                                        prilim=self.prilim,
                                                        certificate=self.certificate)
                except Exception as e:
                    apflog("Error: Cannot reuse googledex?! %s" % (e),level="error" )


def main():
    class Opts:
        def __init__(self):
            self.rank_table = '2025B_ranks_operational'
            self.time_left = '/home/holden/time_left.csv'
            self.test = True
    opt = Opts()
    uco_targets = UCOTargets(opt)
    uco_targets.make_hour_constraints()
    print("Hour constraints:", uco_targets.hour_constraints)
    uco_targets.make_rank_table()
    print("Rank table:", uco_targets.rank_table)
    uco_targets.make_hour_table()
    print("Hour table:", uco_targets.hour_table)
    uco_targets.make_star_table()
    print("Star table:", uco_targets.star_tab[0])
    uco_targets.append_too_column()
    print("Star table:", uco_targets.star_tab[0])

if __name__ == "__main__":
    main()
