import os
import shutil
import datetime

import astropy.io.ascii
import numpy as np

import ParseUCOSched
import SchedulerConsts

try: 
    from apflog import apflog
except ImportError:
    from fake_apflog import *

class UCOTargetTables(object):
    '''
    Class to handle UCO target tables: rank table, hour table, and star table.
    
    '''
    def __init__(self, opt, prilim=0.5):

        self.rank_table_name = opt.rank_table
        self.time_left_name = opt.time_left
        self.debug = opt.test if hasattr(opt, 'test') else False
        self.star_table_name = 'googledex.dat' # historical
        self.star_table = None
        self.stars = None
        self.rank_table = None
        self.rank_table_filename = "rank_table"
        self.hour_table = None
        self.hour_table_filename = "hour_table"
        self.stars = None
        self.halve = opt.halve if hasattr(opt, 'halve') else False
        self.too = None
        self.sheets = None
        self.too_sheets = None
        self.hour_constraints = None

        self.prilim = prilim
        self.certificate = SchedulerConsts.DEFAULT_CERT

        if self.rank_table_name is None:
            apflog("Error: no rank table provided", level='error')
            return

    def __repr__(self):
        return "<UCOTargetTables rank_table=%s star_table=%s>" % \
            (self.rank_table_name, self.star_table_name)


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
        if self.star_table is None or self.rank_table is None:
            return
        if 'too' not in self.rank_table.columns:
            return
        if 'too' in self.star_table.columns:
            return
        too_sheets =  self.rank_table['sheetn'][self.rank_table['too']]
        self.too_sheets = list(too_sheets)

        self.star_table['too'] = np.zeros(len(self.star_table), dtype=bool)
        for sn in too_sheets:
            idxs = self.star_table['sheetn'] == sn
            self.star_table['too'][idxs] = True

    def make_hour_constraints(self):
        '''
        Read hour constraints from file if available.
        
        '''
        if self.rank_table_name is None:
            return

        if self.time_left_name is None:
            return

        if os.path.exists(self.time_left_name):
            try:
                self.hour_constraints = astropy.io.ascii.read(self.time_left_name)
            except Exception as e:
                apflog("Error: Cannot read file of time left %s : %s" % (self.time_left_name, e))


    def make_hour_table(self, obs_datetime=None):
        '''
        Make hour table from rank table and hour constraints.
        
        '''

        req_datetime = obs_datetime if obs_datetime is not None else datetime.datetime.utcnow()

        self.make_hour_constraints()

        if self.rank_table_name is None:
            return

        if self.rank_table is None:
            self.make_rank_table()

        try:
            self.hour_table = ParseUCOSched.make_hour_table(self.rank_table, req_datetime,
                                            hour_constraints=self.hour_constraints)
        except Exception as e:
            apflog("Error: Cannot make hour_table?! %s" % (e),level="error")


    def make_rank_table(self):
        '''
        Get the rank table google sheet if available.
        if not, try to use a backup copy.

        '''

        try:
            self.rank_table = ParseUCOSched.make_rank_table(self.rank_table_name, \
                                outdir=os.getcwd(), outfn=self.rank_table_filename, \
                                hour_constraints=self.hour_constraints, halve_rank=self.halve)
        except Exception as e:
            apflog("Error: Cannot download rank_table?! %s" % (e),level="error")

        if self.rank_table is None:
            # goto backup
            if self.copy_backup(self.rank_table_name):
                try:
                    self.rank_table = ParseUCOSched.make_rank_table(self.rank_table_name,
                                                    outdir=os.getcwd(),
                                                    hour_constraints=self.hour_constraints)
                except Exception as e:
                    apflog("Error: Cannot reuse rank_table?! %s" % (e),level="error")

        if self.rank_table is not None:
            self.sheets = list(self.rank_table['sheetn'][self.rank_table['rank'] > 0])


    def make_star_table(self):
        '''
        Get the star table google sheet if available.
        if not, try to use a backup copy.

        '''

        try:
            self.star_table, _ = ParseUCOSched.parse_UCOSched(self.rank_table, 
                                                        outfn=self.star_table_name,
                                                        outdir=os.getcwd(),
                                                        prilim=self.prilim,
                                                        certificate=self.certificate)
        except Exception as e:
            apflog("Error: Cannot download googledex?! %s" % (e),level="error")

        if self.star_table is None:
            # goto backup
            if self.copy_backup(self.star_table_name):
                try:
                    self.star_table, _ = ParseUCOSched.parse_UCOSched(self.rank_table, 
                                                        outfn=self.star_table_name,
                                                        outdir=os.getcwd(),
                                                        prilim=self.prilim,
                                                        certificate=self.certificate)
                except Exception as e:
                    apflog("Error: Cannot reuse googledex?! %s" % (e),level="error" )
        self.append_too_column()

    def check_files(self, outfn=None):
        '''
        If the star table file is missing, restore it from its backup.

        '''
        if outfn is None:
            outfn = self.star_table_name
        outdir = os.getcwd()
        fullpath = os.path.join(outdir, outfn)
        if os.path.isfile(fullpath):
            return

        # make it so
        backup = fullpath + ".1"
        try:
            shutil.copyfile(backup, fullpath)
        except Exception as e:
            err_str = "Cannot copy %s to %s: %s %s" % (backup, fullpath, type(e), e)
            apflog(err_str, echo=True, level='error')

    def update_from_observed(self, ptime, outfn=None, toofn='too.dat'):
        '''
        Update the local star table file with the observed log,
        set star_table from it, and return the ObservedLog.

        '''
        if outfn is None:
            outfn = self.star_table_name
        observed, self.star_table = ParseUCOSched.update_local_starlist(ptime,\
                                                               outfn=outfn, toofn=toofn, \
                                                                observed_file="observed_targets")
        return observed

    def gen_stars(self):
        '''
        Make the pyephem objects for the current star table.

        '''
        self.stars = ParseUCOSched.gen_stars(self.star_table)
        return self.stars

    def update_hour_table(self, observed, dt, outfn='hour_table', outdir=None):
        '''
        update_hour_table(observed, dt, outfn='hour_table', outdir=None)

        Updates hour_table with history of observations.
        observed is the observed log
        dt is the current datetime
        outfn is the output filename, defaults to hour_table
        outdir is the output directory, defaults to current working directory

        '''
        if self.hour_table is None:
            return None

        if not outdir :
            outdir = os.getcwd()

        outfn = os.path.join(outdir, outfn)

        hours = dict()

        # observed objects have lists as attributes
        # reverse time order, so most recent target observed is first.

        observed.reverse()

        nobj = len(observed.names)
        for i in range(0,nobj):
            own = observed.owners[i]
            if own not in list(hours):
                hours[own] = 0.0

        cur = dt
        for i in range(0,nobj):
            hr, mn = observed.times[i]
            prev = datetime.datetime(dt.year, dt.month, dt.day, hr, mn)
            diff = cur - prev
            hourdiff = diff.days * 24 + diff.seconds / 3600.
            if hourdiff > 0:
                hours[observed.owners[i]] += hourdiff
                cur = prev

        for ky in list(hours.keys()):
            if ky == 'public':
                self.hour_table['cur'][self.hour_table['sheetn'] == 'RECUR_A100'] = hours[ky]
            else:
                self.hour_table['cur'][self.hour_table['sheetn'] == ky] = hours[ky]

        try:
            self.hour_table.write(outfn,format='ascii',overwrite=True)
        except Exception as e:
            apflog("Cannot write table %s: %s %s" % (outfn, type(e), e), level='error', echo=True)

        observed.reverse()

        return self.hour_table

def main():
    class Opts:
        def __init__(self):
            self.rank_table = '2026B_ranks'
            self.time_left = '/home/holden/time_left.csv'
            self.test = True
            self.halve = True
    opt = Opts()
    uco_targets = UCOTargetTables(opt)
    uco_targets.make_hour_table()
    print("Hour table:", uco_targets.hour_table)
    uco_targets.make_star_table()
    print("Star table:", uco_targets.star_table[0])

if __name__ == "__main__":
    main()
