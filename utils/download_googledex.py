from __future__ import print_function
import sys
import os
import argparse
import datetime
import astropy.io.ascii

sys.path.append("../Main")

import ParseUCOSched
import UCOScheduler as ds
import UCOTargets

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Set default options")
    parser.add_argument('-r', '--rank_table', dest='rank_table', default=None, help='Rank table, will make sheet list')
    parser.add_argument('-t', '--time_left', dest='time_left', default=None, help='Use an existing time_left file')
    opt = parser.parse_args()

    if opt.rank_sheet is None:
        print("Needs a rank table")
        sys.exit(0)

    outdir = "."
    outfn = "googledex.dat"
    if os.path.exists(os.path.join(outdir,outfn)):
        os.unlink(os.path.join(outdir,outfn))

    uco_targets = UCOTargets.UCOTargets(opt)

    sheet_list = None

    config = dict()
    config['I2'] = 'Y'
    config['decker']='W'
    config['Bstar']='N'
    config['inst']='levy'

    if opt.time_left is not None and os.path.exists(opt.time_left):
        uco_targets.make_hour_constraints()

    uco_targets.make_rank_table()
    uco_targets.make_hour_table()
    uco_targets.make_star_table()
    uco_targets.append_too_column()

