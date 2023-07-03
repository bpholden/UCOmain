from __future__ import print_function
import sys
import os
import argparse
import datetime
import astropy.io.ascii

sys.path.append("../Main")

import ParseUCOSched
import UCOScheduler as ds

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Set default options")
    parser.add_argument('-s', '--sheet_list', dest='sheetn', default=None, help='List of sheets')
    parser.add_argument('-r', '--rank_sheet', dest='rank_sheet', default=None, help='Rank table, will make sheet list')
    parser.add_argument('-t', '--time_left', dest='time_left', default=None, help='Use an existing time_left file')
    opt = parser.parse_args()

    if opt.sheetn is None and opt.rank_sheet is None:
        print("Need either a list of sheets or a rank table")
        sys.exit(0)

    outdir = "."
    outfn = "googledex.dat"
    if os.path.exists(os.path.join(outdir,outfn)):
        os.unlink(os.path.join(outdir,outfn))


    config = dict()
    config['I2'] = 'Y'
    config['decker']='W'
    config['Bstar']='N'
    config['inst']='levy'

    if opt.time_left is not None and os.path.exists(opt.time_left):
        try:
            hour_constraints = astropy.io.ascii.read(opt.time_left)
        except Exception as e:
            hour_constraints = None
            print("Error: Cannot read file of time left %s : %s" % (opt.time_left,e))
    else:
        hour_constraints = None

    if opt.rank_sheet is not None:
        rank_table = ds.make_rank_table(opt.rank_sheet,hour_constraints=hour_constraints)
        hour_table = ds.make_hour_table(rank_table,datetime.datetime.now(),hour_constraints=hour_constraints)
        sheet_list = list(rank_table['sheetn'][rank_table['rank'] > 0])

    if opt.sheetn is not None:
        sheet_list = opt.sheetn.split(",")

    if sheet_list:
        ParseUCOSched.parse_UCOSched(sheetns=sheet_list,outfn=outfn,outdir=outdir,config=config)
