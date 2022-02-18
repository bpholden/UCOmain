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
    parser.add_argument('sheetn',type=str)
    parser.add_argument('-f', '--frac_sheet', dest='frac_sheet', default=None, help='Download a frac table as well')
    parser.add_argument('-r', '--rank_sheet', dest='rank_sheet', default=None, help='Download a rank table as well')
    parser.add_argument('-t', '--time_left', dest='time_left', default=None, help='Use an existing time_left file')
    opt = parser.parse_args()
    
    sheetns = opt.sheetn.split(",")
    outdir = "."
    outfn = "googledex.dat"
    if os.path.exists(os.path.join(outdir,outfn)):
        os.unlink(os.path.join(outdir,outfn))


    config = dict()
    config['I2'] = 'Y'
    config['decker']='W'
    config['mode']=''
    config['obsblock']=''
    config['Bstar']='N'
    config['owner']='public'
    config['inst']='levy'
    config['raoff'] = None
    config['decoff'] = None 

        
    ParseUCOSched.parseUCOSched(sheetns=sheetns,outfn=outfn,outdir=outdir,config=config)

    if opt.frac_sheet is not None:
        
        if opt.time_left is not None and os.path.exists(opt.time_left):
            try:
                hour_constraints = astropy.io.ascii.read(opt.time_left)
            except Exception as e:
                hour_constraints = None
                print("Error: Cannot read file of time left %s : %s" % (opt.time_left,e))
        else:
            hour_constraints = None

        hour_table = ds.makeHourTable(opt.frac_sheet,datetime.datetime.now(),hour_constraints=hour_constraints)
        
    if opt.rank_sheet is not None:
        rank_table = ds.makeRankTable(opt.rank_sheet)

