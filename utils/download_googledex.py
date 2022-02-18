from __future__ import print_function
import sys
import os
import argparse
import datetime
sys.path.append("../Main")

import ParseUCOSched
import UCOScheduler as ds

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Set default options")
    parser.add_argument('-s', '--sheet_list', dest='sheetn', default=None, help='List of sheets')
    parser.add_argument('-r', '--rank_sheet', dest='rank_sheet', default=None, help='Rank table, will make sheet list')
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
        hour_table = ds.makeHourTable(opt.frac_sheet,datetime.datetime.now())
    if opt.rank_sheet is not None:
        rank_table = ds.makeRankTable(opt.rank_sheet)

