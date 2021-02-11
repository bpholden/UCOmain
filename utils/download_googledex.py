from __future__ import print_function
import sys
import os
import argparse
import datetime
sys.path.append("../master")

import ParseUCOSched
import UCOScheduler_V1 as ds

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Set default options")
    parser.add_argument('sheetn',type=str)
    parser.add_argument('-f', '--frac_sheet', dest='frac_sheet', default=None, help='Download a frac table as well')
    parser.add_argument('-r', '--rank_sheet', dest='rank_sheet', default=None, help='Download a rank table as well')
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

