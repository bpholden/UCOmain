#!/usr/bin/env  /opt/kroot/bin/kpython
from __future__ import print_function

import datetime
import time
import numpy as np
import argparse

import sys
sys.path.append("../Main")
import UCOScheduler as ds
import ParseUCOSched

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Set default options")
    parser.add_argument('desired_stars',type=str)
    parser.add_argument('-s','--sheetns',dest='sheetn',default="RECUR_A100,")

    options = parser.parse_args()    

    outfn = "googledex.dat"
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

    
    desiredstars = options.desired_stars.split(",")
    sheetns = options.sheetn.split(",")

    star_table, stars = ParseUCOSched.parse_UCOSched(sheetns=sheetns,outfn=outfn,config=config)
    now = datetime.datetime.now()
    
    for star in desiredstars:
        aidx, = np.where(star_table["name"] == star)
        idx = aidx[0]
        ret = ds.make_scriptobs_line(star_table[idx],now,decker=star_table['decker'][idx], \
                                                owner=star_table['sheetn'][idx], \
                                                I2=star_table['I2'][idx])
        print(ret)
