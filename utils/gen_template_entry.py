from __future__ import print_function
import sys
sys.path.append("../Main")
from optparse import OptionParser
from datetime import datetime
import time

import numpy as np

import ExposureCalculations
import SchedulerConsts as sc
import UCOScheduler as ds
import ParseUCOSched

if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()    
    if len(args) < 1:
        print ("needs a name")
        sys.exit()

    fp = open("tonight","a+")
    star_table, stars  = ParseUCOSched.parseUCOSched(outfn="googledex.dat")
    bstars = (star_table['Bstar'] == 'Y')|(star_table['Bstar'] == 'y')

    totexptimes = (star_table['texp'] * star_table['nexp'] + 40 * (star_table['nexp']-1))

    dt = datetime.utcfromtimestamp(int(time.time()))
    apf_obs = ds.makeAPFObs(dt)

    for arg in args:
        if np.any(star_table['name'] == arg):
            idx = np.where(star_table['name'] == arg)
            idx = idx[0][0]
            stars[idx].compute(apf_obs)
            res =  ds.makeResult(stars,star_table,totexptimes,star_table['pri'],dt,idx,focval=2,bstar=False)
            
            bidx,bfinidx = ds.findBstars(star_table,idx,bstars)
            bline = ds.makeScriptobsLine(star_table[bstars][bidx],dt,decker="N",I2="Y", owner=res['owner'],focval=2)
            line  = ds.makeScriptobsLine(star_table[idx],dt,decker="N",I2="N", owner=res['owner'],temp=True)
            bfinline = ds.makeScriptobsLine(star_table[bstars][bfinidx],dt,decker="N",I2="Y",owner=res['owner'],focval=0)

            fp.write(bline + "\n")
            oline ="%s #  %s\n" % (line,"pri = %s" % (star_table['pri'][idx])) 
            fp.write(oline)

    fp.write(bline + "\n")
