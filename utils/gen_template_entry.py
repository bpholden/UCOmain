from __future__ import print_function
import sys
sys.path.append("../master")
from optparse import OptionParser
from datetime import datetime
import time

import numpy as np

import ParseGoogledex
import ExposureCalculations
import SchedulerConsts as sc
import UCSCScheduler_V2 as ds

def find_bstar(ras,decs,ra,dec):


    distances = np.sqrt((ra - ras)**2 + (dec - decs)**2)

    min_ind = distances.argmin()
    
    return min_ind

def make_row(star_table,ind,bstar=False):

    row = []

    row.append(star_table[ind, sc.DS_RA])
    row.append( star_table[ind, sc.DS_DEC])
    row.append(star_table[ind, sc.DS_PMRA])
    row.append(star_table[ind, sc.DS_PMDEC])
    row.append(star_table[ind, sc.DS_VMAG])
    row.append(1200)
    row.append(1e9)
    row.append( star_table[ind, sc.DS_APFPRI])
    row.append(0)
    if bstar:
        row.append(2)
    else:
        if star_table[ind, sc.DS_VMAG] > 10:
            row.append(9)
        elif star_table[ind, ds.DS_VMAG] < 8:
            row.append(5)
        else:
            row.append(7)                
    return row


if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()    
    if len(args) < 1:
        print ("needs a name")
        sys.exit()

    fp = open("tonight","a+")
    allnames, star_table, flags, stars  = ParseGoogledex.parseGoogledex()
    bstars = np.array([ True if 'HR' in n else False for n in allnames ], dtype=bool)
    npallnames = np.asarray(allnames)
    dt = datetime.utcfromtimestamp(int(time.time()))
    bline = ""
    for arg in args:
        if arg in allnames:
            i = allnames.index(arg)
            row = make_row(star_table,i)

            bstari = find_bstar(star_table[:,sc.DS_RA][bstars],star_table[:,sc.DS_DEC][bstars],star_table[i,sc.DS_RA],star_table[i,sc.DS_DEC])
            bstarrow = make_row(star_table[bstars],bstari,bstar=True)
            
            line = ds.makeScriptobsLine(allnames[i],row,flags['do'][i],dt,decker="N",I2="N",owner=flags['owner'][i])
            bline = ds.makeScriptobsLine(npallnames[bstars][bstari],bstarrow,'Y',dt,decker="N",I2="Y",owner='public')            
            fp.write(bline + "\n")
            oline ="%s #  %s\n" % (line,"pri = %s" % (star_table[i, sc.DS_APFPRI])) 
            fp.write(oline)

    fp.write(bline + "\n")
