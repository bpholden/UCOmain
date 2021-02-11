#!/usr/bin/env  /opt/kroot/bin/kpython
from __future__ import print_function
import sys
sys.path.append("../master")
#from ExposureCalc import *
from optparse import OptionParser
from datetime import datetime
import time

import numpy as np

import ParseGoogledex
import ExposureCalculations
import SchedulerConsts as sc
import UCSCScheduler_V2 as ds

if __name__ == "__main__":

    THRESHOLD = 60 * 60
    parser = OptionParser()
    parser.add_option("-s","--slowdown",dest="slowdown",default=0.4,type="float")
    parser.add_option("-f","--fwhm",dest="fwhm",default=13,type="float")
    parser.add_option("-e","--el",dest="el",default=70,type="float")
    (options, args) = parser.parse_args()    

    allnames, star_table, flags, stars  = ParseGoogledex.parseGoogledex()

    el = np.zeros_like(star_table[:, sc.DS_BV])
    el += options.el
    fwhm = np.array(options.fwhm)
    i2counts = np.zeros_like(star_table[:, sc.DS_BV])
    totexptimes = np.zeros_like(star_table[:, sc.DS_BV])

    precision = star_table[:, sc.DS_ERR]
    i2counts = ExposureCalculations.getI2_K(precision)
    mstars_inds = np.where(star_table[:, ds.DS_BV] > 1.2)
    i2counts[mstars_inds] = ExposureCalculations.getI2_M(precision[mstars_inds])

    exp_times, exp_counts, i2cnts = ds.calculateUCSCExposureTime(star_table[:, sc.DS_VMAG],precision,el,fwhm,star_table[:, sc.DS_BV],flags['decker'])
    exp_times *= options.slowdown
    totexptimes += ds.computeMaxTimes(exp_times,star_table[:, sc.DS_MAX])

    mxtime = np.zeros_like(star_table[:,sc.DS_MAX])
    mxtime += sc.MAX_EXPTIME
    shorter = (star_table[:,sc.DS_MAX] < sc.MAX_EXPTIME) & (star_table[:,sc.DS_MAX] >0)
    mxtime[shorter] = star_table[:,sc.DS_MAX][shorter]

    star_table[:, sc.DS_EXPT], exps = ds.format_time(totexptimes,i2counts,star_table[:, sc.DS_NSHOTS],star_table[:, sc.DS_MIN],mxtime)
    star_table[:, sc.DS_COUNTS], star_table[:, sc.DS_NSHOTS] = ds.format_expmeter(exp_counts,exps,totexptimes)

    fin_pre = precision
   
    dt = datetime.utcfromtimestamp(int(time.time()))
    for i in range(len(stars)):
        if star_table[i, sc.DS_APFPRI] < 5:
            continue
        row = []

        row.append( star_table[i, sc.DS_RA])
        row.append( star_table[i, sc.DS_DEC])
        row.append( star_table[i, sc.DS_PMRA])
        row.append( star_table[i, sc.DS_PMDEC])
        row.append( star_table[i, sc.DS_VMAG])
        row.append( star_table[i, sc.DS_EXPT])
        row.append( star_table[i, sc.DS_COUNTS])
        row.append( star_table[i, sc.DS_APFPRI])
        row.append(0)
        row.append( star_table[i, sc.DS_NSHOTS])
        line = ds.makeScriptobsLine(allnames[i],row,flags['do'][i],dt,decker=flags['decker'][i],I2=flags['I2'][i],owner=flags['owner'][i])
        print ("%s #  %s" % (line,"pri = %s" % (star_table[i, sc.DS_APFPRI])))

