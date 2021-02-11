from __future__ import print_function
import sys
sys.path.append("../master")
#from ExposureCalc import *

import numpy as np
import pickle
import ephem
import optparse
from datetime import datetime
import re
import os

import UCOScheduler_V1 as ds
import ExposureCalculations as ec
import Generate_Errors as ge
import NightSim as ns
import ParseUCOSched

def compute_simulation(result,curtime,star,apf_obs,slowdowns,fwhms,outfp):
    actel,actaz = ns.compute_el(curtime,star,apf_obs)
    actslow, actfwhm = ns.rand_obs_sample(slowdowns,fwhms)
    actfwhm = ns.gen_seeing_el(actfwhm,actel)
    lastfwhm = actfwhm
    lastslow = actslow
    meterrate = ec.getEXPMeter_Rate(result['VMAG'],result['BV'],actel,actfwhm,result['DECKER'])
    meterrate *= 1 + 0.11*np.random.randn(1)
    meterrate /= actslow
    specrate = ec.getSpec_Rate(result['VMAG'],result['BV'],actel,actfwhm,result['DECKER'])
    specrate *= 1 + 0.11*np.random.randn(1)
    specrate /= actslow
    metertime = result['COUNTS'] / meterrate
    exp_time = result['EXP_TIME']
    barycentertime = curtime
    if metertime < exp_time:
        fexptime = float(metertime)
    else:
        fexptime = float(exp_time)
        
    curtime += (fexptime+40.)/86400
    barycentertime += fexptime/(2.*86400)
    totcounts = fexptime * specrate

    precision, true_error = ge.compute_real_uncertainty(totcounts,result['BV'])
    if actaz < 180:
        actel *= -1.
    outstr = "%s %s %.5f %.1f %.1f %.2f %.2f %.2f %.2f %s" %(result['NAME'] , ephem.Date(curtime), ephem.julian_date(ephem.Date(barycentertime)), fexptime, totcounts,  actel,actaz, actfwhm, actslow, result['owner'])
    print (outstr)
    outfp.write(outstr + "\n")
    return curtime, lastfwhm, lastslow



parser = optparse.OptionParser()
parser.add_option("-d","--date",dest="date",default="today")
parser.add_option("-f","--fixed",dest="fixed",default="")
parser.add_option("-g","--googledex",dest="googledex",default="RECUR_A100,2020B_A000,2020B_A001,2020B_A002,2020B_A003,2020B_A004,2020B_A005,2020B_A006,2020B_A007,2020B_A008,2020B_A009,2020B_A010")
parser.add_option("--frac_table",dest="frac_sheetn",default="2020B_frac")
parser.add_option("--rank_table",dest="rank_sheetn",default="2020B_ranks")
parser.add_option("-i","--infile",dest="infile",default="googledex.dat")
parser.add_option("-o","--outfile",dest="outfile",default=None)
parser.add_option("-b","--bstar",dest="bstar",default=True,action="store_false")
(options, args) = parser.parse_args()    
outdir = "."

if options.date == "today":
    today = datetime.utcnow()
    datestr = "%d/%02d/%02d" % (today.year,today.month,today.day)
else:
    datestr = options.date

    
if options.fixed != "":
    if not os.path.isfile(options.fixed):
        print ("%s is not a file" % (options.fixed))

if not ns.checkdate(datestr):
    print ("%s is not an acceptable date string" % (datestr))
    sys.exit()

if options.outfile == None:
    fdatestr = re.sub("\/","-",datestr)
    outfile = "%s.simout" % (fdatestr )
else:
    outfile = options.outfile
try:
    outfp = open(outfile,"w+")
except Exception as e:
    print ("cannot open file %s for output, %s,  exiting" % (outfile,e))
    sys.exit()

if os.path.exists('hour_table'):
    os.remove('hour_table')

   

rank_table = ds.makeRankTable(options.rank_sheetn)

    
hdrstr = "#starname date time mjd exptime i2counts elevation azimuth fwhm slowdown owner\n"
outfp.write(hdrstr)
        
star_table, stars  = ParseUCOSched.parseUCOSched(sheetns=options.googledex.split(","),outfn=options.infile,outdir=outdir)

fwhms = ns.gen_seeing(val=1.0) # good conditions
slowdowns = ns.gen_clouds(val=0.5) # good conditions

lastslow = 5
lastfwhm = 15
otfn = "observed_targets"
ot = open(otfn,"w")
ot.close()
observing = True
curtime, endtime, apf_obs = ns.sun_times(datestr)
bstar = options.bstar
doTemp = True
tempcount = 0
hour_table = ds.makeHourTable(options.frac_sheetn,curtime.datetime())

while observing:
    curtime = ephem.Date(curtime)

    result = ds.getNext(curtime.datetime(), lastfwhm, lastslow, bstar=bstar, outfn=options.infile,template=doTemp,sheetns=options.googledex.split(","),outdir=outdir,frac_sheet=options.frac_sheetn,rank_sheetn=options.rank_sheetn)
    if result:
        if bstar:
            bstar = False
        if result['isTemp']:
            tempcount = tempcount + 1
        if tempcount == 2:
            doTemp=False # two per night
        curtime += 70./86400 # acquisition time
        (idx,) = np.where(star_table['name'] == result['NAME'])
        idx = idx[0]
        for i in range(0,int(result['NEXP'])):
            (curtime,lastfwhm,lastslow) = compute_simulation(result,curtime,stars[idx],apf_obs,slowdowns,fwhms,outfp)
        ot = open(otfn,"a+")
        for i in range(0,len(result["SCRIPTOBS"])):
            ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
        ot.close()
    else:
        curtime += 2100./86400 # close for lack of target
        lastslow = 5
        lastfwhm = 15
    if curtime > endtime:
        observing = False
        
        
print ("sun rose")
outfp.close()
