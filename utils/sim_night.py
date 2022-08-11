from __future__ import print_function
import sys
sys.path.append("../Main")
#from ExposureCalc import *
import optparse
from datetime import datetime
import re
import os

import numpy as np
import ephem
import astropy.io.ascii

import UCOScheduler as ds

import NightSim 
import ParseUCOSched




parser = optparse.OptionParser()
parser.add_option("-d","--date",dest="date",default="today")
parser.add_option("-f","--fixed",dest="fixed",default="")
parser.add_option("--rank_table",dest="rank_sheetn",default="2022B_ranks")
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

if not NightSim.checkdate(datestr):
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

tleftfn = 'time_left.csv'
if os.path.exists(tleftfn):
    hour_constraints = astropy.io.ascii.read(tleftfn)
else:
    hour_constraints = None
   

rank_table = ds.makeRankTable(options.rank_sheetn)
sheet_list = list(rank_table['sheetn'][rank_table['rank'] > 0])
    
hdrstr = "#starname date time mjd exptime i2counts elevation azimuth fwhm slowdown owner\n"
outfp.write(hdrstr)
        
star_table, stars  = ParseUCOSched.parseUCOSched(sheetns=sheet_list,outfn=options.infile,outdir=outdir,hour_constraints=hour_constraints)

fwhms = NightSim.gen_seeing(val=1.0) # good conditions
slowdowns = NightSim.gen_clouds(val=.6) # typical conditions

lastslow = 5
lastfwhm = 15
otfn = "observed_targets"
ot = open(otfn,"w")
ot.close()
observing = True
curtime, endtime, apf_obs = NightSim.sun_times(datestr)
bstar = options.bstar
doTemp = True
tempcount = 0

hour_table = ds.makeHourTable(rank_table,curtime.datetime(),hour_constraints=hour_constraints)

while observing:
    curtime = ephem.Date(curtime)

    result = ds.getNext(curtime.datetime(), lastfwhm, lastslow, bstar=bstar, outfn=options.infile,template=doTemp,sheetns=sheet_list,outdir=outdir,rank_sheetn=options.rank_sheetn)
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
            (curtime,lastfwhm,lastslow,outstr) = NightSim.compute_simulation(result,curtime,stars[idx],apf_obs,slowdowns,fwhms,result['owner'])
            outfp.write(outstr + "\n")
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
