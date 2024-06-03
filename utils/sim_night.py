from __future__ import print_function
import sys

import optparse
from datetime import datetime, timedelta
import re
import os

import numpy as np
import ephem
import astropy.io.ascii

sys.path.append("../Main")
import UCOScheduler as ds

import NightSim 
import ParseUCOSched

def get_start_time(hr_mn, datestr):
    '''
    get_start_time(hr, mn, datestr)

    Returns a UNIX time stamp for the start time
    '''
    datestr = "%s %s" % (datestr, hr_mn)
    dt = datetime.strptime(datestr,"%Y/%m/%d %H:%M")

    po = timedelta(days=1)

    dt = dt + po
    return float(dt.strftime("%s"))


def parse_options():
    parser = optparse.OptionParser()
    parser.add_option("-d","--date",dest="date",default="today")
    parser.add_option("-f","--fixed",dest="fixed",default="")
    parser.add_option("--rank_table",dest="rank_sheetn",default="2023A_ranks")
    parser.add_option("-i","--infile",dest="infile",default="googledex.dat")
    parser.add_option("-o","--outfile",dest="outfile",default=None)
    parser.add_option("-b","--bstar",dest="bstar",default=True,action="store_false")
    parser.add_option("-s","--start_time",dest="start_time",default=None)
    (options, args) = parser.parse_args()

    return options, args

def find_date(in_date):
    if in_date == "today":
        today = datetime.now()
        po = timedelta(days=1)
        today = today + po
        datestr = "%d/%02d/%02d" % (today.year,today.month,today.day)
    else:
        datestr = in_date
    return datestr

def find_fixed(in_fixed):
    if in_fixed != "":
        if not os.path.isfile(in_fixed):
            print ("%s is not a file" % (in_fixed))
            return None
    else:
        return in_fixed

def make_outfile(outfile, datestr):
    if outfile == None:
        fdatestr = re.sub("\/","-",datestr)
        outfile = "%s.simout" % (fdatestr )

    try:
        outfp = open(outfile,"w+")
    except Exception as e:
        print ("cannot open file %s for output, %s,  exiting" % (outfile,e))
        sys.exit()

    hdrstr = "#starname date time mjd exptime i2counts elevation azimuth fwhm slowdown owner\n"
    outfp.write(hdrstr)

    return outfp


def should_start_list(start_time, cur_time):
    """ Observe.should_start_list()
        should we start a fixed observing list or not? true if start time is None or if w/in + 1 hour - 0.5 hours of start time
    """
    if start_time is None:
        return True
    
    if cur_time > start_time and cur_time - start_time < 3600:
        return True
    if cur_time < start_time and start_time - cur_time < 600:
        return True
    return False


def main():

    options, args = parse_options()
    outdir = "."

    datestr = find_date(options.date)

    start_time = None
    if options.start_time:
        start_time = get_start_time(options.start_time, datestr)

    options.fixed = find_fixed(options.fixed)

    if not NightSim.checkdate(datestr):
        print ("%s is not an acceptable date string" % (datestr))
        sys.exit()

    outfp = make_outfile(options.outfile, datestr)
    

    if os.path.exists('hour_table'):
        os.remove('hour_table')

    tleftfn = 'time_left.csv'
    if os.path.exists(tleftfn):
        hour_constraints = astropy.io.ascii.read(tleftfn)
    else:
        hour_constraints = None
   

    rank_table = ds.make_rank_table(options.rank_sheetn)
    sheet_list = list(rank_table['sheetn'][rank_table['rank'] > 0])
        
    star_table, stars  = ParseUCOSched.parse_UCOSched(sheetns=sheet_list,outfn=options.infile,outdir=outdir,hour_constraints=hour_constraints)

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

    _ = ds.make_hour_table(rank_table,curtime.datetime(),hour_constraints=hour_constraints)

    while observing:
        curtime = ephem.Date(curtime)

        result = ds.get_next(curtime.datetime(), lastfwhm, lastslow, bstar=bstar, \
                             outfn=options.infile,template=doTemp,sheetns=sheet_list,\
                                outdir=outdir,rank_sheetn=options.rank_sheetn,\
                                    start_time=start_time)
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
        
        
    print("Updating star list with final observations")
    curtime = ephem.Date(curtime)
    _, star_table = ParseUCOSched.update_local_starlist(curtime.datetime(),outfn=options.infile,observed_file=otfn)
    print ("sun rose")
    outfp.close()

if __name__ == "__main__":
    main()