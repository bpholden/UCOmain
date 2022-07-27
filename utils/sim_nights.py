from __future__ import print_function
import sys
sys.path.insert(1,"../Main")

import pickle
import optparse
from datetime import datetime, timedelta
import random
import re
import os
import shutil

import ephem
import numpy as np
import astropy

import NightSim
import UCOScheduler as ds
import ParseUCOSched



def read_datefile(datefn):

    try:
        datefile = open(options.datefile)
    except Exception as e:
        print ("cannot open %s, and we died trying: %s" % (args[0], e))
        sys.exit()

    datelist = []
    for line in datefile:
        datestr, = line.split()
        if not NightSim.checkdate(datestr):
            print ("%s is not an acceptable date string" % (datestr))
            sys.exit()
        datelist.append(datestr)
        
    return datelist

def gen_datelist(startstr,endstr):
    datelist = []
    start = datetime.strptime(startstr,"%Y/%m/%d")
    end  = datetime.strptime(endstr,"%Y/%m/%d")

    if end < start:
        print("Start date %s is after end date %s"  %( str(start), str(end)))
        sys.exit(1)

    
    cur = start
    while cur < end:
        breakbeg = datetime(cur.year,12,24)
        breakend = datetime(cur.year+1,1,1)
        td = timedelta(1)
        cur = cur + td
        if cur < end and (cur < breakbeg or cur > breakend):
            datestr = "%d/%02d/%02d" % (cur.year,cur.month,cur.day)
            datelist.append(datestr)
                        
    return datelist


def write_sim_file_results(star_strs,outdir="../SimFiles"):
    outpath = os.path.join(os.curdir,outdir)
    for starname in star_strs.keys():
        ofn = "%s.sim" % (starname)
        ofn = os.path.join(outpath,ofn)
        if os.path.exists(ofn):
            ofp = open(ofn,"a+")
        else:
            ofp = open(ofn,"w")
        for outstr in star_strs[starname]:
            (name,date,time,mjd,expt,i2cnts,actel,actaz,fwhm,slow,owner) = outstr.split()
            outstr = "%s %s %s\n" % (mjd,i2cnts,actel)
            ofp.write(outstr)
        ofp.close()
    

def sim_results(outstr,star_strs,star_dates):
    (name,date,time,mjd,expt,i2cnts,actel,actaz,fwhm,slow,owner) = outstr.split()
    try:
        mjd = float(mjd)
        if name in star_strs.keys():
            star_strs[name].append(outstr)
            star_dates[name].append( mjd )
        else:
            star_strs[name] = [outstr]
            star_dates[name] = [ mjd ]
    except:
        pass
    return


def prep_master(outdir,mastername):
    star_strs = dict()
    star_dates = dict()
    mastername = os.path.join(outdir,mastername)
    newmaster = False
    if not os.path.exists(mastername):
        newmaster = True
    else:
        try:
            masterfp = open(mastername)
        except Exception as e:
            print ("Cannot open file %s for output, %s,  exiting" % (mastername, e))
            sys.exit()
        
        for ln in masterfp:
            sim_results(ln,star_strs,star_dates)
        masterfp.close()
    
    try:
        masterfp = open(mastername,"a+")
    except Exception as e:
        print ("Cannot open file %s for output, %s,  exiting" % (mastername, e))
        sys.exit()

    return masterfp,star_strs, star_dates


def parse_args():
    parser = optparse.OptionParser()
    parser.add_option("-i","--infile",dest="infile",default="googledex.dat")
    parser.add_option("-f","--file",dest="datefile",default="")
    parser.add_option("--seed",dest="seed",default=None)
    parser.add_option("-b","--bstar",dest="bstar",default=True,action="store_false")
    parser.add_option("-o","--outdir",dest="outdir",default=".")        
    parser.add_option("--rank_table",dest="rank_sheetn",default="2022B_ranks")

    parser.add_option("-m","--masterfile",dest="master",default="sim_master.simout")
    (options, args) = parser.parse_args()    

    if len(args) < 2 and options.datefile == "":
        print ("needs either a date pair or an input file which lists the dates")
        sys.exit()    

    if options.datefile != "":
        df = os.path.join(options.outdir,options.datefile)
        datelist = read_datefile(df)
    else:
        datelist = gen_datelist(args[0],args[1])

    if options.seed:
        random.seed(int(options.seed))
        np.random.seed(int(options.seed))

    if not os.path.isdir(options.outdir):
        try:
            os.mkdir(options.outdir)
        except Exception as e:
            print ("cannot make output directory: %s - %s" % (options.outdir,e))
            sys.exit()

    gd = os.path.join(options.outdir,options.infile)
    if not os.path.exists(gd):
        try:    
            shutil.copyfile(options.infile, gd)
        except Exception as e:
            print ("cannot copy %s to %s: %s" % (options.infile,options.outdir,e))
            sys.exit()
            
    return options, datelist


def updateHourConstraints(tleftfn):

    hour_table =  astropy.io.ascii.read('hour_table')
    time_left = astropy.io.ascii.read(tleftfn)

    for sheetn in hour_table['sheetn']:
        if sheetn == 'RECUR_A100':
            runname = 'public'
        else:
            runname = sheetn
            
        used = hour_table['cur'][hour_table['sheetn'] == sheetn]
        time_left['used'][time_left['runname'] == runname] += used
        if runname != 'public':
            time_left['left'][time_left['runname'] == runname] = time_left['alloc'][time_left['runname'] == runname] - time_left['used'][time_left['runname'] == runname]

    time_left.write(tleftfn,format='csv',overwrite=True)
    
###

if __name__ == "__main__":


    options,datelist = parse_args()
    bstar = options.bstar
    masterfp,star_strs,star_dates = prep_master(options.outdir,options.master)

    rank_table = ds.makeRankTable(options.rank_sheetn)
    sheetns = list(rank_table['sheetn'][rank_table['rank'] > 0])
    
    for datestr in datelist:

        if os.path.exists('hour_table'):
            os.remove('hour_table')
            
        tleftfn = 'time_left.csv'
        if os.path.exists(tleftfn):
            hour_constraints = astropy.io.ascii.read(tleftfn)
        else:
            hour_constraints = None
             
        curtime, endtime, apf_obs = NightSim.sun_times(datestr)
        
        hour_table = ds.makeHourTable(rank_table,curtime.datetime(),hour_constraints=hour_constraints)
        
        star_table, stars  = ParseUCOSched.parseUCOSched(sheetns=sheetns,outfn=options.infile,outdir=options.outdir)
    
        fwhms = NightSim.gen_seeing()
        slowdowns = NightSim.gen_clouds()

        doTemp = True
        lastslow = 5
        lastfwhm = 15
        otfn = os.path.join(options.outdir,"observed_targets")
        ot = open(otfn,"w")
        ot.close()
        observing = True
        while observing:

            result = ds.getNext(curtime, lastfwhm, lastslow, bstar=bstar, outfn=options.infile,template=doTemp,
                                    sheetns=sheetns,outdir=options.outdir,rank_sheetn=options.rank_sheetn)
            if result:
                if bstar:
                    bstar = False
                
                curtime += 70./86400 # acquisition time
                (idx,) = np.where(star_table['name'] == result['NAME'])
                idx = idx[0]

                for i in range(0,int(result['NEXP'])):
                    (curtime,lastfwhm,lastslow,outstr) = NightSim.compute_simulation(result,curtime,stars[idx],apf_obs,slowdowns,fwhms,result['owner'])
                    sim_results(outstr,star_strs,star_dates)
                    masterfp.write("%s\n" % (outstr))
                    
                ot = open(otfn,"a+")
                ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
                ot.close()
            else:
                curtime += 2100./86400 # close for lack of target
                lastslow = 5
                lastfwhm = 15
            if curtime > endtime:
                observing = False
        
            curtime = ephem.Date(curtime)
        
        print ("sun rose")
        if hour_constraints:
            updateHourConstraints(tleftfn)

        if os.path.isfile(otfn):
            try:
                os.unlink(otfn)
            except:
                print ("cannot unlink %s" %(otfn))
    if masterfp:
        masterfp.close()
        

# done
