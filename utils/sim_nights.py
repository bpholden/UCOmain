from __future__ import print_function

import argparse
from datetime import datetime, timedelta
import random
import os
import shutil
import sys

import ephem
import numpy as np
import astropy

sys.path.insert(1,"../Main")

import NightSim
import UCOScheduler as ds
import ParseUCOSched
import UCOTargets


def read_datefile(datefn):

    try:
        datefile = open(datefn)
    except Exception as e:
        print ("cannot open %s, and we died trying: %s" % (datefn, e))
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
            ofp = open(ofn,"a+", encoding='utf-8')
        else:
            ofp = open(ofn,"w", encoding='utf-8')
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


def prep_simout(outdir,simoutname):
    star_strs = dict()
    star_dates = dict()
    simoutname = os.path.join(outdir,simoutname)
    if os.path.exists(simoutname):
        try:
            simoutfp = open(simoutname, encoding='utf-8')
        except Exception as e:
            print ("Cannot open file %s for output, %s,  exiting" % (simoutname, e))
            sys.exit()

        for ln in simoutfp:
            sim_results(ln,star_strs,star_dates)
        simoutfp.close()

    try:
        simoutfp = open(simoutname,"a+", encoding='utf-8')
    except Exception as e:
        print ("Cannot open file %s for output, %s,  exiting" % (simoutname, e))
        sys.exit()

    return simoutfp, star_strs, star_dates


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile",dest="infile",default="googledex.dat",
                        help="input googledex file")
    parser.add_argument("-f","--file",dest="datefile",default="",
                        help="input date file")
    parser.add_argument("--seed",dest="seed",default=None,help="random seed")
    parser.add_argument("-b","--bstar",dest="bstar",default=True,
                        action="store_false",help="disable bstar")
    parser.add_argument("-o","--outdir",dest="outdir",default=".",
                        help="output directory")
    parser.add_argument("--rank_table",dest="rank_sheet",
                        default="2026A_rank",help="rank table file")
    parser.add_argument("--tleftfile",dest="time_left",
                        default="time_left.csv",help="time left file")
    parser.add_argument("start_date", nargs='?', default=None, help="start date for simulation (YYYY/MM/DD)")
    parser.add_argument("end_date", nargs='?', default=None, help="end date for simulation (YYYY/MM/DD)")
    parser.add_argument("-s","--simout",dest="simout",
                        default="all_dates.simout",help="simulation output file")
    options = parser.parse_args()

    if options.datefile != "":
        df = os.path.join(options.outdir,options.datefile)
        datelist = read_datefile(df)
    else:
        datelist = gen_datelist(options.start_date, options.end_date)

    if options.seed:
        random.seed(int(options.seed))
        np.random.seed(int(options.seed))

    if not os.path.isdir(options.outdir):
        try:
            os.mkdir(options.outdir)
        except Exception as e:
            print ("cannot make output directory: %s - %s" % (options.outdir, e))
            sys.exit()

    gd = os.path.join(options.outdir,options.infile)
    if not os.path.exists(gd):
        try:    
            shutil.copyfile(options.infile, gd)
        except Exception as e:
            print ("cannot copy %s to %s: %s" % (options.infile,options.outdir, e))
            sys.exit()

    return options, datelist


def update_constraints(googledex_fn):

    star_table = astropy.io.ascii.read(googledex_fn)
    star_table['night_obs'] = 0
    astropy.io.ascii.write(star_table, googledex_fn, format='ecsv', overwrite=True)

    return

def update_hour_constraints(tleftfn):

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
            crun = time_left['runname'] == runname
            time_left['left'][crun] = time_left['alloc'][crun] - time_left['used'][crun]

    time_left.write(tleftfn,format='csv',overwrite=True)
    return
###

def main():

    options,datelist = parse_args()
    bstar = options.bstar
    simoutfp, star_strs, star_dates = prep_simout(options.outdir,options.simout)

    ucotargets = UCOTargets.UCOTargets(options)
    ucotargets.make_rank_table()

    for datestr in datelist:

        if os.path.exists('hour_table'):
            os.remove('hour_table')

        ucotargets.make_hour_table()
        ucotargets.make_star_table()
        stars = ParseUCOSched.gen_stars(ucotargets.star_table)

        curtime, endtime, apf_obs = NightSim.sun_times(datestr)

        fwhms = NightSim.gen_seeing()
        slowdowns = NightSim.gen_clouds()

        do_temp = True
        lastslow = 5
        lastfwhm = 15
        otfn = os.path.join(options.outdir,"observed_targets")
        ot = open(otfn,"w",encoding='utf-8')
        ot.close()
        observing = True
        while observing:

            result = ds.get_next(curtime, lastfwhm, lastslow, ucotargets,\
                                    bstar=bstar, outfn=options.infile,\
                                    do_templates=do_temp,outdir=options.outdir)
            if result:
                if bstar:
                    bstar = False

                curtime += 70./86400 # acquisition time
                (idx,) = np.where(ucotargets.star_table['name'] == result['NAME'])
                idx = idx[0]

                for _ in range(0,int(result['NEXP'])):
                    (curtime,lastfwhm,lastslow,outstr) = \
                        NightSim.compute_simulation(result,curtime, stars[idx],\
                                                    apf_obs, slowdowns, fwhms,\
                                                        result['owner'])
                    sim_results(outstr,star_strs,star_dates)
                    simoutfp.write("%s\n" % (outstr))

                ot = open(otfn,"a+", encoding='utf-8')
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
        if ucotargets.hour_constraints:
            update_hour_constraints(options.time_left)
        update_constraints(os.path.join(options.outdir, options.infile))

        if os.path.isfile(otfn):
            try:
                os.unlink(otfn)
            except:
                print ("cannot unlink %s" %(otfn))
    if simoutfp:
        simoutfp.close()



if __name__ == "__main__":
    main()


# done
