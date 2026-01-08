import sys
sys.path.append("../Main")
#from ExposureCalc import *

import numpy as np
import pickle
import ephem
import optparse
from datetime import datetime
import re
import os
import ExposureCalculations
import Generate_Errors 

def sun_times(datestr):
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = -9.0*np.pi / 180.0
    apf_obs.date = datestr
    sunset = apf_obs.next_setting(ephem.Sun())
    sunrise = apf_obs.next_rising(ephem.Sun())
    return sunset, sunrise, apf_obs
    
def make_obs_sample(fn):
    slow,fwhm = np.loadtxt(fn,unpack=True)
    return slow, fwhm

def gen_seeing(nsize=200,val=-1):
    if val < 0:
        val = np.random.uniform(size=1)
    alpha = 0.52

    if val < 0.9:
        mean = np.random.normal(loc=8.,scale=1.0,size=1)
        rms  = np.random.normal(loc=1.9,scale=1.0,size=1)
    else:
        mean = np.random.normal(loc=19.0,scale=1.0,size=1)
        rms  = np.random.normal(loc=4.5,scale=1.0,size=1)
    real_rms = np.sqrt((1-alpha**2) * rms**2)
    deviates = np.random.normal(loc=0,scale=real_rms,size=nsize)
    for i in np.arange(1,nsize):
        deviates[i] += alpha*deviates[i-1]
    deviates += mean

#    deviates = np.random.normal(loc=mean,scale=rms,size=nsize)
    return deviates

def gen_seeing_el(deviate,el):
    zd = 90 - el
    deviate += (0.0903544076597*zd +  -0.00172591889888*zd*zd + 3.3157238117e-05*zd*zd*zd)
    return deviate

def gen_clouds(nsize=200,val=-1):
    if val < 0:
        val = np.random.uniform(size=1)
    alpha = 0.353
    if val < 0.7:
        mean = np.random.normal(loc=1.,scale=0.1,size=1)
        rms  = np.random.normal(loc=0.05,scale=0.01,size=1)
    elif val > 0.9:
        mean = np.random.normal(loc=2.0,scale=0.4,size=1)
        rms  = np.random.normal(loc=0.5,scale=0.1,size=1)
    else:
        mean = np.random.normal(loc=1.7,scale=0.1,size=1)
        rms  = np.random.normal(loc=0.3,scale=0.05,size=1)

    real_rms = np.sqrt((1-alpha**2) * rms**2)
    deviates = np.random.normal(loc=0,scale=real_rms,size=nsize)
    for i in np.arange(1,nsize):
        deviates[i] += alpha*deviates[i-1]
    deviates += mean
    deviates[deviates < 0.3] = 0.3
    return deviates

def rand_obs_sample(slows,fwhms):
    ls = len(slows) -1
    lf = len(fwhms) -1
    sindx = np.random.randint(0,ls)
    findx = np.random.randint(0,lf)
    return slows[sindx], fwhms[findx]

def compute_el(curtime,star,apf_obs):
    apf_obs.date = curtime
    star.compute(apf_obs)
    actel = np.degrees(star.alt)
    actaz = np.degrees(star.az)
    return actel,actaz


def checkdate(datestr):
    match = re.match("(\d{4})(\/|\-)(\d{1,2})(\/|\-)(\d{1,2})",datestr)

    if not match:
        return False
    if int(match.group(3)) < 1 or int(match.group(3)) > 12:
        return False
    if int(match.group(5)) < 1 or int(match.group(5)) > 31:
        return False

    return True


def compute_simulation(result,curtime,star,apf_obs,slowdowns,fwhms,owner):
    actel,actaz = compute_el(curtime,star,apf_obs)
    actslow, actfwhm = rand_obs_sample(slowdowns,fwhms)
    actfwhm = gen_seeing_el(actfwhm,actel)
    lastfwhm = actfwhm
    lastslow = actslow
    meterrate = ExposureCalculations.getEXPMeter_Rate(result['VMAG'],result['BV'],actel,actfwhm,result['DECKER'])
    meterrate *= 1 + 0.11*np.random.randn(1)
    meterrate /= actslow
    specrate = ExposureCalculations.getSpec_Rate(result['VMAG'],result['BV'],actel,actfwhm,result['DECKER'])
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

    precision, true_error = Generate_Errors.compute_real_uncertainty(totcounts,result['BV'])
    if actaz < 180:
        actel *= -1.
    outstr = "%s %s %.5f %.1f %.1f %.2f %.2f %.2f %.2f %s" %(result['NAME'] , ephem.Date(curtime), ephem.julian_date(ephem.Date(barycentertime)), fexptime, totcounts,  actel,actaz, actfwhm, actslow, owner)
    print (outstr)

    return curtime, lastfwhm, lastslow, outstr



def init_sim_vals():
    vals = dict()

    sim_names = ('name','date','time','jd','etime','I2','El','Az','fwhm','slow','owner')
    for n in sim_names:
        vals[n] = []
    return vals, sim_names


def read_sim_lines(lines,sim_names,sim_vals):
    for l in lines:
        if re.search("\A\#",l):
            continue
        d = l.split()
        if len(d) == 11:
            for i in range(0,len(sim_names)):
                sim_vals[sim_names[i]].append(d[i])

    return

def sum_owner_times(vals):

    owner_tots = dict()
    owner_els = dict()
    owner_nexps = dict()
    for o in vals['owner']:
        if o not in owner_tots.keys():
            owner_tots[o] = 0.
            owner_els[o] = 0.
            owner_nexps[o] = 0
    owner_tots['total'] = 0.
    owner_els['total'] = 0.
    owner_nexps['total'] = 0

    for i in range(0,len(vals['owner'])):
        m = re.search("\ARECUR_A100",vals['owner'][i])
        owner_tots[vals['owner'][i]] += float(vals['etime'][i])
        owner_els[vals['owner'][i]] += float(vals['El'][i])
        owner_nexps[vals['owner'][i]] += 1
        owner_tots['total'] += float(vals['etime'][i])
        owner_els['total'] += float(vals['El'][i])
        owner_nexps['total'] += 1

    for o in owner_els.keys():
        if owner_nexps[o] > 0:
            owner_els[o] /= owner_nexps[o]

    return owner_tots,owner_els,owner_nexps
