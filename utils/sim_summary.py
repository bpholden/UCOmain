from __future__ import print_function
import os
import re
import optparse
import sys

def init_sim_vals():
    vals = dict()

    sim_names = ('name','date','time','jd','etime','I2','El','Az','fwhm','slow','owner')
    for n in sim_names:
        vals[n] = []
    return vals, sim_names

def sum_owner_times(vals):

    owner_tots = dict()
    owner_els = dict()
    owner_nexps = dict()
    for o in vals['owner']:
        if o not in owner_tots.keys():
            owner_tots[o] = 0.
            owner_els[o] = 0.
            owner_nexps[o] = 0
    owner_tots['bstar'] = 0.
    owner_els['bstar'] = 0.
    owner_nexps['bstar'] = 0
            
    for i in range(0,len(vals['owner'])):
        m = re.search("\AHR",vals['name'][i])
        if not m:
            owner_tots[vals['owner'][i]] += float(vals['etime'][i])
            owner_els[vals['owner'][i]] += float(vals['El'][i])
            owner_nexps[vals['owner'][i]] += 1
        else:
            owner_tots['bstar'] = float(vals['etime'][i])
            owner_els['bstar'] += float(vals['El'][i])
            owner_nexps['bstar'] += 1

            
    for o in owner_els.keys():
        if owner_nexps[o] > 0:
            owner_els[o] /= owner_nexps[o]
            
    return owner_tots,owner_els,owner_nexps

def parse_args():
    parser = optparse.OptionParser()
    parser.add_option("-i","--infile",dest="infile",default="sim_master.simout")
    parser.add_option("-g","--gfile",dest="googlefile",default="googledex.dat")
    parser.add_option("-o","--outfile",dest="outfile",default="sim_summary")        
    (options, args) = parser.parse_args()    

    if not os.path.exists(options.infile):
        print ("cannot open %s" % (options.infile))
        sys.exit()
    if not os.path.exists(options.googlefile):
        print ("cannot open %s" % (options.googlefile))
        sys.exit()
            
    return options

###

options = parse_args()
try:
    fp = open(options.infile)
except Exception as e:
    print ("Cannot open %s: %s" % (options.infile,e))
    
simin = fp.read()
fp.close()
lines = simin.split("\n")

sim_vals,sim_names = init_sim_vals()

for l in lines:
    if re.search("\A\#",l):
        continue
    d = l.split()
    if len(d) == 11:
        for i in range(0,len(sim_names)):
            sim_vals[sim_names[i]].append(d[i])

owner_times, owner_els, owner_nexps = sum_owner_times(sim_vals)

for o in owner_times.keys():
    if owner_nexps[o] > 0:
        print ("%s %.2f %.2f %.2f %d" % (o, owner_times[o]/3600., owner_times[o]/owner_nexps[o], owner_els[o], owner_nexps[o]))
