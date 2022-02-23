from __future__ import print_function
import os
import re
import optparse
import sys

import NightSim

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

sim_vals,sim_names = NightSim.init_sim_vals()

NightSim.read_sim_lines(lines,sim_names,sim_vals)


owner_times, owner_els, owner_nexps = NightSim.sum_owner_times(sim_vals)

for o in owner_times.keys():
    if owner_nexps[o] > 0:
        print ("%s %.2f %.2f %.2f %d" % (o, owner_times[o]/3600., owner_times[o]/owner_nexps[o], owner_els[o], owner_nexps[o]))
