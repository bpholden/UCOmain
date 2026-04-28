from __future__ import print_function
import os
import argparse
import sys

import NightSim

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile",dest="infile",default="all_dates.simout")
    parser.add_argument("-g","--gfile",dest="googlefile",default="googledex.dat")
    parser.add_argument("-o","--outfile",dest="outfile",default="sim_summary")
    args = parser.parse_args()

    if not os.path.exists(args.infile):
        print ("cannot open %s" % (args.infile))
        sys.exit()
    if not os.path.exists(args.googlefile):
        print ("cannot open %s" % (args.googlefile))
        sys.exit()

    return args

###
def main():
    options = parse_args()
    try:
        fp = open(options.infile)
    except Exception as e:
        print ("Cannot open %s: %s" % (options.infile,e))
        sys.exit()

    simin = fp.read()
    fp.close()
    lines = simin.split("\n")

    sim_vals, sim_names = NightSim.init_sim_vals()

    NightSim.read_sim_lines(lines, sim_names, sim_vals)

    owner_times, owner_els, owner_nexps = NightSim.sum_owner_times(sim_vals)

    for o, t in owner_times.items():
        if owner_nexps[o] > 0:
            print ("%s %.2f %.2f %.2f %d" %\
                 (o, t/3600., t/owner_nexps[o],\
                  owner_els[o], owner_nexps[o]))

if __name__ == "__main__":
    main()