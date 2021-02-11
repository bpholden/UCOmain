import UCSCScheduler as ds
import numpy as np
import pickle
from datetime import datetime
import time

def compute_jd(dt):

#    zero : 0h May 24, 1968
    zero = datetime(1968,5,24,0,0,0,0)
    delta = dt - zero
    mjd = delta.days + delta.seconds / 86400
    jd = mjd + 2440000.5
    return jd

        
if __name__ == "__main__":

    allnames, star_table, do_flag, stars  = ds.parseGoogledex("","")

    for i in range(len(stars)):
        if star_table[i, ds.DS_APFPRI] < 5:
            continue

        now = datetime.utcfromtimestamp(int(time.time()))
        jd = compute_jd(now)

        delta = jd - star_table[i, ds.DS_LAST]
        if delta > star_table[i, ds.DS_CAD]:
            print "%15s %4.1f %7.1f %d %7.1f " % (allnames[i],star_table[i, ds.DS_APFPRI],star_table[i, ds.DS_LAST],star_table[i, ds.DS_CAD],delta)


