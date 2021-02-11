# Generate the exposure times for all APF_pri > 5 targets
# Start with some assumed values, these may become command line options later

SLOWDOWN = 1.0
AVG_FWHM = 15.0
STAR_EL  = 70.0

import sys
import time
import os

os.chdir("/u/rjhanson/master/")


print os.getcwd()
from UCSCScheduler import *

def exposure_time(star):
    """
     Take in a star (getNext() target form ) and return a total exposure time
    """
    return 200

if __name__ == "__main__":
    print "APF exposure time test"
    with open("/u/rjhanson/master/apf_google_login",'r') as f:
        usr, pwd = pickle.load(f)
    sn, star_table, do_flag, stars = parseGoogledex(usr, pwd)

    times = []
    names = []
    pri = []
    err = []

    for i in range(len(stars)):
        if star_table[i, DS_APFPRI] < 5:
            continue
        precision = star_table[i, DS_ERR]
        if precision < 1.5 and star_table[i, DS_APFPRI] < 10:
            precision = 1.5
        names.append(sn[i])
        err.append(precision)
        pri.append(star_table[i, DS_APFPRI])
        i2counts = getI2(precision)
        exp_counts = getEXPMeter(i2counts, star_table[i, DS_BV])
        exp_time = getEXPTime(i2counts, star_table[i, DS_VMAG], star_table[i, DS_BV], STAR_EL, AVG_FWHM)
        times.append(exp_time)


    import matplotlib.pyplot as plt

    plt.scatter(times, err, c=pri, cmap=plt.get_cmap("jet"), edgecolor='none', size=50)
    plt.xlabel("Exposure Time")
    plt.ylabel("Desired Precision")
    plt.show()
