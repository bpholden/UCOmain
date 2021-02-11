import sys
sys.path.append("../master")
from ExposureCalculations import *
import UCSCScheduler as ds
import numpy as np
import matplotlib.pyplot as plt

THRESHOLD = 3

STAR_EL = 70
AVG_FWHM = 11

def calc_fin_pre(i2counts,exp_time,bmv):
    ratio = exp_time / THRESHOLD
    ni2counts = i2counts / ratio
    if bmv > 1.2:
        A = 4.14
        B = -1.73
    else:
        A = 4.43
        B = -2.26

    log_npre = (np.log10(ni2counts) - A )/ B
    return 10**log_npre,ni2counts

        
if __name__ == "__main__":

    #    for mag in [7,8,9,10]:
    #        for bmv in [1.0, 1.4]:

    cols = ['violet','b','r']
    mags = [4.7,6.53,10.65]
    bmvs = [0.77, 0.68, 0.996]
    for i in range(0,len(mags)):
        uncs = np.linspace(1,10)
        mag = mags[i]
        bmv = bmvs[i]
        if bmv > 1.2:
            stype = "M"
            i2counts = ds.getI2_M(uncs)
        else:
            stype = "K"                
            i2counts = ds.getI2_K(uncs)

        exp_counts = ds.getEXPMeter(i2counts, bmv)
        exp_time = ds.getEXPTime(i2counts, mag, bmv, STAR_EL, AVG_FWHM)

            
#        print "%s %4.1f %3.1f %7.0f %.1g" % (stype,mag,precision,exp_time,exp_counts)

        plt.semilogy(uncs,exp_time,c=cols[i])

        plt.ylabel("Exposure Time (sec)")
        plt.xlabel("Desired Precision m/s")

    plt.ylim(1,2500)
    plt.show()
