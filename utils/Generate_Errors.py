import numpy as np

# needs to take output from sim_night
# generate the skeleton of a vels file for each target
# a file will be MJD error I2 Counts exptime

def compute_precision(i2counts,cols):

    blue = (cols < 1.2)
    red = (cols >= 1.2)
    gerr = .436/np.log(10) # fractional err
    merr = .363/np.log(10)

    A_gk = 4.47
    B_gk = -1.58
    A_m = 4.14
    B_m = -1.73

    if blue:
        mean = (np.log10(i2counts) - A_gk)/B_gk
    else:
        mean = (np.log10(i2counts) - A_m)/B_m                        
#    means = np.zeros_like(i2counts)
#    means[blue] += (np.log10(i2counts[blue]) - A_gk)/B_gk
#    means[red] += (np.log10(i2counts[red]) - A_m)/B_m

    devs = np.random.normal(size=len(i2counts))
    if blue:
        devs *= gerr
    else:
        devs *= merr
    devs += mean
    return 10**devs[0]

def jitter(cols):
    # just use median
    
    # jitter = np.zeros_like(cols)
    # jitter[((cols > 0.4) & (cols < 0.7))] += 2.3 + 17.4*0.02
    # jitter[((cols > 0.7) & (cols < 1.0))] += 2.1 + 4.7 *0.02
    # jitter[((cols > 1.0) & (cols < 1.3))] += 1.6 - 0.003 * 0.3 
    # jitter[((cols > 1.3) & (cols < 1.6))] += 2.1 + 2.7 * 0.5

    if cols < 0.7:
        jitter =  2.3 + 17.4*0.02
    elif cols > 0.7 and cols < 1.0:
        jitter = 2.1 + 4.7 *0.02
    elif cols > 1.0 and cols < 1.3:
        jitter =  1.6 - 0.003 * 0.3 
    else:
        jitter =  2.1 + 2.7 * 0.5
        
    return jitter

def compute_real_uncertainty(i2counts,cols):

    precision = compute_precision(i2counts,cols)
    #    true= np.zeros_like(precision)
    true = 0.0
    true += precision
    floor = 1. # floor is about 1 m/s
    true += floor
    jit = jitter(cols)
    true **= 2
    jit **= 2
    true += jit
    true = np.sqrt(true)
    
    return precision, true
