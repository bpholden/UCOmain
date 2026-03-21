import numpy as np

# needs to take output from sim_night
# generate the skeleton of a vels file for each target
# a file will be MJD error I2 Counts exptime

def compute_precision(i2counts, cols):
    """
    Compute the precision of a measurement based on the I2 counts and the color of the star.
    The precision is computed using the formula:
    precision = 10^(mean + dev)
    where mean is computed from the I2 counts and the color, and dev is a random
    drawn from a normal distribution with a mean and standard deviation from empirical data. The precision is then converted from log space to linear space.
    Parameters:
    i2counts (float): The I2 counts for the observation.
    cols (float): The color of the star (B-V).
    Returns:
    float: The computed precision in m/s.""" 
    blue = cols < 1.2
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
    """
    Compute the jitter of a measurement based on the color of the star.
    The jitter is computed using empirical data that shows a relationship between 
    the color of the star and the expected jitter. 
    The jitter is computed as a piecewise function based on the color of the star.
    Parameters:
    cols (float): The color of the star (B-V).
    Returns:
    float: The computed jitter in m/s.
    """
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
    """
    Compute the real uncertainty of a measurement based on the I2 counts 
    and the color of the star.
    The real uncertainty is computed by combining the precision 
    and the jitter in quadrature.
    Parameters:
    i2counts (float): The I2 counts for the observation.
    cols (float): The color of the star (B-V).
    Returns:
    tuple: A tuple containing the precision and the real uncertainty in m/s.
    """
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
