import math

def make_strs(deg,mn,sec,neg=False):
    """
    make_strs(deg,mn,sec,neg=False)

    Given degrees, minutes, seconds, and a boolean indicating if the value is negative, return
    a tuple of strings representing the values.  The strings are formatted as follows:
    deg: "%d"
    mn: "%d"
    sec: "%.4f"
    """

    sdeg = "%d" % (deg)
    smn = "%d" % (abs(mn))
    ssec = "%.4f" % (abs(sec))
    if neg:
        sdeg = "-" + sdeg
    return sdeg, smn, ssec

def get_RA_rad(hr, mn, sec):
    '''
    Docstring for get_RA_rad
    
    :param hr: hour value
    :param mn: minute value
    :param sec: second value
    :return: tuple (ra in radians, hr string, mn string, sec string)
    '''
    rv = None, "-1", "0", "0"
    try:
        hr = float(hr)
        mn = float(mn)
        sec = float(sec)
        if hr < 0 or hr > 23:
            return rv
        if mn < 0 or mn > 59:
            return rv
        if sec < 0 or sec >= 60:
            return rv
        ra_hours = hr + mn/60. + sec/3600.
        ra_hours *= 15 * math.pi/180.0

        shr, smn, ssec = make_strs(hr,mn,sec)

        return ra_hours, shr, smn, ssec
    except:
        return rv

def get_dec_rad(deg, mn, sec, neg=False):
    '''
    Docstring for get_dec_rad
    
    :param deg: Degree value
    :param mn: Minute value
    :param sec: Second value
    :param neg: Boolean indicating if the value is negative
    :return: tuple (dec in radians, deg string, mn string, sec string)
    '''
    rv = (None, "-90", "0", "0")
    try:
        deg = float(deg)
        mn = float(mn)
        sec = float(sec)
        if deg < -60 or deg > 90:
            return rv
        if mn < 0 or mn > 59:
            return rv
        if sec < 0 or sec >= 60:
            return rv
    except:
        return rv
    if deg < 0:
        neg = True

    if  mn < 0:
        neg = True

    if sec < 0:
        neg = True

    dec = abs(deg) + abs(mn)/60. + abs(sec)/3600.
    dec = dec * math.pi/180.
    if neg:
        dec *= -1

    sdeg, smn, ssec = make_strs(abs(deg),abs(mn),abs(sec),neg=neg)

    return dec, sdeg, smn, ssec
