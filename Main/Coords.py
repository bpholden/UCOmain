import numpy as np

def makeStrs(deg,mn,sec,neg=False):

        
    sdeg = "%d" % (deg)
    smn = "%d" % (abs(mn))
    ssec = "%.4f" % (abs(sec))
    if neg:
        sdeg = "-" + sdeg
    return sdeg, smn, ssec

def getRARad(hr, mn, sec):
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
        ra_hours *= 15 * np.pi/180.0

        shr, smn, ssec = makeStrs(hr,mn,sec)
        
        return ra_hours, shr, smn, ssec
    except:
        return rv

def getDECRad(deg, mn, sec, neg=False):
    rv = (None, "-90", "0", "0")
    try:
        deg = float(deg)
        mn = float(mn)
        sec = float(sec)
        if deg < -60 or deg > 90:
            return rv
        if mn > 59:
            return rv
        if sec >= 60:
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
    dec = dec * np.pi/180.
    if neg:
        dec *= -1

    sdeg, smn, ssec = makeStrs(abs(deg),abs(mn),abs(sec),neg=neg)
    
    return dec, sdeg, smn, ssec

        
def getCoordStr(floatval,isRA=False):

    neg = False
    nround = 2
    if isRA:
        floatval /= 15.
        nround = 3
    if floatval < 0:
        neg = True
    floatval = abs(floatval)
    deghrval = int(floatval)
    minval = (floatval % 1) * 60.0 
    secval = round( (minval % 1) *60.0, nround)

    if neg and deghrval != 0:
        ret = "-" + str(deghrval) + ' '
    else:
        ret = str(deghrval) + ' '
    if neg and deghrval == 0 and minval != 0:
        ret += "-" + str(int(minval)) + ' '
    else:
        ret += str(int(minval)) + ' '
    if neg and deghrval == 0 and minval == 0:
        ret += "-" + str(secval)
    else:
        ret += str(secval)
    return ret


def getLST(date, longitude):
    """Take a datetime and longitude and calculate the Local Sidereal Time."""
    # Assumes date is a datetime object, and that the longitude is formatted as in PyEphem 

    ll = [float(v) for v in longitude.split(':')]
    if ll[0] > 0:
        sign = 1
    else:
        sign = -1
    ut = date.hour + date.minute/60. + date.second/3600.
    lng = ll[0] + sign*ll[1]/60. + sign*ll[2]/3600.
    d  = ephem.julian_date() - 2451545.0
    lst = 100.46 + 0.985647 * d + lng + 15*ut
    return lst % 360.


        
def getElAz(ra, dec, lat, lng, time):
    """Given RA, DEC, Latitude, and a time, returns the corresponding elevation and azimuth angles
       Works with single values, or numpy arrays
       """
    lst = getLST(time, lng)
    ha = ((lst- np.degrees(ra)) % 360.) * np.pi/180.
    el = np.arcsin(np.sin(dec) * np.sin(lat) + \
                   np.cos(dec) * np.cos(lat) * np.cos(ha))
    az = np.arccos( (np.sin(dec) - np.sin(el)*np.sin(lat)) / \
                         (np.cos(el) * np.cos(lat)))
    return (np.degrees(el), np.degrees(az))

