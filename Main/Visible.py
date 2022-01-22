from __future__ import print_function

from datetime import datetime, timedelta
import ephem
import numpy as np
import os
import sys
import time

import SchedulerConsts



def visible(observer, stars, obs_len, pref_min_el=SchedulerConsts.TARGET_ELEVATION_HIGH_MIN, min_el=SchedulerConsts.TARGET_ELEVATION_MIN,
                   max_el=SchedulerConsts.TARGET_ELEVATION_MAX,shiftwest=False):
    """ Args:
            stars: A list of pyephem bodies to evaluate visibility of
            observer: A pyephem observer to use a the visibility reference
            obs_len: A list of observation lengths ( Seconds ). This is the time frame for which visibility is checked
            pref_min_el: Preferred minimum body elevation to be visible ( degrees )            
            min_el: The minimum body elevation to be visible ( degrees ) - only use this if star never goes above preferred limit
            max_el: The maximum body elevation to be visible ( degrees )
        Returns:
            Boolean list representing if body[i] is visible

        Notes: Uses the observer's current date and location
    """
    # Store the previous observer horizon and date since we change these
    prev_horizon = observer.horizon
    cdate = observer.date
    ret = []
    fin_elevations = []
    start_elevations = []
    scaled_elevations = []
    observer.horizon = str(min_el)

    sun = ephem.Sun(observer)
    sun_el = np.degrees(sun.alt)
    sun_az = np.degrees(sun.az)

    bottom_angle = SchedulerConsts.SUNEL_STARTLIM-15 # typically -24 degrees
    
    if sun_el > (bottom_angle) and sun_az > 180 and shiftwest:
        offset = 3*(sun_el - bottom_angle) # note, this is positive
        preferred_angle = (90 - offset)
    else:
        offset = 0.0
        preferred_angle = 90 

    
    # Now loop over each body to check visibility
    for s, dt in zip(stars, obs_len):

        # Is the target visible now?

        observer.date = ephem.Date(cdate)
        s.compute(observer)
        cur_el = np.degrees(s.alt)
        cur_az = np.degrees(s.az)
        start_elevations.append(cur_el)
        
        if cur_el < min_el or cur_el > max_el:
            fin_elevations.append(cur_el)
            scaled_elevations.append(cur_el)
            ret.append(False)
            continue

        if dt > 0:
            # Is the target visible at the end of the observations?
            observer.date = ephem.Date(cdate + dt/86400.)
            s.compute(observer)
            fin_el = np.degrees(s.alt)
            fin_elevations.append(fin_el)
        else:
            fin_elevations.append(cur_el)
            fin_el = cur_el
        
        diff = np.abs(s.a_dec - observer.lat)
        transit_alt = 90.0 - np.degrees(diff)
        se = 90.0 - (transit_alt - cur_el) 

        if offset > 0:
            if cur_az < 180:
                se -= offset
            else:
                se = 90 - np.abs(preferred_angle - se)
        
        scaled_elevations.append(se)
        
        if fin_el < min_el or fin_el > max_el:
            ret.append(False)
            continue
            

        # Does the target remain visible through the observation?
        # The next setting/rising functions throw an exception if the body never sets or rises
        # ex. circumpolar 
        try:
            next_set = observer.next_setting(s)
        except:
            # If it never sets, no need to worry.
            pass
        else:
            # Making the assumption that next_set is a datetime object. Might not be the case
            if next_set < dt:
                # The object will set before the observation finishes
                ret.append(False)
                continue
        #  apflog( "is_visible(): Does the target remain visible through the observation?", echo=True)
        observer.horizon = max_el
        s.compute(observer)

        try:
            next_rise = observer.next_rising(s)
        except:
            # If the body never rises above the max limit no problem
            pass
        else:
            if next_rise < dt:
                # The object rises above the max el before the observation finishes
                ret.append(False)
                continue
        #   apflog( "is_visible(): If the body never rises above the max limit no problem", echo=True)
        observer.horizon = str(pref_min_el)
        s.compute(observer)
        if not s.neverup:
            # will transit above preferred elevation and still rising
            try:
                if ((s.set_time-s.rise_time) > dt/86400.) and cur_el < pref_min_el and np.degrees(s.az) < 180:
                    # this star is currently low on the horizon but will not be above the preferred elevation for the requested exposure time
                    ret.append(False)
                    continue
            except:
                pass
        # Everything seems to be fine, so the target is visible!
        ret.append(True)
#	apflog( "is_visible(): done searching targets", echo=True)
    observer.horizon = prev_horizon
    return ret, np.array(start_elevations), np.array(fin_elevations), np.array(scaled_elevations)



if __name__ == '__main__':
    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = str(SchedulerConsts.TARGET_ELEVATION_MIN)
    apf_obs.date = datetime.utcfromtimestamp(int(time.time()))

    star = ephem.FixedBody()
    star._ra = ephem.hours(":".join(["1", "44", "4.083" ]))
    star._dec = ephem.degrees(":".join(["-15", "56", "14.93"]))
    ret, se, fe, sce = visible(apf_obs,[star],[0.])
    print (ret, se, fe)
    ret, se, fe, sce = visible(apf_obs,[star],[400.])
    print (ret, se, fe, sce)

    
    star = ephem.FixedBody()
    star._ra = ephem.hours(":".join(["1", "44", "4.083" ]))
    star._dec = ephem.degrees(":".join(["-15", "56", "14.93"]))
    ret, se, fe, sce = visible(apf_obs,[star],[400.],)
    print (ret, se, fe, sce)
