from __future__ import print_function

import datetime
import time

import ephem
import numpy as np

import SchedulerConsts

def calc_preferred_angle(shiftwest, sun_el, delta_t):
    """ 
    calc_preferred_angle: Calculate the preferred elevation angle for the telescope
    Args:
        shiftwest: Boolean, True if the telescope is shifted west
        sun_el: The current elevation of the sun ( degrees )
        sun_az: The current azimuth of the sun ( degrees )
        delta_t: The time from the sunset ( seconds )
    """
    bottom_angle = SchedulerConsts.SUNEL_STARTLIM-15 # typically -24 degrees
    offset = 0.0
    preferred_angle = 90

    if shiftwest:
        if sun_el > (bottom_angle):
            offset = 3*(sun_el - bottom_angle) # note, this is positive
            preferred_angle = 90 - offset
        elif delta_t < 3600. and delta_t > 0:
            offset = 45*(1. - delta_t/3600)
            preferred_angle = 90 - offset

    return preferred_angle, offset

def visible(observer, stars, obs_len, pref_min_el=SchedulerConsts.TARGET_ELEVATION_HIGH_MIN, \
                min_el=SchedulerConsts.TARGET_ELEVATION_MIN, \
                   max_el=SchedulerConsts.TARGET_ELEVATION_MAX, shiftwest=False, delta_t=0):
    """ Args:
            stars: A list of pyephem bodies to evaluate visibility of
            observer: A pyephem observer to use a the visibility reference
            obs_len: A list of observation lengths ( Seconds ). 
              This is the time frame for which visibility is checked
            pref_min_el: Preferred minimum body elevation to be visible ( degrees )
            min_el: The minimum body elevation to be visible ( degrees ) 
              - only use this if star never goes above preferred limit
            max_el: The maximum body elevation to be visible ( degrees )
        Returns:
            Boolean list representing if body[i] is visible

        Notes: Uses the observer's current date and location
    """
    # Store the previous observer horizon and date since we change these
    prev_horizon = observer.horizon
    cdate = observer.date
    ret = []
    start_elevations = []
    scaled_elevations = []
    observer.horizon = str(min_el)

    sun = ephem.Sun(observer)
    sun_el = np.degrees(sun.alt)

    preferred_angle, offset = calc_preferred_angle(shiftwest, sun_el, delta_t)

    # Now loop over each body to check visibility
    for star, obs_time in zip(stars, obs_len):

        # Is the target visible now?

        observer.date = ephem.Date(cdate)
        star.compute(observer)
        cur_el = np.degrees(star.alt)
        cur_az = np.degrees(star.az)
        start_elevations.append(cur_el)

        if cur_el > max_el:
            scaled_elevations.append(cur_el)
            ret.append(False)
            continue
        if cur_el < pref_min_el and in_east:
            scaled_elevations.append(cur_el)
            ret.append(False)
            continue
        if cur_el < min_el:
            scaled_elevations.append(cur_el)
            ret.append(False)
            continue

        obs_time_days = obs_time / 86400
        if obs_time > 0:
            # Is the target visible at the end of the observations?
            observer.date = ephem.Date(cdate + obs_time_days)
            star.compute(observer)
            fin_el = np.degrees(star.alt)

            observer.date = ephem.Date(cdate + obs_time_days/2)
            star.compute(observer)
            mid_el = np.degrees(star.alt)

        else:
            fin_el = cur_el
            mid_el = cur_el

        diff = np.abs(star.a_dec - observer.lat)
        transit_alt = 90.0 - np.degrees(diff)
        se = 90.0 - (transit_alt - mid_el)

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
            next_set = observer.next_setting(star)
        except:
            # If it never sets, no need to worry.
            pass
        else:
            # Making the assumption that next_set is a datetime object. Might not be the case
            if next_set < obs_time:
                # The object will set before the observation finishes
                ret.append(False)
                continue

        observer.horizon = max_el
        star.compute(observer)

        try:
            next_rise = observer.next_rising(star)
        except:
            # If the body never rises above the max limit no problem
            pass
        else:
            if next_rise < obs_time:
                # The object rises above the max el before the observation finishes
                ret.append(False)
                continue

        observer.horizon = str(pref_min_el)
        star.compute(observer)
        if not star.neverup:
            # will transit above preferred elevation and still rising
            try:
                if ((star.set_time-star.rise_time) > obs_time_days) \
                    and cur_el < pref_min_el \
                        and np.degrees(star.az) < 180:
                    # this star is currently low on the horizon
                    # but will not be above the preferred elevation
                    # for the requested exposure time
                    ret.append(False)
                    continue
            except:
                ret.append(False)
                continue
        # Everything seems to be fine, so the target is visible!
        ret.append(True)
#	apflog( "is_visible(): done searching targets", echo=True)
    observer.horizon = prev_horizon
    return ret, np.array(start_elevations), np.array(scaled_elevations)

def test_main():
    # This is a test function to check the visibility of a star
    # It will be run when this file is executed
    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = str(SchedulerConsts.TARGET_ELEVATION_MIN)
    apf_obs.date = datetime.datetime.utcfromtimestamp(int(time.time()))

    test_star = ephem.FixedBody()
    test_star._ra = ephem.hours(":".join(["1", "44", "4.083"]))
    test_star._dec = ephem.degrees(":".join(["-15", "56", "14.93"]))
    tret, tse, tsce = visible(apf_obs, [test_star], [0.])
    print(tret, tse, tsce)
    tret, tse, tsce = visible(apf_obs, [test_star], [400.])
    print(tret, tse, tsce)

if __name__ == '__main__':
    test_main()
