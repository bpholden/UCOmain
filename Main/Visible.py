from __future__ import print_function

import datetime

import astroplan
#from astroplan import (AltitudeConstraint, AirmassConstraint,
#                       AtNightConstraint)

import astropy.units
import astropy.time
import astropy.coordinates
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

def visible(observer, stars, obs_len, ptime,
            pref_min_el=SchedulerConsts.TARGET_ELEVATION_HIGH_MIN, \
                min_el=SchedulerConsts.TARGET_ELEVATION_MIN, \
                   max_el=SchedulerConsts.TARGET_ELEVATION_MAX, 
                   shiftwest=False, delta_t=0):
    """ Args:
            stars: A list of astroplan FixedTarget objects to evaluate visibility of
            observer: An astroplan Observer to use as the visibility reference
            obs_len: A list of observation lengths ( astropy.time.TimeDelta ). 
              This is the time frame for which visibility is checked
            ptime: The current observation time ( astropy.time.Time ). 
            pref_min_el: Preferred minimum body elevation to be visible ( degrees )
            min_el: The minimum body elevation to be visible ( degrees ) 
              - only use this if star never goes above preferred limit
            max_el: The maximum body elevation to be visible ( degrees )
        Returns:
            Boolean list representing if body[i] is visible

        Notes: Uses the observer's current date and location
    """
    # Store the previous observer horizon and date since we change these
    ret = []
    start_elevations = []
    scaled_elevations = []

    sun = observer.sun_altaz(ptime)
    sun_el = float(sun.alt.value)

    preferred_angle, offset = calc_preferred_angle(shiftwest, sun_el, delta_t)

    constraints = [astroplan.AltitudeConstraint(min_el*astropy.units.deg, max_el*astropy.units.deg)]
               #AirmassConstraint(5), AtNightConstraint.twilight_civil()]
    # Are targets *always* observable in the time range?


    # Now loop over each body to check visibility
    for star, obs_time in zip(stars, obs_len):

        # Is the target visible now?
        end_time = ptime + obs_time
        always_observable = astroplan.is_always_observable(constraints, observer, star, \
                                                           time_range=(ptime, end_time))
        always_observable = bool(always_observable[0])
        ret.append(always_observable)
        star_pos = observer.altaz(ptime, star)
        cur_el = float(star_pos.alt.value)
        cur_az = float(star_pos.az.value)
        start_elevations.append(cur_el)


        if always_observable is False:
            scaled_elevations.append(cur_el)
            continue

        new_time = ptime + astropy.time.TimeDelta(obs_time.to_value('sec')/2, format='sec')
        if obs_time.to_value('sec') > 0:  # Ensure obs_time is in seconds for comparison
            # mid point elevation
            mid_pos = observer.altaz(new_time, star)
            mid_el = mid_pos.alt.value
        else:
            mid_el = cur_el

        diff = np.abs(star.dec.value - observer.location.lat.value)
        transit_alt = 90.0 - diff
        se = observer.target_hour_angle(new_time, star).to(astropy.units.deg)
        if se.value  < 180:
            in_east = True
            se = 90 - se.value
        else:
            in_east = False
            se = se.value - 270

        if offset > 0:
            if cur_az < 180:
                se -= offset
            else:
                se = 90 - np.abs(preferred_angle - se)

        scaled_elevations.append(se)

        if transit_alt > pref_min_el and in_east and mid_el < pref_min_el:
            # will transit above preferred elevation and still rising
            ret.pop()
            ret.append(False)

    return ret, np.array(start_elevations), np.array(scaled_elevations)

def test_main():
    # This is a test function to check the visibility of a star
    # It will be run when this file is executed
    # Generate a astroplan observer for the APF

    import SunPos

    apf_obs_date = astropy.time.Time(datetime.datetime.now(tz=datetime.timezone.utc), scale='utc')
    apf_obs = SunPos.make_APF_obs()

    test_star_ra = ":".join(["1", "44", "4.083"])
    test_star_dec = ":".join(["-15", "56", "14.93"])
    test_coord = astropy.coordinates.SkyCoord(test_star_ra, test_star_dec,\
                                             unit=(astropy.units.hourangle, astropy.units.deg))
    test_star = astroplan.FixedTarget(coord=test_coord, name="Test Star")

    obstime = astropy.time.TimeDelta(400., format='sec')
    tret, tse, tsce = visible(apf_obs, [test_star], [obstime], apf_obs_date)
    print(tret, tse, tsce)

    test_star_ra = ":".join(["19", "32", "21.5902"])
    test_star_dec = ":".join(["69", "39", "40.2350"])
    test_coord = astropy.coordinates.SkyCoord(test_star_ra, test_star_dec,\
                                             unit=(astropy.units.hourangle, astropy.units.deg))
    test_star = astroplan.FixedTarget(coord=test_coord, name="sig Dra")
    tret, tse, tsce = visible(apf_obs, [test_star], [obstime], apf_obs_date)
    print(tret, tse, tsce)

if __name__ == '__main__':
    test_main()
