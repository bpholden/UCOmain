import time
from skyfield.api import load, wgs84, EarthSatellite

ts = load.timescale()
line1 = '1 51850U 22021A   22185.91365633  .00000107  00000-0  00000-0 0  9998'
line2 = '2 51850   0.0574 112.0706 0002585 105.6714 257.2584  1.00272770  1312'
line1 = '1 29643U 06054A   22187.58585272 -.00000036  00000-0  00000+0 0  9990'
line2 = '2 29643   0.0190  46.1208 0001677  55.6121 282.5292  1.00272324 57090'

satellite = EarthSatellite(line1, line2, 'GOES18', ts)
lat = 37 + 20 / 60 + 33.1 / 3600
long = -121 - 38/60 - 17.7/3600
apf = wgs84.latlon(lat,long, elevation_m=1274)

t = ts.now()
geocentric = satellite.at(t)
print(geocentric.position.km)
lat, lon = wgs84.latlon_of(geocentric)
#print('Latitude:', lat)
#print('Longitude:', lon)

difference = satellite - apf

while True:

    t = ts.now()
    topocentric = difference.at(t)
    alt, az, distance = topocentric.altaz()

    if alt.degrees > 0:
        print(t.tt, alt.degrees, az.degrees)
#    print('Distance: {:.1f} km'.format(distance.km))
    time.sleep(3)
