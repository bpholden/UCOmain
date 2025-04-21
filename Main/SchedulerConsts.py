# Some variables that will soon be moved to a separate file
#TARGET_ELEVATION_MIN = 20 # this elevation is the physical minimum, below this the ADC does not work
TARGET_ELEVATION_MIN = 30 # temporary limit, currently the telescope occasionally fails at lower elevations
TARGET_ELEVATION_HIGH_MIN = 45 # this elevation is the preferred one for stars that will be high in the sky
TARGET_ELEVATION_MAX = 84
# a note about the values above
# the minimum elevation means that if the target is below that value during the observation
# the target will not selected
# the maximum means that if the target is above that value during the observation
# the target will not be selected
# the preferred minimum means that if a target will go above that
# value during the night, it will not be selected until it is above that value

TARGET_EXPOSURE_TIME_MIN = 1200
TARGET_MOON_DIST_MIN = 15
TARGET_MOON_DIST_MAX = 25

SUNEL_ENDLIM = -11.0
SUNEL_STARTLIM = -9.0
SUNEL_HOR = -3.2


# A few constants to make accessing the star table more readable
SEEING_TEMP = 15
SLOWDOWN_MIN = 0.6
SLOWDOWN_TEMP = 2.5
SLOWDOWN_THRESH = 5.0
SLOWDOWN_MAX = 10.0
SLOWDOWN_VMAG_LIM = 11.0
SEEING_THRESH = 30.0

PRI_DELTA = 20

EXP_LIM = 3e9
MAX_PRI = 3
#MAX_PRI = 2

DEFAULT_CERT = 'ucoscheduler-d79e797f0ade.json'
