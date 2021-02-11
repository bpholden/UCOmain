# Some variables that will soon be moved to a separate file
TARGET_ELEVATION_MIN = 20 # this elevation is the physical minimum, below this the ADC does not work
TARGET_ELEVATION_HIGH_MIN = 45 # this elevation is the preferred one for stars that will be high in the sky
TARGET_ELEVATION_MAX = 85
TARGET_EXPOSURE_TIME_MAX = 3 * 60 * 60 # 3 hour
TARGET_EXPOSURE_TIME_MIN = 1240
TARGET_MOON_DIST_MIN = 15
TARGET_MOON_DIST_MAX = 25

SUNEL_ENDLIM = -10.0
SUNEL_STARTLIM = -9.0
SUNEL_HOR = -3.2


# A few constants to make accessing the star table more readable

SLOWDOWN_MIN = 0.6
SLOWDOWN_THRESH = 2.5
SLOWDOWN_MAX = 5.0
SLOWDOWN_VMAG_LIM = 9.0
SEEING_THRESH = 30.0

PRI_DELTA = 20

EXP_LIM = 3e9
MAX_PRI = 10
