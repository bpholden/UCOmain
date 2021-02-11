from __future__ import print_function
import sys
from optparse import OptionParser
from datetime import datetime
import time
import os
import pickle
import re

import numpy as np
import ephem

sys.path.append("../master")
import UCSCScheduler_V2 as ds
from  fake_apflog import *
import Coords
import ParseGoogledex
import SchedulerConsts as sc


if __name__ == "__main__":

    dt = datetime.utcfromtimestamp(int(time.time()))

    parser = OptionParser()
    (options, args) = parser.parse_args()    
    allnames, star_table, flags, stars = ParseGoogledex.parseGoogledex()
    for i in range(len(stars)):
        if star_table[i, sc.DS_APFPRI] > 5 and flags['template'][i] == 'N':
            row = star_table[i,:]
            print ( allnames[i], flags['owner'][i], star_table[i, sc.DS_APFPRI])
