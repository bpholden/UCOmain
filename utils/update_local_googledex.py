#!/usr/bin/env  /opt/kroot/bin/kpython
from __future__ import print_function

import sys
import time
import datetime

sys.path.append("../Main")
import ParseUCOSched

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print("needs a log filename and googledex filename")
        sys.exit()
    obfn = sys.argv[1]
    gdfn = sys.argv[2]

    dt = datetime.datetime.now()
    ParseUCOSched.update_local_starlist(dt,outfn=gdfn,observed_file=obfn)
