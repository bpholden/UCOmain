#!/usr/bin/env  /opt/kroot/bin/kpython
from __future__ import print_function

import sys
sys.path.append("../Main")
import ObservedLog
import ParseUCOSched

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print("needs a observed log and a local copy of the star table")
        sys.exit()
    fn = sys.argv[1]
    outfn = sys.argv[2]
    obslog = ObservedLog.ObservedLog(fn)

    if len(obslog.names) > 0:
        if obslog.sheetns[0] is None:
            sheetns = set(obslog.owners)
        else:
            sheetns = set(obslog.sheetns)

        ParseUCOSched.update_online_sheets(fn,outfn=outfn)
