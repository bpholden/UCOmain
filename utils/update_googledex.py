#!/usr/bin/env  /opt/kroot/bin/kpython
from __future__ import print_function

import sys
sys.path.append("../Main")
import ObservedLog
import ParseUCOSched

if __name__ == "__main__":

    helptext = """
    This script updates the Google Sheets with the latest observed targets and their corresponding data.
    It requires two input files:
    1. The list of observed targets (observed_targets)
    2. A local copy of the star table (googledex.dat)

    The newest saved versions will have a .1 suffix added to the filename.
    """

    if len(sys.argv) <= 2:
        print(helptext)
        sys.exit()
    if '--help' in sys.argv or '-h' in sys.argv:
        print(helptext)
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
