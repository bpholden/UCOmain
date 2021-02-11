#!/usr/bin/env  /opt/kroot/bin/kpython
from __future__ import print_function
import os
import sys
sys.path.append("../master")
import apflog

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("needs a path")
        sys.exit()
    apflog.logpush(os.path.join(sys.argv[1],"observed_targets"))    
    apflog.logpush(os.path.join(sys.argv[1],"robot.log"))    
