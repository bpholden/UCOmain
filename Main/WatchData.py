from __future__ import print_function
from datetime import datetime, timedelta
import glob
import os
import os.path
import sys
import threading
import time

import astropy.io.fits
import astropy.time

import apflog
import ktl

class WatchData(threading.Thread):

    def __init__(self):
        threading.Thread.__init__(self)
        self.setDaemon(True)

        self.apfucam = ktl.Service('apfucam')
        self.outdir = self.apfucam['OUTDIR'].read()
        self.checked_list = {}
        self.signal = True

    def badFile(self,filename):
        hdr = astropy.io.fits.getheader(filename)

        db    = astropy.time.Time(hdr['DATE-BEG'],format='isot')
        de    = astropy.time.Time(hdr['DATE-END'],format='isot')
        midpt = astropy.time.Time(hdr['THEMIDPT'],format='isot')

        if midpt > de:
            rv = True
        elif midpt < db:
            rv = True
        else:
            rv = False

        return rv

    def checkFiles(self):

        file_list = glob.glob(os.path.join(self.outdir,"*.fits"))

        bad = []
        for fn in file_list:
            if fn not in self.checked_list.keys():
                rv = self.badFile(fn)
                self.checked_list[fn] = 'done'
                if rv:
                    bad.append(fn)

        if len(bad) == 0:
            return None

        return bad

    def stop(self):
        self.signal = False
        threading.Thread._Thread__stop(self)
    
    def run(self):
        
        while self.signal:
            bad = self.checkFiles()

            if bad is not None:
#                apflog("%s is bad" % (bad),echo=True,level='alert')
                print("%s is bad" % (bad))
            else:
                print("All good")
#                apflog("All good", echo=True)

            time.sleep(60)

if __name__ == "__main__":

    wd = WatchData()
    wd.start()
    while wd.signal:
        try:
            time.sleep(10)
        except KeyboardInterrupt:
            print("%s has been killed by user." % (wd.name))
            wd.stop()
            sys.exit()
        except:
            print("%s killed by unknown." % (wd.name))
            wd.stop()
            sys.exit()

