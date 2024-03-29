from __future__ import print_function

import re
import os
try:
    from apflog import *
except:
    from fake_apflog import *


class ObservedLog():
    def __init__(self,filename=None):
        self.names = []
        self.times = []
        self.temps = []
        self.owners = []
        self.sheetns = []
        self.filename = filename

        if filename is not None:
            self.read_observed_log()

    def __str__(self):
        return "< ObservedLog %s >" % self.filename

    def __repr__(self):
        return "< ObservedLog %s >" % self.filename

    def parse_key_vals(self,line):
        keyvals = dict()
        ovals = []
        line_list = line.split()
        for ls in line_list:
            split_vals = ls.split('=')
            if len(split_vals) == 2:
                keyvals[split_vals[0]] = split_vals[1]
            else:
                ovals.append(split_vals[0])

        return ovals, keyvals

    def read_observed_log(self):
        """ read_observed_log parses a file to find the object names and times
        ObservedLog.read_observed_log(filename)

        The following attributes in the object are filled:

        names - list of names, must be first column of file called filename
        times - times either as a timestamp in second column or 
            a (hour,minute) tuple from a scriptobs line
        temps - a list of template observations
        owners - a list of owners
        sheetns - a list of names for the sheets
        """

        try:
            f = open(self.filename, 'r')
        except IOError:
            apflog( "Couldn't open %s" % self.filename, level="warn", echo=True)
            return

        for line in f:
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                continue

            ovals, keyvals = self.parse_key_vals(line)
            self.names.append(ovals[0])
            if 'uth' in list(keyvals.keys()):
                self.times.append( ( int(keyvals['uth']), int(keyvals['utm']) ) )
            else:
                self.times.append(float(ovals[1]))

            if 'temp' in list(keyvals.keys()):
                self.temps.append("Y")
            else:
                self.temps.append("N")

            if 'owner' in list(keyvals.keys()):
                self.owners.append(keyvals['owner'])
            else:
                self.owners.append(None)

            if 'coverid' in list(keyvals.keys()):
                self.sheetns.append(keyvals['coversheetid'])
            else:
                if 'owner' in list(keyvals.keys()):
                    self.sheetns.append(keyvals['owner'])
                else:
                    self.sheetns.append(None)

        return

    def reverse(self):

        self.names.reverse()
        self.times.reverse()
        self.temps.reverse()
        self.owners.reverse()
        self.sheetns.reverse()

if __name__ == "__main__":
    fn = 'observed_targets.1'
    if os.path.exists(fn):
        ol = ObservedLog(filename=fn)
        print(ol)
        print(ol.names)
        print(ol.times)
        print(ol.temps)
        print(ol.owners)
        print(ol.sheetns)

        ol.reverse()
        print(ol.names)
        print(ol.times)
        print(ol.temps)
        print(ol.owners)
        print(ol.sheetns)
