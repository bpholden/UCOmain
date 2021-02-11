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
        times - times either as a timestamp in second column or a (hour,minute) tuple from a scriptobs line
        temps - a list of template observations
        owners - a list of owners
        sheetns - a list of names for the sheets
        """

        try:
            f = open(self.filename, 'r')
        except IOError:
            apflog( "Couldn't open %s" % filename,level="warn",echo=True)
            return 
        else: 
            for line in f:
                line = line.strip()
                if len(line) > 0:
                    if line[0] == '#' or line == "":
                        pass
                    else:
                        ovals, keyvals = self.parse_key_vals(line)
                        self.names.append(ovals[0])
                        if 'uth' in keyvals.keys():
                            self.times.append( ( int(keyvals['uth']), int(keyvals['utm']) ) )
                        else:
                            self.times.append(float(ovals[1]))

                        
                        if 'temp' in keyvals.keys():
                            self.temps.append("Y")
                        else:
                            self.temps.append("N")
                            
                        if 'owner' in keyvals.keys():
                            self.owners.append(keyvals['owner'])
                        else:
                            self.owners.append(None)

                        if 'coverid' in keyvals.keys():
                            self.sheetns.append(keyvals['coversheetid'])
                        else:
                            if 'owner' in keyvals.keys():
                                self.sheetns.append(keyvals['owner'])
                            else:
                                self.sheetns.append(None)
            
        self.names.reverse()
        self.times.reverse()
        self.temps.reverse()
        self.owners.reverse()
        self.sheetns.reverse()
        return 
        

def getObserved(filename):
    """ getObserved parses a file to find the object names and times
    names, times, temps = getObserved(filename)
    names - list of names, must be first column of file called filename
    times - times either as a timestamp in second column or a (hour,minute) tuple from a scriptobs line
    temps - a list of template observations

    """
    names = []
    times = []
    temps = []
    nobs = dict()
    try:
        f = open(filename, 'r')
    except IOError:
        apflog( "Couldn't open %s" % filename,level="warn",echo=True)
        return names, times, temps
    else: 
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] == '#' or line == "":
                    pass
                else:
                    ls = line.split()
                    names.append(ls[0])
                    if len(ls) > 15:
                        times.append( (int(ls[14].split('=')[1]), int(ls[15].split('=')[1])) )
                    else:
                        times.append(float(ls[1]))
                    length = len(ls)
                    mtch = re.search("temp\=Y",ls[length-1])
                    if mtch:
                        temps.append("Y")
                    else:
                        temps.append("N")

            
    names.reverse()
    times.reverse()
    temps.reverse()
    return names, times, temps
	

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
