from __future__ import print_function
import sys

from optparse import OptionParser

sys.path.append("../master")
import ParseGoogledex

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-o","--owner",dest="owner",default="S.Vogt")
    parser.add_option("-s","--sheet",dest="sheet",default="2017B")
    (options, args) = parser.parse_args()    


    (names, star_table, flags, stars) = ParseGoogledex.parseGoogledex(sheetns=options.sheet)
     
    for i in range(0,len(names)):

        if options.owner in flags['owner'][i]:
            print (names[i])
