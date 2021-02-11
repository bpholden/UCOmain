import pickle
import csv
from astropy.time import Time
from datetime import datetime


def readin_csv(fn):

    cat = dict()
    with open(fn) as cfp:
        rdr = csv.DictReader(cfp)
        for rw in rdr:
            for k in rw.keys():
                if k not in cat.keys():
                    cat[k] = [rw[k]]
                else:
                    cat[k].append(rw[k])
    cfp.close()
    return cat

def makenames(keylist):
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFpri", "APFcad", "APFnshots", "lastobs", "APFmin", "APFmax", \
                "B-V", "APF Desired Precision", "Close Companion", \
                "APF decker","I2", "owner", "uth","utm","duration"
                ]
    return req_cols

def conv_lastobs(lastobs):
    new = []
    t = Time(datetime.now())
    cdj = t.jd
    for l in lastobs:
        if l == "inf":
            new.append( str(0.0))
        else:
            new.append( "%.2f" % (cdj - float(l)))
    return new

fn ="bj_targets_queue.csv"
ofn ="bj_targets_queue.dat"

cat = readin_csv(fn)
nentries = len(cat['Star Name'])
tab = []

cat['lastobs'] = conv_lastobs(cat['lastobs'])
names = makenames(cat.keys())
tab.append(names)
for n in range(0,nentries):
    srow = []
    for k in names:
        if k == "lastobs":
            cat[k][n] = str(0.)
        srow.append(cat[k][n])
    tab.append(srow)

for n in range(0,nentries):
    row = tab[n]
    for l in range(0,len(row)):
        if type(row[l]) != str:
            print( n,l,type(row[l]))

    
fp = open(ofn,"wb+")
pickle.dump(tab,fp)
fp.close()
