from __future__ import print_function
import sys
import os

from astropy.io import ascii
from astropy.table import Table
import numpy as np

sys.path.append("../master")
import ParseGoogledex
import SchedulerConsts as sc
from ExposureCalculations import getI2_M, getI2_K, getEXPMeter

if __name__ == "__main__":

    req_cols = ["Star Name","APFpri", "RA hr", "RA min", "RA sec", \
                    "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", "B-V", \
                    "texp", "APFnshots", "expcount","I2", "decker","Close Companion",  \
                    "APFcad", "lastobs", "Nobs", "Total Obs", "Template", "owner",\
                    "uth","utm","duration", \
                     "mode", "raoff", "decoff", "obsblock",
                    ]
    

    
    if len(sys.argv) <= 1:
        print ("needs a sheet list")
        sys.exit()
    sheetns = sys.argv[1].split(",")
    for sheetn in sheetns:
        if os.path.exists("googledex.dat"):
            os.unlink("googledex.dat")
        names, star_table, flags, stars = ParseGoogledex.parseGoogledex(sheetns=[sheetn],prilim=-100)

        # translate - make a table

        outdict = dict()
        outdict["Star Name"] = np.asarray(names)
        
        outdict["pmRA"]= star_table[:,sc.DS_PMRA]
        outdict["pmDEC"] = star_table[:,sc.DS_PMDEC]
        outdict["Vmag"] = star_table[:,sc.DS_VMAG]

        outdict["APFpri"] = star_table[:,sc.DS_APFPRI]
        outdict["APFcad"] = star_table[:,sc.DS_CAD]
        outdict["APFnshots"] = star_table[:,sc.DS_NSHOTS]
        low = outdict["APFnshots"] <= 0
        outdict["APFnshots"][low] = 1
        outdict["B-V"] = star_table[:,sc.DS_BV]

        outdict["decker"] = np.asarray(flags["decker"])
        outdict["I2"] = np.asarray(flags["I2"])
        outdict["Close Companion"] = np.asarray(flags["do"])

        outdict["lastobs"] = star_table[:,sc.DS_LAST]
        outdict["Nobs"] = star_table[:,sc.DS_NOB]
        outdict["Total Obs"] = star_table[:,sc.DS_TOT]
    
        outdict["owner"] = np.asarray(flags["owner"])
        outdict["Template"] = np.asarray(flags["template"])

        outdict["uth"] = star_table[:,sc.DS_UTH]
        outdict["utm"] = star_table[:,sc.DS_UTM]
        outdict["duration"] = star_table[:,sc.DS_DUR]
        outdict["mode"] = np.zeros_like(np.asarray(flags["do"]))
        outdict["obsblock"] = np.zeros_like(np.asarray(flags["do"]))                                  
        outdict["raoff"] = np.zeros_like(np.asarray(flags["do"]))                                  
        outdict["decoff"] = np.zeros_like(np.asarray(flags["do"]))                                  


        rah = []
        ram = []
        ras = []
        decd= []
        decm=[]
        decs=[]
        for star in stars:
            crah,cram,cras = str(star._ra).split(":")
            cdecd,cdecm,cdecs = str(star._dec).split(":")
            rah.append(crah)
            ram.append(cram)
            ras.append(cras)
            decd.append(cdecd)
            decm.append(cdecm)
            decs.append(cdecs)
        outdict['RA hr'] = np.asarray(rah)
        outdict['RA min'] = np.asarray(ram)
        outdict['RA sec'] = np.asarray(ras)
        outdict['Dec deg'] = np.asarray(decd)
        outdict['Dec min'] = np.asarray(decm)
        outdict['Dec sec'] = np.asarray(decs)        

        outdict["texp"] = star_table[:,sc.DS_EXPT]

        are_m = outdict['B-V'] > 1.2
        are_k = np.logical_not(are_m)
        
        i2cnts = np.zeros_like(star_table[:,sc.DS_EXPT])
        i2cnts[are_m] = getI2_M(star_table[:,sc.DS_ERR][are_m])
        i2cnts[are_k] = getI2_K(star_table[:,sc.DS_ERR][are_k])

        outdict['expcount'] =  (getEXPMeter(i2cnts,star_table[:,sc.DS_BV]))

        big = outdict['expcount']>1e7
        small = outdict['expcount']<1e7
        
        outdict['expcount'][big] = np.around(outdict['expcount'][big]/1e7)*1e7
        outdict['expcount'][small] = np.around(outdict['expcount'][small]/1e6)*1e6        
        
        outtable = Table(outdict)
        outfn = sheetn + ".csv"
        ascii.write(outtable[req_cols],output=outfn,format='csv')
