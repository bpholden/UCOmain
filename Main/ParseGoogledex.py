from __future__ import print_function
from datetime import datetime, timedelta
import os
import re
import pickle
import sys
import time
import ephem
import numpy as np

import gspread
import json
from oauth2client.service_account import ServiceAccountCredentials

import ObservedLog
import Coords
from SchedulerConsts import MIN_TOTOBS, DS_BV, DS_ERR
import ExposureCalculations as ec

try:
    from apflog import *
except:
    from fake_apflog import *

def checkflag(key,didx,line,regexp,default):
    try:
        match = re.search(regexp,line[didx[key]])
        if match:
            return match.group(1)
        else:
            return default
    except:
        return default


def parseStarname(starname):

    ostarname = starname.strip()
    m= re.search("HD\s+\d+",starname)
    if m:
        ostarname = re.sub("HD\s+","HD",starname)
    m = re.search("\s+",ostarname)
    while m:
        ostarname = re.sub("\s+","_",ostarname)
        m = re.search("\s+",ostarname)
    m = re.search("\+",ostarname)
    while m:
        ostarname = re.sub("\+","p",ostarname)

        
    return ostarname


def int_or_default(value,default=0):
    try:
        attr = int(value)
    except:
        attr = default
    return attr

def float_or_default(value,default=0.0):
    try:
        rv = float(value)
    except:
        rv = default
    return rv



def parseGoogledex(sheetns=["Bstars"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn="googledex.dat",outdir=None,config={'I2': 'Y', 'decker': 'W', 'owner' : '' },force_download=False,prilim=0.5):
    """ parseGoogledex parses google sheets and returns the output as a tuple
    This routine downloads the data if needed and saves the output to a file. If the file exists, it just reads in the file.
    
    names, star_table, do_flag, stars = parseGoogledex(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn="googledex.dat")
    names - a list of stars in the starlist
    star_table - a numpy array
    flags - a dictionary of items on whether or not do="y" needs to be set for scriptobs 
    stars - a list of pyEphem objects 

    """

    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFpri", "APFcad", "APFnshots", "lastobs", "APFmin", "APFmax", \
                "B-V", "APF Desired Precision", "Close Companion", \
                "APF decker","I2", "owner", "uth","utm","duration", "Template",
                "Nobs", "Total Obs"
                ]

    
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    apflog( "Starting Googledex parse",echo=True)    
    if not outdir :
        outdir = os.getcwd()
    if os.path.exists(os.path.join(outdir,outfn)) and force_download is False:
        try:
            f = open(os.path.join(outdir,outfn),'rb')
            full_codex = pickle.load(f)
            f.close()
        except:
            full_codex = makeLocalCopy(req_cols,sheetns=sheetns,certificate=certificate,outfn=os.path.join(outdir,outfn))
    else:
        full_codex = makeLocalCopy(req_cols,sheetns=sheetns,certificate=certificate,outfn=os.path.join(outdir,outfn))

    col_names = full_codex[0]
    codex = full_codex[1:]

    didx = findColumns(col_names,req_cols)
    
    names = []
    star_table = []
    flags = { "do" : [], "decker" : [], "I2" : [], "owner" : [], "template" : [] }
    stars = []
    # Build the star table to return to 
    for ls in codex:
        row = []
        if ls[0] == '':
            continue
        apfpri = float_or_default(ls[didx["APFpri"]])
        nobs = int_or_default(ls[didx["Nobs"]])
        totobs = int_or_default(ls[didx["Total Obs"]],default=-1)

        if totobs > 0 and nobs >= totobs: continue
        if apfpri < prilim: continue
        # Get the star name
        names.append(parseStarname(ls[didx["Star Name"]]))
        
        # Get the RA
        raval = Coords.getRARad(ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]])
        if raval:
            row.append(raval)
        else:
            row.append(-1.)
        # Get the DEC
        decval = Coords.getDECRad(ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]])
        if decval:
            row.append(decval)
        else:
            row.append(-3.14)

        for coln in ("pmRA", "pmDEC"):
            row.append(float_or_default(ls[didx[coln]]))

        # Vmag
        row.append(float_or_default(ls[didx["Vmag"]],default=15.0))
        
        # For now use the old 1e9 count value - these get recalculated 
        row.append(1200.0)
        row.append(1.e9)
        # APFpri
        row.append(apfpri)
        for coln in ["APFcad","APFnshots","lastobs"] :
            row.append(float_or_default(ls[didx[coln]]))

        for coln in [ "B-V", "APF Desired Precision" ]:
            inval = float_or_default(ls[didx[coln]],default=1.0)
            if inval < 0:
                inval = 1.
            if coln is 'B-V' and inval > 2:
                inval = 1
            if coln is 'APF Desired Precision' and inval > 10:
                inval = 10
            row.append(inval)
                    
        for coln in ["uth", "utm"]:
            row.append(int_or_default(ls[didx[coln]]))
                
        # duration:
        row.append(float_or_default(ls[didx["duration"]]))
                
        # APFmin
        row.append(float_or_default(ls[didx["APFmin"]],default=MIN_TOTOBS))
                
        # APFmax
        row.append(float_or_default(ls[didx["APFmax"]]))

        # Nobs
        row.append(nobs)
                
        # Total Obs
        if totobs >= 0:
            row.append(totobs)
        else:
            row.append(0)

        if row[DS_BV] > 1.2:
            i2cnts = ec.getI2_M(row[DS_ERR])
        else:
            i2cnts = ec.getI2_K(row[DS_ERR])
        if i2cnts < 100:
            i2cnts = 100

        row.append(i2cnts)
                
        check = checkflag("Close Companion",didx,ls,"\A(y|Y)","")
        if check == "Y" or check == "y" :
            flags['do'].append(check)
        else:
            flags['do'].append("")
            
        flags['decker'].append(checkflag("APF decker",didx,ls,"\A(W|N|T|S|O|K|L|M|B)",config["decker"]))
        i2select = checkflag("I2",didx,ls,"\A(n|N)",config["I2"])
        flags['I2'].append(i2select.upper())
        tempselect = checkflag("Template",didx,ls,"\A(n|N)",'Y')
        flags['template'].append(tempselect.upper())

        flags['owner'].append(checkflag("owner",didx,ls,"\A(\w?\.?\w+)",config["owner"]))

            
        star_table.append(row)
        star = ephem.FixedBody()
        star.name = ls[0]
        star._ra = ephem.hours(str(":".join([ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]])))
        star._dec = ephem.degrees(str(":".join([ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]])))
        stars.append(star)

    return (names, np.array(star_table), flags, stars)

def updateGoogledexLastobs(filename, sheetns=["2018B"],ctime=None,certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json'):
    """
        Update the online googledex lastobs column assuming things in filename have been observed.
        updateGoogledexLastobs(filename, sheetn="The Googledex",time=None,certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json')

        filename - where the observations are logged

        returns the number of cells updated
    """
#    names, times, temps, owners = ObservedLog.getObserved(filename)
    obslog = ObservedLog.ObservedLog(filename)
    if len(obslog.names) == 0:
        return
    if ctime is None:
        ctime = datetime.utcfromtimestamp(int(time.time()))
    
    nupdates = 0
    for sheetn in sheetns:
        ws = getSpreadsheet(sheetn=sheetn,certificate=certificate)
        
        if ws:
            vals = ws.get_all_values()
        else:
            continue
        col = vals[0].index("lastobs") 
        nobscol = vals[0].index("Nobs")
        tempcol = vals[0].index("Template")
        owncol = vals[0].index("owner")
        wait_time = len(vals)
        time.sleep(wait_time)
    
        for i, v in enumerate(vals):
            # Did we observe this target tonight?
            local_name = parseStarname(v[0])
            if local_name in obslog.names:
                # We observed this target, so update the cell in the worksheet
                # update_cell(row, col, val) - col and row are 1 indexed
                nameidx = obslog.names.index(local_name)
                otime = obslog.times[nameidx]
                taketemp = obslog.temps[nameidx]
                curowner = obslog.owners[nameidx]
                if isinstance(otime,float):
                    t = datetime.utcfromtimestamp(otime)
                else:
                    hr, mn = otime
                    t = datetime(ctime.year, ctime.month, ctime.day, hr, mn)
                jd = float(ephem.julian_date(t))
                try:
                    pastdate = float(v[col])
                    try:
                        n = int(v[nobscol])
                    except:
                        n = 0
                    if jd > pastdate and curowner == v[owncol]:
                        ws.update_cell(i+1, col+1, round(jd, 4) )
                        ws.update_cell(i+1, nobscol+1, n + 1 )
                        nupdates += 2
                except:
                    print (v[0], v[col])
                    ws.update_cell(i+1, col+1, round(jd,4) )
                    ws.update_cell(i+1, nobscol+1, 1 )
                    nupdates += 2
                try:
                   have_temp = v[tempcol]
                   if taketemp == "Y" and have_temp == "N" and curowner == v[owncol]:
                       ws.update_cell(i+1, tempcol+1, "Y")
                       nupdates += 1
                except:
                    apflog( "Error logging template obs for %s" % (v[0]),echo=True,level='error')
                apflog( "Updated %s in %s" % (v[0],sheetn),echo=True)

    return nupdates

def updateLocalGoogledex(intime,googledex_file="googledex.dat", observed_file="observed_targets", frac_table=None):
    """
        Update the local copy of the googledex with the last observed star time.
        updateLocalGoogledex(time,googledex_file="googledex.dat", observed_file="observed_targets")

        opens googledex_file and inputs date of last observation from observed_file
        in principle can use timestamps as well as scriptobs uth and utm values
    """
    # names, times, temps = ObservedLog.getObserved(observed_file)
    obslog = ObservedLog.ObservedLog(observed_file)
    try:
        g = open(googledex_file, 'rb')
        full_codex = pickle.load(g)
        g.close()
    except IOError:
        apflog("googledex file did not exist, so can't be updated",echo=True)
        return obslog.names
    except EOFError:
        apflog("googledex file corrupt, so can't be updated",echo=True)
        return obslog.names

    codex_cols = full_codex[0]

    starNameIdx = codex_cols.index("Star Name")
    lastObsIdx = codex_cols.index("lastobs")
    try:
        nObsIdx = codex_cols.index("Nobs")
    except:
        nObsIdx = -1
    
    for i in range(1, len(full_codex)):
        row = full_codex[i]
        if parseStarname(row[starNameIdx]) in obslog.names:
            # We have observed this star, so lets update the last obs field
            obstime = obslog.times[obslog.names.index(row[starNameIdx])]
            if isinstance(obstime,float):
                t = datetime.utcfromtimestamp(obstime)
            else:
                hr, min = obstime
                if type(intime) != datetime:
                    ctime = datetime.now()
                    td = timedelta(0,3600.*7)
                    intime = ctime + td
                t = datetime(intime.year, intime.month, intime.day, hr, min)


            jd = round(float(ephem.julian_date(t)), 4) 
            apflog( "Updating local googledex star %s from time %s to %s" % (row[starNameIdx], row[lastObsIdx], str(jd)),echo=True)
            row[lastObsIdx] = str(jd)
            if nObsIdx > 0:
                try:
                    row[nObsIdx] = int(row[nObsIdx]) + 1
                except:
                    row[nObsIdx] = 1
            full_codex[i] = row

    with open(googledex_file, 'wb') as f:
        pickle.dump(full_codex, f)
    f.close()
    
    return obslog.names

def makeLocalCopy(req_cols,sheetns=["The Googledex"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn="./googledex.dat"):
    full_codex = []
    # These are the columns we need for scheduling
    hdr_cols = [rc for rc in req_cols]
    hdr_cols.append("Sheetname")
    full_codex.append(hdr_cols)
        
    for sheetn in sheetns:
        worksheet = getSpreadsheet(sheetn=sheetn,certificate=certificate)
        if worksheet:
            cur_codex = worksheet.get_all_values()
            didx = findColumns(cur_codex[0],req_cols)
            
            for row in cur_codex[1:]:
                nrow = []
                for c in req_cols:
                    nrow.append(row[didx[c]])
                nrow.append(sheetn)
                full_codex.append(nrow)

            wait_time = len(nrow)
            time.sleep(wait_time)
        else:
            time.sleep(10)

    f = open(outfn,'wb')
    pickle.dump(full_codex, f)
    f.close()
    return full_codex
    
def getSpreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json'):
    """ Get the spreadsheet from google

    worksheet = getSpreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json')
    worksheet - the worksheet object returned by the gspread module

    sheetn - name of the google sheet, defaults to "The Googledex"
    certificate - certificate used to control access to the google sheet
    
    """
    # this downloads the googledex from the Google Drive
    # the certificate must be available
    # these certificates are generated through the Google Developer Interface
    # the developer must select the correct API for access

    # the certificate has an email associated with it, that email must
    # have the document shared with it to allow access 

    certificate_path = os.path.dirname(__file__)
    
    json_key = json.load(open(os.path.join(certificate_path, certificate)))
    scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']

    credentials = ServiceAccountCredentials.from_json_keyfile_name(os.path.join(certificate_path, certificate), scope)
    try:
        gs = gspread.authorize(credentials)
        apflog("Successfully logged in.", echo=True)
    except:
        apflog("Cannot log into Google API.", echo=True,level='error')
        return None
    apflog("Attempting to Open %s" % (sheetn),echo=True)
    try:
        spreadsheet = gs.open(sheetn)
        apflog("Loaded Main %s" % (sheetn),echo=True)
        worksheet = spreadsheet.get_worksheet(0)
        apflog("Got spreadsheet", echo=True)
    except Exception as e:
        apflog("Cannot Read %s: %s"  % (sheetn, e), echo=True, level='error')
        worksheet = None
    return worksheet


def findColumns(col_names,req_cols,opt_cols=[]):
    """ findColumns finds the indices for the column names in the list of required columns
    indices = findColumns(col_names, req_cols)
    
    indices - a list of indices, each index maps to where in col_names the column is found and in the order of req_cols
    col_names - list of column names to be searched
    req_cols - list of names that should be in the first list
    """
    idx = []
    didx = dict()

    for r in req_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
        else:
            apflog("%s Not found in column names from google spreadsheet" % (r) , level="Warn",echo=True)

    for r in opt_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
            
    # hack to handle an error
    if req_cols[0] == "Star Name" and req_cols[0] not in didx.keys():
        didx[req_cols[0]] = 0
        apflog("Pasting 'Star Name' into column 0 of google spreadsheet" , level="Warn",echo=True)

    return didx

if __name__ == '__main__':


    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFpri", "APFcad", "APFnshots", "lastobs", "APFmin", "APFmax", \
                "B-V", "APF Desired Precision", "Close Companion", \
                "APF decker","I2", "owner", "uth","utm","duration", "Template",
                "Nobs", "Total Obs"
                ]

    sheetnstr="Bstars,A003_PRobertson_2019B,A006_PDalba_2019B,A007_HIsaacson_2019B,A009_MKosiarek_2019B,A011_SKane_2019B,A012_SKane_2019B,A015_AHoward_2019B,A013_ASiemion_2019B,A000_BWelsh_2019B,A001_ICzekala_2019B,A002_ICzekala_2019B,A004_PRobertson_2019B,A007_HIsaacson_2019B,A008_BHolden_2019B,A014_SVogt_2019B,A015_TBrandt_2019B"
    sheetns = sheetnstr.split(",")
    fc = makeLocalCopy(req_cols, sheetns=sheetns)

    names, star_table, flags, stars = parseGoogledex(sheetns=sheetns)
    
    
