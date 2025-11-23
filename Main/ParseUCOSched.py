from __future__ import print_function

import json
import os
import re
import time
import datetime

import astropy
import astropy.io.ascii
import astropy.table
import numpy as np
import ephem
import gspread
from oauth2client.service_account import ServiceAccountCredentials

import Coords
import ObservedLog

from SchedulerConsts import EXP_LIM, MAX_PRI, DEFAULT_CERT

try:
    from apflog import *
except:
    from fake_apflog import *

def check_flag(key,didx,line,regexp,default):
    """ check_flag(key, dict_ind, line, regexp, default)

    - key : the key in the dictionary indices, which key from the input list you want to check
    - didx : dictionary of indices so that keys can be used instead of indices for the line
    - line : a list of entries, map to keys by didx
    - regexp : if the regexp is matched, that value is returned as a string
    - default : value if the regexp fails

    """

    try:
        match = re.search(regexp,line[didx[key]])
        if match:
            return match.group(1)
        else:
            return default
    except:
        return default


def read_star_table(table_filename):
    '''
    star_table = read_star_table(table_filename)

    table_filename - name of the file to read in
    returns an astropy table object
    sets the sheetn and owner columns to empty strings if they are None
    '''
    star_table = astropy.io.ascii.read(table_filename)

#    for coln in ('mode','obsblock','raoff','decoff','sheetn','owner'):
    for coln in ('sheetn','owner'):
        try:
            star_table[coln][star_table[coln] == 'None'] = ''
        except:
            pass

    return star_table



def parse_starname(starname):
    """parse_starname(starname)

    starname - input value which should be the name of the a star, duh
    returns the starname value, doing some clean up to meet
    requirements (no spaces, trim trialing spaces, etc.)
    """
    ostarname = starname.strip()
    m= re.search(r"HD\s+\d+",starname)
    if m:
        ostarname = re.sub(r"HD\s+","HD",starname)
    m = re.search(r"\s+",ostarname)
    while m:
        ostarname = re.sub(r"\s+","_",ostarname)
        m = re.search("\s+",ostarname)
    m = re.search(r"\+",ostarname)
    while m:
        ostarname = re.sub(r"\+","p",ostarname)
        m = re.search(r"\+",ostarname)

    return ostarname


def int_default(value,default=0):
    """
    int_default(value,default=0)
    returns the input value as an integer, and a failure to cast returns default
    """
    try:
        attr = int(value)
    except:
        attr = default
    return attr

def float_default(value,default=0.0):
    """
    float_default(value,default=0.0)
    returns the input value as an float, and a failure to cast returns default
    """
    try:
        rv = float(value)
    except:
        rv = default
    return rv


def get_spreadsheet(sheetn="The Googledex", certificate=DEFAULT_CERT):
    """ Get the spreadsheet from google

    worksheet = get_spreadsheet(sheetn="The Googledex",certificate=DEFAULT_CERT)
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

    certificate_path = os.path.dirname("/usr/local/lick/data/apf/master/")
    if os.path.exists(certificate_path) is False:
        certificate_path = os.path.dirname(__file__)
    finpath = os.path.join(certificate_path, certificate)

    _ = json.load(open(finpath))
    scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']

    credentials = ServiceAccountCredentials.from_json_keyfile_name(os.path.join(certificate_path, certificate), scope)
    try:
        gs = gspread.authorize(credentials)
        apflog("Successfully logged in.", echo=True)
    except Exception as e:
        apflog("Cannot log into Google API.", echo=True,level='error')
        apflog("%s %s" % (type(e), e),echo=True,level='error')
        return None, None
    worksheet = None
    tries = 0
    errlog = None
    while worksheet is None and tries < 3:
        tries = tries + 1
        try:
            spreadsheet = gs.open(sheetn)
            apflog("Loaded %s" % (sheetn),echo=True)
            worksheet = spreadsheet.sheet1
            apflog("Got spreadsheet.sheet1", echo=True)
        except Exception as e:
            errlog = "Cannot Read %s: %s %s" % (sheetn, type(e), e)
            time.sleep(1)
    if worksheet is None:
        apflog(errlog,echo=True,level='error')
        return None, None

    cur_codex = None
    more_sleeping=10.
    while cur_codex is None:
        try:
            cur_codex = worksheet.get_all_values()
        except:
            time.sleep(more_sleeping)
            cur_codex = None

    return worksheet, cur_codex

def retrieve_codex(req_cols,sheetns, certificate=DEFAULT_CERT, sleep=True):
    """retrieve_codex(req_cols,sheetns=["The Googledex"],certificate=DEFAULT_CERT)

    returns the "codex", a list of lists containing all of the columns
    in the req_cols list, source of the data are the Google sheets named
    in sheetns, needs a certificate to authenticate.

    - req_cols : a list of column names in the sheets that are required
    for the final list of lists
    - sheetns : sheets to download from
    - certificate : the thing that allows authentication
    """
    full_codex = []
    # These are the columns we need for scheduling
    full_codex.append(req_cols)
    for sheetn in sheetns:
        wait_time = 0
        worksheet, cur_codex = get_spreadsheet(sheetn=sheetn, certificate=certificate)
        if worksheet:
            if len(cur_codex) <= 0:
                apflog("Worksheet %s exists but is empty, skipping" % (sheetn), \
                       level='error', echo=True)
                continue
            didx = find_columns(cur_codex[0],req_cols)

            for row in cur_codex[1:]:
                nrow = []
                for c in req_cols:
                    if c in list(didx.keys()):
                        nrow.append(row[didx[c]])
                    else:
                        if c == 'sheetn':
                            nrow.append(sheetn)
                        else:
                            nrow.append(None)

                full_codex.append(nrow)
                wait_time += .3
            if sleep and ((sheetns.index(sheetn)+1) < len(sheetns)):
                apflog("Sleeping %.1f seconds to keep Google happy" % (wait_time), \
                       level="info",echo=True)
                time.sleep(wait_time)

    return full_codex


def find_columns(col_names,req_cols):
    """find_columns finds the indices for the column names in the list of

    required columns indices = find_columns(col_names, req_cols)

    indices - a dictionary of indices, each index maps to where in
    col_names the column is found

    col_names - list of column names to be searched
    req_cols - list of names that should be in the first list
    """

    didx = dict()

    for r in req_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
        else:
            apflog("%s Not found in column names from google spreadsheet" % (r) ,\
                   level="Warn",echo=True)

    # hack to handle an error
    if req_cols[0] == "Star Name" and req_cols[0] not in list(didx.keys()):
        didx[req_cols[0]] = 0
        apflog("Pasting 'Star Name' into column 0 of google spreadsheet" , level="Warn",echo=True)

    return didx


def parse_rank_table(sheet_table_name='2022A_ranks',certificate=DEFAULT_CERT):
    '''
    rank_table = parse_rank_table(sheet_table_name='2022A_ranks',certificate=DEFAULT_CERT)

    sheet_table_name - name of the google sheet to download containing the rank table
    certificate - json file for authenticating to access google sheets

    returns an astropy table object of the rank table
    '''
    apflog( "Starting parse of %s" % (sheet_table_name),echo=True)

    sheetns = []
    rank = []
    frac = []
    too = []

    worksheet, cur_codex = get_spreadsheet(sheetn=sheet_table_name,certificate=certificate)
    if worksheet is None:
        return None, None, None, None

    req_cols = ["sheetn", "rank", "frac", "too"]
    didx = find_columns(cur_codex[0], req_cols)

    for row in cur_codex[1:]:
        if row[0] != "":
            sheetns.append(row[didx['sheetn']])
            crank = float_default(row[didx['rank']])
            crank = int(round(crank))
            rank.append(crank)
            cfrac = float_default(row[didx['frac']])
            frac.append(cfrac)
            too.append(row[didx['too']].lower())

    if 'RECUR_A100' not in sheetns:
        sheetns.append('RECUR_A100')
        rank.append(200)
        frac.append(0.05)
        too.append('n')

    return sheetns, rank, frac, too

def init_star_table(col_list):

    """
    star_table = init_star_table(column_list)
    star_table - a Astropy Table object that has the columns needed, most are in column_list
    forces certain columns to be added

    """

    star_table = dict()
    for col in col_list:
        star_table[col] = []
    star_table['name'] = []
    star_table['do'] = []
    star_table['nobs'] = []
    star_table['totobs'] = []
    star_table['ra'] = []
    star_table['dec'] = []
    star_table['nexp'] = []
    star_table['cad'] = []
    star_table['pri'] = []
    star_table['moon'] = []
    star_table['sheetn'] = []
    star_table['only_template'] = []

    return star_table

def parse_codex(config,sheetns=["RECUR_A100"],certificate=DEFAULT_CERT,prilim=1,sleep=True,
                hour_constraints=None):
    '''
    star_table = parse_codex(config,sheetns=["RECUR_A100"],certificate=DEFAULT_CERT,
                             prilim=1,sleep=True,hour_constraints=None)

    config - dictionary of default values for a number of flags
    sheetns - list of google sheet names
    certificate - json file for authenticating to access google sheets
    prilim - limit on priority values, values below this are tossed
    sleep - sleep between downloading sheets
    hour_constraints - astropy table of hour constraints, if None, no constraints are applied
    '''
    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                    "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                    "texp", "I2", "expcount", "decker","Close Companion", \
                    "lastobs", "B-V", \
                    "cad", "pri", "nexp", "count", "binning", \
                    "night_cad", "night_obs", "night_nexp", "DaysNew", \
                    "Template", "Nobs", "Total Obs", "Bstar", "Only Template", \
#                    "mode", "raoff", "decoff",  "obsblock",\
                    "need_cal", "cal_star",
                    'sheetn', 'owner' \
                    ]

    full_codex = retrieve_codex(req_cols, sheetns=sheetns, certificate=certificate, sleep=sleep)
    if full_codex is None or len(full_codex) == 1:
        return None

    # the assumption is that the first row is the column names
    col_names = full_codex[0]
    codex = full_codex[1:]

    didx = find_columns(col_names, req_cols)
    star_table = init_star_table(req_cols)

    if hour_constraints is not None:
        done_names = hour_constraints['runname'][hour_constraints['left'] < 0]
    else:
        done_names = []

    # Go line by line through the list of lists that represents the google sheet
    for ls in codex:
        if ls[0] == '':
            continue
        # Get the priority, we get this early to avoid parsing the rest of the line
        # if the priority is too low
        if "pri" in didx and ls[didx["pri"]] is not None:
            apfpri = int_default(ls[didx["pri"]],default=-1)
        else:
            apfpri = int_default(ls[didx["APFpri"]],default=-1)

        # if the nobs >= total obs, we are done with this target
        nobs = int_default(ls[didx["Nobs"]])
        totobs = int_default(ls[didx["Total Obs"]],default=-1)

        # if the sheet name is in the done list, we will skip this target
        csheetn = check_flag("sheetn",didx,ls,r"\A(.*)", 'RECUR_A100')
        if 'owner' in didx:
            owner = check_flag("owner",didx,ls,r"\A(.*)", csheetn)

        # these are all the conditions that will cause us to skip this target
        if totobs > 0 and nobs >= totobs: continue
        if apfpri < prilim: continue
        if csheetn in done_names: continue

        if apfpri > MAX_PRI: apfpri = MAX_PRI

        # star names that are allowed in Google sheets may not work
        # in the actual star list sent to scroptobs
        # so we do some cleanup here
        name = parse_starname(ls[didx["Star Name"]])

        # Get the RA
        raval, rahr, ramin, rasec = Coords.get_RA_rad(ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]])
        if raval is None:
            # alarm
            apflog("Error in RA coordinates for %s: %s" %(name, csheetn),level='warn',echo=True)
            continue

        # Get the DEC
        decval, decdeg, decmin, decsec = Coords.get_dec_rad(ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]])
        if decval is None:
            # alarm
            apflog("Error in Dec coordinates for %s: %s" %(name, csheetn),level='warn',echo=True)
            continue

        # why are we doing this you may ask?
        # we use Google sheets which cannot have -0 for a value
        # but if we pass a value like 00:-16:00 to eostele, it generates
        # an incorrect declination value
        # so, we move the - to the front of the sexagesimal string
        # the radian values above are only used for the scheduler, we still
        # command the telescope in the raw units


        if name and raval and decval:
            star_table["name"].append(name)

            star_table['ra'].append(raval)
            star_table['RA hr'].append(rahr)
            star_table['RA min'].append(ramin)
            star_table['RA sec'].append(rasec)

            star_table['dec'].append(decval)
            star_table["Dec deg"].append(decdeg)
            star_table["Dec min"].append(decmin)
            star_table["Dec sec"].append(decsec)
        else:
            continue

#        mode = check_flag("mode",didx,ls,"\A(b|B|a|A|c|C)",config["mode"])
#        if type(mode) == str:
#            mode = mode.upper()
#        star_table['mode'].append(mode)
#        star_table['raoff'].append(check_flag("raoff",didx,ls,"\A((\+|\-)?\d+\.?\d*)",config["raoff"]))
#        star_table['decoff'].append(check_flag("decoff",didx,ls,"\A((\+|\-)?\d+\.?\d*)",config["decoff"]))

        for coln in ("pmRA", "pmDEC"):
            star_table[coln].append(float_default(ls[didx[coln]]))


        star_table['Vmag'].append(float_default(ls[didx["Vmag"]], default=10.0))
        star_table['texp'].append(float_default(ls[didx["texp"]], default=1200))
        expcount = float_default(ls[didx["expcount"]],default=1e9)
        expcount = min(expcount,EXP_LIM)

        star_table['expcount'].append(expcount)
        if "nexp" in didx and ls[didx["nexp"]] is not None:
            star_table['nexp'].append(int_default(ls[didx["nexp"]], default=1))
        elif "count" in didx and ls[didx['count']] is not None:
            star_table['nexp'].append(int_default(ls[didx["count"]], default=1))
        else:
            star_table['nexp'].append(int_default(ls[didx["APFnshots"]], default=1))

        # most programs do not have binning, so we default to 1,1
        if "binning" in didx and ls[didx["binning"]] is not None:
            binp = int_default(ls[didx["binning"]], default=1)
            if binp != 1 and binp != 2 and binp != 4 and binp != 8:
                binp=1
            bin_str = "%d,%d" % (binp, binp)
            star_table['binning'].append(bin_str)
        else:
            star_table['binning'].append("1,1")

        # some targets want dark time, but many targets do not
        # the external representation is the number of days new, which is
        # is what the coversheet has, internal representation is the fraction
        # of the moon that is illuminated
        if "DaysNew" in didx and ls[didx["DaysNew"]] is not None:
            days_from_new = float_default(ls[didx['DaysNew']], default=15.0)
            star_table['moon'].append(days_from_new / 15.0)
        else:
            star_table['moon'].append(1.0)

        # scheduler specific, cadence is the time between exposures
        # in nights, there are two different allowed column headings
        if "cad" in didx and ls[didx['cad']] is not None:
            cad_value = float_default(ls[didx["cad"]], default=0.7)
        else:
            cad_value = float_default(ls[didx["APFcad"]], default=0.7)

        if cad_value > 0:
            star_table['cad'].append(cad_value)
        else:
            star_table['cad'].append(0.7)

        # night_cad is how often per night a target can be observed
        # night_obs is how many times it has been observed in the current night
        night_cad = float_default(ls[didx["night_cad"]], default=-1.0)
        night_nexp = 1
        if night_cad > 0:
            night_cad /= 60*24
            night_nexp = float_default(ls[didx["night_nexp"]], default=2)
        star_table['night_cad'].append(night_cad)
        star_table['night_obs'].append(0)
        star_table['night_nexp'].append(night_nexp)

        # cal_star and need_cal are flags for whether or not a target is a calibrator
        # or if the star needs to have a calibrator observed the same night
        cal_star_val = check_flag("cal_star", didx, ls, r"\A(y|Y)","N")
        need_cal_val = check_flag("need_cal", didx, ls, r"\A(y|Y)","N")
        star_table['cal_star'].append(cal_star_val.upper())
        star_table['need_cal'].append(need_cal_val.upper())

        star_table['pri'].append(apfpri)
        star_table["lastobs"].append(float_default(ls[didx["lastobs"]], default=0))

        inval = float_default(ls[didx["B-V"]], default=0.7)
        if inval < 0:
            inval = 1.
        if coln == 'B-V' and inval > 2:
            inval = 1
        star_table['B-V'].append(inval)

        # Nobs - number of observations
        star_table['nobs'].append(nobs)

        # Total Obs
        if totobs >= 0:
            star_table['totobs'].append(totobs)
        else:
            star_table['totobs'].append(0)

        # another case where the column name is not consistent with the scriptobs name
        check = check_flag("Close Companion", didx, ls, r"\A(y|Y)","")
        if "do" in didx:
            check = check_flag("do", didx, ls, r"\A(y|Y)","")
        if check == "Y" or check == "y" :
            star_table['do'].append(check)
        else:
            star_table['do'].append("")

        # details about the spectrometer configuration
        decker_name = check_flag("decker", didx, ls, r"\A(W|N|T|S|O|K|L|M|B)", config["decker"])
        star_table['decker'].append(decker_name)
        i2select = check_flag("I2", didx, ls, r"\A(n|N)", config["I2"])
        star_table['I2'].append(i2select.upper())
        tempselect = check_flag("Template",didx,ls,r"\A(n|N)",'Y')
        star_table['Template'].append(tempselect.upper())
        tempselect = check_flag("Only Template",didx,ls,r"\A(y|Y)",'N')
        star_table['only_template'].append(tempselect.upper())
#        star_table['obsblock'].append(check_flag("obsblock",didx,ls,"\A(\w+)",config["obsblock"]))
#        star_table['inst'].append(check_flag("inst",didx,ls,"(levy|darts)",config['inst']).lower())

        # need to check raoff and decoff values and alarm on failure

        # a Bstar is a specific calibration star, so has its own flag
        if 'Bstar' in didx:
            star_table['Bstar'].append(check_flag('Bstar', didx, ls, "(Y|y)", 'N'))
            star_table['sheetn'].append(csheetn)
        else:
            if 'RECUR_A100' in csheetn :
                star_table['Bstar'].append("Y")
                star_table['sheetn'].append('RECUR_A100')
            else:
                star_table['Bstar'].append("N")
                star_table['sheetn'].append(csheetn)

    badkeylist = []
    for k in list(star_table.keys()):
        if len(star_table[k]) == 0:
            badkeylist.append(k)
    for k in badkeylist:
        del star_table[k]

    # This just reorders the columns
    # This way the ascii table has these columns in front to make finding targets by specific programs easier
    star_table_names = list(star_table.keys())
    for n in  ('Dec sec', 'Dec min', 'Dec deg', 'RA sec', 'RA min', 'RA hr', 'sheetn', 'name'):
        if n in star_table_names:
            star_table_names.remove(n)
            star_table_names = [n] + star_table_names

    star_table = astropy.table.Table(star_table, names=star_table_names)

    return star_table

def gen_stars(star_table):
    """
    pyephem_objs = gen_stars(star_table)

    given a star_table returned by parse_codex (or init_star_table) returns
    a list of pyephem objects for every object in the table

    Inputs star_table - astropy Table that must have the RA and Dec in
    sexigresimal format with each column for each part of the
    coordinates separate
    """
    stars = []
    if 'name' in star_table.colnames:
        for i in range(0,len(star_table['name'])):
            star = ephem.FixedBody()
            star.name = star_table['name'][i]
            star._ra = ephem.hours(":".join([str(star_table["RA hr"][i]), str(star_table["RA min"][i]), str(star_table["RA sec"][i])]))
            star._dec = ephem.degrees(":".join([str(star_table["Dec deg"][i]), str(star_table["Dec min"][i]), str(star_table["Dec sec"][i])]))
            stars.append(star)

    return stars



def parse_UCOSched(sheetns=["RECUR_A100"], certificate=DEFAULT_CERT, outfn="sched.dat",
                   outdir=None, config={'I2': 'Y', 'decker': 'W',  'Bstar' : 'N' },
                   force_download=False, prilim=0.5, hour_constraints=None):
    """ parse_UCOSched parses google sheets and returns the output as a tuple
    This routine downloads the data if needed and saves the output to a file.
    If the file exists, it just reads in the file.

    star_table, stars = parse_UCOSched(sheetns=["RECUR_A100"],certificate='cert.json',outfn="sched.dat",outdir=None,config={'I2': 'Y', 'decker': 'W', 'Bstar' : 'N' },force_download=False,prilim=0.5)

    star_table - an astropy table
    stars - a list of pyEphem objects

    Inputs:
    sheetns - list of google sheet names
    certificate - json file for authenticating to access google sheets
    outfn - output file name, will read this in if it already exists
            instead of downloading sheets if force_download is False
    outdir - output directory for outfn, defaults to ./
    config - default values for a number of flags
    force_download - force the google sheets to be downloaded even if outfn already exists
    prilim - limit on priority values, values below this are tossed

    """


    stars = None
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    if not outdir :
        outdir = os.getcwd()
    outfn = os.path.join(outdir,outfn)
    if os.path.exists(outfn) and force_download is False:
        try:
            star_table = read_star_table(outfn)
        except:
            star_table  = parse_codex(config, sheetns=sheetns, certificate=certificate, prilim=prilim, 
                                      hour_constraints=hour_constraints)

    else:
        star_table = parse_codex(config, sheetns=sheetns, certificate=certificate, prilim=prilim, 
                                 hour_constraints=hour_constraints)

    if star_table is None:
        return (star_table, stars)

    stars = gen_stars(star_table)

    if len(stars) > 0:
        astropy.io.ascii.write(star_table, outfn, format='ecsv', overwrite=True)
    else:
        star_table = None
        stars = None

    return (star_table, stars)


def parse_TOO(too_sheetns=None, outfn='googledex.dat', outdir=None, certificate=DEFAULT_CERT, prilim=0.5):
    '''
    star_table = parse_TOO(too_sheetns=None, outfn='googledex.dat', outdir=None, certificate=DEFAULT_CERT, prilim=0.5)

    too_sheetns - list of google sheet names
    outfn - output file name, will read this in if it already exists
    outdir - output directory for outfn, defaults to ./
    certificate - json file for authenticating to access google sheets
    prilim - limit on priority values, values below this are tossed

    writes the output to outfn

    this is a wrapper for parse_codex that updates the googledex with the ToO targets
    the too_sheetns are the names of the google sheets that contain the ToO targets,
    should not be all sheet lists as that will take too long to download
    '''
    if not outdir :
        outdir = os.getcwd()

    if too_sheetns is None:
        return

    outfn = os.path.join(outdir,outfn)
    if os.path.exists(outfn) is False:
        return

    try:
        star_table = read_star_table(outfn)
    except:
        return

#    config={'I2': 'Y', 'decker': 'W', 'owner' : '', 'mode' : '', 'obsblock' : '', 'Bstar' : 'N' , 'raoff' : None, 'decoff' : None }
    config={'I2': 'N', 'decker': 'W', 'owner' : '',  'Bstar' : 'N'  }
    too_table = parse_codex(config, sheetns=too_sheetns, certificate=certificate, \
                            prilim=prilim, sleep=False)

    for n in too_sheetns:
        cur = star_table['sheetn'] == n
        if np.any(cur):
            # entries already exist, we will delete them
            # this ensures that all changes are propogated
            star_table.remove_rows(cur)
    # now just append new values
    star_table = astropy.table.vstack([star_table, too_table])

    # write to the googledex
    astropy.io.ascii.write(star_table,outfn, format='ecsv', overwrite=True)




def update_local_starlist(intime, observed_file="observed_targets", outfn='parsesched.dat', toofn='too.dat', outdir=None):
    """
        Update the local copy of the googledex with the last observed star time.

        update_local_starlist(time,googledex_file="googledex.dat", observed_file="observed_targets")

        opens googledex_file and inputs date of last observation from observed_file
        in principle can use timestamps as well as scriptobs uth and utm values
    """

    if not outdir :
        outdir = os.getcwd()

    obslog = ObservedLog.ObservedLog(filename=os.path.join(outdir, observed_file))

    outfn = os.path.join(outdir, outfn)
    if os.path.exists(outfn):
        star_table = read_star_table(outfn)
    else:
        return obslog, None

    toofn = os.path.join(outdir, toofn)
    if os.path.exists(toofn):
        too_table = read_star_table(toofn)
    else:
        too_table = None

    for name in obslog.names:

        prev = 0
        for _ in range(0,obslog.names.count(name)):

            obstime, _, owner, prev = log_values(name, obslog, prev)

            if owner == 'public':
                owner = 'RECUR_A100'

            if isinstance(obstime,float):
                t = datetime.datetime.utcfromtimestamp(obstime)
            else:
                hr, mn = obstime
                if type(intime) != datetime.datetime:
                    intime = datetime.datetime.utcnow()
                t = datetime.datetime(intime.year, intime.month, intime.day, hr, mn)

            jd = round(float(ephem.julian_date(t)), 4)

            selection = (star_table['name'] == name) & (star_table['sheetn'] == owner)
            if any(selection):
                if np.any(jd > star_table['lastobs'][selection]):
                    star_table['lastobs'][selection] = jd
                    star_table['nobs'][selection] += 1
                    star_table['night_obs'][selection] += 1
                    log_str = "Updating local googledex star %s " % (name)
                    log_str += "in program %s to %.4f" % (owner, jd)
                    apflog(log_str, echo=True)
            elif too_table is not None:
                selection = (too_table['name'] == name) & (too_table['sheetn'] == owner)
                if any(selection) and jd > too_table['lastobs'][selection]:
                    log_str =  "Updating ToO target %s " % (name)
                    log_str += "from time %.4f to %.4f" % (too_table['lastobs'][selection], jd)
                    apflog(log_str, echo=True)
                    too_table['lastobs'][selection] = jd
                    too_table['nobs'][selection] += 1

    astropy.io.ascii.write(star_table,outfn, format='ecsv', overwrite=True)
    if too_table is not None:
        astropy.io.ascii.write(too_table,toofn, format='ecsv', overwrite=True)
        star_table = astropy.table.vstack(too_table,star_table)

    return obslog, star_table

def observed_JD(star_table_row,otime,ctime):
    '''
    observed_JD(star_table_row,otime,ctime)

    star_table_row - row from the googledex
    otime - observation time as a time stamp
    ctime - current time as a time stamp

    returns the julian date of the last observation
    '''
    jd = None

    # ideally get JD from the local table -
    # this will always be the value of the last observation
    if  star_table_row is not None:
        if len(star_table_row['lastobs']) > 0:
            jd = float(star_table_row['lastobs'][0])

    # if the above fails, we should be able to use the observing log
    # but this is JUST the UT hour and minute, not the day so we have to use the otime
    # value to calculate the full JD
    if jd is None:
        if isinstance(otime,float):
            t = datetime.datetime.utcfromtimestamp(otime)
        else:
            hr, mn = otime
            t = datetime.datetime(ctime.year, ctime.month, ctime.day, hr, mn)
            jd = float(ephem.julian_date(t))

    return jd

def log_values(local_name, obslog, prev):
    '''
    log_values(local_name, obslog, prev)

    local_name - name of the target in the local googledex
    obslog - ObservedLog object
    prev - previous index in the ObservedLog object

    returns the observation time, temperature, owner,
        and the index of the next observation matching local_name
    '''
    nameidx = obslog.names.index(local_name,prev)

    prev = nameidx+1
    # observation details
    otime = obslog.times[nameidx]
    taketemp_val = obslog.temps[nameidx]
    taketemp = False
    if taketemp_val == 'Y':
        taketemp = True
    curowner = obslog.owners[nameidx]
    if curowner == 'public':
        curowner = 'RECUR_A100'

    return otime, taketemp, curowner, prev

def update_a_sheet(sheetn, obslog, star_table, ctime):
    '''
    update_a_sheet(worksheet_vals, obslog)

    worksheet_vals - list of lists from the google sheet
    obslog - ObservedLog object

    returns the number of cells updated
    '''

    if ctime is None:
        ctime = datetime.datetime.utcfromtimestamp(int(time.time()))

    worksheet, worksheet_vals = get_spreadsheet(sheetn=sheetn,certificate=DEFAULT_CERT)
    if worksheet is None:
        return 0
    # Google does not like too many requests at once
    time.sleep(len(worksheet_vals))

    nupdates = 0
    n_temps = 0

    didx = find_columns(worksheet_vals[0], ['Star Name', 'lastobs', 'Nobs', 'Template', 'night_obs','pri','Total Obs'])

    for i, v in enumerate(worksheet_vals):
        # Starting at the top of vals is important.
        # the i is the row in the list of lists.
        # v is the current row in the list.
        # The columns that are updated are assigned above

        # Did we observe this target tonight?
        local_name = parse_starname(v[didx['Star Name']])

        if local_name in obslog.names:
            # We observed this target, so update the cell in the worksheet
            prev = 0
            for n_appear in range(0,obslog.names.count(local_name)):
                # a target can be observed more than once a night

                otime, taketemp, curowner, prev = log_values(local_name, obslog, prev)

                try:
                    c_row = star_table['name'] == local_name
                    c_row = c_row & (star_table['sheetn'] == sheetn)
                    star_table_row = star_table[c_row]
                except:
                    star_table_row = None

                jd = observed_JD(star_table_row, otime, ctime)
                # this will skip a done target
                # this is done in case a target was observed
                # and the target appears in the sheet more than once
                if v[didx['pri']] is not None:
                    pri = int_default(v[didx['pri']])
                    if pri < 0:
                        continue
                if v[didx['Total Obs']] is not None:
                    total_obs = int_default(v[didx['Total Obs']])
                    if total_obs > 0 and int_default(v[didx['Nobs']]) >= total_obs:
                        continue

                try:
                    pastdate = float(v[didx['lastobs']])
                except:
                    pastdate = 0.0

                try:
                    nobs = int(v[didx['Nobs']])
                except:
                    nobs = 0

                if didx['Template'] > 0:
                    try:
                        have_temp = v[didx['Template']]
                        if taketemp and have_temp == "N" and curowner == sheetn:
                            worksheet.update_cell(i+1, didx['Template']+1, "Y")
                            nupdates += 1
                            n_temps += 1
                            apflog( "Updated %s to having a template in %s" % (v[didx['Star Name']],sheetn),echo=True)
                    except:
                        apflog( "Error logging template obs for %s" % (v[didx['Star Name']]),echo=True,level='error')
                else:
                    pass

                new_nobs = nobs + n_appear + 1
                # update_cell(row, col, val) - col and row are 1 indexed
                try:
                    if round(jd, 3) > pastdate and curowner == sheetn and not taketemp:
                        worksheet.update_cell(i+1, didx['lastobs']+1, round(jd, 3) )
                        worksheet.update_cell(i+1, didx['Nobs']+1, new_nobs )
                        if didx['night_obs'] >= 0:
                            worksheet.update_cell(i+1, didx['night_obs']+1, n_appear+1)
                            nupdates += 1
                        log_str = "Updated %s from %.4f " % (v[didx['Star Name']],pastdate)
                        log_str += "to %.4f and %d in %s"  % (round(jd, 3),new_nobs,sheetn)
                        nupdates += 2
                        apflog(log_str, echo=True)
                except:
                    apflog("Error updating %s in %s" % (v[didx['Star Name']],sheetn),echo=True,level='error')
                    apflog("Cannot update %s to %.4f and %d in %s" % (v[didx['Star Name']],round(jd,3),1,sheetn),echo=True)

    return nupdates, n_temps


def update_online_sheets(observed_file,ctime=None,certificate=DEFAULT_CERT,outfn='parsesched.dat',outdir=None):
    """
        Update the online googledex  assuming things in filename have been observed.
        update_online_sheets(filename, sheetn="The Googledex",time=None,certificate=DEFAULT_CERT)

        filename - where the observations are logged
        ctime - current time as a time stamp
        certificate - required for authentication
        outfn - the local copy of the star list
        outdir - directory of data files
        returns the number of cells updated
    """

    if not outdir :
        outdir = os.getcwd()

    obslog = ObservedLog.ObservedLog(filename=os.path.join(outdir, observed_file))
    if len(obslog.names) == 0:
        return

    outfn = os.path.join(outdir,outfn)
    star_table = read_star_table(outfn)

    # OK the following code is going to rely on the fact that owner
    # is the sheet name
    # this will need to be updated when we start using coverid
    needed_sheetns = list( set(obslog.owners))

    if 'public' in needed_sheetns:
        needed_sheetns.remove('public')
        needed_sheetns.append('RECUR_A100')

    nupdates = 0
    n_temps = 0
    for sheetn in needed_sheetns:

        try:
            c_updates, c_temps = update_a_sheet(sheetn, obslog, star_table, ctime)
            nupdates += c_updates
            n_temps += c_temps
        except:
            apflog("Error updating sheet %s" % (sheetn),echo=True,level='error')
            continue

    return nupdates, n_temps
