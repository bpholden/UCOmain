# test_UCOScheduler.py
# smoke tests for UCOScheduler, run as a script from a directory holding
# rank_table, googledex.dat and time_left.csv
from __future__ import print_function

import time
import datetime

import numpy as np

import ParseUCOSched
import Observability
import ScriptobsLine
import UCOTargets
from UCOScheduler import get_next

try:
    import ktl
except:
    pass

def test_basic_ops(ucotargets):
    """
    test_basic_ops()
    """

    # Test the basic operations of the scheduler
    # This is a test function to see if the basic operations work
    # It will not be run in production

    try:
        ktl.write('apftask', 'SCRIPTOBS_LINE_RESULT', 3, binary=True)
    except:
        pass

    # For some test input what would the best target be?
    OTFN = "observed_targets"
    ot = open(OTFN, "w")
    starttime = time.time()
    result = get_next(starttime, 7.99, 0.4, ucotargets, bstar=True, \
                      do_templates=False)
    while len(result['SCRIPTOBS']) > 0:
        ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
    ot.close()

    for i in range(5):

        result = get_next(starttime, 7.99, 0.4, ucotargets, bstar=False, \
                         do_templates=False)
        #result = smartList("tst_targets", time.time(), 13.5, 2.4)

        if result is None:
            print("Get None target")

        while len(result["SCRIPTOBS"]) > 0:
            ot = open(OTFN, "a+")
            while len(result['SCRIPTOBS']) > 0:
                ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
            ot.close()
            starttime += result["TOTEXP_TIME"]

    print("Done")
    ot.close()

    return starttime

def test_failure(starttime, ucotargets):
    '''
    test_failure(starttime, ucotargets)
    starttime - time to start the test
    ucotargets - UCOTargets object
    '''
    print("Testing a failure")
    try:
        ktl.write('apftask', 'SCRIPTOBS_LINE_RESULT', 2, binary=True)
    except:
        pass
    result = get_next(starttime, 7.99, 0.4, ucotargets, bstar=False, \
                     do_templates=True, )
    print(result)
    print("Nonsensical start time")
    result = get_next(starttime, 7.99, 0.4, ucotargets, bstar=True, \
                     do_templates=True, start_time=1)
    print(result)
    return

def test_templates(ucotargets):
    """
    test_templates(ucotargets)
    ucotargets - UCOTargets object
    """
    print("Testing templates")
    t_dt = datetime.datetime.now()
    tstar_table, _ = ParseUCOSched.parse_UCOSched(ucotargets.rank_table, \
                                                     outfn='googledex.dat', outdir=".", \
                                                        config=ScriptobsLine.config_defaults('public'))
    tidx, = np.asarray(tstar_table['name'] == '185144').nonzero()
    tidx = tidx[0]
    tbstars = (tstar_table['Bstar'] == 'Y')|(tstar_table['Bstar'] == 'y')
    tbidx, tbfinidx = Observability.find_Bstars(tstar_table, tidx, tbstars)
    decker = "N"
    tline  = ScriptobsLine.make_scriptobs_line(tstar_table[tidx], t_dt, \
                                decker=decker, I2="N", owner='public', temp=True)
    if "decker=W" in tline:
        decker = "W"
    tbline = ScriptobsLine.make_scriptobs_line(tstar_table[tbstars][tbidx], t_dt, \
                                decker=decker, I2="Y", owner='public', focval=2)

    tbfinline = ScriptobsLine.make_scriptobs_line(tstar_table[tbstars][tbfinidx], t_dt, \
                                   decker=decker, I2="Y", owner='public', focval=0)
    temp_res= []
    temp_res.append(tbfinline + " # temp=Y end")
    temp_res.append(tline + " # temp=Y")
    temp_res.append(tbline + " # temp=Y")
    out_r = [print(r) for r in temp_res]
    print(out_r)
    print("Done")

def test_main():
    """
    test_main()
    """

    # Test the basic operations of the scheduler
    # This is a test function to see if the basic operations work
    # It will not be run in production

    RANK_TABLEN='2025B_ranks_operational'

    class Opt:
        def __init__(self):
            self.test = True
            self.time_left = "time_left.csv"
            self.rank_table = RANK_TABLEN

    uco_targets = UCOTargets.UCOTargets(Opt())

    # this calls make_rank_table
    uco_targets.make_hour_constraints()
    # this calls make_hour_table 
    uco_targets.make_hour_table()

    starttime = test_basic_ops(uco_targets)

    test_failure(starttime, uco_targets)
    test_templates(uco_targets)

if __name__ == '__main__':

    test_main()
