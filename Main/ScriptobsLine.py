# ScriptobsLine.py
# pure string generation for scriptobs lines, split out of UCOScheduler.py
from __future__ import print_function

import numpy as np

import SchedulerConsts

def make_scriptobs_line(star_table_row, t, decker="W", I2="Y", owner='public', focval=0, coverid='', temp=False):
    """ given a name, a row in a star table and a do_flag, will generate
    a scriptobs line as a string:
    line = make_scriptobs_line(star_table_row, t, decker="W",I2="Y")

    star_table_row -contains all of the data needed for the line except
    t - a datetime object, this is used to fill in the uth and utm fields
    decker - one character field for the decker, defaults to "W"
    I2 - one character field for whether or not the Iodine cell is in, must be "Y" or "N"
    temp - a boolean for whether or not this is a template observation
    """

    # Add the RA as three elements, HR, MIN, SEC
    rastr = "%s %s %s " % (star_table_row['RA hr'],
                           star_table_row['RA min'],
                           star_table_row['RA sec'])

    # Add the DEC as three elements, DEG, MIN, SEC
    decstr = "%s %s %s " % (star_table_row['Dec deg'],
                            star_table_row['Dec min'],
                            star_table_row['Dec sec'])
    # Start with the target name
    ret = "%s %s %s 2000 " % (str(star_table_row['name']),rastr, decstr)

    # Proper motion RA and DEC
    ret += 'pmra=%.4f ' % (star_table_row['pmRA'])
    ret += 'pmdec=%.4f ' % (star_table_row['pmDEC'])
    # V Mag
    ret += 'vmag=%.2f ' % (star_table_row['Vmag'])

    # T Exp
    if temp:
        ret += 'texp=1200 '
    else:
        ret += 'texp=%d ' % int(star_table_row['texp'])

    # I2
    if temp:
        I2 = 'N'
    ret += 'I2=%s ' % (I2)
    # lamp
    ret += 'lamp=none '
    # start time
    ret += 'uth=%02d utm=%02d ' % (int(t.hour),int(t.minute))

    # Exp Count
    if star_table_row['expcount'] > SchedulerConsts.EXP_LIM:
        ret += 'expcount=%.3g ' % (SchedulerConsts.EXP_LIM)
    elif temp:
        ret += 'expcount=%.3g ' % (1e9)
    else:
        ret += 'expcount=%.3g ' % (star_table_row['expcount'])
    # Decker
    if temp and star_table_row['Vmag'] > 9:
        decker = 'W'
    elif temp and star_table_row['Vmag'] <= 9:
        decker = 'N'
    ret += 'decker=%s ' % (decker)
    # do flag
    if star_table_row['do']:
        ret += 'do=Y '
    else:
        ret += 'do= '
    # Count
    if temp:
        count = num_template_exp(star_table_row['Vmag'])
    else:
        count = int(star_table_row['nexp'])

    ret += 'count=%d ' % (count)

    ret += 'foc=%d ' % (int(focval))

    if owner != '':
        if owner == 'RECUR_A100':
            owner = 'public'
        ret += 'owner=%s ' % str(owner)

    if coverid != '':
        ret += 'coverid=%s ' % str(coverid)

    ret += 'binning=%s ' % str(star_table_row['binning'])

#    if star_table_row['mode'] != None:
#        if star_table_row['mode'] == BLANK:
#            ret += ' blank=Y'
#        elif star_table_row['mode'] == ACQUIRE:
#            ret += ' guide=Y'
#    else:
#        ret += ''

#    raoff  = star_table_row['raoff']
#    decoff = star_table_row['decoff']
#    if raoff == 'None':
#        raoff = ''
#    if decoff == 'None':
#        decoff = ''
#    if raoff is not '' and decoff is not '':
#        ret += ' raoff=' + str(raoff) + ' decoff=' + str(decoff)

    return str(ret)

def num_template_exp(vmag):
    '''
    num_template_exp(vmag)

    vmag - V magnitude of target
    count - number of exposures for a template observation

    '''
    count = 7

    if vmag > 10:
        count = 9

    elif vmag  < 8:
        count = 5

    return count

def config_defaults(owner):
    '''
    config_defaults(owner)
    owner - string, owner of the targets

    config - dictionary of default values for the config
    '''

    config = dict()
    config['I2'] = 'Y'
    config['decker'] = 'W'
    config['mode'] = ''
    config['obsblock'] = ''
    config['Bstar'] = 'N'
    config['owner'] = owner
    config['inst'] = 'levy'
    config['raoff'] = ''
    config['decoff'] = ''

    return config

# make_obs_block has no live caller: the obsblock path in
# UCOScheduler.make_result is commented out. It is kept here because
# re-enabling obsblocks is a plausible future want.
def make_obs_block(star_table, idx, dt, focval):
    '''

    make_obs_block(star_table, idx, dt, focval)

    star_table - astropy table of targets
    idx - index of target in star_table
    dt - datetime object
    focval - focus value

    rv - list of scriptobs lines for an obsblock
    '''

    rv = []

    cur_obsblock = star_table['obsblock'][idx]

    allinblock = star_table['obsblock'] == cur_obsblock
    allinblock = allinblock & (star_table['sheetn'] == star_table['sheetn'][idx])

    if np.any(star_table['mode'][allinblock] == SchedulerConsts.FIRST):
        first = star_table['mode'][allinblock] == SchedulerConsts.FIRST
    elif np.any(star_table['mode'][allinblock] == SchedulerConsts.ACQUIRE):
        first = star_table['mode'][allinblock] == SchedulerConsts.ACQUIRE
    else:
        first = None

    if np.any(star_table['mode'][allinblock] == SchedulerConsts.LAST):
        last = star_table['mode'][allinblock] == SchedulerConsts.LAST
    else:
        last = None

    rest = star_table['mode'][allinblock] != SchedulerConsts.FIRST
    rest = rest & (star_table['mode'][allinblock] != SchedulerConsts.ACQUIRE)
    rest = rest & (star_table['mode'][allinblock] != SchedulerConsts.LAST)
    rest_idxs, = np.where(rest)


    if np.any(first):
        first_idxs, = np.where(first)
        for idx in first_idxs:
            scriptobs_line = make_scriptobs_line(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                                owner=star_table['sheetn'][allinblock][idx], \
                                                I2=star_table['I2'][allinblock][idx], focval=focval)
            rv.append(scriptobs_line)

    for idx in rest_idxs:
        scriptobs_line = make_scriptobs_line(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                               owner=star_table['sheetn'][allinblock][idx], \
                                               I2=star_table['I2'][allinblock][idx], focval=focval)
        rv.append(scriptobs_line)

    if np.any(last):
        last_idxs, = np.where(last)
        for idx in last_idxs:
            scriptobs_line = make_scriptobs_line(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                                owner=star_table['sheetn'][allinblock][idx], \
                                                I2=star_table['I2'][allinblock][idx], focval=focval)
            rv.append(scriptobs_line)


    rv.reverse()
    rv[0] += ' # obsblock=%s end' % (cur_obsblock)
    return rv
