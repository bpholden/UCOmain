# Target.py
# the target returned by UCOScheduler.get_next
from __future__ import print_function

# old result-dict key -> attribute, in the order the dict was built
_KEYS = [
    ('RA', 'ra'),
    ('DEC', 'dec'),
    ('PM_RA', 'pm_ra'),
    ('PM_DEC', 'pm_dec'),
    ('VMAG', 'vmag'),
    ('BV', 'bv'),
    ('COUNTS', 'counts'),
    ('EXP_TIME', 'exp_time'),
    ('NEXP', 'nexp'),
    ('TOTEXP_TIME', 'totexp_time'),
    ('NAME', 'name'),
    ('PRI', 'pri'),
    ('DECKER', 'decker'),
    ('I2', 'i2'),
    ('BINNING', 'binning'),
    ('isTemp', 'is_temp'),
    ('isBstar', 'is_bstar'),
    ('isTOO', 'is_too'),
    ('mode', 'mode'),
    ('owner', 'owner'),
    ('SCRIPTOBS', 'scriptobs'),
    ('template_conditions_met', 'template_conditions_met'),
]


class Target(object):
    '''
    A target selected by UCOScheduler.get_next.

    scriptobs - the scriptobs lines for this target, in pop order
    (the last line in the list is sent first).
    Values taken from the star table keep their numpy types.

    '''
    def __init__(self, name, owner, ra, dec, pm_ra, pm_dec, vmag, bv, counts,
                 exp_time, nexp, totexp_time, pri, decker, i2, binning,
                 is_bstar=False, is_too=False, mode=''):
        self.ra = ra
        self.dec = dec
        self.pm_ra = pm_ra
        self.pm_dec = pm_dec
        self.vmag = vmag
        self.bv = bv
        self.counts = counts
        self.exp_time = exp_time
        self.nexp = nexp
        self.totexp_time = totexp_time
        self.name = name
        self.pri = pri
        self.decker = decker
        self.i2 = i2
        self.binning = binning
        self.is_temp = False
        self.is_bstar = is_bstar
        self.is_too = is_too
        self.mode = mode
        self.owner = owner
        self.scriptobs = []
        self.template_conditions_met = False

    @classmethod
    def from_star_table(cls, star_table, idx, star, totexptime, pri, bstar=False):
        '''
        Build a Target from row idx of star_table.
        star is the computed ephem body for that row.

        '''
        return cls(name=star_table['name'][idx],
                   owner=star_table['sheetn'][idx],
                   ra=star.a_ra,
                   dec=star.a_dec,
                   pm_ra=star_table['pmRA'][idx],
                   pm_dec=star_table['pmDEC'][idx],
                   vmag=star_table['Vmag'][idx],
                   bv=star_table['B-V'][idx],
                   counts=star_table['expcount'][idx],
                   exp_time=star_table['texp'][idx],
                   nexp=star_table['nexp'][idx],
                   totexp_time=totexptime,
                   pri=pri,
                   decker=star_table['decker'][idx],
                   i2=star_table['I2'][idx],
                   binning=star_table['binning'][idx],
                   is_bstar=bstar,
                   is_too=star_table['too'][idx])

    def set_template_lines(self, lines, decker):
        '''
        Replace the scriptobs lines with a template sequence.

        '''
        self.scriptobs = list(lines)
        self.is_temp = True
        self.decker = decker

    def to_dict(self):
        '''
        The old result dict, same keys, values and order.

        '''
        return dict((key, getattr(self, attr)) for key, attr in _KEYS)

    def __repr__(self):
        return repr(self.to_dict())
