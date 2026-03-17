from __future__ import print_function
import sys
import argparse
from datetime import datetime
import time

import numpy as np

sys.path.append("../Main")
import ParseUCOSched


if __name__ == "__main__":

    dt = datetime.utcfromtimestamp(int(time.time()))

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g","--googledex",
        metavar="googledex.dat",
        type=str,
        default="googledex.dat",
        help="Path to the googledex.dat file"
    )
    parser.add_argument(
        "-r","--rank_table",
        metavar="rank_table.dat",
        type=str,
        default="2026A_ranks",
        help="Name of the rank table"
    )
    args = parser.parse_args()

    rank_table = ParseUCOSched.make_rank_table(sheet_table_name=args.rank_table)
    star_table, stars = ParseUCOSched.parse_UCOSched(rank_table, outfn=args.googledex)


    need = (star_table['Template'] == 'N') & (star_table['I2'] == 'Y')
    need = need & (star_table['Bstar'] == 'N')
    need = need & (star_table['pri'] > 0)
    names = star_table['name'][need]
    for i in range(len(names)):
        ostr = names[i].strip() + " "
        ostr += star_table['sheetn'][need][i].strip() + " "
        ostr += str(star_table['Vmag'][need][i]) + " "
        ra = star_table['RA hr'][need][i] + ":" + star_table['RA min'][need][i] + ":" + star_table['RA sec'][need][i]
        dec = star_table['Dec deg'][need][i] + ":" + star_table['Dec min'][need][i] + ":" + star_table['Dec sec'][need][i]
        ostr += ra + " " + dec
        print(ostr)
