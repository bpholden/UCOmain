#!/usr/bin/env python

with open('stars_APF','r') as f:
    new_stars = f.readlines()

with open('/u/user/devel_scripts/ucscapf/stars_APF','r') as f:
    old_stars = f.readlines()

new_names = [v.split()[0] for v in new_stars]

with open("Stars_Test",'w') as f:
    for l in old_stars:
        if l.split()[0] in new_names:
            f.write(l)

 
