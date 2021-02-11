#!/usr/bin/python
import csv
import pickle

fn = "googledex.csv"
googleReader = csv.reader(open('googledex.csv', 'rb'))

gd = []

for r in googleReader:
    gd.append(r)

gdat =open("googledex.dat","w")
pickle.dump(gd,gdat)

