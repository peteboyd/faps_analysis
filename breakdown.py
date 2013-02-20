#!/usr/bin/env python
import sys
import scipy
from scipy import stats
from itertools import izip
import math

if len(sys.argv) == 1:
    print "Usage: %prog [csv file]"
    sys.exit(1)
file = sys.argv[1]

def clean(name):
    if name.startswith("mmol/g-"):
        name = name[7:]
    if name.endswith(".report"):
        name = name[:-7]
    if name.endswith(".csv"):
        name = name[:-4]
    if name.endswith(".merged"):
        name = name[:-7]
    return name

def order(x):
    x = [float(i) for i in x]
    L = len(x)
    rangeL = range(L)
    z = izip(x, rangeL)
    z = izip(z, rangeL)
    D = sorted(z)
    D.reverse()
    return [d[1] for d in D]

def rank(x):
    x = [float(i) for i in x]
    L = len(x)
    ordering = order(x)
    ranks = [0] * len(x)
    for i in range(L):
        ranks[ordering[i]] = i
    ranks.reverse()
    return ranks

def grab_rep_key(dic):
    for i in dic.keys():
        if "REPEAT" in i:
            return i
    return None

def grab_percentage(L, percent):
    # assumes that the percentage will be greater than 1%
    if percent > 1.:
        percent = float(percent)/100.
    return int(math.floor(L*percent))

def compute_add_on(L, segment):
    excess = L % segment
    if excess != 0:
        print("the additional %i will be distributed evenly over all sets"%
            (L%segment))
        # TODO: include extra to the final set.
        distrib = excess / L 
        print("percentage of MOFs in each set: %7.4f"%(
            float(segment + distrib) / float(L)))
    else:
        distrib = 0
    return distrib

def strip_key(key):
    name = clean(key)
    name = [i for i in name.split(base) if i]
    if not name:
        name = "EGULP.param.1"
    else:
        name = name[0]
    if name.startswith("."):
        name = name[1:]
    if name.endswith("."):
        name = name[:-1]
    return name

base = clean(file)
filestream = open(file, "r")
headers = filestream.readline().lstrip("#").strip().split(",")
csv = {head: [] for head in headers}
to_csv = {i: head for i, head in enumerate(headers)}
for line in filestream:
    line = line.rstrip("#").strip().split(",")
    [csv[to_csv[i]].append(val) for i, val in enumerate(line)]
filestream.close()

rep_key = grab_rep_key(csv)
# re-order the columns
if rep_key is None:
    print("ERROR: could not find the REPEAT column "
           +"to properly rank the structures")
    sys.exit(1)
rankings = rank(csv[rep_key])
ordering = order(csv[rep_key])

# take the ten percents
segment = grab_percentage(len(rankings), 20)
L = len(rankings)
distrib = compute_add_on(L, segment)
# calculate the pearson and spearmans with repeat only..
start = 0
print("Ranking of MOFs based on REPEAT uptake...")
print("MOFname, REPEAt uptake, original order")
for ord in ordering:
    print "%s,%f,%i"%(csv['MOFname'][ord], float(csv[rep_key][ord]), ord)
for rankx in range(L / segment):
    # note, debug to make sure all this indexing is correct
    finish = start + segment + distrib
    # check to see if finish is close to the end.
    if (segment+distrib) > len(ordering[finish:]):
        finish += len(ordering[finish:])

    #order = ordering[start:finish-1]
    order = ordering[0: finish]
    #percentile = (float(start)/float(L)*100, float(finish-1)/float(L)*100)
    percentile = (float(1)/float(L)*100, float(finish-1)/float(L)*100)
    keys = [i for i in csv.keys() if i is not rep_key]
    keys.pop(keys.index("MOFname"))
    for key in keys:
        name = strip_key(key)
        #print("Comparing REPEAT with %s"%(name))
        outcsv = open("%s_correlations.csv"%(name), "a")
        if rankx == 0:
            outcsv.writelines("%s,%s,%s,%s,%s\n"%("order", "pearson", "spearman", "rms", "rmsd"))
        rep_vals = [float(csv[rep_key][i]) for i in order]
        other_vals = [float(csv[key][i]) for i in order]
        # remove zero entries
        for (rep, other) in zip(reversed(rep_vals), reversed(other_vals)):
            if rep == 0. or other == 0.:
                rep_vals.pop(rep_vals.index(rep))
                other_vals.pop(other_vals.index(other))
        print("%s, %s"%(rep_key, key))
        diff = []
        for rep_val, other_val in zip(rep_vals, other_vals):
            diff.append(other_val - rep_val)
            print rep_val, other_val
        rms = scipy.sqrt(scipy.mean([i*i for i in diff]))
        rmsd = scipy.std(diff)
        rho, pval = stats.spearmanr(other_vals, rep_vals)
        pears, pval2 = stats.pearsonr(other_vals, rep_vals)
        outcsv.writelines("%i,%f,%f,%f,%f\n"%(rankx, pears, rho, rms, rmsd))
        outcsv.close()
        #print("Pearson: %7.5f"%pears)
        #print("Spearman: %7.5f"%rho)
    start = finish

