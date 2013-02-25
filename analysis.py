#!/usr/bin/env python

"""analysis.py

grab data from different sources and perform selections on that data based on 
certain criteria.

"""
import os
import sys
import random
from math import cos, sin, sqrt, pi, exp
import numpy as np
from scipy import stats
from optparse import OptionParser
import ConfigParser
from StringIO import StringIO
import subprocess
import zipfile
import uuid
import itertools

DEG2RAD = pi / 180.0
ATOM_NUM = [
    "ZERO", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uuo"]
LOOKUPDIR = "/shared_scratch/pboyd/OUTCIF/FinalCif"
WORKDIR = ""
ORG_MAX = 10 
ORG_PAIR_MAX = ORG_MAX / 2
FNL_MAX = 20
FNL_PAIR_MAX = FNL_MAX / 5
# following to be decided at run-time
DATA_MAX = 0
MET_MAX = 0
TOP_MAX = 0

class CSV(dict):
    """
    Reads in a .csv file for data parsing.

    """
    def __init__(self, filename, _MOFNAME=True):
        self._columns = {"MOF":"MOFname", "uptake":"mmol/g",
                      "temperature":"T/K", "pressure":"p/bar"}
        self.filename = filename
        head_read = open(filename, "r")
        self.headings = head_read.readline().lstrip("#").split(",")
        head_read.close()
        if _MOFNAME:
            self._parse_by_mofname()
        else:
            self._parse_by_heading()

    def obtain_data(self, column, _TYPE="float", **kwargs):
        """return the value of the data in column, based on values of 
        other columns assigned in the kwargs.

        """
        # create a set of lists for the data we care about
        trunc = []
        trunc_keys = {} 
        # check to see if the columns are actually in the csv file
        for ind, key in enumerate([column] + kwargs.keys()):
            try:
                rightkey = self._columns[key]
            except KeyError:
                rightkey = key
            if rightkey not in self.headings:
                print("Could not find the column %s in the csv file %s "%
                        (rightkey, self.filename) + "returning...")
                return 0. if _TYPE is "float" else None
            else:
                trunc.append(self[rightkey])
                trunc_keys[ind] = key
                if key == column:
                    colind = ind
    
        for entry in itertools.izip_longest(*trunc):
            # tie an entry list index to column + kwargs keys
            kwargs_id =[i for i in range(len(entry)) if trunc_keys[i] in 
                    kwargs.keys()]
            if all([entry[i] == kwargs[trunc_keys[i]] for i in kwargs_id]):
                # grab the entry for the column
                col = entry[colind]
                return float(col) if _TYPE is "float" else col

        print("Didn't find the data point requested in the csv file %s"%
                self.filename)
        return 0. if _TYPE is "float" else None

    def _parse_by_heading(self):
        """The CSV dictionary will store data to heading keys"""
        filestream = open(self.filename, "r")
        # burn the header, as it's already read
        burn = filestream.readline()
        for line in filestream:
            line = line.lstrip("#").split(",")
            for ind, entry in enumerate(line):
                try:
                    entry = float(entry)
                except ValueError:
                    #probably a string
                    pass
                self.setdefault(self.headings[ind], []).append(entry)
        filestream.close()

    def _parse_by_mofname(self):
        """The CSV dictionary will have MOFnames as headings and contain
        sub-dictionaries for additional row data.

        """
        filestream = open(self.filename, "r")
        try:
            mofind = self.headings.index(self._columns["MOF"])
        except ValueError:
            print("ERROR: the csv file %s does not have %s as a column! "%
                    (self.filename, self._columns["MOF"]) + 
                    "EXITING ...")
            sys.exit(0)

        try:
            uptind = self.headings.index(self._columns["uptake"])
        except ValueError:
            print("WARNING: the csv file %s does not have %s as a column"%
                    (self.filename, self._columns["uptake"]) +
                    " the uptake will be reported as 0.0 mmol/g")
        burn = filestream.readline()
        for line in filestream:
            line = line.lstrip("#").split(",")
            mofname = line[mofind].strip()
            mofname = clean(mofname)
            try:
                uptake = line[uptind]
            except UnboundLocalError:
                uptake = 0.
            self.setdefault(mofname, {})["mmol/g"] = float(uptake)
        filestream.close()


class FunctionalGroups(dict):
    """
    Reads in a .sqlout file and returns a dictionary containing mofnames and
    their functionalizations.

    """

    def __init__(self, filename):
        """Read in all the mofs and store in the dictionary."""
        filestream = open(filename, 'r')
        for line in filestream:
            line = line.split("|")
            mof = "%s.sym.%s"%(line[0], line[1])
            # use a dictionary to sort the functionalizations
            groups = line[2]
            dic = {}
            if len(groups) > 0:
                [dic.setdefault(i.split("@")[0], []).append(i.split("@")[1])
                        for i in groups.split(".")]
                if len(dic.keys()) == 1:
                    dic[None] = []
            # check if the entry already exists!
            if self._check_duplicate(mof):
                if self[mof] == dic:
                    pass
                    # duplicate
                    #print "Duplicate found %s"%(mof)
                else:
                    print "Warning: duplicate found for %s"%(mof) + \
                            " but with different functionalization!"
            else:
                self[mof] = dic
        filestream.close()

    def _check_duplicate(self, entry):
        """Return true if the key exists, otherwise, false."""
        if self.has_key(entry):
            return True
        return False

class MOFlist(list):
    """Returns a list of MOFs parsed from a simple text file."""

    def __init__(self, filename):
        filestream = open(filename, "r")
        header = filestream.readline()
        mofind = header.split(',').index("MOFname")
        for line in filestream:
            line = line.strip()
            line = line.split(",")
            line = line[mofind]
            if "sym" in line:
                line.lstrip("#")
                line = clean(line)
                self.append(line)

        filestream.close()

class Selector(object):
    """
    Take in a dictionary of a mof, complete with data, and select a set 
    according to some criteria.

    """
    metal_indices = {
        "Zn" : (1, 3),
        "Cu" : (2,),
        "Co" : (4,),
        "Cd" : (5,),
        "Mn" : (6,),
        "Zr" : (7,),
        "In" : (8, 11),
        "V" : (9,),
        "Ba" : (11,),
        "Ni" : (12,)
        }

    bad_organics = [
            9,
            15,
            10,
            (2, 11),
            (20, 23),
            (24, 27),
            (15, 25)
            ]

    functional_groups = [
            "H",
            "F",
            "Cl",
            "Br",
            "I",
            "NH2",
            "Me",
            "Et",
            "Pr",
            "NO2",
            "HCO",
            "COOH",
            "OH",
            "CN",
            "Ph",
            "OMe",
            "OEt",
            "OPr",
            "NHMe",
            "SO3H"
            ]

    def __init__(self, mof_dic, metals=None, topologies=[],
                 ignore='used_mofs'):
        # THESE global variables are adjusted below.
        global TOP_MAX
        global MET_MAX
        self.mof_dic = mof_dic.copy()
        used_mofs = MOFlist(ignore)
        # remove mofs in the dictionary which have already been used
        for mof in used_mofs:
            try:
                self.mof_dic.pop(mof)
            except KeyError:
                pass
        # assert if a metal is selected
        self.metalind = []
        if metals is not None:
           for metal in metals:
               # try to extract a metal index from the list of metals
               try:
                   metind = int(metal)
                   if not metind in [ind for sublist in 
                                        self.metal_indices.values()
                                        for ind in sublist]:
                       print("ERROR: metal index %i is not in the database!"
                               %metind)
                       sys.exit(1)
                   self.metalind.append(metind)
               except ValueError:
                   try:
                       indices = self.metal_indices[metal]
                       for ind in indices:
                           self.metalind.append(ind)
                   except KeyError:
                       print("ERROR: metal %s is not in the database!"%(metal))
                       sys.exit(1)
        # all metals are fair game
        else:
           for metalind in self.metal_indices.values():
               for i in metalind:
                   self.metalind.append(i)
        # set upper bounds for the topologies and metal indices
        if topologies:
            TOP_MAX = DATA_MAX / len(topologies) + 1
        MET_MAX = DATA_MAX / len(self.metalind) + 1
        self.topologies = topologies
        top_bool = True if topologies else False
        self.trim_undesired(_TOPOLOGY=top_bool)

    def trim_undesired(self, _METAL=True, _ORGANIC=True, _TOPOLOGY=True):
        """Reduces the number of mofs in the mof_dic dictionary based
        on requested metal indices, topologies, and bad organic species

        """
        temp_moflist = self.mof_dic.keys()
        for mof in temp_moflist:
            met, org1, org2, top, fnl = parse_mof_data(mof)
            pop = False
            if _METAL:
                if met() not in self.metalind:
                    pop = True
            if _TOPOLOGY:
                if top() not in self.topologies:
                    pop = True
            if _ORGANIC:
                o1 = org1()
                o2 = org2()
                orgpair = tuple(sorted([o1, o2]))
                if (o1 in self.bad_organics):
                    pop = True
                elif (o2 in self.bad_organics):
                    pop = True
                elif (orgpair in self.bad_organics):
                    pop = True
            if pop:
                self.mof_dic.pop(mof)

    def trim_non_existing(self):
        """Reduces the number of mofs in the mof_dic based on if the actual
        .cif file for the MOF exist in a specified directory. If it does not,
        the mof is deleted from the dictionary.

        """
        dir = LOOKUPDIR
        print "Checking %s for all mofs..."%(dir)
        temp_moflist = self.mof_dic.keys()
        for mof in temp_moflist:
            # re-name mof to the out directory standards.
            met, org1, org2, top, fnl = parse_mof_data(mof)
            newmofname = "str_m" + str(met()) + "_o" + str(org1()) + \
                    "_o" + str(org2()) + "_f0_" + top() + ".sym." + \
                    str(fnl()) + ".out.cif"
            if not os.path.isfile(dir + "/" + newmofname):
                self.mof_dic.pop(mof)
                #print "%s  does not exist!"%(mof)

    def top_select(self, exclude=[], inclusive=[], partial=[], gridmax=None):
        """Order the mofs by top ranked structures.  store in a dataset
        dictionary and write to a csv file.

        """
        dataset = {}
        moflist = self.gen_moflist(inclusive, partial, exclude)
        if inclusive or partial:
            for group in inclusive + partial:
                # obtain list of mofs containing group,

                partial_list = self.isolate_group(moflist, group)
                ranked_list = self.rank_by_uptake(partial_list)
                for mof in ranked_list:
                    ngrid = self.grid_points(mof, gridmax)
                    ngrid_test = (ngrid > 0 and (ngrid <= gridmax if
                                  gridmax is not None else True))
                    if ngrid_test:
                        groups = self.mof_dic[mof]['functional_groups'].keys()
                        dataset.setdefault(tuple(groups), []).append(mof)
            for key, value in dataset.items():
                dataset[key] = value[:FNL_MAX]
        else:
            ranked_list = self.rank_by_uptake(moflist)
            for mof in ranked_list[:DATA_MAX]:
                groups = self.mof_dic[mof]['functional_groups'].keys()
                ngrid = self.grid_points(mof, gridmax)
                ngrid_test = (ngrid > 0 and (ngrid <= gridmax if
                              gridmax is not None else True))
                if ngrid_test:
                    dataset.setdefault(tuple(groups), []).append(mof)

        self.write_dataset(dataset, gridmax=None)

    def rank_by_uptake(self, mof_list):
        order = []
        for mof in mof_list:
            # grab uptake
            try:
                uptake = self.mof_dic[mof]['mmol/g']
            except KeyError:
                uptake = 0.

            order.append(tuple([uptake, mof]))
        order = [i for i in reversed(sorted(order))]
        order = [i[1] for i in order[:]]
        return order

    def isolate_group(self, moflist, group):
        plist = []
        for mof in moflist:
            groups = self.mof_dic[mof]['functional_groups'].keys()
            if group in groups:
                plist.append(mof)
        return plist

    def set_maxima(self, exclude, inclusive, partial, oexclude, oinclude,
                   topologies):
        """Set maximum values for the global variables based on the
        values set by the arrays.

        """
        global ORG_MAX
        global ORG_PAIR_MAX
        global FNL_MAX
        global FNL_PAIR_MAX
        global MET_MAX
        global TOP_MAX
        # set upper bounds for the topologies and metal indices.
        if topologies:
            TOP_MAX = DATA_MAX / len(topologies) + 1
        MET_MAX = DATA_MAX / len(self.metalind) + 1


    def random_select(self, exclude=[], inclusive=[], partial=[],
                      weight=None, gridmax=None):
        """Select a list of MOFs randomly 
        'exclude' defines functional groups which will not be represented

        'inclusive' defines functional groups which will be represented as the
            only functional group in the mof

        'partial' defines functional groups which can be paired with other
            functional groups excluding those in 'exclude'

        """
        print "Performing random selection.. values used:\n" + \
              "inclusive functional groups: ", inclusive, "\n" + \
              "exclude functional groups: ", exclude, "\n" + \
              "partial functional groups: ", partial, "\n" + \
              "weight uptake (gaussian): %s\n"%weight + \
              "max grid points: ",gridmax, "\n" + \
              "Total number of MOFs: %i\n"%DATA_MAX + \
              "Maximum per functional group: %i\n"%FNL_MAX + \
              "Maximum per functional group pair: %i\n"%FNL_PAIR_MAX + \
              "Maximum per organic linker: %i\n"%ORG_MAX + \
              "Maximum per organic linker pair: %i\n"%ORG_PAIR_MAX + \
              "Maximum per metal index: %i\n"%MET_MAX + \
              "Maximum per topology: %i\n"%TOP_MAX + \
              "Gridmax is set to ", gridmax
        dataset = {}
        organics_count, fnl_groups_count = {}, {}
        top_count, met_index_count = {}, {}
        mofcount = 0
        if inclusive and partial:
            # check if any co-exist in these lists
            pair = [set(i) for i in itertools.product(inclusive, partial)
                    if len(set(i)) == 1]
            if pair:
                for i in pair:
                    group = i[0]
                    print "Warning duplicate %s found in inclusive"%(group) + \
                          " partial lists. Appending to inclusive only..."
                    partial.pop(partial.index(group))
        # generate a list of valid mofs which obey the inclusive, partial and 
        # exclude lists.
        moflist = self.gen_moflist(inclusive, partial, exclude)
        done = False
        while not done:
            try:
                mof = random.choice(moflist)
                moflist.pop(moflist.index(mof))

            except IndexError:
                print "Sampled all MOFs without completing list! Writing " + \
                      "output file anyways.."
                self.write_dataset(dataset, gridmax)
                return
            met, org1, org2, top, fnl = parse_mof_data(mof)
            uptake = self.mof_dic[mof]['mmol/g']
            try:
                (fnl_grp1, fnl_grp2) = self.mof_dic[mof]['functional_groups']
            except ValueError:
                # doesn't find two functional groups in the dictionary
                fnl_grp1, fnl_grp2 = None, None
            except KeyError:
                # doesn't find the key [mof] or ['functional_groups'] in their
                # associated dictionaries
                fnl_grp1, fnl_grp2 = None, None
            org_max = self.check_dictionary_counts(organics_count,
                                                organic1=org1(),
                                                organic2=org2())
            fnl_max = self.check_dictionary_counts(fnl_groups_count, 
                                                   fnl_group1=fnl_grp1, 
                                                   fnl_group2=fnl_grp2)
            top_max = self.check_dictionary_counts(top_count,
                                                   topology=top())
            met_max = self.check_dictionary_counts(met_index_count,
                                                   metal=met())
            if weight is not None:
                upt_wght = self.weight_by_gaussian(uptake)
            else:
                upt_wght = True
            ngrid = self.grid_points(mof, gridmax)
            ngrid_test = (ngrid > 0 and (ngrid <= gridmax if 
                         gridmax is not None else True))
            # check to see if only the weighted uptake failed.
            # If so, put the mof back in the pool for selection later..
            # I did this because the lists are not completing.
            if not upt_wght and not org_max and not fnl_max and not top_max \
                    and not met_max and ngrid_test:
                # put the mof back in the random selection pool.
                moflist.append(mof)

            if ngrid_test and not org_max and not fnl_max and not top_max and \
                    not met_max and upt_wght:
                # increment counts 
                self.increment_dictionary_counts(organics_count,
                                                 org1(), 
                                                 org2())
                # keep track of functional group counts
                self.increment_dictionary_counts(fnl_groups_count,
                                                 fnl_grp1,
                                                 fnl_grp2)
                # keep track of metal index counts
                self.increment_dictionary_counts(met_index_count,
                                                 met())
                # keep track of topology counts
                self.increment_dictionary_counts(top_count,
                                                 top())
                # append to list
                dataset.setdefault(tuple((fnl_grp1, fnl_grp2)), []).append(mof)
                self.mof_dic[mof]['ngrid'] = ngrid
                mofcount += 1
            if mofcount >= DATA_MAX:
                done = True

        self.write_dataset(dataset, gridmax)

    def check_dictionary_counts(self, dictionary, **kwargs):

        if 'fnl_group1' in kwargs.keys():
            maximum = FNL_MAX
            combine_max = FNL_PAIR_MAX
        elif 'organic1' in kwargs.keys():
            maximum = ORG_MAX
            combine_max = ORG_PAIR_MAX
        elif 'metal' in kwargs.keys():
            maximum = MET_MAX
            combine_max = MET_MAX
        elif 'topology' in kwargs.keys():
            maximum = TOP_MAX
            combine_max = TOP_MAX
        # check the combination of arguments
        # only currently valid for organic and functional group dictionaries.
        combine = sorted(kwargs.keys())
        try:
            combine.pop(combine.index(None))
        except ValueError:
            pass
        if len(combine) > 1:
            dictionary.setdefault(tuple(combine), 0)
            if dictionary[tuple(combine)] >= combine_max:
                return True 
        # check the individual value counts
        for key, value in kwargs.items():
            dictionary.setdefault(value, 0)
            if dictionary[value] >= maximum:
                return True 
        return False

    def increment_dictionary_counts(self, dictionary, *args):
        # increment individual entries.
        for arg in args:
            if arg is not None:
                dictionary.setdefault(arg, 0)
                dictionary[arg] += 1
        # increment combinations as well...
        if len(set(args)) > 1 and None not in args:
            dictionary.setdefault(tuple(args), 0)
            dictionary[tuple(args)] += 1

    def weight_by_gaussian(self, uptake, a=1, b=5, c=1.3):
        """Weight according to a gaussian distribution based on uptake.
        Defaults are a distribution amplitude of 1, centered around
        5 mmol/g, smeared by 1.5 (width)"""
        gauss = gaussian(a, b, c)
        # probability of accepting the gaussian weight determined
        # by a random number between 0 and a.
        if gauss(uptake) < (random.random()*a):
            return False
        return True

    def grid_points(self, mofname, gridmax):
        """Determine the maximum number of grid points needed for the esp."""
        # have to source the correct file.
        mofname = clean(mofname)
        dirmof = (LOOKUPDIR + '/' + 
                  mofname + '.out.cif')
        ngrid = -1
        if os.path.isfile(dirmof) and gridmax is not None:
            from_cif = CifFile(dirmof)
            ngrid = GrabGridPoints(from_cif.cell,
                                   from_cif.atom_symbols,
                                   from_cif.cart_coordinates)
        elif gridmax is None:
            ngrid = 1

        return ngrid 

    def write_dataset(self, dataset, gridmax):
        """Writes the data to a file."""
        os.chdir(WORKDIR)
        basename = ""
        if self.metalind:
            metal_titles = []
            for met in self.metalind:
                for key, value in self.metal_indices.items():
                    if met in value:
                        metal_titles.append(key)
            for metal in set(metal_titles):
                basename += "%s_"%metal
        basename += "dataset"
        if self.topologies:
            for top in self.topologies:
                basename += "_%s"%top
        count = 0
        filename = create_csv_filename(basename) 
        print("Writing dataset to %s.csv..."%(filename))
        outstream = open(filename+".csv", "w")
        header="MOFname,mmol/g,functional_group1,functional_group2"
        if gridmax is not None:
            header += ",ngrid\n"
        else:
            header += "\n"
        outstream.writelines(header)
        for fnl, mofs in dataset.iteritems():
            if len(fnl) == 1:
                print "Functional group: %s"%fnl
                groups = tuple([fnl[0], None])
            elif len(fnl) == 2:
                print "Functional groups: %s, %s"%(fnl)
                groups = fnl
            for mof in mofs:
                print "     " + mof
                try:
                    uptake = self.mof_dic[mof]['mmol/g']
                except KeyError:
                    uptake = 0.
                line = "%s,%f,%s,%s"%(mof, uptake, groups[0], groups[1])
                if gridmax is not None:
                    try:
                        ngrid = self.mof_dic[mof]['ngrid']
                    except KeyError:
                        ngrid = 0
                    line += ",%i\n"%(ngrid)
                else:
                    line += "\n"
                outstream.writelines(line)
        outstream.close()
        print("Done.")

    def gen_moflist(self, inclusive, partial, exclude):
        """Returns a list of MOFs which are chosen based on the 
        discrimination input lists.

        """
        moflist = []
        if not inclusive and not partial and not exclude:
            return self.mof_dic.keys()

        for key, value in self.mof_dic.iteritems():
            try:
                (group1, group2) = value['functional_groups'].keys()
            except KeyError:
                value['functional_groups'] = {"None1":[], "None2":[]}
                (group1, group2) = (None, None)
            except ValueError:
                group1, group2 = None, None

            if group1 in exclude or group2 in exclude:
                # do not append
                pass
            elif (group1 in partial or group2 in partial) and \
                (group1 not in inclusive and group2 not in inclusive) \
                and partial:
                # append
                moflist.append(key)
            elif ((group1 in inclusive and group2 is None) or
                  (group1 is None and group2 in inclusive)) and inclusive:
                # append
                moflist.append(key)
            elif not inclusive and not partial:
                # append if no restrictions made on inclusive and partial
                moflist.append(key)

        return moflist

class GrabNewData(object):
    """Takes a mof dictionary and adds new uptake data to it."""
    def __init__(self, mofs, basedir, extended=False):
        self.mofs = mofs
        # determines if hydrogen replacements are shown in the output
        self.extended = extended
        # base directory containing all the faps job info.
        self.basedir = basedir
    def grab_data(self, temp=298.0, press=0.15):
        """Descends into directories and grabs appropriate data."""
        os.chdir(self.basedir)
        all_list = os.listdir('.')
        directories = [i for i in all_list if os.path.isdir(i)]
        for mof in self.mofs.keys():
            if mof in directories:
                os.chdir(mof)
                if os.path.isfile(mof + "-CO2.csv"):
                    data = CSV(mof + "-CO2.csv", _MOFNAME=False)
                    new_uptake = data.obtain_data("mmol/g",
                                                  temperature=temp,
                                                  pressure=press)
                else:
                    print "ERROR: could not find %s-CO2.csv"%(mof)
                    new_uptake = 0.
                os.chdir('..')
            else:
                print "ERROR: could not find %s in the directory %s"%(
                        mof, self.basedir)
                new_uptake = 0. 
            self.mofs[mof]['new_uptake'] = new_uptake
        os.chdir(WORKDIR)

    def write_data(self, filename="default.report.csv"):
        """Write all MOF data to a csv file."""
        os.chdir(WORKDIR)
        basename = clean(filename)
        count = 0
        # make sure there's no overwriting, goes up to 99
        filename = create_csv_filename(basename)
        print("Writing report to %s.csv..."%(filename))
        outstream = open(filename + ".csv", "w")
        if self.extended:
            header = "MOFname,mmol/g,Functional_grp1," +\
                    "Functional_grp2,grp1_replacements,grp2_replacements\n"
        else:
            header = "MOFname,mmol/g,Functional_grp1," +\
                    "Functional_grp2\n"

        outstream.writelines(header)
        for mof in self.mofs.keys():
            csv_uptake = self.mofs[mof]['mmol/g']
            new_uptake = self.mofs[mof]['new_uptake']
            try:
                replaced_groups = self.mofs[mof]['functional_groups']
            except KeyError:
                replaced_groups = {"None1": [], "None2": []}
            names = replaced_groups.keys()
            names.sort()
            if len(names) == 2:
                fnl_grp1 = names[0]
                fnl_grp2 = names[1]
                rep_1 = " ".join(replaced_groups[fnl_grp1])
                rep_2 = " ".join(replaced_groups[fnl_grp2])

            elif len(names) == 1:
                fnl_grp1 = names[0]
                fnl_grp2 = None
                rep_1 = " ".join(replaced_groups[fnl_grp1])
                rep_2 = ["N/A"]

            if self.extended:
                fmt = "%s,%f,%s,%s,%s,%s\n"
                line = fmt%(mof, new_uptake,
                            fnl_grp1, fnl_grp2, rep_1, rep_2)
            else:
                fmt = "%s,%f,%s,%s\n"
                line = fmt%(mof, new_uptake, fnl_grp1,
                            fnl_grp2)
            outstream.writelines(line)
        print("Done.")
        outstream.close()

class Options(object):
    """Read in the options from the config file."""
    def __init__(self, name=None):
        self.job = ConfigParser.SafeConfigParser()
        self.defaults = ConfigParser.SafeConfigParser()
        self._set_paths()
        
    def _set_paths(self):
        if __name__ != '__main__':
            self.script_dir = os.path.dirname(__file__)
        else:
            self.script_dir = os.path.abspath(sys.path[0])
        self.job_dir = os.getcwd()

    def load_defaults(self):
        default_path = os.path.join(self.script_dir, 'defaults.ini')
        print default_path
        try:
            filetemp = open(default_path, 'r')
            default = filetemp.read()
            filetemp.close()
            if not '[defaults]' in default.lower():
                default = '[defaults]\n' + default
            default = StringIO(default)
        except IOError:
            print("Error loading defaults.ini")
            default = StringIO('[defaults]\n')
        self.defaults.readfp(default)

class CommandLine(object):
    """Parse command line options and communicate directives to the program."""

    def __init__(self):
        self.commands = {}
        self.command_options()

    def command_options(self):
        usage = "%prog [options]\n" + \
                "%prog -r -d Cu_dataset/NO_CHG -q cusql.sqlout" +\
                " -c Cu_dataset.csv\n" + \
                "%prog -s -M Cu -i F,Cl,Br,I,SO3H,NHMe -q cusql.sqlout" +\
                " -G 150000 -I used_mofs -c combined.csv " + \
                "-L /scratch/tdaff/FinalCif\n" + \
                "%prog -t -M Cu,Zn -x F,Cl,Br,I,SO3H,NHMe -q allsql.sqlout" +\
                " -G 150000 -I used_mofs -c combined.csv -N 300 -F 20\n" + \
                "%prog -e cif_whole_error -c combined.csv -q allsql.sqlout"
        parser = OptionParser(usage=usage)
        parser.add_option("-r", "--report", action="store_true",
                          dest="report",
                          help="issue a report based on some calculated data.")
        parser.add_option("-s", "--dataset", action="store_true",
                          dest="dataset",
                          help="generate a dataset of randomly selected MOFs.")
        parser.add_option("-e", "--extract", action="store", type="string",
                          dest="extract",
                          help="write a report on a list of MOFs which " + \
                               "contains counting of functional groups etc.")
        parser.add_option("-t", "--topranked", action="store_true",
                          dest="top",
                          help="create a dataset of top ranked structures.")
        parser.add_option("-c", "--csvfile", action="store", type="string",
                          dest="csvfilename",
                          help="location of csv file with existing data in it.")
        parser.add_option("-d", "--directory", action="store",
                          dest="dirname",
                          help="location of directory containing new data.")
        parser.add_option("-q", "--sqlfile", action="store", type="string",
                          dest="sqlname",
                          help="location of sql file containing mof "+
                               "functionalizations.")
        parser.add_option("-I", "--ignore", action="store", type="string",
                           dest="ignorefile", default='used_mofs',
                           help="location of list of MOFs to ignore in sampling.")
        parser.add_option("-x", "--excludelist", 
                          type="string", dest="exclude",
                          action="callback",
                          callback=self.parse_commas, default=[],
                          help="comma (,) delimited list of functional groups"+
                               " to exclude from the sampling")
        parser.add_option("-i", "--inclusivelist", 
                          type="string", dest="inclusive", 
                          action="callback",
                          callback=self.parse_commas, default=[],
                          help="comma (,) delimited list of functional groups"+
                               " to sample, excluding all others.")
        parser.add_option("-p", "--partiallist", 
                          type="string", dest="partial",
                          action="callback", default=[],
                          callback=self.parse_commas,
                          help="comma (,) delimited list of functional groups"+
                               " to include in the sampling.")
        parser.add_option("-L", "--lookupdir", action="store", type="string",
                          dest="lookup", 
                          default="/shared_scratch/pboyd/OUTCIF/FinalCif",
                          help="lookup directory with all the output cifs.")
        parser.add_option("-M", "--metal", dest="metals", type="string",
                          action="callback", 
                          callback=self.parse_commas, default=[],
                          help="specify metals to generate dataset. Current" +
                          " options are: Zn, Cu, Co, Cd, Mn, Zr, In, V, Ba or Ni" + 
                          " or you can enter the associated indices") 
        parser.add_option("-T", "--topologies", dest="topologies", type="string",
                          action="callback", 
                          callback=self.parse_commas, default=[],
                          help="specify topologies to include in the dataset." +
                          " delimited by commas") 
        parser.add_option("-G", "--gridpoints", action="store", type="int",
                          dest="maxgridpts",
                          help="specify the max number of grid points to "+\
                               "allow the structures to be chosen")
        parser.add_option("-N", "--mofcount", action="store", type="int",
                          dest="nummofs", default="120",
                          help="Number of MOFs to include in the set.")
        parser.add_option("-F", "--fnlmax", action="store", type="int",
                          dest="fnlmax", default="20",
                          help="Maximum number of MOFs with a particular "+\
                                "functional group.")
        parser.add_option("-C", "--combine", type="string",
                          dest="combine", action="callback",
                          callback=self.parse_commas,
                          help="comma (,) delimited list of csv files to " + \
                          "merge into a larger file.")
        parser.add_option("-W", "--weight",
                          action="store_true",
                          help="Make the random selected data a gaussian " + \
                            "weight on the CO2 uptake.")

        (local_options, local_args) = parser.parse_args()
        self.options = local_options

    def parse_commas(self, option, opt, value, parser):
        setattr(parser.values, option.dest, value.split(','))

class GrabGridPoints(int):
    """Extracts an integer of grid points from an egulp calculation."""
    code_loc = "/home/pboyd/bin"
    # directory to submit jobs and remove later.
    dir = "/home/pboyd/.temp"
    config_file = 'tempconfig.ini'
    param_file = 'tempparam.dat'
    config_lines = ("build_grid  1\n" +
                   "build_grid_from_scratch 1 none 0.25 0.25 0.25 2.0\n" +
                   "save_grid 0 temp\n" +
                   "calculate_pot_diff 0\n" +
                   "skip_everything 1\n" +
                   "point_charges_present 0\n" +
                   "include_pceq 0\n" +
                   "imethod 0\n")
    param_lines = ("1\n    1   4.52800000   6.94520000")
 
    def __new__(cls, cell, atoms, coordinates, value=0):
        i = int.__new__(cls, value)
        i.cell = cell
        # make the geo_file difficult to duplicate
        i.geo_file = str(uuid.uuid1())
        i.atoms = atoms
        i.nums = [ATOM_NUM.index(j) for j in atoms]
        i.coordinates = coordinates
        i._write_geo_file()
        i = i._submit_egulp_job()
        return i

    def _write_geo_file(self):
        """Sets up the .geo file to run e-gulp."""
        filestream = open(self.dir + "/" + self.geo_file, "w")
        filestream.writelines("temp\n")
        a, b, c = self.cell[:]
        filestream.writelines(('%15.12f %15.12f %15.12f \n'%(tuple(a))))
        filestream.writelines(('%15.12f %15.12f %15.12f \n'%(tuple(b))))
        filestream.writelines(('%15.12f %15.12f %15.12f \n'%(tuple(c))))
        filestream.writelines("%i\n"%(len(self.atoms)))
        for atom, (x, y, z) in zip(self.nums, self.coordinates):
            filestream.writelines("%6d %12.7f %12.7f  %12.7f %12.7f\n"%
                                    (atom, x, y, z, 0.))
        filestream.close()

    def _submit_egulp_job(self):
        """Submits the EGULP calculation to determine the number of grid
        points to calculate.

        """
        os.chdir(self.dir)
        if not os.path.isfile(self.config_file):
            configstream = open(self.config_file, 'w')
            configstream.writelines(self.config_lines)
            configstream.close()
        if not os.path.isfile(self.param_file):
            parameterstream = open(self.param_file, 'w')
            parameterstream.writelines(self.param_lines)
            parameterstream.close()
        code_exe = self.code_loc + "/egulppot"
        job = subprocess.Popen([code_exe, self.geo_file, 
                                self.param_file, self.config_file],
                                shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        comm = job.communicate()
        if comm[1] is not '':
            print "ERROR in EGULP calculation!"
            return -1
        else:
            os.remove(self.geo_file)
            return self._grab_ngrid_points(comm[0].split('\n')) 

    def _grab_ngrid_points(self, lines):
        """Grab the number of grid points from the output file."""
        search_string = "total number of simulation points"
        for line in reversed(lines):
            if search_string in line:
                parse = line.lstrip(search_string)
                return int(parse.strip())
        return -1
 
class CifFile(object):
    """This class will grab data from the cif file such as the
    cartesian atom coordinates and the cell vectors.
    
    """

    def __init__(self, moffilename):
        """Read in the mof file and parse the data to write a .geo file
        for the egulp method.

        """
        self.mofdata = parse_cif(moffilename)
        self.cell = self._grab_cell_vectors()
        self.atom_symbols = self._grab_atom_names()
        self.cart_coordinates = self._grab_cartesians()

    def _grab_cell_vectors(self):
        """Select the lattice parameters from the mofdata and convert to a
        3x3 array of cell vectors.

        """
        a = float(self.mofdata['_cell_length_a'])
        b = float(self.mofdata['_cell_length_b'])
        c = float(self.mofdata['_cell_length_c'])
        alpha = DEG2RAD * float(self.mofdata['_cell_angle_alpha'])
        beta = DEG2RAD * float(self.mofdata['_cell_angle_beta'])
        gamma = DEG2RAD * float(self.mofdata['_cell_angle_gamma'])
        veca = a * np.array([1., 0., 0.])
        vecb = np.array([b * cos(gamma), b * sin(gamma), 0.])
        c_x = c * cos(beta)
        c_y = c * (cos(alpha) - cos(gamma)*cos(beta)) / sin(gamma) 
        vecc = np.array([c_x, c_y,
                         sqrt(c**2 - c_x**2 - c_y**2)])
        return np.array([veca, vecb, vecc])

    def _grab_atom_names(self):
        """Retun a tuple of the atom names."""
        return tuple(self.mofdata['_atom_site_type_symbol'])

    def _grab_cartesians(self):
        """Return a tuple of coordinates of all the atoms."""
        coordinates = []
        for x, y, z in zip(self.mofdata['_atom_site_fract_x'],\
                       self.mofdata['_atom_site_fract_y'],\
                       self.mofdata['_atom_site_fract_z']):
            scaled_pos = np.array([float(x), float(y), float(z)])
            coordinates.append(tuple(np.dot(scaled_pos, self.cell)))
        return tuple(coordinates)

def parse_cif(moffilename):
    """Return a dictionary containing parsed information."""
    cifstream = open(moffilename, "r")
    cifdata = {}
    loopindices = {}
    loopcount = 0
    loopread = False
    for line in cifstream:
        line = line.strip()
        if not loopread and line.startswith("_"):
            parsedline = line.split()
            cifdata[parsedline[0]] = parsedline[1]
        if loopread and line.startswith("_"):
            loopindices[loopcount] = line
            cifdata[line] = []
            loopcount += 1
        # append data to the appropriate heading
        elif loopread:
            parsedline = line.split()
            # pretty crappy check to see if we're on a line of data
            if len(parsedline) == loopcount:
                for ind, entry in enumerate(parsedline):
                    cifdataentry = loopindices[ind]
                    cifdata[cifdataentry].append(entry)
            else:
                loopread = False

        elif "loop_" in line:
            loopindices = {}
            loopcount = 0
            loopread = True
    cifstream.close()
    return cifdata 


def parse_mof_data(mofname):
    mofname = clean(mofname)
    mofparse = mofname.split("_")
    data = mofparse[5].split(".")
    def metal_index():
        return int(mofparse[1].lstrip("m"))

    def org1_index():
        return int(mofparse[2].lstrip("o"))

    def org2_index():
        return int(mofparse[3].lstrip("o"))

    def topology():
        return data[0]

    def fnl_number():
        try:
            num = int(data[2])
        except IndexError:
            num = 0
        return num 

    return metal_index, org1_index, org2_index, topology, fnl_number


def pair_csv_fnl(mof, fnl):
    """Join csv and fnl to one dictionary."""
    for name in mof.keys():
        try:
            mof[name]["functional_groups"] = fnl[name]

        except KeyError:
            mof[name]["functional_groups"] = {"None1":[], "None2":[]}
            #print "Warning: could not find data for %s"%name

    del fnl

def create_a_dataset(csvfile=None, sqlfile=None, metals=[],
                     inclusive=[], exclude=[], partial=[],
                     gridmax=None,
                     ignore=None,
                     topologies=[],
                     weight=None):
    """Create a dataset with pre-defined number of MOFs. Change in code
    if needed.

    """
    mofs = CSV(csvfile)
    fnl = FunctionalGroups(sqlfile)
    pair_csv_fnl(mofs, fnl)
    sel = Selector(mofs, metals=metals, ignore=ignore, topologies=topologies)
    sel.trim_non_existing()
    sel.random_select(exclude=exclude, inclusive=inclusive, partial=partial, 
                      gridmax=gridmax, weight=weight)

def write_report(directory=None, sqlfile=None, csvfile=None):
    """Write a report on some new uptake data located in 'directory' based
    on some old data in the 'csvfile'

    """
    mofs = CSV(csvfile)
    fnl = FunctionalGroups(sqlfile)
    pair_csv_fnl(mofs, fnl)
    data = GrabNewData(mofs, basedir=directory, extended=False)
    data.grab_data()
    base_list = directory.split("/")
    base_list = [i for i in base_list if i]
    # take the last two entries in the directory as the file name
    basename = ".".join(base_list[-2:])
    dir = WORKDIR 
    data.write_data(filename=dir + "/" + basename + ".report.csv")

def generate_top_structures(csvfile=None, sqlfile=None, 
                            exclude=[], inclusive=[], partial=[],
                            metals=[],
                            gridmax=None,
                            ignore=None):
    
    mofs = CSV(csvfile)
    fnl = FunctionalGroups(sqlfile)
    pair_csv_fnl(mofs, fnl)
    sel = Selector(mofs, metals=metals, ignore=ignore)
    sel.trim_non_existing()
    sel.top_select(exclude=exclude, inclusive=inclusive, 
                   partial=partial, gridmax=gridmax)

def combine_csvs(*args):
    """Combine csv files to a single, unified csv.  Compute correlation
    coefficients for the uptake columns, if the data match.

    """
    basenames = [clean(i.split("/")[-1]) for i in args]
    csvs = [CSV(i) for i in args]
    merged = {}
    # re-name the uptake column, this will be the only columns which
    # are not over-written
    for ind, csv in enumerate(csvs):
        for key, dic in csv.items():
            uptake = dic.pop('mmol/g')
            dic['mmol/g-%s'%basenames[ind]] = uptake
            merged.setdefault(key, {})
            merged[key].update(dic)
    # combine uptake columns
    corr = [[] for i in range(len(basenames))]
    for mof, value in merged.items():
        if all(['mmol/g-%s'%i in value.keys() for i in basenames]):
            if all([value['mmol/g-%s'%i] > 0. for i in basenames]):
                [corr[basenames.index(i)].append(value['mmol/g-%s'%i]) for 
                        i in basenames]
            else:
                print("MOF %s not included in correlation calculations "%mof
                        + "because of 0 value for uptake")
        else:
            print("MOF %s is not included in correlation calculations "%mof +
                    "because it couldn't be found in one or more csv files.")
    for upt1, upt2 in itertools.combinations(corr, 2):
        ind1 = corr.index(upt1)
        ind2 = corr.index(upt2)
        print("Correlation coefficients for %s, and %s:"%(
                 basenames[ind1], basenames[ind2]))
        rho, pval = stats.spearmanr(upt1, upt2)
        print("Spearman rank: %7.5f"%rho)
        pears, p2 = stats.pearsonr(upt1, upt2)
        print("Pearson correlation: %7.5f"%pears)

    names = {}
    [names.setdefault(i.split('.')[0], 0) for i in basenames]
    name = names.keys()[0]
    name += ".merged"
    write_csv(name, merged)

def write_csv(basename, dic):
    """Writes a specific csv where the dictionary contains a mof name
    followed by a dictionary of values.

    """
    filename = create_csv_filename(basename)
    print("Writing to %s.csv"%filename)
    outstream = open(filename + ".csv", "w")
    lead = []
    headings = {}
    # set up the headings first, then append data
    for mof, data in dic.items():
        [headings.setdefault(i, []) for i in data.keys()]
    # now append data
    for mof, data in dic.items():
        lead.append(mof)
        for head in headings.keys():
            try:
                headings[head].append(data[head])
            except KeyError:
                headings[head].append("")
    outstream.writelines("MOFname," + ",".join(headings.keys()) + "\n")
    for ind, entries in enumerate(itertools.izip(*headings.values())):
        line = "%s,"%lead[ind]
        line += ",".join([str(i) for i in entries])
        line += "\n"
        outstream.writelines(line)
    print("Done.")
    outstream.close()

def extract_info(file_name=None, sqlfile=None, csvfile=None):
    """Extract as much information from a list of mofnames, and provide
    an ammended report with the findings.

    """
    all_mofs = CSV(csvfile)
    all_fnls = FunctionalGroups(sqlfile)
    pair_csv_fnl(all_mofs, all_fnls)
    mofs = MOFlist(file_name)
    mof_dat = {}
    fnl_count = {}
    fnl_pair = {}
    organic_count = {}
    metal_count = {}
    org_pair_count = {}
    full_mof_count = {}
    for mof in mofs:
        mof = clean(mof)
        try:
            data = all_mofs[mof]
            mof_dat[mof] = data
        except KeyError:
            print("ERROR: could not find data for %s"%mof)

    def analyse_data():
        for mof in mof_dat.keys():
            met, org1, org2, top, junk = parse_mof_data(mof)
            stuff = mof_dat[mof]
            fnls = stuff['functional_groups'].keys()
            if len(fnls) == 1:
                fnl1 = fnls[0]
                fnl2 = None
                fnl_count.setdefault(fnl1, 0) 
                fnl_count[fnl1] += 1
            elif len(fnls) == 2:
                fnl1, fnl2 = fnls
                fnl_count.setdefault(fnl1, 0)
                fnl_count[fnl1] += 1
                fnl_count.setdefault(fnl2, 0)
                fnl_count[fnl2] += 1
                pair = sorted([fnl1, fnl2])
                fnl_pair.setdefault(tuple(pair), 0) 
                fnl_pair[tuple(pair)] += 1
            metal_count.setdefault(met(), 0)
            metal_count[met()] += 1
            if org1() == org2():
                organic_count.setdefault(org1(), 0) 
                organic_count[org1()] += 1
            else:
                organic_count.setdefault(org1(), 0)
                organic_count[org1()] += 1
                organic_count.setdefault(org2(), 0) 
                organic_count[org2()] += 1
                pair = sorted([org1(), org2()])
                org_pair_count.setdefault(tuple(pair), 0)
                org_pair_count[tuple(pair)] += 1
            fullmof = "str_m%i_o%i_o%i_%s"%(met(),org1(),org2(),top())
            full_mof_count.setdefault(fullmof, 0) 
            full_mof_count[fullmof] += 1

    def write_overall_csv():
        os.chdir(WORKDIR)
        basename = clean(filename)
        count = 0
        files = []
        print("Writing overall reports...")
        # first file prints out individual functional group stats
        fnlfile = basename + ".fnl_grps"
        # create original filename
        fnlfile = create_csv_filename(fnlfile)
        print("Writing %s.csv"%fnlfile)
        # store filename for archiving later
        files.append("%s.csv"%fnlfile)
        outstream = open(fnlfile + ".csv", "w")
        # write data to file
        outstream.writelines("#functional_group,count\n")
        for fnl, count in fnl_count.items():
            outstream.writelines("%s,%i\n"%(fnl,fnl_count[fnl]))
        outstream.close()
        # second file prints out pair functional group stats
        fnlpfile = basename + ".fnl_pairs"
        fnlpfile = create_csv_filename(fnlpfile)
        print("Writing %s.csv"%fnlpfile)
        files.append("%s.csv"%fnlpfile)
        outstream = open(fnlpfile + ".csv", "w")
        outstream.writelines("#functional_group1,functional_group2,count\n")
        for pair, count in fnl_pair.items():
            outstream.writelines("%s,%s,%i\n"%(pair[0], pair[1], count))
        outstream.close()
        # third file prints out individual organic group stats
        orgfile = basename + ".org"
        orgfile = create_csv_filename(orgfile)
        print("Writing %s.csv"%orgfile)
        files.append("%s.csv"%orgfile)
        outstream = open(orgfile + ".csv", "w")
        outstream.writelines("#organic_index,count\n")
        for org, count in organic_count.items():
            outstream.writelines("%i,%i\n"%(org, count))
        outstream.close()
        # fourth file prints out organic pair stats
        orgpfile = basename + ".org_pairs"
        orgpfile = create_csv_filename(orgpfile)
        print("Writing %s.csv"%orgpfile)
        files.append("%s.csv"%orgpfile)
        outstream = open(orgpfile + ".csv", "w")
        outstream.writelines("#organic_index1,organic_index2,count\n")
        for orgs, count in org_pair_count.items():
            outstream.writelines("%i,%i,%i\n"%(orgs[0], orgs[1], count))
        outstream.close()
        # fifth file prints out metal index stats
        metfile = basename + ".metal"
        metfile = create_csv_filename(metfile)
        print("Writing %s.csv"%metfile)
        files.append("%s.csv"%metfile)
        outstream = open(metfile + ".csv", "w")
        outstream.writelines("#metal_index,count\n")
        for met, count in metal_count.items():
            outstream.writelines("%i,%i\n"%(met, count))
        outstream.close()
        # sixth file prints out whole mof stats
        moffile = basename + ".mofname"
        moffile = create_csv_filename(moffile)
        print("Writing %s.csv"%moffile)
        files.append("%s.csv"%moffile)
        outstream = open(moffile + ".csv", "w")
        outstream.writelines("#MOFname,count\n")
        for mof, count in full_mof_count.items():
            outstream.writelines("%s,%i\n"%(mof, count))
        outstream.close()
        print("Done.")

        # Archive data
        zipname = basename + ".overall_reports"
        zipname = create_csv_filename(zipname, extension=".zip")
        print("Archiving to %s.zip ..."%zipname)
        zip = zipfile.ZipFile(zipname + ".zip", "w")
        for f in files:
            zip.write("%s"%f)
        zip.close()
        print("Done.")
        # Clear filenames
        for f in files:
            os.remove(f)

    def write_specific_csv():
        os.chdir(WORKDIR)
        basename = clean(file_name)
        count = 0
        filename = basename + ".specific_report"
        filename = create_csv_filename(filename)
        outstream = open(filename + ".csv", "w")
        print ("Writing specific report to %s.csv"%filename)
        header = "MOFname,functional_group1,functional_group2," +\
                "mof_occur,metal_occur,org1_occur,org2_occur," +\
                "org_pair_occur,fnl1_occur,fnl2_occur," +\
                "fnl_pair_occur\n"
        outstream.writelines(header)
        for mof, vals in mof_dat.items():
            met, org1, org2, top, fnl = parse_mof_data(mof)
            fnls = vals['functional_groups'].keys()
            # grab all the data for each
            fnls = sorted(fnls)
            if len(fnls) == 1:
               fnl1 = fnls[0]
               fnl2 = None
               fnl1_occur = fnl_count[fnl1]
               fnl_pair_occur = fnl_count[fnl1]
               fnl2_occur = 0
            elif len(fnls) == 2:
               fnl1, fnl2 = fnls
               fnl1_occur = fnl_count[fnl1]
               fnl2_occur = fnl_count[fnl2]
               pair = sorted([fnl1, fnl2])
               try:
                   fnl_pair_occur = fnl_pair[tuple(pair)]
               except KeyError:
                   fnl_pair_occur = 0
            met_occur = metal_count[met()]
            orgs = sorted([org1(), org2()])
            if len(set(orgs)) == 1:
                org1_occur = organic_count[org1()]
                org2_occur = 0
                org_pair_occur = organic_count[org1()]
            elif len(set(orgs)) == 2:
                org1_occur = organic_count[org1()]
                org2_occur = organic_count[org2()]
                org_pair_occur = org_pair_count[tuple(orgs)]

            mof_abbrev = "str_m%i_o%i_o%i_%s"%(met(),org1(),org2(),top())
            mof_occur = full_mof_count[mof_abbrev]
            outstream.writelines("%s,%s,%s,%i,%i,%i,%i,%i,%i,%i,%i\n"
                    %(mof,fnl1,fnl2,mof_occur,met_occur,org1_occur,
                        org2_occur,org_pair_occur,fnl1_occur,fnl2_occur,
                        fnl_pair_occur))
        outstream.close()
        print("Done.")
    return analyse_data, write_overall_csv, write_specific_csv

def create_csv_filename(basename, extension=".csv"):
    count = 0
    filename = basename
    # make sure there's no overwriting, goes up to 99
    while os.path.isfile(filename + extension):
        count += 1
        filename = basename + ".%02d"%(count)
    return filename

def clean(name):
    if name.endswith('.cif'):
        name = name[:-4]
    elif name.endswith('.niss'):
        name = name[:-5]
    elif name.endswith('-CO2.csv'):
        name = name[:-8]
    elif name.endswith('.flog'):
        name = name[:-5]
    elif name.endswith('.out.cif'):
        name = name[:-8]
    elif name.endswith('.tar'):
        name = name[:-4]
    elif name.endswith('.db'):
        name = name[:-3]
    elif name.endswith('.faplog'):
        name = name[:-7]
    elif name.endswith('.db.bak'):
        name = name[:-7]
    elif name.endswith('.csv'):
        name = name[:-4]
    return name

def gaussian(a, b, c):
    def g(x):
        return a*exp(-(x - b)**2/(2*c**2))
    return g

def main():
    global LOOKUPDIR
    global WORKDIR
    global DATA_MAX
    global FNL_MAX
    options = Options()
    options.load_defaults()
    cmd = CommandLine()
    LOOKUPDIR = cmd.options.lookup
    WORKDIR = os.getcwd()
    DATA_MAX = cmd.options.nummofs
    FNL_MAX = cmd.options.fnlmax
    if cmd.options.dataset:
        if cmd.options.csvfilename:
            test = os.path.isfile(cmd.options.csvfilename)
            if not test:
                print("ERROR: could not find the .csv file")
                sys.exit()
        else:
            print("ERROR: .csv filename not set in the command line")
            sys.exit()
        if cmd.options.sqlname:
            test = os.path.isfile(cmd.options.sqlname)
            if not test:
                print("ERROR: could not find the sql file")
                sys.exit()
        else:
            print("ERROR: sql file not set in the command line")
        if cmd.options.maxgridpts:
            print("Excluding MOFs with gridpoints exceeding %i"%
                    (cmd.options.maxgridpts))

        create_a_dataset(csvfile=cmd.options.csvfilename,
                         sqlfile=cmd.options.sqlname,
                         metals=cmd.options.metals,
                         inclusive=cmd.options.inclusive,
                         exclude=cmd.options.exclude,
                         partial=cmd.options.partial,
                         gridmax=cmd.options.maxgridpts,
                         ignore=cmd.options.ignorefile,
                         topologies=cmd.options.topologies,
                         weight=cmd.options.weight)
                         
    if cmd.options.report:
        if cmd.options.dirname:
            test = os.path.isdir(cmd.options.dirname)
            if not test:
                print("ERROR: could not find the directory containing MOFs")
        else:
            print("ERROR: directory not set in the command line")
            sys.exit()
        if cmd.options.csvfilename:
            test = os.path.isfile(cmd.options.csvfilename)
            if not test:
                print("ERROR: could not find the .csv file")
                sys.exit()
        else:
            print("ERROR: .csv filename not set in the command line")
            sys.exit()
        if cmd.options.sqlname:
            test = os.path.isfile(cmd.options.sqlname)
            if not test:
                print("ERROR: could not find the sql file")
                sys.exit()
        else:
            print("ERROR: sql file not set in the command line")
            sys.exit()

        write_report(directory=cmd.options.dirname,
                     sqlfile=cmd.options.sqlname,
                     csvfile=cmd.options.csvfilename)

    if cmd.options.extract:
        if cmd.options.csvfilename:
            test = os.path.isfile(cmd.options.csvfilename)
            if not test:
                print("ERROR: could not find the .csv file")
                sys.exit()
        else:
            print("ERROR: .csv filename not set in the command line")
            sys.exit()
        if cmd.options.sqlname:
            test = os.path.isfile(cmd.options.sqlname)
            if not test:
                print("ERROR: could not find the sql file")
                sys.exit()
        else:
            print("ERROR: sql file not set in the command line")
        if not os.path.isfile(cmd.options.extract):
            print("ERROR: could not find the mof file")
            sys.exit()

        an, over, spec = extract_info(
                     file_name=cmd.options.extract,
                     sqlfile=cmd.options.sqlname,
                     csvfile=cmd.options.csvfilename)
        an()
        over()
        spec()

    if cmd.options.top:
        if cmd.options.csvfilename:
            test = os.path.isfile(cmd.options.csvfilename)
            if not test:
                print("ERROR: could not find the .csv file")
                sys.exit()
        else:
            print("ERROR: .csv filename not set in the command line")
            sys.exit()
        if cmd.options.sqlname:
            test = os.path.isfile(cmd.options.sqlname)
            if not test:
                print("ERROR: could not find the sql file")
                sys.exit()

        generate_top_structures(csvfile=cmd.options.csvfilename,
                                sqlfile=cmd.options.sqlname,
                                exclude=cmd.options.exclude,
                                inclusive=cmd.options.inclusive,
                                partial=cmd.options.partial,
                                metals=cmd.options.metals,
                                gridmax=cmd.options.maxgridpts,
                                ignore=cmd.options.ignorefile
                                )
    if cmd.options.combine:
        combine_csvs(*cmd.options.combine)

if __name__ == "__main__":
    main()
