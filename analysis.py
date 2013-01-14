#!/usr/bin/env python

"""analysis.py

grab data from different sources and perform selections on that data based on 
certain criteria.

"""
import os
import sys
import random
import math
# include readcsv stuff from read_csv.py
# read in csv file
# read in sql file
# match mofs with their functionalization
# remove mofs without actual structures.
# perform random selection stuff
# create zip file containing directory structure and .cif files
# create a report on the selection.

# include a dictionary containing mofnames as the key and uptakes as the value.

class CSV(dict):
    """
    Reads in a .csv file for data parsing.

    """
    def __init__(self, filename):
        self._mof_name_column = "MOFname"
        self._uptake_column = "mmol/g"
        self._temp_column = "T/K"
        self._press_column = "p/bar"
        self.filename = filename

    def read_from_csv_multiple(self):
        """Reads in a file, stores data in memory."""
        filestream = open(self.filename, "r")
        headings = filestream.readline().strip("#").split(",")
        mofind = headings.index(self._mof_name_column)
        uptind = headings.index(self._uptake_column)
        # read in the first line of the file as the headings
        # currently assumes only one uptake value at a statepoint
        for line in filestream:
            parsed = line.split(",")
            mofname = parsed[mofind].rstrip("-CO2.csv")
            uptake = parsed[uptind]
            self.setdefault(mofname, {})["old_uptake"] = float(uptake)

        filestream.close()

    def read_from_csv_single(self, T=None, P=None):
        filestream = open(self.filename, "r")
        # pull out the uptake in mmol/g from a given T and P
        headings = filestream.readline().strip("#").split(",")
        pressind = headings.index(self._press_column)
        tempind = headings.index(self._temp_column)
        uptkind = headings.index(self._uptake_column)

        for line in filestream:
            parsed = line.split(",")
            if T and P:
                if (T == float(parsed[tempind])) and\
                        (P == float(parsed[pressind])):
                    uptake = float(parsed[uptkind])
                    self['new_uptake'] = uptake
        try:
            self['new_uptake']
        except KeyError:
            self['new_uptake'] = None
            mofname = self.filename.rstrip("-CO2.csv")
            print "ERROR: Could not find the uptake values for %s"%(mofname)
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
            # check if the entry already exists!
            if self._check_duplicate(mof):
                if self[mof] == dic:
                    # duplicate
                    print "Duplicate found %s"%(mof)
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
        for line in filename:
            if "sym" in line:
                line.lstrip("#")
                line.rstrip("-CO2.csv")
                self.append(line)

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
        "Ba" : (10,),
        "Ni" : (12,)
        }

    bad_organics = [
            29,
            28,
            9,
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

    def __init__(self, mof_dic, metal=None, weight=False):
        self.mof_dic = mof_dic
        self.weight = weight
        used_mofs = MOFlist("used_mofs")
        # remove mofs in the dictionary which have already been used
        for mof in used_mofs:
            try:
                self.mof_dic.pop(mof)
            except KeyError:
                pass
        # assert if a metal is selected
        if metal:
           try:
               self.metalind = self.metal_indices[metal]
               self.trim_metals()

           except KeyError:
               self.metalind = None

    def trim_metals(self):
        """Reduces the number of mofs in the mof_dic dictionary based
        on requested metal indices.

        """
        temp_moflist = self.mof_dic.keys()
        # delete all metals not with the correct index
        for mof in temp_moflist:
            # obtain the metal index
            met, org1, org2, top, fnl = parse_mof_data(mof)
            ind = met()

            if ind not in self.metalind:
                self.mof_dic.pop(mof)

    def trim_organics(self):
        """Reduces the number of mofs in the mof_dic dictionary based
        on bad organic building units.  These building units, or 
        combinations of building units tend to break VASP.

        """

        temp_moflist = self.mof_dic.keys()
        for mof in temp_moflist:
            # obtain organic indices
            met, org1, org2, top, fnl = parse_mof_data(mof)
            o1 = org1()
            o2 = org2()
            orgpair = tuple(sorted([o1, o2]))
            if (o1 in self.bad_organics) or (o2 in self.bad_organics) or\
                    (orgpair in self.bad_organics):
                self.mof_dic.pop(mof)

    def trim_non_existing(self, dir="/shared_scratch/pboyd/OUTCIF/FinalCif"):
        """Reduces the number of mofs in the mof_dic based on if the actual
        .cif file for the MOF exist in a specified directory. If it does not,
        the mof is deleted from the dictionary.

        """
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
                print "%s  does not exist!"%(mof)

    def random_select(self, exclude=[], inclusive=[], partial=[]):
        """Select a list of MOFs randomly based on low uptake.
        'exclude' defines functional groups which will not be represented

        'inclusive' defines functional groups which will be represented as the
            only functional group in the mof

        'partial' defines functional groups which can be paired with other
            functional groups excluding those in 'exclude'

        """
        chosen = []
        dataset = {}
        inclusive_max = 100
        partial_max = 100
        group_max = 20
        if inclusive and partial:
            # check if any co-exist in these lists
            pair = [set(i) for i in itertools.product(inclusive, partial)
                    if len(set(i)) == 1]
            if pair:
                for i in pair:
                    group = i[0]
                    print "Warning duplicate %s found in inclusive"%(group) + \
                          "partial lists. Appending to inclusive only..."
                    partial.pop(partial.index(group))

        if inclusive:
            # To reduce the time spent on randomly selecting MOFs, organize 
            # into categories based on the inclusive functional groups
            for group in inclusive:
                moflist = self.gen_moflist(group, exclude)
                # do a random selection, determine uptake, check against
                # already selected mofs, etc.
                done = False
                group_count = 0
                while not done:
                    mof = random.choice(moflist)
                    skip = False
                    if mof in chosen:
                        skip = True
                    try:
                        uptake = self.mof_dic[mof]['old_uptake']
                    except KeyError:
                        uptake = 0.0
                    if self.weight:
                        # perform a weighted check against the uptake.
                        if math.exp(-uptake/4.) < random.random():
                            skip = True

                    if group_count >= group_max:
                        skip = True
                        done = True

                    if not skip:
                        group_count += 1
                        chosen.append(mof)
                        dataset.setdefault(tuple((group,)), []).append(mof)

        elif exclude and not inclusive:
            fnl_dic = {}
            mofcount = 0
            inclusive = self.functional_groups[:]
            for group in exclude:
                if group in inclusive:
                    id = inclusive.index(group)
                    inclusive.pop(id)

            [fnl_dic.setdefault(i, 0) for i in inclusive]
            moflist = self.mof_dic.keys() 
            done = False
            while not done:
                mof = random.choice(moflist)
                skip = False

                try:
                    groups = self.mof_dic[mof]['functional_groups'].keys()
                except KeyError:
                    skip = True

                if mof in chosen:
                    skip = True
                try:
                    uptake = self.mof_dic[mof]['old_uptake']
                except KeyError:
                    uptake = 0.0
                # perform a weighted check against the uptake.
                if self.weight:
                    if math.exp(-uptake/4.) < random.random():
                        skip = True
                for i in groups:
                    if i in exclude:
                        skip = True
                    elif fnl_dic[i] > group_max:
                        skip = True

                if not skip:
                    for i in groups:
                        fnl_dic[i] += 1
                    chosen.append(mof)
                    mofcount += 1
                    dataset.setdefault(tuple(groups), []).append(mof)
                if mofcount >= inclusive_max:
                    done = True

        self.write_dataset(dataset)
        # TODO(pboyd): add partial considerations.

    def write_dataset(self, dataset):
        """Writes the data to a file."""

        filename = ""
        if self.metalind:
            for key, value in self.metal_indices.items():
                if value == self.metalind:
                    filename += "%s_"%key
        filename += "dataset"
        count = 0
        while os.path.isfile(filename+".csv"):
            count += 1
            filename += "%02d"%(count)
        outstream = open(filename+".csv", "w")
        header="MOFname,mmol/g,functional_group1,functional_group2\n"
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
                uptake = self.mof_dic[mof]['old_uptake']
                outstream.writelines("%s,%f,%s,%s\n"%(mof, uptake, 
                                                      groups[0], groups[1]))
        outstream.close()

    def gen_moflist(self, group, exclude, incl=True):
        """Returns a list of MOFs based on a functional group."""
        moflist = []
        for key, value in self.mof_dic.iteritems():
            try:
                groups = set(value['functional_groups'].keys())
            except KeyError:
                value['functional_groups'] = {('None'):[]}
                groups = ('None')
            if incl:
                if len(groups) == 1 and list(groups)[0] == group:
                    moflist.append(key)
            else:
                append = True
                for grp in list(groups):
                    if grp in exclude:
                        append = False
                if append:
                    moflist.append(key)

        return moflist

class GrabNewData(object):
    """Takes a mof dictionary and adds new uptake data to it."""
    def __init__(self, mofs, basedir, extended=False):
        self.pwd = os.getcwd()
        self.mofs = mofs
        # determines if hydrogen replacements are shown in the output
        self.extended = extended
        # base directory containing all the faps job info.
        self.basedir = basedir

    def grab_data(self):
        """Descends into directories and grabs appropriate data."""
        os.chdir(self.basedir)
        all_list = os.listdir('.')
        directories = [i for i in all_list if os.path.isdir(i)]
        for mof in self.mofs.keys():
            if mof in directories:
                os.chdir(mof)
                if os.path.isfile(mof + "-CO2.csv"):
                    data = CSV(mof + "-CO2.csv")
                    data.read_from_csv_single(T=298.0, P=0.15)
                    new_uptake = data['new_uptake']
                else:
                    print "ERROR: could not find %s-CO2.csv"%(mof)
                    new_uptake = 0.
                os.chdir('..')
            else:
                print "ERROR: could not find %s in the directory %s"%(
                        mof, self.basedir)
                new_uptake = 0. 

            self.mofs[mof]['new_uptake'] = new_uptake
        os.chdir(self.pwd)

    def write_data(self, filename="default.report.csv"):
        """Write all MOF data to a csv file."""
        os.chdir(self.pwd)
        basename = filename.rstrip('.csv')
        done = False
        count = 0
        # make sure there's no overwriting, goes up to 99
        while os.path.isfile(filename):
            count += 1
            filename = basename + "%02d.csv"%(count)
        outstream = open(filename, "w")
        if self.extended:
            header = "MOFname,old_uptake,new_uptake,Functional_grp1" +\
                    "Functional_grp2,grp1_replacements,grp2_replacements\n"
        else:
            header = "MOFname,old_uptake,new_uptake,Functional_grp1" +\
                    "Functional_grp2\n"

        outstream.writelines(header)
        for mof in self.mofs.keys():
            old_uptake = self.mofs[mof]['old_uptake']
            new_uptake = self.mofs[mof]['new_uptake']
            try:
                replaced_groups = self.mofs[mof]['functional_groups']
            except KeyError:
                replaced_groups = {None: ["N/A"], None: ["N/A"]}
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
                fmt = "%s,%f,%f,%s,%s,%s,%s\n"
                line = fmt%(mof, old_uptake, new_uptake,
                            fnl_grp1, fnl_grp2, rep_1, rep_2)
            else:
                fmt = "%s,%f,%f,%s,%s\n"
                line = fmt%(mof, old_uptake, new_uptake, fnl_grp1,
                            fnl_grp2)
            outstream.writelines(line)
        outstream.close()

def parse_mof_data(mofname):
    mofname = mofname.rstrip("-CO2.csv")
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
            mof[name]["functional_groups"] = {None:["N/A"], None:["N/A"]}
            print "Warning: could not find data for %s"%name

    del fnl

def create_a_dataset(csvfile=None, sqlfile=None, metal=None):
    """Create a dataset with pre-defined number of MOFs. Change in code
    if needed.

    """
    mof = CSV(csvfile)
    mof.read_from_csv_multiple()
    fnl = FunctionalGroups(sqlfile)
    pair_csv_fnl(mof, fnl)
    sel = Selector(mof, metal=metal)
    sel.trim_organics()
    sel.trim_non_existing()
    #sel.random_select(exclude=["F", "Cl", "Br", "I", "SO3H", 
    #                                    "NO2", "HCO", "NH2"])
    sel.random_select(inclusive=["F", "Cl", "Br", "I", "SO3H"])

def write_report(directory=None, sqlfile=None, csvfile=None):
    """Write a report on some new uptake data located in 'directory' based
    on some old data in the 'csvfile'

    """
    mofs = CSV(csvfile)
    mofs.read_from_csv_multiple()
    fnl = FunctionalGroups(sqlfile)
    pair_csv_fnl(mofs, fnl)
    data = GrabNewData(mofs, basedir=directory, extended=False)
    data.grab_data()
    data.write_data(filename=directory + ".report.csv")

def main():
    create_a_dataset(csvfile="combined.csv", 
                     sqlfile="vsql.sqlout",
                     metal="V")

    #write_report(directory='wilmer_qeq_gulp_optimized',
    #             sqlfile='allsql.sqlout',
    #             csvfile='top_ranked.csv')

if __name__ == "__main__":
    main()
