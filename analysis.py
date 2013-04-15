#!/usr/bin/env python

"""analysis.py

grab data from different sources and perform selections on that data based on 
certain criteria.

"""
import os
import sys
import random
from math import cos, sin, sqrt, pi, exp, floor
import numpy as np
import scipy
from scipy import stats
from options import Options
import subprocess
import zipfile
import uuid
import itertools
from logging import info, debug, warning, error, critical
sys.path.append("/home/pboyd/codes_in_development/faps")
from faps import PyNiss, Structure, Atom, Cell, Guest, Symmetry 
import pickle

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

class CSV(dict):
    """
    Reads in a .csv file for data parsing.

    """
    def __init__(self, filename, _MOFNAME=True):
        self._columns = {"MOF":"MOFname", "uptake":"mmol/g",
                      "temperature":"T/K", "pressure":"p/bar",
                      "heat of adsorption":"hoa/kcal/mol"}
        self.filename = filename
        if not os.path.isfile(filename):
            error("Could not find the file: %s"%filename)
            sys.exit(1)
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
                warning("Could not find the column %s in the csv file %s "%
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

        warning("Didn't find the data point requested in the csv file %s"%
                self.filename)
        return 0. if _TYPE is "float" else None

    def _parse_by_heading(self):
        """The CSV dictionary will store data to heading keys"""
        filestream = open(self.filename, "r")
        # burn the header, as it's already read
        # if the file is empty append zeroes..
        if self._line_count(self.filename) <= 1: 
            for ind in range(len(self.headings)):
                self.setdefault(self.headings[ind], []).append(0.)
            filestream.close()
            return
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

    def _line_count(self, filename):
        with open(filename) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def _parse_by_mofname(self):
        """The CSV dictionary will have MOFnames as headings and contain
        sub-dictionaries for additional row data.

        """
        filestream = open(self.filename, "r")
        try:
            mofind = self.headings.index(self._columns["MOF"])
        except ValueError:
            error("the csv file %s does not have %s as a column! "%
                    (self.filename, self._columns["MOF"]) + 
                    "EXITING ...")
            sys.exit(0)

        try:
            uptind = self.headings.index(self._columns["uptake"])
        except ValueError:
            warning("the csv file %s does not have %s as a column"%
                    (self.filename, self._columns["uptake"]) +
                    " the qst will be reported as 0.0 kcal/mol")
        try:
            hoaind = self.headings.index(self._columns["heat of adsorption"])
        except ValueError:
            warning("the csv file %s does not have %s as a column"%
                    (self.filename, self._columns["heat of adsorption"]) +
                    " the qst will be reported as 0.0 kcal/mol")
        burn = filestream.readline()
        for line in filestream:
            line = line.strip()
            if line:
                line = line.lstrip("#").split(",")
                mofname = line[mofind].strip()
                mofname = clean(mofname)
                try:
                    uptake = line[uptind]
                except UnboundLocalError:
                    uptake = 0.
                self.setdefault(mofname, {})["mmol/g"] = float(uptake)
                try:
                    hoa = line[hoaind]
                except UnboundLocalError:
                    hoa = 0.
                self.setdefault(mofname, {})["hoa"] = float(hoa)
        filestream.close()


class FunctionalGroups(dict):
    """
    Reads in a .sqlout file and returns a dictionary containing mofnames and
    their functionalizations.

    """

    def __init__(self, filename):
        """Read in all the mofs and store in the dictionary."""
        if not os.path.isfile(filename):
            error("could not find the file: %s"%(filename))
            sys.exit(1)
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
                elif len(dic.keys()) == 0:
                    dic[None] = []
                    dic[False] = []
            # check if the entry already exists!
            if self._check_duplicate(mof):
                if self[mof] == dic:
                    # duplicate
                    debug("Duplicate found %s"%(mof))
                    #pass
                else:
                    warning("duplicate found for %s"%(mof) +
                            " but with different functionalization!")
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
        "Ba" : (10,),
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

    def __init__(self, options, mof_dic):
        self.options = options
        self.mof_dic = mof_dic.copy()
        self.dataset = {'master_list':[]}
        self.ignore = MOFlist(self.options.ignore_list)
        self._assign_metalind()
        if self.options.max_gridpoints == 0:
            self.options.max_gridpoints = None
        self.trim_undesired()

    def _print_info(self):
        info("Inclusive functional groups: %s"%
                (', '.join(self.options.fnl_include)))
        info("Exclude functional groups: %s"%
                (', '.join(self.options.fnl_exclude)))
        info("Partial functional groups: %s"%
                (', '.join(self.options.fnl_partial)))
        info("Weight uptake (gaussian): %s"%self.options.gaussian)
        info("Max grid points: %i"%self.options.max_gridpoints)
        info("Total number of MOFs: %i"%self.options.total_mofs)
        info("Maximum per functional group: %i"%self.options.functional_max)
        info("Maximum per organic linker: %i"%self.options.organic_max)
        info("Maximum per metal index: %i"%self.options.metal_max)
        info("Maximum per topology: %i"%self.options.topology_max)

    def _assign_metalind(self):
        # assert if a metal is selected
        self.metalind = []
        if self.options.metals:
           for metal in self.options.metals:
               # try to extract a metal index from the list of metals
               try:
                   metind = int(metal)
                   if not metind in [ind for sublist in 
                                        self.metal_indices.values()
                                        for ind in sublist]:
                       error("Metal index %i is not in the database!"
                               %metind)
                       sys.exit(1)
                   self.metalind.append(metind)
               except ValueError:
                   try:
                       indices = self.metal_indices[metal]
                       for ind in indices:
                           self.metalind.append(ind)
                   except KeyError:
                       error("ERROR: metal %s is not in the database!"%(metal))
                       sys.exit(1)
        # all metals are fair game
        else:
           for metalind in self.metal_indices.values():
               for i in metalind:
                   self.metalind.append(i)

    def _assign_maxima(self):
        # topologies
        if self.options.topology_max == 0 and len(self.options.topologies) != 0:
            self.options.topology_max = (self.options.total_mofs /
                                         len(self.options.topologies))
        elif self.options.topology_max == 0 and len(self.options.topologies) == 0:
            self.options.topology_max = self.options.total_mofs
        # metal
        if self.options.metal_max == 0 and len(self.metalind) != 0:
            self.options.metal_max = (self.options.total_mofs /
                                      len(self.metalind))
        # organic
        if self.options.organic_max == 0 and len(self.options.org_include) != 0:
            self.options.organic_max = (self.options.total_mofs /
                                        len(self.options.org_include))
        elif self.options.organic_max == 0:
            self.options.organic_max = self.options.total_mofs / 30

        # functional
        if self.options.functional_max == 0:
            self.options.functional_max = (self.options.total_mofs / 20)


    def trim_undesired(self):
        """Trims the MOFs with bad organic linkers from the database.
        This is temporary and should be removed.
        """
        temp_moflist = self.mof_dic.keys()
        for mof in temp_moflist:
            met, org1, org2, top, fnl = parse_mof_data(mof)
            pop = False
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

    def non_existing(self, mof):
        """Checks a MOF against a directory to see if the cif exists."""
        dir = self.options.lookup 
        if not os.path.isdir(self.options.lookup):
            error("Could not find the .cif directory, exiting..")
            sys.exit(1)
        # re-name mof to the out directory standards.
        met, org1, org2, top, fnl = parse_mof_data(mof)
        newmofname = "str_m" + str(met()) + "_o" + str(org1()) + \
                    "_o" + str(org2()) + "_f0_" + top() + ".sym." + \
                    str(fnl()) + ".out.cif"
        if not os.path.isfile(dir + "/" + newmofname):
            return True
        return False

    def top_select(self):
        """Order the mofs by top ranked structures.  store in a dataset
        dictionary and write to a csv file.

        """
        mofnames = self.mof_dic.keys()
        info("Size of the list of MOFs to sample from: %i"%(
            len(mofnames)))
        for mof in mofnames:
            if self._valid_mof(mof):
                self._bin_mof(mof)
                self.dataset['master_list'].append(mof)
        # rank bins by uptake
        limit = self.options.total_mofs if self.options.total_mofs else None
        for bin, dic in self.dataset.items():
            if bin == 'master_list':
                self.dataset[bin] = self._rank_by_uptake(dic)
                self.write_dataset(self.dataset[bin][:limit],
                                   basename="top_ranked") 
            else:
                for key, val in dic.items():
                    self.dataset[bin][key] = self._rank_by_uptake(val)
                    if isinstance(key, tuple):
                        name = "-".join([str(i) for i in key])
                    else:
                        name = str(key)
                    self.write_dataset(self.dataset[bin][key][:limit],
                        basename="%s_%s_top_ranked"%(bin,name))

    def _rank_by_uptake(self, list):
        """Rank by mmol/g uptake."""
        order = []
        for mof in list:
            try:
                uptake = self.mof_dic[mof]['mmol/g']
            except KeyError:
                uptake = 0.

            order.append(tuple([uptake, mof]))
        order = [i for i in reversed(sorted(order))]
        return [i[1] for i in order[:]]

    def _bin_mof(self, mof):
        """Bin mof into different categories"""
        met, org1, org2, top, fnl = parse_mof_data(mof)
        met = met()
        org1 = org1()
        org2 = org2()
        top = top()
        try:
            (fnl_grp1, fnl_grp2) = self.mof_dic[mof]['functional_groups']
            if not fnl_grp1:
                fnl_grp1 = "None"
            if not fnl_grp2:
                fnl_grp2 = "None"
        except ValueError:
            # doesn't find two functional groups in the dictionary
            fnl_grp1, fnl_grp2 = "None", "None"

        org_pair = tuple(sorted([org1, org2]))
        fnl_pair = tuple(sorted([fnl_grp1, fnl_grp2]))
        mofname = "str_m%i_o%i_o%i_%s"%(met, org_pair[0], org_pair[1], top)
        self.dataset.setdefault('mof_id',{}).setdefault(mofname,[]).append(mof)
        self.dataset.setdefault('fnl_p',{}).setdefault(fnl_pair,[]).append(mof)
        self.dataset.setdefault('org_p',{}).setdefault(org_pair,[]).append(mof)
        self.dataset.setdefault('org',{}).setdefault(org1,[]).append(mof)
        self.dataset.setdefault('org',{}).setdefault(org2,[]).append(mof)
        self.dataset.setdefault('met',{}).setdefault(met, []).append(mof)
        if fnl_grp1 == fnl_grp2:
            self.dataset.setdefault('fnl',{}).setdefault(fnl_grp1,[]).append(mof)
        else:
            self.dataset.setdefault('fnl',{}).setdefault(fnl_grp1,[]).append(mof)
            self.dataset.setdefault('fnl',{}).setdefault(fnl_grp2,[]).append(mof)

    def _valid_mof(self, mof):
        """Discriminate a mof based on whatever criteria is set in the
        input file.
        Current checks in order are:

        1  - MOF is not in the ignore list
        2  - MOF exists in the lookup directory
        3  - the reported uptake (mmol/g) is above the cutoff
        4  - the reported uptake (mmol/g) samples a gaussian dist.
        5  - Functional groups are not in the EXCLUDE list
        6  - Functional groups are in the PARTIAL list
        7  - Functional group is in the INCLUDE list
        8  - Organic linkers are not in the EXCLUDE list
        9  - Organic linkers are in the INCLUDE list (read: PARTIAL)
        10 - Metal indices are in the metal INCLUDE list
        11 - The topology is in the topology INCLUDE list
        12 - ESP gridpoints is below a given maximum
        """
        # grab the data for the mof
        met, org1, org2, top, fnl = parse_mof_data(mof)
        org1, org2 = org1(), org2()
        met = met()
        top = top()
        try:
            (fnl1, fnl2) = self.mof_dic[mof]['functional_groups']
        except ValueError:
            # doesn't find two functional groups in the dictionary
            fnl1, fnl2 = None, None
            # return false in the case that discriminating against
            # functional group was requested. False due to the  
            # assumption that the functional group listing for this mof 
            # could not be found.
            if self.options.fnl_include or self.options.fnl_partial:
                return False
 
        uptake = self.mof_dic[mof]['mmol/g']

        # 1 - check if the MOF exists in the ignore list 
        if mof in self.ignore:
            debug("%s was found in the ignore list"%(mof))
            return False

        # 2 - check if the MOF exists in the directory of .cif files
        if self.options.lookup:
            if self.non_existing(mof):
                debug("%s  does not exist!"%(mof))
                return False

        # 3 - check if the uptake reported in mmol/g is above the cutoff
        if uptake < self.options.uptake_cutoff:
            debug("%s has an uptake, %4.2f. This is lower than the cutoff"%
                    (mof, uptake))
            return False

        # 4 - determine if the uptake fits in the gaussian distribution
        if self.options.gaussian:
            # NOTE: this creates a gaussian each instance it's called - 
            # I should implement a pre-computed gaussian curve to reduce 
            # computational expense
            if not weight_by_gaussian(uptake):
                debug("%s has an uptake, %4.2f not sampled by the gaussian."%
                        (mof, uptake))
                return False

        # 5 - check if the functional groups are in the EXCLUDE list
        if (fnl1 in self.options.fnl_exclude) or \
                (fnl2 in self.options.fnl_exclude):
            debug("%s has functional groups in the exclude list"%mof)
            return False

        # 6 - check if the functional groups are in the PARTIAL list
        if self.options.fnl_partial:
            if (fnl1 not in self.options.fnl_partial) and \
                (fnl2 not in self.options.fnl_partial):
                debug("%s has functional groups not in the partial list"%mof)
                return False

        # 7 - check if the functional group is in the INCLUDE list
        if self.options.fnl_include:
            if fnl1 and not fnl2:
                if fnl1 not in self.options.fnl_include:
                    debug("%s has a functional group not in the include list"
                            %(mof))
                    return False
            elif fnl2 and not fnl1:
                if fnl2 not in self.options.fnl_include:
                    debug("%s has a functional group not in the include list"
                            %(mof))
                    return False
            else:
                debug("%s contains two functional groups, the fnl_include" +
                        " list restricts only one per MOF"%(mof))
                return False

        # 8 - check if the organic linkers are in the EXCLUDE list
        if (org1 in self.options.org_exclude) or \
                (org2 in self.options.org_exclude):
            debug("%s has organic linkers in the exclude list"%mof)
            return False

        # 9 - check if the organic linkers are in the INCLUDE list
        if self.options.org_include:
            if (org1 not in self.options.org_include) and \
                (org2 not in self.options.org_include):
                debug("%s has organic linkers not in the include list"%mof)
                return False

        # 10 - check if the metal index is in the list of METALS.
        if met not in self.metalind:
            debug("%s metal linker is not in the metal list"%(mof))
            return False

        # 11 - check if the topology is in the topology list.
        if self.options.topologies:
            if top not in self.options.topologies:
                debug("%s topology is not in the topology include list"%(mof))
                return False

        # 12 - determine if the number of esp grid points is below
        # a maximum (if requested)
        # TODO(pboyd): include some estimating scheme, this takes
        # way to long for the entire database.
        if self.options.report_ngrid:
            try:
                ngrid = self.mof_dic[mof]['ngrid']
            except KeyError:
                ngrid = self.grid_points(mof)
                self.mof_dic[mof]['ngrid'] = ngrid

            if ngrid <= 0:
                debug("gridpoint calculation had errors, ignoring %s"%mof)
                return False

            if ngrid > self.options.max_gridpoints:
                debug("%s contains %i gridpoints.  The max is %i"%
                        (mof,ngrid,self.options.max_gridpoints))
                return False

        # all tests passed, return True.
        return True

    def random_select(self):
        """Select a list of MOFs randomly 

        """
        # check if total_mofs is allocated, if not, then assign default of 100
        if not self.options.total_mofs:
            self.options.total_mofs = 100
        self._assign_maxima()
        self._print_info()
        organics_count, fnl_groups_count = {}, {}
        top_count, met_index_count = {}, {}
        mofcount = 0
        # generate a list of valid mofs which obey the inclusive, partial and 
        # exclude lists.
        moflist = self.mof_dic.keys()
        done = False
        info("Size of the list of MOFs to sample from: %i"%(len(moflist))) 
        while not done:
            try:
                mof = random.choice(moflist)
                moflist.pop(moflist.index(mof))
            except IndexError:
                warning("Sampled all MOFs without completing list! Writing " + 
                      "output file anyways..")
                self.write_dataset(dataset)
                return
            if self._valid_mof(mof):
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
            
                if not org_max and not fnl_max and not top_max and \
                        not met_max:
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
                    self.dataset['master_list'].append(mof)
                    mofcount += 1
            if mofcount >= self.options.total_mofs:
                done = True
        self.write_dataset(self.dataset['master_list'])

    def check_dictionary_counts(self, dictionary, **kwargs):

        if 'fnl_group1' in kwargs.keys():
            maximum = self.options.functional_max
            combine_max = self.options.functional_max / 2
        elif 'organic1' in kwargs.keys():
            maximum = self.options.organic_max 
            combine_max = self.options.organic_max / 2
        elif 'metal' in kwargs.keys():
            maximum = self.options.metal_max
            combine_max = self.options.metal_max
        elif 'topology' in kwargs.keys():
            maximum = self.options.topology_max
            combine_max = self.options.topology_max
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
        5 mmol/g, smeared by 1.3 (width)"""
        gauss = gaussian(a, b, c)
        # probability of accepting the gaussian weight determined
        # by a random number between 0 and a.
        if gauss(uptake) < (random.random()*a):
            return False
        return True

    def grid_points(self, mofname):
        """Determine the maximum number of grid points needed for the esp."""
        # have to source the correct file.
        mofname = clean(mofname)
        dirmof = os.path.join(self.options.lookup, mofname+'.out.cif') 
        ngrid = -1
        if os.path.isfile(dirmof):
            from_cif = CifFile(dirmof)
            ngrid = GrabGridPoints(from_cif.cell,
                                   from_cif.atom_symbols,
                                   from_cif.cart_coordinates)
        return ngrid 

    def write_dataset(self, dataset, basename=None):
        """Writes the data to a file."""
        os.chdir(self.options.job_dir)
        if basename is None:
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
            if self.options.topologies:
                for top in self.options.topologies:
                    basename += "_%s"%top
        count = 0
        order_by_fnl = {}
        filename = create_csv_filename(basename) 
        info("Writing dataset to %s.csv ..."%(filename))
        outstream = open(filename+".csv", "w")
        header="MOFname,mmol/g,hoa/kcal/mol,functional_group1,functional_group2"
        if (self.options.max_gridpoints is not None) or\
                (self.options.report_ngrid):
            header += ",ngrid\n"
        else:
            header += "\n"
        outstream.writelines(header)
        for mof in dataset:
            try:
                fnl1, fnl2 = self.mof_dic[mof]['functional_groups'].keys()
            except (KeyError, ValueError):
                fnl1, fnl2 = None, None
            if not fnl1:
                fnl1 = "None"
            if not fnl2:
                fnl2 = "None"
            order_by_fnl.setdefault((fnl1,fnl2),[]).append(mof)
            try:
                uptake = self.mof_dic[mof]['mmol/g']
            except KeyError:
                uptake = 0.
            try:
                hoa = self.mof_dic[mof]['hoa']
            except KeyError:
                hoa = 0.
            line = "%s,%f,%f,%s,%s"%(mof, uptake, hoa, fnl1, fnl2)
            if (self.options.report_ngrid):
                try:
                    ngrid = self.mof_dic[mof]['ngrid']
                except KeyError:
                    ngrid = self.grid_points(mof)
                line += ",%i\n"%(ngrid)
            else:
                line += "\n"
            outstream.writelines(line)
        outstream.close()
        # write this stuff to the terminal (log file)
        # this was legacy from an older version of the code.
        for fnl, mofs in order_by_fnl.items():
            debug("Functional groups: %s, %s"%(fnl[0], fnl[1]))
            for mof in mofs:
                debug("     " + mof)

        info("Done.")

class GrabNewData(object):
    """Takes a mof dictionary and adds new uptake data to it."""
    def __init__(self, options, mofs, extended = False):
        self.mofs = mofs
        self.options = options
        # determines if hydrogen replacements are shown in the output
        self.extended = extended
        # base directory containing all the faps job info.
        if self.options.directory:
            if os.path.exists(self.options.directory):
                self.basedir = self.options.directory
            else:
                error("The directory specified in the input file "+
                        "could not be found: %s"%(self.options.directory))
                sys.exit(1)
        else:
            error("No directory with new data specified in the "+
                    "input file!")
            sys.exit(1)

    def grab_data(self, temp=298.0, press=0.15):
        """Descends into directories and grabs appropriate data."""
        os.chdir(self.basedir)
        all_list = os.listdir('.')
        directories = [i for i in all_list if os.path.isdir(i)]
        for mof in self.mofs.keys():
            if mof in directories:
                info("grabbing data from %s"%mof)
                os.chdir(mof)
                if os.path.isfile(mof + "-CO2.csv"):
                    data = CSV(mof + "-CO2.csv", _MOFNAME=False)
                    new_uptake = data.obtain_data("mmol/g",
                                                  temperature=temp,
                                                  pressure=press)
                    new_hoa = data.obtain_data("heat of adsorption",
                                               temperature=temp,
                                               pressure=press)
                else:
                    warning("could not find %s-CO2.csv"%(mof))
                    new_uptake = 0.
                    new_hoa = 0.
                os.chdir('..')
            else:
                warning("could not find %s in the directory %s"%(
                        mof, self.basedir))
                new_uptake = 0.
                new_hoa = 0.
            self.mofs[mof]['new_uptake'] = new_uptake
            self.mofs[mof]['new_hoa'] = new_hoa
        os.chdir(self.options.job_dir)

    def write_data(self, filename="default.report.csv"):
        """Write all MOF data to a csv file."""
        os.chdir(self.options.job_dir)
        basename = clean(filename)
        count = 0
        # make sure there's no overwriting, goes up to 99
        filename = create_csv_filename(basename)
        info("Writing report to %s.csv ..."%(os.path.basename(filename)))
        outstream = open(filename + ".csv", "w")
        header = "MOFname,mmol/g,hoa/kcal/mol,functional_group1,functional_group2"
        if self.extended:
            header += ",grp1_replacements,grp2_replacements"
        if self.options.report_ngrid:
            header += ",ngrid"
        header += "\n"
        outstream.writelines(header)
        for mof in self.mofs.keys():
            new_uptake = self.mofs[mof]['new_uptake']
            new_hoa = self.mofs[mof]['new_hoa']
            try:
                replaced_groups = self.mofs[mof]['functional_groups']
            except KeyError:
                replaced_groups = {False:[], None:[]}
            names = replaced_groups.keys()
            names.sort()
            fnl_grp1 = names[0]
            fnl_grp2 = names[1]
            rep_1 = " ".join(replaced_groups[fnl_grp1])
            rep_2 = " ".join(replaced_groups[fnl_grp2])
            line = "%s,%f,%f,%s,%s"%(mof, new_uptake, new_hoa, 
                                     fnl_grp1, fnl_grp2)
            if self.extended:
                line += "%s,%s"%(rep_1, rep_2)
            if self.options.report_ngrid:
                line += ",%i"%(self.grid_points(mof))
            line += "\n"
            outstream.writelines(line)
        info("Done.")
        outstream.close()

    def grid_points(self, mofname):
        """Determine the maximum number of grid points needed for the esp."""
        # have to source the correct file.
        mofname = clean(mofname)
        dirmof = os.path.join(self.options.lookup, mofname+'.out.cif') 
        ngrid = -1
        if os.path.isfile(dirmof):
            from_cif = CifFile(dirmof)
            ngrid = GrabGridPoints(from_cif.cell,
                                   from_cif.atom_symbols,
                                   from_cif.cart_coordinates)
        return ngrid 


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
            warning("ERROR in EGULP calculation!")
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

class JobHandler(object):
    """Reads the data in options and directs the program to the
    appropriate job.

    """

    def __init__(self, options):
        self.options = options
        self._check_inconsistencies()
        self._direct_job()

    def _check_job_type(self):
        if not self.options.combine and not self.options.cmd_opts.combine and\
                not self.options.dataset and not self.options.report and \
                not self.options.extract and not self.options.report and not \
                self.options.comparison and not self.options.csv_info_file:
            error("no job type requested!")
            sys.exit(1)

    def _check_inconsistencies(self):
        """General function to check if requested job directives conflict."""
        # check if the necessary boolean was requested for the gridpoints.
        if self.options.max_gridpoints > 0 and not self.options.report_ngrid:
            error("The gridpoint max cutoff was set to %i but report_ngrid was"
                    +" not enabled in the input file, set report_ngrid = True!")
            sys.exit(1)
        
        # check if the number of gridpoints max is really low.
        if self.options.report_ngrid and self.options.dataset and \
                self.options.max_gridpoints < 50:
            warning("The max number of gridpoints was set to %i, "%(
                self.options.max_gridpoints) + "this is likely too low"
                + " set max_gridpoints to > 50,000.  Will continue anyways..")
        if self.options.dataset and not self.options.ignore_list:
            ignore = str(uuid.uuid1())
            warning("No ignore_list of MOFs chosen, will create a "+
                    "blank file %s.ignore to use. Remove "%ignore +
                    "after program termination.")
            self.options.ignore_list = "%s.ignore"%ignore
            ignorestream = open("%s.ignore"%ignore, "w")
            ignorestream.writelines("MOFname,")
            ignorestream.close()

    def _direct_job(self):
        self._check_job_type()
        if self.options.combine or self.options.cmd_opts.combine:
            comb_list = self.options.combine + self.options.cmd_opts.combine
            info("Combining csv files: %s"%(', '.join(comb_list)))
            self.combine_csvs(*set(comb_list))
            return
        # if combine is requested, skip all this
        if self.options.dataset or self.options.report or \
                self.options.extract or self.options.report:
            # both need a list of MOFs with functional groups
            if not self.options.csv_file:
                error("no .csv file specified in the input")
                sys.exit(1)
            if not self.options.sql_file:
                error("no .sqlout file specified in the input")
                sys.exit(1)

            self.mofs = CSV(self.options.csv_file)
            self.fnl_dic = FunctionalGroups(self.options.sql_file)
            pair_mofs_fnl(self.mofs, self.fnl_dic)
        
        if self.options.extract:
            info("Extracting info from csv file %s"%(self.options.csv_file))
            self.extract_report()

        elif self.options.comparison:
            if not self.options.compare_csv_1 and not \
                    self.options.compare_csv_2:
                error("Please enter the two .csv files to compare in the "+
                        "input file.")
                sys.exit(1)
            if not self.options.compare_location_1 and not \
                    self.options.compare_location_2:
                error("Please enter the location of the job runs for both "+
                        "csv files.")
                sys.exit(1)
            if not self.options.sql_file:
                error("no .sqlout file specified in the input")
                sys.exit(1)
            info("Comparing two .csv files %s, and %s"%(
                os.path.basename(self.options.compare_csv_1), 
                os.path.basename(self.options.compare_csv_2)))
            self.run_comparison()

        elif self.options.dataset:
            if self.options.top_ranked:
                info("Generating top ranked structures..")
                self.generate_top_structures()
            else:
                info("Generating randomly distributed structures..")
                self.create_a_dataset()

        elif self.options.report:
            info("Report requested")
            self.write_report()

        elif self.options.csv_info_file:
            info("Compiling info for %s"%(self.options.csv_info_file))
            job = GetInfo(self.options)
            info("Found %i MOFs in the file."%(len(job.moflist)))
            job.write_csv_file()

    def extract_report(self):
        """Write a bunch of reports based on the stats of the MOFs in 
        the csv file.

        """
        ext = Extract(self.options, self.mofs)
        ext.write_overall_csv()
        ext.write_specific_csv()

    def create_a_dataset(self):
        """Create a dataset with pre-defined number of MOFs. Change in code
        if needed.
    
        """
        sel = Selector(self.options, self.mofs)
        sel.random_select()

    def run_comparison(self):
        """Run a comparison between two .csv files.. generally used for large
        amounts of data."""
        job = Comparison(self.options)
        job.bin_data()
        job.write_output_data()

    def write_report(self):
        """Write a report on some new uptake data located in 'directory' based
        on some old data in the 'csvfile'

        """
        data = GrabNewData(self.options, self.mofs)
        data.grab_data()
        # first strip the .csv file of all pre_defined directories
        basename = clean(os.path.basename(self.options.csv_file))
        # then create the base file name
        filename = '.'.join([basename,
                         os.path.basename(self.options.directory)])
        data.write_data(filename)

    def generate_top_structures(self):
        sel = Selector(self.options, self.mofs)
        sel.top_select()

    def combine_csvs(self, *args):
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
                hoa = dic.pop('hoa')
                dic['mmol/g-%s'%basenames[ind]] = uptake
                dic['hoa/kcal/mol-%s'%basenames[ind]] = hoa
                merged.setdefault(key, {})
                merged[key].update(dic)
        # combine uptake columns and hoa columns
        upt_corr = [[] for i in range(len(basenames))]
        hoa_corr = [[] for i in range(len(basenames))]
        for mof, value in merged.items():
            if all(['mmol/g-%s'%i in value.keys() for i in basenames]):
                if all([value['mmol/g-%s'%i] > 0. for i in basenames]):
                    [upt_corr[basenames.index(i)].append(value['mmol/g-%s'%i]) for 
                            i in basenames]
                else:
                    warning("MOF %s not included in correlation calculations "%mof
                        + "because of 0 value for uptake")
                if all([value['hoa/kcal/mol-%s'%i] > 0. for i in basenames]):
                    [hoa_corr[basenames.index(i)].append(value['hoa/kcal/mol-%s'%i]) for 
                            i in basenames]
                else:
                    warning("MOF %s not included in correlation calculations "%mof
                        + "because of 0 value for heat of adsorption")

            else:
                warning("MOF %s is not included in correlation calculations "%mof +
                    "because it couldn't be found in one or more csv files.")
        for upt1, upt2 in itertools.combinations(upt_corr, 2):
            if upt1 and upt2:
                ind1 = upt_corr.index(upt1)
                ind2 = upt_corr.index(upt2)
                info("Uptake correlation coefficients for %s, and %s:"%(
                     basenames[ind1], basenames[ind2]))
                rho, pval = stats.spearmanr(upt1, upt2)
                info("Spearman rank: %7.5f"%rho)
                pears, p2 = stats.pearsonr(upt1, upt2)
                info("Pearson correlation: %7.5f"%pears)
        for hoa1, hoa2 in itertools.combinations(hoa_corr, 2):
            if hoa1 and hoa2:
                ind1 = hoa_corr.index(hoa1)
                ind2 = hoa_corr.index(hoa2)
                info("Hoa correlation coefficients for %s, and %s:"%(
                     basenames[ind1], basenames[ind2]))
                rho, pval = stats.spearmanr(hoa1, hoa2)
                info("Spearman rank: %7.5f"%rho)
                pears, p2 = stats.pearsonr(hoa1, hoa2)
                info("Pearson correlation: %7.5f"%pears)

        names = {}
        [names.setdefault(i.split('.')[0], 0) for i in basenames]
        name = names.keys()[0]
        name += ".merged"
        write_csv(name, merged)

class GetInfo(object):
    """Gets as much information as I could think of from a list of MOFs
    provided by the user.

    """
    def __init__(self, options):
        self.options = options
        self.fnl_dic =  FunctionalGroups(self.options.sql_file)
        self.mof_dic = CSV(self.options.csv_file)
        self.moflist = self._get_list_of_mofs()
        pair_mofs_fnl(self.mof_dic, self.fnl_dic)

    def _get_list_of_mofs(self):
       return MOFlist(self.options.csv_info_file) 

    def write_csv_file(self):
        filename = create_csv_filename(clean(self.options.csv_info_file))
        info("Writing info file as %s.csv"%(filename))
        csvstream = open(filename+".csv", "w")
        header = "MOFname,mmol/g,hoa/kcal/mol,functional_group1,functional_group2\n"
        csvstream.writelines(header)
        for mof in self.moflist:
            try:
                mmol_g = self.mof_dic[mof]['mmol/g']
            except KeyError:
                mmol_g = 0.
            try:
                hoa = self.mof_dic[mof]['hoa']
            except KeyError:
                hoa = 0.
            try:
                fnl_1 = self.mof_dic[mof]['functional_groups'].keys()[0]
                fnl_2 = self.mof_dic[mof]['functional_groups'].keys()[1]
            except KeyError:
                fnl_1 = None
                fnl_2 = None
            fnl_1 = fnl_1 if fnl_1 else None
            fnl_2 = fnl_2 if fnl_2 else None
            csvstream.writelines("%s,%f,%f,%s,%s\n"%(mof,mmol_g,hoa,fnl_1,fnl_2))
        csvstream.close()

class Extract(object):
    """Extract as much information from a list of mofnames, and provide
    an ammended report with the findings.

    NOTE: the file writing is really ugly... if I have time, repair.
    """
    def __init__(self, options, mofs):
        self.options = options
        self.mofs = mofs.copy() 
        self.fnl_count = {}
        self.fnl_pair = {}
        self.organic_count = {}
        self.metal_count = {}
        self.org_pair_count = {}
        self.full_mof_count = {}
        self._analyse_data()

    def _increment_counts(self, dic, *args):
        for arg in args:
            dic.setdefault(arg, 0)
            dic[arg] += 1

    def _grid_points(self, mofname):
        """Determine the maximum number of grid points needed for the esp."""
        # have to source the correct file.
        mofname = clean(mofname)
        dirmof = os.path.join(self.options.lookup, mofname+'.out.cif') 
        ngrid = -1
        from_cif = CifFile(dirmof)
        ngrid = GrabGridPoints(from_cif.cell,
                               from_cif.atom_symbols,
                               from_cif.cart_coordinates)
        return ngrid

    def _analyse_data(self):
        for mof, info in self.mofs.items():
            met, org1, org2, top, junk = parse_mof_data(mof) 
            fnls = info['functional_groups'].keys()
            fnl_pair = sorted([i for i in fnls])
            org_pair = sorted([org1(), org2()])
            [self._increment_counts(self.fnl_count, i) for i in 
                    fnl_pair]
            if len(fnl_pair) > 1:
                self._increment_counts(self.fnl_pair, tuple(fnl_pair))
            [self._increment_counts(self.organic_count, i) for i in 
                     org_pair]
            if len(org_pair) > 1:
                self._increment_counts(self.org_pair_count, tuple(org_pair))
            self._increment_counts(self.metal_count, met())

            fullmof = "str_m%i_o%i_o%i_%s"%(met(),org1(),org2(),top())
            self._increment_counts(self.full_mof_count, fullmof)

    def write_overall_csv(self):

        os.chdir(self.options.job_dir)
        basename = clean(self.options.csv_file)
        count = 0
        files = []
        info("Writing overall reports...")
        # first file prints out individual functional group stats
        fnlfile = basename + ".fnl_grps"
        # create original filename
        fnlfile = create_csv_filename(fnlfile)
        info("Writing %s.csv"%fnlfile)
        # store filename for archiving later
        files.append("%s.csv"%fnlfile)
        outstream = open(fnlfile + ".csv", "w")
        # write data to file
        outstream.writelines("#functional_group,count\n")
        for fnl, count in self.fnl_count.items():
            outstream.writelines("%s,%i\n"%(fnl,count))
        outstream.close()
        # second file prints out pair functional group stats
        fnlpfile = basename + ".fnl_pairs"
        fnlpfile = create_csv_filename(fnlpfile)
        info("Writing %s.csv"%fnlpfile)
        files.append("%s.csv"%fnlpfile)
        outstream = open(fnlpfile + ".csv", "w")
        outstream.writelines("#functional_groups,count\n")
        for pair, count in self.fnl_pair.items():
            outstream.writelines("(%s %s),%i\n"%(pair[0], pair[1], count))
        outstream.close()
        # third file prints out individual organic group stats
        orgfile = basename + ".org"
        orgfile = create_csv_filename(orgfile)
        info("Writing %s.csv"%orgfile)
        files.append("%s.csv"%orgfile)
        outstream = open(orgfile + ".csv", "w")
        outstream.writelines("#organic_index,count\n")
        for org, count in self.organic_count.items():
            outstream.writelines("%i,%i\n"%(org, count))
        outstream.close()
        # fourth file prints out organic pair stats
        orgpfile = basename + ".org_pairs"
        orgpfile = create_csv_filename(orgpfile)
        info("Writing %s.csv"%orgpfile)
        files.append("%s.csv"%orgpfile)
        outstream = open(orgpfile + ".csv", "w")
        outstream.writelines("#organic_indices,count\n")
        for orgs, count in self.org_pair_count.items():
            outstream.writelines("(o%i o%i),%i\n"%(orgs[0], orgs[1], count))
        outstream.close()
        # fifth file prints out metal index stats
        metfile = basename + ".metal"
        metfile = create_csv_filename(metfile)
        info("Writing %s.csv"%metfile)
        files.append("%s.csv"%metfile)
        outstream = open(metfile + ".csv", "w")
        outstream.writelines("#metal_index,count\n")
        for met, count in self.metal_count.items():
            outstream.writelines("%i,%i\n"%(met, count))
        outstream.close()
        # sixth file prints out whole mof stats
        moffile = basename + ".mofname"
        moffile = create_csv_filename(moffile)
        info("Writing %s.csv"%moffile)
        files.append("%s.csv"%moffile)
        outstream = open(moffile + ".csv", "w")
        outstream.writelines("#MOFname,count\n")
        for mof, count in self.full_mof_count.items():
            outstream.writelines("%s,%i\n"%(mof, count))
        outstream.close()
        info("Done.")

        # Archive data
        zipname = basename + ".overall_reports"
        zipname = create_csv_filename(zipname, extension=".zip")
        info("Archiving to %s.zip ..."%zipname)
        zip = zipfile.ZipFile(zipname + ".zip", "w")
        for f in files:
            zip.write("%s"%f)
        zip.close()
        info("Done.")
        # Clear filenames
        for f in files:
            os.remove(f)

    def write_specific_csv(self):
        os.chdir(self.options.job_dir)
        basename = clean(self.options.csv_file)
        count = 0
        filename = basename + ".specific_report"
        filename = create_csv_filename(filename)
        outstream = open(filename + ".csv", "w")
        info("Writing specific report to %s.csv"%filename)
        header = "MOFname,functional_group1,functional_group2," +\
                "mof_occur,metal_occur,org1_occur,org2_occur," +\
                "org_pair_occur,fnl1_occur,fnl2_occur," +\
                "fnl_pair_occur"

        if self.options.report_ngrid:
            header += ",ngrid"
        header += "\n"
        outstream.writelines(header)
        for mof, vals in self.mofs.items():
            met, org1, org2, top, fnl = parse_mof_data(mof)
            fnls = vals['functional_groups'].keys()
            # grab all the data for each
            fnls = sorted([i for i in fnls if i])
            if len(fnls) == 1:
               fnl1 = fnls[0]
               fnl2 = None
               fnl1_occur = self.fnl_count[fnl1]
               fnl_pair_occur = self.fnl_count[fnl1]
               fnl2_occur = 0
            elif len(fnls) == 2:
               fnl1, fnl2 = fnls
               fnl1_occur = self.fnl_count[fnl1]
               fnl2_occur = self.fnl_count[fnl2]
               pair = sorted([fnl1, fnl2])
               try:
                   fnl_pair_occur = self.fnl_pair[tuple(pair)]
               except KeyError:
                   fnl_pair_occur = 0
            met_occur = self.metal_count[met()]
            orgs = sorted([org1(), org2()])
            if len(set(orgs)) == 1:
                org1_occur = self.organic_count[org1()]
                org2_occur = 0
                org_pair_occur = self.organic_count[org1()]
            elif len(set(orgs)) == 2:
                org1_occur = self.organic_count[org1()]
                org2_occur = self.organic_count[org2()]
                org_pair_occur = self.org_pair_count[tuple(orgs)]

            mof_abbrev = "str_m%i_o%i_o%i_%s"%(met(),org1(),org2(),top())
            mof_occur = self.full_mof_count[mof_abbrev]
            line = ("%s,%s,%s,%i,%i,%i,%i,%i,%i,%i,%i"
                    %(mof,fnl1,fnl2,mof_occur,met_occur,org1_occur,
                        org2_occur,org_pair_occur,fnl1_occur,fnl2_occur,
                        fnl_pair_occur))
            if self.options.report_ngrid:
                line += ",%i"%(self._grid_points(mof))
            line += "\n"
            outstream.writelines(line)
        outstream.close()
        info("Done.")


class Comparison(object):
    """Compares two runs on a database with different charge set parameters.
    The user defined value will be compared and binned by a user defined
    deviation interval, For each interval the charges assigned to each
    atom will be reported, as well as the population of functional groups and
    building units.

    """

    def __init__(self, options):

        self.options = options
        self.mof_dic = {}
        self.fnl_dic = FunctionalGroups(self.options.sql_file)
        self._join_dictionaries()
        self._determine_min_max()
        self._init_bin_dics()

    def _join_dictionaries(self):
        mof_dic1 = CSV(self.options.compare_csv_1)
        mof_dic2 = CSV(self.options.compare_csv_2)
        for mof in mof_dic1.keys():
            mof = clean(mof)
            try:
                mof_dic2[mof]
                mofname = clean(mof)
                self.mof_dic[mofname] = {}
                # add the first entries
                for key, val in mof_dic1[mof].items():
                    self.mof_dic[mofname][key+".1"] = val
                # add the second entries
                for key, val in mof_dic2[mof].items():
                    self.mof_dic[mofname][key+".2"] = val
                mof_dic2.pop(mof)
            except KeyError:
                warning("No %s found in the file %s!"%(mof,
                            self.options.compare_csv_2))
                mof_dic1.pop(mof)
        for key in mof_dic2.keys():
            warning("No %s found in the file %s!"%(key, 
                          self.options.compare_csv_1))
        pair_mofs_fnl(self.mof_dic, self.fnl_dic)
        del self.fnl_dic
        del mof_dic1
        del mof_dic2

    def _determine_min_max(self):

        self.min = 1
        self.max = 0
        for mof in self.mof_dic.keys():
            diff = self._get_difference(mof)
            if diff < self.min:
                self.min = diff
            if diff > self.max:
                self.max = diff
        info("Bin intervals: %3.5f"%((self.max - self.min)/
            float(self.options.bin_interval)))

    def _init_bin_dics(self):
        self.atomic_charge_difference = {}
        self.atomic_charge_1 = {}
        self.atomic_charge_2 = {}
        self.mofname_bin = {}
        self.met_bin = {}
        self.org_pair_bin = {}
        self.org_bin = {}
        self.fnl_pair_bin = {}
        self.fnl_bin = {}
        for i in range(self.options.bin_interval):
            self.atomic_charge_difference.setdefault(i, {}) 
            self.atomic_charge_1.setdefault(i, {})
            self.atomic_charge_2.setdefault(i, {}) 
            self.mofname_bin.setdefault(i, [])
            self.met_bin.setdefault(i, {})
            self.org_pair_bin.setdefault(i, {})
            self.org_bin.setdefault(i, {})
            self.fnl_pair_bin.setdefault(i, {})
            self.fnl_bin.setdefault(i, {})

    def bin_data(self):
        self.mof_count = 0
        for mof in self.mof_dic.keys():
            if self._validate_data(mof):
                self.mof_count += 1 
                diff = self._get_difference(mof)
                bin_num = self._determine_bin(diff)
                charge_dict1 = self._get_charges(mof, 1)
                charge_dict2 = self._get_charges(mof, 2)
                for atom, (type1, charge1) in charge_dict1.items():
                    try:
                        (type2, charge2) = charge_dict2[atom]
                        difference = (charge1 - charge2)
                        # store the difference, and the charges for each atom
                        # type in a separate dictionary
                        self.atomic_charge_difference[bin_num].setdefault(type1,
                            []).append(difference)
                        self.atomic_charge_1[bin_num].setdefault(type1, 
                            []).append(charge1)
                        self.atomic_charge_2[bin_num].setdefault(type2, 
                            []).append(charge2)
                    except KeyError:
                        warning("Could not find %s in the second charge set"%(atom))
                met, o1, o2, top, junk = parse_mof_data(mof)
                try:
                    (fnl1, fnl2) = self.mof_dic[mof]["functional_groups"].keys()
                except ValueError:
                    (fnl1, fnl2) = (None, None)
                #YEUCH, ugly stuff next
                self.met_bin[bin_num].setdefault(met(), 0)
                self.met_bin[bin_num][met()] += 1
                self.mofname_bin[bin_num].append(mof)
                opair = tuple(sorted([o1(), o2()]))
                self.org_pair_bin[bin_num].setdefault(opair, 0)
                self.org_pair_bin[bin_num][opair] += 1
                self.org_bin[bin_num].setdefault(o1(), 0)
                self.org_bin[bin_num].setdefault(o2(), 0)
                if o1() == o2():
                    self.org_bin[bin_num][o1()] += 1
                else:
                    self.org_bin[bin_num][o1()] += 1
                    self.org_bin[bin_num][o2()] += 1
                fpair = tuple(sorted([fnl1, fnl2]))
                self.fnl_pair_bin[bin_num].setdefault(fpair, 0)
                self.fnl_bin[bin_num].setdefault(fnl1, 0)
                self.fnl_bin[bin_num].setdefault(fnl2, 0)
                self.fnl_pair_bin[bin_num][fpair] += 1
                self.fnl_bin[bin_num][fnl1] += 1
                self.fnl_bin[bin_num][fnl2] += 1
        del self.mof_dic

    def write_output_data(self):
        """This is ugly."""
        # create directories in the current directory for:
        # charge difference, charges of 1, charges of 2,
        # mof list, functional, building units
        info("Writing .csv files")
        os.makedirs(os.path.join(self.options.job_dir, "charge_1"))
        os.makedirs(os.path.join(self.options.job_dir, "charge_2"))
        os.makedirs(os.path.join(self.options.job_dir, "charge_diff"))
        os.makedirs(os.path.join(self.options.job_dir, "mof"))
        os.makedirs(os.path.join(self.options.job_dir, "functional"))
        os.makedirs(os.path.join(self.options.job_dir, "building_units"))
        for bin_int in range(self.options.bin_interval):
            #CHARGE 1
            basepath = os.path.join(self.options.job_dir, "charge_1")
            csvfile = open(os.path.join(basepath, "%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,%s,%s\n"%("Atom_type","mean_chg","stdev"))
            for atomtype in self.atomic_charge_1[bin_int].keys():
                mean_chge = scipy.mean(self.atomic_charge_1
                                       [bin_int][atomtype])
                stdev_chge = scipy.std(self.atomic_charge_1
                                       [bin_int][atomtype])
                csvfile.writelines("%s,%f,%f\n"%
                        (atomtype, mean_chge, stdev_chge))
            csvfile.close()
            #CHARGE 2
            basepath = os.path.join(self.options.job_dir, "charge_2")
            csvfile = open(os.path.join(basepath, "%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,%s,%s\n"%("Atom_type","mean_chg","stdev"))
            for atomtype in self.atomic_charge_2[bin_int].keys():
                mean_chge = scipy.mean(self.atomic_charge_2
                                       [bin_int][atomtype])
                stdev_chge = scipy.std(self.atomic_charge_2
                                       [bin_int][atomtype])
                csvfile.writelines("%s,%f,%f\n"%
                        (atomtype, mean_chge, stdev_chge))
            csvfile.close()
            #CHARGE DIFFERENCE
            basepath = os.path.join(self.options.job_dir, "charge_diff")
            csvfile = open(os.path.join(basepath, "%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,%s,%s\n"%("Atom_type","mean_diff","stdev"))
            for atomtype in self.atomic_charge_2[bin_int].keys():
                mean_chge = scipy.mean(self.atomic_charge_difference
                                       [bin_int][atomtype])
                stdev_chge = scipy.std(self.atomic_charge_difference
                                       [bin_int][atomtype])
                csvfile.writelines("%s,%f,%f\n"%
                        (atomtype, mean_chge, stdev_chge))
            csvfile.close()
            #MOF STORAGE
            basepath = os.path.join(self.options.job_dir, "mof")
            csvfile = open(os.path.join(basepath, "%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,\n"%("MOFname"))
            for mof in self.mofname_bin[bin_int]:
                csvfile.writelines("%s,\n"%mof)
            csvfile.close()
            # FUNCTIONAL GROUPS
            basepath = os.path.join(self.options.job_dir, "functional")
            csvfile = open(os.path.join(basepath, "fnl_%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,%s\n"%("Functional_group","rep"))
            for fnl,count in self.fnl_bin[bin_int].items():
                percent = float(count) / float(self.mof_count)
                csvfile.writelines("%s,%f\n"%(fnl,percent)) 
            csvfile.close()
            csvfile = open(os.path.join(basepath, "fnl_pair_%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,%s\n"%("Functional_pair","rep"))
            for fnlpair, count in self.fnl_pair_bin[bin_int].items():
                fnlpair = list(fnlpair)
                percent = float(count) / float(self.mof_count)
                for i in range(2):
                    if not fnlpair[i]:
                        fnlpair[i] = "None"
                csvfile.writelines("%s,%f\n"%(" ".join(fnlpair),percent))
            csvfile.close()
            # BUILDING UNITS
            basepath = os.path.join(self.options.job_dir, "building_units")
            csvfile = open(os.path.join(basepath, "org_pair_%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,%s\n"%("Organic_pair","rep"))
            for orgpair, count in self.org_pair_bin[bin_int].items():
                percent = float(count) / float(self.mof_count)
                orgpair = [str(i) for i in orgpair]
                csvfile.writelines("%s,%f\n"%(" ".join(orgpair),percent))
            csvfile.close()
            csvfile = open(os.path.join(basepath, "org_%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,%s\n"%("Organic_ind","rep"))
            for org, count in self.org_bin[bin_int].items():
                percent = float(count) / float(self.mof_count)
                csvfile.writelines("%i,%f\n"%(org,percent))
            csvfile.close()
            csvfile = open(os.path.join(basepath, "met_%i.csv"%(bin_int)),'w')
            csvfile.writelines("%s,%s\n"%("Metal_ind","rep"))
            for met, count in self.met_bin[bin_int].items():
                percent = float(count) / float(self.mof_count)
                csvfile.writelines("%i,%f\n"%(met,percent))
            csvfile.close()

    def _determine_bin(self, val):
        interval = self.max - self.min
        bin_freq = interval / float(self.options.bin_interval)
        if val == self.max:
            return (self.options.bin_interval - 1)
        return int(floor((val - self.min)/bin_freq))

    def _validate_data(self, mof):
        try:
            self.mof_dic[mof][self.options.comparison_value + '.1']
            pass1 = True
        except KeyError:
            warning("No %s value for %s in the first csv file!"%(
                self.options.comparison_value, mof))
            pass1 = False
        try:
            self.mof_dic[mof][self.options.comparison_value + '.2']
            pass2 = True
        except KeyError:
            warning("No %s value for %s in the second csv file!"%(
                self.options.comparison_value, mof))
            pass2 = False
        return (pass1 and pass2)

    def _get_difference(self, mof):
        """Evaluate the difference between the desired value."""
        return (self.mof_dic[mof][self.options.comparison_value + '.1'] -
                self.mof_dic[mof][self.options.comparison_value + '.2'])

    def _get_charges(self, mof, run=1):
        """return the charges assigned to each atom in the MOF."""
        if run == 1:
            basedir = self.options.compare_location_1
        elif run == 2:
            basedir = self.options.compare_location_2
        filename = self._locate(basedir, "%s.niss"%(mof))
        if filename:
            niss = open(filename,'rb')
            sim = pickle.load(niss)
            niss.close()
            chargedic = {}
            for atom in sim.structure.atoms:
                chargedic[atom.site] = (atom.type, atom.charge)
            return chargedic 
        else:
            warning("Couldn't find the charges file for %s"%mof)
            return {} 

    def _locate(self, basedir, file):
        """Crawl a directory for a particular file, return the location."""
        for root, dirs, files in os.walk(basedir):
            if file in files:
                return os.path.join(root,file)


class CifFile(object):
    """This class will grab data from the cif file such as the
    cartesian atom coordinates and the cell vectors.
    
    """

    def __init__(self, moffilename):
        """Read in the mof file and parse the data to write a .geo file
        for the egulp method.

        """
        self.mofdata = self.parse_cif(moffilename)
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

    def parse_cif(self, moffilename):
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


def pair_mofs_fnl(mof, fnl):
    """Join csv and fnl to one dictionary."""
    for name in mof.keys():
        try:
            mof[name]["functional_groups"] = fnl[name]

        except KeyError:
            mof[name]["functional_groups"] = {False:[], None:[]}
            debug("could not find data for %s"%name)
    del fnl

def write_csv(basename, dic):
    """Writes a specific csv where the dictionary contains a mof name
    followed by a dictionary of values.

    """
    filename = create_csv_filename(basename)
    info("Writing to %s.csv"%filename)
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
    info("Done.")
    outstream.close()

def create_csv_filename(basename, extension=".csv"):
    count = 0
    filename = basename
    # make sure there's no overwriting, goes up to 99
    while os.path.isfile(filename + extension):
        count += 1
        filename = basename + ".%02d"%(count)
    return filename

def clean(name):
    if name.startswith('./run_x'):
        name = name[10:]
    if name.endswith('.cif'):
        name = name[:-4]
    elif name.endswith('.niss'):
        name = name[:-5]
    elif name.endswith('.out-CO2.csv'):
        name = name[:-12]
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
    options = Options()
    job = JobHandler(options)

if __name__ == "__main__":
    main()
