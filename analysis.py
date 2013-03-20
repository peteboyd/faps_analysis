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
from options import Options
import subprocess
import zipfile
import uuid
import itertools
from logging import info, debug, warning, error, critical

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
                      "temperature":"T/K", "pressure":"p/bar"}
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
            error("the csv file %s does not have %s as a column! "%
                    (self.filename, self._columns["MOF"]) + 
                    "EXITING ...")
            sys.exit(0)

        try:
            uptind = self.headings.index(self._columns["uptake"])
        except ValueError:
            warning("the csv file %s does not have %s as a column"%
                    (self.filename, self._columns["uptake"]) +
                    " the uptake will be reported as 0.0 mmol/g")
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
        used_mofs = MOFlist(self.options.ignore_list)
        # remove mofs in the dictionary which have already been used
        for mof in used_mofs:
            try:
                self.mof_dic.pop(mof)
            except KeyError:
                pass
        self._assign_metalind()
        self._assign_maxima()
        top_bool = True if self.options.topologies else False
        self._print_info()
        if self.options.max_gridpoints == 0:
            self.options.max_gridpoints = None
        self.trim_undesired(_TOPOLOGY=top_bool)

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
                if top() not in self.options.topologies:
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
        dir = self.options.lookup 
        info("Checking %s for all mofs..."%(dir))
        if not os.path.isdir(self.options.lookup):
            error("Problem finding the directory, exiting..")
            sys.exit(1)
        temp_moflist = self.mof_dic.keys()
        for mof in temp_moflist:
            # re-name mof to the out directory standards.
            met, org1, org2, top, fnl = parse_mof_data(mof)
            newmofname = "str_m" + str(met()) + "_o" + str(org1()) + \
                    "_o" + str(org2()) + "_f0_" + top() + ".sym." + \
                    str(fnl()) + ".out.cif"
            if not os.path.isfile(dir + "/" + newmofname):
                self.mof_dic.pop(mof)
                #warning("%s  does not exist!"%(mof))

    def top_select(self):
        """Order the mofs by top ranked structures.  store in a dataset
        dictionary and write to a csv file.

        """
        dataset = {}
        moflist = self._gen_moflist()
        info("Size of the list of MOFs to sample from: %i"%(len(moflist)))
        _used_mofs = []
        if self.options.fnl_include or self.options.fnl_partial:
            for group in self.options.fnl_include + self.options.fnl_partial:
                # obtain list of mofs containing group,
                partial_list = self.isolate_group(moflist, group)
                ranked_list = self.rank_by_uptake(partial_list)
                for mof in ranked_list:
                    if self.options.max_gridpoints:
                        try:
                            ngrid = self.mof_dic[mof]['ngrid']
                        except KeyError:
                            ngrid = self.grid_points(mof)
                            self.mof_dic[mof]['ngrid'] = ngrid
                        ngrid_test = (ngrid > 0 and (ngrid <= 
                                  self.options.max_gridpoints))
                    else:
                        ngrid = True
                    if self.options.uptake_cutoff:
                        # grab the uptake from the original dictionary.
                        ads = self.mof_dic[mof]['mmol/g']
                        uptake = True if (ads >= self.options.uptake_cutoff)\
                                else False
                    else:
                        uptake = True
                    if ngrid_test and uptake and mof not in _used_mofs:
                        groups = self.mof_dic[mof]['functional_groups'].keys()
                        dataset.setdefault(tuple(groups), []).append(mof)
                        _used_mofs.append(mof)

            # if uptake cutoff is requested, then keep all above the
            # cutoff, otherwise respect the functional_max parameter
            if not self.options.uptake_cutoff:
                for key, value in dataset.items():
                    dataset[key] = value[:self.options.functional_max]
        else:
            ranked_list = self.rank_by_uptake(moflist)
            rankcount = 0
            for mof in ranked_list:
                rankcount += 1
                groups = self.mof_dic[mof]['functional_groups'].keys()
                if self.options.max_gridpoints:
                    try:
                        ngrid = self.mof_dic[mof]['ngrid']
                    except KeyError:
                        ngrid = self.grid_points(mof)
                        self.mof_dic[mof]['ngrid'] = ngrid
                    ngrid_test = (ngrid > 0 and (ngrid <= 
                              self.options.max_gridpoints))
                else:
                    ngrid_test = True
                # uptake_cutoff will override the max number of mofs set
                # in the input, this requires two booleans, max and uptake
                # to ensure the final dataset is correct.
                if self.options.uptake_cutoff:
                    # grab the uptake from the original dictionary.
                    ads = self.mof_dic[mof]['mmol/g']
                    uptake = True if (ads >= self.options.uptake_cutoff)\
                             else False
                    max = False
                else:
                    uptake = True
                    max = False if (rankcount <= self.options.total_mofs)\
                            else True 
                if ngrid_test and uptake and not max and mof not in _used_mofs:
                    dataset.setdefault(tuple(groups), []).append(mof)
                    _used_mofs.append(mof)

        self.write_dataset(dataset)

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

    def random_select(self):
        """Select a list of MOFs randomly 

        """
        dataset = {}
        organics_count, fnl_groups_count = {}, {}
        top_count, met_index_count = {}, {}
        mofcount = 0
        # generate a list of valid mofs which obey the inclusive, partial and 
        # exclude lists.
        moflist = self._gen_moflist()
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
            if self.options.gaussian:
                upt_wght = self.weight_by_gaussian(uptake)
            else:
                upt_wght = True
            if self.options.max_gridpoints:
                try:
                    ngrid = self.mof_dic[mof]['ngrid']
                except KeyError:
                    ngrid = self.grid_points(mof)
                    self.mof_dic[mof]['ngrid'] = ngrid
                ngrid_test = (ngrid > 0 and 
                    (ngrid <= self.options.max_gridpoints))
            else:
                ngrid_test = True
            # check to see if only the weighted uptake failed.
            # If so, put the mof back in the pool for selection later..
            # I did this because the lists are not completing.
            if not upt_wght and not org_max and not fnl_max and not top_max \
                    and not met_max and ngrid_test:
                # put the mof back in the random selection pool.
                moflist.append(mof)

            # debug to figure out why full lists are not being made
            if not ngrid_test or org_max or fnl_max or top_max or met_max or \
                    not upt_wght:
                debug("Tests for mof %s"%mof)
                debug("Grid point test: %s"%ngrid_test)
                debug("Weight on uptake test: %s"%upt_wght)
                debug("Max Organic test: %s"%org_max)
                debug("Max Functional group test: %s"%fnl_max)
                debug("Max topology test: %s"%top_max)
                debug("Max metal test: %s"%met_max)
            
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
                mofcount += 1
            if mofcount >= self.options.total_mofs:
                done = True
        self.write_dataset(dataset)

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
        5 mmol/g, smeared by 1.5 (width)"""
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

    def write_dataset(self, dataset):
        """Writes the data to a file."""
        os.chdir(self.options.job_dir)
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
        filename = create_csv_filename(basename) 
        info("Writing dataset to %s.csv ..."%(filename))
        outstream = open(filename+".csv", "w")
        header="MOFname,mmol/g,functional_group1,functional_group2"
        if (self.options.max_gridpoints is not None) or\
                (self.options.report_ngrid):
            header += ",ngrid\n"
        else:
            header += "\n"
        outstream.writelines(header)
        for fnl, mofs in dataset.iteritems():
            if len(fnl) == 1:
                debug("Functional group: %s"%fnl)
                groups = tuple([fnl[0], None])
            elif len(fnl) == 2:
                debug("Functional groups: %s, %s"%(fnl))
                groups = fnl
            for mof in mofs:
                debug("     " + mof)
                try:
                    uptake = self.mof_dic[mof]['mmol/g']
                except KeyError:
                    uptake = 0.
                line = "%s,%f,%s,%s"%(mof, uptake, groups[0], groups[1])
                if (self.options.max_gridpoints is not None) or \
                        (self.options.report_ngrid):
                    try:
                        ngrid = self.mof_dic[mof]['ngrid']
                    except KeyError:
                        ngrid = self.grid_points(mof)
                    line += ",%i\n"%(ngrid)
                else:
                    line += "\n"
                outstream.writelines(line)
        outstream.close()
        info("Done.")

    def _gen_moflist(self):
        """Returns a list of MOFs which are chosen based on the 
        discrimination input lists.

        """
        moflist = []
        inclusive = self.options.fnl_include
        partial = self.options.fnl_partial
        exclude = self.options.fnl_exclude
        if not inclusive and not partial and not exclude:
            return self.mof_dic.keys()

        for key, value in self.mof_dic.iteritems():
            try:
                (group1, group2) = value['functional_groups'].keys()
            except KeyError:
                value['functional_groups'] = {False:[], None:[]}
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
                else:
                    warning("could not find %s-CO2.csv"%(mof))
                    new_uptake = 0.
                os.chdir('..')
            else:
                warning("could not find %s in the directory %s"%(
                        mof, self.basedir))
                new_uptake = 0. 
            self.mofs[mof]['new_uptake'] = new_uptake
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
        header = "MOFname,mmol/g,Functional_grp1,Functional_grp2"
        if self.extended:
            header += ",grp1_replacements,grp2_replacements"
        if self.options.report_ngrid:
            header += ",ngrid"
        header += "\n"
        outstream.writelines(header)
        for mof in self.mofs.keys():
            new_uptake = self.mofs[mof]['new_uptake']
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
            line = "%s,%f,%s,%s"%(mof, new_uptake, fnl_grp1, fnl_grp2)
            if self.extended:
                line += "%s,%s"%(rep_1, rep_2)
            if self.options.report_ngrid:
                line += ",%i"%(grid_points(mof))
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
        self._direct_job()

    def _check_if_ok(self):
        if not self.options.combine and not self.options.cmd_opts.combine and\
                not self.options.dataset and not self.options.report and \
                not self.options.extract and not self.options.report:
            error("no job type requested!")
            sys.exit(1)

    def _direct_job(self):
        self._check_if_ok()
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
        sel.trim_non_existing()
        sel.random_select()

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
        sel.trim_non_existing()
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
                    warning("MOF %s not included in correlation calculations "%mof
                        + "because of 0 value for uptake")
            else:
                warning("MOF %s is not included in correlation calculations "%mof +
                    "because it couldn't be found in one or more csv files.")
        for upt1, upt2 in itertools.combinations(corr, 2):
            ind1 = corr.index(upt1)
            ind2 = corr.index(upt2)
            info("Correlation coefficients for %s, and %s:"%(
                 basenames[ind1], basenames[ind2]))
            rho, pval = stats.spearmanr(upt1, upt2)
            info("Spearman rank: %7.5f"%rho)
            pears, p2 = stats.pearsonr(upt1, upt2)
            info("Pearson correlation: %7.5f"%pears)

        names = {}
        [names.setdefault(i.split('.')[0], 0) for i in basenames]
        name = names.keys()[0]
        name += ".merged"
        write_csv(name, merged)

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
