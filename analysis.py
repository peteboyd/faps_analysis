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
from optparse import OptionParser
import subprocess

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
global LOOKUPDIR

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
        headings = filestream.readline().strip().lstrip("#").split(",")
        mofind = headings.index(self._mof_name_column)
        try:
            uptind = headings.index(self._uptake_column)
        except ValueError:
            pass
        # read in the first line of the file as the headings
        # currently assumes only one uptake value at a statepoint
        for line in filestream:
            parsed = line.split(",")
            mofname = parsed[mofind].strip()
            if mofname.endswith("-CO2.csv"):
                mofname = mofname[:-8]
            try:
                uptake = parsed[uptind]
            except UnboundLocalError:
                uptake = 0.
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
            if mofname.endswith("-CO2.csv"):
                mofname = mofname[:-8]
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
        for line in filestream:
            line = line.strip()
            line = line.split()
            line = line[0]
            if "sym" in line:
                line.lstrip("#")
                if line.endswith("-CO2.csv"):
                    line = line[:-8]
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
            29,
            28,
            9,
            15,
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

    def random_select(self, exclude=[], inclusive=[], partial=[], gridmax=None):
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
                        if exp(-uptake/4.) < random.random():
                            skip = True

                    if group_count >= group_max:
                        skip = True
                        done = True

                    if not skip:
                        # special case for gridpoints, since their
                        # determination is so expensive
                        if gridmax:
                            # have to source the correct file.
                            dirmof = (LOOKUPDIR + '/' + 
                                      mof[:-4] + '.out.cif')
                            from_cif = CifFile(dirmof)
                            ngrid = GrabGridPoints(from_cif.cell,
                                    from_cif.atom_symbols,
                                    from_cif.cart_coordinates)
                            if ngrid > 0 and ngrid <= gridmax:
                                group_count += 1
                                chosen.append(mof)
                                dataset.setdefault(tuple((group,)), []).\
                                        append(mof)
                                self.mof_dic[mof]['ngrid'] = ngrid

                        else:
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
                    if exp(-uptake/4.) < random.random():
                        skip = True
                for i in groups:
                    if i in exclude:
                        skip = True
                    elif fnl_dic[i] > group_max:
                        skip = True

                if not skip:
                    # special case for gridpoints, since their
                    # determination is so expensive
                    if gridmax:
                        # have to source the correct file.
                        dirmof = (LOOKUPDIR + '/' + 
                                  mof[:-4] + '.out.cif')
                        from_cif = CifFile(dirmof)
                        ngrid = GrabGridPoints(from_cif.cell,
                                from_cif.atom_symbols,
                                from_cif.cart_coordinates)
                        if ngrid > 0 and ngrid <= gridmax:
                            for i in groups:
                                fnl_dic[i] += 1
                            chosen.append(mof)
                            mofcount += 1
                            dataset.setdefault(tuple(groups), []).append(mof)
                            self.mof_dic[mof]['ngrid'] = ngrid
                    else:
                        for i in groups:
                            fnl_dic[i] += 1
                        chosen.append(mof)
                        mofcount += 1
                        dataset.setdefault(tuple(groups), []).append(mof)

                if mofcount >= inclusive_max:
                    done = True

        self.write_dataset(dataset, gridmax)
        # TODO(pboyd): add partial considerations.

    def write_dataset(self, dataset, gridmax):
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
        header="MOFname,mmol/g,functional_group1,functional_group2"
        if gridmax:
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
                uptake = self.mof_dic[mof]['old_uptake']
                line = "%s,%f,%s,%s"%(mof, uptake, groups[0], groups[1])
                if gridmax:
                    try:
                        ngrid = self.mof_dic[mof]['ngrid']
                    except KeyError:
                        ngrid = 0
                    line += ",%i\n"%(ngrid)
                else:
                    line += "\n"
                outstream.writelines(line)
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
        if filename.endswith(".csv"):
            basename = filename[:-4]
        else:
            basename = filename
        done = False
        count = 0
        # make sure there's no overwriting, goes up to 99
        while os.path.isfile(filename):
            count += 1
            filename = basename + "%02d.csv"%(count)
        outstream = open(filename, "w")
        if self.extended:
            header = "MOFname,old_uptake,new_uptake,Functional_grp1," +\
                    "Functional_grp2,grp1_replacements,grp2_replacements\n"
        else:
            header = "MOFname,old_uptake,new_uptake,Functional_grp1," +\
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

class CommandLine(object):
    """Parse command line options and communicate directives to the program."""

    def __init__(self):
        self.commands = {}
        self.command_options()

    def command_options(self):
        usage = "%prog [options]\n" + \
                "%prog -r -M Cu -d Cu_dataset/NO_CHG -q cusql.sqlout" +\
                " -c Cu_dataset_q0.csv"
        parser = OptionParser(usage=usage)
        parser.add_option("-r", "--report", action="store_true",
                          dest="report",
                          help="issue a report based on some calculated data.")
        parser.add_option("-s", "--dataset", action="store_true",
                          dest="dataset",
                          help="generate a dataset of randomly selected MOFs.")
        parser.add_option("-c", "--csvfile", action="store", type="string",
                          dest="csvfilename",
                          help="location of csv file with existing data in it.")
        parser.add_option("-d", "--directory", action="store",
                          dest="dirname",
                          help="location of directory containing .cif files.")
        parser.add_option("-q", "--sqlfile", action="store", type="string",
                          dest="sqlname",
                          help="location of sql file containing mof "+
                               "functionalizations.")
        parser.add_option("-i", "--ignore", action="store", type="string",
                           dest="ignorefile",
                           help="location of list of MOFs to ignore in sampling.")
        parser.add_option("-x", "--excludelist", action="store", type="string",
                          dest="exclude",
                          help="comma (,) delimited list of functional groups"+
                               " to exclude from the sampling")
        parser.add_option("-n", "--inclusivelist", action="store", 
                          type="string",
                          dest="inclusive",
                          help="comma (,) delimited list of functional groups"+
                               " to sample, excluding all others.")
        parser.add_option("-L", "--lookupdir", action="store", type="string",
                          dest="lookup", 
                          default="/shared_scratch/pboyd/OUTCIF/FinalCif",
                          help="lookup directory with all the output cifs.")
        parser.add_option("-M", "--metal", action="store", type="string",
                          dest="metal",
                          help="specify metals to generate dataset. Current" +
                          " options are: Zn, Cu, Co, Cd, Mn, Zr, In, V, Ba or Ni") 

        parser.add_option("-G", "--gridpoints", action="store", type="int",
                          dest="maxgridpts",
                          help="specify the max number of grid points to "+\
                               "allow the structures to be chosen")

        (local_options, local_args) = parser.parse_args()
        self.options = local_options

class GrabGridPoints(int):
    """Extracts an integer of grid points from an egulp calculation."""
    code_loc = "/home/pboyd/bin"
    # directory to submit jobs and remove later.
    dir = "/home/pboyd/.temp"
    geo_file = 'temp.geo'
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
            return self._grab_ngrid_points(comm[0].split('\n')) 

    def _grab_ngrid_points(self, lines):
        """Grab the number of grid points from the output file."""
        search_string = "total number of simulation points"
        for line in reversed(lines):
            if search_string in line:
                parse = line.lstrip(search_string)
                return int(parse.strip())
 
class CifFile(object):
    """This class will grab data from the cif file such as the
    cartesian atom coordinates and the cell vectors.
    
    """

    def __init__(self, moffilename, code_loc="/home/pboyd"):
        """Read in the mof file and parse the data to write a .geo file
        for the egulp method.

        """
        self.mofdata = parse_cif(moffilename)
        self.code_loc = code_loc
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
    if mofname.endswith("-CO2.csv"):
        mofname = mofname[:-8]
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

def create_a_dataset(csvfile=None, sqlfile=None, metal=None,
                     inclusive=None, exclude=None, gridmax=None):
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
    if exclude and inclusive:
        sel.random_select(inclusive=inclusive, exclude=exclude, gridmax=gridmax)
    elif exclude and inclusive is None:
        sel.random_select(exclude=exclude, gridmax=gridmax)
    elif inclusive and exclude is None:
        sel.random_select(inclusive=inclusive, gridmax=gridmax)
    #sel.random_select(exclude=["F", "Cl", "Br", "I", "SO3H", 
    #                                    "NO2", "HCO", "NH2"])
    #sel.random_select(inclusive=["F", "Cl", "Br", "I", "SO3H"])

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
    basename = ".".join(directory.split("/"))
    dir = os.getcwd() 
    data.write_data(filename=dir + "/" + basename + ".report.csv")

def main():
    cmd = CommandLine()
    LOOKUPDIR = cmd.options.lookup
    if cmd.options.dataset:
        inclusive, exclude = None, None
        if cmd.options.exclude:
            exclude = cmd.options.exclude.split(",")
        if cmd.options.inclusive:
            inclusive = cmd.options.inclusive.split(",")
        if not cmd.options.exclude and not cmd.options.inclusive:
            raise Error ("No exclusions or inclusions?")
        if cmd.options.csvfilename:
            test = os.path.isfile(cmd.options.csvfilename)
            if not test:
                print ("ERROR: could not find the .csv file")
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
                         metal=cmd.options.metal,
                         inclusive=inclusive,
                         exclude=exclude
                         gridmax=cmd.options.maxgridpts)
                         
    if cmd.options.report:
        if cmd.options.dirname:
            test = os.path.isdir(cmd.options.dirname)
            if not test:
                print ("ERROR: could not find the directory containing MOFs")
        else:
            print("ERROR: directory not set in the command line")
            sys.exit()
        if cmd.options.csvfilename:
            test = os.path.isfile(cmd.options.csvfilename)
            if not test:
                print ("ERROR: could not find the .csv file")
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

        write_report(directory=cmd.options.dirname,
                     sqlfile=cmd.options.sqlname,
                     csvfile=cmd.options.csvfilename)

if __name__ == "__main__":
    main()
