#!/usr/bin/env python

import sys
import os
from shutil import copyfile
import subprocess
from optparse import OptionParser
import stat

#pwd = "/work/nrap682/pboyd"
pwd = os.getcwd()
WQEQ_DIR = "/home/pboyd/lib/WilmerQEq"
faps_dir = "/home/pboyd/codes_in_development/faps"
faps_tools = "/home/pboyd/codes_in_development/faps/tools"

just_repeat = """dft_code = vasp
dispersion = False
charge_method = repeat
no_gcmc = True
no_properties = True
optim_cell = False
optim_h = False
optim_all = False
quiet = True
"""

faps_vasp_opt = """dft_code = vasp
dispersion = False
charge_method = repeat
guests = CO2
mc_eq_steps = 2000000
mc_prod_steps = 10000000
mc_numguests_freq = 50000
mc_probability_plot = False
mc_temperature = 298.0
mc_pressure = 0.15
no_properties = False
optim_cell = True 
optim_all = True 
quiet = True
"""

faps_job_uff_qeq = """charge_method = gulp
no_dft = True
guests = CO2
mc_eq_steps = 2000000
mc_prod_steps = 10000000
mc_numguests_freq = 50000
mc_probability_plot = False
mc_temperature = 298.0
mc_pressure = 0.15
no_dft = True
optim_all = False
optim_cell = False
optim_h = False
quiet = True
qeq_parameters = 
    Cu    4.20000000    4.22000000
"""

faps_job_vasp_gcmc_sp="""charge_method = repeat
dft_code = vasp
find_maxima = False 
guests = CO2
mc_eq_steps = 2000000
mc_prod_steps = 10000000
mc_numguests_freq = 50000
mc_probability_plot = False 
mc_temperature = 298.0
mc_pressure = 0.15
optim_all = False
optim_cell = False 
optim_h = False
quiet = True
"""

faps_job_noq_gcmc="""find_maxima = False 
guests = CO2
mc_eq_steps = 2000000
mc_prod_steps = 10000000
mc_numguests_freq = 50000
mc_probability_plot = False 
mc_temperature = 298.0
mc_pressure = 0.15
no_charges = True
no_dft = True
optim_all = False
optim_cell = False 
optim_h = False
quiet = True
"""

faps_job_egulp_param3="""find_maxima = False 
guests = CO2
mc_eq_steps = 2000000
mc_prod_steps = 10000000
mc_numguests_freq = 50000
mc_probability_plot = False 
mc_temperature = 298.0
mc_pressure = 0.15
charge_method = egulp
qeq_parameters = 
    Cu    5.28660000    3.28889297
    Zn    3.47871363    4.14128701
     C    5.448105      5.53192500
     N    6.30432       6.03174000
     O   10.70614       8.57850000
     F    6.364800     11.136000
    Cl    5.824000      7.272000
    Br    5.760000      8.760000
     I    5.426667      5.720000
     S    3.080533      4.336800
no_dft = True
optim_all = False
optim_cell = False 
optim_h = False
quiet = True
"""

faps_job_egulp_param4d="""find_maxima = False 
guests = CO2
mc_eq_steps = 2000000
mc_prod_steps = 10000000
mc_numguests_freq = 50000
mc_probability_plot = False 
mc_temperature = 298.0
mc_pressure = 0.15
charge_method = egulp
qeq_parameters = 
    Cu    5.42900000    3.46800000
    Zn    3.70100000    4.46300000
     C    5.431000      5.85700000
     N    6.68800       6.62200000
     O    8.71400       8.56800000
     F    6.416000     11.131000
    Cl    5.821000      7.273000
    Br    5.692000      8.760000
     I    5.431000      5.720000
     S    3.369000      5.092000
no_dft = True
optim_all = False
optim_cell = False 
optim_h = False
quiet = True
"""

faps_job_egulp_param4dt2="""find_maxima = False 
guests = CO2
mc_eq_steps = 2000000
mc_prod_steps = 10000000
mc_numguests_freq = 50000
mc_probability_plot = False 
mc_temperature = 298.0
mc_pressure = 0.15
charge_method = egulp
qeq_parameters = 
    Cu    5.429000      3.468000
    Zn    3.701000      4.463000
     C    5.431000      5.857000
     N    6.688000      6.622000
     O    8.714000      8.568000
     F    6.416000     11.131000
    Cl    5.821000      7.273000
    Br    5.692000      8.760000
     I    5.431000      5.720000
     S    3.369000      5.092000
    800   8.714000      8.568000
    801  10.597000      9.744000
    802   7.968000     10.323000
no_dft = True
optim_all = False
optim_cell = False 
optim_h = False
egulp_typed_atoms = True
quiet = True
"""

faps_job_egulp_param4dt3="""find_maxima = False 
guests = CO2
mc_eq_steps = 2000000
mc_prod_steps = 10000000
mc_numguests_freq = 50000
mc_probability_plot = False 
mc_temperature = 298.0
mc_pressure = 0.15
charge_method = egulp
qeq_parameters = 
    Cu    5.429000      3.468000
    Zn    3.701000      4.463000
     C    5.431000      5.857000
     N    6.688000      6.622000
     O    8.714000      8.568000
     F    6.416000     11.131000
    Cl    5.821000      7.273000
    Br    5.692000      8.760000
     I    5.431000      5.720000
     S    3.369000      5.092000
    800   8.714000      8.568000
    801  10.528000      9.543000
    802   8.086000     10.187000
   1001   4.035000      6.722000
no_dft = True
optim_all = False
optim_cell = False 
optim_h = False
egulp_typed_atoms = True
quiet = True
"""

class CommandLine(object):
    """Parse command line options and communicate directives to the program."""

    def __init__(self):
        self.commands = {}
        self.command_options()

    def command_options(self):
        usage = "%prog -f [filename] [options]"
        if len(sys.argv) == 1:
            print usage
            sys.exit(0)
        parser = OptionParser(usage=usage)
        parser.add_option("-Z", "--zerochg", action="store_true",
                          dest="CHARGE_ZERO",
                          help="Create a series of submission directories " +
                          "containing mofs with zero charges.")
        parser.add_option("-W", "--wilmer", action="store_true",
                          dest="CHARGE_WILMER",
                          help="Create submission directories where the MOFs "+
                          "have partial atomic charges assigned from Wilmer's " +
                          "QEq method.")
        parser.add_option("--egulp3", action="store_true",
                          dest="CHARGE_EGULP_3",
                          help="Create submission directories where the MOFS "+
                               "have partial atomic charges assigned from "+
                               "EGULP parameter set 3.")
        parser.add_option("--egulp4d", action="store_true",
                          dest="CHARGE_EGULP_4d",
                          help="Create submission directories with EGULP "+
                               "parameter set 4d.")
        parser.add_option("--egulp4dt2", action="store_true",
                          dest="CHARGE_EGULP_4dt2",
                          help="Create submission directories with EGULP "+
                               "parameter set 4d.t2.")
        parser.add_option("--egulp4dt3", action="store_true",
                          dest="CHARGE_EGULP_4dt3",
                          help="Create submission directories with EGULP "+
                               "parameter set 4dt3.")
        parser.add_option("--vasp_opt", action="store_true",
                          dest="VASP_OPT",
                          help="Create submission directories to optimize " +
                               "the MOFs at the DFT level.")
        parser.add_option("-U", "--uff", action="store_true", 
                          dest="CHARGE_UFF",
                          help="Create submission directories where the MOFs "+
                               "have PAtC's assigned from the UFF "+
                               "parameterization.")
        parser.add_option("-f", "--file", action="store", type="string",
                          dest="file_to_read", 
                          help=".csv file containing all the MOF names for submission.")
        parser.add_option("-d", "--dir", action="store", type="string",
                          dest="lookup_dir", 
                          default="/shared_scratch/tdaff/GROIN_20130307/FinalCif",
                          help="Directory containing all the .cif files to "+
                          "copy over to the submission directory")
        (local_options, local_args) = parser.parse_args()
        self.options = local_options

    def parse_commas(self, option, opt, value, parser):
        setattr(parser.values, option.dest, value.split(','))

def faps_sub():
    line = \
"""#!/bin/bash
faps_dir="/home/pboyd/codes_in_development/faps"

for dir in `ls -d -- */`; do
    mof=${dir%"/"}
    echo $mof
    cd $mof
    $faps_dir/faps.py -s ${mof}
    cd ..
done"""
    file = open('faps_sub.sh', 'w')
    file.writelines(line)
    file.close()
    os.chmod('faps_sub.sh', stat.S_IRWXU)

def zipper_chg_set(basename):
    line = \
"""#!/bin/bash
zipname="%s.zip"
read -a dirs <<< `ls -d -- */`
for dir in ${dirs[@]}; do
    mof=${dir%%'/'}
    if [ -d ${mof}/faps_${mof}_repeat ]; then
        cp ${mof}/faps_${mof}_repeat/faps-${mof}.out ${mof}/${mof}_REPEAT.out
        cp ${mof}/faps_${mof}_vasp/CONTCAR ${mof}/${mof}_CONTCAR
        cp ${mof}/faps_${mof}_vasp/OUTCAR ${mof}/${mof}_OUTCAR
        zip $zipname ${mof}/${mof}.cif
        zip $zipname ${mof}/${mof}_REPEAT.out
        zip $zipname ${mof}/${mof}_CONTCAR
        zip $zipname ${mof}/${mof}_OUTCAR
        rm ${mof}/${mof}_REPEAT.out
        rm ${mof}/${mof}_CONTCAR
        rm ${mof}/${mof}_OUTCAR
    fi
done
"""%(basename)
    file = open('zipper_chg.sh', 'w')
    file.writelines(line)
    file.close()
    os.chmod('zipper_chg.sh', stat.S_IRWXU)

def zipper_fap_cif(basename):
    line = \
"""#!/bin/bash
zipname="%s.zip"
read -a dirs <<< `ls -d -- */`
zip $zipname faps_sub.sh
zip $zipname zipper_chg.sh
for dir in ${dirs[@]}; do
    mof=${dir%%'/'}
    zip $zipname ${mof}/${mof}.cif
    zip $zipname ${mof}/${mof}.fap
done
"""%(basename)
    file = open('zipper.sh', 'w')
    file.writelines(line)
    file.close()
    os.chmod('zipper.sh', stat.S_IRWXU)

def reset_charges_to_zero(cif_file):
    """Reset all the charges reported in a faps-cif file to zero."""
    cifstream = open(cif_file, "r")
    site_loop = "_atom_site_label"
    charge_loop = "_atom_type_partial_charge"
    charge_loop2 = "_atom_type_parital_charge"
    check, scan, read = False, False, False
    lines = ""
    count = 0
    for line in cifstream:
        line = line.strip()

        if check:
            if site_loop in line:
                lines += line + "\n"
                scan = True
                check = False
            else:
                lines += line + "\n"
                check = False

        elif "loop_" in line:
            lines += line + "\n"
            check = True

        elif scan:
            if "_" == line[0]:
                lines += line + "\n"
                count += 1
                if charge_loop in line or charge_loop2 in line:
                    charge_ind = count
            else:
                scan = False
                read = True
        elif not read:
            lines += line + "\n"

        if read:

            try:
                charge_ind
                newl = line.split()
                try:
                    charge = float(newl[charge_ind])
                    #print "Replacing charge: %9.4f, with %3.2f"%(charge, 0.)
                    for id, entry in enumerate(newl):
                        try:
                            newl[id] = float(entry)
                        except ValueError:
                            pass
                    newl[charge_ind] = 0.
                    newl = tuple(newl)
                    newline = "%-5s%-6s%-6s%8.6f %8.6f %8.6f %3.1f\n"%(newl)
                    lines += newline

                except IndexError:
                    lines += line + "\n"
                    read = False
            except NameError:
                print "ERROR - could not find the loop for %s"%(charge_loop)
                sys.exit()

    cifstream.close()
    return lines
def charge_zero_cif(structname):
    # for zero charge assignment
    newmof = reset_charges_to_zero(structname + ".cif")
    os.remove(structname + ".cif")
    newone = open(structname + ".cif", "w")
    newone.writelines(newmof)
    newone.close()

def wilmer_cif(structname, submitdir):
    os.chdir(submitdir)
    copyfile(WQEQ_DIR + "/ionizationdata.dat", 
             submitdir + "/ionizationdata.dat")
    copyfile(WQEQ_DIR + "/chargecenters.dat", 
             submitdir + "/chargecenters.dat")
    job = subprocess.Popen(["%s/WQEq.x"%WQEQ_DIR, "%s.cif"%structname],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    # check for error in job submission
    job_comm = job.communicate()
    if job_comm[1]:
        print("ERROR in wilmer job submission! you should check "+
                "the cif file for %s"%structname)
        sys.exit(1)
    try:
        os.rename("%s.Wilmer_EQeq.cif"%structname, 
                  "%s.cif"%structname)
    except OSError:
        print("could not find the file %s.Wilmer_EQeq.cif"%structname)

    os.remove(submitdir + "/ionizationdata.dat")
    os.remove(submitdir + "/chargecenters.dat")
    return 

def clean(x):
    if x.endswith("-CO2.csv"):
        x = x[:-8]
    if x.endswith(".csv"):
        x = x[:-4]
    return x

def gen_submit_dir(cmd, local_dir, basefile):
    if cmd.options.CHARGE_ZERO:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "NO_CHG"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "NO_CHG")
    elif cmd.options.CHARGE_EGULP_3:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "EGULP.param.3"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "EGULP.param.3")
    elif cmd.options.CHARGE_EGULP_4d:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "EGULP.param.4d"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "EGULP.param.4d")
    elif cmd.options.CHARGE_EGULP_4dt2:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "EGULP.param.4d.t2"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "EGULP.param.4d.t2")
    elif cmd.options.CHARGE_EGULP_4dt3:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "EGULP.param.4d.t3"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "EGULP.param.4d.t3")
    elif cmd.options.CHARGE_WILMER:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "WILMER_QEQ"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "WILMER_QEQ")
    elif cmd.options.CHARGE_UFF:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "UFF_QEQ"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "UFF_QEQ")
    elif cmd.options.VASP_OPT:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "VASP_OPT"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "VASP_OPT")

    else:
        try:
            os.makedirs(os.path.join(local_dir, basefile, "REPEAT"))
        except OSError:
            print("Directory already exists!")
        submit_dir = os.path.join(local_dir, basefile, "REPEAT") 
    return submit_dir

def main():
    cmd = CommandLine()
    filestream = open(cmd.options.file_to_read, "r")
    basefile = cmd.options.file_to_read.split("/")[-1]
    basefile = clean(basefile)
    submit_dir = gen_submit_dir(cmd, pwd, basefile)

    print("Creating submission directory: %s"%submit_dir)
    os.chdir("%s"%(submit_dir))
    # create the submission script
    faps_sub()
    # create the zipper script
    if not any([cmd.options.CHARGE_EGULP_3, cmd.options.CHARGE_WILMER,
        cmd.options.CHARGE_ZERO, cmd.options.CHARGE_UFF]):
        zipname = os.path.join(submit_dir, '.'.join(submit_dir.split('/')[-2:]))
        zipper_fap_cif(zipname)
        zipper_chg_set(zipname)
    else:
        zipname = os.path.join(submit_dir, '.'.join(submit_dir.split('/')[-2:]))
        zipper_fap_cif(zipname)

    for line in filestream:
        os.chdir("%s"%(submit_dir))
        line = line.split(",")[0]
        line = clean(line)
        if ".sym." in line:
            #parse the line to get the base stuff
            stuff = line.split("_")
            metind = int(stuff[1].lstrip("m"))
            orgind1 = int(stuff[2].lstrip("o"))
            orgind2 = int(stuff[3].lstrip("o"))
            otherstuff = stuff[5].split(".")
            topology = otherstuff[0]
            fnl = int(otherstuff[2])
            # now get the mof from the lookup directory 
            base_name = "str_m%i_o%i_o%i_f0_%s"%(metind,
                                                 orgind1,
                                                 orgind2,
                                                 topology)
            dir_name = "met%i%s"%(metind, topology)
            structname = base_name + ".sym.%i"%(fnl)
            # create the directory where the new job will be run
            os.makedirs("%s/%s"%(submit_dir, structname))
            # move the new cif file to the new directory
            copyfile("%s/%s.out.cif"%(cmd.options.lookup_dir,structname),
                     "%s/%s/%s.cif"%(submit_dir, structname, structname))
            # change to the new directory
            os.chdir("%s/%s"%(submit_dir, structname))
            if cmd.options.CHARGE_ZERO:
                charge_zero_cif(structname)
                faplines = faps_job_noq_gcmc
            elif cmd.options.CHARGE_EGULP_3:
                faplines = faps_job_egulp_param3
            elif cmd.options.CHARGE_EGULP_4d:
                faplines = faps_job_egulp_param4d
            elif cmd.options.CHARGE_EGULP_4dt2:
                faplines = faps_job_egulp_param4dt2
            elif cmd.options.CHARGE_EGULP_4dt3:
                faplines = faps_job_egulp_param4dt3
            elif cmd.options.CHARGE_WILMER:
                faplines = faps_job_noq_gcmc
                wilmer_cif(structname, "%s/%s"%(submit_dir,structname))
            elif cmd.options.CHARGE_UFF:
                faplines = faps_job_uff_qeq
            elif cmd.options.VASP_OPT:
                faplines = faps_vasp_opt
            else:
                faplines = faps_job_vasp_gcmc_sp
            # create the job info file
            fapfile = open("%s/%s/%s.fap"%(submit_dir, structname, structname), "w")
            fapfile.writelines(faplines)
            fapfile.close()
    filestream.close()

if __name__ == "__main__":
    main()


