#!/usr/bin/env python
import logging
from logging import debug, error, info
import os
import sys
import copy
import textwrap
from optparse import OptionParser
import ConfigParser
from StringIO import StringIO

class Options(object):
    """Read in the options from the config file."""
    def __init__(self):
        # read in from the command line first
        self._command_options()
        self._init_logging()
        self.job = ConfigParser.SafeConfigParser()
        self.defaults = ConfigParser.SafeConfigParser()
        self._set_paths()
        self._load_defaults()
        self._load_job()
        self._set_attr()

    def _set_paths(self):
        if __name__ != '__main__':
            self.script_dir = os.path.dirname(__file__)
        else:
            self.script_dir = os.path.abspath(sys.path[0])
        self.job_dir = os.getcwd()

    def _load_defaults(self):
        """Load data from the default.ini in the code path."""
        default_path = os.path.join(self.script_dir, 'defaults.ini')
        try:
            filetemp = open(default_path, 'r')
            default = filetemp.read()
            filetemp.close()
            if not '[defaults]' in default.lower():
                default = '[defaults]\n' + default
            default = StringIO(default)
        except IOError:
            error("Error loading defaults.ini")
            default = StringIO('[defaults]\n')
        self.defaults.readfp(default)

    def _load_job(self):
        """Load data from the local job name."""
        if self.cmd_opts.input_file is not None:
            job_path = os.path.join(self.job_dir, self.cmd_opts.input_file)
            try:
                filetemp = open(job_path, 'r')
                job = filetemp.read()
                filetemp.close()
                if not '[job]' in job.lower():
                    job = '[job]\n' + job
                job = StringIO(job)
            except IOError:
                error("Error loading %s"%(self.cmd_opts.input_file))
                job = StringIO('[job]\n')
        else:
            job = StringIO('[job]\n')
        self.job.readfp(job)

    def _command_options(self):
        """Load data from the command line."""

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
        parser.add_option("-f", "--file", action="store", type="string",
                          dest="input_file", default="None",
                          help="Specify a job-specific input file.")
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
                          dest="weight",
                          help="Make the random selected data a gaussian " + \
                            "weight on the CO2 uptake.")
        parser.add_option("-S", "--silent", action="store_true",
                          help="Print nothing to the console.")
        parser.add_option("-Q", "--quiet", action="store_true",
                          help="Print only warnings and errors.")
        parser.add_option("-V", "--verbose", action="store_true",
                          help="Print everything to the console.")

        (local_options, local_args) = parser.parse_args()
        self.cmd_opts = local_options

    def _init_logging(self):
        """Initiate the logging"""
        if self.cmd_opts.silent:
            stdout_level = logging.CRITICAL
            file_level = logging.INFO
        elif self.cmd_opts.quiet:
            stdout_level = logging.ERROR
            file_level = logging.INFO
        elif self.cmd_opts.verbose:
            stdout_level = logging.DEBUG
            file_level = logging.DEBUG
        else:
            stdout_level = logging.INFO
            file_level = logging.INFO

        logging.basicConfig(level=file_level,
                            format='[%(asctime)s] %(levelname)s %(message)s',
                            datefmt='%Y%m%d %H:%M:%S',
                            filename="log.out",
                            filemode='a')

        logging.addLevelName(10, '--')
        logging.addLevelName(20, '>>')
        logging.addLevelName(30, '**')
        logging.addLevelName(40, '!!')
        logging.addLevelName(50, 'XX')

        console = ColouredConsoleHandler(sys.stdout)
        console.setLevel(stdout_level)
        formatter = logging.Formatter('%(levelname)s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

    def parse_commas(self, option, opt, value, parser):
        setattr(parser.values, option.dest, value.split(','))

    def _set_attr(self):
        """Sets attributes to the base class. default options are over-written
        by job-specific options.

        """
        for key, value in self.defaults.items('defaults'):
            setattr(self, key, value)
        for key, value in self.job.items('job'):
            setattr(self, key, value)

class ColouredConsoleHandler(logging.StreamHandler):
    """Makes colourised and wrapped output for the console."""
    def emit(self, record):
        """Colourise and emit a record."""
        # Need to make a actual copy of the record
        # to prevent altering the message for other loggers
        myrecord = copy.copy(record)
        levelno = myrecord.levelno
        if levelno >= 50:  # CRITICAL / FATAL
            front = '\033[30;41m'  # black/red
            text = '\033[30;41m'  # black/red
        elif levelno >= 40:  # ERROR
            front = '\033[30;41m'  # black/red
            text = '\033[1;31m'  # bright red
        elif levelno >= 30:  # WARNING
            front = '\033[30;43m'  # black/yellow
            text = '\033[1;33m'  # bright yellow
        elif levelno >= 20:  # INFO
            front = '\033[30;42m'  # black/green
            text = '\033[1m'  # bright
        elif levelno >= 10:  # DEBUG
            front = '\033[30;46m'  # black/cyan
            text = '\033[0m'  # normal
        else:  # NOTSET and anything else
            front = '\033[0m'  # normal
            text = '\033[0m'  # normal

        myrecord.levelname = '%s%s\033[0m' % (front, myrecord.levelname)
        myrecord.msg = textwrap.fill(
            myrecord.msg, initial_indent=text, width=76,
            subsequent_indent='\033[0m   %s' % text) + '\033[0m'
        logging.StreamHandler.emit(self, myrecord)
