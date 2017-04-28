#!/bin/env python

import aipy as a, numpy as np, capo.omni as co
import sys, optparse, os, warnings, subprocess, shutil
import glob, re

### Define functions for handling redundancies.

warnings.filterwarnings("ignore")


class params:
   """ Holds parameters of the analysis, to be read from a config file  """
   PREFIX='BH'
   instrument='hera'   #Needed for beam normalization
   pols='I'
   chans='95_115'
   NBOOT=100
   cal=""
   OUT='.'
   in_files=None
   EVEN_DATAPATH=OUT+'/even/'
   ODD_DATAPATH=OUT+'/odd/'
   WINDOW='blackman-harris'
   length=14
   auto_redundancy=1
   min=8              ## Minimum number for a redundant group to be included

   def __init__(self):
       in_files = []    # New file list


def read_config(cfg_file):
   if not os.path.isfile(cfg_file):
        raise ValueError("Invalid config file: " + cfg_file)
   opts=params()
   for line in open(cfg_file).readlines():
    if re.match(r'^export', line):
        line = " ".join(line.split(' ')[1:]).strip()
        line = line.split("#")[0]       #In case there are comments
        key, val = line.split("=")
        val = val.strip()
        key = key.strip()
        if val.startswith("'") or val.startswith("\""): val = val[1:]
        if val.endswith("'") or val.endswith("\""): val = val[:-1]
        if hasattr(opts, key):
              setattr(opts,key,val)
   if hasattr(opts, 'chans'):
         if not type(opts.chans) is list: opts.chans = opts.chans.split(" ")   #Make into a list
   if hasattr(opts, 'pols'):
         if not type(opts.pols) is list:  opts.pols =  opts.pols.split(" ")   #Make into a list
   if type(opts.auto_redundancy) is str: opts.auto_redundancy = int(opts.auto_redundancy)
   return opts



if __name__ == '__main__':
      # If being run directly, take in a list of base file names and a config.

    o = optparse.OptionParser()
    o.set_usage('paper_pipeline.py -c <config file> input_file(s)')  #Miriad files only!
    
    o.add_option('-c', dest='config', help='Config file')
    o.add_option('-s', dest='skiptoplot', action='store_true', help='Skip to plot')
#    o.add_option('-b', '--batch', action='store_true', help='Run as a slurm batch job')
    
    opt,args = o.parse_args(sys.argv[1:])

    if opt.config is None:
       raise AttributeError("Config file needed")
    config = opt.config
    
    if len(args) == 0:
       raise AttributeError("Please provide some files to process")
    
    opts = read_config(config)
    
    #Locate cal file and move nearby, if necessary
#    seek=False
    seek = not os.path.isfile(opts.cal+".py")
#    try:
#       exec('import %s' % opts.cal)
#    except ImportError as :
#       seek=True
    
    if seek:
        cmd = "find -L ~/SIM_SETUP -name "+opts.cal+".py -type 'f'"
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)   #Careful with shell=True!
        finds = p.stdout.readlines()
        if len(finds) == 0:
              print "Cal file "+opts.cal+" not found. Exiting"
              sys.exit()
        else:
              found = finds[0].strip()
              shutil.copyfile(found,'.')                                #Copy the cal file locally
    
#    if opts.auto_redundancy == 1:
#         s = build_redundant('count', opts, cal=opts.cal, restore=opts.cal+"_reds")
#         seps = map(str, range(s))
#         print "Seps: ", seps
    
    files = {}
    for f in args:
        files[f] = f
    
    ## Make stokes?
#    stokefiles = [ r+"A" for r in files.values() ]
    stokefiles = []
    for k,f in files.iteritems():
       fil = f+"P"
       files[k] = fil
       if (not os.path.exists(fil)) and (fil not in stokefiles):
           stokefiles.append(f)
    #print stokefiles
    
    if len(stokefiles)>0:
        print "Make stokes"
        cmd = 'python ~/capo/ael/mk_stokes.py --stokes="I" '+' '.join(stokefiles)
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)   #Careful with shell=True!
        rcode = p.wait()
        if not rcode==0:
            print p.stderr.readlines()
            sys.exit()
    
    ## Make separations directories
    
    basedir = os.path.abspath(opts.OUT)

    
    ### I now have a directory of files filled with files for each separation.
    
    for chan in opts.chans:
        chandir = basedir+"/"+opts.PREFIX+"/"+chan
        if not os.path.exists(chandir): os.makedirs(chandir)
        for pol in opts.pols:
           print "Starting work on pol "+pol
           poldir = chandir+"/"+pol
           if not os.path.exists(poldir): os.makedirs(poldir)
           outdir = poldir #+"/sep_"+s    #Cut the sepdirs option
#           if not os.path.exists(outdir): os.makedirs(outdir)
           prev_dir = os.getcwd()
#           os.chdir(sepdirs[s])
           cur_files = files.values()
#           cur_files = glob.glob("*P")
           cmd = "~/capo/pspec_pipeline/pspec_oqe_simple_hera.py " \
                 +" -C "+opts.cal+" -b "+opts.NBOOT+" -a cross -c " +chan \
                 +" --weight=I --length="+opts.length+ " --instrument="+opts.instrument \
                 +" --window="+opts.WINDOW + ' -p '+pol \
                 +" --output="+outdir+" "+' '.join(cur_files)

#           if opt.batch:
#                 cmd = "sbatch --mem=20G -t 2:00:00 -n 10 "+cmd

           print cmd
           p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)   #Careful with shell=True!
#           if opt.batch:
#                rcode = p.wait()
#                lines = " ".join(p.stdout.readlines()).strip()
#                pid   = lines.split(" ")[-1]
#                cmd = '~/extra_scripts/wait_for_slurm.sh '+pid
#                p0 = subprocess.Popen(cmd, shell=True)
#                rcode = p0.wait()
#           else:
           rcode = p.wait()
           if not rcode==0:
                print "\n".join(p.stderr.readlines())
                sys.exit()
           bootstrap_files = glob.glob(outdir+"/pspec_boot*.npz")
           cmd = "python ~/capo/pspec_pipeline/pspec_simple_2d_to_1d.py " \
                 +" --output="+outdir+" --nboot="+opts.NBOOT+" "+' '.join(bootstrap_files)
           print cmd
           p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)   #Careful with shell=True!
           rcode = p.wait()
           if not rcode==0:
                print p.stderr.readlines()
                sys.exit()
           shutil.move(outdir+'/pspec_pk_k3pk.npz', outdir+"/"+opts.PREFIX+"_"+chan+".npz")


