#!/bin/env python

#SBATCH -J paper_pipeline
#SBATCH -t 3:00:00
#SBATCH -n 16
#SBATCH --mem=15G

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
   NGPS=5
   sep=-1
   cal=""
   obsfile=""        ## Load the beam normalization from an obs save file.
   OUT='.'
   in_files=None
   WINDOW='blackman-harris'
   length="14"
   auto_redundancy=1
   min=8              ## Minimum number for a redundant group to be included

   def __init__(self):
       in_files = []    # New file list


def read_config(cfg_file):
   if not os.path.isfile(cfg_file):
        raise ValueError("Invalid config file: " + cfg_file)
   cfg=params()
   for line in open(cfg_file).readlines():
    if re.match(r'^export', line):
        line = " ".join(line.split(' ')[1:]).strip()
        line = line.split("#")[0]       #In case there are comments
        key, val = line.split("=")
        val = val.strip()
        key = key.strip()
        if val.startswith("'") or val.startswith("\""): val = val[1:]
        if val.endswith("'") or val.endswith("\""): val = val[:-1]
        if hasattr(cfg, key):
              setattr(cfg,key,val)
   if hasattr(cfg, 'chans'):
         if not type(cfg.chans) is list: cfg.chans = cfg.chans.split(" ")   #Make into a list
   if hasattr(cfg, 'pols'):
         if not type(cfg.pols) is list:  cfg.pols =  cfg.pols.split(" ")   #Make into a list
   if type(cfg.auto_redundancy) is str: cfg.auto_redundancy = int(cfg.auto_redundancy)
   return cfg



if __name__ == '__main__':
      # If being run directly, take in a list of base file names and a config.

    o = optparse.OptionParser()
    o.set_usage('paper_pipeline.py -c <config file> input_file(s)')  #Miriad files only!
    
    o.add_option('-c', dest='config', help='Config file')
    o.add_option('-s', dest='skiptoplot', action='store_true', help='Skip to plot')
    o.add_option('--append', help='Append string to prefix in config file')
    o.add_option('--reconfig', help='Modify the configuration on the fly with comma-separated settings')
 
    opts,args = o.parse_args(sys.argv[1:])

    if opts.config is None:
       raise AttributeError("Config file needed")
    config = opts.config
    
    if len(args) == 0:
       raise AttributeError("Please provide some files to process")
   
    cfg = read_config(config)
    if not opts.append is None:
       cfg.PREFIX += "_" + opts.append
    
    if not opts.reconfig is None:
       cmds = opts.reconfig.split(',')
       for c in cmds:
           if not '=' in c: continue
           k,v = c.split('=')
           if hasattr(cfg,k): setattr(cfg,k,v)
           print k,v

    if hasattr(cfg, 'chans'):
         if not type(cfg.chans) is list: cfg.chans = cfg.chans.split(" ")   #Make into a list
#    import IPython; IPython.embed()
#    sys.exit()

    #Locate cal file and move nearby, if necessary
#    seek=False
    seek = not os.path.isfile(cfg.cal+".py")
#    try:
#       exec('import %s' % cfg.cal)
#    except ImportError as :
#       seek=True
    
    if seek:
        cmd = "find -L ~/SIM_SETUP -name "+cfg.cal+".py -type 'f'"
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)   #Careful with shell=True!
        finds = p.stdout.readlines()
        if len(finds) == 0:
              print "Cal file "+cfg.cal+" not found. Exiting"
              sys.exit()
        else:
              found = finds[0].strip()
              shutil.copy(found,'.')                                #Copy the cal file locally
    
#    if cfg.auto_redundancy == 1:
#         s = build_redundant('count', cfg, cal=cfg.cal, restore=cfg.cal+"_reds")
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
    
    basedir = os.path.abspath(cfg.OUT)

    
    ### I now have a directory of files filled with files for each separation.
    
    for chan in cfg.chans:
        chandir = basedir+"/"+cfg.PREFIX+"/"+chan
        if not os.path.exists(chandir): os.makedirs(chandir)
        for pol in cfg.pols:
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
                 +" -C "+cfg.cal+" -b "+cfg.NBOOT+" -c " +chan \
                 +" --weight=I --length="+str(cfg.length)+ " --instrument="+cfg.instrument \
                +" --window="+cfg.WINDOW + ' -p '+pol + ' --NGPS='+str(cfg.NGPS) + ' --sep='+str(cfg.sep)
           if not cfg.obsfile == "": cmd += " --obspath="+str(cfg.obsfile)
           cmd += " --output="+outdir+" "+' '.join(cur_files)

#           if opt.batch:
#                 cmd = "sbatch --mem=20G -t 2:00:00 -n 10 "+cmd

           p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)   #Careful with shell=True!
           for c in iter(lambda: p.stdout.readline(), ""):
                 print c.strip()
           stdout, stderr = p.communicate()
           rcode = p.returncode
           if not rcode==0:
                print stderr
                sys.exit()
           bootstrap_files = glob.glob(outdir+"/pspec_boot*.npz")
           cmd = "python ~/capo/pspec_pipeline/pspec_simple_2d_to_1d.py " \
                 +" --output="+outdir+" --nboot="+cfg.NBOOT+" "+' '.join(bootstrap_files)
           p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)   #Careful with shell=True!
           rcode = p.wait()
           if not rcode==0:
                print p.stderr.readlines()
                sys.exit()
           shutil.move(outdir+'/pspec_pk_k3pk.npz', outdir+"/"+cfg.PREFIX+"_"+chan+".npz")


