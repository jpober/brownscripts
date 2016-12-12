Readme file for PAPER64 Foreground Subtraction and Analysis Scripts

The following scripts are relevant for the use of Foreground Subtraction being applied by FHD to PAPE64 Observation Data:

FHDPipePrep.sh:
Slurm script used to batch process many observation miriad files. Requires the following necessary scripts
to function.

   -fix_anttable.py: Fixes the antenna labels and appends a 'T' to the miriad file.
   -add_uvws.py: Adds uvw coords. which may have been removed or improperly added before and appends a 'U' to the miriad file.
   -miriad2uvfits.py: Takes in a miriad file and outputs an FHD acceptable uvfits file.
   
FgndSubtract.sh:
Master slurm script which batch processes PSA64 observations. Must be used in conjunction with the paper_psa64.cfg config. script,
which is required by FHD to pass it important parameters. This is where parameters such as subtraction catalog, catalog sources, etc
are given, and not all are required as they can be defaulted in FHD already.

batch_sav2miriad.sh:
Slurm script to batch process output FHD dirty and model visibilities. The output files after all of the internal processing scripts
should give you a dirty, residual, model, and filtered model miriad visibilities. This script relies on the following python scripts.

   -save2miriad.py: Converts all IDL sav files in the FHD output directories to miriad files. Output from this script should just be 
   		    a Dirty and Model miriad file with a 'HP' or 'SP' suffix appended and a 'P' prefix added to denote its been through
		    FHD. The 'P' added means it has been converted from linear 'xx'/'yy' polarizations to Stokes 'I'.
   -reduce_seps.sh (Optional): Reduces the number of baseline separations in the miriad files. This can save an appreciable amount of
   		   	       space if the analysis only requires e.g. 30m baselines.
   -ModelDelayFilter.py: Applies a wideband filter to delays outside of the horizon limit for a specific baseline type. This is necessary
   		       to remove high delay noise introduced by FHD, and to reduce the amount of this noise that can be subtracted into
		       the residual visibility.
   -uv_sub.py: Used to subtract the high delay filtered model visibilities from the dirty visibilities, which outputs the residual 
   	       visibilities.

#########

The following scripts are relevant for the analysis of the foreground subtracted observations.

mk_pspec_ratio.sh: Modified script from the original PSA64 data analysis. This script handles the bootstrapping of the LST Binned and 
		   separated miriad observations. It relies on the following modified bootstrapping scripts.
   -pspec_psa64.cfg: Config. file for parameters related to the PSPEC sampling, channels, and separations.
   -jrk_pspec_cov_boot_v003.py: Bootstrap samples the LST binned observations of both the dirty and residual observations in parallel.
   				The output is N bootstrap samples of the ratio of the power spectrum of Residual/Dirty.
   -jrk_pspec_cov_v002.py: Continued bootstrap sampling over the separations relevant for the power spectrum analysis, output is an npz
   			   file which contains the folded(as of 10Dec2016 this is wrong) and unfolded ratio'd PSPECs per k mode.


PlotRatioPSPEC.py: Takes output npz from the mk_pspec_ratio.sh process, and plots the ratio'd PSPECs. Has the ability to compare over multiple
		   channel ranges (redshifts).

power_check.py: !!!NEEDS TO BE REWRITTEN!!! Checks percentage of power removed from the dirty visibilities compared to the residuals. Intended
		to be used in the analysis over power removal per baseline length. Outputs npz files with necessary ratio of sampled(or summed) data.

PlotPowerCompare.py: !!!NEEDS TO BE REWRITTEN!!! Imports npz files from power_check.py output for analysis.