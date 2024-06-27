#!/bin/sh
################################################################################
################################################################################
## JOB_LOOP.SH
##        Version 3, 19.08.23
##        Thomas Bissinger
################################################################################
## HOW TO RUN
## To call this script, just go to the folder it is stored in (via cd) and then
## type
## >> ./job_loop.sh MODE MODEL STATE
## in the command line, where MODE is to be replaced by either the option integ
## or create. They stand for two different tasks this script can perform, namely
## create                  creates input (directories and files)
## submit                  submits jobs to cluster, runs MolDyPoP
## matlab                  submits jobs to cluster, runs matlab
## run_local               locally runs jobs. Only recommended for fixing small systems
## submit_individual       submits jobs individually, i.e. running through all runnumbers. To fix single failed jobs.
## 
## The command MODEL decides which model is to be used
## 
## The command STATE decides whether the equilibration or sampling integration
## is performed
################################################################################
################################################################################


################################################################################
## SETTING VARIABLES
## First part of this script just sets parameters, starting with reading the
## input, then setting system variables and simulation parameters
################################################################################
## Checks for correct input
if [[ ! $3 || ( $1 != "create" && $1 != "submit" && $1 != "matlab" && $1 != "run_local" && $1 != "submit_individual" ) ]]; then
  echo 'CALLING SCRIPT JOB_LOOP.SH. Required input
$1:  mode            run mode
     Possibilities:
     "create"                creates inputs
     "submit"                submits jobs to cluster, runs MolDyPoP
     "matlab"                submits jobs to cluster, runs matlab
     "run_local"             locally runs jobs. Only recommended for fixing small systems
     "submit_individual"     submits jobs individually, i.e. running through all runnumbers. To fix single failed jobs.
$2:  state           State in which to simulate
     Possibilities:
     "eq"                    equilibration run is performed
     "integ"                 integration run is performed
$3:  model           model used
     Possibilities:
     "xy"                    XY model
     "mxy"                   mobile XY model [only model for which this is optimized]
     "fmxy"                  frozen mobile xy model (disordered XY model)
'
  exit
fi
mode=$1
state=$2
model=$3

## Variable sim_id gives name to the simulation's output files (will be
## assigned to the variable job_id in the infile.in). May depend on the
## computation mode, but right now it doesn't
sim_id=$state



## Basic directory structure.
## localdir        usually the $(pwd) of the folder where computations are done
## bindir          where specific binaries are stored
## basedir         basic directory for sim results
## submitscript    script to be submitted to the cluster
## globallog       log file contianing information for the submitted jobs
## executable      the executable that runs the simulation. Different for different models
localdir=/data/scc/thobi/00_MolDyPop_Final/Example; # SET TO YOUR SETUP
bindir=$localdir/bin;
submitscript=$bindir/submitscript_runs;
globallog=$localdir/logs/submitlog.log;
executable=$bindir/MolDyPoP

## Parameter values for the sims
## T_VALS           values of T for the loop (noise strength)
## SQRTN_VALS       values for sqrtN for the loop (system size)
## N_T_VALS         number of elements of T_VALS
## N_SQRTN_VALS     number of elements of SQRTN_VALS
#T_VALS=(.11 .17 .18 .185 .19 .195 .20 .24)
SQRTN_VALS=(16 32 64)
N_SQRTN_VALS=${#SQRTN_VALS[@]}
basedir="$localdir/$state"_data"/$model";
eq_basedir="$localdir/eq_data/$model";
if [[ ( $model == "mxy" ) || ( $model == "fmxy" ) ]]; then
  L_VALS=(9.25 18.5 37)
  qmax_vals=( 5.60 2.80 1.40 )
  ANNEAL_STEP_VALS=( 1 2 3)
  TMAX_VALS=(5e2 5e2 5e2)
  if [[ "$state" == "eq" ]]; then
    EQ_TMAX_VALS=(3000 6000 9000)
    T_VALS=.01
    EQ_TPRINTSTEP=.005
    echo ${T_VALS[@]}
  else
    EQ_TMAX_VALS=(500 500 500)
    T_VALS=(.05 .09 .13 .15 .17 .19 .21 .23 .25)
    EQ_TPRINTSTEP="NaN"
  fi
  
elif [[ $model == "xy" ]]; then
  L_VALS=(16 32 64 128 256)
  samp_timestep=(1 2 4 8 10)
  qmax_vals=( 1.50 1.50)
  EQ_TMAX_VALS=(3000 6000)
  TMAX_VALS=( 1e3 2e3)
  ANNEAL_STEP_VALS=( 1 2)
  if [[ "$state" == "eq" ]]; then
    T_VALS=0.01
    EQ_TPRINTSTEP=.05
  else
    T_VALS=(.10 .30 .50 .70 .80 .90 1.00 1.10 1.20 1.30 1.40)
  fi

else
  echo "ERROR   Unknown model $model. Exiting"
  exit
fi
N_T_VALS=${#T_VALS[@]}

## Sets runmin and runmax for each system size. Runs go from runmin to runmax
## runmin_vals      number of first run (for each system size)
## runmax_vals      number of last run (for each system size)
runmin_vals=(1 1 1)
runmax_vals=(100 100 100)

## Number of bins for storing spatial correlations. Can vary with system size.
## rbin_vals        number of bins in r (for example required for SCF functions)
## qbin_vals        number of bins in q (for example required for TCF functions like gxx)
rbin_vals=( 256 256 256)
qbin_vals=( 4 8  16  32 32)

## Chooses values for sqrtN. Choice with 0 and N_SQRTN_VALS gives all of them.
## Individual or range choices possible, too
istart=2
imax=$((N_SQRTN_VALS))

################################################################################
## BASIC LOOP
## First loop (i) goes over all system sizes
## Second loop (j) goes over all T/noise strength
################################################################################
for (( i=istart; i<imax; i++ ))
do
  sqrtN=${SQRTN_VALS[$i]}
  echo $sqrtN
  rbin=${rbin_vals[$i]}
  qbin=${qbin_vals[$i]}
  qmax=${qmax_vals[$i]}
  L=${L_VALS[$i]}
  runmin=${runmin_vals[$i]}
  runmax=${runmax_vals[$i]}
  anneal_step=${ANNEAL_STEP_VALS[$i]}
  for (( j=0; j<N_T_VALS; j++ ))
  do
    T=${T_VALS[$j]}
    workdir="$basedir/sqrtN_$sqrtN/T_$T"

    ############################################################################
    ## MODE CREATE
    ## Creates directory structure and inputs for the simulation to run on the
    ## cluster.
    ############################################################################
    if [[ $mode == "create" ]]; then
      ## Sets folders and files for current working directory (curworkdir)
      curworkdir=$workdir/run_1
      inputdir=$curworkdir/input
      outputdir=$curworkdir/output
      curinfile=$inputdir/infile.in

      ## Creates a directory if non exists and fills it with the necessary
      ## inputs
      if [[ ! -d $curworkdir ]]; then
    	echo "Creating $curworkdir"
    	mkdir -p $curworkdir
    	mkdir $inputdir
    	mkdir $outputdir
    	cp $bindir/basic_infile.in $curinfile
    	cp $bindir/basic_config_integ.in $inputdir/samp_config.in
        cp $bindir/basic_config_eq.in $inputdir/eq_samp_config.in

        ## Changes options in current infile
        infile_change_option.sh "$curinfile" "system" "$model"
        infile_change_option.sh "$curinfile" "mode" "$state"
    	infile_change_option.sh "$curinfile" "kT" "$T"
        infile_change_option.sh "$curinfile" "sqrtN" "$sqrtN"
        infile_change_option.sh "$curinfile" "qmax" "$qmax"
        infile_change_option.sh "$curinfile" "L" "$L"
        infile_change_option.sh "$curinfile" "Tmax" "${TMAX_VALS[$i]}"
        infile_change_option.sh "$curinfile" "eq_Tmax" "${EQ_TMAX_VALS[$i]}"
	if [[ "$state" == "eq" ]]; then
	        infile_change_option.sh "$curinfile" "eq_Tprintstep" "$EQ_TPRINTSTEP"
        	infile_change_option.sh "$curinfile" "eq_anneal_step" "$anneal_step"
	fi
	
      fi
    elif [[ "$mode" == "submit" ]]; then
      ############################################################################
      ## MODE SUBMIT
      ## Submits jobs to cluster
      ############################################################################
      ## The following lines set the options for the scc job. They may depend
      ## on system size. For more dTils, cf. also
      ## http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
      ## -l h_vmem   virtual memory. Larger jobs required more memory
      ## -q          which queue to run on. Larger jobs should run on scc,
      ##             all others can run on whatever nodes are free
      ## -tc         number of jobs to run simultaneously
      ## -t          tasks to be run. Setting this like below means that the
      ##             script will run the tasks from runmin to runmax
  	if [[ ($sqrtN -ge 256) ]]; then
  		qspec="-l h_vmem=15G -q scc"
    	else
    		qspec="-l h_vmem=12G"
    	fi
      	qspec="$qspec -tc 500"
      	qspec="$qspec -t $runmin-$runmax"

      ## Here the actual submission of the job happens. qsub submits the script
      ## denoted by $submitscript to the queue. All commands before $submitscript
      ## (that is the options stored in qspec) are options for the submission.
      ## The commands after $submitscript are options that are passed on to the
      ## submitscript, for dTils on those, look there.
      if [[ "$runmax" -ge "$runmin" ]]; then
	cd $workdir
       	qsub $qspec $submitscript "$mode" "$model" "$sim_id" "$T" "$sqrtN"
        echo "qsub $qspec $submitscript \"$mode\" \"$model\" \"$sim_id\" \"$T\" \"$sqrtN\""
      fi


    elif [[ "$mode" == "matlab" ]]; then
      ############################################################################
      ## MODE MATLAB
      ## Submits jobs to cluster, runs matlab script
      ############################################################################
      ## The following lines set the options for the scc job. They may depend
      ## on system size. For more dTils, cf. also
      ## http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
      ## -l h_vmem   virtual memory. Larger jobs required more memory
      ## -q          which queue to run on. Larger jobs should run on scc,
      ##             all others can run on whatever nodes are free
      ## -tc         number of jobs to run simultaneously
      ## -t          tasks to be run. Setting this like below means that the
      ##             script will run the tasks from runmin to runmax
    	qspec="-l h_vmem=50G"
      qspec="$qspec -l h_rt=00:20:00"
      #qspec="$qspec -tc 125"
      #qspec="$qspec -t $runmin-$runmax"

      ## COMMENTED OUT FOR NOW
      ## Here the actual submission of the job happens. qsub submits the script
      ## denoted by $submitscript to the queue. All commands before $submitscript
      ## (that is the options stored in qspec) are options for the submission.
      ## The commands after $submitscript are options that are passed on to the
      ## submitscript, for details on those, look there.
      # nsuccess=$(ls $workdir/run_*/output/sampling_output_$sim_id.mat | wc -l)
      # echo $workdir $nsuccess
      #if [[ "$nsuccess" != "$runmax" ]]; then
      #	qsub $qspec $submitscript "$mode" "$model" "$sim_id" "$T" "$sqrtN"
      #fi
      #if [[ ! -e $workdir/samp_Dynamics.mat ]]; then
      #  echo $workdir
        qsub $qspec $submitscript "$mode" "$model" "$sim_id" "$T" "$sqrtN" 
      #fi
    elif [[ $mode == "run_local" || $mode == "submit_individual" ]]; then
      for ((runnr=runmin; runnr <= runmax; runnr++)); do
        ## Sets folders and files for current working directory (curworkdir)
        curworkdir=$workdir/run_$runnr
        inputdir=$curworkdir/input
        outputdir=$curworkdir/output
        curinfile=$inputdir/infile.in

	exitvalue=$(grep -F "Exit value:" "$outputdir/data_$sim_id.out" | cut -d ' ' -f 3);
        if [[ "$exitvalue" != "1" ]]; then
          if [[ $mode == "run_local" ]]; then
	    mkdir -p $curworkdir
	    mkdir $outputdir
	    cp -r "$workdir/run_1/input" $curworkdir
            cd $curworkdir
            pwd
	    if [[ "$state" == integ ]]; then
		    eq_outputdir="$eqbasedir/sqrtN_$sqrtN/T_$T/run_$runnr/output"
		    eqfile=$(ls $eq_outputdir/eq_snapshot_T_* | grep -F "$T""00")
	            infile_change_option.sh "$curinfile" "init_file" "$eqfile"
	    fi

            infile_change_option.sh "$curinfile" "randomseed" "$RANDOM"
            $executable
            if [[ -e $curworkdir/submitcheck.log ]]; then
              job_id=$(cat $curworkdir/submitcheck.log | cut -d ' ' -f 9 | cut -d '_' -f 2)
              rm $curworkdir/submitcheck.log
              qdel $job_id.$runnr
            fi
          elif [[ $mode == "submit_individual" ]]; then
            echo $curworkdir
            rm $curworkdir/submitcheck.log
            if [[ ($sqrtN -ge 128) ]]; then
        			qspec="-l h_vmem=15G -q scc"
          	else
          		qspec="-l h_vmem=12G"
          	fi
          	qsub $qspec -t $runnr $submitscript "submit" "$model" "$sim_id" "$T" "$sqrtN" 
            echo "qsub $qspec -t $runnr $submitscript \"submit\" \"$model\" \"$sim_id\" \"$T\" \"$sqrtN\""
          fi
        fi
      done
    fi
  done
done
