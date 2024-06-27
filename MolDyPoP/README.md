\mainpage MolDyPoP Main page

\section intro_sec Introduction

Welcome to MolDyPoP, the simulation environment for Molecular Dynamics for Polar
Particles. You will find a full documentation of the MolDyPoP package on these
pages. Additionally, there is a quick tutorial on how to set yourself up for a
simulation run.

If you truly wish to understand the code, you will have to be able to dig a
little bit. The core MolDyPoP package is written in C++ and this documentation
provides you with an exhaustive explanation of what the individual C++ files
are meant to do, and how the classes and namespaces work.

However, calls to MolDyPoP typically happen via the console and you will need to
submit MolDyPoP calculations to the SCC cluster, which in itself is not part of
the C++ routine and requires knowledge of some shell coding. Some basic scripts
that may be of value to you have been compiled here - to properly understand
them requires an understanding of the bash-shell.

Finally, the results produced in MolDyPoP runs are MATLAB-executable files.
The MATLAB-scripts I used in data analysis are therefore also provided here.

Note that neither the bash nor the MATLAB scripts are fully commented for use,
the commenting effort has been focused on giving a proper explanation of the
workings of MolDyPoP. I will give quick explanations of what the other scripts
do, and with some digging in bash and experimenting with MATLAB syntax you
should be able to understand those scripts fairly quickly.


\section struct_sec Structure

You will find four folders in this implementation.

1. A folder "MolDyPoP". Contains the MolDyPoP C++ code and its documentation.

2. A folder "bin". Contains binary files and shell scripts that may be useful.

3. A folder "Example". This folder contains a fully set up calculation that can
be executed to obtain results (see below). It contains some useful scripts.

4. A folder "MatlabScripts". Contains various MATLAB scritps, some of which will
be explained below. The others may be of interest. A selection of relevant
MATLAB scripts is contained in the Folder Example/Matlab
<table>
  <caption id="multi_row">Folders</caption>
    <tr><th> Folder name   	<th>Content
    <tr><td> MolDyPoP				<td>Contains the MolDyPoP C++ code and its documentation.
    <tr><td> bin  		      <td>Contains binary files and shell scripts that may be useful.
    <tr><td> Example        <td>Contains a fully set up calculation that can be executed
                                to obtain results (see below). It contains some useful scripts.
    <tr><td> MatlabScripts  <td>Contains various MATLAB scritps, some of which will
                                be explained below. The others may be of interest. A
                                selection of relevant MATLAB scripts is contained
                                in the Folder Example/Matlab
</table>



\section Example The Example Code

Here's a quick guide on how to get started with the example code.

1. Go to the folder Example. There are five subfolders:
<table>
  <caption id="multi_row">Folders</caption>
    <tr><th> Folder name   	<th>Content
    <tr><td> bin				    <td>Contains scripts and executables and further data
                                required for running simulations
    <tr><td> Matlab		      <td>Contains Matlab scripts required for simulation
                                and simulation analysis
    <tr><td> logs           <td>Stores log files
    <tr><td> scratchdir     <td>Will contain all scratch files
    <tr><td> eq_data        <td>Will store data for equilibration
    <tr><td> integ_data     <td>Will store data for integration simulation runs
</table>
 All folders except for bin and Matlab should be empty at this point.

2. Go to the bin-Folder. Relevant files here:
<table>
  <caption id="multi_row">Folders</caption>
    <tr><th> Folder name   	<th>Content
    <tr><td> MolDyPoP				<td>MolDyPoP executable
    <tr><td> job_loop.sh		<td>This file does most of the work with the folder
                                architecture.
    <tr><td> submitscript_runs <td>Will be submitted to the scc-Cluster, manages
                                everything that happens on the cluster.
    <tr><td> basic_*        <td>Files that form the basis for calculation inputs.
                                They are used by job_loop.sh and submitscript_runs.
</table>
Both job_loop.sh and submitscript_runs are written in bash.

3. Browse through job_loop.sh to understand how it works. Adjust the variables
for the folder architecture in job_loop.sh. Do the same for submitsript_runs
(the script is a bit more complicated and less commented. Make sure you understand
where files go. You also should carefully read and adjust the first few commented
lines -- a line starting with #$ will be read in SCC. These first lines determine
where the output and error files for the calculations go).

4. Run "./job_loop.sh create eq mxy". This will create a folder structure and input
files in the eq_data folder. It will create files for different system sizes
and different temperature values for the mxy model.

5. Run "./job_loop.sh submit eq mxy". This submits jobs to the scc-cluster. Now,
the scc will calculate equilibrated data. (The maximum system size is
set to N = (64)^2 = 4048 particles -- this may be a bit large, computations will
take rather long. If you want smaller system sizes or faster calculations, you
can play around with the parameters in job_loop.sh). Attention: From experiene,
some of the jobs will not be executed correctly (I claim it's the scc, but I
might be biased). If a job did not terminate successfully, you can either
resubmit with the command "./job_loop.sh submit eq mxy". This will restart calculations
that have not been properly started, will not start calculations in folders where
there are still calculations running - however if a job ends abruptly (without
a possibility of cleaning up), there is no failsafe for that. You can also directly
resubmit jobs by using "./job_loop.sh submit_individual eq mxy". This may however
temper with currently running scripts. If only singular files did not work properly,
you can also consider running them on your local machine by using
"./job_loop.sh run_local eq mxy".

6. Now, after some time you should have fully equilibrated data. Then, you can call
"./job_loop.sh create integ mxy". This will create a folder structure and input
files in the integ_data folder. It will create files for different system sizes
and different temperature values for the mxy model.

7. Finally, a call of "./job_loop.sh submit integ mxy" is in order. This will submit
the sampling runs to the cluster. These generate data. Again, you may need
to use the tips above if not all calculations run properly.

8. You can now call "./job_loop.sh matlab integ mxy" and "./job_loop.sh matlab eq mxy"
in order to run matlab scripts that extract and collect data. These will for example
be stored in files such as "integ_data/sqrtN_16/T_.13/samp_integ.mat".

9. If all data is created correctly, you can now use two simple plot scripts by
going to the folder "Matlab/Plots". There, the magnetization can be plotted
in "Plot_Magnetization.m" and in "Plot_TCF.m"



\section Further Further Introduction

You can now play around with the example code. It is not optimized to run fast,
it does not feature anything but the MXY model.

Also, you can check more of the data sets created.

To understand what you are actually doing, I recommend starting by looking at
job_loop.sh. Then you can check the inner workings of submitscript_runs. You may
also be interested in understanding more about the meaning of specific input values.
Short summaries may be found at the bottom of basic_infile.in, but a better
overview can be found in this MolDyPoP documentation in the documentation of the
class definition of parameters in parameters.h. To understand the output values,
look at the documentation in sampler.h. Finally, to understand more about the
MATLAB code, you may want to look at the MATLAB code created in submitscript_runs
and the files related to it.

Enjoy your time.
