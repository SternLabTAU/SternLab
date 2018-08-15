**Here is a set of commands that you might want to have in your .rc file:**
------

**module loads**
* To load R: module load R/R302
* To load python 3.4: module load python/anaconda3-4.0.0
* To load python 3.6: module load python/anaconda_python-3.6.1
* To allow multi-threading using CPU: module load rocks-openmpi
* To allow multi-threading using GPU: module load openmpi-x86_64
* Colors control: PS1='\[\033[1;36m\]\u\[\033[1;31m\]@\[\033[1;32m\]\h:\[\033[1;35m\]\w\[\033[1;31m\]\$\[\033[0m\] '

**Cool qstat aliases:**
* To get job description: qstat -f <JobID>
* To get a jog array progress by each array entry: qstat -f <JobID[]>
* To get all computational nodes for all jobs running: for job in `qstat | grep username | cut -f1 -d.`; do echo  $job; qstat -f $job | grep host;  done
* **-a** Display all jobs in any status (running, queued, held)
* **-r** Display all running or suspended jobs
* **-n** Display the execution hosts of the running jobs
* **-i** Display all queued, held or waiting jobs
* **-u** username  Display jobs that belong to the specified user
* **-s**  Display any comment added by the administrator or scheduler. This option is typically used to find clues of why a job has not started running.
* **-f job_id**  Display detailed information about a specific job
* **-xf job_id**,
* **-xu user_id** Display status information for finished jobs (within the past 7 days). This option is only available in the newer version of PBS.
<br/><br/>

**Predefined variables:**
------
```
PBS_HOME          The home dirrectory Where you ran qsub command
$PBS_O_LOGNAME     Your UID
$PBS_O_HOST        Your host or node name from where you ran qsub
$PBS_O_QUEUE       The queue name from which you submitted your job
$PBS_QUEUE         The queue name where your job is running
$PBS_JOBID         The job ID
$PBS_JOBNAME       The job name
$PBS_ENVIRONMENT   The PBS enviromental variables
$PBS_O_WORKDIR     The working directory from where you ran qsub
```
**Here are a few useful qalter commands:**
------

ENTER HERE

<br/><br/>
	
**Our nodes:**
------
**Address:** power8.tau.ac.il

**Our nodes:**
* compute-0-163 
* compute-0-164
* compute-0-165
* compute-0-166
* compute-0-167
* compute-0-168
* compute-0-169
* compute-0-170
* compute-0-171


**Please note: users that at some point during they're studies had a user created for them under computer science directory must request Danny to change their defalut directory in power. The password must be reset as well.**

**Syntax to submit a job:** qsub \<job_file_path\>

**Job syntax: bash. Example:**
```
#!/bin/bash

#PBS -S /bin/bash
#PBS -r y
#PBS -q adis
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N SomeCoolName
#PBS -e /jobs_output/
#PBS -o /jobs_output/
#PBS -l select=ncpus=1:mem=4gb 
<commands>
```
Unfotunately, there is no way in PBS to name the log files of the job and also redirect them to a designated directory. So, currently you have two options:

Either name the log files, but have them created in the directory from which the jobs where submitted. Do this by using:
```
#PBS -e /jobs_output/$JOB_NAME.$JOB_ID.ER
#PBS -o /jobs_output/$JOB_NAME.$JOB_ID.OU
```
Or name the log files by the ID of the job, and have them created in your designated directory. Do this by using:
```
#PBS -e /jobs_output/
#PBS -o /jobs_output/
```
The best solution is to add your username to 'get_cmdfile_dir' under 'pbs_jobs.py'. This way you can create a folder in which all job outputs are
saved by date and job alias. 

**How to enter cluster interactive mode:** Here, you have two options:
Enter by submitting a job to cluster interative mode: 
```
qsub /sternadi/home/volume1/shared/InteraciveScripts/clustInteractive.cmd
```
Where clustInteractive.cmd consists of:
```##############################################################
# Interactive Job Script 0.2
#
# This gives you a command prompt running in the kerenh queue.
# Resources:
# - Bash shell
# - 4 GB RAM
# - 1 CPU core
#####################################
#!/bin/bash
#PBS -S /bin/bash
#PBS -q adis
#PBS -N clustInteractive
#PBS -l select=ncpus=1:mem=4gb
#PBS -I
```
This option is better, due to a number of reasons:
* The session is considered, and thus, monitored by the PBS queue system, and you can limit the number of CPUs and memory that it can use ahead (with the ```-I``` argument).
* Unlike the direct connection to a node, the session is not terminated when the connection is interrupted. This is ahuge advantage, because say you run a script that submits jobs in bulk, and then repeatedly check which are done - you no longer need to keep the connection alive to assure that the script finishes its task. If the connection to the cluster is interrupted (because you went home / went to class or any other reason), you can revive the screen of the session. This is done by simply executing the command ```screen``` as soon as the session begins and before interrupting the connection, execute ```ctrl +a```. Dvory also provided addional information on how to utilize this approach:

Please see screen main useful commands below:
* screen                                          # to create a screen in a machine
* ctl +a, d                                       # to detach a screen – go out from a screen – but make it resumeable
* screen –ls                                      # to see all the screens and their IDs
* screen –r <id>                                  # to attach a screen
* ctl+A ?                                         # shows all options in the screen
* clt+A, x                                        # lock the screen for user
* ctl+a, k                                        # kill the screen
* ctl+a, <tab>                                    # resume to last screen
* screen –L                                       # Make terminal logging : to write to a file all the corresponding in the terminal of that screen
* screen -X -S <sessionid> kill                   # kill the screen outside
Also, If you press the arrow down key stroke, you may tell whether you are on a screen or not (if the display blinks – then you are in a screen)


**To check which modules are available:** modules avail


**Tips from HPC team:**
* Use intel's compilers (icc for c and icpc for c++) instead of gcc and g++ (more effective). The icc and icpc equivalent to gcc > 6.0 are available in /powerapps/share/intel/parallel_studio_xe_2018/bin/ 
* The default shell script language in power8 is bash, so in order to set environment vairables, use ```export <VAR_NAME>=<VAR_VALUE>``` (bash compatible) instead of ```setenv <VAR_NAME> <VAR_VALUE>``` (csh / tcsh compatible)

**Useful PDB commands**
* To delete all your jobs: qselect -u <username> | xargs qdel
* To check which nodes your jobs are going to: qstat -1 -n | grep <username> (if usernae fails, try to trim it)


