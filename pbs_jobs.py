#! /usr/local/python_anaconda/bin/python3.4

import os
import subprocess
from time import sleep

def create_pbs_cmd(cmdfile, alias, tnum, gmem, cmds, dir = "", load_python=True):
    with open(cmdfile, 'w') as o:
        o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
        o.write("#PBS -q adis\n")
        o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
        o.write("#PBS -N "+ alias+"\n")
        if (gmem):
            mem=gmem*1000
            o.write("#PBS -l mem="+str(mem)+"mb\n")
       # if jnum != 1:
       #     o.write("#PBS -J 1-"+str(jnum)+"\n\n")
       # #o.write("#PBS -J 3-4 \n")
        if dir != "":
            o.write("ls -land %s\n" % dir)
        o.write("id\n")
        o.write("hostname\n")
        if load_python:
            o.write("module load python/anaconda_python-3.4.0\n")
        o.write("'echo %s\n'" % cmds)
        o.write("\n")
        o.write(cmds)
    o.close()

def submit(cmdfile):
    cmd = "/opt/pbs/bin/qsub " + cmdfile
    result = os.popen(cmd).read()
    return result.split(".")[0]


def check_pbs(job_id):
    """
    :param job_id: The PBS job id
    :return: "Done!", when the job is done
    """
    status = "Running..."
    try:
        process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
        while process != "":
            process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
            sleep(0.05)
    except (subprocess.CalledProcessError):
        process = ""
    if process == "":
        status = "Done!"
    return status
