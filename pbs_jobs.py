#! /usr/local/python_anaconda/bin/python3.4

import os
import subprocess
from time import sleep
import getpass
import datetime

def create_pbs_cmd(cmdfile, alias, gmem=2, cmds="", dir = "", load_python=True, jnum=False):
    with open(cmdfile, 'w') as o:
        o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
        o.write("#PBS -q adis\n")
        o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
        o.write("#PBS -N "+ alias+"\n")
        if alias in cmdfile and datetime.datetime.today().strftime('%Y-%m') in cmdfile:
            o.write("#PBS -o %s\n" % "/".join(cmdfile.split("/")[:-1]))
        if (gmem):
            mem=gmem*1000
            o.write("#PBS -l mem="+str(mem)+"mb\n")
        if jnum:
           if jnum != 1:
               o.write("#PBS -J 1-"+str(jnum)+"\n\n")
       # #o.write("#PBS -J 3-4 \n")
        if dir != "":
            o.write("ls -land %s\n" % dir)
        o.write("id\n")
        o.write("date\n")
        o.write("hostname\n")
        if load_python:
            o.write("module load python/anaconda_python-3.4.0\n")
        
        o.write("\n")
        o.write(cmds)
        o.write("\n")
        o.write("date\n")
    o.close()


def create_array_pbs_cmd(cmdfile, alias, jnum, gmem=7, cmds="", dir="", load_python=True):
    with open(cmdfile, 'w') as o:
        o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
        o.write("#PBS -q adis\n")
        o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
        o.write("#PBS -N " + alias + "\n")
        if alias in cmdfile and datetime.datetime.today().strftime('%Y-%m') in cmdfile:
            o.write("#PBS -o %s\n" % "/".join(cmdfile.split("/")[:-1]))
        if (gmem):
            mem = gmem * 1000
            o.write("#PBS -l mem=" + str(mem) + "mb\n")
        if jnum == 1:
            raise Exception('Jnum parameter should be larger then 1\n')
        else:
            o.write("#PBS -J 1-" + str(jnum) + "\n\n")
            # #o.write("#PBS -J 3-4 \n")
        if dir != "":
            o.write("ls -land %s\n" % dir)
        o.write("id\n")
        o.write("date\n")
        o.write("hostname\n")
        if load_python:
            o.write("module load python/anaconda_python-3.4.0\n")

        o.write("\n")
        o.write(cmds)
        o.write("\n")
        o.write("date\n")
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
        print("")
    except (subprocess.CalledProcessError):
        process = ""
    if process == "":
        status = "Done"
    return status


def get_cmdfile_dir(cmdfile, alias):
    username = getpass.getuser()
    lab_users_dic = {"taliakustin":"/sternadi/home/volume1/taliakustin/temp", "daniellem1":"/sternadi/home/volume1/daniellem1/temp"}
    if username in lab_users_dic.keys():
        tmp_dir = lab_users_dic[username]
        if not os.path.exists(tmp_dir):
            os.system("mkdir %s" % tmp_dir)
        date = datetime.datetime.today().strftime('%Y-%m')
        tmp_dir = tmp_dir + "/%s" % date
        if not os.path.exists(tmp_dir):
            os.system("mkdir %s" % tmp_dir)
        tmp_dir = tmp_dir + "/%s" % alias
        if not os.path.exists(tmp_dir):
            os.system("mkdir %s" % tmp_dir)
        cmdfile = tmp_dir + "/" + cmdfile
    return cmdfile
