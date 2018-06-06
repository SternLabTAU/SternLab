
"""
@Author: odedkushnir

"""

import os
import pbs_runners



srr = "SRR349755"
ref = "JN562723.SRR349755.fasta"
folder = "/sternadi/home/volume3/okushnir/SRP/%s/fastq" % (srr)

# for d in folders:
output_dir = "/sternadi/home/volume3/okushnir/SRP/%s/q30_consensus_1e-03" % (srr)
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
cmd = "python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py -i %s -o %s -r /sternadi/home/volume3/" \
      "okushnir/SRP/%s -NGS_or_Cirseq 1 -rep 1  -q 30 -ev 0.001" % (folder, output_dir, ref)
pbs_runners.script_runner(cmd, alias="pipeline_d")