#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename, check_dirname
import os


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-p", "--pipeline_path", dest="pipeline_path", default=None, help="run other pipeline script path")
    parser.add_option("-t", "--type", dest="pipeline_type", default=None, help="cirseq or NGS - C or N")
    parser.add_option("-i", "--input", dest="input", help="input directory with fastq.gz files")
    parser.add_option("-o", "--o", dest="output", help="output directory")
    parser.add_option("-r", "--reference", dest="reference", help="reference genome seq (fasta)")
    parser.add_option("-e", "--e", dest="error", default = None, help="error output file (with path)")
    parser.add_option("s", "--start", dest="start", default=1, type="int", help="start at step number. default: 1")
    parser.add_option("--end", dest="end", default=4, type="int", help="end at step number, default: 4")
    parser.add_option("--type_of_input_file", dest="type_of_input_file", default="f", help="type of input files, optional f if fastq and not zipped files")
    parser.add_option("-g", "--gaps", dest="gaps", default="Y", help="refer to gaps? Y/N default Y")
    parser.add_option("--NGS_or_Cirseq", dest="NGS_or_Cirseq", type="int", help="NGS/Cirseq? type 2 for Cirseq (default;min num repeats=2) or 1 for NGS (min num repeats=1)")
    parser.add_option("-q", "--q_score", dest="q_score", default=23, type="int", help="Q-score cutoff, default =23")
    parser.add_option("-b", "--blast", dest="blast_id", type="int", help="% id for blast, default=85")


    cirseq_pipeline_path = "/sternadi/home/volume1/shared/tools/pipeline/v5.2_power7/pipeline_runner.v5.1.pl"
    NGS_pipeline_path = "/sternadi/home/volume1/shared/tools/pipelineNGS/5.1_power7/pipeline_runner.v5.1.1.pl "

    (options, args) = parser.parse_args()
    # chooses which pipeline path to run
    pipeline_path = options.pipeline_path
    type = options.type
    NGS_or_Cirseq = options.NGS_or_Cirseq
    if pipeline_path == None and type == None:
        raise Exception("Need to specify the type of pipeline wanted or an alternative pipeline running path")
    elif pipeline_path != None and type == None:
        pipeline_path = check_filename(pipeline_path)
        if NGS_or_Cirseq in [1, 2]:
            type = NGS_or_Cirseq
        else:
            raise Exception("Need to specify if it's a cirseq or NGS pipeline")
    elif pipeline_path == None and type == "C":
        pipeline_path = cirseq_pipeline_path
        type = 2
    elif pipeline_path == None and type == "N":
        pipeline_path = NGS_pipeline_path
        type = 1
    else:
        raise Exception("Can't specify also a pipeline path and also a pipeline type - need to specify only one")

    input = options.input
    input = check_dirname(input)
    output = options.output
    output = check_dirname(output)
    reference = options.reference
    reference = check_filename(reference)
    error = options.error
    if error == None:
        error = output + "/error.txt"
    error = check_filename(error, Truefile=False)

    start = options.start
    end = options.end
    if start not in [1,2,3]:
        raise Exception("Not a valid start step - needs to be between 1:3")
    if end not in [2,3,4]:
        raise Exception("Not a valid end step - needs to be between 2:4")

    type_of_input_file = options.type_of_input_file
    gaps = options.gaps
    if gaps not in ["Y", "N"]:
        raise Exception("Not a valid gap - must be Y or N")
    q_score = options.q_score
    blast_id = options.blast_id

    path_to_save_pipeline_command = output + "/pipeline_command.txt"

    if type == 1: #NGS
        cmd = "perl %s %s %s %s %s %i %i %s %s %i %i %i"% (pipeline_path, input, output, reference, error,
                                                           start, end, type_of_input_file, gaps, NGS_or_Cirseq,
                                                           q_score, blast_id)
    elif type == 2: #cirseq
        cmd = "perl %s %s %s %s %s %i %i %s %s %i %i"% (pipeline_path, input, output, reference, error,
                                                           start, end, type_of_input_file, gaps,
                                                           q_score, blast_id)

    os.system(cmd)

    o = open(path_to_save_pipeline_command, "w")
    o.write(cmd)
    o.close()

    print("Ran pipeline")


if __name__ == "__main__":
    main()
