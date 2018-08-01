#! /usr/local/python_anaconda/bin/python3.4

import re


def filter_fastq_by_read_id(fastq_path, filtered_fastq_path, read_id_list):
    '''
    Filters fastq file, keeps only reads in the given list of read ids.
    :param fastq_path: input fastq path.
    :param filtered_dastq_path: path to create filtered fastq in.
    :param read_id_list: list of read ids to keep.
    '''
    new_fastq = ''
    with open(fastq_path) as f:
        old_fastq_str = f.read()
    for read_id in read_id_list:
        p = re.compile('(@' + read_id + '\n.*?\n)@M', re.DOTALL)
        read_data = p.findall(old_fastq_str)[0]
        new_fastq += read_data
    with open(filtered_fastq_path, 'w') as f:
        f.write(new_fastq)
    return

def filter_out_fastq_by_read_id(fastq_path, filtered_fastq_path, read_id_list):
    '''
    Filters fastq file, keeps only reads *not* in the given list of read ids.
    :param fastq_path: input fastq path.
    :param filtered_dastq_path: path to create filtered fastq in.
    :param read_id_list: list of read ids to filter out.
    '''
    with open(fastq_path) as f:
        fastq_str = f.read()
    for read_id in read_id_list:
        p = re.compile('(@' + read_id + '\n.*?\n)@M', re.DOTALL)
        read_data = p.findall(fastq_str)[0]
        fastq_str = fastq_str.replace(read_data, '')
    with open(filtered_fastq_path, 'w') as f:
        f.write(fastq_str)
    return