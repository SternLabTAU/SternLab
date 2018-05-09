process goes like this:
split fastq -> run alignment on cluster -> merge outputs -> analyse primer-ids distribution

Linux command for splitting original fastq file(s)
split -l 200000 --numeric-suffixes=1 --suffix-length=3 ${fastq_filename} ${output_prefix_with_dot_in_the_end}

Linux command for rename splitted files to remove extra "0" in the filenames (after CD'ing into the directory)
for file in `ls`; do a=`echo ${file} | cut -f2 -d'.'`; b=$((10#$a)); mv $file ${output_prefix_with_dot_in_the_end}${b}; done

run on cluster barcode_aligner.py with the barcode fasta
concatenate output to a single file
run extract_primer_id.py on the concatenated file
