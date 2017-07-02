#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
import glob
from file_utilities import check_filename, check_dirname
from freqs_utilities import merge_freqs_files, filter_freqs_for_regression_analysis, add_mutation_to_freq_file
from pandas_utilities import merge_dfs



def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="freqs_dir", help="dir with freqs file")


    (options, args) = parser.parse_args()
    freqs_dir = options.freqs_dir
    freqs_dir = check_dirname(freqs_dir)

    freqs_files = glob.glob(freqs_dir + "/*.freqs")

    with_mutation_files = glob.glob(freqs_dir + "/*_with_mutations.csv")
    if with_mutation_files == []:
        print("Adding mutation to freqs files")
        for freqs_file in freqs_files:
            print(freqs_file)
            output = freqs_file.split(".freqs")[0] + "_with_mutations.csv"
            add_mutation_to_freq_file(output, freqs_file = freqs_file)
            with_mutation_files.append(output)


    segment_files = glob.glob(freqs_dir + "/Segment_[0-9].csv") + glob.glob(freqs_dir + "/Segment_[0-9][0-9].csv")
    if segment_files == []:
        print("Merging segment from different passages")
        for s in range(1,11):
            specific_segment_mutation_files = glob.glob(freqs_dir + "/P*-S%s_with_mutations.csv" % s)
            segment_file, segment_freqs = merge_freqs_files(specific_segment_mutation_files, freqs_dir + "/Segment_%i.csv" % s)
            segment_files.append(segment_file)


    filtered_files = glob.glob(freqs_dir + "/Segment_[0-9]_filtered.csv") + \
                    glob.glob(freqs_dir + "/Segment_[0-9][0-9]_filtered.csv")

    print("Filtering positions from segment csvs")
    for segment_file in segment_files:
        output = segment_file.split(".csv")[0] + "_filtered.csv"
        filtered_file, filtered_ferqs = filter_freqs_for_regression_analysis(output, freqs_file=segment_file)
        filtered_files.append(filtered_file)


    merge_dfs(filtered_files, freqs_dir + "/All_segments_filtered_for_regression.csv")





if __name__ == "__main__":
    main()