import pandas as pd
import matplotlib.pyplot as plt
import glob
from optparse import OptionParser

def main():
    # for Cluster
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--tmp_dir", dest="tmp dir", help="tmp dir")
    parser.add_option("-o", "--output_dir", dest="output dir", help="output dir if the plots")
    (options, args) = parser.parse_args()

    tmp_dir = options.tmp_dir
    out_dir = options.output_dir

    reads_edges_count(tmp_dir, out_dir)


def reads_edges_count(tmp_dir, out_dir):
    input_files=glob.glob(tmp_dir +"/*.blast")
    counter_list=[]
    for file in input_files:
        data_frame = pd.read_csv(file, sep='\t', names = ["name", "num1", "num2", "start", "stop","strand", "somthing", "mutations"])

        counter_start=0
        for i in data_frame['start']:
            if i >450 and i<490:
                counter_start+=1
            if i< 50 and i>3 :
                counter_start+=1
            else:
                continue
        for i in data_frame['stop']:
            if i >450 and i<490:
                counter_start+=1
            if i< 50 and i>3 :
                counter_start+=1
            else:
                continue
        counter_list.append(counter_start/len(data_frame['start']))
        plt.hist(data_frame["start"], bins=50)
        plt.hist(data_frame["stop"], bins=50)
        plt.savefig(out_dir+'/plot.png')
    sum_lst=sum(counter_list)
    print(sum_lst/len(counter_list))


if __name__ == "__main__":
    main()





