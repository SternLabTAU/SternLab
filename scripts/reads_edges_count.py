import pandas as pd
import matplotlib.pyplot as plt
import glob




input_file=r'RCseq_no_ss10.part10.fasta.blast'
data_frame = pd.read_csv(input_file, sep='\t', names = ["name", "num1", "num2", "start", "stop","strand", "somthing", "mutations"])
print(data_frame)

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




print(counter_start/len(data_frame['start']))


plt.hist(data_frame["start"], bins=50)
plt.hist(data_frame["stop"], bins=50)

#plt.ylim(0,100)
plt.show()



