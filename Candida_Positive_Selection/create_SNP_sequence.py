import pandas as pd
from optparse import OptionParser


def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-p", "--patient1", dest="PATIENT1", help="Path to patient1's data")
    parser.add_option("-o", "--output", dest="DEST", help="Output path")

    (options, args) = parser.parse_args()

    DEST = options.DEST
    PATIENT1 = options.PATIENT1

    create_SNP_sequence(DEST, PATIENT1)

def create_SNP_sequence(DEST, PATIENT1):
    """
    Creates a new fasta file that contains the SNP sequence for all the samples.
    :param DEST: Output file path
    :param PATIENT1: Patient1's data file
    """
    data = pd.read_csv(PATIENT1)
    data.columns = ["chr", "pos", "ref"] + list(data.columns)[3:]

    sequences = ["" for i in range(1,19)]

    new_letters = ["" for i in range(1,19)]
    for line in data.iterrows():

        new_letters[0] = line[1]["ref"]
        for sample in [x for x in range(1,18) if x != 10]:
            new_letters[sample] = line[1][str(sample)]

        if valid_line(new_letters):
            sequences[0] += new_letters[0]
            for sample in [x for x in range(1, 18) if x != 10]:
                if new_letters[sample] == "-":
                    sequences[sample] += new_letters[0]
                else:
                    sequences[sample] += new_letters[sample].upper()

    # write into new file file
    new_file = open(DEST, "w")
    for sample in (x for x in range(0, 18) if x != 10):
        new_file.write(">"+str(sample)+"\n")
        new_file.write(sequences[sample])
        new_file.write("\n")

    new_file.close()


def valid_line(lst):
    """
    Recives a list of all the letters in a specific SNP and decided if it should be included.
    If all the letters in all the samples are the same - return False
    If one of the letters is not leagal ("N,N" etc) - retrun False
    Else, return True
    """
    identical = True
    for letter in lst[1:]:
        if len(letter) != 1 and letter != "":
            return False
        elif letter != lst[1] and letter != "":
            identical = False
    if identical:
        return False
    else:
        return True


if __name__ == "__main__":
    main()