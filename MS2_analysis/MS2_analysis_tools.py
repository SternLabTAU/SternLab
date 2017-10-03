from FITS.create_fits_input import remove_del_ins
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
primers = list(range(1, 20)) + list(range(1291, 1304)) + list(range(1179, 1200)) + list(range(2270, 2288)) +\
                list(range(2167, 2188)) + list(range(3548, 3570))
problematic = [17, 18, 19, 20, 21, 22, 23, 183, 188, 189, 190, 196, 274, 317, 364, 452, 762, 2719, 3117, 3133, 3139,
                   3143, 3146, 3150, 3401, 3539, 3542]

def main():
    path = r"C:\Users\Owner\Google Drive\Lab\Analysis\ms2_20170420\Freqs\freqs_170430"
    df = pd.read_csv(path, sep='\t')

    df = df[~df['Pos'].isin(primers + problematic)]

    test_selective_swipe(df, 18, 37, 'A')
    test_selective_swipe(df, 18, 37, 'B')
    test_selective_swipe(df, 18, 41, 'A')
    test_selective_swipe(df, 18, 41, 'B')
    test_selective_swipe(df, 20, 37, 'B')
    test_selective_swipe(df, 20, 37, 'A')
    test_selective_swipe(df, 20, 41, 'B')
    test_selective_swipe(df, 20, 41, 'A')

def test_selective_swipe(df, passage, degree, replica):
    df = df[(df['Degree'] == degree) & (df['Time'] == passage) & (df['Replica'] == replica)]
    df = remove_del_ins(df)
    df["Mutation"] = df['Ref'] + df['Base']
    df = df[(df['Mutation'] == 'AG') | (df['Mutation'] == 'GA') | (df['Mutation'] == 'CT') | (df['Mutation'] == 'TC')]
    df = df[df["Freq"] >= 0.0004]

    assert(len(set(df['Pos'].values)) == df.shape[0])

    plt.scatter(x=df['Pos'].values, y=df['Freq'].values, color='green', alpha=0.5)
    plt.title("Minor transition allele frequency as a function of position in the genome\n p{}{}{}".format(passage, degree, replica))
    plt.yscale("log")
    plt.ylim(0.0001, 1)
    plt.show()

if __name__ == "__main__":
    main()