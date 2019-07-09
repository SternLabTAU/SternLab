import pandas as pd
import glob
import sys
import numpy as np
from Bio import SeqIO

# Starting assumptions:
# iterate on gene
        # produce 3 reading frames (generate aa sequence on wide interval)- according to primary variant
        # align (no indels) 3 primary RFs to annotation\ref to determine RF + starting point + endpoint
            # choose RF by minimal %AA variance
            # mark AA mutations- AtoB + position

# Simple version:
#   - Only Major muts\variants
#   - No indels


def get_nuc_seq_from_fasta(ref_fasta, gene_start_index, gene_nuc_length):
    return list(SeqIO.parse(ref_fasta, "fasta"))[0]


def protein_mutations_loop():
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files_ET86_2/**/*.freqs')
    mutations = pd.DataFrame(columns=['sample_id', 'pos_major'])

    pol_start_index = 1456
    pol_nuc_gene_length = 2500

    # iterate on genes
    # TODO- resolve comparison method- can Seq be transformed to dataframe?
    ref_gene_nuc_seq = get_nuc_seq_from_fasta(ref_fasta,pol_start_index, pol_nuc_gene_length)
    ref_gene_aa_seq = produce_aa_seq(ref_gene_nuc_seq)

    # pol sanity check
    pol_aa_comparison = "FFRETLAFQQGKAREFPSEQTRANSPTRESQTRANSPTTRELQVRGSNTFSEAGAERQGSLNFPQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEINLPGKWKPKMIGGIGGFIKVRQYDQIIIEICGKKAIGTVLVGPTPVNIIGRNMLTQLGRTLNFPISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALTAICEEMEQEGKISRIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEGFRKYTAFTIPSTNNETPGIRYQYNVLPQGWKGSPPIFQSSMPQILEPFRAPNPEIVIYQYMDDLYVGSDLEIGQHRAPIEELREHLLKWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIQLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLCKLLRGAKALTDIVTLTEEAELELAENREILKEPVHGVFYDPSKDLIAEIQKQGNDQWTFQFYQEPFKNLKTGKFAKRGTAHTNDVKQLTAVVQKIALESIVIWGKTPKFRLPIQKETWEAWWTDYWQATWIPEWEFVNTPPLVKLWYQLEKEPIAGVETFYVDGAANRETKIGKAGYVTDRGRQKIVSLTETTNQKTELQAIQLALQDSGSEVNIVTDSQYALGIILAQPDKSESEIVNQIIEQLISKERVYLSWVPAHKGIGGNEQVDKLVSSGIRKVLFLDGIDKAQEEHEKYHSNWRAMANEFNIPPVVPKEIVACCDKCQLKGEAIHGQVNCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVRVIHTDNGSNFTSNAVKAACWWAGIQQEFGIPYNPQSQGVVESMNKELKKIIGQVREQAEHLKTAVQMAVFIHNFKRRGGIGGYSAGERIIDIIASDIQTKELQNQILKIQNFRVYYRDSRDPIWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGADCVAGRQDED"
    ref_sanity_check = compare(ref_gene_aa_seq, pol_aa_comparison)
    print(ref_sanity_check)

    for freq_file in freq_files:
        significant_mutes = major_protein_mutations_per_sample(freq_file, ref_gene_nuc_seq, pol_start_index, pol_nuc_gene_length)

    return significant_mutes


def compare(ref_gene_aa_seq, aa_seq):
    significant_mutes = []

    # TODO- transform for loop to pandas operation
    for idx, ref_aa in enumerate(ref_gene_aa_seq):
        if ref_aa != aa_seq[idx]:
            significant_mutes.append((idx, ref_aa + aa_seq[idx]))

    return significant_mutes


def get_nuc_seq_from_freq_file_wo_dels(freq_file, ref_gene_start_index, ref_gene_length, ref_seq_for_deletions):
    df = get_nuc_seq_from_freq_file_with_dels(freq_file, ref_gene_length, ref_gene_start_index)
    # for deletions- replace with Ref base (in order to preserve all positions)
    if ref_seq_for_deletions != None:
        df.set_index('Pos').join(ref_seq_for_deletions) # TODO- complete this line (creation of Ref_for_del)
        df['Base'] = np.where(df['Base'] == '-', df['Ref_for_del'], df['Base'])
    else:
        df['Base'] = np.where(df['Base'] == '-', df['Ref'], df['Base'])


def get_nuc_seq_from_freq_file_with_dels(freq_file, ref_gene_length, ref_gene_start_index):
    df = pd.read_csv(freq_file, sep='\t')
    # filter only major variant
    df = df.loc[df['Rank'] == 0][['Pos', 'Base']]
    # filter relevant interval
    df = df.loc[(df['Pos'] > ref_gene_start_index) & (df['Pos'] < ref_gene_start_index + ref_gene_length)]
    return df


def produce_aa_seq(nuc_seq):



def major_protein_mutations_per_sample(freq_file, ref_gene_nuc_seq, ref_gene_start_index, ref_gene_length):
    # Two challenges-
    # 1. Determine significant mutations
    # 2. Determine coordinates- start + stop of each gene per sample (start with one)
    #    -> this is actually not such a challenge when you are coordinated to pipeline?

    significant_mutes = []
    min_mutations = sys.maxsize
    rf = 0

    # When using a freq files aligned to annotated Reference- reading frame is actually already known
    for i in [0]:
        # nuc_seq = get_nuc_seq_from_freq_file_wo_dels(freq_file, ref_gene_start_index + i, ref_gene_length) # assume major variant, no indels
        # aa_seq = produce_aa_seq(nuc_seq)
        # rf_significant_mutes = compare(ref_gene_aa_seq, aa_seq)

        nuc_seq_with_dels = get_nuc_seq_from_freq_file_with_dels(freq_file, ref_gene_start_index + i, ref_gene_length) # assume major variant, no indels
        rf_significant_mutes = compare_by_nucs(ref_gene_nuc_seq, nuc_seq_with_dels)

        if len(rf_significant_mutes) < min_mutations:
            significant_mutes = rf_significant_mutes
            min_mutations = len(significant_mutes)
            rf = i

        # TODO- extract minor muts too
        # TODO- manual examine known NRTI muts

        # TODO- extract advanced muts?
        # TODO- indels needs to be done manually. But significant indel areas can be marked\shown somehow

    return significant_mutes

def extract_muts_statistic():
    # how many indels\subs etc. interesting?

def compare_muts_between_samples():
    # Check diff from previous

    # Check diff from ancestor

if __name__ == "__main__":
    protein_mutations_loop()

