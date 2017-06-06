#! /usr/local/python_anaconda/bin/python3.4




def unalign_seq(seq, gap = "-"):
    """
    unalign single seq
    :param seq:  seq in seqIO format
    :param gap: gap type (default: - )
    :return: seq without gaps
    """
    seq.seq = seq.seq.ungap("-")
    return seq
