"""script to analyze Tracer data:
1) load fasta file;
2) filter sequencing error;
    2.1) exclude low coverage (normalize, threshold);
    2.2) exclude sequences without indels (same sequence length as reference);
    2.3) estimate Hamming distance between each pair of sequence with same length, collapse very close sequences
         (n mismatches)
3) perform pairwise alignment with reference sgRNA sequence;
4) Prepare a table of sequence alignment and statistics: indel length, pam presence

usage:
    %s < in.fa > out.txt
"""

# data frame handling package
import pandas as pd

# package to work with system inputs and outputs
import sys

# module to work with biological sequences
from Bio import pairwise2

# Import format_alignment method
from Bio.pairwise2 import format_alignment

# Regular expression operations
import re

# Input sequence to align to
left = "TGTGAACCGCATCGAGCT"
right = "GAAGGGTAAGAGCTATGCTGGAA"

# Hamming distance function
def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


# pairwise alignment function

def pw_aligner(x, y):
    """ x = sequence of interest, left = reference upstream cut site, right =  reference downstream cut site"""

    # align with left
    ox = pairwise2.align.localms(x, left, 1, 0, -1, -1, one_alignment_only = True, penalize_end_gaps = False)

    # find number of matches
    for a in ox:
        ox_form = format_alignment(*a)
    ox_match = ox_form.split('\n')[1].count('|')

    # align with right
    oy = pairwise2.align.localms(x[::-1], right[::-1], 1, 0, -1, -1, one_alignment_only = True, penalize_end_gaps = False)
    # find number of matches
    for a in oy:
        oy_form = format_alignment(*a)
    oy_match = oy_form.split('\n')[1].count('|')


    # create an output
    out = print("  ".join(str(x) for x in [ox_match, oy_match]))


if __name__ == "__main__":

    # print usage message
    # __doc__ %= sys.argv[0]
    file = open(
        "/Users/ogorodnikov/Box Sync/Ogorodnikov/LabNoteBook/05_tracer/Dry/oak180813_tracer0001_TdTtest/Collapsed/Cas9_Cre_S2_L001_R1_001.fa")
    with open(sys.argv[1], 'r') as file:

        # read as data frame
        fastaDFin = pd.read_csv(file, sep='\t', header=None)

        # split IDs and reads
        fastaDict = {'ID': fastaDFin[::2][0].str.split("-", 0).str[0].tolist(),
                     'Reads': pd.to_numeric(fastaDFin[::2][0].str.split("-", 0).str[1]).tolist(),
                     'Seq': fastaDFin[1::2][0].tolist()}

        # make data frame
        fastaDictDF = pd.DataFrame(fastaDict)

        # count length of sequence
        fastaDictDF['Length'] = fastaDictDF['Seq'].str.len()

        # normalize reads to million
        fastaDictDF['Norm'] = fastaDictDF['Reads'] * 1e6 / fastaDictDF['Reads'].sum()

        # filter out all below 5 reads/million
        fastaDictDFfilter = fastaDictDF[fastaDictDF.Norm > 5]

        # make results DF
        results = pd.DataFrame()

        # loop over every unique value and extract DFs for each
        for uniqLenght in fastaDictDFfilter['Length'].unique():

            # extrsact relevant DF
            fastaDF_rel = fastaDictDFfilter[fastaDictDFfilter.Length == uniqLenght]

            ## iterative filtering and collapsing of reads

            # create an intermediate results DF
            res = pd.DataFrame()

            # iterate over each sequence
            for sequence in fastaDF_rel['Seq']:
                # quantify Hamming distance
                fastaDF_rel_dist = fastaDF_rel
                fastaDF_rel_dist['Hamming'] = fastaDF_rel_dist['Seq'].apply(hamming_distance, args=(sequence,))

                # combine all reads with distance is <= 2
                hF = fastaDF_rel_dist[fastaDF_rel_dist.Hamming <= 2]

                # collapse reads in hF to the read with maximum alignments
                maxVal = hF['Reads'].max()
                # add all reads together for collapsed read
                resPre = hF[hF.Reads == maxVal]
                resPre['Reads'] = hF['Reads'].sum()

                # append to res DF
                res = res.append(resPre)

                # exclude collapsed sequences from fastaDF_rel (which we loop over)
                fastaDF_rel = fastaDF_rel[-fastaDF_rel['Seq'].isin(hF['Seq'])]

            results = results.append(res)

    # delete 41 nt sequences (we look for indels)
    result_toContinue = results[results.Length != 41].iloc[:,0:3]




for a in ox:
    print(format_alignment(*a))

print(contents)





#output variable
outFile = ""

#import input file (fasta)
if len(sys.argv) > 1:
    print(sys.argv)
    print(__doc__)
    sys.exit()

    records = sum(1 for _ in sys.stdin) / 4
    print("records in file", records, file = sys.stderr)


outFile = open("demofile.txt", "a")



### Some changes