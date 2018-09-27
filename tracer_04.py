"""Program to analyze Tracer data:
1) load fasta file;
2) filter sequencing error;
    2.1) exclude low coverage (normalize, threshold);
    2.2) exclude sequences without indels (same sequence length as reference);
    2.3) estimate Hamming distance between each pair of sequence with same length, collapse very close sequences
         (n mismatches)
3) perform pairwise alignment with reference sgRNA sequence;
4) Returns a table of sequence alignment and statistics: ID, number of reads, matches up- and downstream the cut,
        InDel length
        position of PAM (in downstream portion), and alignment scheme

"""

# data frame handling package
import pandas as pd

# package to work with system inputs and outputs
import sys

# module to work with biological sequences
from Bio import pairwise2

# Import format_alignment method
from Bio.pairwise2 import format_alignment


# Help message
import argparse
from argparse import RawTextHelpFormatter


parser = argparse.ArgumentParser(
    description='''Program to analyze Tracer data:
1) load fasta file;
2) filter sequencing error;
\t2.1) exclude low coverage (normalize, threshold);
\t2.2) exclude sequences without indels (same sequence length as reference);
\t2.3) estimate Hamming distance between each pair of sequence with same length, collapse very close sequences
         (n mismatches)
3) perform pairwise alignment with reference sgRNA sequence;
4) Returns a table of sequence alignment and statistics: ID, number of reads, matches up- and downstream the cut, 
        InDel length
        position of PAM (in downstream portion), and alignment scheme''', formatter_class=RawTextHelpFormatter,
    epilog="""Trace wisely""")
parser.add_argument('-i', help='Input fasta file with aligned reads in format: >ID-reads sequence. As an output'
                                    'of fastx_collapser')
parser.add_argument('-o', help='output text file')
parser.add_argument('--left', help='sequence upstream SpCas9 cut, which is 3 nt before PAM')
parser.add_argument('--right', help='sequence upstream SpCas9 cut')
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)



# Input sequence to align to
left = sys.argv[3]
right = sys.argv[4]

# Hamming distance function
def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


# pairwise alignment function

def pw_aligner(x, left, right):
    """ x = sequence of interest, left = reference upstream cut site, right =  reference downstream cut site"""
    """ Return the field with ID, number of reads, matches up- and downstream the cut, InDel length
        position of PAM (in downstream portion), and alignment scheme"""

    # align with left
    ox = pairwise2.align.localms(x, left, 1, 0, -1, -1, one_alignment_only = True, penalize_end_gaps = False)

    # find number of matches
    for a in ox:
        ox_form = format_alignment(*a)
    ox_match = ox_form.split('\n')[1].count('|')

    # exclude aligned sequence to align to the right
    x_cut_index = ox_form.split('\n')[1].rfind('|')
    x_cut = ox_form.split('\n')[0][x_cut_index + 1:]

    # align with right
    oy = pairwise2.align.localms(x_cut, right, 1, 0, -1, -1, one_alignment_only = True, penalize_end_gaps = False)
    # find number of matches
    for a in oy:
        oy_form = format_alignment(*a)
    oy_match = oy_form.split('\n')[1].count('|')

    # exclude aligned sequence to align to assay indel
    y_cut_index = oy_form.split('\n')[1].find('||')
    y_cut = oy_form.split('\n')[0][:y_cut_index]



    # for print output
    x_query_out = ox_form.split('\n')[0][:x_cut_index + 1]
    x_code_out = ox_form.split('\n')[1][:x_cut_index + 1]
    x_ref_out = ox_form.split('\n')[2][:x_cut_index + 1]

    # for print output
    y_query_out = oy_form.split('\n')[0][y_cut_index:]
    y_code_out = oy_form.split('\n')[1][y_cut_index:]
    y_ref_out = oy_form.split('\n')[2][y_cut_index:]

    # print indel
    indel_query_out = y_cut
    indel_code_out = ' ' * y_cut_index
    indel_ref_out = '*' * y_cut_index

    # create an output
    out = ('ID {}\tReads {}\tL_match {}\tR_match {}\tInDel {}\tPAM {}\n'
                '{}\t{}\t{}\n'
                '{}\t{}\t{}\n'
                '{}\t{}\t{}\n'.format(
        result_toContinue[result_toContinue.Seq == x].ID.item(),
        result_toContinue[result_toContinue.Seq == x].Reads.item(),
        ox_match,
        oy_match,
        indel_query_out.__len__(),
        y_query_out.find('GG'),
        x_query_out,
        indel_query_out,
        y_query_out,
        x_code_out,
        indel_code_out,
        y_code_out,
        x_ref_out,
        indel_ref_out,
        y_ref_out
        ))


    return(print(out))

# record std.out
original = sys.stdout

if __name__ == "__main__":

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
            result_toContinue = results[results.Length != 41].iloc[:, 0:3]

            # open output file
            outPut = open(sys.argv[2], 'w')
            sys.stdout = outPut

            # make an output
            scheme_result = result_toContinue['Seq'].apply(pw_aligner, args=(left, right,))

            sys.stdout = original




