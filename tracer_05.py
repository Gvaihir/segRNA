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
argsP = parser.parse_args()


# Input sequence to align to
left = argsP.left
right = argsP.right

# Hamming distance function
#def hamming_distance(s1, s2):
 #   """Return the Hamming distance between equal-length sequences"""
  #  if len(s1) != len(s2):
   #     raise ValueError("Undefined for sequences of unequal length")
    #return sum(el1 != el2 for el1, el2 in zip(s1, s2))


# pairwise alignment function

def pw_aligner(x, left, right):
    """ x = sequence of interest, left = reference upstream cut site, right =  reference downstream cut site"""
    """ Return the field with ID, number of reads, matches up- and downstream the cut, InDel length
        position of PAM (in downstream portion), and alignment scheme"""

    # align with left
    ox = pairwise2.align.localms(x, left, 1, 0, -2, -1, one_alignment_only = True, penalize_end_gaps = False)

    # find number of matches
    for a in ox:
        ox_form = format_alignment(*a)
    ox_match = ox_form.split('\n')[1].count('|')

    # exclude aligned sequence to align to the right
    x_cut_index = max(ox_form.split('\n')[1].rfind(i) for i in '||')
    x_cut = ox_form.split('\n')[0][x_cut_index + 1:]

    # align with right
    oy = pairwise2.align.localms(x_cut, right, 1, 0, -2, -1, one_alignment_only = True, penalize_end_gaps = False)
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
    indel_ref_out = oy_form.split('\n')[2][:y_cut_index]

    # count Insertions and deletions
    insert = len(y_cut)-y_cut.count('-')
    deletion = len(right) - len(oy_form.split('\n')[0][y_cut_index:])

    # create an output
    out = ('ID {}\tReads {}\tL_match {}\tR_match {}\tInsert {}\tDel {}\tPAM {}\n'
                '{}\t{}\t{}\n'
                '{}\t{}\t{}\n'
                '{}\t{}\t{}\n'.format(
        result_toContinue[result_toContinue.Seq == x].ID.item(),
        result_toContinue[result_toContinue.Seq == x].Reads.item(),
        ox_match,
        oy_match,
        insert,
        deletion,
        y_query_out.find('GGTAAGA'),
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

    with open(argsP.i, 'r') as file:

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
        result_toContinue_pre = fastaDictDFfilter.copy()


        result_toContinue = result_toContinue_pre.sort_values(by=['Reads'], ascending=False)
        # open output file
        outPut = open(argsP.o, 'w')
        sys.stdout = outPut

        # make an output
        scheme_result = result_toContinue['Seq'].apply(pw_aligner, args=(left, right,))

        sys.stdout = original




