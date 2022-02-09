#! /Users/zhaohuanan/opt/anaconda3/envs/py27/bin/python
# _*_ coding: UTF-8 _*_


# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    ---
    Author: ZHAO Herman

    Version-01:
      2020-09-13 
            1. multi-plot of target-seq
            2. modify: -i, now, the input param should be a list spite by comma (-i sample1.bmat,sample2.bmat,...)
            3. add new param: -l, it is the label of each sample (-l B-1,B-2,M1-1,M1-2,M2-1,M2-2)
            4. add new param: --to_base, what base to plot (A/G/C/T/Ins/Del or A,G, --to_base T --to_base A,T)
            5. add new param: --count_ratio, plot count info, ratio info or count and ratio info (--count_ratio all) 
            6. add new param: --output_matrix
    Version-02:
      2020-11-26
            1. fix align bugs
    E-Mail: hermanzhaozzzz@gmail.com     
    """
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = \
    """

    There will be a lot of parameters to set,

    but default paramters can handle output figure in most cases.

    """

PIPELINE = \
    """

    1. try to find all PAM and align sgRNA related to the PAM position

    2. if no appropriate alignment, try to run alignment with reverse complementary seq

    3. get the plot region

    4. make plot paramters and data

    5. plot

    """
# Learning Part END-----------------------------------------------------------

import argparse
import gzip
import os
import string
import time

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns

import Bio
from Bio import Align
from Bio import SeqIO

np.set_printoptions(suppress=True)


###############################################################################
# function part
###############################################################################
def hex_to_rgb(hex_color):
    """
    INPUT
        <hex_color> 
            Color format  like #FFFFAA
    
    RETURN
        <rgb_color>    
            RGB tuple like (255, 255, 170)
    """
    value = hex_color.lstrip('#')
    value_len = len(value)
    return( tuple(int(value[index: index + 2], 16) for index in range(0, value_len, 2)) )


def rgb_to_hex(rgb_color):
    """
    INPUT
        <rgb_color>
            Color format like (255, 255, 170)
    
    RETURN
        <hex_color>
            Color format  like #FFFFAA
    """
    rgb_color = ('#%02x%02x%02x' % rgb_color).upper()
    return(rgb_color)


def make_color_list(low_color_RGB, high_color_RGB, length_out = 20, back_format="Hex"):
    """
    INPUT
        <low_color_RGB> <high_color_RGB>
            Format like (210, 179, 150), tuple, list, or np.array
        
        <back_format>
            Hex OR RGB
        
    RETURN
        <color_list>
    """
    low_color = np.array(low_color_RGB)
    high_color = np.array(high_color_RGB)
    
    color_list = []
    for index in range(0, length_out + 1):
        rgb_color = [abs(i) for i in low_color + (high_color - low_color) // length_out * index]
        if back_format == "Hex":
            color_list.append(rgb_to_hex(tuple(rgb_color)))
        else:
            color_list.append(tuple(rgb_color))
            
    return(color_list)


def map_color(value_vec, breaks, color_list):
    """
    INPUT:
        <value_vec>
            np.array or a list of values
            
        <breaks>
            A sorted value list, which can split all num into len(color_list) intervals. 
            e.g. [0.01, 0.1, 0.5, 1] make all real num into 5 intervals, (-Inf,0.01], (0.01,0.1], (0.1, 0.5],  (0.5, 1], (1, +Inf] 
        
        <color_list>
            A hex-format color list, which have to match with breaks
    
    RETURN
        <value_color_vec>
            A list map the value_vec with breaks 
    """
    value_idx_list = []
    
    for value in value_vec:
        match_state = False
        for index, break_value in enumerate(breaks):
            if value <= break_value:
                value_idx_list.append(index)
                match_state = True
                break
        
        if not match_state:
            value_idx_list.append(index+1)
    
    return( tuple(color_list[col_idx] for col_idx in value_idx_list))


def find_seq_PAM_index(query_seq):
    """
    <INPUT>
        query_seq
    
    <HELP>
        e.g. query_seq = "AAAGAGAG"
        return index list = [1,3,5], which means search NAG, NGG at the same time
    """
    pam_index_list = []
    
    for index, base in enumerate(query_seq[:-2]):
        if query_seq[index+1:index+3] == "AG":
            pam_index_list.append((index, "NAG"))
            
        elif query_seq[index+1:index+3] == "GG":
            pam_index_list.append((index, "NGG"))
            
    return(pam_index_list)


def analysis_align_obj(alignment, reverse_state = False):
    """
    INPUT:
        <alignment obj>
        
    OUTPUT:
        <info> 
            1. match count 
            2. mismatch count 
            3. gap count 
            4. alignment.score
            
    HELP:
        2019-11-15 fix-1
            The gap count should be the num of gap contain in sgRNA alignment region.
            
            e.g.
            
            AGTGGTAAGAAGAAGACGAGACATAATGAG
            ------||||||||||||||X|----||||
            ------AAGAAGAAGACGAGCC----TGAG
            
            gap count should be 4, rather than 10.
        
        2019-11-15 fix-2
            add return info, start_index, end_index, now the retrun list will be 
            
            return_list = [
                match_count,
                mismatch_count,
                gap_count,
                alignment.score,
                start_index,
                end_index
            ]
            
            The <start_index> and <end_index> are index related to sgRNA alignment string
        
    """
    
    # define params 
    match_count = 0
    mismatch_count = 0
    gap_count = 0

    alignment_list =  str(alignment).split("\n")
    query_length = len(alignment.query)    
    
    if reverse_state:
        target_list = alignment_list[0][::-1]
        info_list = alignment_list[1][::-1]
        query_list = alignment_list[2][::-1]

    else:
        target_list = alignment_list[0]
        info_list = alignment_list[1]
        query_list = alignment_list[2]        
        
    count_state = False
    count_query_base = 0
    
    ref_align_start_index = 0
    ref_align_end_index = len(alignment.target) - 1
    
    # counting 
    for index, info_base in enumerate(info_list):
        if not count_state:
            if query_list[index] != "-":
                count_state = True
                ref_align_start_index = index
                
        if not count_state:
            continue
        
        else:
            if info_base == "|":
                match_count += 1
                count_query_base += 1

            elif info_base == "X":
                mismatch_count += 1
                count_query_base += 1

            elif info_base == "-":
                gap_count += 1
                if query_list[index] != "-":
                    count_query_base += 1 
            
        if count_query_base >= query_length:
            ref_align_end_index = index
            break
                
    return_list = [
        match_count,
        mismatch_count,
        gap_count,
        alignment.score,
        ref_align_start_index,
        ref_align_end_index
    ]
    
    return(return_list)


def sign_value(x):
    """
    HELP
        sign function
    """
    if x < 0:
        return(-1)
    elif x == 0:
        return(0)
    else:
        return(1)    


def cmp_align_list(align_a, align_b):
    """
    INPUT
        like [17, 3, 0, 73.0, 1, 22 'AAGAAGAAGACGAGTCTGCA', '||||||||||||||X|||XX', 'AAGAAGAAGACGAGCCTGAG']
    
    HELP
        compare function for align_list 
    """
    align_alphabet = {"-":0, "X":1, "|":2}
    sort_index_list = [3, 2, 1, 0]
    sort_rev_state_list = [True, False, False, True]    
    
    for order_index, align_index in enumerate(sort_index_list):
        if (align_a[align_index]  - align_b[align_index]) == 0:
                continue
        else:
            if sort_rev_state_list[order_index]:
                return( -1 * sign_value(align_a[align_index]  - align_b[align_index]))
            else:
                return( sign_value(align_a[align_index]  - align_b[align_index]))
    
    for index, char_a in enumerate(align_a[7]):
        if index <= (len(align_b[7]) - 1):
            value_a = align_alphabet[char_a]
            value_b = align_alphabet[align_b[7][index]]

            if value_a > value_b:
                return 1
            elif value_a < value_b:
                return -1 

    return 0


def run_sgRNA_alignment(align_ref_seq, align_sgRNA, sgRNA_aligner, extend_len=3):
    """
    INPUT
        <align_ref_seq>
            
        <align_sgRNA> 
            sgRNA seq without PAM
            
        <possible_sgRNA_region>
    
    RETURN
        <final_align_res_list>
    """
    align_sgRNA_rev = align_sgRNA[::-1]
    
    # find all PAM
    PAM_info_list = find_seq_PAM_index(align_ref_seq)
    if len(PAM_info_list) == 0:
        return([])
    
    # forward alignment 
    final_align_res_list = []

    for PAM_start_idx, PAM_type in PAM_info_list:

        # filter 5' end PAM 
        if (PAM_start_idx - extend_len) < len(align_sgRNA):
            continue

        # select PAM and sgRNA in the possible region 
        region_seq_start = PAM_start_idx-len(align_sgRNA) - extend_len
        region_seq_end = PAM_start_idx + 3 

        # alignment part 
        region_seq = align_ref_seq[region_seq_start : PAM_start_idx]
        align_res = sgRNA_aligner.align(region_seq[::-1], align_sgRNA_rev)

        # parse alignment
        ## if contain multiple alignment result with the same score, keep the best one;
        ## sort reason score -> gap -> mismatch -> match 
        align_res_list = []
        for align in align_res:
            align_analysis_res = analysis_align_obj(align, reverse_state=True)
            align_info_list = [ x[::-1] for x in str(align).strip().split("\n")]
            align_analysis_res += align_info_list
            align_analysis_res += [PAM_start_idx, PAM_type]
            align_res_list.append(align_analysis_res)

        if len(align_res_list) == 1:
            final_align_res_list.append(align_res_list[0])

        else:
            align_res_list.sort(cmp=cmp_align_list)
            final_align_res_list.append(align_res_list[0])

    # sort final alignment 
    final_align_res_list.sort(cmp=cmp_align_list)
    
    return(final_align_res_list)


def run_no_PAM_sgRNA_alignment_no_chop(align_ref_seq, align_sgRNA_full, no_PAM_sgRNA_aligner):
    """
    INPUT
        <align_ref_seq>
            
        <align_sgRNA> 
            sgRNA seq without PAM
            
        <no_PAM_sgRNA_aligner>
            An obj from BioPython pairwise alignment
    
    RETURN
        <final_align_res_list>
    """

    # alignment part 
    align_res = no_PAM_sgRNA_aligner.align(align_ref_seq, align_sgRNA_full)

    # parse alignment
    ## if contain multiple alignment result with the same score, keep the best one;
    ## sort reason score -> gap -> mismatch -> match 
    align_res_list = []
    for align in align_res:
        align_analysis_res_temp = analysis_align_obj(align, reverse_state=False)
        align_analysis_res = align_analysis_res_temp[:]
        align_analysis_res += str(align).strip().split("\n")

        PAM_start_index = align_analysis_res[5] - 2
        PAM_type = ref_seq[PAM_start_index : PAM_start_index+3]

        align_analysis_res += [PAM_start_index, PAM_type]
        align_res_list.append(align_analysis_res)

    if len(align_res_list) == 1:
        return (align_res_list)

    else:
        align_res_list.sort(cmp=cmp_align_list)
        return (align_res_list)
###############################################################################
# default sgRNA
###############################################################################

default_sgRNA_dict = {
    "VEGFA": "GACCCCCTCCACCCCGCCTCCGG",
    "EMX1": "GAGTCCGAGCAGAAGAAGAAGGG",
    "HEK3": "GGCCCAGACTGAGCACGTGATGG",
    "HEK4": "GGCACTGCGGCTGGAGGTGGGGG",
    "RNF2": "GTCATCTTAGTCATTACCTGAGG",
    "VEGFA-on-target": "GACCCCCTCCACCCCGCCTCCGG",
    "EMX1-on-target": "GAGTCCGAGCAGAAGAAGAAGGG",
    "HEK3-on-target": "GGCCCAGACTGAGCACGTGATGG",
    "HEK4-on-target": "GGCACTGCGGCTGGAGGTGGGGG",
    "RNF2-on-target": "GTCATCTTAGTCATTACCTGAGG"
}

###############################################################################
# main part
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="convert mpileup file to info file")

    parser.add_argument("-i", "--input_bmat",
                        help=".bmat file of Target-Seq data", required=True)

    parser.add_argument("--sgRNA",
                        help="sgRNA sequence with PAM (NGG/NAG) sequence", default=None)
    parser.add_argument("-l", "--label",
                        help="label for panels, linked by the comma", required=True)
    parser.add_argument("--to_base",
                        help="ref base to A/G/C/T/Ins/Del", default='A,G,C,T,Ins,Del')
    parser.add_argument("--count_ratio",
                        help="plot count, ratio or all", default='all')
    parser.add_argument("-o", "--out_figure",
                        help="Output figure filename", required=True)
    
    parser.add_argument("--output_matrix",
                        help="Output count and ratio matrix", default=False)
    parser.add_argument("--plot_heatmap",
                        help="Plot heatmap filename", required=True)
    parser.add_argument("--mut_direction",
                        help="[from Base][to Base], [from Base][to Base] split by comma,", default='CT,GA')

    parser.add_argument("--out_figure_format",
                        help="Support 'pdf' and 'png' result, default=pdf", default="pdf")

    parser.add_argument("--out_figure_dpi",
                        help="Out figure dpi, default=100", default="100")

    parser.add_argument("--show_indel",
                        help="If show indel info in the out figure, default=True", default="True")

    parser.add_argument("--show_index",
                        help="If show index info in the out figure, default=True", default="True")

    parser.add_argument("--block_ref",
                        help="Not show color of reference sites in the out figure, default=True", default="True")

    parser.add_argument("--box_border",
                        help="If show box border in the out figure, default=False", default="False")

    parser.add_argument("--box_space",
                        help="Set space size between two boxes, default=0.03", default="0.03")

    parser.add_argument("--min_color",
                        help="min color to plot heatmap with RGB format, default=250,239,230", default="250,239,230")

    parser.add_argument("--max_color",
                        help="max color to plot heatmap with RGB format, default=154,104,57", default="166,117,71")

    parser.add_argument("--min_ratio",
                        help="Lower than this ratio plot as white color, default=0.001", default="0.001")

    parser.add_argument("--max_ratio",
                        help="Higher than this ratio plot as max color, default=0.99", default="0.99")

    parser.add_argument("--region_extend_length",
                        help="From the middle site to extend <region_extend_length> bp at both side, default=25",
                        default="25")

    parser.add_argument("--align_settings",
                        help="Set <align_match_score>  <align_mismatch_score> <align_gap_open_score> <align_gap_extension_score>, default=5,-4,-24,-8",
                        default="5,-4,-24,-8")

    parser.add_argument("--align_min_score",
                        help="If alignment score lower than this, consider as no appropriate alignment, default=15",
                        default="15")

    ARGS = parser.parse_args()

    ###############################################################################
    # read parameters and load file
    ###############################################################################
    # set sgRNA info
    if ARGS.sgRNA in default_sgRNA_dict.keys():
        sgRNA_full_length = default_sgRNA_dict[ARGS.sgRNA]
    else:
        sgRNA_full_length = str(ARGS.sgRNA).upper()

    sgRNA_seq = sgRNA_full_length[:-3]
    sgRNA_PAM_info = sgRNA_full_length[-3:]

    # set plot region
    region_extend_length = int(ARGS.region_extend_length)

    # set min align score
    align_min_score = float(ARGS.align_min_score)

    # ---------------------------------------------------------------->>>>>
    # load .bmat file
    # ---------------------------------------------------------------->>>>>
    label_panel = ARGS.label.replace(' ', '').split(',')
    #
    ls_bmat = [i.strip() for i in ARGS.input_bmat.split(',')]
    ls_bmat_table = [pd.read_csv(path_bmat, sep='\t') for path_bmat in ls_bmat]
    for index, bmat in enumerate(ls_bmat_table):
        bmat['label'] = label_panel[index]
        bmat.drop('chr_name', axis=1, inplace=True)

    df_tmp = ls_bmat_table[0]
    bmat_table = ls_bmat_table[0]

    for df_bmat in ls_bmat_table[1:]:
        df_tmp = pd.merge(df_tmp, df_bmat, on='chr_index', how='outer')

    ls_columns = ['chr_index']
    for label in label_panel:
        str_tmp = ' '.join(
            ['ref_base_', 'A_', 'G_', 'C_', 'T_', 'del_count_', 'insert_count_', 'ambiguous_count_', 'deletion_',
             'insertion_', 'ambiguous_', 'mut_num_', 'label_ '])
        ls_tmp = str_tmp.replace('_ ', '_{label} ').format(label=label).strip().split(' ')
        ls_columns.extend(ls_tmp)

    df_tmp.columns = ls_columns

    # define ref_seq as the referencing sequence
    if len(ls_bmat_table) == 1:
        ref_seq = "".join(bmat_table['ref_base_' + label_panel[0]].tolist())

    elif len(ls_bmat_table) > 1:
        ref_seq = "".join(bmat_table.ref_base)
    else:
        raise ValueError('ref info gets wrong!')

    df_bmat_all = df_tmp.copy()
    # ---------------------------------------------------------------->>>>>
    # set alignment
    # ---------------------------------------------------------------->>>>>
    align_match, align_mismatch, align_gap_open, align_gap_extension = map(int, str(ARGS.align_settings).split(","))
    
    seq_aligner = Align.PairwiseAligner()
    seq_aligner.match = align_match
    seq_aligner.mismatch = align_mismatch
    seq_aligner.open_gap_score = align_gap_open
    seq_aligner.extend_gap_score = align_gap_extension
    seq_aligner.query_left_gap_score = 0
    seq_aligner.query_right_gap_score = 0
    seq_aligner.mode = "global"

    print("-" * 80)
    print(str(seq_aligner))
    print("-" * 80)

#     sgRNA_aligner = Align.PairwiseAligner()
#     sgRNA_aligner.match = align_match
#     sgRNA_aligner.mismatch = align_mismatch
#     sgRNA_aligner.open_gap_score = align_gap_open
#     sgRNA_aligner.extend_gap_score = align_gap_extension
#     sgRNA_aligner.query_left_gap_score = align_gap_open
#     sgRNA_aligner.query_right_gap_score = 0
#     sgRNA_aligner.mode = "global"

#     print("-" * 80)
#     print(str(sgRNA_aligner))
#     print("-" * 80)

#     no_PAM_sgRNA_aligner = Align.PairwiseAligner()
#     no_PAM_sgRNA_aligner.match = align_match
#     no_PAM_sgRNA_aligner.mismatch = align_mismatch
#     no_PAM_sgRNA_aligner.open_gap_score = align_gap_open
#     no_PAM_sgRNA_aligner.extend_gap_score = align_gap_extension
#     no_PAM_sgRNA_aligner.query_left_gap_score = 0
#     no_PAM_sgRNA_aligner.query_right_gap_score = 0
#     no_PAM_sgRNA_aligner.mode = "global"

#     print("-" * 80)
#     print(str(no_PAM_sgRNA_aligner))
#     print("-" * 80)

    # ---------------------------------------------------------------->>>>>
    # alignment
    # ---------------------------------------------------------------->>>>>
    sgRNA_align_extend_len = 3

    # make ref seq
    ref_seq_BioPy = Bio.Seq.Seq(ref_seq, Bio.Alphabet.IUPAC.unambiguous_dna)
    ref_seq_rc = str(ref_seq_BioPy.reverse_complement())

    # PAM fwd alignment
#     final_align_fwd = run_sgRNA_alignment(ref_seq, sgRNA_seq, sgRNA_aligner, sgRNA_align_extend_len)
    final_align_fwd = run_no_PAM_sgRNA_alignment_no_chop(ref_seq, sgRNA_full_length, seq_aligner)


    # PAM rev alignment
#     final_align_rev = run_sgRNA_alignment(ref_seq_rc, sgRNA_seq, sgRNA_aligner, sgRNA_align_extend_len)
    final_align_rev = run_no_PAM_sgRNA_alignment_no_chop(ref_seq_rc, sgRNA_full_length, seq_aligner)

    # fwd alignment
    print("Forward best alignment:")
    print(final_align_fwd[0][6][final_align_fwd[0][4]: final_align_fwd[0][5] + 1])
    print(final_align_fwd[0][7][final_align_fwd[0][4]: final_align_fwd[0][5] + 1])
    print(final_align_fwd[0][8][final_align_fwd[0][4]: final_align_fwd[0][5] + 1])
    print(final_align_fwd[0][3])
    print("-" * 80)

    # rev alignment
    print("Reverse best alignment:")
    print(final_align_rev[0][6][final_align_rev[0][4]: final_align_rev[0][5] + 1])
    print(final_align_rev[0][7][final_align_rev[0][4]: final_align_rev[0][5] + 1])
    print(final_align_rev[0][8][final_align_rev[0][4]: final_align_rev[0][5] + 1])
    print(final_align_rev[0][3])
    print("-" * 80)

    # make alignment info
    sgRNA_align = [""] * len(ref_seq)
    sgRNA_align_insert = [""] * len(ref_seq)

    final_align_direction = None
    DNA_rev_cmp_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "-": "-"}

    # define align direction and final align res
    if final_align_fwd[0][3] >= final_align_rev[0][3]:
        if final_align_fwd[0][3] >= align_min_score:
            final_align_direction = "Forward PAM Alignment"
            final_align = final_align_fwd[0]

    elif final_align_rev[0][3] > final_align_fwd[0][3]:
        if final_align_rev[0][3] >= align_min_score:
            final_align_direction = "Reverse PAM Alignment"
            final_align = final_align_rev[0]

    if final_align_direction == None:
#         final_align_direction = "No PAM Alignment"
#         final_align = run_no_PAM_sgRNA_alignment_no_chop(ref_seq, sgRNA_full_length, no_PAM_sgRNA_aligner)[0]
        raise IOError("Alignment Error!")

        
    # make sgRNA alignment info
    final_align_info = final_align[7][final_align[4]: final_align[5] + 1]
    final_align_ref = final_align[6][final_align[4]: final_align[5] + 1]
    final_align_sgRNA = final_align[8][final_align[4]: final_align[5] + 1]
    final_align_ref_gap_count = final_align_ref.count("-")
    final_align_sgRNA_gap_count = final_align_sgRNA.count("-")

    ref_align_gap_count = 0
    ref_del_str = ""

    if (final_align_direction == "Forward PAM Alignment") or (final_align_direction == "Reverse PAM Alignment"):
        sgRNA_start = final_align[9] - len(sgRNA_seq) - final_align_sgRNA_gap_count + final_align_ref_gap_count

        for align_index, align_ref in enumerate(final_align_ref):
            if align_ref != "-":
                sgRNA_align[sgRNA_start + align_index - ref_align_gap_count] = ref_del_str + final_align_sgRNA[
                    align_index]
                ref_del_str = ""
            else:
                ref_align_gap_count += 1
                ref_del_str += final_align_sgRNA[align_index]
                sgRNA_align_insert[sgRNA_start + align_index] = final_align_sgRNA[align_index]

        # add PAM info
        if final_align_direction == "Reverse Alignment":
            sgRNA_align = sgRNA_align[::-1]
            sgRNA_align_insert = sgRNA_align_insert[::-1]
#         PAM_ref_start_index = sgRNA_start + align_index + 1 - ref_align_gap_count
#         sgRNA_align[PAM_ref_start_index: PAM_ref_start_index + 3] = sgRNA_full_length[-3:]

#         if final_align_direction == "Reverse PAM Alignment":
#             sgRNA_align = sgRNA_align[::-1]
#             sgRNA_align_insert = sgRNA_align_insert[::-1]

#     elif final_align_direction == "No PAM Alignment":
#         sgRNA_start = final_align[4]

#         for align_index, align_ref in enumerate(final_align_ref):
#             if align_ref != "-":
#                 sgRNA_align[sgRNA_start + align_index - ref_align_gap_count] = final_align_sgRNA[align_index]
#             else:
#                 ref_align_gap_count += 1
#                 sgRNA_align_insert[sgRNA_start + align_index] = final_align_sgRNA[align_index]

    # set possible_sgRNA_region
    sgRNA_align_start = final_align[4]    
    sgRNA_align_end = final_align[5]
    
    
    possible_sgRNA_region_start = max(sgRNA_align_start - region_extend_length, 0)
    possible_sgRNA_region_end = min(sgRNA_align_end + region_extend_length, len(ref_seq) - 1)
    possible_sgRNA_region = [possible_sgRNA_region_start, possible_sgRNA_region_end]
#     if final_align_direction == None:
#         possible_sgRNA_region = [len(ref_seq) // 2 - region_extend_length, len(ref_seq) // 2 + region_extend_length]

#     else:
#         if final_align_direction == "No PAM Alignment":
#             sgRNA_align_start = final_align[4]
#             sgRNA_align_end = final_align[5]

#         else:
#             if final_align_direction == "Forward PAM Alignment":
#                 print 1
#                 sgRNA_align_end = final_align[9] + 3
#                 sgRNA_align_start = sgRNA_align_end - len(sgRNA_full_length)

#             elif final_align_direction == "Reverse PAM Alignment":
#                 print 2
#                 sgRNA_align_start = len(ref_seq) - (final_align[9] + 3)
#                 sgRNA_align_end = sgRNA_align_start + len(sgRNA_full_length)

#         possible_sgRNA_region_start = max(sgRNA_align_start - region_extend_length, 0)
#         possible_sgRNA_region_end = min(sgRNA_align_end + region_extend_length, len(ref_seq) - 1)
#         possible_sgRNA_region = [possible_sgRNA_region_start, possible_sgRNA_region_end]
    # select bmat_table
    bmat_table_select = df_bmat_all[possible_sgRNA_region[0]: possible_sgRNA_region[1]]
    # ---------------------------------------------------------------->>>>>
    # make plot
    # ---------------------------------------------------------------->>>>>

    # --------------------------------------------------->>>>>
    # set color
    # --------------------------------------------------->>>>>
    # show indel
    indel_plot_state = eval(ARGS.show_indel)
    index_plot_state = eval(ARGS.show_index)
    box_border_plot_state = eval(ARGS.box_border)

    # set panel size
    panel_box_width = 0.4
    panel_box_heigth = 0.4
    panel_space = 0.05
    panel_box_space = float(ARGS.box_space)

    # color part
    base_color_dict = {"A": "#04E3E3", "T": "#F9B874", "C": "#B9E76B", "G": "#F53798", "N": "#AAAAAA", "-": "#AAAAAA"}

    # make color breaks
    color_break_num = 20
    break_step = 1.0 / color_break_num
    min_color_value = float(ARGS.min_ratio)
    max_color_value = float(ARGS.max_ratio)
    color_break = np.round(np.arange(min_color_value, max_color_value, break_step), 5)

    # make color list
    # low_color = (210, 179, 150)
    # low_color = (240, 225, 212)
    # low_color = (250, 239, 230)
    # high_color = (154, 104, 57)
    # low_color = (217, 231, 245)
    # high_color = (9, 42, 96)

    low_color = tuple(map(int, ARGS.min_color.split(",")))
    high_color = tuple(map(int, ARGS.max_color.split(",")))

    try:
        color_list = make_color_list(low_color, high_color, len(color_break) - 1, "Hex")
        color_list = ["#FFFFFF"] + color_list
    except:
        print low_color, high_color
        print color_break

    # get plot info
    total_box_count = possible_sgRNA_region[1] - possible_sgRNA_region[0]

    # calculate base info and fix zero
    ls_base_sum_count = []
    ls_total_sum_count = []
    for label in label_panel:
        base_sum_count = bmat_table_select[["A_%s" % label, "G_%s" % label, "C_%s" % label, "T_%s" % label]].apply(
            lambda x: x.sum(), axis=1)
        # print base_sum_count
        total_sum_count = bmat_table_select[["A_%s" % label, "G_%s" % label, "C_%s" % label, "T_%s" % label,
                                             "del_count_%s" % label, "insert_count_%s" % label]].apply(
            lambda x: x.sum(), axis=1)

        base_sum_count[base_sum_count == 0] = 1
        
        # base_sum_count[np.where(base_sum_count == 0)[0]] = 1
        # total_sum_count[np.where(total_sum_count == 0)[0]] = 1
        base_sum_count = base_sum_count + 1
        total_sum_count = total_sum_count + 1
        
        ls_base_sum_count.append(base_sum_count)
        ls_total_sum_count.append(total_sum_count)
    # make plot size
    plot_data_list = None
    panel_space_coef = None
    panel_height_coef = None
    # check to_base:  like AGCTInsDel or combine like A,G,T,Ins
    ls_base_raw = ['A', 'G', 'C', 'T', 'Ins', 'Del']
    if ARGS.to_base in ls_base_raw:
        ls_to_base_raw = [ARGS.to_base]
    elif ',' in ARGS.to_base:
        ls_to_base_raw = [i.strip() for i in ARGS.to_base.strip().split(',')]
        for base in ls_to_base_raw:
            if base in ls_base_raw:
                pass
            else:
                print 'It seems like that [to_base] param has something wrong.'
    else:
        print 'It seems like that [to_base] param has something wrong.'
    # to_base = ARGS.to_base # A G C T Del Ins all
    print 'Split the to_base params:'
    print '    ', ls_to_base_raw
    count_ratio = ARGS.count_ratio  # count ratio all
    print 'Plot count and/or ratio:'
    print '    ', count_ratio

    # panel_height_coef = [0.5, 0.9, 0.9] + [0.5] * 6 + [0.5] * 6
    if count_ratio == 'all':

        panel_height_coef = [0.5, 0.9, 0.9] + [0.5] * len(ls_bmat_table) * len(ls_to_base_raw) * 2
    else:
        panel_height_coef = [0.5, 0.9, 0.9] + [0.5] * len(ls_bmat_table) * len(ls_to_base_raw)
    # 根据bmat表格个数来控制方格高度
    # panel_space_coef = [1, 1, 1] + [0.3] * 3 + [1, 0.3, 1] + [0.3] * 3 + [1, 0.3]
    if count_ratio == 'all':
        panel_space_coef = [1, 1, 1] + ([0.3] * (len(ls_bmat_table) - 1) + [4]) * len(ls_to_base_raw) * 2
    else:
        panel_space_coef = [1, 1, 1] + ([0.3] * (len(ls_bmat_table) - 1) + [4]) * len(ls_to_base_raw)
    # 更正heatmap align错误
#     plot_heatmap_index = bmat_table_select.chr_index
    plot_data_list = [
        ["Index", bmat_table_select.chr_index],
        ["On-target", sgRNA_align[possible_sgRNA_region[0]: possible_sgRNA_region[1]]],
        ["Ref", bmat_table_select['ref_base_%s' % label_panel[0]]]
    ]
    if count_ratio == 'count' or count_ratio == 'all':
        for to_base in ls_to_base_raw:
            if to_base == 'A':
                for label in label_panel:
                    plot_data_list.append(
                        ["{label}: to A".format(label=label),
                         np.array(bmat_table_select["A_{label}".format(label=label)])]
                    )
            if to_base == 'G':
                for label in label_panel:
                    plot_data_list.append(
                        ["{label}: to G".format(label=label),
                         np.array(bmat_table_select["G_{label}".format(label=label)])]
                    )
            if to_base == 'C':
                for label in label_panel:
                    plot_data_list.append(
                        ["{label}: to C".format(label=label),
                         np.array(bmat_table_select["C_{label}".format(label=label)])]
                    )
            if to_base == 'T':
                for label in label_panel:
                    plot_data_list.append(
                        ["{label}: to T".format(label=label),
                         np.array(bmat_table_select["T_{label}".format(label=label)])]
                    )
            if to_base == 'Del':
                for label in label_panel:
                    plot_data_list.append(
                        ["{label}: to Del".format(label=label),
                         np.array(bmat_table_select["del_count_{label}".format(label=label)])]
                    )
            if to_base == 'Ins':
                for label in label_panel:
                    plot_data_list.append(
                        ["{label}: to Ins".format(label=label),
                         np.array(bmat_table_select["insert_count_{label}".format(label=label)])]
                    )
    if count_ratio == 'ratio' or count_ratio == 'all':
        for to_base in ls_to_base_raw:
            if to_base == 'A':
                for index, label in enumerate(label_panel):
                    plot_data_list.append(
                        ["{label}: to A(%)".format(label=label),
                         np.array(bmat_table_select["A_{label}".format(label=label)] / ls_base_sum_count[index])]
                    )
            if to_base == 'G':
                for index, label in enumerate(label_panel):
                    plot_data_list.append(
                        ["{label}: to G(%)".format(label=label),
                         np.array(bmat_table_select["G_{label}".format(label=label)] / ls_base_sum_count[index])]
                    )
            if to_base == 'C':
                for index, label in enumerate(label_panel):
                    plot_data_list.append(
                        ["{label}: to C(%)".format(label=label),
                         np.array(bmat_table_select["C_{label}".format(label=label)] / ls_base_sum_count[index])]
                    )
            if to_base == 'T':
                for index, label in enumerate(label_panel):
                    plot_data_list.append(
                        ["{label}: to T(%)".format(label=label),
                         np.array(bmat_table_select["T_{label}".format(label=label)] / ls_base_sum_count[index])]
                    )
            if to_base == 'Del':
                for index, label in enumerate(label_panel):
                    plot_data_list.append(
                        ["{label}: to Del(%)".format(label=label), np.array(
                            bmat_table_select["del_count_{label}".format(label=label)] / ls_total_sum_count[index])]
                    )
            if to_base == 'Ins':
                for index, label in enumerate(label_panel):
                    plot_data_list.append(
                        ["{label}: to Ins(%)".format(label=label), np.array(
                            bmat_table_select["insert_count_{label}".format(label=label)] / ls_total_sum_count[index])]
                    )

    # get box and space info
    box_height_list = np.array(panel_height_coef) * panel_box_heigth
    panel_space_list = np.array(panel_space_coef) * panel_space
    # for i, j in enumerate(panel_space_coef):
    #     print i,j
    # print panel_space
    # calculate figure total width and height
    figure_width = total_box_count * panel_box_width + (total_box_count - 1) * panel_box_space + panel_box_width * 2
    figure_height = sum(box_height_list) + sum(panel_space_list)
    # make all box_x
    box_x_vec = np.arange(0, figure_width + panel_box_width, panel_box_width + panel_box_space)
    box_x_vec = box_x_vec[:(len(ref_seq) + 1)]

    # make box border
    if box_border_plot_state:
        box_edgecolor = "#AAAAAA"
        box_linestyle = "-"
        box_linewidth = 2
    else:
        box_edgecolor = "#FFFFFF"
        box_linestyle = "None"
        box_linewidth = 0

    # make box_y initialize
    current_y = 0

    # ---------------------------------------------------------------->>>>>>>>
    # plot region
    # ---------------------------------------------------------------->>>>>>>>

    # make plot

    # output_matrix state
    print 'Output matrix state: True'

    bl_matrix = True
    df_matrix = pd.DataFrame(np.zeros((len(panel_height_coef), len(plot_data_list[0][1]))))

    ls_row_name = []
    for panel_index in range(len(panel_height_coef)):
        ls_row_name.append(plot_data_list[panel_index][0])

        for index, box_value in enumerate(plot_data_list[panel_index][1]):
            row = panel_index
            col = index
            df_matrix.iloc[row, col] = box_value
    df_matrix.index = ls_row_name

    if bool(ARGS.output_matrix):
        bmat_table_select = bmat_table_select.copy()
        bmat_table_select['On-target'] = df_matrix.loc['On-target', :].tolist()
        bmat_table_select.to_csv(ARGS.output_matrix, sep='\t', index=None)

    # plot heatmap

    heatmap_state = bool(ARGS.plot_heatmap)
    if heatmap_state:
        heatmap_path = ARGS.plot_heatmap.strip()
        # 预设参数
        num_extend = int(ARGS.region_extend_length)
        input_mut_direction = ARGS.mut_direction
        # rename columns of df_matrix
        df_matrix.columns = range(1, df_matrix.shape[1] + 1)
        # print df_matrix

        # 解析参数
        ls_input_mut_direction = [base for base in input_mut_direction.replace(' ', '').split(',')]
        ls_input_mut_direction

        dt_base = {'A': [], 'G': [], 'C': [], 'T': []}
        for info in ls_input_mut_direction:
            dt_base[info[0]] += info[1]
        print '[mut_direction] for heatmap: %s' % dt_base
        print '[label_panel] for heatmap: %s' % label_panel
        print '[region_extend_length] for heatmap: %s' % region_extend_length
        print '[block_ref] for Target-seq multiplot: %s' % ARGS.block_ref

        ls_col_not_null = df_matrix.columns[-(df_matrix.loc['On-target', :] == '')].tolist()
        ls_col_all = df_matrix.columns.tolist()
        ls_plot_range = list(range(ls_col_not_null[0] - num_extend, ls_col_not_null[-1] + num_extend))

        try:
            df_matrix = df_matrix.loc[:, ls_plot_range]
            # 有不在的就不为空，不为空则判断成立
            # print 'testTTTTTTTTTTTTTTT'

        except:
            raise ValueError('The param [num_extend] is too large or <0, it must be a integer>=0')
        ls_on_target = df_matrix.loc['On-target', :].tolist()
        ls_ref = df_matrix.loc['Ref', :].tolist()

        ls_ratio = []
        for label in label_panel:
            df_sample = df_matrix[[label == i.split(': to')[0] and '(%)' not in i for i in df_matrix.index]]
            df_sample.index = ['A', 'G', 'C', 'T', 'Ins', 'Del']
            ls_A = df_sample.T['A'].isnull().tolist()
            ls_G = df_sample.T['G'].isnull().tolist()
            ls_C = df_sample.T['C'].isnull().tolist()
            ls_T = df_sample.T['T'].isnull().tolist()
            ls_bool = []

            for i in range(len(ls_A)):
                if ls_A[i] & ls_C[i] & ls_G[i] & ls_T[i]:
                    ls_bool.append(False)
                else:
                    ls_bool.append(True)
            ls_bl_select_df_sample = ls_bool
            df_sample = df_sample.loc[:, ls_bl_select_df_sample]
            df_sample = df_sample.loc[['A', 'G', 'C', 'T'], :].T.copy()
            df_sample = df_sample.applymap(int)
            df_sample['sum'] = df_sample['A'] + df_sample['G'] + df_sample['C'] + df_sample['T']
            df_sample['Ref'] = np.array(ls_ref)[ls_bl_select_df_sample]
            df_sample['To_base'] = df_sample['Ref'].map(lambda x: dt_base[x])
            df_sample['Mut_count'] = 0

            for recode in df_sample.iterrows():
                index_row = recode[0]
                for to_base in recode[1]['To_base']:
                    df_sample.loc[index_row, 'Mut_count'] += recode[1][to_base]

            df_sample['Mut_ratio'] = df_sample['Mut_count'] / df_sample['sum']
            ls_sample_ratio = df_sample['Mut_ratio'].tolist()
            ls_ratio.append(ls_sample_ratio)

        df_ratio_all = pd.DataFrame(ls_ratio).T
        df_ratio_all.columns = label_panel
        

        try:
            df_ratio_all.index = df_matrix.T['Index'][ls_bl_select_df_sample].map(float).map(int).tolist()
            df_ratio_all['On-target'] = df_matrix.T['On-target'][ls_bl_select_df_sample].tolist()
            df_ratio_all['Ref'] = df_matrix.T['Ref'][ls_bl_select_df_sample].tolist()
        except ValueError:
            start_idx = df_matrix.T['Index'][ls_bl_select_df_sample].map(float).map(int).tolist()[0]
            end_idx = start_idx + len(df_ratio_all.index.tolist())
            df_ratio_all.index = np.arange(start_idx, end_idx)
            df_ratio_all['On-target'] = df_matrix.T['On-target'].tolist()
            df_ratio_all['Ref'] = df_matrix.T['Ref'].tolist()


        
        # df_ratio_all.index = df_matrix.T['Index'].map(float).map(int).tolist()
        # df_ratio_all['On-target'] = df_matrix.T['On-target'].tolist()
        # df_ratio_all['Ref'] = df_matrix.T['Ref'].tolist()
        
        df_ratio_all['On-Target'] = 0
        df_ratio_all['Reference'] = 0
        df_ratio_all = df_ratio_all[['On-target', 'Ref', 'On-Target', 'Reference'] + label_panel].T.copy()

        base_color_dict = {"A": "#04E3E3", "T": "#F9B874", "C": "#B9E76B", "G": "#F53798", "N": "#FFFFFF"}


        def plot_agct(x):
            if x == '':
                return base_color_dict['N']  # white
            elif x in ['A', 'G', 'C', 'T']:
                return base_color_dict[x]
            else:
                return '#AAAAAA'


        df_onTarget_Ref = df_ratio_all.iloc[:2, :].fillna(' ')
        df_onTarget_Ref_color = df_ratio_all.iloc[:2, :].fillna('').applymap(plot_agct)
        df_plot = df_ratio_all.iloc[2:, :].copy()



    def map_hex_for_matrix(x):
        # 判断value范围并返回颜色的Hex值
        for value in ls_break:
            if x < value:
                break
            else:
                #             print('x:%s, value:%s' % (x, value))
                continue
        return ls_color[ls_break.index(value) - 1]


    base_color_dict = {"A": "#04E3E3", "T": "#F9B874", "C": "#B9E76B", "G": "#F53798", "N": "#FFFFFF"}


    def plot_agct(x):
        if x == ' ':
            return base_color_dict['N']  # white
        elif x in ['A', 'G', 'C', 'T']:
            return base_color_dict[x]
        else:
            return '#FFFFFF'


#     ls_break = [0, 0.05, 0.1, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100]
#     ls_break = [0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
#     ls_break = [0, 0.003, 0.006, 0.009, 0.012, 0.015, 0.018, 0.022, 0.026, 0.030, 0.035, 0.040, 0.05, 0.06, 0.07, 0.08, 1.00, 2.00, 3.00]
    df_plot_rec = df_plot.applymap(float) * 100
    ls_max = []
    for i in df_plot_rec.values.tolist():
        ls_max.extend(i)
    ls_break = list(np.arange(0,max(ls_max),max(ls_max)/100))
    ls_color_middle = make_color_list(hex_to_rgb('#87BDDB'), hex_to_rgb('#0A306A'),length_out=80, back_format="Hex")
    ls_color_top = ['#0A306A'] * 20
    ls_color_bottom = ['#EFEFEF'] * 5 + ['#DDEAF6'] * 5 + ['#A9CEE4'] * 5 + ['#87BDDB'] * 5
    ls_color = ls_color_bottom + ls_color_middle + ls_color_top
#     print(df_plot_rec.values.tolist())
    print(ls_color)
    print(ls_break)
#     ls_color = [
#         '#EFEFEF',
#         '#DDEAF6',
#         '#A9CEE4',
#         '#87BDDB',
#         '#6AABD3',
#         '#57A0CD',
#         '#4691C6',  
#         '#3A88C0',
#         '#2D7EBA',
#         '#2272B4',
#         '#1C67AD',
#         '#1A60A7',
#         '#1757A0',
#         '#164F99',
#         '#12458B',
#         '#0F3F82',
#         '#0D3776',
#         '#0A306A'
#     ]
    # data for test
    # df_plot = pd.read_csv('./df_plot.csv', index_col=0)
    # df_onTarget_Ref_color = pd.read_csv('./df_onTarget_Ref_color.csv', index_col=0)
    # df_onTarget_Ref = pd.read_csv('./df_onTarget_Ref.csv', index_col=0)
    df_onTarget_Ref_tmp = df_onTarget_Ref.T
    df_onTarget_Ref_tmp['On-target'][df_onTarget_Ref_tmp['On-target'].isnull()] = ' '
    df_onTarget_Ref = df_onTarget_Ref_tmp.T
    # heatmap_path = './jupyter-lab_test_heatmap.pdf'
    # num_extend = 10
    # input_mut_direction = 'GA,CT'
    df_plot_rec = df_plot.applymap(float) * 100
    # heatmap颜色，去map_hex_for_matrix函数中调整
    df_plot_rec_cmap = df_plot_rec.applymap(map_hex_for_matrix)
    df_plot_rec.loc['On-Target', :] = df_onTarget_Ref.loc['On-target', :]
    df_plot_rec.loc['Reference', :] = df_onTarget_Ref.loc['Ref', :]
    df_plot_rec_cmap.loc['On-Target', :] = df_onTarget_Ref.loc['On-target', :].map(plot_agct)
    df_plot_rec_cmap.loc['Reference', :] = df_onTarget_Ref.loc['Ref', :].map(plot_agct)



    figure_width_heatmap = df_plot_rec.shape[1]
    figure_height_heatmap = max(df_plot_rec.shape[0], len(ls_break)*0.8)
    scale = 1.1  # 1.1
    fig_heatmap = plt.figure(figsize=(figure_width_heatmap * scale, figure_height_heatmap * scale))
    ax = fig_heatmap.add_subplot(111, aspect="equal")

    plt.xlim([0, figure_width_heatmap + 5])
    plt.ylim([-figure_height_heatmap - 0.5, 0])
    plt.axis("off")
    
    # export heatmap values
    print("export heatmap reference table...")
    df_plot_rec.to_csv(heatmap_path.replace(".pdf",".csv").replace(".PDF",".csv"))

    for row in range(df_plot_rec.shape[0]):
        row_name = df_plot_rec.index[row]

        # add panel name
        ax.text(
            x=-1,
            y=0.55 - row - 1.1,
            s=row_name,
            horizontalalignment='right',
            verticalalignment='center',
            fontsize=34 / 1.1 * scale,
            fontname="DejaVu Sans",
            alpha=1,
        )

        # add on-target and ref
        if row_name in ['On-Target', 'Reference']:
            for col in range(df_plot.shape[1]):
                # plot Rectangle
                site_x = col
                site_y = - row + 0.05 - 1.05
                ax.add_patch(matplotlib.patches.Rectangle(
                    (site_x, site_y),
                    width=1,
                    height=1,
                    linestyle='-',
                    fill=True,
                    facecolor=df_plot_rec_cmap.iloc[row, col],
                    edgecolor='#AAAAAA' if df_plot_rec_cmap.iloc[row, col] != '#FFFFFF' else '#FFFFFF',
                    linewidth=3.5

                ))
                # plot text
                site_x = col
                site_y = -row - 0.55
                ax.text(
                    x=site_x + 0.5,
                    y=site_y,
                    s=df_plot_rec.iloc[row, col],
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=34 / 1 * scale,
                    fontname="DejaVu Sans",
                    alpha=1
                )


        # add seq panels
        else:
            for col in range(df_plot.shape[1]):
                # plot Rectangle
                site_x = col
                site_y = -row - 1.05
                ax.add_patch(matplotlib.patches.Rectangle(
                    (site_x, site_y),
                    width=1,
                    height=1,
                    linestyle='-',
                    fill=True,
                    facecolor=df_plot_rec_cmap.iloc[row, col],
                    edgecolor='#FFFFFF',
                    linewidth=3

                ))

    # add cbar
    # start from this site
    site_x = df_plot.shape[1] + 1
    site_y = -figure_height_heatmap * 0.1 - 1
    step_scale = 0.1
    print '[cbar_scale]: %s' % step_scale


    for color in ls_color:
        site_y += step_scale
        ax.add_patch(matplotlib.patches.Rectangle(
            (site_x, site_y),
#             width=step_scale,
            width=1,
            height=step_scale,
#             linestyle='-',
            fill=True,
            facecolor=color,
#             edgecolor='#AAAAAA',
#             linewidth=3 * step_scale

        ))

    # add cbar text
    # start from this site
    site_x = df_plot.shape[1] + 1
    site_y = -figure_height_heatmap * 0.1
    for idx,label in enumerate(ls_break):
        site_y += step_scale
        if idx%10==0:
            ax.text(
                x=site_x + 1.1,
                y=site_y - 1,
                s=round(label,4),
                horizontalalignment='left',
                verticalalignment='center',
                fontsize=34 / 1.5,
    #             fontsize=34 / 1.5 * step_scale,
                fontname="DejaVu Sans",
                alpha=1
            )
    # cax = plt.gcf().axes[-1]
    # cax.tick_params(labelsize=44, direction='in', top=False, bottom=False, left=False, right=False)

    fig_heatmap.savefig(heatmap_path, bbox_inches='tight')  # 减少边缘空白
    # plt.show()



    # set new figure
    fig = plt.figure(figsize=(figure_width * 1.1, figure_height * 1.1))
    ax = fig.add_subplot(111, aspect="equal")
    plt.xlim([0, figure_width])
    plt.ylim([-figure_height, 0])
    plt.axis("off")

    text_list = []
    patches = []
    # print panel_height_coef
    for panel_index in range(len(panel_height_coef)):
        # print 'panel index:', panel_index
        # panel name
        panel_name = plot_data_list[panel_index][0]
        panel_name_x = box_x_vec[0]
        # print panel_name
        # print panel_name_x
        panel_name_y = current_y - box_height_list[panel_index] * 0.5
        # print panel_name_y
        text_list.append((panel_name_x, panel_name_y, panel_name, 10))

        # plot panel box
        # plot Index行
        if panel_name == "Index":
            # don't draw box, only add text
            for index, box_value in enumerate(plot_data_list[panel_index][1]):
                box_x = box_x_vec[index + 1]
                text_list.append(
                    (box_x + panel_box_width * 0.5, current_y - box_height_list[panel_index] * 0.5, str(box_value), 6))

            # make next panel_y
            current_y = current_y - (box_height_list[panel_index] + panel_space_list[panel_index])
        # plot On-target和Ref行
        elif panel_name in ["On-target", "Ref"]:
            # if panel_name = Ref, form a new list to store plot info for checking the block_ref param
            if panel_name == 'Ref':
                plot_data_list_ref = plot_data_list[panel_index][1]
            for index, box_value in enumerate(plot_data_list[panel_index][1]):
                if box_value == "":
                    box_fill = False
                    box_color = "#FFFFFF"

                else:
                    if "Reverse" in final_align_direction:
                        if panel_name == "Ref":
                            box_value = "".join([DNA_rev_cmp_dict.get(x) for x in box_value])
                        else:
                            ### fix bug
                            for x in box_value:
                                box_value = DNA_rev_cmp_dict.get(x)

                        box_color = base_color_dict.get(box_value[0])

                    else:
                        box_color = base_color_dict.get(box_value[-1])

                    if not box_color:
                        box_fill = False
                        box_color = "#FFFFFF"
                    else:
                        box_fill = True

                box_x = box_x_vec[index + 1]

                patches.append(Rectangle(
                    xy=(box_x, current_y - box_height_list[panel_index]),
                    width=panel_box_width,
                    height=box_height_list[panel_index],
                    fill=box_fill,
                    alpha=1,
                    linestyle=box_linestyle,
                    linewidth=box_linewidth,
                    edgecolor=box_edgecolor,
                    facecolor=box_color)
                )

                # text
                text_list.append(
                    (box_x + 0.5 * panel_box_width, current_y - 0.5 * box_height_list[panel_index], str(box_value), 16))

            # make next panel_y
            current_y = current_y - (box_height_list[panel_index] + panel_space_list[panel_index])
        # 用来plot Sample行的
        # 这里plot counts
        elif (panel_name.split('to ')[-1] in ls_base_raw):
            # print 'panel name:', panel_name.split('to ')[-1].replace('(%)', '')
            # print panel_name, 'check panel name'
            # 改进了这里total 和 base的值的引用
            if count_ratio == 'all':
                ls_base_sum_count_all = ls_base_sum_count * 2 * len(ls_to_base_raw)
                ls_total_sum_count_all = ls_total_sum_count * 2 * len(ls_to_base_raw)
                total_sum_count = ls_total_sum_count_all[panel_index - 3]
                base_sum_count = ls_base_sum_count_all[panel_index - 3]
            else:
                ls_base_sum_count_all = ls_base_sum_count * len(ls_to_base_raw)
                ls_total_sum_count_all = ls_total_sum_count * len(ls_to_base_raw)
                total_sum_count = ls_total_sum_count_all[panel_index - 3]
                base_sum_count = ls_base_sum_count_all[panel_index - 3]

            panel_name = panel_name.split('to ')[-1]

            if panel_name in ["Del", "Ins"]:
                box_ratio = plot_data_list[panel_index][1] / total_sum_count
            else:
                box_ratio = plot_data_list[panel_index][1] / base_sum_count

            box_color_list = map_color(box_ratio, color_break, color_list)
            ls_base_sum = base_sum_count.tolist()
            for index, box_value in enumerate(plot_data_list[panel_index][1]):
                # if count_ratio == 'all':
                box_color = box_color_list[index]
                # check block_ref state: count plot
                if ARGS.block_ref == 'True' or ARGS.block_ref == True:
                    if plot_data_list_ref.iloc[index] == panel_name:
                        # if box_value >= ls_base_sum[index] * 0.8:
                        # box_fill = False
                        box_color = "#FFFFFF"

                box_x = box_x_vec[index + 1]
                patches.append(Rectangle(
                    xy=(box_x, current_y - box_height_list[panel_index]),
                    width=panel_box_width,
                    height=box_height_list[panel_index],
                    fill=True,
                    alpha=1,
                    linestyle=box_linestyle,
                    linewidth=box_linewidth,
                    edgecolor=box_edgecolor,
                    facecolor=box_color)
                )
                # text
                text_list.append(
                    (box_x + 0.5 * panel_box_width, current_y - 0.5 * box_height_list[panel_index], str(box_value), 6))
                # print box_value
            # make next panel_y
            current_y = current_y - (box_height_list[panel_index] + panel_space_list[panel_index])

        # 这里plot ratio(剩下的都是base(%))
        # 注意各个参数的length和对应关系,最后一个BUG就在这儿![已解决]
        else:
            box_color_list = map_color(plot_data_list[panel_index][1], color_break, color_list)
            # print box_color_list
            panel_name = panel_name.split('to ')[-1].replace('(%)', '')

            for index, box_value in enumerate(plot_data_list[panel_index][1]):
                # ratio
                box_color = box_color_list[index]

                # check block_ref state: ratio plot
                if ARGS.block_ref == 'True' or ARGS.block_ref == True:
                    if plot_data_list_ref.iloc[index] == panel_name:
                        # if box_value >= 0.85:
                        # box_fill = False
                        box_color = "#FFFFFF"

                box_x = box_x_vec[index + 1]
                patches.append(Rectangle(
                    xy=(box_x, current_y - box_height_list[panel_index]),
                    width=panel_box_width,
                    height=box_height_list[panel_index],
                    fill=True,
                    alpha=1,
                    linestyle=box_linestyle,
                    linewidth=box_linewidth,
                    edgecolor=box_edgecolor,
                    facecolor=box_color)
                )

                if '(%)' in plot_data_list[panel_index][0]:
                    # text
                    # ref
                    # print 'test', box_value
                    # if ARGS.block_ref == True:
                    #     if box_value >= 0.99:
                    #         box_value = 0

                    text_list.append((box_x + 0.5 * panel_box_width, current_y - 0.5 * box_height_list[panel_index],
                                      round(box_value * 100, 4), 6))

                else:
                    text_list.append((box_x + 0.5 * panel_box_width, current_y - 0.5 * box_height_list[panel_index],
                                      box_value, 6))

            # make next panel_y
            if panel_index < len(panel_space_list):
                current_y = current_y - (box_height_list[panel_index] + panel_space_list[panel_index])

    # plot box
    ax.add_collection(PatchCollection(patches, match_original=True))

    # add text
    for text_x, text_y, text_info, text_fontsize in text_list:
        if " to " in str(text_info):
            plt.text(
                x=text_x+0.3,
                y=text_y,
                s=text_info,
                horizontalalignment='right',
                verticalalignment='center',
                fontsize=text_fontsize,
                fontname="DejaVu Sans"
            )
        else:
            plt.text(
            x=text_x,
            y=text_y,
            s=text_info,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=text_fontsize,
            fontname="DejaVu Sans"
        )

    # output plot
    if str(ARGS.out_figure_format) == "PNG":
        fig.savefig(fname=ARGS.out_figure, bbox_inches='tight', dpi=int(ARGS.out_figure_dpi), format="png")

    else:
        fig.savefig(fname=ARGS.out_figure, bbox_inches='tight', dpi=int(ARGS.out_figure_dpi), format="pdf")
# plt.show()


