#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_



# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: MENG Howard

    Version-01:
      2019-08-19 merge barcode contained target-seq fastq

    Version-02:
      2019-09-23 demultiplex fastq using primer and than merge by barcode

    Version-03:
      2019-10-10 merge barcode using demultiplex fastq as input file

    E-Mail: meng_howard@126.com
    """
# Version information END ----------------------------------------------------

# Learning Part START --------------------------------------------------------
LEARNING_PART = \
    """
        Input R1 and R2 FASTQ files and merge those target-seq FASTQ files
    """
# Learning Part END-----------------------------------------------------------

import argparse
import gzip
import sys
import os
import string
import time
import math


###############################################################################
# function part
###############################################################################
class FastqRead(object):
    '''
        Define fastq read object
    '''

    def __init__(self, fastq_read_list, phred=33):
        self.head = fastq_read_list[0]
        self.sequence = fastq_read_list[1]
        self.info = fastq_read_list[2]
        self.quality = fastq_read_list[3]
        self.phred = phred
        self.length = len(self.sequence)
        read_id_raw = fastq_read_list[0]
        
        if ':' in read_id_raw:
            # Illumina
            self.read_id = fastq_read_list[0].split(" ")[0]
        else:
            # MGI reads
            self.read_id = fastq_read_list[0].split("/")[0]
        # elif read_id_raw[-2] == '/':
        #     # MGI reads
        #     self.read_id = fastq_read_list[0].split("/")[0]
        # else:
            # raise OError("read1 file not match with read2 file! (Or parse read id failed, MGI or Illumina?)")

    def trimmer(self, start=0, end=None):
        # 返回trim的fastq read对象
        if (start < 0) or (start >= self.length):
            start = 0

        if (end >= self.length) or (end <= start):
            end = None

        trim_sequence = self.sequence[start:end]
        trim_quality = self.quality[start:end]
        trim_info = self.info
        trim_read = [self.head, trim_sequence, trim_info, trim_quality]
        return (FastqRead(trim_read, phred=self.phred))

    def get_phred(self):
        self.phred_score = [ord(x) - self.phred for x in self.quality]
        return (self.phred_score)

    def write_format(self):
        # 返回一个str可以用于直接写入文件
        write_str = "{0}\n{1}\n{2}\n{3}\n".format(
            self.head,
            self.sequence,
            self.info,
            self.quality)

        return (write_str)


class PrimerInfo(object):
    """
    For target-seq primer info
    """

    def __init__(self, primer_line_list):
        self.chr_name = primer_line_list[0]
        self.region_start = int(primer_line_list[1])
        self.region_end = int(primer_line_list[2])
        self.region_id = primer_line_list[3]
        self.left_barcode_length = int(primer_line_list[4])
        self.right_barcode_length = int(primer_line_list[5])
        self.upstream_primer = primer_line_list[6].upper()
        self.downstream_primer = primer_line_list[7].upper()
        self.target_strand = primer_line_list[8]

def back_final_seq(temp_read_obj_list, min_read_count_cutoff, min_base_ratio_cutoff):
    """
    <HELP>
        min_read_count_cutoff, the reads with same barcode lower than this number will filtered
        min_base_ratio_cutoff, if base ratio lower than min_base_ratio_cutoff, will filtered

    <RETURN>
        [final_seq, final_phred] OR None
    """
    # set cutoff
    min_base_count_cutoff = int(math.ceil(min_read_count_cutoff * min_base_ratio_cutoff))
    median_index = len(temp_read_obj_list) // 2

    if len(temp_read_obj_list) >= min_read_count_cutoff:
        read_seq_mat = zip(*[fq_obj.sequence for fq_obj in temp_read_obj_list])
        read_phred_mat = zip(*[fq_obj.quality for fq_obj in temp_read_obj_list])

        read_final_seq = ""
        read_final_phred = ""

        for index, base_list in enumerate(read_seq_mat):
            test_base = base_list[0]
            pass_min_base_count_cutoff = False

            if base_list.count(test_base) >= min_base_count_cutoff:
                read_final_seq += test_base
                read_final_phred += sorted(read_phred_mat[index])[median_index]
                pass_min_base_count_cutoff = True

            else:
                for test_base in ["A", "T", "C", "G", "N"]:
                    if test_base != base_list[0]:
                        if base_list.count(test_base) >= min_base_count_cutoff:
                            read_final_seq += test_base
                            read_final_phred += sorted(read_phred_mat[index])[median_index]
                            pass_min_base_count_cutoff = True

            if not pass_min_base_count_cutoff:
                return (None)

        # make back fastq list
        return ([read_final_seq, read_final_phred])

    else:
        return (None)


###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="convert mpileup file to info file")

    parser.add_argument("-i", "--InputDir",
                        help="Input dir", required=True)

    parser.add_argument("-p", "--PrimerInfoFile",
                        help="Primer info file the format like: \
         <chr_name> <primer_region_start> <primer_region_end> \
         <region_index> <left_barcode_length> <right_barcode_length> \
         <left_primer> <right_primer> <target_strand> TAB split", required=True)

    parser.add_argument("-o", "--OutputDir",
                        help="Output dir, default=Input dir", default=None)

    parser.add_argument("-r", "--ReportFile",
                        help="Report file of processing information, default=stdout", default="stdout")

    parser.add_argument("--Read1BarcodeLen",
                        help="Read1 barcode length, default=4", default="4")

    parser.add_argument("--Read2BarcodeLen",
                        help="Read2 barcode length, default=6", default="6")

    parser.add_argument("--MinMergeReadNumCutoff",
                        help="A group of reads have to larger than this cutoff, default=3", default="3")

    parser.add_argument("--MinMergeReadRatioCutoff",
                        help="A group of reads' identity have to larger than this cutoff, default=0.7", default="0.7")

    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    full_cmd = """python target_seq_merge_fq_V03.py \\
    -i {input_dir} \\
    -p {primer_file} \\
    -o {output_dir} \\
    -r {report_file} \\
    --Read1BarcodeLen {Read1BarcodeLen} \\
    --Read2BarcodeLen {Read2BarcodeLen} \\
    --MinMergeReadNumCutoff {MinMergeReadNumCutoff} \\
    --MinMergeReadRatioCutoff {MinMergeReadRatioCutoff}\\
    """.format(
        input_dir=ARGS.InputDir,
        primer_file=ARGS.PrimerInfoFile,
        output_dir=ARGS.OutputDir,
        report_file=ARGS.ReportFile,
        Read1BarcodeLen=ARGS.Read1BarcodeLen,
        Read2BarcodeLen=ARGS.Read2BarcodeLen,
        MinMergeReadNumCutoff=ARGS.MinMergeReadNumCutoff,
        MinMergeReadRatioCutoff=ARGS.MinMergeReadRatioCutoff
    )

    # log
    sys.stderr.write("\n" + "-" * 80 + "\n")
    sys.stderr.write(full_cmd + "\n")
    sys.stderr.write("-" * 80 + "\n")

    # ------------------------------------->>>>>
    # file part
    # ------------------------------------->>>>>
    base_dir = ARGS.InputDir
    if ARGS.OutputDir == None:
        out_base_dir = os.path.abspath(ARGS.InputDir)
    else:
        out_base_dir = os.path.abspath(ARGS.OutputDir)
        if not os.path.exists(out_base_dir):
            os.makedirs(out_base_dir)
    if not os.path.exists(out_base_dir+'/merge.fastq'):
        os.makedirs(out_base_dir+'/merge.fastq')

    primer_info_filepath = ARGS.PrimerInfoFile

    report_file = None
    if ARGS.ReportFile == "stdout":
        report_file = sys.stdout
    else:
        report_file = open(ARGS.ReportFile, "w")

    # ------------------------------------->>>>>
    # setting and cutoff
    # ------------------------------------->>>>>
    read_1_barcode_len = int(ARGS.Read1BarcodeLen)
    read_2_barcode_len = int(ARGS.Read2BarcodeLen)
    min_read_count_cutoff = int(ARGS.MinMergeReadNumCutoff)
    min_base_ratio_cutoff = float(ARGS.MinMergeReadRatioCutoff)
    # -------------------------------------------------------------------->>>>>
    # load the primer info
    # -------------------------------------------------------------------->>>>>
    # log
    sys.stderr.write("Load primer info... %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    primer_dict = {}
    with open(primer_info_filepath, "r") as primer_in_file:
        header = primer_in_file.readline()
        for line in primer_in_file:
            line_list = line.strip().split("\t")
            primer_obj = PrimerInfo(line_list)
            if (primer_obj.right_barcode_length > 0) and (primer_obj.left_barcode_length > 0):
                if not primer_dict.get(primer_obj.region_id):
                    primer_dict[primer_obj.region_id] = primer_obj

    # log
    sys.stderr.write("Load primer info... Done! %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # output report file header
    report_file_header = "case_name\tR1_barcode\tR2_barcode\traw_fq_count\tif_final_output"
    report_file.write(report_file_header + "\n")

    # output
    all_out_count = 0

    # load all reads for each primer into memory
    for primer_key in primer_dict:
        primer_temp_read_dict = {}
        line_count = 0
        final_out_count = 0
        in_fq_demultiplex_R1_filebase = 'demultiplex.fastq/{base_name}_demultiplex_R1.fastq'.format(base_name=primer_key)
        in_fq_demultiplex_R2_filebase = 'demultiplex.fastq/{base_name}_demultiplex_R2.fastq'.format(base_name=primer_key)

        in_fq_demultiplex_R1_filename = os.path.join(os.path.abspath(base_dir), in_fq_demultiplex_R1_filebase)
        in_fq_demultiplex_R2_filename = os.path.join(os.path.abspath(base_dir), in_fq_demultiplex_R2_filebase)

        if in_fq_demultiplex_R1_filename[-3:] == ".gz":
            try:
                in_fq_demultiplex_1 = gzip.open(in_fq_demultiplex_R1_filename, "r")
            except IOError:
                continue
        else:
            try:
                in_fq_demultiplex_1 = open(in_fq_demultiplex_R1_filename, "r")
            except IOError:
                continue
                
        if in_fq_demultiplex_R2_filename[-3:] == ".gz":
            in_fq_demultiplex_2 = gzip.open(in_fq_demultiplex_R2_filename, "r")
        else:
            in_fq_demultiplex_2 = open(in_fq_demultiplex_R2_filename, "r")

        while True:
            # load fastq file
            fq_1_list = []
            fq_2_list = []
            for index in range(4):
                line_count += 1
                fq_1_line = in_fq_demultiplex_1.readline()
                fq_2_line = in_fq_demultiplex_2.readline()

                if fq_1_line and fq_2_line:
                    fq_1_list.append(fq_1_line.strip())
                    fq_2_list.append(fq_2_line.strip())
                else:
                    break
            # EOF
            if (len(fq_1_list) != 4) or (len(fq_2_list) != 4):
                break

            # init the fastq obj
            fq_1_obj = FastqRead(fq_1_list, phred=33)
            fq_2_obj = FastqRead(fq_2_list, phred=33)

            if fq_1_obj.read_id != fq_2_obj.read_id:
                # DBUG
                sys.stderr.write("error read pair: r1: %s, r2: %s" % (fq_1_obj.read_id, fq_2_obj.read_id))
                raise IOError("read1 file not match with read2 file!")

            # get barcode sequence
            fq_1_barcode = fq_1_obj.sequence[:primer_dict[primer_key].left_barcode_length]
            fq_2_barcode = fq_2_obj.sequence[:primer_dict[primer_key].right_barcode_length]

            temp_read_barcode_key = "{left_barcode}.{right_barcode}".format(
                left_barcode=fq_1_barcode,
                right_barcode=fq_2_barcode
            )

            # ------------------------------------------------->>>>>>>>
            # no query key then get new file
            # ------------------------------------------------->>>>>>>>
            if primer_temp_read_dict.get(temp_read_barcode_key) == None:
                primer_temp_read_dict[temp_read_barcode_key] = [[], [], 0]

            # ------------------------------------------------->>>>>>>>
            # write output fastq file
            # ------------------------------------------------->>>>>>>>
            fq_1_obj_rm_barcode = fq_1_obj.trimmer(primer_dict[primer_key].left_barcode_length)
            fq_2_obj_rm_barcode = fq_2_obj.trimmer(primer_dict[primer_key].right_barcode_length)

            fq_1_obj_rm_barcode.head = fq_1_obj_rm_barcode.read_id + ":" + fq_1_obj.sequence[
                                                                           :primer_dict[primer_key].left_barcode_length] + ":" + fq_2_obj.sequence[
                                                                                                        :primer_dict[primer_key].right_barcode_length] + " " + \
                                       fq_1_obj_rm_barcode.head.split(" ")[1]
            fq_2_obj_rm_barcode.head = fq_1_obj_rm_barcode.read_id + ":" + fq_1_obj.sequence[
                                                                           :primer_dict[primer_key].left_barcode_length] + ":" + fq_2_obj.sequence[
                                                                                                        :primer_dict[primer_key].right_barcode_length] + " " + \
                                       fq_2_obj_rm_barcode.head.split(" ")[1]

            primer_temp_read_dict[temp_read_barcode_key][0].append(fq_1_obj_rm_barcode)
            primer_temp_read_dict[temp_read_barcode_key][1].append(fq_2_obj_rm_barcode)
            # count
            primer_temp_read_dict[temp_read_barcode_key][2] += 1
        # log
        sys.stderr.write("Load FASTQ reads for primer\t%s\tDone! Total reads count:\t%d\t%s\n" % (
        primer_key, line_count // 4, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        # make output filename
        out_fq_R1_filebase = 'merge.fastq/{base_name}_merge_barcode_R1.fastq'.format(base_name=primer_key)
        out_fq_R2_filebase = 'merge.fastq/{base_name}_merge_barcode_R2.fastq'.format(base_name=primer_key)

        out_fq_R1_filename = os.path.join(os.path.abspath(out_base_dir), out_fq_R1_filebase)
        out_fq_R2_filename = os.path.join(os.path.abspath(out_base_dir), out_fq_R2_filebase)

        out_fq_R1 = open(out_fq_R1_filename, "w")
        out_fq_R2 = open(out_fq_R2_filename, "w")

        for barcode_info in primer_temp_read_dict.keys():
            if primer_temp_read_dict[barcode_info][2] >= min_read_count_cutoff:
                temp_read_list = primer_temp_read_dict[barcode_info]
                final_report_state = False

                # make output fastq str list
                merge_R1_info = back_final_seq(temp_read_list[0], min_read_count_cutoff, min_base_ratio_cutoff)
                merge_R2_info = back_final_seq(temp_read_list[1], min_read_count_cutoff, min_base_ratio_cutoff)

                if (merge_R1_info != None) and (merge_R2_info != None):
                    final_out_count += 1
                    final_report_state = True

                    if final_out_count % 10000 == 0:
                        # log
                        sys.stderr.write("Merge FASTQ reads count for primer\t%s\t:\t%d\t%s\n" % (
                            primer_key, final_out_count, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

                    out_head = "@{base_name}:{R1_R2_barcode}:{raw_count}".format(
                        base_name=primer_key,
                        R1_R2_barcode=barcode_info.replace(".", ":"),
                        raw_count=len(temp_read_list[0])
                    )

                    out_head_R1 = out_head + " 1"
                    out_head_R2 = out_head + " 2"

                    out_R1_list = [out_head_R1, merge_R1_info[0], "+", merge_R1_info[1]]
                    out_R2_list = [out_head_R2, merge_R2_info[0], "+", merge_R2_info[1]]
                    out_fq_R1.write("\n".join(map(str, out_R1_list)) + "\n")
                    out_fq_R2.write("\n".join(map(str, out_R2_list)) + "\n")

                    # report region
                report_list = [
                    primer_key,
                    barcode_info.split(".")[0],
                    barcode_info.split(".")[1],
                    len(temp_read_list[0]),
                    final_report_state
                ]

                report_file.write("\t".join(map(str, report_list)) + "\n")
        out_fq_R1.close()
        out_fq_R2.close()
        all_out_count += final_out_count
        # log
        sys.stderr.write("Merge FASTQ reads for primer\t%s\tdone! Final out reads count:\t%d\t%s\n" % (
            primer_key, final_out_count, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        # free of memory
        del primer_temp_read_dict

    # log
    sys.stderr.write("Merge FASTQ done! Final out reads count:\t%d\t%s\n" % (
    all_out_count, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

# 2019-10-10