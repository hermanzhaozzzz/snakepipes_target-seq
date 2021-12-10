#! /home/zhaohuanan/anaconda3/envs/py27/bin/python
# _*_ coding: UTF-8 _*_


# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: MENG Howard

    Version-01:
      2019-10-10 demultiplex fastq using primer

    E-Mail: meng_howard@126.com
    """
# Version information END ----------------------------------------------------

# Learning Part START --------------------------------------------------------
LEARNING_PART = \
    """
        Input R1 and R2 FASTQ files and demultiplex those target-seq FASTQ files by primer
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


def get_primer_region_index(query_seq, primer_dict, fastq_state="read1", identify_mismatch_num=1, identify_match_num=9):
    """
    <HELP>
        back reads belong to which primer_dict key

    <INPUT>
        primer_dict, key= primer_name, value= pirmer_obj

    <RETURN>
        primer_dict key OR None
    """

    query_seq = query_seq.upper()
    primer_key = None
    primer_seq = None

    for primer_key in primer_dict.keys():
        if fastq_state == "read1":
            primer_seq = primer_dict[primer_key].upstream_primer

        elif fastq_state == "read2":
            primer_seq = primer_dict[primer_key].downstream_primer

        else:
            raise IOError("fastq_state should be read1 or read2")

        mismatch_count = 0
        match_count = 0
        for index in range(min(len(query_seq), len(primer_seq))):
            if query_seq[index] != primer_seq[index]:
                mismatch_count += 1
            else:
                match_count += 1

            if mismatch_count > identify_mismatch_num:
                break

        if (mismatch_count <= identify_mismatch_num) and (match_count >= identify_match_num):
            return (primer_key, mismatch_count, match_count)

    return (None)

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="demultiplex target-seq FASTQ files")

    parser.add_argument("-1", "--FastqRead1",
                        help="Input read1 FASTQ file, support .gz suffix", required=True)

    parser.add_argument("-2", "--FastqRead2",
                        help="Input read2 FASTQ file, support .gz suffix", required=True)

    parser.add_argument("-p", "--PrimerInfoFile",
                        help="Primer info file the format like: \
         <chr_name> <primer_region_start> <primer_region_end> \
         <region_index> <left_barcode_length> <right_barcode_length> \
         <left_primer> <right_primer> <target_strand> TAB split", required=True)

    parser.add_argument("-o", "--OutputDir",
                        help="Output dir, default=Input FASTQ file dir", default=None)

    parser.add_argument("--Read1BarcodeLen",
                        help="Read1 barcode length, default=4", default="4")

    parser.add_argument("--Read2BarcodeLen",
                        help="Read2 barcode length, default=6", default="6")

    parser.add_argument("--IdentifyLength",
                        help="Identify length after barcode sequence, default=10", default="10")

    parser.add_argument("--IdentifyMismatchCutoff",
                        help="Identify length after barcode sequence mismatch cutoff default=1", default="1")

    parser.add_argument("--MaxReadNumLoadOneTime",
                        help="The max number of reads load into into memory each time, default=1000000",
                        default="1000000")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    full_cmd = """python target_seq_demultiplex_fastq_V01.py \\
    -1 {read1_file} \\
    -2 {read2_file} \\
    -p {primer_file} \\
    -o {output_dir} \\
    --Read1BarcodeLen {Read1BarcodeLen} \\
    --Read2BarcodeLen {Read2BarcodeLen} \\
    --IdentifyLength {IdentifyLength} \\
    --IdentifyMismatchCutoff {IdentifyMismatchCutoff} \\
    --MaxReadNumLoadOneTime {MaxReadNumLoadOneTime}
    """.format(
        read1_file=ARGS.FastqRead1,
        read2_file=ARGS.FastqRead2,
        primer_file=ARGS.PrimerInfoFile,
        output_dir=ARGS.OutputDir,
        Read1BarcodeLen=ARGS.Read1BarcodeLen,
        Read2BarcodeLen=ARGS.Read2BarcodeLen,
        IdentifyLength=ARGS.IdentifyLength,
        IdentifyMismatchCutoff=ARGS.IdentifyMismatchCutoff,
        MaxReadNumLoadOneTime=ARGS.MaxReadNumLoadOneTime
    )

    # log
    sys.stderr.write("\n" + "-" * 80 + "\n")
    sys.stderr.write(full_cmd + "\n")
    sys.stderr.write("-" * 80 + "\n")

    # ------------------------------------->>>>>
    # file part
    # ------------------------------------->>>>>
    in_fq_1_filename = ARGS.FastqRead1
    in_fq_2_filename = ARGS.FastqRead2

    base_dir = None
    if ARGS.OutputDir == None:
        base_dir = os.path.abspath(os.path.dirname(in_fq_1_filename))
    else:
        base_dir = os.path.abspath(ARGS.OutputDir)
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
    if not os.path.exists(base_dir+'/demultiplex.fastq'):
        os.makedirs(base_dir+'/demultiplex.fastq')

    primer_info_filepath = ARGS.PrimerInfoFile

    # ------------------------------------->>>>>
    # setting and cutoff
    # ------------------------------------->>>>>
    read_1_barcode_len = int(ARGS.Read1BarcodeLen)
    read_2_barcode_len = int(ARGS.Read2BarcodeLen)
    identify_length_after_barcode = int(ARGS.IdentifyLength)
    identify_mismatch_num = int(ARGS.IdentifyMismatchCutoff)
    identify_match_num = identify_length_after_barcode - identify_mismatch_num - 1
    max_read_load_each_time = int(ARGS.MaxReadNumLoadOneTime)
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

    # -------------------------------------------------------------------->>>>>
    # load all FASTQ reads into dict
    # -------------------------------------------------------------------->>>>>
    # log
    sys.stderr.write(
        "Load all FASTQ file into memory ... %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # open FASTQ file
    in_fq_1 = None
    in_fq_2 = None

    if in_fq_1_filename[-3:] == ".gz":
        in_fq_1 = gzip.open(in_fq_1_filename, "r")
    else:
        in_fq_1 = open(in_fq_1_filename, "r")

    if in_fq_2_filename[-3:] == ".gz":
        in_fq_2 = gzip.open(in_fq_2_filename, "r")
    else:
        in_fq_2 = open(in_fq_2_filename, "r")

    # load all FASTQ file into memory
    temp_read_dict = {}
    found_primer_dict = {}
    line_count = 0

    while True:
        primer_found_flag = False
        # load fastq file
        fq_1_list = []
        fq_2_list = []
        for index in range(4):
            line_count += 1
            fq_1_line = in_fq_1.readline()
            fq_2_line = in_fq_2.readline()

            if fq_1_line and fq_2_line:
                fq_1_list.append(fq_1_line.strip())
                fq_2_list.append(fq_2_line.strip())
            else:
                break

        # log part START
        # Write temp reads to file and clear memory
        if (line_count // 4) % max_read_load_each_time == 0:
            sys.stderr.write("Load FASTQ reads count:\t%d\t%s\n" % (
            line_count // 4, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
            for primer_key in temp_read_dict.keys():
                # make output demultiplex filename
                out_fq_demultiplex_R1_filebase = 'demultiplex.fastq/{base_name}_demultiplex_R1.fastq'.format(base_name=primer_key)
                out_fq_demultiplex_R2_filebase = 'demultiplex.fastq/{base_name}_demultiplex_R2.fastq'.format(base_name=primer_key)

                out_fq_demultiplex_R1_filename = os.path.join(os.path.abspath(base_dir), out_fq_demultiplex_R1_filebase)
                out_fq_demultiplex_R2_filename = os.path.join(os.path.abspath(base_dir), out_fq_demultiplex_R2_filebase)

                out_fq_demultiplex_R1 = open(out_fq_demultiplex_R1_filename, "a")
                out_fq_demultiplex_R2 = open(out_fq_demultiplex_R2_filename, "a")
                for fq_obj in temp_read_dict[primer_key]['R1']:
                    out_fq_demultiplex_R1.write(fq_obj.write_format())
                for fq_obj in temp_read_dict[primer_key]['R2']:
                    out_fq_demultiplex_R2.write(fq_obj.write_format())
                out_fq_demultiplex_R1.close()
                out_fq_demultiplex_R2.close()
            temp_read_dict = {}
        # EOF
        if (len(fq_1_list) != 4) or (len(fq_2_list) != 4):
            break

        # init the fastq obj
        fq_1_obj = FastqRead(fq_1_list, phred=33)
        fq_2_obj = FastqRead(fq_2_list, phred=33)

        if fq_1_obj.read_id != fq_2_obj.read_id:
            raise IOError("read1 file not match with read2 file!")

        # get adjacent sequence
        fq_1_adjseq = fq_1_obj.sequence[read_1_barcode_len:read_1_barcode_len + identify_length_after_barcode]
        fq_2_adjseq = fq_2_obj.sequence[read_2_barcode_len:read_2_barcode_len + identify_length_after_barcode]

        fq_1_primer_key = get_primer_region_index(
            query_seq=fq_1_adjseq,
            primer_dict=primer_dict,
            fastq_state="read1",
            identify_mismatch_num=identify_mismatch_num,
            identify_match_num=identify_match_num)

        fq_2_primer_key = get_primer_region_index(
            query_seq=fq_2_adjseq,
            primer_dict=primer_dict,
            fastq_state="read2",
            identify_mismatch_num=identify_mismatch_num,
            identify_match_num=identify_match_num)

        # demultiplex fastq using primer:
        if (fq_1_primer_key != None) and (fq_2_primer_key != None):
            if fq_1_primer_key == fq_2_primer_key:
                primer_found_flag = True
                # ------------------------------------------------->>>>>>>>
                # make new dict
                # ------------------------------------------------->>>>>>>>
                if temp_read_dict.get(fq_1_primer_key[0]) == None:
                    temp_read_dict[fq_1_primer_key[0]] = {'R1': [], 'R2': []}
                temp_read_dict[fq_1_primer_key[0]]['R1'].append(fq_1_obj)
                temp_read_dict[fq_1_primer_key[0]]['R2'].append(fq_2_obj)
        #change barcode length and try again
        if not primer_found_flag:
            fq_1_adjseq = fq_1_obj.sequence[read_2_barcode_len:read_2_barcode_len + identify_length_after_barcode]
            fq_2_adjseq = fq_2_obj.sequence[read_1_barcode_len:read_1_barcode_len + identify_length_after_barcode]

            fq_1_primer_key = get_primer_region_index(
                query_seq=fq_1_adjseq,
                primer_dict=primer_dict,
                fastq_state="read1",
                identify_mismatch_num=identify_mismatch_num,
                identify_match_num=identify_match_num)

            fq_2_primer_key = get_primer_region_index(
                query_seq=fq_2_adjseq,
                primer_dict=primer_dict,
                fastq_state="read2",
                identify_mismatch_num=identify_mismatch_num,
                identify_match_num=identify_match_num)

            # demultiplex fastq using primer:
            if (fq_1_primer_key != None) and (fq_2_primer_key != None):
                if fq_1_primer_key == fq_2_primer_key:
                    Flag = True
                    # ------------------------------------------------->>>>>>>>
                    # make new dict
                    # ------------------------------------------------->>>>>>>>
                    if temp_read_dict.get(fq_1_primer_key[0]) == None:
                        temp_read_dict[fq_1_primer_key[0]] = {'R1': [], 'R2': []}
                    temp_read_dict[fq_1_primer_key[0]]['R1'].append(fq_1_obj)
                    temp_read_dict[fq_1_primer_key[0]]['R2'].append(fq_2_obj)

    # Write remaining reads to file and clear memory
    for primer_key in temp_read_dict.keys():
        # make output demultiplex filename
        out_fq_demultiplex_R1_filebase = 'demultiplex.fastq/{base_name}_demultiplex_R1.fastq'.format(base_name=primer_key)
        out_fq_demultiplex_R2_filebase = 'demultiplex.fastq/{base_name}_demultiplex_R2.fastq'.format(base_name=primer_key)

        out_fq_demultiplex_R1_filename = os.path.join(os.path.abspath(base_dir), out_fq_demultiplex_R1_filebase)
        out_fq_demultiplex_R2_filename = os.path.join(os.path.abspath(base_dir), out_fq_demultiplex_R2_filebase)

        out_fq_demultiplex_R1 = open(out_fq_demultiplex_R1_filename, "a")
        out_fq_demultiplex_R2 = open(out_fq_demultiplex_R2_filename, "a")
        for fq_obj in temp_read_dict[primer_key]['R1']:
            out_fq_demultiplex_R1.write(fq_obj.write_format())
        for fq_obj in temp_read_dict[primer_key]['R2']:
            out_fq_demultiplex_R2.write(fq_obj.write_format())
        out_fq_demultiplex_R1.close()
        out_fq_demultiplex_R2.close()
    temp_read_dict = {}

    in_fq_1.close()
    in_fq_2.close()
    # log
    sys.stderr.write("Load FASTQ reads Done! Total reads count:\t%d\t%s\n" % (
    line_count // 4, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))


# 2019-10-10