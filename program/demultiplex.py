import gzip
import os
import sys
import time

import fire


class FastqRead(object):
    """Define fastq read object"""

    def __init__(self, fastq_read_list, phred=33):
        self.head = fastq_read_list[0]
        self.sequence = fastq_read_list[1]
        self.info = fastq_read_list[2]
        self.quality = fastq_read_list[3]
        self.phred = phred
        self.length = len(self.sequence)
        read_id_raw = fastq_read_list[0]

        if ":" in read_id_raw:
            # Illumina
            self.read_id = fastq_read_list[0].split(" ")[0]
        else:
            # MGI reads
            self.read_id = fastq_read_list[0].split("/")[0]

    def trimmer(self, start=0, end=None):
        if (start < 0) or (start >= self.length):
            start = 0
        if (end is None) or (end >= self.length) or (end <= start):
            end = None
        trim_sequence = self.sequence[start:end]
        trim_quality = self.quality[start:end]
        trim_info = self.info
        trim_read = [self.head, trim_sequence, trim_info, trim_quality]
        return FastqRead(trim_read, phred=self.phred)

    def get_phred(self):
        self.phred_score = [ord(x) - self.phred for x in self.quality]
        return self.phred_score

    def write_format(self):
        write_str = "{0}\n{1}\n{2}\n{3}\n".format(
            self.head, self.sequence, self.info, self.quality
        )
        return write_str


class PrimerInfo(object):
    """For target-seq primer info"""

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


def main(
    fastq_read1: str,
    fastq_read2: str,
    primer_info_file: str,
    output_dir: str = None,
    read1_barcode_len: int = 4,
    read2_barcode_len: int = 6,
    identify_length: int = 10,
    identify_mismatch_cutoff: int = 1,
    max_read_num_load: int = 1000000,
):
    """
    Target-seq数据分选主函数

    Args:
        fastq_read1:       输入read1 FASTQ路径（支持.gz）
        fastq_read2:       输入read2 FASTQ路径（支持.gz）
        primer_info_file:  引物信息文件路径
        output_dir:        输出目录（默认输入文件所在目录）
        read1_barcode_len: Read1 Barcode长度（默认4）
        read2_barcode_len: Read2 Barcode长度（默认6）
        identify_length:   引物识别区域长度（默认10）
        identify_mismatch_cutoff: 允许的最大错配数（默认1）
        max_read_num_load: 单次加载最大reads数（默认1000000）
    """
    # 参数预处理
    identify_match_num = identify_length - identify_mismatch_cutoff - 1
    base_dir = _prepare_output_dir(output_dir, fastq_read1)

    # 加载引物信息
    primer_dict = load_primer_info(primer_info_file)

    # 处理Fastq文件
    process_fastq_files(
        fastq_read1,
        fastq_read2,
        primer_dict,
        base_dir,
        read1_barcode_len,
        read2_barcode_len,
        identify_length,
        identify_mismatch_cutoff,
        identify_match_num,
        max_read_num_load,
    )


def _prepare_output_dir(output_dir, input_path):
    """处理输出目录"""
    if not output_dir:
        base_dir = os.path.abspath(os.path.dirname(input_path))
    else:
        base_dir = os.path.abspath(output_dir)

    os.makedirs(os.path.join(base_dir, "demultiplex.fastq"), exist_ok=True)
    return base_dir


def load_primer_info(primer_file):
    """加载引物信息文件"""
    primer_dict = {}
    with open(primer_file, "r") as f:
        next(f)  # 跳过标题行
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 9:
                primer = PrimerInfo(parts)
                primer_dict[primer.region_id] = primer
    return primer_dict


def process_fastq_files(
    r1_path,
    r2_path,
    primer_dict,
    base_dir,
    r1_barcode_len,
    r2_barcode_len,
    identify_len,
    mismatch_cutoff,
    match_num,
    max_reads,
):
    """核心处理逻辑"""
    # 打开文件
    with (
        gzip.open(r1_path, "rt") if r1_path.endswith(".gz") else open(r1_path) as fq1,
        gzip.open(r2_path, "rt") if r2_path.endswith(".gz") else open(r2_path) as fq2,
    ):
        temp_data = {}
        line_count = 0

        while True:
            # 读取批次数据
            batch = read_batch(fq1, fq2, max_reads)
            if not batch:
                break

            # 处理每个read pair
            for fq1_obj, fq2_obj in batch:
                line_count += 1
                process_read_pair(
                    fq1_obj,
                    fq2_obj,
                    temp_data,
                    r1_barcode_len,
                    r2_barcode_len,
                    identify_len,
                    mismatch_cutoff,
                    match_num,
                    primer_dict,
                )

            # 定期保存数据
            save_data(temp_data, base_dir)
            log_progress(line_count)

        # 保存剩余数据
        save_data(temp_data, base_dir)
        log_final(line_count)


def read_batch(fq1, fq2, max_reads):
    """读取一批reads"""
    batch = []
    for _ in range(max_reads):
        r1_lines = [fq1.readline().strip() for _ in range(4)]
        r2_lines = [fq2.readline().strip() for _ in range(4)]

        if not all(r1_lines) or not all(r2_lines):
            break

        batch.append((FastqRead(r1_lines), FastqRead(r2_lines)))
    return batch


def get_primer_region_index(
    query_seq,
    primer_dict,
    fastq_state="read1",
    identify_mismatch_num=1,
    identify_match_num=9,
):
    """
    <HELP>
        back reads belong to which primer_dict key
    <INPUT>
        primer_dict, key= primer_name, value= pirmer_obj
    <RETURN>
        primer_dict key OR None
    """
    query_seq = query_seq.upper()
    for primer_key in primer_dict:
        if fastq_state == "read1":
            primer_seq = primer_dict[primer_key].upstream_primer
        elif fastq_state == "read2":
            primer_seq = primer_dict[primer_key].downstream_primer
        else:
            raise ValueError("fastq_state should be read1 or read2")
        mismatch_count = 0
        match_count = 0
        for index in range(min(len(query_seq), len(primer_seq))):
            if query_seq[index] != primer_seq[index]:
                mismatch_count += 1
            else:
                match_count += 1
            if mismatch_count > identify_mismatch_num:
                break
        if (mismatch_count <= identify_mismatch_num) and (
            match_count >= identify_match_num
        ):
            return (primer_key, mismatch_count, match_count)
    return (None, 0, 0)


def process_read_pair(
    r1, r2, temp_data, r1_len, r2_len, ident_len, mismatch, match_num, primers
):
    """处理单个read pair"""
    # 提取识别序列
    r1_seq = r1.sequence[r1_len : r1_len + ident_len]
    r2_seq = r2.sequence[r2_len : r2_len + ident_len]

    # 引物识别
    r1_primer = get_primer_region_index(r1_seq, primers, "read1", mismatch, match_num)
    r2_primer = get_primer_region_index(r2_seq, primers, "read2", mismatch, match_num)

    # 有效分选
    if r1_primer[0] and r2_primer[0] and (r1_primer[0] == r2_primer[0]):
        primer_id = r1_primer[0]
        if primer_id not in temp_data:
            temp_data[primer_id] = {"R1": [], "R2": []}
        temp_data[primer_id]["R1"].append(r1)
        temp_data[primer_id]["R2"].append(r2)


def save_data(data_dict, base_dir):
    """保存分选数据"""
    for pid, data in data_dict.items():
        r1_path = os.path.join(base_dir, f"demultiplex.fastq/{pid}_R1.fastq")
        r2_path = os.path.join(base_dir, f"demultiplex.fastq/{pid}_R2.fastq")

        with open(r1_path, "a") as f1, open(r2_path, "a") as f2:
            for read in data["R1"]:
                f1.write(read.write_format())
            for read in data["R2"]:
                f2.write(read.write_format())

    data_dict.clear()


def log_progress(count):
    """记录处理进度"""
    sys.stderr.write(f"\rProcessed {count} read pairs...")
    sys.stderr.flush()


def log_final(total):
    """最终日志"""
    sys.stderr.write(f"\nTotal processed: {total} read pairs\n")
    sys.stderr.write(f"Completed at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")


if __name__ == "__main__":
    fire.Fire(main)
