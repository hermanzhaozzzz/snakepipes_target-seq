#!/usr/bin/env python3
"""
TargetSeq FASTQ 合并工具 (单区域版)
功能：根据引物信息去除barcode，按指定阈值合并相同barcode的reads

1. cutoff=0
    功能：仅去除barcode，不进行合并
    处理流程：
    直接修剪R1左侧和R2右侧的barcode
    所有原始reads会被保留（仅修剪后的形式）
    输出结果：
    output_r1 和 output_r2 中的reads数量与输入文件相同
    报告文件 (report_path) 为空（无统计信息）
    适用场景：快速获取去除barcode的原始数据，不关心合并
2. cutoff=3
    功能：合并相同barcode的reads，保留出现次数≥3的barcode
    处理流程：
    收集阶段：将所有reads按barcode分组
    合并阶段：对每个barcode组：
    若组内reads数≥3，合并生成一条代表该组的read
    否则丢弃该组所有reads
    输出结果：
    输出文件仅包含满足条件的barcode合并后的reads
    报告文件列出所有barcode及其计数，标记是否通过 (Passed=Y/N)
    适用场景：平衡数据量和准确性，过滤低频错误barcode
3. cutoff=5
    逻辑类似cutoff=3，但条件更严格
    效果：
    仅保留出现次数≥5的barcode组
    输出数据量比cutoff=3时更少
    合并后的reads可靠性更高（基于更多原始reads）
    适用场景：对数据准确性要求更高，允许丢失部分低频barcode
4. cutoff=10
    最严格条件，仅保留高频barcode
    效果：
    输出数据量大幅减少
    合并结果代表高度一致的序列（适合高精度需求）
    可能丢失真实存在的低频变异
    适用场景：研究高度保守区域或要求极低错误率
"""

import gzip
import os
from typing import List, Optional

import fire


class FastqRead:
    """FASTQ 记录处理类（保持原样）"""

    def __init__(self, lines: List[str], phred: int = 33):
        if len(lines) != 4:
            raise ValueError("Invalid FASTQ record")
        self.header = lines[0].strip()
        self.sequence = lines[1].strip()
        self.plus = lines[2].strip()
        self.quality = lines[3].strip()
        self.phred = phred
        if " " in self.header:
            self.read_id = self.header.split(" ")[0][1:]
        elif "/" in self.header:
            self.read_id = self.header.split("/")[0][1:]
        else:
            self.read_id = self.header[1:]

    def trim(self, start: int, end: Optional[int] = None) -> "FastqRead":
        new_seq = self.sequence[start:end]
        new_qual = self.quality[start:end]
        return FastqRead([f"@{self.read_id}", new_seq, self.plus, new_qual])

    def write(self) -> str:
        return f"{self.header}\n{self.sequence}\n{self.plus}\n{self.quality}\n"


class PrimerInfo:
    """引物信息处理类（保持原样）"""

    def __init__(self, fields: List[str]):
        if len(fields) < 9:
            raise ValueError("Invalid primer line")
        self.chr = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.region_id = fields[3]
        self.left_barcode_len = int(fields[4])
        self.right_barcode_len = int(fields[5])
        self.left_primer = fields[6].upper()
        self.right_primer = fields[7].upper()
        self.strand = fields[8]


# =============================================================================
# 核心逻辑
# =============================================================================
def process_single_region(
    r1_path: str,
    r2_path: str,
    primer_path: str,
    region_name: str,
    output_r1: str,
    output_r2: str,
    report_path: str,
    cutoff: int,
):
    """处理单个区域的主函数"""
    # 加载引物信息
    primer = load_primer_info(primer_path, region_name)

    # 创建输出目录
    os.makedirs(os.path.dirname(output_r1), exist_ok=True)
    os.makedirs(os.path.dirname(output_r2), exist_ok=True)

    # Cutoff为0时的快速处理模式
    if cutoff == 0:
        simple_trim(
            r1_path,
            r2_path,
            output_r1,
            output_r2,
            primer.left_barcode_len,
            primer.right_barcode_len,
        )
        with open(report_path, "w") as report:
            report.write("")
        return

    # 正常合并模式
    merge_reads_with_validation(
        r1_path,
        r2_path,
        output_r1,
        output_r2,
        primer.left_barcode_len,
        primer.right_barcode_len,
        report_path,
        cutoff,
    )


def load_primer_info(primer_path: str, region_name: str) -> PrimerInfo:
    """从文件加载指定区域的引物信息"""
    with open(primer_path) as f:
        header = f.readline()  # 跳过标题行
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[3] == region_name:
                return PrimerInfo(fields)
    raise ValueError(f"Region '{region_name}' not found in primer file")


def simple_trim(
    r1_in: str, r2_in: str, r1_out: str, r2_out: str, trim_left: int, trim_right: int
):
    """简单去除barcode模式（cutoff=0时使用）"""
    # 自动处理gzip压缩
    open_in = gzip.open if r1_in.endswith(".gz") else open
    mode_in = "rt" if r1_in.endswith(".gz") else "r"

    open_out = gzip.open if r1_out.endswith(".gz") else open
    mode_out = "wt" if r1_out.endswith(".gz") else "w"

    with (
        open_in(r1_in, mode_in) as f1,
        open_in(r2_in, mode_in) as f2,
        open_out(r1_out, mode_out) as o1,
        open_out(r2_out, mode_out) as o2,
    ):
        while True:
            # 读取R1
            r1_lines = [f1.readline() for _ in range(4)]
            if not all(r1_lines):
                break

            # 读取R2
            r2_lines = [f2.readline() for _ in range(4)]
            if not all(r2_lines):
                break

            # 处理R1
            read1 = FastqRead(r1_lines)
            trimmed_r1 = read1.trim(trim_left)
            o1.write(trimmed_r1.write())

            # 处理R2
            read2 = FastqRead(r2_lines)
            trimmed_r2 = read2.trim(trim_right)
            o2.write(trimmed_r2.write())


def merge_reads_with_validation(
    r1_in: str,
    r2_in: str,
    r1_out: str,
    r2_out: str,
    trim_left: int,
    trim_right: int,
    report_path: str,
    cutoff: int,
):
    """合并reads主逻辑（cutoff>0时使用）"""
    # 初始化数据结构
    barcode_groups = {}

    # 第一次遍历：收集相同barcode的reads
    open_in = gzip.open if r1_in.endswith(".gz") else open
    mode_in = "rt" if r1_in.endswith(".gz") else "r"

    with open_in(r1_in, mode_in) as f1, open_in(r2_in, mode_in) as f2:
        while True:
            r1_lines = [f1.readline() for _ in range(4)]
            if not all(r1_lines):
                break

            r2_lines = [f2.readline() for _ in range(4)]
            if not all(r2_lines):
                break

            # 处理R1
            read1 = FastqRead(r1_lines)
            trimmed_r1 = read1.trim(trim_left)
            barcode1 = read1.sequence[:trim_left]

            # 处理R2
            read2 = FastqRead(r2_lines)
            trimmed_r2 = read2.trim(trim_right)
            barcode2 = read2.sequence[:trim_right]

            # 合并barcode标识
            barcode_key = f"{barcode1}:{barcode2}"

            if barcode_key not in barcode_groups:
                barcode_groups[barcode_key] = {
                    "r1_reads": [],
                    "r2_reads": [],
                    "count": 0,
                }

            barcode_groups[barcode_key]["r1_reads"].append(trimmed_r1)
            barcode_groups[barcode_key]["r2_reads"].append(trimmed_r2)
            barcode_groups[barcode_key]["count"] += 1

    # 第二次遍历：合并并输出
    open_out = gzip.open if r1_out.endswith(".gz") else open
    mode_out = "wt" if r1_out.endswith(".gz") else "w"

    with (
        open_out(r1_out, mode_out) as o1,
        open_out(r2_out, mode_out) as o2,
        open(report_path, "w") as report,
    ):
        report.write("Barcode\tTotalReads\tPassed\n")

        for bc_key, group in barcode_groups.items():
            count = group["count"]
            passed = "Y" if count >= cutoff else "N"
            report.write(f"{bc_key}\t{count}\t{passed}\n")

            if count < cutoff:
                continue
            # 修改点3：同步处理R1/R2合并结果
            merged_r1 = merge_single_reads(group["r1_reads"])
            merged_r2 = merge_single_reads(group["r2_reads"])

            # 关键修改：仅当两者都成功时写入
            if merged_r1 and merged_r2:
                header = construct_merged_header(bc_key, count)
                o1.write(f"{header}/1\n{merged_r1[0]}\n+\n{merged_r1[1]}\n")
                o2.write(f"{header}/2\n{merged_r2[0]}\n+\n{merged_r2[1]}\n")


def merge_single_reads(reads: List[FastqRead], min_ratio=0.6) -> Optional[tuple]:
    """改进后的合并函数，移除了不必要的长度检查"""
    try:
        seq_matrix = [list(r.sequence) for r in reads]
        qual_matrix = [list(r.quality) for r in reads]

        merged_seq = []
        for pos in zip(*seq_matrix):
            bases = {}
            for b in pos:
                bases[b] = bases.get(b, 0) + 1
            for b in ["A", "T", "C", "G", "N"]:
                if bases.get(b, 0) / len(pos) >= min_ratio:
                    merged_seq.append(b)
                    break
            else:
                return None

        merged_qual = []
        for qpos in zip(*qual_matrix):
            merged_qual.append(sorted(qpos)[len(qpos) // 2])

        return ("".join(merged_seq), "".join(merged_qual))
    except Exception as e:
        print(f"Error merging reads: {e}")
        # 处理异常情况
        # 例如：返回None或记录错误信息
        return None


def construct_merged_header(barcode: str, count: int) -> str:
    """统一生成头信息前缀"""
    return f"@{barcode}_x{count}"


# =============================================================================
# 命令行接口
# =============================================================================
def main(
    input_r1: str,
    input_r2: str,
    primer_info: str,
    region_name: str,
    output_r1: str,
    output_r2: str,
    report_file: str,
    cutoff: int = 3,
):
    """
    主函数：处理命令行参数并调用处理函数
    Args:
    input_r1       : Read1输入文件路径
    input_r2       : Read2输入文件路径
    primer_info    : 引物信息文件路径
    region_name         : 目标区域名称
    output_r1      : Read1输出路径
    output_r2      : Read2输出路径
    report_file    : 合并报告文件路径（默认：merge_report.txt）
    cutoff         : 合并阈值（默认3，0表示仅去除barcode）
    """

    process_single_region(
        input_r1,
        input_r2,
        primer_info,
        region_name,
        output_r1,
        output_r2,
        report_file,
        cutoff,
    )


if __name__ == "__main__":
    fire.Fire(main)
