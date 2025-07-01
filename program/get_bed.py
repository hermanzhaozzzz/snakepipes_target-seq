# getbed.py

import fire


def generate_bed_files(
    primer_table: str, region_name: str, output_region_bed: str, output_focus_bed: str
):
    """
    生成区域和sgRNA的BED文件
    Args:
        primer_table: 引物信息表路径
        output_region_bed: 输出region bed文件路径
        output_focus_bed: 输出sgRNA bed文件路径
        region_name: 目标区域名称
    """
    region_name = str(region_name)
    with open(primer_table) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 13:
                continue
            if fields[3] != region_name:
                continue

            # 主区域BED
            chr_main = fields[0]
            start_main = str(int(fields[1]) - 1)
            end_main = fields[2]
            strand_main = fields[8]

            # sgRNA区域
            chr_sgrna = fields[10]
            start_sgrna = str(int(fields[11]) - 1)
            end_sgrna = fields[12]
            strand_sgrna = fields[8]

            # 写入主BED
            with open(output_region_bed, "w") as fout:
                fout.write(
                    f"{chr_main}\t{start_main}\t{end_main}\t{region_name}\t0\t{strand_main}\n"
                )

            # 写入sgRNA BED
            with open(output_focus_bed, "w") as fout:
                fout.write(
                    f"{chr_sgrna}\t{start_sgrna}\t{end_sgrna}\t{region_name}_sgRNA\t0\t{strand_sgrna}\n"
                )
            break


if __name__ == "__main__":
    fire.Fire(generate_bed_files)
