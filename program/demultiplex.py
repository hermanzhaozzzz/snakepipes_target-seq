#!/usr/bin/env python3
import gzip
import os
import time
from collections import defaultdict
import fire

class Primer:
    def __init__(self, line):
        fields = line.strip().split('\t')
        self.region_index = fields[3]
        self.left_barcode_length = int(fields[4])
        self.right_barcode_length = int(fields[5])
        self.left_primer = fields[6].upper()
        self.right_primer = fields[7].upper()
        
        # Precompute primer properties
        self.left_len = len(self.left_primer)
        self.right_len = len(self.right_primer)
        self.left_prefix = self.left_primer[:4]  # First 4 bases for quick check
        self.right_prefix = self.right_primer[:4]  # First 4 bases for quick check

class FastqWriter:
    def __init__(self, output_dir, buffer_size=100000, max_reads=5000000):
        self.max_reads = max_reads
        self.output_dir = output_dir
        self.buffer_size = buffer_size
        self.buffers = defaultdict(lambda: {'R1': [], 'R2': []})
        self.handles = {}
        self.counts = defaultdict(int)  # Track read counts per primer
        
        os.makedirs(output_dir, exist_ok=True)

    def add_read(self, primer_id, r1, r2):
        # Skip if already reached limit
        if self.counts[primer_id] >= self.max_reads:
            return
        
        # Check if adding this read would exceed limit
        if self.counts[primer_id] + 1 > self.max_reads:
            return
        
        self.buffers[primer_id]['R1'].append(r1)
        self.buffers[primer_id]['R2'].append(r2)
        self.counts[primer_id] += 1  # Increment immediately to ensure atomicity
        
        if len(self.buffers[primer_id]['R1']) >= self.buffer_size:
            self._flush_buffer(primer_id)

    def _flush_buffer(self, primer_id):
        if primer_id not in self.handles:
            r1_path = os.path.join(self.output_dir, f"{primer_id}_R1.fastq.gz")
            r2_path = os.path.join(self.output_dir, f"{primer_id}_R2.fastq.gz")
            self.handles[primer_id] = {
                'R1': gzip.open(r1_path, 'at'),
                'R2': gzip.open(r2_path, 'at')
            }
            
        h = self.handles[primer_id]
        for read in self.buffers[primer_id]['R1']:
            h['R1'].write('\n'.join(read) + '\n')
        for read in self.buffers[primer_id]['R2']:
            h['R2'].write('\n'.join(read) + '\n')
        
        self.buffers[primer_id]['R1'].clear()
        self.buffers[primer_id]['R2'].clear()

    def close_all(self):
        for pid in list(self.buffers.keys()):
            self._flush_buffer(pid)
        for h in self.handles.values():
            h['R1'].close()
            h['R2'].close()

    def is_full(self, primer_id):
        return self.counts.get(primer_id, 0) >= self.max_reads

def fastq_reader(file_path):
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header: break
            seq = f.readline().strip()
            _ = f.readline()  # +
            qual = f.readline().strip()
            yield (header, seq, qual)

def match_primer(seq, primer_seq, primer_prefix, max_mismatch):
    # Quick prefix check
    prefix = seq[:4]
    mismatches = sum(c1 != c2 for c1, c2 in zip(prefix, primer_prefix))
    if mismatches > max_mismatch:
        return False
    
    # Full sequence check
    mismatches = 0
    for c1, c2 in zip(seq, primer_seq):
        if c1 != c2:
            mismatches += 1
            if mismatches > max_mismatch:
                return False
    return True

def process_pairs(primers, r1_path, r2_path, output_dir, max_mismatch=1, buffer_size=100000, max_reads=5000000):
    writer = FastqWriter(output_dir, buffer_size=buffer_size, max_reads=max_reads)
    r1_gen = fastq_reader(r1_path)
    r2_gen = fastq_reader(r2_path)
    
    full_primers = set()
    
    for (h1, s1, q1), (h2, s2, q2) in zip(r1_gen, r2_gen):
        # Early exit if all primers are full
        if len(full_primers) == len(primers):
            break
            
        s1, s2 = s1.upper(), s2.upper()
        found = False
        
        for primer in primers:
            if primer.region_index in full_primers:
                continue
                
            # Extract primer regions
            r1_start = primer.left_barcode_length
            r1_end = r1_start + primer.left_len
            if r1_end > len(s1): continue
            
            r2_start = primer.right_barcode_length
            r2_end = r2_start + primer.right_len
            if r2_end > len(s2): continue
            
            r1_sub = s1[r1_start:r1_end]
            r2_sub = s2[r2_start:r2_end]
            
            # Check match
            if (match_primer(r1_sub, primer.left_primer, primer.left_prefix, max_mismatch) and
                match_primer(r2_sub, primer.right_primer, primer.right_prefix, max_mismatch)):
                writer.add_read(primer.region_index, 
                              [h1, s1, '+', q1],
                              [h2, s2, '+', q2])
                
                # Check if this primer just became full
                if writer.is_full(primer.region_index):
                    full_primers.add(primer.region_index)
                
                found = True
                break  # First match policy
        
        if not found:
            # Optional: Handle unmatched reads
            pass
    
    writer.close_all()

def main(primer_table: str, r1: str, r2: str, output_dir: str, mismatch=1, buffer_size=100000, max_reads=5000000):
    """
    Demultiplex paired-end fastq files based on primer sequences.

    Args:
        primer_table (str): Path to the primer table file.
        r1 (str): Path to the forward read fastq file.
        r2 (str): Path to the reverse read fastq file.
        output_dir (str): Directory to save demultiplexed fastq files.
        mismatch (int): Maximum number of mismatches allowed for primer matching.
        buffer_size (int): Number of reads to buffer before writing to disk.
        max_reads (int): Maximum number of reads per primer region.
    """
    start_time = time.time()
    
    # Load primers
    with open(primer_table) as f:
        headers = f.readline()  # Skip header
        primers = [Primer(line) for line in f]
    
    print(f"Loaded {len(primers)} primers")
    process_pairs(primers, r1, r2, output_dir, mismatch, buffer_size=buffer_size, max_reads=max_reads)
    print(f"Completed in {time.time()-start_time:.2f} seconds")

if __name__ == '__main__':
    fire.Fire(main)
