#!/usr/bin/env python3
# coding=utf-8
"""
Temp file not for usage.
"""
from static.biopython_use import biopython_use

biopython_use(
    result_loc='static', large_job=False, input_file='Course4_dataset_v04_mod.fastq', database='refseq_protein',
    file_format='fastq', index=None,
    print_results=False, program='blastx'
)
