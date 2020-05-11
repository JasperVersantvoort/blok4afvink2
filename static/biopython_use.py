#!/usr/bin/env python3
# coding=utf-8
"""
TEMP file, for figuring shit out.
"""
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.NCBIWWW import qblast


# print(help(qblast))


def data_printer(blast_record, e_value_thresh=0.04):
    """
    :param blast_record: NCBIXML.parse(result_handle)
    :param e_value_thresh:
    :return: None
    """
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < e_value_thresh:
                print("****Alignment****")
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print("e value:", hsp.expect)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")
    return


if __name__ == '__main__':
    record_dict = SeqIO.index("/home/s_k_tiger/Downloads/Course4_dataset_v04_mod.fastq", format="fastq")
    record = record_dict['HWI-M02942:21:000000000-ACNW4:1:1101:19810:2642_1']
    result_handle = qblast(
        program='blastn', database='nt', sequence=record.format("fasta"), ncbi_gi=True, hitlist_size=10, megablast=False
    )
    with open("/home/s_k_tiger/Downloads/my_blast.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
        result_handle.close()
    result_handle = open("/home/s_k_tiger/Downloads/my_blast.xml")
    blast_records = NCBIXML.parse(result_handle)
    data_printer(blast_record=blast_records, e_value_thresh=0.04)
