#!/usr/bin/env python3
# coding=utf-8
"""
TEMP file, for figuring shit out.
"""
import time  # https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen

# Bio = Biopython but Pycharm doesn't know.
# noinspection PyPackageRequirements
from Bio import SeqIO
# noinspection PyPackageRequirements
from Bio.Blast import NCBIXML
# noinspection PyPackageRequirements
from Bio.Blast.NCBIWWW import qblast


# print(help(qblast))


def data_printer(blast_record, e_value_thresh=0.04):
    # noinspection SpellCheckingInspection
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


def bio_blaster(
        input_file: str, file_format: str, output_file: str, index: str = None,
        program='blastn', database='nt', gi_format=True, size=10
):
    """Blasting sort of automated. TODO: testing
    :param input_file: The file path which file contains a/the sequence that is to be blasted.
    :param file_format: The format of the input_file.
    :param output_file: Where to store the results.
    :param index: If the file contains multiple seqs and you want to analise 1 place its ref code here.
    :param program: What blast program to use.
    :param database: Which database to use.
    :param gi_format: Whether to request the gi_format.
    :param size: The amount of results to request.
    :return: 'result_handle_bak', not sure if this data can be used.
    """
    record_dict = SeqIO.index(input_file, format=file_format)
    if index is not None:
        record = record_dict[index]
    else:
        record = record_dict
    del record_dict
    if program == 'blastn':
        result_handle = qblast(
            program='blastn', database=database, sequence=record.format("fasta"),
            ncbi_gi=gi_format, hitlist_size=size, megablast=False
        )
    else:
        result_handle = qblast(
            program=program, database=database, sequence=record.format("fasta"),
            ncbi_gi=gi_format, hitlist_size=size
        )
    result_handle_bak = result_handle
    with open(output_file, "w") as out_handle:
        out_handle.write(result_handle.read())
        result_handle.close()
    return result_handle_bak


large_job = False
if __name__ == '__main__':
    time.sleep(0.4)  # post no more than three URL requests per second
    if large_job:  # Limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays.
        job_time = time.gmtime(time.time())
        if not job_time[6] < 5 and not job_time[3] < 21:
            # noinspection SpellCheckingInspection
            data_bak = bio_blaster(
                input_file="/home/s_k_tiger/Downloads/Course4_dataset_v04_mod.fastq", file_format="fastq",
                output_file="/home/s_k_tiger/Downloads/my_blast.xml",
                index='HWI-M02942:21:000000000-ACNW4:1:1101:19810:2642_1'
            )
        del job_time
    results = open("/home/s_k_tiger/Downloads/my_blast.xml")
    blast_records = NCBIXML.parse(results)
    data_printer(blast_record=blast_records, e_value_thresh=0.04)
