#!/usr/bin/env python3
# coding=utf-8
"""
For doing blasting.
"""
# Bio = Biopython but Pycharm doesn't know.
# noinspection PyPackageRequirements
from Bio import SeqIO
# noinspection PyPackageRequirements
from Bio.Blast import NCBIXML
# noinspection PyPackageRequirements
from Bio.Blast.NCBIWWW import qblast
# print(help(qblast))


def data_printer(blast_record, e_value_thresh=0.04, line_len=80):
    # noinspection SpellCheckingInspection
    """
            :param line_len: The max length of the sequence lines. (Min 4)
            :param blast_record: NCBIXML.parse(result_handle).
            :param e_value_thresh: Treshold of e_value above which the results are not printed.
            :return: None.
            """
    if line_len > 3:
        line_len -= 3
    else:
        print("Error in 'data_printer()': param 'line_len' is too low. (min 4)")
        return
    for record in blast_record:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_value_thresh:
                    print("****Alignment****")
                    print("sequence:", alignment.title)
                    print("length:", alignment.length)
                    print("e value:", hsp.expect)
                    if len(hsp.query) > line_len:
                        print(hsp.query[0:line_len] + "...")
                        print(hsp.match[0:line_len] + "...")
                        print(hsp.sbjct[0:line_len] + "...\n")
                    else:
                        print(f'{hsp.query}\n{hsp.match}\n{hsp.sbjct}\n')
    return


def bio_blaster(
        input_file: str, file_format: str, output_file: str, index: str = None,
        program='blastn', database='nt', gi_format=True, size=10
):
    """Blasting sort of automated.
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


def biopython_use(
        result_loc: str, get_data=False, large_job=True, input_file: str = None, file_format: str = None,
        index: str = None, print_results=False, e_value_thresh=0.04
):
    """Getting and/or printing results.
    :param result_loc: Where to store/get the results.
    :param get_data: Whether to blast or not.
    :param large_job: Check whether you are allowed to do large jobs.
    :param input_file: The file path which file contains a/the sequence that is to be blasted.
    :param file_format: The format of the input_file.
    :param index: If the file contains multiple seqs and you want to analise 1 place its ref code here.
    :param print_results: Whether to print the results.
    :param e_value_thresh: if print_results is True; Threshold of e_value above which the results are not printed.
    :return: None.
    """
    if get_data:
        # https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
        from time import time, sleep, gmtime
        sleep(0.4)  # post no more than three URL requests per second
        if large_job:  # Limit large jobs to either weekends or in between 9:00 PM to 5:00 AM (EST).
            job_time = gmtime(time())
            if not job_time[6] < 5 and not job_time[3] < 21:
                # noinspection SpellCheckingInspection
                bio_blaster(
                    input_file=input_file, file_format=file_format,
                    output_file=result_loc,
                    index=index
                )
            del job_time
        else:
            bio_blaster(
                input_file=input_file, file_format=file_format,
                output_file=result_loc,
                index=index
            )
    results = open(result_loc)
    blast_records = NCBIXML.parse(results)
    if print_results:
        data_printer(blast_record=blast_records, e_value_thresh=e_value_thresh)
    return
