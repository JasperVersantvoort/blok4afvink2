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


def data_printer(blast_record, e_value_thresh=0.04, line_len=80) -> str:
    # noinspection SpellCheckingInspection
    """prints the results.
    :param line_len: The max length of the sequence lines. (Min 4)
    :param blast_record: NCBIXML.parse(result_handle).
    :param e_value_thresh: Treshold of e_value above which the results are not printed.
    :return: str - A textblock
    """
    alignments = ''
    if line_len > 3:
        line_len -= 3
    else:
        alignments = "Error in 'data_printer()': param 'line_len' is too low. (min 4)"
        return alignments
    for record in blast_record:     # All 1 records
        for alignment in record.alignments:     # Each alignment
            for hsp in alignment.hsps:  # hsp = ???
                if hsp.expect < e_value_thresh:
                    alignments += '****Alignment****\n'
                    alignments += "sequence: " + alignment.title + '\n'
                    alignments += "length: " + alignment.length + '\n'
                    alignments += "e value: " + hsp.expect + '\n'
                    if len(hsp.query) > line_len:
                        alignments += hsp.query[0:line_len] + "...\n"
                        alignments += hsp.match[0:line_len] + "...\n"
                        alignments += hsp.sbjct[0:line_len] + "...\n\n"
                    else:
                        alignments += f'{hsp.query}\n{hsp.match}\n{hsp.sbjct}\n\n'
    return alignments


def bio_blaster(
        input_file: str, file_format: str, output_file: str, index: str = None,
        program='blastn', database='nt', gi_format=True, size=10
) -> str:
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
    from time import sleep
    record_dict = SeqIO.index(input_file, format=file_format)
    if index is not None:
        record = [record_dict[index].format("fasta")]
    else:
        record = []
        for seq in record_dict:
            record.append(record_dict[seq].format("fasta"))
    del record_dict
    result_handle_store = ''
    for req in record:
        if index is None:
            sleep(0.4)
        if program == 'blastn':
            result_handle = qblast(     # blastn has the option 'megablast' which the others do not.
                program='blastn', database=database, sequence=req,
                ncbi_gi=gi_format, hitlist_size=size, megablast=False
            )
        else:
            result_handle = qblast(     # The actual blasting
                program=program, database=database, sequence=req,
                ncbi_gi=gi_format, hitlist_size=size
            )
        result_handle_store += result_handle.read()
        result_handle.close()   # The result handle is like an open file and must be closed.
    with open(output_file + '.xml', "w") as out_handle:  # Saving the results
        out_handle.write(result_handle_store)
    return result_handle_store


def biopython_use(
        result_loc: str, get_data=False, large_job=True, input_file: str = None, file_format: str = None,
        index: str = None, print_results=False, e_value_thresh=0.04
) -> None:
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
    if print_results:
        results = open(result_loc)  # reopening the results for reading.
        blast_records = NCBIXML.parse(results)
        print(data_printer(blast_record=blast_records, e_value_thresh=e_value_thresh))
    return
