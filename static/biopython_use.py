#!/usr/bin/env python3
# coding=utf-8
# noinspection SpellCheckingInspection
"""For doing blasting.
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.NCBIWWW import qblast
from time import time, sleep, gmtime
"""
from typing import Generator, Union, Optional


def data_printer(blast_record: Generator, e_value_thresh: Union[float, int] = 0.04, line_len: int = 80) -> str:
    # noinspection SpellCheckingInspection
    """prints the results. Adapted from https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec125
    :param line_len: The max length of the sequence lines. (Min 4)
    :param blast_record: Generator = NCBIXML.parse(results)
    :param e_value_thresh: Treshold of e_value above which the results are not printed.
    :return: str - A textblock
    """
    assert e_value_thresh > 0
    assert line_len > 3
    alignments = ''
    line_len -= 3
    try:
        for record in blast_record:     # All 1 records
            for alignment in record.alignments:     # Each alignment
                for hsp in alignment.hsps:
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
    except TypeError or AttributeError:
        alignments = "Error in 'data_printer()': param 'blast_record' is malformed."
    return alignments


def bio_blaster(
        input_file: str, file_format: str, output_file: str, index: Optional[str] = None,
        program: str = 'blastn', database: str = 'nr', gi_format: bool = True, size: int = 10
) -> None:
    """Blasting sort of automated.
    :param input_file: The file path which file contains a/the sequence that is to be blasted.
    :param file_format: The format of the input_file.
    :param output_file: Where to store the results.
    :param index: If the file contains multiple seqs and you want to analise 1 place its ref code here.
    :param program: What blast program to use.
    :param database: Which database to use.
    :param gi_format: Whether to request the gi_format.
    :param size: The amount of results to request.
    :return: None.
    """
    # noinspection SpellCheckingInspection
    assert file_format in (  # https://biopython.org/wiki/SeqIO
        'abi', 'abi-trim', 'ace', 'cif-atom', 'cif-seqres', 'clustal', 'embl', 'fasta', 'fasta-2line', 'fastq-sanger',
        'fastq', 'fastq-solexa', 'fastq-illumina', 'gck', 'genbank', 'gb', 'ig', 'imgt', 'nexus', 'pdb-seqres',
        'pdb-atom', 'phd', 'phylip', 'pir', 'seqxml', 'sff', 'sff-trim', 'snapgene', 'stockholm', 'swiss', 'tab',
        'qual', 'uniprot-xml', 'xdna'
    )
    # noinspection SpellCheckingInspection
    assert program in (  # https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc92
        'blastn', 'blastp', 'blastx', 'tblast', 'tblastx'
    )
    if program == 'blastp' or program == 'blastx':  # Checking if the program supports the database.
        assert database in (
            'nr', 'refseq_protein', 'landmark', 'swissprot', 'pataa', 'pdb', 'env_nr', 'tsa_nr'
        )
    else:
        assert database in (
            'nr', 'nt', 'refseq_rna', 'refseq_representative_genomes', 'refseq_genomes', 'wgs', 'est', 'SRA', 'TSA',
            'HTGS', 'pat', 'pdb', 'RefSeq_Gene', 'gss', 'dbsts'
        )
    assert size > 0
    if output_file[-4:] == '.xml':
        output_file = output_file[:-4]
    # Bio = Biopython but Pycharm doesn't know.
    # noinspection PyPackageRequirements
    from Bio import SeqIO
    # noinspection PyPackageRequirements
    from Bio.Blast.NCBIWWW import qblast
    record_dict = SeqIO.index(input_file, format=file_format)
    record: str = ''
    if index is not None:   # Only one of the contained sequences.
        assert index in record_dict.keys()
        record = record_dict[index].format("fasta")
    else:   # All contained sequences
        for seq in record_dict:
            record += record_dict[seq].format("fasta") + '\n'
    del record_dict
    result_handle = qblast(     # The actual blasting see print(help(qblast))
        program=program, database=database, sequence=record, ncbi_gi=gi_format, hitlist_size=size, megablast=False
    )
    with open(output_file + '.xml', "w") as out_handle:  # Saving the results
        out_handle.write(result_handle.read())
    result_handle.close()   # The result handle is like an open file and must be closed.
    return


def biopython_use(
        result_loc: str, large_job: bool = True, input_file: Optional[str] = None, program: str = 'blastn',
        database: str = 'nr', file_format: Optional[str] = None, index: Optional[str] = None,
        print_results: bool = False, e_value_thresh: Union[float, int] = 0.04
) -> None:
    """Getting and/or printing results.
    :param program: What BLAST program to use.
    :param result_loc: Where to store/get the results.
    :param large_job: Check whether you are allowed to do large jobs.
    :param input_file: The file path which file contains a/the sequence that is to be blasted.
    :param database: The database to use.
    :param file_format: The format of the input_file.
    :param index: If the file contains multiple seqs and you want to analise 1 place its ref code here.
    :param print_results: Whether to print the results.
    :param e_value_thresh: if print_results is True; Threshold of e_value above which the results are not printed.
    :return: None.
    """
    from os.path import isdir
    assert e_value_thresh > 0 or not print_results
    assert isdir(result_loc)
    if input_file is not None:
        from time import time, sleep, gmtime
        if result_loc[-1:] == '/':
            result_loc = result_loc[:-1]
        result_loc += '/results5.xml'
        # https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
        sleep(0.4)  # post no more than three URL requests per second
        if large_job:  # Limit large jobs to either weekends or in between 9:00 PM to 5:00 AM (EST).
            job_time = gmtime(time())
            if not job_time[6] < 5 and not job_time[3] < 21:
                # noinspection SpellCheckingInspection
                bio_blaster(
                    input_file=input_file, file_format=file_format,
                    output_file=result_loc, index=index, program=program, database=database
                )
            del job_time
        else:
            bio_blaster(
                input_file=input_file, file_format=file_format, database=database,
                output_file=result_loc, index=index, program=program
            )
    if print_results:
        # Bio = Biopython but Pycharm doesn't know.
        # noinspection PyPackageRequirements
        from Bio.Blast import NCBIXML   # only importing if necicairy.
        with open(result_loc, 'r') as results:   # reopening the results for reading.
            for blast_records in NCBIXML.parse(results):
                print(data_printer(blast_record=blast_records, e_value_thresh=e_value_thresh))
    return
