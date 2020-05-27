# Bin1d groep 3
# XML data in database zetten
import mysql.connector
from xml.etree import ElementTree


def connector(blast_version, header, hit_id, acc, perc, tscore,
                          evalue, defen, qseq ):

    conn = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        db='mlfrg',
        user='mlfrg@hannl-hlo-bioinformatica-mysqlsrv',
        password='chocolade45')
    cursor = conn.cursor()
    cursor.execute(
        "insert into research_sequence (sequence,score, header, `read`, plusorminux) values ('"qseq"', 0, "header", 1, 'plus');")
    rows = cursor.fetchall()
    cursor.close()
    conn.close()


def open_xml(xml_file):
    # Parse door xml file
    par = ElementTree.parse(xml_file)

    #queries = par.findall("./BlastOutput_iterations/Iteration/Iteration_hits")
    count = 1
    iteration = par.findall(path="./BlastOutput_iterations/Iteration")
    blast_version = par.find("./BlastOutput_program").text
    for queries in iteration:
        header = queries.find('Iteration_query-def').text
        for hits in queries:

            for hit in hits:
                hit_id = hit.find("Hit_id")
                defen = hit.find("Hit_def")
                acc = hit.find("Hit_accession")
                perc = hit.find("Hit_hsps/Hsp/Hsp_identity")
                tscore = hit.find("Hit_hsps/Hsp/Hsp_score")
                evalue = hit.find("Hit_hsps/Hsp/Hsp_evalue")
                qseq = hit.find("Hit_hsps/Hsp/Hsp_qseq")
                if hit_id is not None:
                    hit_id = hit_id.text
                if acc is not None:
                    acc = acc.text
                if perc is not None:
                    perc = perc.text
                if tscore is not None:
                    tscore = tscore.text
                if evalue is not None:
                    evalue = evalue.text
                if defen is not None:
                    defen = defen.text
                if qseq is not None:
                    qseq = qseq.text

                connector(blast_version, header, hit_id, acc, perc, tscore, evalue, defen, qseq)







                # evalue = hit.find("Hit_hsps/Hsp/Hsp_evalue").text
                # print("E value = ", evalue)


def main():
    xml_file = "Tester.xml"
    open_xml(xml_file)


main()
