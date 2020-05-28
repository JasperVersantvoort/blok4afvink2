# Bin1d groep 3
# XML data in database zetten
import mysql.connector
from xml.etree import ElementTree


def connector(blast_version, header, hit_id, acc, perc, tscore,
              evalue, defen, qseq):
    conn = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        db='mlfrg',
        user='mlfrg@hannl-hlo-bioinformatica-mysqlsrv',
        password='chocolade45')
    cursor = conn.cursor()

    check_header_execute = "select header from research_sequence " \
                           "where header = '" + header + "';"
    cursor.execute(check_header_execute)
    check = cursor.fetchall()

    if len(check) == 0:
        print("Nieuwe sequentie header: ", header)
        if qseq is not None and header is not None:
            execute_res_seq = "insert into research_sequence " \
                              "(sequence,score, header, `read`, plusorminux) " \
                              "values ('" + qseq + "',0, '" + header + "'," + \
                              header[-1] + ", 'null')"

            cursor.execute(execute_res_seq)
            conn.commit()
            cursor.close()
            cursor = conn.cursor()

    if perc is not None and acc is not None and evalue is not None and blast_version is not None and defen is not None:
        print(header, "hit", hit_id)
        res_seq_id_execute = "select id from research_sequence " \
                             "where header = '" + header + "'"
        cursor.execute(res_seq_id_execute)
        res_seq_id = cursor.fetchall()
        execute_result = "insert into results (percent_identity, acc_code, e_value,total, max, research_sequence_id, blast_version, description)" \
                         " values (" + perc + ",'" + acc + "'," + evalue + ",0,0," + str(
            res_seq_id[0][
                0]) + ",'" + blast_version + "','" + defen + "');"
        cursor.execute(execute_result)
        conn.commit()

    cursor.close()
    conn.close()


def open_xml(xml_file):
    # Parse door xml file
    par = ElementTree.parse(xml_file)

    # queries = par.findall("./BlastOutput_iterations/Iteration/Iteration_hits")
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

                connector(blast_version, header, hit_id, acc, perc, tscore,
                          evalue, defen, qseq)

                # evalue = hit.find("Hit_hsps/Hsp/Hsp_evalue").text
                # print("E value = ", evalue)


def main():
    xml_file = "Tester.xml"
    open_xml(xml_file)


main()
