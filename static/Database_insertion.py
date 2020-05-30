# Bin1d groep 3
# XML data in database zetten
import mysql.connector
from xml.etree import ElementTree


def connector(blast_version, header, hit_id, acc, ident,
              evalue, defen, qseq, hit_num):
    """
    Connect met de mlfrg, en zet de resultaten in de database

    blast_version = blast die gebruikt
    header = header van sequentie
    hit_id = het id van die hit
    acc = de accessie code van de hit
    ident = de identitie van de hit
    evalue = de e value van de hit
    defen = de defenitie van de hit
    qseq = het gebruikte stuk sequentie voor de hit
    hit_num = Het nummer van de hit

    return: De restultaten in de mlfrg database
    """

    conn = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        db='mlfrg',
        user='mlfrg@hannl-hlo-bioinformatica-mysqlsrv',
        password='chocolade45')
    cursor = conn.cursor()

    # Checkt welke blast is uitgevoerd,
    # als hij deze niet kent voegt hij die toe aan de database
    check_blast_execute = \
        "select blast_version from blast_version " \
        "where blast_version = '" + blast_version + "';"
    cursor.execute(check_blast_execute)
    blast_check = cursor.fetchall()

    if len(blast_check) == 0 and blast_version is not None:
        execute_blast = "insert into blast_version (blast_version) " \
                        "values ('" + blast_version + "')"
        cursor.execute(execute_blast)
        conn.commit()
        cursor.close()
        cursor = conn.cursor()

    # Kijkt of er een nieuwe header is en voegt deze toe
    check_header_execute = "select header from research_sequence " \
                           "where header = '" + header + "';"
    cursor.execute(check_header_execute)
    check = cursor.fetchall()

    if len(check) == 0:
        print("Nieuwe sequentie header: ", header)
        if qseq is not None and header is not None:
            execute_res_seq = "insert into research_sequence " \
                              "(sequence,score, header, `read`," \
                              " plusorminux) " \
                              "values ('" + qseq + "',null, '" + \
                              header + "'," + \
                              header[-1] + ", null)"

            cursor.execute(execute_res_seq)
            conn.commit()
            cursor.close()
            cursor = conn.cursor()

    # kijkt of er een nieuw organism is.
    # Wanneer deze nieuw is voegt hij deze toe aan de database
    if defen is not None:
        check_organism_execute = "select id from organism " \
                                 "where name = '" + defen.split(',')[0] + "';"
        cursor.execute(check_organism_execute)
        check_organism = cursor.fetchall()

        if len(check_organism) == 0:
            org_id = defen.split(',')[0]
            execute_organism = "insert into organism(name, results_id) " \
                               "values ('" + defen.split(',')[0] + "',0)"
            cursor.execute(execute_organism)
            conn.commit()
            cursor.close()
            cursor = conn.cursor()
        else:
            org_id = check_organism[0][0]

    if ident is not None and acc is not None and evalue is not None \
            and blast_version is not None and defen is not None:
        # kijkt of de data aan de eisen voldoet om in de database te mogen
        if int(hit_num) < 10 and float(evalue) < 0.00005 and int(ident) > 40:
            print(header, "hit", hit_id)
            res_seq_id_execute = "select id from research_sequence " \
                                 "where header = '" + header + "'"
            cursor.execute(res_seq_id_execute)
            res_seq_id = cursor.fetchall()
            execute_result = "insert into results (percent_identity," \
                             " acc_code, e_value,total, max," \
                             " research_sequence_id, description," \
                             " organism_id,blast_version_blast_version)" \
                             " values (" + ident + ",'" + acc + "',"\
                             + evalue + ",null,null," + str(
                res_seq_id[0][0]) + ",'" + defen + "'," + str(
                org_id) + ",'" + blast_version + "');"
            cursor.execute(execute_result)
            conn.commit()

    cursor.close()
    conn.close()


def open_xml(xml_file):
    # Parse door xml file
    par = ElementTree.parse(xml_file)

    # queries = par.findall("./BlastOutput_iterations/Iteration/Iteration_hits")
    iteration = par.findall(path="./BlastOutput_iterations/Iteration")
    blast_version = par.find("./BlastOutput_program").text
    for queries in iteration:
        header = queries.find('Iteration_query-def').text
        for hits in queries:

            for hit in hits:
                hit_id = hit.find("Hit_id")
                defen = hit.find("Hit_def")
                hit_num = hit.find("Hit_num")
                acc = hit.find("Hit_accession")
                ident = hit.find("Hit_hsps/Hsp/Hsp_identity")
                tscore = hit.find("Hit_hsps/Hsp/Hsp_score")
                evalue = hit.find("Hit_hsps/Hsp/Hsp_evalue")
                qseq = hit.find("Hit_hsps/Hsp/Hsp_qseq")

                if hit_id is not None:
                    hit_id = hit_id.text
                if acc is not None:
                    acc = acc.text
                if ident is not None:
                    ident = ident.text
                if tscore is not None:
                    tscore = tscore.text
                if evalue is not None:
                    evalue = evalue.text
                if defen is not None:
                    defen = defen.text
                if qseq is not None:
                    qseq = qseq.text
                if hit_num is not None:
                    hit_num = hit_num.text

                connector(blast_version, header, hit_id, acc, ident,
                          evalue, defen, qseq, hit_num)


def main():
    xml_file = "Tester.xml"
    open_xml(xml_file)


main()
