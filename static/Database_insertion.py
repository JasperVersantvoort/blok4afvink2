# Bin1d groep 3
# XML data in database zetten
import mysql.connector
from xml.etree import ElementTree


def connector():
    conn = None
    # try:
    conn = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        db='mlfrg',
        user='mlfrg@hannl-hlo-bioinformatica-mysqlsrv',
        password='chocolade45')
    cursor = conn.cursor()
    cursor.execute("select * from results")
    rows = cursor.fetchall()
    print(rows)
    cursor.close()
    conn.close()


def open_xml(xml_file):


    # Parse door xml file
    par = ElementTree.parse(xml_file)
    # hits = par.findall("/BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    hits = par.findall("./BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    print(hits)

    # hits = find_hits()
    # for regel in open_xml:
    #     if regel.startswith('<!DOCTYPE'):


def main():
    xml_file = "results.xml"
    test = connector()
    open_xml(xml_file)


main()
