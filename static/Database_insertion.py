#!/usr/bin/env python3
# coding=utf-8
"""
Putting xml files into the database.
"""


def connector(
        host: str, db: str, user: str, password: str, command: str, fetch: bool = False, commit: bool = False
) -> str:
    """For the database operations
    :param host: The host of the database.
    :param db: The database to connect to.
    :param user: The user to log in as.
    :param password: Your password.
    :param command: What command to execute.
    :param fetch: Whether to collect the output.
    :param commit: Whether to commit any made changes.
    :return: If fetch is True the collected output else ''.
    """
    import mysql.connector
    try:
        conn = mysql.connector.connect(host=host, user=user, passwd=password, db=db)
        del password
        cursor = conn.cursor()
        cursor.execute(command)
        if commit:
            conn.commit()
        if fetch:
            fetch = cursor.fetchall()
        else:
            fetch = ''
    except mysql.connector.Error as error:
        print(error)
        return ''
    cursor.close()
    conn.close()
    return fetch


def xml_to_database(xml_file: str) -> None:
    """TODO add database credentials and check contents.
    :param xml_file: The file whose data is to be exported.
    :return: None
    """
    # Bio = Biopython but Pycharm doesn't know.
    # noinspection PyPackageRequirements
    from Bio.Blast import NCBIXML
    # noinspection SqlNoDataSourceInspection
    codes = connector(
        host='', db='', user='', password='',
        command='select Hit_id from data', fetch=True
    )
    with open(xml_file, 'r') as xml:
        blast_record = NCBIXML.parse(xml)
        try:
            for record in blast_record:
                for alignment in record.alignments:     # Each alignment
                    if alignment.Hit_id not in codes:
                        for hsp in alignment.hsps:
                            # noinspection SqlNoDataSourceInspection
                            connector(
                                host='', db='', user='', password='',
                                command='INSERT INTO data '
                                        '(Hit_id, title, length, "e-value", query, match, subject)'
                                        ' VALUES '
                                        f'({alignment.Hit_id}, {alignment.title}, {alignment.length}, {hsp.expect}, '
                                        f'{hsp.query}, {hsp.match}, {hsp.sbjct})', commit=True
                            )
        except TypeError or AttributeError:
            print("Did you use correct credentials?")
    return
