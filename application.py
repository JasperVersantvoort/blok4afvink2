# Course 4, groep d3 in Azure
# Jasper Versantvoort, Yuri Wit, Lars Weijenborg en Christel van Haren

import mysql.connector
from flask import Flask, render_template, request

app = Flask(__name__)


@app.route('/')
def website():
    """
    Zorgt voor de homepagina op de website met informatie over het
    project.
    :return: De basis html pagina
    """
    return render_template("project.html")


@app.route('/project.html')
def website_home():
    """
    Zorgt ervoor dat je ook weer terug kunt komen op de homepagina.
    :return: de basis html pagina
    """
    return render_template("project.html")


@app.route('/resultaten.html', methods=["POST", "GET"])
def resultaten():
    """

    :return:
    """
    conn = mysql.connector.connect(host='hannl-hlo-bioinformatica'
                                        '-mysqlsrv.mysql.database'
                                        '.azure.com',
                                   user='mlfrg@hannl-hlo'
                                        '-bioinformatica-mysqlsrv',
                                   password='chocolade45', db='mlfrg')
    cursor = conn.cursor()

    aantal = "select o.name, count(organism_id) as Occurence from  " \
            "results join organism o on results.organism_id = o.id " \
            "group by organism_id order by Occurence DESC limit 10;"
    print(aantal)
    cursor.execute(aantal)
    rows = cursor.fetchall()
    cursor.close()
    conn.close()
    return rows


@app.route('/blast.html', methods=["POST", "GET"])
def blast():
    """
    Geeft een pagina weer waar je blast kunt toevoegen aan de database
    :return: De html pagina met blast toevoegen aan de database
    """
    if request.method == "GET":
        blast_version = request.args.get("blast", "None")
        header = request.args.get("header", "None")
        hit_id = request.args.get("hit_id", "None")
        acc = request.args.get("acc", "None")
        organisme = request.args.get("organisme", "None")
        ident = request.args.get("ident", "None")
        # kijkt of de ident aan de eisen voeldoet en of het een getal is
        try:
            ident = int(ident)
            if ident < 40:
                return render_template("blast.html",
                                       reactie="Te lage identity")
            else:
                ident = str(ident)
        except ValueError:
            return render_template("blast.html",
                                   reactie="Geen geldige identity")
        # kijkt of de evalue aan de eisen voeldoet en of het een getal is
        evalue = request.args.get("evalue", "None")
        try:
            evalue = float(evalue)
            if evalue >= 0.00001:
                return render_template("blast.html",
                                       reactie=
                                       "Te hoge e-value")
            else:
                evalue = str(evalue)
        except ValueError:
            return render_template("blast.html",
                                   reactie=
                                   "Geen geldige e-value")

        defen = request.args.get("defenitie", "None")
        qseq = request.args.get("qseq", "None")
        if blast_version == "None" or header == "None" or hit_id == "None" \
                or acc == "None" or organisme == "None" or organisme == "None" \
                or ident == "None" or evalue == "None" or defen == "None" \
                or qseq == "None":
            # print("if dus minimaal een is None")
            print(blast_version, header, hit_id)
            return render_template("blast.html",
                                   reactie=
                                   "Het invoeren is mislukt")
        else:
            print("Alles heeft een waarde")
            connector(blast_version, header, hit_id, acc, organisme, ident,
                      evalue, defen, qseq)
            return render_template("blast.html",
                                   reactie=
                                   "De Blast resultaten zijn toegevoegd!")


@app.route('/database.html', methods=["POST", "GET"])
def site_database():
    """
    Zorgt ervoor dat je kunt zoeken in de database em dat je op de
    e-value kunt filteren met de POST methode.
    :return: Geeft de database html pagina weer
    """
    if request.method == "POST":
        zoeken = request.form.get("zoek", "None")
        print("Het zoekwoord is", zoeken)
        if zoeken == "":
            zoeken = "None"
        e_value = str(request.form.get("evalue", '1'))
        rows = database(zoeken, e_value)
        return render_template("database.html", database=rows,
                               zoek=zoeken, evalue=e_value)
    else:
        e_value = str(request.form.get("evalue", '1'))
        rows = database("None", e_value)
        return render_template("database.html", database=rows,
                               zoek="None", evalue=e_value)


def database(zoek, e_value):
    """
    Geeft de description van de database en kan deze filteren op een
    zoekwoord en filteren op de e-value.
    :param e_value: filtert op evalue
    :param zoek: Het ingegeven zoekwoord
    :return: Een lijst met de juiste discriptions
    """

    conn = mysql.connector.connect(host='hannl-hlo-bioinformatica'
                                        '-mysqlsrv.mysql.database'
                                        '.azure.com',
                                   user='mlfrg@hannl-hlo'
                                        '-bioinformatica-mysqlsrv',
                                   password='chocolade45', db='mlfrg')
    cursor = conn.cursor()

    if zoek == "None":
        print("Het zoekwoord is:", zoek)
        query = "select acc_code, name, " \
                "description, e_value from results " \
                "join organism on " \
                "organism_id = organism.id " \
                "where e_value <= " + e_value + " order " \
                                                "by name "
        # print (query)
        cursor.execute(query)
        rows = cursor.fetchall()
        cursor.close()
        conn.close()
        return rows
    else:
        query = "select acc_code, name, description,e_value from " \
                "results join organism on organism_id = organism.id  " \
                "where (name like '%" + zoek + "%' or acc_code like  " \
                                               "'%" + zoek + "%' or description like '%" + zoek + \
                "%')  and e_value <= " + e_value + " order by name "
        # print(query)
        cursor.execute(query)
        rows = cursor.fetchall()
        cursor.close()
        conn.close()
        return rows


def connector(blast_version, header, hit_id, acc, organisme, ident,
              evalue, defen, qseq):
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

    return: De restultaten in de mlfrg database
    """

    conn = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        db='mlfrg',
        user='mlfrg@hannl-hlo-bioinformatica-mysqlsrv',
        password='chocolade45')
    cursor = conn.cursor()
    try:
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
            print(execute_blast)
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
            # Als er in organism een ' staat dan wordt die weg gehaald
            # Anders geeft dit fouten in de query
            if "'" in organisme:
                organism = organisme.replace("'", "")
            check_organism_execute = "select id from organism " \
                                     "where name = '" + organisme + "';"
            print("Organisme: ", check_organism_execute)
            cursor.execute(check_organism_execute)
            check_organism = cursor.fetchall()

            if len(check_organism) == 0:
                execute_organism = "insert into organism(name, results_id) " \
                                   "values ('" + organisme + "',0)"
                print(execute_organism)
                cursor.execute(execute_organism)
                conn.commit()
                cursor.close()
                cursor = conn.cursor()
                check_organism_execute = "select id from organism " \
                                         "where name = '" + organisme + "';"
                cursor.execute(check_organism_execute)
                check_organism = cursor.fetchall()
                org_id = check_organism[0][0]
            else:
                org_id = check_organism[0][0]

        # Als er in defen een ' staat dan wordt die weg gehaald
        # Anders geeft dit fouten in de query
        if "'" in defen:
            defen = str(defen).replace("'", "")

        res_seq_id_execute = "select id from research_sequence " \
                             "where header = '" + header + "'"
        cursor.execute(res_seq_id_execute)
        res_seq_id = cursor.fetchall()
        execute_result = "insert into results (percent_identity," \
                         " acc_code, e_value,total, max," \
                         " research_sequence_id, description," \
                         " organism_id,blast_version_blast_version)" \
                         " values (" + ident + ",'" + acc + "'," \
                         + evalue + ",null,null," + str(
            res_seq_id[0][0]) + ",'" + defen + "'," + str(
            org_id) + ",'" + blast_version + "');"
        print("Results: ", execute_result)
        cursor.execute(execute_result)
        conn.commit()
        cursor.close()
        conn.close()
    except IndexError:
        print("Index Error")
        pass


if __name__ == '__main__':
    app.run()
