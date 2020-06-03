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


@app.route('/resultaten.html')
def resultaten():
    """
    Geeft een pagina weer met de belangrijkste resultaten
    :return: De html pagina met resultaten
    """
    return render_template("resultaten.html")


@app.route('/database.html', methods=["POST", "GET"])
def site_database():
    """
    Zorgt ervoor dat je kunt zoeken in de database met de POST methode.
    :return: Geeft de database html pagina weer
    """
    if request.method == "POST":
        zoeken = request.form.get("zoek", "None")
        rows = database(zoeken)
        return render_template("database.html", database=rows,
                               zoek=zoeken)
    else:
        rows = database("None")
        return render_template("database.html", database=rows,
                               zoek="None")


def database(zoek):
    """
    Geeft de description van de database en kan deze filteren op een
    zoekwoord
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
        cursor.execute("select acc_code, name, "
                       "description, e_value from results "
                       "join organism on "
                       "organism_id = organism.id order by name")
        rows = cursor.fetchall()
        cursor.close()
        conn.close()
        return rows
    else:
        query = "select acc_code, name, description,e_value from " \
                "results join organism on organism_id = organism.id  " \
                "where name like '%" + zoek + "%' or acc_code like " \
                                              "'%" + zoek + "%' or description like '%" + zoek + "%'" \
                                                                                                 "order by name"
        cursor.execute(query)
        rows = cursor.fetchall()
        cursor.close()
        conn.close()
        return rows


def evalue_nul():
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
    cursor.execute("select acc_code, name, description,e_value from "
                   "results join organism on organism_id =  "
                   "organism.id where e_value = 0")
    rij_1 = cursor.fetchall()
    cursor.close()
    conn.close()
    return rij_1


if __name__ == '__main__':
    app.run()
