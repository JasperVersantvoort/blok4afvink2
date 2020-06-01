import mysql.connector
from flask import Flask, render_template, request

app = Flask(__name__)


@app.route('/')
def website():
    return render_template("project.html")


@app.route('/project.html')
def website_home():
    return render_template("project.html")


@app.route('/over ons.html')
def site_over_ons():
    return render_template("over ons.html")


@app.route('/database.html', methods=["POST", "GET"])
def site_database():
    if request.method == "POST":
        zoeken = request.form.get("zoek", "")
        rows = database(zoeken)

        return render_template("database.html", database=rows,
                               zoek=zoeken)

    else:
        rows = database("None")
        return render_template("database.html", database=rows,
                               zoek="None")


def database(zoek):
    """ haalt de description uit de ensembldb database
     en filtert deze op het zoekwoord

    :param zoek: Het ingegeven zoekwoord
    :return: Een lijst met de juiste discriptions
    """
    print("zoek woord is:", zoek)
    conn = mysql.connector.connect(host='hannl-hlo-bioinformatica'
                                        '-mysqlsrv.mysql.database.azure.com',
                                   user='mlfrg@hannl-hlo-bioinformatica-mysqlsrv',
                                   password = 'chocolade45', db='mlfrg')
    cursor = conn.cursor()
    cursor.execute("select description from results")
    rows = cursor.fetchall()
    des = []
    for row in rows:
        if str(row) != "(None,)":
            if zoek.upper() in str(row).upper() or zoek is 'None':
                des.append(row)

    cursor.close()
    conn.close()
    return des

if __name__ == '__main__':
    app.run()
