from flask import Flask, render_template

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


@app.route('/database.html')
def site_database():
    return render_template("database.html")

if __name__ == '__main__':
    app.run()
