from flask import Flask, render_template

app = Flask(__name__)


@app.route('/')
def website():
    return render_template("../blok4afvink2_2/templates/project.html")


if __name__ == '__main__':
    app.run()
