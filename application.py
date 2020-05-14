#!/usr/bin/env python3
# coding=utf-8
"""TODO: fil this in
"""
from flask import Flask, render_template, request

app = Flask(__name__)


@app.route('/')
def dna_naar_eiwit():
    """

    :return: Webpagina van dnanaareiwit.html weergeven
    """
    seq = request.args.get("seq", '')
    # seq.lower, zodat de sequentie zowel groot als klein ingevoerd kan worden
    eiwit = eiwitje(seq.lower())

    return render_template("dnanaareiwit.html", seq=eiwit)


def eiwitje(seq: str) -> str:
    """

    :param seq: De ingevoerde string van de webpagina
    :return: De eiwit sequentie, of een tekst die aangeeft
    wat fout is met de sequentie
    """
    code = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
            'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
            'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
            'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
            'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
            'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
            'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
            'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
            'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
            'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
            'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
            'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
            'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
            'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
            'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
            'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
            }
    eiwit = ''
    try:
        if len(seq) == 0:
            eiwit = "Geef goede sequentie"

        elif len(seq)/3 != float(len(seq)/3):
            eiwit = "Lengte sequentie is niet deelbaar door 3"
        else:
            dna = True
            for s in seq.lower():
                if s not in ["a", "t", "g", "c"]:
                    dna = False
                    eiwit = "sequentie bestaat niet uit a, g, c of t"

            if dna:
                try:
                    for i in range(0, len(seq), 3):
                        eiwit += (code[seq[i:(i + 3)]])
                except IndexError:
                    eiwit = "geen goede dna sequentie"
    except TypeError:
        eiwit = "geef DNA sequentie -- typeError"
    return eiwit


if __name__ == '__main__':
    app.run()
