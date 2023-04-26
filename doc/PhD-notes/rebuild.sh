#! /bin/bash
set -e

shopt -s extglob
rm phd_report.!(tex)

fn=phd_report
ftex=${fn}.tex
pdflatex $ftex
makeglossaries $fn
pdflatex $ftex
bibtex $fn
pdflatex $ftex
