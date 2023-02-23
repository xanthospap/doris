#! /usr/bin/bash

fname=phd_report

shopt -s extglob
rm -v phd_report!(.tex) 2>/dev/null
shopt -u extglob

pdflatex $fname
biber $fname
pdflatex $fname
makeglossaries $fname
makeindex ${fname}.nlo -s nomencl.ist -o ${fname}.nls
pdflatex $fname
pdflatex $fname
