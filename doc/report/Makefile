.PHONY: all clean

all: report.pdf

clean:
	rm report.pdf

report.pdf: *.tex *.bib
	pdflatex report
	bibtex report
	pdflatex report
	pdflatex report
