## Makefile to use knitr in the vignette


all: knit 

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif
RPROG=	$(R_HOME)/bin/R
RSCRIPT=$(R_HOME)/bin/Rscript

knit:
	$(RSCRIPT) -e ".libPaths(c('~/R', .libPaths())); require(knitr); knit( 'knit_pmc_tutorial.Rnw', 'pmc_tutorial.Rnw');"
	$(RPROG) CMD Sweave pmc_tutorial.Rnw 
	pdflatex pmc_tutorial.tex
	bibtool -s -d -x pmc_tutorial.aux ~/Documents/Mendeley/bib/library.bib > pmc_tutorial.bib
	bibtex pmc_tutorial 
	pdflatex pmc_tutorial.tex
	pdflatex pmc_tutorial.tex
	rm -f *.bbl *.blg *.aux *.out *.log *.spl *.tikz 
	cp *.bib *.pdf pmc_tutorial.Rnw ../../../vignettes/

pdfclean:
	rm -f *.pdf

clean:
	rm -f index.html
	rm -f *.bbl *.blg *.aux *.out *.log *.spl *.tikz *tikzDictionary 
	rm -f *.nex

cleancache:
	rm -rf cache


