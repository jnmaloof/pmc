## Makefile to use knitr in the vignette


all: pmc_tutorial.pdf

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif
RPROG=	$(R_HOME)/bin/R
RSCRIPT=$(R_HOME)/bin/Rscript

pmc_tutorial.pdf: pmc_tutorial/pmc_tutorial.Rnw
	$(RPROG) CMD Sweave pmc_tutorial.Rnw 
	#Compile in latex, but keep aux file till we make bibliography
	$(RSCRIPT) -e "tools::texi2dvi( 'pmc_tutorial.tex', pdf = TRUE, clean = FALSE )"	
	bibtex pmc_tutorial # make the bibliography 
	$(RSCRIPT) -e "tools::texi2dvi( 'pmc_tutorial.tex', pdf = TRUE, clean = TRUE )" # compile and cleanup

knit:
	cd pmc_tutorial
	$(RSCRIPT) --vanilla -e ".libPaths(c('~/R', .libPaths())); require(knitr); knit( 'pmc_tutorial.Rnw', '../pmc_tutorial.Rnw');"
	cd ..
	$(RPROG) CMD Sweave pmc_tutorial.Rnw 
	#Compile in latex, but keep aux file till we make bibliography
	$(RSCRIPT) -e "tools::texi2dvi( 'pmc_tutorial.tex', pdf = TRUE, clean = FALSE )"
	bibtool -s -d -x pmc_tutorial.aux ~/Documents/Mendeley/bib/library.bib pmc_tutorial.bib
	bibtex pmc_tutorial # make the bibliography 
	$(RSCRIPT) -e "tools::texi2dvi( 'pmc_tutorial.tex', pdf = TRUE, clean = TRUE )" # compile and cleanup

pdfclean:
	rm -f *.pdf

clean:
	rm -f index.html
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl


