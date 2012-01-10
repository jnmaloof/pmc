## Makefile to use knitr in the vignette


all: pmc_tutorial.pdf

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif
RPROG=	$(R_HOME)/bin/R
RSCRIPT=$(R_HOME)/bin/Rscript

pmc_tutorial.pdf: pmc_tutorial/pmc_tutorial.Rnw
	cp -f pmc_tutorial/pmc_tutorial.Rnw .
	$(RSCRIPT) --vanilla -e "require(knitr); knit( 'pmc_tutorial.Rnw');"
	#Compile in latex, but keep aux file till we make bibliography
	$(RSCRIPT) -e "tools::texi2dvi( 'pmc_tutorial.tex', pdf = TRUE, clean = FALSE )" 
	bibtex pmc_tutorial # make the bibliography 
	$(RSCRIPT) -e "tools::texi2dvi( 'pmc_tutorial.tex', pdf = TRUE, clean = TRUE )" # compile and cleanup
#	cp pmc_tutorial/pmc_tutorial-fake.Rnw pmc_tutorial.Rnw

pdfclean:
	rm -f *.pdf

clean:
	rm -f index.html
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl


