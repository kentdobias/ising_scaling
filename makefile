
DOC=essential-ising
FIGS=fig-susmag

all: ${FIGS:%=figs/%.tex}
	rubber $(DOC).tex
	dvipdf $(DOC).dvi

figs/%.tex: figs/%.gplot
	gnuplot $< > $@

clean:
	rubber --clean $(DOC)
	rm -f $(DOC).pdf
	rm -f $(DOC)Notes.bib
	rm -f ${FIGS:%=figs/%.tex}

