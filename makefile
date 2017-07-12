
DOC=essential-ising
FIGS=fig-susmag
DATA=data_square-4096
FUNCS=fig-mag_scaling-func fig-sus_scaling-func

all: ${DOC}.pdf

%.tex: figs/%.gplot ${DATA:%=data/%.dat} ${FUNCS:%=figs/%.dat}
	gnuplot $< > $@

${DOC}.pdf: ${DOC}.tex ${DOC}.bib ${FIGS:%=%.tex}
	rubber $(DOC).tex
	dvipdf $(DOC).dvi

clean:
	rubber --clean $(DOC)
	rm -f $(DOC).pdf
	rm -f $(DOC)Notes.bib
	rm -f ${FIGS:%=%.tex}

