
DOC=essential-ising
FIGS=fig-susmag
DATA=data_square-4096
FUNCS=fig-mag_scaling-func fig-sus_scaling-func

all: ${FIGS:%=figs/%.tex}
	rubber $(DOC).tex
	dvipdf $(DOC).dvi

figs/%.tex: figs/%.gplot ${DATA:%=data/%.dat} ${FUNCS:%=figs/%.dat}
	gnuplot $< > $@

clean:
	rubber --clean $(DOC)
	rm -f $(DOC).pdf
	rm -f $(DOC)Notes.bib
	rm -f ${FIGS:%=figs/%.tex}

