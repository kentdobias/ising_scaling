
DOC=essential-ising
FIGS=fig-sus fig-mag

all: ${FIGS:%=figs/%.tex}
	rubber $(DOC).tex
	dvipdf $(DOC).dvi

figs/%.tex: figs/%.gplot
	gnuplot $< > $@

clean:
	rubber --clean $(DOC)

