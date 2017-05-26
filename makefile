
DOC=essential-ising.tex
FIGS=scaling_func

all: ${FIGS:%=figs/%.tex}
	rubber --pdf $(DOC)

figs/%.tex: figs/%.gplot
	gnuplot $< > $@

clean:
	rubber --clean $(DOC)

