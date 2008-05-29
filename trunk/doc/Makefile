# Makefile to build PDF and html version of manual.
# This requires lyx and inkscape
DPI=300
LYX=lyx
PDFLATEX=pdflatex
BIBTEX=bibtex
MAKEIDX=makeindex
SVG2PNG=inkscape -D -d $(DPI) -e
SVG2EPS=inkscape -D -E 
SVG2PDF=inkscape -D -A
PNG2PDF=convert
LATEX2HTML=tth -i -e1 

%.png : %.svg
	$(SVG2PNG) $@ $<
%.pdf : %.svg
	$(SVG2PDF) $@ $<
%.eps : %.svg
	$(SVG2EPS) $@ $<

%.bbl : %.tex
	$(BIBTEX) $<
	
all: pdf html

pdf: manual.pdf

manual.pdf: manual.tex atom_numbering.png angle_defs.png dihedral_defs.png
	$(PDFLATEX) manual
	$(BIBTEX) manual
	$(MAKEIDX) manual
	$(PDFLATEX) manual
	$(PDFLATEX) manual
	
manual.tex: manual.lyx atom_numbering.eps angle_defs.eps dihedral_defs.eps
	$(LYX) --export pdflatex manual.lyx
	rm -f *.pdf rm *.bbl *.aux

# angle_defs.pdf:
# 	$(SVG2PNG) angle_defs.png angle_defs.svg
# 	$(PNG2PDF) angle_defs.png angle_defs.pdf
	
	
html: manual.tex atom_numbering.png angle_defs.png dihedral_defs.png
	$(LATEX2HTML) manual.tex

.PHONY: clean

clean:
	rm -f *.tex *.pdf *.png *.eps *.aux *.idx *.log *.out *.toc *.html *.tid *.tin *.tms *.blg *.bbl *.ind *.ilg
