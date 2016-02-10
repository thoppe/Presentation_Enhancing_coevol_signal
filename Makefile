name = biophys_poster
CLEAN_TARGETS = *.log *.aux *~ $(name).pdf *.bbl *.blg

all:
	aspell check -t $(name).tex
	make build
	make bibliography
	make build
	evince $(name).pdf


bibliography:
	bibtex $(name).aux

edit:
	emacs $(name).tex &

build:
	pdflatex $(name).tex

clean:
	rm -vf $(CLEAN_TARGETS)

commit:
	make all
	git commit -a 
