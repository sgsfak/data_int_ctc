SERIES = 16443 18672 19783 22820 27562 31364 9574
GSE = $(patsubst %,data/gse%.rda,$(SERIES))

COMPARISONS = c1 c2 c3
PDFS = $(patsubst %, %.pdf, $(COMPARISONS))
MDS = $(patsubst %, %.md, $(COMPARISONS))


all: $(GSE) $(PDFS) results.pdf

$(GSE): R/get_data.R
$(GSE): download_files
	
.SECONDARY: download_files $(MDS)
.PHONY: clean cacheclean
download_files:
	Rscript --vanilla -e 'source("R/get_data.R");get_data()'

.SUFFIXES: .pdf .md
.md.pdf:
	pandoc -V geometry:margin=1in -s $< -o $@

%.md : %.Rmd
	Rscript -e "library(knitr); knit(\"$<\")"

%.html : %.Rmd
	Rscript -e "library(knitr); knit2html(\"$<\")"

clean:
	rm -rf $(PDFS)

cacheclean: clean
	rm -rf cache
