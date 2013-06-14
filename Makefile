SERIES = 16443 18672 19783 22820 27562 31364 9574
GSE = $(patsubst %,data/gse%.rda,$(SERIES))

all: $(GSE) c1.pdf c2.pdf c3.pdf

$(GSE): R/get_data.R
$(GSE): download_files
	
.SECONDARY: download_files
download_files:
	Rscript --vanilla -e 'source("R/get_data.R");get_data()'

.SUFFIXES: .pdf .md
.md.pdf:
	pandoc -V geometry:margin=1in -s $< -o $@

%.md : %.Rmd
	Rscript -e "library(knitr); knit(\"$<\")"
