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

c123_common_siggenes.txt: c1_siggenes.txt c2_siggenes.txt c3_siggenes.txt
	cat $^ | sort | uniq -c | awk '$$1==3{print $$2}' > $@

c124_common_siggenes.txt: c1_siggenes.txt c2_siggenes.txt c4_siggenes.txt
	cat $^ | sort | uniq -c | awk '$$1==3{print $$2}' > $@

union_c123_c124_siggenes.txt: c123_common_siggenes.txt c124_common_siggenes.txt
	cat $^ | sort | uniq > $@

clean:
	rm -rf $(PDFS)

cacheclean: clean
	rm -rf cache
