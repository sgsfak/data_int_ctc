## Data integration ##

This is a git repository for the code of our multi-platform data integration
work for the  identification of biomarkers with relation to Circulating Tumor
Cells (CTCs)

The R code can be found in the `R` directory. The `Makefile` can be used to
produce various documentation of the analysis that exists in the `.pdf` files.
You need [GNU make](http://www.gnu.org/software/make/) for that:
    
    make

and also [pandoc](http://johnmacfarlane.net/pandoc/) for creating the final pdf
files (through Markdown).

The code exists in the `R` folder and in the "source" `.Rmd` files. We are
using the excellent [knitr](http://yihui.name/knitr/) package for
transformation of the `Rmd` files to pdf. Therefore as the result of the `make`
command a number of PDF files are produced -- one for each `Rmd` file.

### Important note (13/06/2013)

There seems that one of the datasets
([GSE27562](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27562)) was
updated just a couple of days ago (in June 7th). The result of this is that now
the C2 comparison yields 103 overexpressed genes instead of 79 that was
initially (in the paper) described!

