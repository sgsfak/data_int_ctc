library(siggenes)
library(inSilicoMerging)


load("data/gse12763.rda")
gse12763$Study = "GSE12763"
load("data/gse15852.rda")
gse15852$Study = "GSE15852"
load("data/gse27562.rda")
gse27562$Study = "GSE27562"
load("data/gse16443.rda")
gse16443$Origin = "Normal Blood"

removeBatch <- TRUE
if (removeBatch) {
    km = kmeans(t(exprs(gse15852)), 2)
    gse15852_cl1 = gse15852[,km$cluster==1]
    gse15852_cl2 = gse15852[,km$cluster==2]
    gse15852_m = merge(list(gse15852_cl2, gse15852_cl1), method="COMBAT")
    gse15852 = gse15852_m
}

load("data/gse22820.rda")
load("data/gse19783.rda")
load("data/gse31364.rda")
gse31364$Origin = "Cancer Tissue"
gse22820$Origin = "Cancer Tissue"
gse19783$Origin = "Cancer Tissue"

gse27562$Origin = 'Normal Blood'
gse22820$Disease ='Cancer'
d=  list(gse15852=gse15852, gse12763 = gse12763, gse27562=gse27562, gse22820=gse22820, gse31364=gse31364, gse19783=gse19783, gse16443=gse16443)
mall = merge(d, method="COMBAT")

# -----------------------------
normal = gse15852$characteristics_ch1.2 == 'grade: normal'
s_gse15852 = sampleNames(gse15852)[!normal]
s_gse12763 = sampleNames(gse12763)
normal <- pData(gse27562)$characteristics_ch1=='phenotype: Normal'
s_gse27562 = sampleNames(gse27562)[normal]

idx = gse19783$characteristics_ch1.1 == 'disease state: Breast Primary Tumor'
gse19783_samples = sampleNames(gse19783[,idx])
gse22820_samples = sampleNames(gse22820[,gse22820$characteristics_ch1!='disease state: normal breast tissue'])
gse31364_samples = sampleNames(gse31364)

gse16443_samples = sampleNames(gse16443[,gse16443$Disease == 'Control'])

## The sample names to use for the comparison
s = c(s_gse15852, s_gse12763, s_gse27562, gse22820_samples, gse31364_samples, gse19783_samples, gse16443_samples)
c3.data=mall[,s]



## alternative/old approach: prefilter and then merge:
# d1=  list(gse15852=gse15852[,s_gse15852], gse12763 = gse12763[,s_gse12763], gse27562=gse27562[,s_gse27562])
# m = merge(d1, method="COMBAT")

# ------ FIND DIFF EXPRESIION ----------
findDiffExp = FALSE
if (findDiffExp) {
    cl = ifelse(c3.data$Origin == 'Normal Blood', 0, 1)
    sam.out = sam(exprs(c3.data), cl, B=500, rand=0xDEAD)
    delta = findDelta(sam.out, fdr=0.05)[2,1]
    plot(sam.out, delta)

    sam.sum <- summary(sam.out, delta)
    w <- which(sam.sum@mat.sig$d.value > 0)
    num.genes.over <- length(w)
    siggenes.over <- list.siggenes(sam.out, delta)[w]
}




