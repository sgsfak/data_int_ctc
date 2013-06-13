#------------------------------
# get_data: loads the data from GEO and builds a 'datasets' list
#------------------------------
require(GEOquery)
require(hgug4112a.db)
require(impute)
require(genefilter)
require(hgu133a.db)
require(ggplot2)

needs.log.tr <- function(gset) {
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
       (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    return (LogC)
}

log.tr <- function(gset) {
    if (needs.log.tr(gset)) {
        ex <- exprs(gset)
        ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex)
    }
    return (gset);
}

get_c1_data <- function() {
    gset<-getGEO(GEO="GSE22820", destdir="data")
    if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
    gse22820<-gset[[idx]]
    gse22820 <- log.tr(gse22820);
    pData(gse22820)$Study <- "GSE22820"
    
    #GPL6480 platform is Agilent-014850 Whole Human Genome Microarray 4x44K G4112F
    # See http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480
    ###gse22820@annotation <- 'hgug4112a' # See https://stat.ethz.ch/pipermail/bioconductor/2008-November/025062.html
    gse22820 <- nsFilter(gse22820, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)$eset
    names <- as.character(mget(featureNames(featureData(gse22820)), GPL6480ENTREZID))
    featureNames(featureData(gse22820)) <- names
    rownames(exprs(gse22820)) <- names
    
    save(gse22820, file="data/gse22820.rda")
    
    ## --------------------------------------------
    
    gset<-getGEO(GEO="GSE19783", destdir="data")
    if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
    gse19783<-gset[[idx]]
    pData(gse19783)$Study <- "GSE19783"
    
    ##gse19783@annotation <- 'hgug4112a' # See https://stat.ethz.ch/pipermail/bioconductor/2008-November/025062.html
    gse19783 <- nsFilter(gse19783, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)$eset
    names <- as.character(mget(featureNames(featureData(gse19783)), GPL6480ENTREZID))
    featureNames(featureData(gse19783)) <- names
    rownames(exprs(gse19783)) <- names
    
    save(gse19783, file="data/gse19783.rda")
    ## --------------------------------------------
    
    gset<-getGEO(GEO="GSE31364", destdir="data")
    gse31364<-gset[[1]]
    pData(gse31364)$Study = "GSE31364"
    
    ## Probably has missing values
    if (sum(is.na(exprs(gse31364))) > 0) {
        a <- impute.knn(exprs(gse31364))
        exprs(gse31364) <- a$data
    }
    
    # remove unknown probes, etc
    a <- nsFilter(gse31364, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)
    gse31364 <- a$eset
    names <- as.character(mget(featureNames(featureData(gse31364)), GPL14378ENTREZID))
    featureNames(featureData(gse31364)) <- names
    rownames(exprs(gse31364)) <- names
    
    save(gse31364, file="data/gse31364.rda")
    ## --------------------------------------------
    gset<-getGEO(GEO="GSE18672", destdir="data")
    gse18672 <- gset[[1]]
    # idx <- gse18672$characteristics_ch1.1 != 'clinical characteristic/disease status: breast tumor'
    # gse18672 <- gse18672[,idx]
    # gse18672$Disease <- 'Control'
    gse18672$Study <- "GSE18672"
    annotation(gse18672) <- "hgug4112a"
    
    ## Probably has missing values
    if (sum(is.na(exprs(gse18672))) > 0) {
        a <- impute.knn(exprs(gse18672))
        exprs(gse18672) <- a$data
    }
    
    # remove unknown probes, etc
    a <- nsFilter(gse18672, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)
    gse18672 <- a$eset
    names <- as.character(mget(featureNames(featureData(gse18672)), hgug4112aENTREZID))
    featureNames(featureData(gse18672)) <- names
    rownames(exprs(gse18672)) <- names
    
    save(gse18672, file="data/gse18672.rda")
    ## Let's see some graphs
    ## compare the normal samples from the two sets
    ## that have "control" cases
    
    #gse1_control <- as.vector(exprs(gse22820[,gse22820$Disease=='Control']))
    #gse4_control <- as.vector(exprs(gse18672[,gse18672$Disease=='Control']))
    #g <- c(rep.int("gse22820", length(gse1_control)), rep.int("gse18672", length(gse4_control)));
    #ggplot(data.frame(x=c(gse1_control, gse4_control), g=g), aes(x=x, fill=g)) + geom_density()
    
    
    
    ## --------------------------------------------
    gset<-getGEO(GEO="GSE33526", destdir="data")
    gse33526 <- gset[[1]];
    gse33526$bmi <- as.numeric(sub("bmi: ", "", gse33526$characteristics_ch2.1))
    ## gse33526 <- gse33526[,gse33526$bmi < 25]
    annotation(gse33526) <- "hgug4112a"
    
    
    #id_probe <- fData(gse33526)[,c("ID", "NAME")]
    #c<-id_probe[featureNames(gse33526),]
    
    ee <- aggregate(exprs(gse33526), list(NAME=fData(gse33526)$NAME), median)
    rownames(ee) <- ee$NAME
    ee$NAME <- NULL
    
    
    
    ## --------------------------------------------
    
    gset <- getGEO(GEO="GSE9574", destdir="data");
    gse9574 <- gset[[1]];
    ## idx <- gse9574$characteristics_ch1 == 'disease-state: normal (reduction mammoplasty)'
    ## gse9574 <- gse9574[,idx]
    ## gse9574$Disease <- 'Control'
    gse9574$Study <- "GSE9574"
    annotation(gse9574) <- "hgu133a"
    
    gse9574 <- log.tr(gse9574)
    
    a <- nsFilter(gse9574, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)
    gse9574 <- a$eset
    names <- as.character(mget(featureNames(featureData(gse9574)), hgu133aENTREZID))
    featureNames(featureData(gse9574)) <- names
    rownames(exprs(gse9574)) <- names
    
    
    save(gse9574, file="data/gse9574.rda")
}

get_c2_data <- function() {
    gset <- getGEO('GSE27562', destdir="data")
    
    gse1 <- gset[[1]]
    pData(gse1)$Disease <- ifelse(pData(gse1)$characteristics_ch1=='phenotype: Normal', "Control", "Cancer")
    gse1$Study <- 'GSE27562'
    annotation(gse1) <- "hgu133plus2"
    
    gse1 <- nsFilter(gse1, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)$eset
    names <- as.character(mget(featureNames(featureData(gse1)), hgu133plus2ENTREZID))
    featureNames(featureData(gse1)) <- names
    rownames(exprs(gse1)) <- names
    
    gse27562 <- gse1
    
    save(gse27562, file="data/gse27562_c2.rda")
    
    gset <- getGEO('GSE16443', destdir="data")
    eset <- gset[[1]]
    gse2 <- eset[!is.na(fData(eset)$GENE), ]
    
    
    exclude <- grepl("Breast-feeding|week of mens cycle|Pregnant", pData(gse2)$characteristics_ch1.7)
    gse2 <- gse2[, ! exclude]
    
    e <- exprs(gse2)
    i<-apply(e, 1, IQR)
    d <- data.frame(probe=rownames(e), gene=fData(gse2)$GENE, iqr = i)
    dp <- aggregate(iqr~ gene, d, max)
    dd <- base::merge(d, dp)
    
    ix <- featureNames(gse2) %in% dd$probe
    gse2 <- gse2[ix,]
    featureNames(featureData(gse2)) <- dd$gene
    rownames(exprs(gse2)) <- dd$gene
    gse2$Study <- "GSE16443"
    pData(gse2)$Disease <- ifelse(pData(gse2)$characteristics_ch1.1 == "status: Healthy", "Control", "Cancer")
    
    gse16443 <- gse2
    save(gse16443, file="data/gse16443.rda")
}

get_c3_data <- function() {
    gset <- getGEO("GSE15852", destdir="data")
    gse15852 <- log.tr(gset[[1]])
    gse15852$Origin <- 'Cancer Tissue'
    gse15852$Study <- 'GSE15852'
    annotation(gse15852) <- "hgu133a"
    # normal <- gse15852$characteristics_ch1.2 == 'grade: normal'
    # gse15852_can <- gse15852[,!normal]
    
    gsef <- nsFilter(gse15852, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)$eset
    names <- as.character(mget(featureNames(featureData(gsef)), hgu133aENTREZID))
    featureNames(featureData(gsef)) <- names
    rownames(exprs(gsef)) <- names
    gse15852 <- gsef
    
    save(gse15852, file="data/gse15852.rda")
    
    
    gset <- getGEO("GSE12763", destdir="data")
    gse12763 <- log.tr(gset[[1]])
    gse12763$Origin <- 'Cancer Tissue'
    gse12763$Study <- 'GSE12763'
    annotation(gse12763) <- "hgu133plus2"
    gsef <- nsFilter(gse12763, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)$eset
    names <- as.character(mget(featureNames(featureData(gsef)), hgu133plus2ENTREZID))
    featureNames(featureData(gsef)) <- names
    rownames(exprs(gsef)) <- names
    gse12763 <- gsef
    
    
    save(gse12763, file="data/gse12763.rda")
    
    gset <- getGEO('GSE27562', destdir="data")
    gse27562 <- gset[[1]]
    gse27562$Origin <- 'Normal Blood'
    gse27562$Study <- 'GSE27562'
    annotation(gse27562) <- "hgu133plus2"
    # normal <- pData(gse27562)$characteristics_ch1=='phenotype: Normal'
    # gse27562 = gse27562[,normal]
    gsef <- nsFilter(gse27562, require.entrez=TRUE,remove.dupEntrez=TRUE,
                 var.func=IQR, var.cutoff=0.5, var.filter=FALSE)$eset
    names <- as.character(mget(featureNames(featureData(gsef)), hgu133plus2ENTREZID))
    featureNames(featureData(gsef)) <- names
    rownames(exprs(gsef)) <- names
    gse27562 <- gsef
    
    save(gse27562, file="data/gse27562.rda")
}


get_data <- function() {
    get_c1_data()
    get_c2_data()
    get_c3_data()
}
