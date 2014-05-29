#===================================================================
#       Setup and install annotation packages for the platforms used
# See http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/SQLForge.pdf
#===================================================================
require(AnnotationDbi)
require(AnnotationForge)
require(human.db0)

## Annotation for GPL14378
## 
tmpout = tempdir()
makeDBPackage("HUMANCHIP_DB",
        affy=FALSE,
        prefix="GPL14378",
        fileName="GPL14378.csv",
        baseMapType="gbNRef",
        outputDir = tmpout,
        version="1.0.0",
        manufacturer = "Agilent",
        chipName = "Agendia_human_DiscoverPrint_v1",
        manufacturerUrl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14378")
install.packages(paste(tmpout, "GPL14378.db", sep="/"), repos=NULL, type="source")

## Annotation for GPL6480
##
tmpout = tempdir()
makeDBPackage("HUMANCHIP_DB",
        affy=FALSE,
        prefix="GPL6480",
        fileName="GPL6480.csv",
        baseMapType="gb",
        outputDir = tmpout,
        version="1.0.0",
        manufacturer = "Agilent",
        chipName = "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)",
        manufacturerUrl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480")

install.packages(paste(tmpout, "GPL6480.db", sep="/"), repos=NULL, type="source")

#----- Old stuff
createGPL14378AnnDbi <- function() {
    # Try to create an SQLForge annot file
    # See http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/SQLForge.pdf
    require(AnnotationDbi)
    tmpout = tempdir()
    myMeta = c("DBSCHEMA"="HUMANCHIP_DB",
               "ORGANISM"="Homo sapiens",
               "SPECIES"="Human",
               "MANUFACTURER"="Agilent",
               "CHIPNAME"="Agendia human DiscoverPrint v1",
               "MANUFACTURERURL"="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14378")
    populateDB("HUMANCHIP_DB", affy = FALSE, prefix = "GPL14378",
               fileName = "GPL14378.csv", metaDataSrc = myMeta,
               baseMapType = "gbNRef", outputDir = tmpout)

    seed <- new("AnnDbPkgSeed",
                Package = "GPL14378.db",
                Version = "1.0.0",
                PkgTemplate = "HUMANCHIP.DB",
                AnnObjPrefix = "GPL14378")
    makeAnnDbPkg(seed,
                 file.path(tmpout, "GPL14378.sqlite"),
                 dest_dir = tmpout)
     makeDBPackage("HUMANCHIP_DB",
                    affy=FALSE,
                    prefix="GPL14378",
                    fileName="GPL14378.csv",
                    baseMapType="gbNRef",
                    outputDir = tmpout,
                    version="1.0.0",
                    manufacturer = "Agilent",
                    chipName = "Agendia_human_DiscoverPrint_v1",
                    manufacturerUrl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14378")
    ## 
}

    ## makeDBPackage("HUMANCHIP_DB",
    ##                affy=FALSE,
    ##                prefix="GPL6480",
    ##                fileName="GPL6480.csv",
    ##                baseMapType="gb",
    ##                outputDir = tmpout,
    ##                version="1.0.0",
    ##                manufacturer = "Agilent",
    ##                chipName = " 	Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)",
    ##                manufacturerUrl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480")
