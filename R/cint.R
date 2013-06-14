require(org.Hs.eg.db)

c1 <- read.table("c1_siggenes.txt",col.names="g")$g
c2 <- read.table("c2_siggenes.txt",col.names="g")$g
c3 <- read.table("c3_siggenes.txt",col.names="g")$g

common <- intersect(c1, intersect(c2, c3))
cl <- mget(as.character(common), org.Hs.egSYMBOL)

# df <- data.frame(do.call(rbind, cl))
# df

m <- sapply(cl, function(x) x[[1]])
df <- data.frame(ezgene=names(m), symbol=m)
# df <- data.frame(symbol=sapply(cl, function(x) x[[1]]))
write.table(df, file="c1c2c3_commongenes.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
# print(df, row.names=FALSE)
