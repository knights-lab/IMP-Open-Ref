require(reshape2) 
require(optparse)

option_list <- list(
 		make_option(c("-i", "--inputfile"), type="character", 
 			help="blast format file"),
        make_option(c("-o", "--outputfile"), type="character", 
 			help="output distance matrix")
        )

 	opts <- parse_args(OptionParser(option_list=option_list), 
 		args=commandArgs(trailing=TRUE))

data <- read.table(opts$inputfile, header=F, sep="\t")
colnames(data) <- c("query.label", "target.label", "percent.identity", "alignment.length", "num.mismatches", "num.gaps", "query.start", "query.end", "target.start", "target.end")
n <- max(table(data[,2]))

# pick the min % identity or else it'll result in non-symmetric matrix
d1 <- data[with(data, order(query.label, target.label)), 1:3]
d2 <- data[with(data, order(target.label, query.label)), 1:3]
dd <- cbind(d1[,1:2], p1=d1$percent.identity,p2=d2$percent.identity)
dd <- transform(dd, min = pmin(p1,p2))

dm <- dcast(dd, query.label~target.label, value.var="min")
rownames(dm) <- dm[,1]
dm <- dm[,-1]
dm <- (100-dm)/100 # convert percent id to dissimilarity
dm[is.na(dm)] <- 1 # set NAs as being most distant

cat("\t", file=opts$outputfile)
write.table(dm, file=opts$outputfile, quote=F, sep="\t", append=T)



