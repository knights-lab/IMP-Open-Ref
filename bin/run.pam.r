library('optparse')
library('cluster')

option_list <- list(
 		make_option(c("-i", "--inputfile"), type="character", 
 			help="dissimilarity matrix"),
        make_option(c("-k", "--k"), type="numeric", 
 			help="num clusters"),
        make_option(c("-o", "--outputfile"), type="character", 
 			help="file to save pam object")
        )

 	opts <- parse_args(OptionParser(option_list=option_list), 
 		args=commandArgs(trailing=TRUE))
    inputfp <- opts$inputfile
    k <- opts$k
    outputfp <- opts$outputfile
dm <- read.table(opts$inputfile, header=T, sep="\t", row=1, check.names=F)

pam.obj <- pam(dm, k)

save(pam.obj, file = outputfp) # save this just in case

medoid.ix <- pam.obj$id.med # INDICES of medoid representatives

medoid.ids <- rownames(pam.obj$data)[medoid.ix]

write.table(medoid.ids, "medoid.ids.txt", row.names=F, quote=F, col.names=F)