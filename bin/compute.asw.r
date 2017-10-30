library('optparse')
library(fpc)
option_list <- list(
 		make_option(c("-i", "--inputfile"), type="character", 
 			help="dissimilarity matrix"),
        make_option(c("-s", "--start"), type="numeric", 
 			help="num clusters start"),
        make_option(c("-e", "--end"), type="numeric", 
 			help="num clusters end")
        )

 	opts <- parse_args(OptionParser(option_list=option_list), 
 		args=commandArgs(trailing=TRUE))
    inputfp <- opts$inputfile
    kstart <- opts$start
    kend <- opts$end
dm <- read.table(opts$inputfile, header=T, sep="\t", row=1, check.names=F)
pk <- pamk(data=dm, krange=kstart:kend, diss=TRUE)
write.table(cbind(num.clusters=kstart:kend, asw=pk$crit[-1]), file=paste0(kstart, "_", kend, "_", inputfp), quote=FALSE, sep="\t", row.names=FALSE)