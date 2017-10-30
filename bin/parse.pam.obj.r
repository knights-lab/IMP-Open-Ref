# used locally only for loading pam object that was stored when doing clustering on teraminx

require('cluster')

setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/Embalmer-Open-Ref")

load(file = "pam.111.obj")

medoid.ix <- pam.obj$id.med # INDICES of medoid representatives

medoid.ids <- rownames(pam.obj$data)[medoid.ix]

write.table(medoid.ids, "medoid.ids.txt", row.names=F, quote=F, col.names=F)