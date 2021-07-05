savefile <- "/home/zack_ely/denovo_2/clpca_denov.Rdata"
save2file <- "/home/zack_ely/denovo_2/varinfo_denov.Rdata"
save3file <- "/home/zack_ely/denovo_2/pwpca_denov.Rdata"
library(scde)
library(htmltools)
pol <- read.table(file = "prep_scde.csv", row.names = 1, header=TRUE, sep=",")
cd <- clean.counts(pol, min.lib.size = 100)

knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = 8, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 8, plot = TRUE)
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

library(org.Mm.eg.db)
#translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Mm.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids 
# convert GO lists from ids to gene names
go.env <- eapply(org.Mm.egGO2ALLEGS, function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment

pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 8)
##df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)

clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 8, plot = TRUE)
#df <- pagoda.top.aspects(clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)

save(clpca, file=savefile)
save(varinfo, file=save2file)
save(pwpca, file=save3file)
#df <- pagoda.top.aspects(clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)

# get full info on the top aspects
#tam <- pagoda.top.aspects(clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
# determine overall cell clustering
#hc <- pagoda.cluster.cells(tam, varinfo)

#tamr <- pagoda.reduce.loading.redundancy(tam, clpca)

#tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)

#col.cols <- rbind(groups = cutree(hc, 3))
# compile a browsable app, showing top three clusters with the top color bar
#pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20))

#app <- make.pagoda.app(tamr2, tam, varinfo, go.env, clpca, col.cols = col.cols, cell.clustering = hc, title = "NPCs")

# show app in the browser (port 1468)
#show.app(app, "denovo", browse = TRUE, port = 1468) 

#save.image(file = "denovo.RData", version = NULL, ascii = FALSE, compress = !ascii, safe = TRUE)

#save.html(app, pagoda_go_output.html, background = "white", libdir = "data_denovo")

