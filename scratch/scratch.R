# check null distirbution for LRscore
gexp <- as.matrix(SummarizedExperiment::assays(spe_brain)[["logcounts"]])
rowAnnots <- SingleCellExperiment::rowData(spe_brain)

q <- apply(gexp,2,quantile,0.99)
gg <- log(1+sweep(gexp,2,q/median(q),"/"))

ligand <- "App"
receptor <- "Lrp1"
lig_rec <- "App_Lrp1"

expression_min_prop <- 0.05


spot_dist <- calc_spot_dist(spe_brain)

n_cutoff <- ncol(gexp) * expression_min_prop
subset_idx <- rowSums(gexp > 0) > n_cutoff
gexp <- gexp[subset_idx, ]
rowAnnots <- subset(rowAnnots, subset = subset_idx)

# ligand and receptor expression values
lig_ids <- which(rowAnnots$symbol == ligand)
rec_ids <- which(rowAnnots$symbol == receptor)


if (length(lig_ids) > 1) {
  lig_exp <- colMeans(gexp[lig_ids, ])
} else {
  lig_exp <- gexp[lig_ids, ]
}

if (length(rec_ids) > 1) {
  rec_exp <- colMeans(gexp[rec_ids, ])
} else {
  rec_exp <- gexp[rec_ids, ]
}

sqrt.prod_all <- sqrt(lig_exp * rec_exp)

sqrt.prod <-
  sqrt(lig_exp[spot_dist$src] * rec_exp[spot_dist$dst])

