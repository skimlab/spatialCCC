#' Return LRdb table
#'
#' get_LRdb() returns the full size table while get_LRdb_small()
#'   returns a smaller table, useful for testing.
#'   See [LRdb_human] and [LRdb_mouse] for more detail about the database.
#'
#' @param species string, either "human" (default) or "mouse"
#' @param n_samples the number of LR-pairs to sample; if 0, return full table.
#'
#' @return LRdb table in data frame
#'
#'
#' @export
#'
#' @examples
#' LRdb <- get_LRdb(species = "human")
#'
get_LRdb <- function(species = c("human", "mouse"), n_samples = 0) {
  if (length(species) > 1) {
    species <- species[1]
  }

  species_option <- c("human", "mouse")

  if (length(species) != 1 ||
      is.na(pmatch(species, species_option))) {
    cli::cli_abort(message = "{.var species} should be a string, either 'human' or 'mouse'")
  }

  species <- species_option[pmatch(species, species_option)]

  if (species == "human")
    LRdb <- LRdb_human
  else
    LRdb <-LRdb_mouse

  if (n_samples > 0)
    LRdb <- dplyr::slice_sample(LRdb, n = n_samples)

  LRdb
}

#' Shuffle LRdb data frame
#'
#' Shuffle ligand and receptor pairing, utilized mainly to generate
#' NULL distribution of LRscore, which later can be used to compute
#' p-value.
#'
#' @param lrdb LRdb table in data frame
#' @param reverse LR pairing is reversed if TRUE.
#'        If FALSE, LR pairing is shuffled.
#'
#' @return LRdb table in data frame
#'
#' @export
shuffle_LRdb <- function(lrdb, reverse = TRUE) {
  if (reverse) {
    indices <- rev(1:nrow(lrdb))
  } else {
    indices <- sample(1:nrow(lrdb))
  }

  lrdb$receptor_gene_symbol <- lrdb$receptor_gene_symbol[indices]
  lrdb$receptor_gene_id <- lrdb$receptor_gene_id[indices]
  lrdb$receptor_ensembl_protein_id <- lrdb$receptor_ensembl_protein_id[indices]
  lrdb$receptor_ensembl_gene_id <- lrdb$receptor_ensembl_gene_id[indices]
  lrdb$evidence <- "random"

  lrdb %>%
    mutate(LR = paste(ligand_gene_symbol, receptor_gene_symbol, sep = "_"))
}

#'
#' @noRd
#'
get_LRdb_small <- function(species = c("human", "mouse")) {
  LRdb <-
    get_LRdb(species)

  LRdb_small <-
    LRdb %>%
    dplyr::slice_sample(n = 100)

  LR_pair_list <- c(
    "Apoe_Sdc1",
    "App_Dcc",
    "Cxcl12_Ackr3",
    "Jam3_Itgam",
    "Cx3cl1_Itga4")

  rbind(LRdb_small, dplyr::filter(LRdb, .data$LR %in% LR_pair_list)) %>%
    dplyr::distinct()
}
