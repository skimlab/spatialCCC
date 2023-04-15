#' Return LRdb table
#'
#' get_LRdb() returns the full size table while get_LRdb_small()
#'   returns a smaller table, useful for testing.
#'   See [LRdb_human] and [LRdb_mouse] for more detail about the database.
#'
#' @param species string, either "human" (default) or "mouse"
#'
#'
#' @return LRdb table in data frame
#'
#'
#' @export
#'
#' @examples
#' LRdb <- get_LRdb(species = "human")
#'
get_LRdb <- function(species = c("human", "mouse")) {
  if (length(species) > 1) {
    species <- species[1]
  }

  species_option <- c("human", "mouse")

  if (length(species) != 1 ||
      is.na(pmatch(species, species_option))) {
    cli::cli_abort(message = "{.var species} should be a string, either 'human' or 'mouse'")
  }

  species <- species_option[pmatch(species, species_option)]

  if (species == "human") {
    LRdb_human
  } else {
    LRdb_mouse
  }
}

#'
#' @rdname get_LRdb
#'
#' @export
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

  rbind(LRdb_small, dplyr::filter(LRdb, LR %in% LR_pair_list)) %>%
    dplyr::distinct()
}
