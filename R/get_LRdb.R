#' Return LRdb table
#'
#' @param species string, either "human" (default) or "mouse"
#'
#' @return LRdb table in data frame
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
