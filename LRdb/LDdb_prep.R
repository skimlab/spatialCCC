LRdb_human <-
  readr::read_tsv("CellTalkDB/human_lr_pair.txt", show_col_types = FALSE) %>%
  dplyr::rename(LR = "lr_pair")

LRdb_mouse <-
  readr::read_tsv("CellTalkDB/mouse_lr_pair.txt", show_col_types = FALSE) %>%
  dplyr::rename(LR = "lr_pair")
