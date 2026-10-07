# Export New Zealand district cross-vote tables from ei.Datasets, collapsed to 5 x 5:
# rows = party vote, columns = party of the candidate voted for; National, Labour, Green, NZ First, rest (incl. informal).
# Usage: Rscript scripts/export_nz.R <path to ei_NZ_YYYY.RData> <year> <out csv>
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
load(args[1]); x <- get(paste0("ei_NZ_", args[2]))
grp <- function(s) case_when(
  str_detect(s, regex("National Party", ignore_case = TRUE)) ~ "National",
  str_detect(s, regex("Labour Party", ignore_case = TRUE)) ~ "Labour",
  str_detect(s, regex("Green Party", ignore_case = TRUE)) ~ "Green",
  str_detect(s, regex("New Zealand First", ignore_case = TRUE)) ~ "NZ First",
  TRUE ~ "rest")
lev <- c("National", "Labour", "Green", "NZ First", "rest")
out <- map_dfr(seq_len(nrow(x)), function(i) {
  m <- x$District_cross_votes[[i]]
  m |> rename(party_vote = 1) |>
    pivot_longer(-party_vote, names_to = "cand", values_to = "votes") |>
    mutate(row = grp(party_vote),
           col = grp(coalesce(str_match(cand, "\\(([^()]*(\\([^()]*\\))?[^()]*)\\)\\s*$")[, 2], "rest"))) |>
    group_by(row, col) |> summarise(votes = sum(votes), .groups = "drop") |>
    complete(row = lev, col = lev, fill = list(votes = 0)) |>
    mutate(area = x$Number_of_district[i], district = x$District[i],
           row_no = match(row, lev), col_no = match(col, lev))
}) |> select(area, district, row, col, row_no, col_no, votes) |> arrange(area, row_no, col_no)
write_csv(out, args[3])
cat(n_distinct(out$area), "districts, total", sum(out$votes), "\n")
