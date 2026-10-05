# Combine the per-task results into one data set.   Rscript scripts/combine_sim.R
library(here)
files <- list.files(here("results", "sim_tables"), pattern = "^task_.*\\.rds$", full.names = TRUE)
res   <- lapply(files, readRDS)
dat   <- do.call(rbind, lapply(res, `[[`, "rows"))
rownames(dat) <- NULL
write.csv(dat, here("results", "sim_tables.csv"), row.names = FALSE)
# ranks of the true log cells, one row per (task, parameterisation, cell)
rk <- do.call(rbind, lapply(res, function(x) {
  do.call(rbind, lapply(names(x$ranks), function(nm)
    data.frame(task = x$rows$task[1], parameterisation = nm, cell = seq_along(x$ranks[[nm]]), rank = x$ranks[[nm]])))
}))
saveRDS(rk, here("results", "sim_tables_ranks.rds"))
cat(length(files), "tasks;", nrow(dat), "fits;", sum(!is.na(dat$error)), "errors\n")
