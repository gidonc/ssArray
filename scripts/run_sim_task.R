# One task of the simulated-tables evaluation: simulate one table, fit every parameterisation to it.
#   Rscript scripts/run_sim_task.R 1          (task number; under SLURM it is taken from SLURM_ARRAY_TASK_ID)
# Set SIM_FAST=1 for a quick smoke test (100 + 100 draws).
suppressPackageStartupMessages({ library(cmdstanr); library(here) })
source(here("R", "single_area.R"))
source(here("R", "sim_tables.R"))

arg  <- commandArgs(TRUE)
task <- as.integer(if (length(arg)) arg[1] else Sys.getenv("SLURM_ARRAY_TASK_ID"))
fast <- nzchar(Sys.getenv("SIM_FAST"))
cfg  <- sim_design()[task, ]
stopifnot(nrow(cfg) == 1)
cat(sprintf("task %d: sigma_b %g, eps %g, rep %d\n", task, cfg$sigma_b, cfg$eps, cfg$rep))

mod <- cmdstan_model(here("stan", "single_area.stan"), include_paths = here("stan"))
set.seed(100000 + task)
sim <- sim_table(cfg$sigma_b, cfg$eps)

rows <- list(); ranks <- list()
for (nm in names(PARAMS)) {
  cat("  fitting", nm, "\n")
  r <- sim_fit_safe(mod, sim, param = PARAMS[[nm]], eps = cfg$eps, sigma_b = cfg$sigma_b, seed = task,
                    warmup = if (fast) 100 else 500, draws = if (fast) 100 else 500)
  rows[[nm]]  <- cbind(cfg, parameterisation = nm, r$row)
  ranks[[nm]] <- r$ranks
}
out <- here("results", "sim_tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(rows = do.call(rbind, rows), ranks = ranks, table = sim$T),
        file.path(out, sprintf("task_%04d.rds", task)))
cat("done\n")
