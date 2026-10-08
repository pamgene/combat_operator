# Effect of defect D1 on the batch spread: per batch, the SD of each peptide across the
# batch's samples, averaged over peptides (the spread ComBat's scale parameter acts on), for the raw
# data, the old one-step L/S correction (bug D1) and the new one (fixed), no reference.
#
#   docker run --rm -v "$PWD:/src" --entrypoint Rscript combat_operator:dev /src/dev/spread_per_batch.R [out_dir]
# Output: spread_per_batch.csv in out_dir (default dev/compare).

if (dir.exists("/src")) setwd("/src")
source("dev/mock_ctx.R")
source("R/functions.R")
source("R/pgcombat.R")
New = pgCombat
old_env = new.env(); eval(parse(text = readLines("dev/reference/pgcombat_shiny_ecfbdb7.R", warn = FALSE)), envir = old_env)
Old = old_env$pgCombat
if (!file.exists("dev/data/simulated.csv")) source("dev/simulate_data.R")

rows = list()
for (dsn in names(datasets)) {
  ds = datasets[[dsn]]
  ctx = mock_ctx(read_dataset(ds), ds)
  X = build_matrix(ctx$data$.ri, ctx$data$.ci, ctx$data$.y, nrow(ctx$rselect()), nrow(ctx$cselect()))
  b = factor(sample_batches(ctx$data$.ci, batch_labels(ctx$select(ctx$colors)), ncol(X)))
  sd_batch = function(M) sapply(levels(b), function(l) mean(apply(M[, b == l, drop = FALSE], 1, sd)))
  old = Old$new()$fit(X, b)$Xc
  new = New$new()$fit(X, b)$Xc
  rows[[length(rows) + 1]] = data.frame(data = dsn, batch = levels(b),
                                        sd_raw = sd_batch(X), sd_old_LS = sd_batch(old), sd_new_LS = sd_batch(new))
}
res = do.call(rbind, rows)
out_dir = if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "dev/compare"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(out_dir, "spread_per_batch.csv"), row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
