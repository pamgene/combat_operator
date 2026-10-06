# Probe of defect D2: the stop rule of the empirical-Bayes iteration (it.sol) in the
# old ComBat code (dev/reference/pgcombat_shiny_ecfbdb7.R, identical to sva::ComBat):
#
#   change <- max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old)
#
# The denominator has no abs(): a negative g.old gives a negative ratio, which max()
# ignores. A g.old of exactly 0 gives Inf (or NaN for 0/0).
#
#   docker run --rm -v "$PWD:/src" --entrypoint Rscript combat_operator:dev /src/dev/probe_d2_convergence.R [out_dir]
#
# Fits the L/S model with and without a reference batch on every data set and records,
# per batch: iterations, the stop-rule value after iteration 1, the smallest |gamma.hat|,
# and how many peptides have a negative gamma.hat (ignored by the stop rule).
# Output: d2_convergence.csv in out_dir (default dev/compare).

if (dir.exists("/src")) setwd("/src")
source("dev/mock_ctx.R")
if (!file.exists("dev/data/simulated.csv")) source("dev/simulate_data.R")

src = readLines("dev/reference/pgcombat_shiny_ecfbdb7.R", warn = FALSE)
i = grep("count <- count + 1", src, fixed = TRUE)
src = append(src, '        if (count == 1) assign("first_change", change, envir = .GlobalEnv)', i)
i = grep("adjust <- rbind(g.new, d.new)", src, fixed = TRUE)
src = append(src, paste0('      assign("log_it", rbind(get0("log_it", .GlobalEnv), data.frame(iterations = count, ',
                         'change_after_iteration_1 = get("first_change", .GlobalEnv), min_abs_gamma_hat = min(abs(g.hat)), ',
                         'n_negative_gamma_hat = sum(g.hat < 0), n_peptides = length(g.hat))), envir = .GlobalEnv)'), i - 1)
env = new.env(); eval(parse(text = src), envir = env)

rows = list()
for (dsn in names(datasets)) for (use_ref in c(FALSE, TRUE)) {
  ds = datasets[[dsn]]
  ctx = mock_ctx(read_dataset(ds), ds)
  X = with(ctx$data, tapply(.y, list(.ri, .ci), c))
  b = factor(with(ctx$data, tapply(ctx$data[[ds$colours]], .ci, `[`, 1)))
  ref = if (use_ref) levels(b)[1] else NULL
  log_it <<- NULL
  env$pgCombat$new()$fit(X, b, ref.batch = ref, mean.only = FALSE)
  rows[[length(rows) + 1]] = cbind(data = dsn, reference = if (use_ref) ref else "-",
                                   batch = levels(b), is_reference = levels(b) %in% ref, log_it)
}
res = do.call(rbind, rows)
out_dir = if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "dev/compare"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(out_dir, "d2_convergence.csv"), row.names = FALSE)
options(width = 200)
print(res, row.names = FALSE, digits = 4)
