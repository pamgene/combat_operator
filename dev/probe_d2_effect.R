# Defect D2, effect on the results: the stop rule of it.sol in R/pgcombat.R
#
#   change <- max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old)   (conv 1e-4)
#
# ignores peptides with a negative g.old (negative ratio). Here the same fit is run
# with a strict rule: relative change on absolute values (floored at 1e-12 against
# 0/0 for the reference batch), conv 1e-10, at most 10000 iterations. Reported per
# case: max abs difference of the corrected values (strict - current).
#
#   docker run --rm -v "$PWD:/src" --entrypoint Rscript combat_operator:dev /src/dev/probe_d2_effect.R [out_dir]
# Output: d2_effect.csv in out_dir (default dev/compare).

if (dir.exists("/src")) setwd("/src")
source("dev/mock_ctx.R")
source("R/functions.R")
if (!file.exists("dev/data/simulated.csv")) source("dev/simulate_data.R")

src = readLines("R/pgcombat.R", warn = FALSE)
load_pgcombat = function(src) { e = new.env(); eval(parse(text = src), envir = e); e$pgCombat }
Current = load_pgcombat(src)
i = grep("max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old)", src, fixed = TRUE)
j = grep("conv =", src, fixed = TRUE)
k = grep("count <- count + 1", src, fixed = TRUE)
stopifnot(length(i) == 1, length(j) == 1, length(k) == 1)
strict = src
strict[i] = sub("max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old)",
                "max(abs(g.new - g.old) / pmax(abs(g.old), 1e-12), abs(d.new - d.old) / pmax(abs(d.old), 1e-12))",
                strict[i], fixed = TRUE)
strict[j] = sub("conv =\\s*\\.0001", "conv = 1e-10", strict[j])
strict[j] = sub("conv =\\s*$", "conv = 1e-10 + 0 *", strict[j])
strict[k] = paste(strict[k], '; if (count > 10000) stop("it.sol did not converge")')
Strict = load_pgcombat(strict)

rows = list()
for (dsn in names(datasets)) {
  ds = datasets[[dsn]]
  ctx = mock_ctx(read_dataset(ds), ds)
  X = build_matrix(ctx$data$.ri, ctx$data$.ci, ctx$data$.y, nrow(ctx$rselect()), nrow(ctx$cselect()))
  b = factor(sample_batches(ctx$data$.ci, batch_labels(ctx$select(ctx$colors)), ncol(X)))
  for (ref in list(NULL, levels(b)[1])) {
    a = Current$new()$fit(X, b, ref.batch = ref)
    s = Strict$new()$fit(X, b, ref.batch = ref)
    rows[[length(rows) + 1]] = data.frame(
      data = dsn, model = "L/S", reference = if (is.null(ref)) "-" else ref,
      max_abs_diff_CmbCor = max(abs(s$Xc - a$Xc)),
      max_abs_diff_gamma_star = max(abs(s$gamma.star - a$gamma.star)),
      max_abs_diff_delta_star = max(abs(s$delta.star - a$delta.star)),
      value_range = diff(range(X)))
  }
}
res = do.call(rbind, rows)
out_dir = if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "dev/compare"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(out_dir, "d2_effect.csv"), row.names = FALSE)
options(width = 200)
print(res, row.names = FALSE, digits = 3)
