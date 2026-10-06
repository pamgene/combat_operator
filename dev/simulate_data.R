# Simulated test data with a known batch effect -> dev/data/simulated.csv
#
#   Rscript dev/simulate_data.R   (from the repo root)
#
# 150 peptides x 24 samples, 3 batches (A, B, C) of 8 samples.
# - every value ~ N(8, 1)
# - batch shift per peptide: s_p ~ N(0, 1); A + 0, B + 1 * s_p, C - 1 * s_p
# - batch C spread doubled around 8 (scale effect)
# - condition per batch: Control (A: 3, B: 2, C: 3), T1, T2 (the rest);
#   T1 raises peptides 1-20 by 1.5 (biological effect)
set.seed(1)
n_pep = 150
batch = rep(c("A", "B", "C"), each = 8)
cond = c(rep("Control", 3), rep("T1", 3), rep("T2", 2),
         rep("Control", 2), rep("T1", 3), rep("T2", 3),
         rep("Control", 3), rep("T1", 2), rep("T2", 3))
X = matrix(rnorm(n_pep * 24, 8, 1), n_pep)
X = X + outer(rnorm(n_pep), c(A = 0, B = 1, C = -1)[batch])
X[, batch == "C"] = 8 + (X[, batch == "C"] - 8) * 2
X[1:20, cond == "T1"] = X[1:20, cond == "T1"] + 1.5

d = data.frame(
  sim.peptide = rep(sprintf("pep%03d", 1:n_pep), 24),
  sim.sample = rep(sprintf("s%02d", 1:24), each = n_pep),
  sim.batch = rep(batch, each = n_pep),
  sim.condition = rep(cond, each = n_pep),
  sim.value = as.vector(X)
)
write.csv(d, "dev/data/simulated.csv", row.names = FALSE)
cat("wrote dev/data/simulated.csv:", nrow(d), "rows\n")
