# ComBat operator — domain glossary

| Term | Meaning |
|---|---|
| Batch | Technical grouping to correct for (e.g. Run, Barcode) = combination of all colour factors, written `value1;value2` in colour order |
| Sample | One crosstab column, e.g. one array (Barcode + Array/Row) |
| Reference batch | One batch whose values stay unchanged; the other batches are mapped onto it (ComBat `ref.batch`) |
| Fit condition | Samples of one or more conditions (defined by column factors, e.g. Grouping = REF); the Combat model is fitted on these only and applied to all samples |
| Fit samples | Samples used to fit the model: all samples, or the fit condition |
| Combat model | Internal `pgCombat` fit: per-peptide L, S and per-batch gamma*, delta* |
| Fit / apply | Fit = estimate the Combat model from the fit samples; apply = correct samples with an estimated model |
| CmbCor | Combat-corrected value per peptide × sample (main output) |
| L / L/S model | Mean-only (level) / mean + scale (level and spread) correction |
| One-step / two-step (old Shiny operator) | One-step: fit and correct in one run. Two-step: save the model in one step, apply it in another; replaced by the fit condition |
