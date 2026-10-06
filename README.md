# ComBat Operator

> **In development.** This operator replaces the Shiny-based [`dascombat_shiny_fit_operator`](https://github.com/pamgene/dascombat_shiny_fit_operator). The current version reproduces the old operator's results exactly, including its known defects; see the [phase A report](docs/reports/phaseA_equivalence/report.md). The defects are being fixed next, and a full user guide (which setup and model to use when) will follow.

##### Description

Removes technical batch effects (e.g. between runs or barcodes) from normalized peptide data with ComBat. Returns the corrected values plus a PCA plot before and after the correction.

##### Input projection

| Crosstab | Content |
|---|---|
| y-axis | normalized values (log- or VSN-transformed, not S100) |
| rows | peptide IDs |
| columns | factors defining one sample per column (e.g. Barcode and Array/Row), plus the factors used by the fit condition (e.g. Grouping) |
| colors | technical batch factor(s), e.g. Run or Barcode; several factors are combined into one batch |

##### Settings

| Setting | Meaning |
|---|---|
| `ModelType` | `L/S` corrects the level and the spread of each batch, `L` only the level |
| `UseReferenceBatch`, `ReferenceBatch` | keep one batch unchanged and map the other batches onto it. With several colour factors give one value per factor, separated by `;` (e.g. `run01;650097104`) |
| `UseFitCondition`, `FitConditionFactors`, `FitConditionValues` | learn the correction only from the samples of a condition (e.g. `Grouping` = `REF`) and apply it to all samples. Several factors: `Supergroup;Test Condition` = `Control;DMSO`; several conditions: `Control;DMSO \| Control;Vehicle`. (These parameters replace the two-step shiny ComBat.) |

##### Output

| Output | Content |
|---|---|
| `CmbCor` | ComBat-corrected value per peptide and sample |
| PCA plot (PNG) | PC1 vs PC2 coloured by batch, before (left) and after (right) correction; fit-condition samples drawn as triangles |
