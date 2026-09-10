# AMPATH tuberculosis and pregnancy analysis

A single R workflow reconstructs the manuscript analyses from supplied AMPATH extracts: eligibility, pregnancy and TB definitions, cohort and self-controlled case series (SCCS) models, sensitivity analyses, pathway analyses, tables, and figures.

This repository contains code and documentation only. Data, calculated results, manuscript files, source-file fingerprints, fitted models, local configuration, and historical audit records are not distributed.

## Run the complete analysis

Requirements: R 4.5.3, the packages in `environment.yml`, a C compiler, and the `zip` command. The tested environment is Apple Silicon macOS. `environment-osx-arm64.lock.txt` records its exact package builds; use `environment.yml` to resolve an environment on another platform. The native bootstrap helper currently targets macOS/Linux.

With an existing micromamba installation:

```sh
micromamba create -n ampath-manuscript-reproduction -f environment.yml
micromamba activate ampath-manuscript-reproduction
```

Copy `config/paths.example.R` to `config/paths.local.R`. Set paths to four authorized extracts, kept outside the project:

| Configuration entry | Required input |
|---|---|
| `ampath_paths$demographics` | SAS demographics and treatment/history extract |
| `ampath_paths$pregnancies` | SAS pregnancy/nonpregnancy intervals |
| `ampath_paths$visits` | Longitudinal visits CSV |
| `ampath_outcome_path` | Pregnancy-outcome CSV with the same patient/date row order as visits |

Then run:

```sh
Rscript --vanilla run.R
```

This runs all modules in dependency order and checks the generated results. No saved results or study-value reference files are required. Production runs from R, including compilation of the bundled C helper; Python is not required.

The single entry point resolves paths from its own location. `Rscript --vanilla run.R --list` displays the module sequence. To resume a local analysis, use module numbers, for example `--from=14 --through=18`. Resuming requires the preceding local intermediates and manifests; a new installation must start at module 1.

## Project layout

| Path | Purpose |
|---|---|
| `run.R` | Complete analysis entry point |
| `R/` | Reusable analysis, presentation, and sharing functions |
| `scripts/` | Ordered analysis modules |
| `config/paths.example.R` | Placeholder input configuration |
| `docs/METHODS.md` | Analysis rules and code map |
| `docs/VARIABLES.md` | Main derived-variable definitions |
| `src/` | Bootstrap compatibility helper and third-party notices |
| `tests/check_code.R` | Syntax and synthetic checks; no study data |
| `tools/check_results.R` | Consistency checks on locally generated results |
| `tools/check_share.R` | Review of approved sharing files |
| `tools/export_code.R` | Code-only ZIP exporter |

Results are written locally under `outputs/`; participant-level intermediates and fitted models are under `work/`. Both entire directories are ignored by Git and excluded from code exports. The workflow produces three main tables, eleven supplementary tables, and four figures, plus unrounded statistics and model diagnostics. Editable tables are under `outputs/stage4_documents/`, individual tables under `outputs/stage4_tables/`, and figures under `outputs/stage4_figures/`.

Reproducibility starts with the supplied extracts. Recreating the earlier derivation of the pregnancy-interval extract additionally requires its original generating program and upstream inputs.

## Checks

```sh
Rscript --vanilla tests/check_code.R
Rscript --vanilla tools/check_share.R
```

After a full local analysis, `tools/check_results.R` verifies model convergence and rank, event and eligibility totals, input immutability, required figure files, and agreement between Word table cells and their CSV sources. These are consistency checks, not comparisons with distributed study results. The consolidated workflow was also rebuilt privately from the original extracts and compared with the previously reviewed analyses.

The GitHub workflow runs code parsing, synthetic tests, and the sharing check without access to any study extracts. It uses the documented [checkout](https://github.com/actions/checkout) and [R setup](https://github.com/r-lib/actions/tree/v2/setup-r) actions. Running these checks locally does not certify a future GitHub run.

## Share the code

```sh
Rscript --vanilla tools/export_code.R
```

The exporter writes `release/ampath-tb-pregnancy-code-only.zip`. It copies only the exact paths listed in `share-files.txt`, reruns the sharing check on the copied files, and verifies the ZIP inventory. Git history and local results are never copied. Adding a new source file requires review and an explicit addition to that list.

For GitHub, use this code-only project, or unpack that ZIP into a new repository. Do not upload the parent research workspace or an earlier results-inclusive archive. Git ignores all data/results folders, local configuration, common dataset and document formats, credentials files, and compiled libraries. The sharing check also rejects unapproved tracked files, including files added with Git's force option.

The public repository can be shared without granting access to the study extracts. Third-party licensing for the native helper is retained in `src/NUMPY_LICENSE.txt`; that notice applies to the adapted third-party code.
