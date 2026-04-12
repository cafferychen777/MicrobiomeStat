# MicrobiomeStat

MicrobiomeStat is an R package for statistical analysis and visualization of microbiome and multi-omics data.

## Repository Layout

- `R/`: package implementation.
- `man/`: generated Rd documentation.
- `tests/testthat/`: regression and unit tests.
- `vignettes/`: end-user long-form documentation.
- `data/`: packaged example datasets.
- `inst/`: installed assets, templates, external example files, and pkgdown source assets.
- `dev/`: developer-only scripts and benchmarks that are tracked in Git but excluded from package builds.
- `docs/`: generated pkgdown site committed for GitHub Pages.

## Organizing Principles

- Keep package-shipping content in standard R package directories only.
- Keep reusable developer automation and reproducibility helpers under `dev/`.
- Keep large benchmark inputs outside the Git repository in the sibling workspace directory `../datasets/`.
- Keep throwaway local exploration outside the package tree in `../scratch/`.

## Benchmarks

The Hackathon benchmark harness lives in `dev/hackathon/`. Those scripts resolve paths relative to themselves and expect the following workspace layout:

```text
MicrobiomeStat/
├── repo/
└── datasets/
    └── benchmark_data/
        ├── Hackathon/
        └── results/
```

See `dev/README.md` and `dev/hackathon/README.md` for the developer-facing workflow.
