# csthelpers

## Overview
- R package for climate stress testing helper functions, weighting methods, and plotting utilities. Scope covers package source, tests, and lightweight development helpers in this repo; package data under `data/` and `inst/extdata/` is reference material and out of scope unless a task explicitly targets it.

## Repo Map
- `R/`: package source files
- `tests/testthat/`: unit tests
- `man/`: generated Rd files; update only when roxygen changes require regeneration
- `inst/extcode/`, `inst/examples/`: example scripts and exploratory usage
- `data/`, `inst/extdata/`: read-only reference/package data unless explicitly asked
- `TEMP/`: disposable scratch output for quick plots/checks; empty it before finishing
- `.agents/skills/`, `.claude/skills/`: externally managed skill links; do not edit

## Key Commands
```bash
# Tests
R -q -e "devtools::test()"

# Single test file during iteration
R -q -e "testthat::test_file('tests/testthat/test-file.R')"

# Lint
R -q -e "lintr::lint_package()"

# Regenerate docs when roxygen comments change
R -q -e "devtools::document()"
```

## Design philosophy

This is a **prototype-level repository**, not a production system. Optimize for clarity and traceability over generality.

- **Keep it simple.** Prefer flat, direct implementations. Avoid unnecessary abstraction layers, wrapper functions, config objects, registries, or factories.
- **No over-engineering.** Do not generalize for hypothetical future use cases. Solve the current problem with the minimum viable structure.
- **Readable by default.** Code should be easy to step through, explain to a collaborator, and modify without deep context. Favor explicit control flow over clever indirection.
- **Bias toward deletion.** Before adding a new layer, ask: can this be merged, flattened, or removed?
- **Function signatures stay lean.** Short argument lists, direct parameter passing, no nested config objects unless clearly justified.


## Conventions
- Keep edits surgical and follow existing patterns in nearby `R/` and `tests/testthat/` files.
- Do not change exported function names, arguments, defaults, or return structures unless explicitly requested.
- Do not add new package dependencies without explicit need.
- Keep files UTF-8.
- Prefer explicit base R control flow over `purrr`-style rewrites.
- Update tests when behavior changes; regenerate `man/` and `NAMESPACE` only when roxygen changes require it.

## Workflow
- Prefer targeted test files while iterating; run `devtools::test()` before finishing.
- Run `lintr::lint_package()` when changes touch style-sensitive or newly added R code.
- Do not reformat unrelated files or make broad refactors while handling focused tasks.
- If you create scratch plots or other temporary artifacts in `TEMP/`, remove them before finishing.
- Summarize files changed, commands run, results, and any unresolved issues at the end.

## Hard Constraints
- IMPORTANT: Do not modify `data/` or `inst/extdata/` unless the task explicitly requires it.
- IMPORTANT: Do not edit `.agents/skills/` or `.claude/skills/`; they are managed externally.
- Do not edit generated files such as `man/` or `NAMESPACE` unless the task requires regenerated documentation.
- Do not change CI, deployment, secrets, or other repo infrastructure unless explicitly requested.

## References
- Package metadata: `DESCRIPTION`
- Test entrypoint: `tests/testthat.R`
- Lint config: `.lintr.R`
- Repo-specific package guidance: `inst/csthelpers_instructions.txt`
