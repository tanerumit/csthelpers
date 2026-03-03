# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`csthelpers` is an R package providing utilities for climate stress testing (CST) and hydrological analysis. It contains functions for drought indices, environmental flow metrics, climate response surface visualization, and GCM scenario analysis.

## Build and Development Commands

```bash
# Check package (includes tests)
R CMD check .

# Build package
R CMD build .

# Install locally
R CMD INSTALL .

# Run tests
Rscript -e "testthat::test_local()"

# Run a single test file
Rscript -e "testthat::test_file('tests/testthat/test-compute_spei.R')"

# Generate documentation (roxygen2)
Rscript -e "devtools::document()"

# Load package during development
Rscript -e "devtools::load_all()"
```

## Architecture

### Core Function Categories

**Climate/Drought Indices:**
- `compute_spei()` - Wrapper around SPEI package for Standardized Precipitation Evapotranspiration Index
- `max_dry_spell()` - Maximum consecutive dry days using run-length encoding
- `compute_dpi()` - Drought Probability Index (internal helper in max_dry_spell.R)

**Environmental Flow (E-flow) Analysis:**
- `calculate_eflow_metrics()` - Computes 70 eco-hydrological indicators (HF/MF/LF/RFC/IF groups)
- `calculate_erfa_class()` - Environmental Flow Risk Assessment classification comparing baseline vs scenario

**Climate Response Surfaces:**
- `climateSurface()` - Generates 2D contour plots of system response to temperature/precipitation changes
- `climateSurfaceBatch()` - Batch version for multiple surfaces
- `GCMplausibilityRange()` - Computes plausibility ranges from GCM scenario ellipses

**Visualization:**
- `sos_radialplot()` - Faceted radial/spider plots for comparing indicators across scenarios
- `plot_basin_map()` - Basin maps with rivers and points using sf/ggplot2

**Utilities:**
- `findNN()` - Brute-force nearest neighbor search
- `select_maximin()` - Space-filling subset selection via maximin algorithm
- `find_sequences_below()` - Find consecutive sequences below threshold using RLE
- `slice_event_window()` - Extract fixed-length windows anchored to events
- `noleap_date_sequence()` - Generate date sequences excluding Feb 29

### Key Dependencies

The package heavily relies on tidyverse ecosystem (dplyr, ggplot2, tidyr, tibble, lubridate), plus sf for spatial operations and zoo for rolling statistics. Optional dependencies include SPEI, fields, and sp for specific functions.

### Code Patterns

- Functions use tidyverse-style NSE with `.data[[var]]` pronoun for column access
- Most visualization functions return ggplot2 objects for further customization
- E-flow metrics use run-length encoding (`rle()`) for sequence detection
- GCM functions rely on ggplot2 internals via `ggplot_build()` for ellipse extraction

## Testing

Tests use testthat (edition 3). Currently only `compute_spei` has tests. Tests skip when optional packages (like SPEI) are not installed.
