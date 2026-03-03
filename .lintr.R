linters <- lintr::linters_with_defaults(
  line_length_linter    = lintr::line_length_linter(120),
  indentation_linter    = lintr::indentation_linter(2),
  object_name_linter    = lintr::object_name_linter(styles = "snake_case"),
  commented_code_linter = lintr::commented_code_linter(),
  object_usage_linter   = lintr::object_usage_linter()
)

exclusions <- list(
  "renv/",
  "inst/extdata/",
  "docs/",
  "R/legacy_.*\\.R$"
)
