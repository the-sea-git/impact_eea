# GEV Fitting Analysis

This R script performs Generalized Extreme Value (GEV) fitting analysis.

## Libraries
library(evd)  # Required for GEV functions

## Function to fit GEV
fit_gev <- function(data) {
  # Fit the GEV distribution
  gev_fit <- fgev(data)
  return(gev_fit)
}

## Example Data
set.seed(123)
example_data <- rgev(100, loc=0, scale=1, shape=0)

## Fit GEV to Example Data
gev_model <- fit_gev(example_data)

# Print GEV parameters
print(summary(gev_model))
