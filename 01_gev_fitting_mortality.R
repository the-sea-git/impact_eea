# GEV Fitting for Extreme Event Attribution 

# This R script performs Generalized Extreme Value (GEV) fitting analysis (including bootstrapping to obtain uncertainty ranges)

library(extRemes)

# ---------------------------------------------------------
# Configuration for an Example (Two simulations out of 1000)
# ---------------------------------------------------------

set.seed(123)

YEA1_GEV <- 1950
YEA2_GEV <- 2022

nYEA_GEV <- YEA2_GEV - YEA1_GEV + 1
vYEA <- YEA1_GEV:YEA2_GEV

iYEA1_GEV <- 1
iYEA2_GEV <- nYEA_GEV

nSIM <- 2

# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

# Use project-relative paths
load("data/blockmaxima_mortality.RData")

# Choose main variable from loaded file (attributable fraction)
BLOCK_MAX <- amax_AF

# Regions
vREG <- rownames(BLOCK_MAX)
nREG <- length(vREG)

# ---------------------------------------------------------
# Load GMST
# ---------------------------------------------------------

TS_GMST_raw <- read.csv(
  file = "data/gmst_noaa.csv",
  header = TRUE,
  skip = 2
)

# Restrict years
gmst_subset <- TS_GMST_raw[
  TS_GMST_raw$Year >= YEA1_GEV &
    TS_GMST_raw$Year <= YEA2_GEV,
]

# Baseline adjustment (1880–1910)
baseline <- mean(
  TS_GMST_raw$Lowess.5[
    TS_GMST_raw$Year >= 1880 &
      TS_GMST_raw$Year <= 1910
  ],
  na.rm = TRUE
)

TS_GMST <- gmst_subset$Lowess.5 - baseline

# ---------------------------------------------------------
# Output arrays
# ---------------------------------------------------------

loc0_gev  <- array(NA, dim = c(nREG, nSIM + 1))
loc1_gev  <- array(NA, dim = c(nREG, nSIM + 1))
scale_gev <- array(NA, dim = c(nREG, nSIM + 1))
shape_gev <- array(NA, dim = c(nREG, nSIM + 1))

p_2022  <- array(NA, dim = c(nREG, nSIM + 1))
p_1800  <- array(NA, dim = c(nREG, nSIM + 1))
pr_2022 <- array(NA, dim = c(nREG, nSIM + 1))

# ---------------------------------------------------------
# Main loop
# ---------------------------------------------------------

for (iREG in seq_len(nREG)) {

  message(
    sprintf(
      "Region %s/%s: %s",
      iREG,
      nREG,
      vREG[iREG]
    )
  )

  for (iSIM in 0:nSIM) {

    # -----------------------------------------------------
    # Prepare data
    # -----------------------------------------------------

    GEV_DATA <- data.frame(
      MAXI = BLOCK_MAX[iREG, iYEA1_GEV:iYEA2_GEV, 1 + iSIM],
      GMST = TS_GMST[iYEA1_GEV:iYEA2_GEV]
    )

    # Bootstrap sample
    if (iSIM == 0) {

      GEV_DATA_BOOT <- GEV_DATA

    } else {

      GEV_DATA_BOOT <- GEV_DATA[
        sample(nYEA_GEV, nYEA_GEV, replace = TRUE),
      ]
    }

    # Skip problematic regions
    if (
      median(GEV_DATA_BOOT$MAXI, na.rm = TRUE) <= 0.5 ||
      any(is.na(GEV_DATA_BOOT$MAXI))
    ) {
      next
    }

    # -----------------------------------------------------
    # Fit GEV model
    # -----------------------------------------------------

    fit_gmle <- tryCatch(

      fevd(
        MAXI,
        data = GEV_DATA_BOOT,
        type = "GEV",
        method = "GMLE",
        location.fun = ~ GMST,
        use.phi = TRUE
      ),

      error = function(e) {
        message(
          sprintf(
            "Fit failed for region %s, sim %s",
            vREG[iREG],
            iSIM
          )
        )
        return(NULL)
      }
    )

    # Skip failed fits
    if (is.null(fit_gmle)) {
      next
    }

    # -----------------------------------------------------
    # Store parameters
    # -----------------------------------------------------

    loc0_gev[iREG, iSIM + 1]  <- fit_gmle$results$par[1]
    loc1_gev[iREG, iSIM + 1]  <- fit_gmle$results$par[2]
    scale_gev[iREG, iSIM + 1] <- fit_gmle$results$par[3]
    shape_gev[iREG, iSIM + 1] <- fit_gmle$results$par[4]

   # -----------------------------------------------------
# Probabilities for 2022 / current warming level
# -----------------------------------------------------

# Index corresponding to year 2022
i2022 <- which(vYEA == 2022)

# Event value in 2022
event_value <- GEV_DATA$MAXI[i2022]

# GMST value in 2022
gmst_2022 <- TS_GMST[i2022]

p_2022[iREG, iSIM + 1] <- 1 - pevd(
  event_value,
  loc = loc0_gev[iREG, iSIM + 1] +
    loc1_gev[iREG, iSIM + 1] * gmst_2022,
  scale = exp(scale_gev[iREG, iSIM + 1]),
  shape = shape_gev[iREG, iSIM + 1]
)

p_1800[iREG, iSIM + 1] <- 1 - pevd(
  event_value,
  loc = loc0_gev[iREG, iSIM + 1],
  scale = exp(scale_gev[iREG, iSIM + 1]),
  shape = shape_gev[iREG, iSIM + 1]
)

pr_2022[iREG, iSIM + 1] <-
  p_2022[iREG, iSIM + 1] /
  p_1800[iREG, iSIM + 1]

# ---------------------------------------------------------
# Save outputs
# ---------------------------------------------------------

save(
  vREG,
  loc0_gev,
  loc1_gev,
  scale_gev,
  shape_gev,
  p_1800,
  pr_2022,
  p_2022,
  file = "output/output_mort.RData"
)
