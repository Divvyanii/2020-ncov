# Script to generate quarantine effect plots
# Run this after main_model.R has been executed

# Load plotting function
source("R/plot_quarantine_effect.R")

# Generate quick comparison plot (faster, single run per scenario)
# This compares different quarantine effectiveness levels: 0%, 30%, 50%, 70%
cat("Generating quick quarantine comparison plot...\n")
plot_quarantine_quick(
  quarantine_levels = c(0, 0.3, 0.5, 0.7),
  nn = 1e3,  # number of particles
  dt = 0.25,
  filename = "quarantine_effect_quick"
)

cat("\nDone! Plot saved to plots/quarantine_effect_quick.pdf\n")
cat("To generate a more detailed plot with bootstrap replicates, use:\n")
cat("  plot_quarantine_effect(quarantine_levels = c(0, 0.3, 0.5, 0.7), rep_plot = 20)\n")
cat("Note: This will take longer as it runs multiple bootstrap replicates.\n")

