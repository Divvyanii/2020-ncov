# Run analysis--------------------------------------------------------------



# - - -
# Run bootstrap SMC 
run_fits(rep_plot=100, # number of repeats
         nn=2e3, #number of particles
         dt=0.25,
         filename="1"
)

# Output plots
plot_outputs(filename="1") # Figure 2

# plot_dispersion(filename="1") # Figure 3

# - - -
# Plot quarantine effect (optional)
# Uncomment to generate quarantine comparison plots
# source("R/plot_quarantine_effect.R")
# plot_quarantine_quick(quarantine_levels = c(0, 0.3, 0.5, 0.7), nn=1e3, dt=0.25, filename="quarantine_effect")


# Run models --------------------------------------------------------------



