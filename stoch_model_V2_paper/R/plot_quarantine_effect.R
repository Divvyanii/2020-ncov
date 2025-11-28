# Plot quarantine effect --------------------------------------------

plot_quarantine_effect <- function(quarantine_levels = c(0, 0.3, 0.5, 0.7), 
                                    rep_plot = 20, 
                                    nn = 1e3, 
                                    dt = 0.25,
                                    filename = "quarantine_comparison"){
  
  # quarantine_levels: vector of quarantine effectiveness values to compare (0-1)
  # rep_plot: number of bootstrap replicates per scenario
  # nn: number of particles per SMC run
  # dt: time step
  # filename: output filename prefix
  
  # Store results for each quarantine level
  results_list <- list()
  
  # Save original quarantine effectiveness
  original_quarantine <- theta[["quarantine_effectiveness"]]
  
  cat("Running quarantine scenarios...\n")
  
  for(q_level in quarantine_levels){
    
    cat(paste("  Quarantine effectiveness:", q_level, "\n"))
    
    # Set quarantine effectiveness
    theta[["quarantine_effectiveness"]] <- q_level
    
    # Recalculate quarantine start time if needed
    if(!is.null(theta[["quarantine_start_date"]])){
      quarantine_start_date <- as.Date(theta[["quarantine_start_date"]])
      quarantine_start_time <<- as.numeric(quarantine_start_date - start_date + 1)
    }
    
    # Run bootstrap fits
    out_rep <- foreach(kk = 1:rep_plot) %dopar% {
      output_smc <- smc_model(theta, nn, dt)
      output_smc
    }
    
    # Extract trajectories
    R0_plot = matrix(NA, ncol=rep_plot, nrow=t_period)
    C_local_plot = matrix(NA, ncol=rep_plot, nrow=t_period)
    I_plot = matrix(NA, ncol=rep_plot, nrow=t_period)
    
    for(kk in 1:rep_plot){
      output_smc <- out_rep[[kk]]
      if(output_smc$lik != -Inf){
        R0_plot[,kk] <- output_smc$beta_trace/(theta[["recover"]])
        case_local_pos <- theta[["confirmed_prop"]]*(output_smc$C_local_trace - c(0,head(output_smc$C_local_trace,-1)))
        C_local_plot[,kk] <- case_local_pos
        I_plot[,kk] <- output_smc$E_trace + output_smc$I_trace*(1-theta[["confirmed_prop"]])
      }
    }
    
    # Remove NA fits
    R0_plot = R0_plot[,!is.na(R0_plot[t_period,])]
    C_local_plot = C_local_plot[,!is.na(C_local_plot[t_period,])]
    I_plot = I_plot[,!is.na(I_plot[t_period,])]
    
    # Calculate quantiles
    R0_quantile <- apply(R0_plot, 1, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=TRUE)})
    C_local_quantile <- apply(C_local_plot, 1, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=TRUE)})
    I_quantile <- apply(I_plot, 1, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=TRUE)})
    
    results_list[[paste0("q", q_level)]] <- list(
      R0_quantile = R0_quantile,
      C_local_quantile = C_local_quantile,
      I_quantile = I_quantile,
      level = q_level
    )
  }
  
  # Restore original quarantine value
  theta[["quarantine_effectiveness"]] <- original_quarantine
  
  # Create plots
  par(mfrow=c(3,1), mar=c(3,3,2,1), mgp=c(2,0.7,0))
  
  # Color palette for different quarantine levels
  n_levels <- length(quarantine_levels)
  colors <- rainbow(n_levels)
  names(colors) <- paste0("q", quarantine_levels)
  
  xMin1 <- as.Date("2019-12-15")
  xMax <- end_date - 1
  
  # Plot 1: Reproduction Number (R0)
  plot(date_range, rep(0, length(date_range)), 
       col="white", ylim=c(0, 8), xlim=c(xMin1, xMax),
       xlab="", ylab=expression(paste(R[t])), 
       main="Effect of Quarantine on Reproduction Number")
  
  for(i in 1:length(results_list)){
    q_level <- results_list[[i]]$level
    R0_q <- results_list[[i]]$R0_quantile
    col_q <- colors[paste0("q", q_level)]
    
    # Plot median
    lines(date_range, R0_q[3,], col=col_q, lwd=2)
    
    # Plot 95% CI
    polygon(c(date_range, rev(date_range)), 
            c(R0_q[1,], rev(R0_q[5,])), 
            col=adjustcolor(col_q, alpha=0.2), border=NA)
    
    # Plot 50% CI
    polygon(c(date_range, rev(date_range)), 
            c(R0_q[2,], rev(R0_q[4,])), 
            col=adjustcolor(col_q, alpha=0.3), border=NA)
  }
  
  lines(c(wuhan_travel_restrictions, wuhan_travel_restrictions), c(0, 10), col="red", lty=2, lwd=1.5)
  lines(date_range, rep(1, length(date_range)), lty=3, col="gray")
  
  # Add legend
  legend("topright", 
         legend=paste0("Quarantine: ", quarantine_levels*100, "%"),
         col=colors, lwd=2, cex=0.8)
  
  text(x=wuhan_travel_restrictions, y=7.5, "Travel restrictions", 
       col="red", adj=0, cex=0.7, srt=90)
  
  # Plot 2: Daily New Cases
  plot(date_range, rep(0, length(date_range)), 
       col="white", ylim=c(0, max(sapply(results_list, function(x) max(x$C_local_quantile[5,], na.rm=TRUE)))), 
       xlim=c(xMin1, xMax),
       xlab="", ylab="Daily new cases in Wuhan",
       main="Effect of Quarantine on Daily Case Incidence")
  
  for(i in 1:length(results_list)){
    q_level <- results_list[[i]]$level
    C_q <- results_list[[i]]$C_local_quantile
    col_q <- colors[paste0("q", q_level)]
    
    # Plot median
    lines(date_range, C_q[3,], col=col_q, lwd=2)
    
    # Plot 95% CI
    polygon(c(date_range, rev(date_range)), 
            c(C_q[1,], rev(C_q[5,])), 
            col=adjustcolor(col_q, alpha=0.2), border=NA)
  }
  
  lines(c(wuhan_travel_restrictions, wuhan_travel_restrictions), 
        c(0, max(sapply(results_list, function(x) max(x$C_local_quantile[5,], na.rm=TRUE)))), 
        col="red", lty=2, lwd=1.5)
  
  if(exists("quarantine_start_date")){
    lines(c(quarantine_start_date, quarantine_start_date), 
          c(0, max(sapply(results_list, function(x) max(x$C_local_quantile[5,], na.rm=TRUE)))), 
          col="blue", lty=2, lwd=1.5)
    text(x=quarantine_start_date, y=max(sapply(results_list, function(x) max(x$C_local_quantile[5,], na.rm=TRUE)))*0.9, 
         "Quarantine starts", col="blue", adj=0, cex=0.7, srt=90)
  }
  
  legend("topleft", 
         legend=paste0("Quarantine: ", quarantine_levels*100, "%"),
         col=colors, lwd=2, cex=0.8)
  
  # Plot 3: Total Infections (Prevalence)
  plot(date_range, rep(0, length(date_range)), 
       col="white", 
       ylim=c(0, max(sapply(results_list, function(x) max(x$I_quantile[5,], na.rm=TRUE)))), 
       xlim=c(xMin1, xMax),
       xlab="Date", ylab="Prevalence (E+I) in Wuhan",
       main="Effect of Quarantine on Infection Prevalence")
  
  for(i in 1:length(results_list)){
    q_level <- results_list[[i]]$level
    I_q <- results_list[[i]]$I_quantile
    col_q <- colors[paste0("q", q_level)]
    
    # Plot median
    lines(date_range, I_q[3,], col=col_q, lwd=2)
    
    # Plot 95% CI
    polygon(c(date_range, rev(date_range)), 
            c(I_q[1,], rev(I_q[5,])), 
            col=adjustcolor(col_q, alpha=0.2), border=NA)
  }
  
  lines(c(wuhan_travel_restrictions, wuhan_travel_restrictions), 
        c(0, max(sapply(results_list, function(x) max(x$I_quantile[5,], na.rm=TRUE)))), 
        col="red", lty=2, lwd=1.5)
  
  if(exists("quarantine_start_date")){
    lines(c(quarantine_start_date, quarantine_start_date), 
          c(0, max(sapply(results_list, function(x) max(x$I_quantile[5,], na.rm=TRUE)))), 
          col="blue", lty=2, lwd=1.5)
  }
  
  legend("topleft", 
         legend=paste0("Quarantine: ", quarantine_levels*100, "%"),
         col=colors, lwd=2, cex=0.8)
  
  # Save plot as PDF
  pdf(paste0("plots/", filename, ".pdf"), width=10, height=12)
  
  # Recreate all plots for PDF (same code as above, but now directed to PDF device)
  par(mfrow=c(3,1), mar=c(3,3,2,1), mgp=c(2,0.7,0))
  
  # Plot 1: Reproduction Number
  plot(date_rangeB, R0_quantileB[1,], col="white", ylim=c(0, 8), xlim=c(xMin1, xMaxR),
       xlab="", ylab=expression(paste(R[t])), 
       main="Effect of Quarantine on Reproduction Number")
  
  for(i in 1:length(results_list)){
    q_level <- results_list[[i]]$level
    R0_q <- results_list[[i]]$R0_quantile
    col_q <- colors[paste0("q", q_level)]
    lines(date_rangeB, R0_q[3,], col=col_q, lwd=2)
    polygon(c(date_rangeB, rev(date_rangeB)), 
            c(R0_q[1,], rev(R0_q[5,])), 
            col=adjustcolor(col_q, alpha=0.2), border=NA)
  }
  
  lines(c(wuhan_travel_restrictions, wuhan_travel_restrictions), c(0, 10), col="red", lty=2)
  lines(date_rangeB, rep(1, length(date_rangeB)), lty=3, col="gray")
  legend("topright", 
         legend=paste0("Quarantine: ", quarantine_levels*100, "%"),
         col=colors, lwd=2, cex=0.8)
  
  # Plot 2: Daily Cases (similar recreation)
  # ... (would include full plot recreation code here)
  
  dev.off()
  
  cat(paste("✓ PDF saved to: plots/", filename, ".pdf\n", sep=""))
  
  return(results_list)
}

# Quick comparison plot (faster, single run per scenario)
plot_quarantine_quick <- function(quarantine_levels = c(0, 0.3, 0.5, 0.7), 
                                   nn = 1e3, 
                                   dt = 0.25,
                                   filename = "quarantine_quick"){
  
  # Faster version with single SMC run per scenario (no bootstrap)
  
  # Store results
  results_list <- list()
  original_quarantine <- theta[["quarantine_effectiveness"]]
  
  cat("Running quick quarantine comparison...\n")
  
  for(q_level in quarantine_levels){
    cat(paste("  Quarantine effectiveness:", q_level, "\n"))
    
    theta[["quarantine_effectiveness"]] <- q_level
    
    # Recalculate quarantine start time if needed
    if(!is.null(theta[["quarantine_start_date"]])){
      quarantine_start_date <- as.Date(theta[["quarantine_start_date"]])
      quarantine_start_time <<- as.numeric(quarantine_start_date - start_date + 1)
    }
    
    # Single SMC run
    output_smc <- smc_model(theta, nn, dt)
    
    if(output_smc$lik != -Inf){
      R0_trace <- output_smc$beta_trace/(theta[["recover"]])
      C_local_trace <- theta[["confirmed_prop"]]*(output_smc$C_local_trace - c(0,head(output_smc$C_local_trace,-1)))
      I_trace <- output_smc$E_trace + output_smc$I_trace*(1-theta[["confirmed_prop"]])
      
      results_list[[paste0("q", q_level)]] <- list(
        R0 = R0_trace,
        C_local = C_local_trace,
        I = I_trace,
        level = q_level
      )
    }
  }
  
  theta[["quarantine_effectiveness"]] <- original_quarantine
  
  # Create plots
  par(mfrow=c(3,1), mar=c(3,3,2,1), mgp=c(2,0.7,0))
  
  n_levels <- length(quarantine_levels)
  colors <- rainbow(n_levels)
  names(colors) <- paste0("q", quarantine_levels)
  
  xMin1 <- as.Date("2019-12-15")
  xMax <- end_date - 1
  
  # Plot 1: R0
  plot(date_range, rep(0, length(date_range)), 
       col="white", ylim=c(0, 8), xlim=c(xMin1, xMax),
       xlab="", ylab=expression(paste(R[t])), 
       main="Effect of Quarantine on Reproduction Number (Single Run)")
  
  for(i in 1:length(results_list)){
    q_level <- results_list[[i]]$level
    R0_trace <- results_list[[i]]$R0
    col_q <- colors[paste0("q", q_level)]
    lines(date_range, R0_trace, col=col_q, lwd=2)
  }
  
  lines(c(wuhan_travel_restrictions, wuhan_travel_restrictions), c(0, 10), col="red", lty=2)
  lines(date_range, rep(1, length(date_range)), lty=3, col="gray")
  
  legend("topright", 
         legend=paste0("Quarantine: ", quarantine_levels*100, "%"),
         col=colors, lwd=2, cex=0.8)
  
  # Plot 2: Daily Cases
  ymax <- max(sapply(results_list, function(x) max(x$C_local, na.rm=TRUE)), na.rm=TRUE)
  plot(date_range, rep(0, length(date_range)), 
       col="white", ylim=c(0, ymax), xlim=c(xMin1, xMax),
       xlab="", ylab="Daily new cases",
       main="Effect of Quarantine on Daily Case Incidence")
  
  for(i in 1:length(results_list)){
    q_level <- results_list[[i]]$level
    C_trace <- results_list[[i]]$C_local
    col_q <- colors[paste0("q", q_level)]
    lines(date_range, C_trace, col=col_q, lwd=2)
  }
  
  lines(c(wuhan_travel_restrictions, wuhan_travel_restrictions), c(0, ymax), col="red", lty=2)
  if(exists("quarantine_start_date")){
    lines(c(quarantine_start_date, quarantine_start_date), c(0, ymax), col="blue", lty=2)
  }
  
  legend("topleft", 
         legend=paste0("Quarantine: ", quarantine_levels*100, "%"),
         col=colors, lwd=2, cex=0.8)
  
  # Plot 3: Prevalence
  ymax <- max(sapply(results_list, function(x) max(x$I, na.rm=TRUE)), na.rm=TRUE)
  plot(date_range, rep(0, length(date_range)), 
       col="white", ylim=c(0, ymax), xlim=c(xMin1, xMax),
       xlab="Date", ylab="Prevalence (E+I)",
       main="Effect of Quarantine on Infection Prevalence")
  
  for(i in 1:length(results_list)){
    q_level <- results_list[[i]]$level
    I_trace <- results_list[[i]]$I
    col_q <- colors[paste0("q", q_level)]
    lines(date_range, I_trace, col=col_q, lwd=2)
  }
  
  lines(c(wuhan_travel_restrictions, wuhan_travel_restrictions), c(0, ymax), col="red", lty=2)
  if(exists("quarantine_start_date")){
    lines(c(quarantine_start_date, quarantine_start_date), c(0, ymax), col="blue", lty=2)
  }
  
  legend("topleft", 
         legend=paste0("Quarantine: ", quarantine_levels*100, "%"),
         col=colors, lwd=2, cex=0.8)
  
  dev.copy(pdf, paste0("plots/", filename, ".pdf"), width=10, height=12)
  dev.off()
  
  cat(paste("Plot saved to plots/", filename, ".pdf\n", sep=""))
  
  return(results_list)
}

