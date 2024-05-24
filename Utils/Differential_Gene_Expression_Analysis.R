#!/usr/bin/env Rscript

# Rscript to define functions necessary to run scripts
# @Author: Clara Savay (clara.savary@yahoo.fr, OmixAnalytics)
# Date: April 16th, 2024


# Elbow algorithm
elbow_pcs <- function(object, ndims = 50, graph = "pca"){
  
  # Perform Elbow plot for ranking principle components based on the percentage of variance explained by each one
  ElbowPlot(object, ndims = ndims)
  
  ## First metric:
  # Determine percent of variation associated with each PC
  pct <- object[[graph]]@stdev / sum(object[[graph]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  ## Second metric:
  # Identify the PC where the percent change in variation between consecutive PCs is less than 0.1%:
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  ## Create a dataframe with values
  plot_df <- data.frame(
    pct  = pct,
    cumu = cumu,
    rank = 1:length(pct)
  )
  # Elbow plot to visualize 
  p <- ggplot(
    plot_df,
    aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()
  print(p)
  
  message(glue::glue("Number of Principal Components (npcs) selected: ", pcs))
  
  return(pcs)
  
}



