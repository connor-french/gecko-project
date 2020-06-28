require(raster)
require(dplyr)
require(ggplot2)
require(reshape2)
require(forcats)

raster_correlation <- function(folder = getwd(), raster_suffix = ".tif", cutoff = 0.75, factor_vec = NULL) {
  
  #get all of the raster files for correlation analysis
  files <- list.files(folder, pattern = paste0(raster_suffix, "$"), full.names = TRUE)
  
  raster_stack <- raster::stack(files)
  
  #calculate the pairwise pearson correlation coefficient
  cors <- raster::layerStats(raster_stack, "pearson", na.rm = TRUE)
  
  cor_values <- cors$`pearson correlation coefficient`
  
  #reorder the correlation matrix for visualization
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  cor_values <- reorder_cormat(cor_values)
  
  #get the lower triangle of the correlation matrix. The top triangle is redundant
  cor_values[upper.tri(cor_values)] <- NA
  
  #convert to a tibble for plotting
  cor_lower <- reshape2::melt(cor_values) %>% 
    as_tibble() 
  
  #plot a heatmap of the correlation matrix
  cor_plot <- ggplot(data = cor_lower, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", 
                         mid = "white", 
                         high = "red", 
                         na.value = "transparent",
                         name = "Pearson\nCorrelation") + 
    labs(x = NULL, y = NULL) +
    theme_minimal()
  
  #return a tibble of the raster pairs with the highest correlation coefficient
  cor_lower_filtered <- na.omit(cor_lower) %>% 
    filter(abs(value) > cutoff, Var1 != Var2) %>% 
    mutate(Var1 = forcats::fct_relevel(Var1, factor_vec),
           Var2 = forcats::fct_relevel(Var2, factor_vec)) %>% 
    arrange(Var1, Var2)
  
  list(cor_heatmap = cor_plot, top_cor = cor_lower_filtered, all_cor = cor_lower)
}


