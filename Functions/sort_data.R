##%######################################################%##
#                                                          #
####             Function to organise data              ####
####          for predictions for forest paper          ####
#                                                          #
##%######################################################%##


# note that currently all data above is needed first
sort_data <- function(modout = NULL, moddata = NULL, scalers = NULL, from = NULL,
                      to = NULL, vals = vals, variable = NULL, fac = NULL, n = NULL, logval = FALSE){
  
  
  
  # 1. Create table of values to predicted
  
  # basic table of median values and reference factors
  pred_tab <- data.frame(landcovers.5k = median(moddata$landcovers.5k),
                         homogen = median(moddata$homogen),
                         fert.total_log = median(moddata$fert.total_log),
                         percNH = median(moddata$percNH),
                         Hansen_mindist_log =  median(moddata$Hansen_mindist_log),
                         Forest_biome = "Temperate Broadleaf & Mixed Forests",
                         Use_intensity = "Minimal use",
                         Predominant_land_use = "Primary vegetation",
                         Tropical = "Temperate",
                         Species_richness = 0,
                         logAbun = 0, 
                         RCAR_110km = 0)
  
  
  # organise factor levels
  # check levels of factor variables
  levels(pred_tab$Predominant_land_use) <- levels(modout$data$Predominant_land_use)
  levels(pred_tab$Use_intensity) <- levels(modout$data$Use_intensity) 
  levels(pred_tab$Forest_biome) <- levels(modout$data$Forest_biome) 
  
  if(!is.null(modout$data$Tropical)){
  levels(pred_tab$Tropical) <- levels(modout$data$Tropical) 
  }
  
  
  # 2. make changes to table to specify values to predict
  # number of values to predict
  
  
  # transform the values
  vals_trans <- sapply(vals, 
                       FUN = rescale, 
                       centre = scalers[scalers$variable == variable, 'centre'], 
                       scale = scalers[scalers$variable == variable, 'scale'],
                       logval = logval)
  
  
  # add reps of pred_tab and add the new variable vals
  pred_tab2 <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))
  
  # replace the variable of interest with the new ones
  pred_tab2[,variable] <- vals_trans
  
  if(is.null(fac)){
  return(pred_tab2)
    }
  
  #### if looking at interactions ####
  
  if(!is.null(fac) & !is.null(n)){
    
    # add more rows for the dif factor levels
    
    if(n >=2){
      pred_tab3 <- do.call("rbind", replicate(n, pred_tab2, simplify = FALSE))
      pred_tab3[, fac][(nrow(pred_tab3)/n + 1):(nrow(pred_tab3)/n + 1000)] <- levels(moddata[, fac])[2]
      
    }
    
    if(n >=3){
      
      pred_tab3[, fac][(nrow(pred_tab3)/n + 1001):(nrow(pred_tab3)/n + 2000)] <- levels(moddata[, fac])[3]
      
      
    }
    
  }
  
  return(pred_tab3)
}
