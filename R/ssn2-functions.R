#new/updated SSN functions to run with SSN2 package

#' Run checks on SSN models
#'
#' @param model the fitted glmssn model
#' @param response the variable of interest
#'
#' @import SSN2
#' @import ggplot2
#' @importFrom caret R2
#' @importFrom patchwork wrap_plots
#'
#' @return a list of checks (1) is the spatial residuals with just the fixed effects
#' (2) is the a histogram of the standardized residuals
#' (3) is a table of potential outlier values based on residuals and cook's distance
#' (4) is the cross validation plot, potential outliers are in red and numbered with rid
#' (5) is plots of the explanatory variables vs the standarized residuals (currently only works with continous variables)
#' (6) is the R-squared value for the cross validation plot
#' (7) is the unexplained variance (the nugget) after accounting for the spatial autocorrelation and covariates
#'
check_model2 <- function(model, response){
  '_resid.stand_' <- exp_var <- residuals <- observe <- predict <- pid <- NULL
  #create list to hold plots
  plots <- list()
  dataset <- model$ssn.object

  #get residuals
  resid <- SSN2::residuals.glmssn(model)
  obs_resid <- getSSNdata.frame(resid)

  #see residuals
  plots[[1]] <-  plot_ssn_resid(model, obs_resid)
  plots[[2]] <- ggplot2::ggplot(data=obs_resid, aes(x=`_resid.stand_`)) + ggplot2::geom_histogram(bins=15)

  #get residuals with explanatory variables
  model_sum <- str_split2(as.character(model$args$formula)[3], pattern="[ + ]")
  model_sum <- subset(model_sum, model_sum != "")
  resid_plots <- list()
  for(x in 1:length(model_sum)){
    var <- model_sum[x]
    index <- which(colnames(obs_resid) == var)
    resid_data <- data.frame(residuals=obs_resid$`_resid.stand_`, exp_var = obs_resid[,index])
    plot <- ggplot2::ggplot(data=resid_data, aes(x=exp_var, y=residuals)) + ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept=0, color="red", linetype="dashed") +
      ggplot2::labs(x=var)
    resid_plots[[x]] <- plot
  }
  rplots[[5]] <- patchwork::wrap_plots(resid_plots)


  #check for potential outliers
  #get leverage outliers
  model_sum <- model$estimates
  model_sum2 <- model$args
  p <- nrow(model_sum$betahat) + length(model_sum2$CorModels) + 1
  pot_outlier <- subset(obs_resid, abs(obs_resid$`_resid.student_`) >2 | obs_resid$`_CooksD_` > 0.5)
  col_num <- which(colnames(pot_outlier) %in% c("_resid.student_", "_CooksD_", "pid", "SITE","_leverage_",response,
                                                "_resid.standard_", "_resid_", "_fit_"))
  pot_outlier <- pot_outlier[,col_num]
  pot_outlier$cooks_flag <- F
  pot_outlier$cooks_flag[(pot_outlier$`_CooksD_` > mean(obs_resid$`_CooksD_`, na.rm=T)*3) | (pot_outlier$`_CooksD_` > 1)] <- T
  pot_outlier$resid_flag <- F
  pot_outlier$resid_flag[pot_outlier$`_resid.student_` > 3] <- T

  plots[[3]] <- pot_outlier

  #see the model fit
  #get observations
  obs <- SSN2::getSSNdata.frame(dataset)
  cv.out <- SSN2::CrossValidationSSN(model)
  cols <- c(which(colnames(obs) == "pid"), which(colnames(obs) == response))
  cv.out <- merge(cv.out, obs[,cols], by="pid")
  colnames(cv.out)<- c("pid","predict", "se","observe")

  outliers <- subset(cv.out, cv.out$pid %in% pot_outlier$pid)
  outliers <- merge(outliers, pot_outlier[,1:2], by="pid")
  plots[[4]] <- ggplot2::ggplot(cv.out, aes(x=observe, y=predict)) +
    ggplot2::geom_point(alpha=0.5) +
    ggplot2::geom_abline (slope=1, linetype = "dashed", color="Red") +
    ggplot2::geom_point(data=outliers, aes(x=observe, y=predict), color="red3") +
    ggplot2::geom_text(data=outliers, aes(x=observe, y=predict, label=SITE),
                       vjust=1.1,hjust="inward", size=3) +
    labs(x=expression(bold("Observed DOC ("*mg*L^{-1}*")")),
         y=expression(bold("Predicted DOC ("*mg*L^{-1}*")"))) +
    theme_pub() + theme(axis.title = element_text(face="bold"))

  #get R2
  R2 <- caret::R2(cv.out$predict, cv.out$observe)
  plots[[6]] <- R2

  #get nugget
  var_comp <- SSN2::varcomp(model)
  plots[[7]] <- var_comp[which(var_comp$VarComp=="Nugget"),2]

  return(plots)
}

#' plot SSN data (updated for SSN2)
#'
#' takes either a glmssn fit or a glmssn.predict as an input and returns a list of the
#' items needed to make a nice clean map of the data. If you use multiple datasets it
#' will set them on a single scale for comparison accross the maps.
#'
#' @import SSN2
#' @import dplyr
#' @importFrom sf read_sf
#' @importFrom raster extent
#' @importFrom BAMMtools getJenksBreaks
#' @importFrom pals kovesi.rainbow cubicl ocean.haline parula
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#'
#' @param ssn_obj a list or individual ssn_glm or ssn_lm object
#' @param sites a vector of sites names, should match the ssn.ssn file names
#' @param type a vector of the type of data being plotted, should be either "obs" or "pred"
#' @param response the variable to plot as a character
#' @param size the size to plot to the variable, a vector the same length as models, either a number or "se
#' @param scale Either "continuous" or "binned"
#' @param nbins the number of bins to use in plotting
#' @param palette a character with the color palette to use. Currently 'parula', 'ocean.haline', 'cubicl', and 'kovesi.rainbow' from the 'pals' package are supported (and custom 'katie_pal')
#' @param prec an integer for the number of significant figures used for binning the intensity
#' @param bin_method either "jenks" or "quantile"
#'
#' @return a list of the items needed to make a nice clean map of the data
#'
plot_ssn2 <- function(ssn_obj, sites=c("sites", "preds"),
                     type=c("obs", "pred"),
                     response =c("DOC_mgL", "DOC_mgL"),
                     size=c(2, "se"),
                     scale="binned", nbins=8, palette="cubicl",
                     prec=2, bin_method="jenks"){
  stopifnot(length(ssn_obj)==length(sites),
            length(ssn_obj)==length(type),
            length(ssn_obj)==length(size),
            scale %in% c("continuous", "binned"),
            class(nbins) == "numeric",
            bin_method %in% c("jenks", "quantile"))

  val <- NULL
  #functions to make breaks and labels
  .place_val <- function(x) {
    stopifnot(class(x)=="numeric")
    place <- 0
    x <- abs(x)
    if(x == 1){
      place <- 1
    }
    else if(abs(x) > 1){
      while(x > 1){
        x <- x/10
        place <- place + 1
      }
    }else{
      while(x < 1){
        x <- x*10
        place <- place - 1
      }}
    place

  } #finds place value of max to determine scale
  .ceiling_dec <- function(x, level=-.place_val(x)) round(x + 5.0001*10^(-level-1), level) #ceiling function with decimal
  .floor_dec <- function(x, level=-.place_val(x)) round(x - 5.0001*10^(-level-1), level) #ceiling function with decimal

  #reformat data
  clean_data <- list()
  #get min and max of both groups for plotting
  min <- vector()
  max <- vector()
  vals <- vector()
  for(x in 1:length(type)){
    if(type[x] == "obs"){
      df <- SSN2::ssn_get_data(ssn_obj[[x]])
      plot_var_index <- which(colnames(df) == response[x])
      data <- df[,c(which(colnames(df)=="pid"),plot_var_index)]
      data <- stats::na.omit(data)
      data$size <- size[x]
      data <- st_drop_geometry(data)
      colnames(data) <- c("pid", "val", "size")
      clean_data[[x]] <- data
      min <- c(min, min(data$val))
      max <- c(max, max(data$val))
      vals <- c(vals, data$val)
    }else if(type[x]== "pred"){
      pred_points <- SSN2::ssn_get_data(ssn_obj[[x]], name=sites[x])
      pred_points <- st_drop_geometry(pred_points)
      data <- data.frame(pid=pred_points[,which(colnames(pred_points)=="pid")],
                         val=pred_points$preds, size=pred_points$se)
      clean_data[[x]] <- data
      min <- c(min, min(data$val))
      max <- c(max, max(data$val))
      vals <- c(vals, data$val)
    }else{stop("please choose pred or obs")}
  }

  min <- min(min, na.rm=T)
  max <- max(max, na.rm=T)
  vals <- stats::na.omit(vals) #to prevent errors

  #get streams
  if(length(which(type=="obs"))!= 0){
    loc <- ssn_obj[[which(type=="obs")[1]]]$path
  }else{
    loc <- ssn_obj[[which(type=="pred")[1]]][["ssn.object"]]$path
  }
  streams <- sf::read_sf(paste(loc, "edges.shp", sep="/"))

  if(scale=="binned"){
    #get data bins
    #get bins
    if(bin_method == "jenks"){
      breaks <- BAMMtools::getJenksBreaks(vals, nbins+1)
    }else if(bin_method == "quantile"){
      breaks <- unique(quantile(vals, prob = seq(0, 1, 1/(nbins))))
    }else{
      stop("please choose either jenks or quantile for binning method")
    }

    #add zero to breaks
    #if endpoint is less than zero, add that, otherwise zero
    start <- min(.floor_dec(min(vals), prec), 0)
    #round make cutoff nicer
    breaks <- .ceiling_dec(breaks, prec)
    breaks <- c(start, breaks[-1])

    #get color palette
    if(palette == "parula"){
      colors <- pals::parula(n=nbins)
    }else if(palette == "ocean.haline"){
      colors <- pals::ocean.haline(n=nbins)
    }else if(palette == "cubicl"){
      colors <- pals::cubicl(n=nbins)
    }else if(palette == "kovesi.rainbow"){
      colors <- pals::kovesi.rainbow(n=nbins)
    }else if(palette == "katie_pal"){
      katie <- grDevices::colorRampPalette(c("#d9ead3","#274e13"))
      colors <- katie(nbins)
    }else{
      stop("Please choose 'parula', 'ocean.haline', 'cubehelix', or 'kovesi.rainbow'")
    }

    #create labels
    labs_df <- data.frame(start=breaks)
    labs_df$end <- c(labs_df$start[2:nrow(labs_df)], NA)
    labs_df <- labs_df[-nrow(labs_df),]
    if(max(vals) < 1){
      labs_df$start <- sprintf(paste("%.", (max(nchar(labs_df$start))-2), "f", sep=""), labs_df$start)
      labs_df$end <- sprintf(paste("%.", (max(nchar(labs_df$end), na.rm=T)-2), "f", sep=""), labs_df$end)
    }
    labs_df$label <- paste(as.character(labs_df$start), as.character(labs_df$end), sep="-")
    labs <- labs_df$label
    labs[1] <- paste("<", labs_df$end[1], sep="")
    breaks <- c(as.numeric(labs_df$start),labs_df$end[nrow(labs_df)])

    #ensure the colors go to the right break
    labs_df$break_group <- paste("(",labs_df$start, ",",labs_df$end, "]", sep="")
    names(colors) <- labs_df$break_group
  }

  #get plotting info
  plots <- list()
  for(i in 1:length(sites)){
    #get data
    df <- clean_data[[i]]

    site_data <- sf::read_sf(paste(loc, "/", sites[i], ".shp", sep=""))

    #merge sites with data to plot
    plot_data <- merge(site_data, df, by="pid")

    #crop to area of interest
    area<- raster::extent(plot_data)
    x_lim <- c(area[1]-(area[2]-area[1])*.1, area[2]+(area[2]-area[1])*.1)
    y_lim <- c(area[3]-(area[4]-area[3])*.1, area[4]+(area[4]-area[3])*.1)

    #format SE sizes
    if(length(unique(df$size))> 1){
      #get bins
      if(bin_method == "jenks"){
        size_breaks <- BAMMtools::getJenksBreaks(df$size, 6)
      }else if(bin_method == "quantile"){
        size_breaks <- unique(quantile(df$size, prob = seq(0, 1, 1/(6))))
      }else{
        stop("please choose either jenks or quantile for binning method")
      }

      #round make cutoff nicer
      size_breaks <- .floor_dec(size_breaks, prec)

      #create labels
      labs_df <- data.frame(start=size_breaks)
      labs_df$end <- c(labs_df$start[2:nrow(labs_df)], NA)
      labs_df <- labs_df[-nrow(labs_df),]
      if(max(plot_data$val) < 1){
        labs_df$start <- sprintf(paste("%.", (max(nchar(labs_df$start))-2), "f", sep=""), labs_df$start)
        labs_df$end <- sprintf(paste("%.", (max(nchar(labs_df$end), na.rm=T)-2), "f", sep=""), labs_df$end)
      }
      labs_df$label <- paste(as.character(labs_df$start), as.character(labs_df$end), sep="-")
      size_labs <- labs_df$label
      size_breaks <- as.numeric(labs_df$start)
    }

    if(scale == "binned"){
      #bin data
      plot_data <- plot_data %>% mutate(bins = cut(val, breaks = c(breaks)))

      #fix labels
      #cut_labs <- labs[which(names(colors) %in% plot_data$bins)]

      #binned
      if(length(unique(df$size))==1){
        plot_spec <- list(streams=streams, data=plot_data, colors=colors, size=as.numeric(df$size[1]),color_labs=labs,  size_breaks=NA, size_labs=NA,x_extent=x_lim, y_extent=y_lim)

      }else{
        plot_spec <- list(streams=streams, data=plot_data, colors=colors, size=NA, color_labs=labs, size_breaks=size_breaks, size_labs=size_labs, x_extent=x_lim, y_extent=y_lim)
      }

    }else{
      #get palette
      if(palette == "parula"){
        colors <- pals::parula(100)
      }else if(palette == "ocean.haline"){
        colors <- pals::ocean.haline(100)
      }else if(palette == "cubicl"){
        colors <- pals::cubicl(100)
      }else if(palette == "kovesi.rainbow"){
        colors <- pals::kovesi.rainbow(100)
      }else if(palette == "katie_pal"){
        katie <- grDevices::colorRampPalette(c("#d9ead3","#274e13"))
        colors <- katie(100)
      }else{
        stop("Please choose 'parula', 'ocean.haline', 'cubehelix', or 'kovesi.rainbow'")
      }

      #continuous
      if(length(unique(df$size))==1){
        plot_spec <- list(streams=streams, data=plot_data, color=colors, color_min=min, color_max=max, size_breaks=NA, size_labs=NA,x_extent=x_lim, y_extent=y_lim)

      }else{
        plot_spec <- list(streams=streams, data=plot_data, colors=colors, color_min=min, color_max=max, size_breaks=size_breaks, size_labs=size_labs, x_extent=x_lim, y_extent=y_lim)

      }

    }

    plots[[i]] <- plot_spec
  }

  return(plots)

}

#' Plot glmssn residuals in ggplot2
#'
#' Given a glmssn object, it will extract the residuals and plot as an
#' ggplot
#'
#' @import SSN2
#' @importFrom sf read_sf
#' @importFrom ggspatial annotation_scale
#' @importFrom viridis scale_color_viridis
#'
#' @param model the glmssn object
#' @param sites the shapefile name for the sites in the .ssn object
#' @param size the column used to determine the size of the points
#'
#' @return a ggplot object of residuals on the stream network
#'
plot_ssn_resid2 <- function(model, sites="sites.shp", size="_resid_"){
  color_col <- NULL

  dataset <- model$ssn.object
  loc <- dataset@path
  obs <- SSN2::getSSNdata.frame(dataset)
  streams <- sf::read_sf(paste(loc, "edges.shp", sep="/"))
  sites <- sf::read_sf(paste(loc, "sites.shp", sep="/"))

  #get residuals
  resid <- SSN2::residuals.glmssn(model)
  plot_df <- SSN2::getSSNdata.frame(resid)

  #merge sites with data to plot
  plot_data <- merge(sites, plot_df, by="pid")
  #remove NA's
  plot_data <- plot_data[!is.na(plot_data$`_fit_`),]

  #rename things to plot
  colnames(plot_data)[which(colnames(plot_data) == "_resid_")] <- "color_col"
  colnames(plot_data)[which(colnames(plot_data) == size)] <- "size_col"

  ggplot2::ggplot() + ggplot2::geom_sf(data = streams, color="gray50")+
    ggplot2::theme_bw() + ggplot2::theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    ggspatial::annotation_scale() + ggplot2::geom_sf(data=plot_data, aes(color=color_col)) +
    viridis::scale_color_viridis(option = "C")+ ggplot2::labs(color="Residiuals")
}



#' Test multiple autocorrelation variance structures (fixed with SSN2)
#'
#' @param formula the best fitting varaibles from glmssn model
#' @param dataset the ssn object being modeled
#' @param corModels all the variance structures you want tested (tail up, tail down, euclidian)
#' @param addfunccol the additive function column, **currently broken, only uses afvArea
#' @param random a formula (~Random), with the random variables you want to add
#'
#' @import SSN2
#' @importFrom utils txtProgressBar setTxtProgressBar combn
#' @importFrom stringr str_replace_all
#'
#' @return a printout of the fits of the different models
#'
#' @examples
#' \dontrun{
#' best_var(DOC_mgL ~ PRECIP + BFI + AWC + AREA_km2 + sev, fall_ssn)}
#'
best_var_2 <- function(formula, dataset,
                     tailup_type = "exponential",
                     taildown_type = "exponential",
                     euclid_type = "exponential",
                     addfunccol = "afvArea",
                     random=NULL){
  corModels <- paste(c("tailup_","taildown_", "euclid_"), c(tailup_type,taildown_type,euclid_type), sep="")
  possible_mods <- do.call("c", lapply(seq_along(corModels), function(i) utils::combn(corModels, i, FUN = list)))

  #add progress bar
  pb <- utils::txtProgressBar(min = 0,max = length(possible_mods), style = 3,char = "=")

  models <- list()
  #run through all models
  for(x in 1:length(possible_mods)){
    cor <- possible_mods[[x]]
    cor <- data.frame(type=str_split2(cor, "_", piece=1), model=str_split2(cor, "_", piece=2))
    full_cor <- data.frame(type=c("tailup","taildown","euclid"))
    full_cor <- merge(full_cor, cor, by="type", all=T)
    full_cor$model[is.na(full_cor$model)==T] <- "none"
    full_cor <- full_cor %>% arrange(factor(type, levels = c("tailup","taildown","euclid")))

    model <- SSN2::ssn_lm(formula, dataset,
                          tailup_type = full_cor$model[1],
                          taildown_type = full_cor$model[2],
                          euclid_type = full_cor$model[3],
                          random= random,
                    additive = "afvArea")
    models[[x]] <- model
    utils::setTxtProgressBar(pb,x)
  }
  close(pb)
  #compare the fits
  fits <- as.data.frame(t(sapply(models, SSN2::glance)))
  fits2<- as.data.frame(t(sapply(models, SSN2::loocv)))
  fits <- cbind(fits,fits2)
  best_AIC <- models[[which(unlist(fits$AICc) == min(unlist(fits$AICc)))]]$formula
  best_RMSPE <- fits$Variance_Components[fits$RMSPE == min(fits$RMSPE)]

  fits$easy_names <-  sapply(1:nrow(fits), function(x){
    name <- str_split2(fits$Variance_Components[x], "[+]")
    name <- name[1:(length(name)-1)]
    name <- stringr::str_replace_all(name, " ", "")
    name <- paste(name, collapse=",")})

  print(fits[,c(3,5,8,11:13)])
  cat("The best variance models are:",
      paste("AIC: (", which(fits$Variance_Components == best_AIC), ") ", best_AIC, sep=""),
      paste("RMSPE: (",which(fits$Variance_Components == best_RMSPE), ") ", best_RMSPE, sep=""), sep="\n")
}

#' Test all potential types of variance structures
#'
#' @param formula the best fitting varaibles from glmssn model
#' @param dataset the ssn object being modeled
#' @param corstr the types of autocorrelation to include
#' @param cortype the types of autocorrelation models
#' @param addfunccol additive function column
#' @param random a character vector of random effects (the column name)

#' @import dplyr
#' @import SSN2
#' @importFrom tidyr expand
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stringr str_replace_all
#'
#' @return a df with the model fits

#' @examples
#' \dontrun{
#' test <- best_var_str(DOC_mgL ~ PRECIP + BFI + AWC + AREA_km2 + sev, fall_ssn,
#' corstr=c("tailup","taildown"))}
#'
best_var_str2 <- function(formula, dataset,
                         corstr=c("tailup", "taildown","Euclid"),
                         cortype=c("Exponential", "LinearSill","Spherical",
                                   "Mariah", "Epanech","Cauchy", "Gaussian"),
                         addfunccol = "afvArea",
                         random=NULL, parallel=T, ncores=NULL){
  stopifnot(corstr %in% c("tailup", "taildown","Euclid"),
            cortype %in% c("Exponential", "LinearSill","Spherical","Mariah", "Epanech","Cauchy" , "Gaussian"))
  cor1 <- cor2 <- cor3 <- NULL

  #remove models that have incompatible model/type
  invalid <- c(paste(c("Epanech", "LinearSill", "Mariah"),".Euclid", sep=""),
               paste(c("Cauchy", "Gaussian"), ".tailup", sep=""),
               paste(c("Cauchy", "Gaussian"), ".taildown", sep=""))

  #for each type of corstr, get all options
  if(length(corstr) == 1){
    corModels <- data.frame(cor1=paste(cortype, corstr, sep="."))
  }else if(length(corstr) == 2){
    corModels1 <- paste(cortype, corstr[1], sep=".")
    corModels2 <- paste(cortype, corstr[2], sep=".")
    corModels <- data.frame(cor1=corModels1, cor2=corModels2)
    corModels <- corModels %>% tidyr::expand(cor1, cor2)
  }else if(length(corstr) == 3){
    corModels1 <- paste(cortype, corstr[1], sep=".")
    corModels2 <- paste(cortype, corstr[2], sep=".")
    corModels3 <- paste(cortype, corstr[3], sep=".")
    corModels <- data.frame(cor1=corModels1, cor2=corModels2, cor3=corModels3)
    corModels <- corModels %>% tidyr::expand(cor1, cor2, cor3)
  }

  possible_mods <- do.call("c", lapply(seq_along(corModels), function(i) utils::combn(corModels, i, FUN = list)))




  corModels <- subset(corModels, !(corModels$cor1 %in% invalid))
  if(ncol(corModels) > 1){
    corModels <- subset(corModels, !(corModels$cor2 %in% invalid))}
  if(ncol(corModels) > 2){
    corModels <- subset(corModels, !(corModels$cor3 %in% invalid))}

  #make df into characters
  corModels <- data.frame(lapply(corModels, as.character), stringsAsFactors=FALSE)

  #try all the models
  #add progress bar
  pb <- utils::txtProgressBar(min = 0,max = nrow(corModels), style = 3,char = "=")

  #set up parallel processing
  if(is.null(ncore)){nCores <- parallelly::availableCores()}
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)

  models <- list()
  #run through all models
  if(parallel == T){
    foreach (i=1:nrow(corModels), .combine=rbind) %dopar% {
      model <- SSN2::glmssn(formula, dataset,
                            CorModels = c(unname(as.character(corModels[x,])),random),
                            addfunccol= addfunccol)
      models[[x]] <- model}
    parallel::stopCluster(cl) #close out, this is the last use

  }else{
    for(x in 1:nrow(corModels)){
      model <- SSN2::glmssn(formula, dataset,
                            CorModels = c(unname(as.character(corModels[x,])),random),
                            addfunccol= addfunccol)
      models[[x]] <- model
      utils::setTxtProgressBar(pb,x)
    }
    close(pb)
  }

  #compare the fits
  fits <- SSN2::InfoCritCompare(models)

  fits$easy_names <-  sapply(1:nrow(fits), function(x){
    name <- str_split2(fits$Variance_Components[x], "[+]")
    name <- name[1:(length(name)-1)]
    name <- stringr::str_replace_all(name, " ", "")
    name <- paste(name, collapse=",")})
  fits <- fits[,c(3,5,8,11:14)]
  best_AIC <- fits$Variance_Components[fits$AIC == min(fits$AIC)]
  best_RMSPE <- fits$Variance_Components[fits$RMSPE == min(fits$RMSPE)]

  cat("The best variance models are:",
      paste("AIC: (", which(fits$Variance_Components == best_AIC), ") ", best_AIC, sep=""),
      paste("RMSPE: (",which(fits$Variance_Components == best_RMSPE), ") ", best_RMSPE, sep=""), sep="\n")

  return(fits)

}

#' Test all potential types of variance structures
#'
#' @param formula the best fitting variables from glmssn model
#' @param dataset the ssn object being modeled
#' @param corstr the types of autocorrelation to include
#' @param cortype the types of autocorrelation models
#' @param addfunccol additive function column
#' @param random a formula (~Random), with the random variables you want to add

#' @import dplyr
#' @import SSN2
#' @importFrom tidyr expand
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stringr str_replace_all
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallelly availableCores
#'
#' @return a df with the model fits
#' @export
#' @examples
#' \dontrun{
#' test <- best_var_str(DOC_mgL ~ PRECIP + BFI + AWC + AREA_km2 + sev, fall_ssn,
#' corstr=c("tailup","taildown"))}
#'
best_var_str_ssn2 <- function(formula, dataset,
                              tailup_type = c("none", "linear", "spherical", "exponential","mariah","epa"),
                              taildown_type = c("none", "linear", "spherical", "exponential","mariah","epa"),
                              euclid_type = c("none", "spherical", "exponential", "gaussian","cosine", "cubic","pentaspherical",
                                              "wave", "jbessel", "gravity","rquad","magnetic"),
                              addfunccol = "afvArea",
                              random=NULL, parallel=T, nCores=NULL){
  stopifnot(tailup_type %in% c("linear", "spherical", "exponential","mariah","epa", "none"),
            taildown_type %in% c("linear", "spherical", "exponential","mariah","epa", "none"),
            euclid_type %in% c("spherical", "exponential", "gaussian","cosine", "cubic","pentaspherical",
                               "wave", "jbessel", "gravity","rquad","magnetic", "none"))

  #make a matrix with all options
  corModels <- data.frame(tailup=NA,taildown=NA, euclid=NA)
  for(x in tailup_type){
    for(y in taildown_type){
      for(z in euclid_type){
        run <- data.frame(tailup=x, taildown=y, euclid=z)
        corModels <- rbind(corModels, run)}}}
  corModels <- corModels[-1,]

  #try all the models
  #add progress bar
  pb <- utils::txtProgressBar(min = 0,max = nrow(corModels), style = 3,char = "=")

  if(parallel==T){
    #set up parallel processing
    if(is.null(nCores)){nCores <- parallelly::availableCores()}
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
  }

  models <- list()
  #run through all models
  if(parallel==T){
    foreach (x=1:nrow(corModels), .combine=rbind) %dopar%{
      model <- SSN2::ssn_lm(formula, dataset,
                            tailup_type = corModels[x,1],
                            taildown_type = corModels[x,2],
                            euclid_type = corModels[x,3],
                            random= random,
                            additive = "afvArea")
      models[[x]] <- model
    }
    parallel::stopCluster(cl) #close out, this is the last use
  }else{
    for(x in 1:nrow(corModels)){
      model <- SSN2::ssn_lm(formula, dataset,
                            tailup_type = corModels[x,1],
                            taildown_type = corModels[x,2],
                            euclid_type = corModels[x,3],
                            random= random,
                            additive = "afvArea")
      models[[x]] <- model
      utils::setTxtProgressBar(pb,x)
    }
    close(pb)
  }

  #compare the fits
  fits <- as.data.frame(t(sapply(models, SSN2::glance)))
  fits2<- as.data.frame(t(sapply(models, SSN2::loocv)))
  fits <- cbind(fits,fits2)
  best_AIC <- which(unlist(fits$AICc) == min(unlist(fits$AICc)))
  best_RMSPE <- which(unlist(fits$RMSPE) == min(unlist(fits$RMSPE)))
  best_R2 <- which(unlist(fits$pseudo.r.squared) == min(unlist(fits$pseudo.r.squared)))

  fits$formula <- paste(corModels$tailup, corModels$taildown, corModels$euclid, sep=" + ")

  cat("The best variance models are:",
      paste("AIC:", fits$formula[best_AIC], sep=""),
      paste("RMSPE:",fits$formula[best_RMSPE], sep=""),
      paste("R2:",fits$formula[best_R2] , sep=""),sep="\n")

  return(fits)

}

#' Plot and report variable importance based on standarized linear regression coefficients
#'
#' @param model the fitted glmssn object
#'
#' @import ggplot2
#' @importFrom stats reorder
#'
#' @return A list, with the first object being a table of the coefficients and the second
#' being a plot of the coefficients. Includes the standard error bars. If you did model
#' selection these are likely larger than is shown, but you can still show, just make it
#' very clear in the paper that these are only error from the final model.
#'
var_imp2 <- function(model){
  nice_name <- Estimate <- direction <- lower <- higher <- NULL

  nice_names <- data.frame(
    vars=c("MIN_TEM","MAX_TEM","DNBR","SOIL_OM","AWC","LOG_AREA","PRECIP",
           "ARIDITY","Sample_Time","ELEV","CLAY","TWI","PH","BFI","FOREST","SLOPE"),
    nice_name=c("Minimum Temperature", "Maximum Temperature","Burn Severity (dNBR)",
                "Soil Organic Matter","Available Water Capacity","log of Basin Area",
                "Annual Precipitation","Aridity Index",
                "Sample Time","Elevation","Soil Clay Percentage","Topographic Wetness Index",
                "Soil pH", "Baseflow Index","Percent Forested","Slope"))


  coeff <- summary(model)
  vals <- coeff[["fixed.effects.estimates"]]
  vals <- subset(vals, vals$FactorLevel != "(Intercept)")
  vals <- vals[order(abs(vals$Estimate), decreasing=T),]
  vals$order <- 1:nrow(vals)
  vals$percentage <- abs(vals$Estimate) / sum(abs(vals$Estimate)) * 100
  vals$direction <- "positive"
  vals$direction[vals$Estimate <0 ] <- "negative"
  vals$direction <- factor(vals$direction, levels=c("positive", "negative"), ordered=T)
  vals$lower <- NA
  vals$higher <- NA
  for(x in 1:nrow(vals)){
    ci1 <- vals$Estimate[x] + vals$std.err[x]
    ci2 <- vals$Estimate[x] - vals$std.err[x]
    vals$lower[x] <- min(c(ci1, ci2))
    vals$higher[x] <- max(c(ci1,ci2))
  }

  vals <- merge(vals, nice_names, by.x="FactorLevel", by.y="vars")
  plot <- ggplot2::ggplot(vals, aes(y=stats::reorder(nice_name,abs(Estimate)), x=Estimate,
                           fill=direction)) +
    ggplot2::geom_bar(stat="identity") +
    theme_pub() + ggplot2::labs(x= "Standardized Coefficient", y="Parameter",
                       fill="Coefficient \n Direction") +
    ggplot2::theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
          panel.grid.major.y = element_blank(), axis.text.y=element_text(size=14)) +
    ggplot2::scale_fill_manual(values=c("#ef8a62","#67a9cf")) + theme(legend.position = "none") +
    ggplot2::theme(        plot.margin = margin(5.5, 10, 5.5, 5.5, "pt")) +
    ggplot2::geom_errorbar(aes(xmin=lower,
                      xmax=higher), width=.2,
                  position=position_dodge(.9))
  return(list(vals, plot))
}


#' Calculate impact of burn severity
#'
#' uses the dNBR thresholds from the Holiday Farm fire to extract the impact of
#' burn severity on DOC (mg/L)
#'
#' @importFrom stats confint
#' @param model the fitted glmssn object
#'
#' @return a table showing the range and 95% CI of burn severity for each severity level
#'
burn_impact2 <- function(model){
  coeff <- summary(model)
  vals <- coeff[["fixed.effects.estimates"]]
  slope <- vals$Estimate[vals$FactorLevel == "DNBR"]
  lowerslope <- vals$Estimate[vals$FactorLevel == "DNBR"] + 2*vals$std.err[vals$FactorLevel == "DNBR"]
  upperslope <- vals$Estimate[vals$FactorLevel == "DNBR"] - 2*vals$std.err[vals$FactorLevel == "DNBR"]

  unburned <- c(lowerslope,slope,upperslope) *40
  low <- c(lowerslope,slope,upperslope)*320
  moderate <-c(lowerslope,slope,upperslope)*660
  high <- c(lowerslope,slope,upperslope)*772
  fire <- data.frame(severity_group = c("Unburned","Low","Moderate","High"),
                     min_change_L95 = c(0,unburned[1],low[1], moderate[1]),
                     min_change = c(0,unburned[2],low[2], moderate[2]),
                     min_change_U95 = c(0,unburned[3],low[3], moderate[3]),
                     max_change_L95 = c(unburned[1], low[1], moderate[1], high[1]),
                     max_change = c(unburned[2], low[2], moderate[2], high[2]),
                     max_change_U95 = c(unburned[3], low[3], moderate[3], high[3]))
  return(fire)

}

#' Perform double selection using a p-value limit to determine which variable to include
#'
#' @param focal a character of the column name of the variable you're interested in
#' @param response a character of the column name of the response variable
#' @param p_limit variable with p-value below this limit will be included in the final model
#' @param df the dataframe with the data to fit in it
#' @param covar a vector of column names of the variables to include as potential covariates
#'
#' @importFrom stats lm as.formula
#' @return a formula for the final model
#'
ds_select2 <- function(focal="DNBR", response="DOC_mgL", p_limit=0.1, df, covar){
  focal_col <- which(colnames(df) %in% focal) #get column(s) of focal
  covar_col <- which(colnames(df) %in% covar) #get column(s) of covariates
  response_col <- which(colnames(df) %in% response) #get column of response
  fit_df <- df[,c(response_col, focal_col, covar_col)]

  modeld0 <- stats::lm(stats::as.formula(paste(response, "~", paste(covar, collapse=" + "))), fit_df)
  selected0 <- as.data.frame(summary(modeld0)$coefficients)
  selected0 <- subset(selected0, selected0$`Pr(>|t|)` < p_limit)
  selected0 <- row.names(selected0)

  #step 2: fit model with sev ~ covariates
  modeld1 <- stats::lm(stats::as.formula(paste(focal, "~", paste(covar, collapse=" + "))), fit_df)
  selected1 <- as.data.frame(summary(modeld1)$coefficients)
  selected1 <- subset(selected1, selected1$`Pr(>|t|)` < p_limit)
  selected1 <- row.names(selected1)

  #step 3: final model includes significant variables for both
  selected=sort(union(selected0,selected1))
  selected <- selected[selected != "(Intercept)"]
  form <- stats::as.formula(paste(response, "~", paste(selected, collapse=" + "), "+", focal))
  form  }
