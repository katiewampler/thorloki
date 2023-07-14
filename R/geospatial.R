#' Clean MTBS wildfire files
#'
#' MTBS downloads come in multiple zipped folders and the file names are long and meaningless. This
#' will pull all files in a single folder and rename files to be more user friendly.
#'
#'@importFrom sf st_read
#'@importFrom stringr str_detect str_replace_all str_split regex

#' @param wd the location with the raw MTBS files
#' @param save_loc the location to save all the files
#' @export
#' @examples
#' \dontrun{
#' wd <- "Z:/3_GIS-GPS/Wildfire/MTBS/McKenzie/"
#' save_loc <- "Z:/3_GIS-GPS/Wildfire/MTBS/McKenzie"
#' clean_mtbs(wd, save_loc)}

clean_mtbs <- function(wd, save_loc){
  stopifnot(file.exists(wd), file.exists(save_loc))

  #get zipped files of fires
  files <- list.files(wd, recursive = T)
  files <- files[stringr::str_detect(files, ".zip")]

  lapply(1:length(files), function(f){

    x <- files[f]

    #get place to put zip files
    exdir <- paste(wd, unlist(stringr::str_split(x, "/"))[1], sep="/")

    #upzip file
    unzip(paste(wd, x, sep="/"), exdir=exdir)

    #find shape file to extract fire name
    zip_files <- list.files(exdir)
    shp_file <- zip_files[stringr::str_detect(zip_files, "burn_bndy.shp")]

    #read file to get fire name
    name <- sf::st_read(paste(exdir, shp_file, sep="/"))
    filename <- name$Incid_Name

    #format name for file
    filename <- tolower(filename)
    filename <- stringr::str_replace_all(filename, " ", "_")
    filename <- gsub("[][()]", "", filename)

    #rename files
    #get dates of data
    dates <- str_split2(shp_file, "_")[2:3]

    #get header, remove files from other fires
    header <- str_split2(shp_file, "_")[1]
    zip_files <- zip_files[str_detect(zip_files, stringr::regex(header, ignore_case=T))]

    #if both, remove completely
    new_files <- gsub(paste("_", paste(dates, collapse ="_"), sep=""), "", zip_files)

    #if just one use pre or post
    new_files <- gsub(dates[1], "pre", new_files)
    new_files <- gsub(dates[2], "post", new_files)

    #replace long name with fire name
    new_files <- gsub(header, filename, new_files)

    #replace metadata name
    new_files <- gsub(toupper(header), filename, new_files)

    #rename files
    file.rename(paste(exdir, zip_files, sep="/"), paste(save_loc, new_files, sep="/"))
  })
}

#' Summarise raster dataset after clipping to a basin shape
#'
#' @importFrom graphics hist
#' @importFrom stats quantile
#' @importFrom utils unzip
#' @importFrom raster compareCRS mask crop freq cellStats
#' @importFrom sf st_transform
#' @importFrom terra crs
#' @import ggplot2
#'
#' @param raster object of class raster
#' @param basin object of class sf, a shapefile used to clip the raster optional
#' @param type either 'categorical' or 'numeric'
#' @param return either 'p' for plot or 't' for table (see details for more info)
#' @param data_type either 'character' or 'wildfire' (see details for more info)
#' @return if 'type' is numeric will provide mean, median, min, max, if categorical will return table with percent as a decimal for each category
#'
#' @details If data_type is 'character' it will return a table with either percentage of total in the basin
#' or if it's a numeric it will give quantiles to describe the range of values.
#' If the data_type is 'wildfire' it will return the percent burned at each severity
#' with 1 being unburned to low severity and 4 being high severity. If return is 'p', it will
#' provide a bar plot showing the relative counts of each group, 't' returns a table with
#' that data.
#'
summary_raster <- function(raster, basin=NULL, type="numeric", return ="t", data_type="character"){
  stopifnot(type %in% c("numeric", "categorical")| return %in% c("t", "p")| data_type %in% c("character", "wildfire")|
              class(raster) == "RasterLayer" | class(basin)[1] == c("sf"))
 if(is.null(basin) == F){
   if(raster::compareCRS(raster, basin) == F){
     basin_pj <- sf::st_transform(basin, terra::crs(raster))
   }else{basin_pj <- basin}

   #clip to basin file
   raster_crop <- raster::mask(raster, basin_pj)
   raster_crop <- raster::crop(raster_crop, basin_pj)
 }else{
   raster_crop <- raster
 }

  if(data_type == "character"){
    if(type=="categorical"){
      raster_sum <- as.data.frame(raster::freq(raster_crop))
      raster_sum <- raster_sum[order(raster_sum$count, decreasing=T),]
      raster_sum <- raster_sum[(is.na(raster_sum$value) == F),]
      raster_sum$count <- raster_sum$count / (sum(raster_sum$count))
      raster_plot <- raster_sum[1:10,]
      plot <- ggplot(raster_plot, aes(x=as.factor(value), y=count)) + geom_bar(stat="identity")
    }else{
      raster_sum <- quantile(raster_crop, probs=c(0.005, 0.01, 0.05, 0.25,0.5,0.75,0.95,0.99, 0.995))
      raster_sum <- c(raster::cellStats(raster_crop, min), as.numeric(raster_sum), cellStats(raster_crop, max), cellStats(raster_crop, mean))
      raster_sum <- data.frame(stat=c("min", "0.5%", "1%", "5%", "25%", "median", "75%","95%","99%","99.5%", "max", "mean"), raster_sum)
      plot <- hist(raster_crop,breaks=100)
    }
    if(return == "t"){
      raster_sum
    }else{
      plot
    }
  }else{
    raster_sum <- as.data.frame(freq(raster_crop))
    raster_sum <- subset(raster_sum, raster_sum$value %in% c(1:4))
    raster_sum$percentage <- round(raster_sum$count / sum(raster_sum$count),2)

    raster_sum
  }
}

#' Convert the projection geospatial data
#'
#' Converts a raster or shapefile into the same coordinate reference system as
#' another raster or shapefile.Helpful when trying to use multiple datasets.
#'
#' @importFrom terra res
#' @importFrom sf st_crs st_transform
#' @importFrom raster projectRaster crs
#' @param input the shapefile or raster you want to convert to a new coordinate reference system
#' @param goal the shapefile or raster with the coordinate reference system desired
#' @param type either 'numeric' or 'categorical' depending on if your raster is discrete or continuous values
#' @param res the resolution of the projected raster, if not specified with default to 30m
#' @return the input object with the coordinate reference system of the goal object
#' @export
#'
#' @examples
#' \dontrun{
#' convert_crs(basin, DEM)
#' }
convert_crs <-  function(input, goal, type="numeric", res=NULL){
  stopifnot(class(input)[1] %in% c("sf", "RasterLayer"),
            class(goal)[1] %in% c("sf", "RasterLayer"),
            type %in% c("numeric", "categorical"))

  #check first to see if a conversion is needed
  if(raster::compareCRS(goal, input) == T){
    message("CRS already matches that of your goal object")
    return(input)
  }else{
    if(class(input)[1] == "RasterLayer"){ #rasters need to be transformed differently
      #figure out the resolution of the projected raster
      if(is.null(res)){
        if(class(goal)[1] == "RasterLayer"){
          res <- terra::res(goal)[1]
        }else{
          #if a shapefile get the units of that to know how to set 30m
          unit <- sf::st_crs(goal, parameters = TRUE)$units_gdal
          res <- ifelse(unit == "degree", 0.0003280119,30)
        }
      }
      method <- ifelse(type == "numeric", "bilinear", "ngb") #set summary type
      #project raster
      input_prj <- raster::projectRaster(input, crs=raster::crs(goal), method=method, res=res)
    }else{
      #project shapefile
      input_prj <- sf::st_transform(input, raster::crs(goal))
    }
    return(input_prj)
  }}


#' Clip raster to shapefile
#'
#' Function will ensure the raster and shapefile are in the same coordinate reference
#' system, then will clip to the shapefile boundry.
#'
#' @importFrom raster rasterToPoints
#' @importFrom sf st_buffer
#' @param raster the raster you want to clip
#' @param sf the shapefile you want to use to clip the raster
#' @param type either 'numeric' or 'categorical' depending on if your raster is discrete or continuous values
#' @param res the resolution of the projected raster, if not specified with default to 30m
#' @param return either 'df' or 'raster' to specify the form of the returned raster
#'
#' @return if 'return' is df it will return the raster as a dataframe suitable for
#' plotting in ggplot2. If 'return' is raster it will return the raster as a
#' rasterLayer object
#' @export
#'
clean_raster <- function(raster, sf, type="numeric",
                         res=NULL, return="df"){
  stopifnot(class(raster) == "RasterLayer", class(sf)[1] == c("sf"),
              type %in% c("numeric", "categorical"),
              return %in% c("df", "raster"))
  #get units of raster to set buffer
  unit <- sf::st_crs(sf, parameters = TRUE)$units_gdal
  buffer <- ifelse(unit == "degree", 0.1,5000)
  res <- ifelse(unit == "degree", 0.0003280119,30)

  method <- ifelse(type == "numeric", "bilinear", "ngb") #set summary type

  #ensure it's in the right projection, if yes just clips to basin and formats as df
  if(compareCRS(raster, sf) == T){
    raster_crop <- raster::crop(raster, sf)
    raster_crop <- raster::mask(raster_crop, sf)
    raster_df <- as.data.frame(raster::rasterToPoints(raster_crop))
    colnames(raster_df) <- c("x", "y", "val")
  }else{
    #if not it will clip to smaller area before projecting to save time
    sf_prj <- convert_crs(sf::st_buffer(sf, dist=buffer), raster)
    raster_crop <- raster::crop(raster, sf_prj)
    raster_crop <- raster::mask(raster_crop, sf_prj)

    #reproject raster
    raster_prj <- raster::projectRaster(raster_crop, crs=crs(sf), method=method, res=res)

    #clip to actual basin area
    raster_crop <- raster::crop(raster_prj, sf)
    raster_crop <- raster::mask(raster_crop, sf)
    #format as df
    raster_df <- as.data.frame(raster::rasterToPoints(raster_crop))
    colnames(raster_df) <- c("x", "y", "val")
  }
  if(return == "df"){
    raster_df
  } else{raster_crop}
}
