date <- "2022_11_01"
group <- "Group 1"
wd <- paste("Z:/1_Research/2_PNNL_EWEB/Photos", date, group, sep="/")
datasheet <- read.csv(paste("Z:/1_Research/2_PNNL_EWEB/Field Notes/", date, "/Master_Sample_List_", date, ".csv", sep=""))

#' Clean field photo names
#'
#' Uses the times sites were visited from the datasheet and the metadata from
#' the photos to rename confusing/undescriptive names from phones and tablets
#' to the site it was taken at with the data and time the photo was taken.
#'
#' @importFrom stringr str_detect str_sub str_replace_all str_split
#' @import dplyr
#' @importFrom exifr read_exif
#' @importFrom data.table as.ITime
#' @importFrom stats sd

#'
#' @param wd the file path where the photos are located written as a string
#' @param datasheet a dataframe of the datasheet
#' @param group string with the name of the group (see details)
#'
#' @description The group name is used to subset the datasheet if there are
#' multiple groups on the datasheet. This should be a string.
#' Datasheet should have the following columns: Group (if grouping is used),
#' Sample_Time_PST, Site_Num
#'
#' @returns The photos in the same place just renamed with the site and time.
#' @export
#'
#' @examples
#' \dontrun{
#' date <- "2022_11_01"
#' group <- "Group 1"
#' wd <- paste("Z:/1_Research/2_PNNL_EWEB/Photos", date, group, sep="/")
#' datasheet <- read.csv(paste("Z:/1_Research/2_PNNL_EWEB/Field Notes/",
#' date, "/Master_Sample_List_", date, ".csv", sep=""))
#' clean_photo_names(wd, datasheet, group="Group1")}

clean_photo_names <- function(wd, datasheet, group=NA){
  stopifnot(c("Sample_Time_PST", "Site_Num") %in% colnames(datasheet), file.exists(wd))
  files <- list.files(wd)
  files <- files[stringr::str_detect(files, pattern=".jpg|.jpeg")]
  avenza <- stringr::str_detect(files, pattern="_IMG") #won't pull date time from avenza photos

  #read datasheet
  if(is.na(group) == F){
    stopifnot("Group" %in% colnames(datasheet))
    df <- subset(datasheet, datasheet$Group == group)}
  df$Sample_Time_PST <- data.table::as.ITime(df$Sample_Time_PST)

  #get metadata
  meta <- exifr::read_exif(paste(wd,files, sep="/"), tags=c("FileName", "DateTimeOriginal"))

  #check on time zones
  min_photo <-  data.table::as.ITime(stringr::str_sub(min(meta$DateTimeOriginal, na.rm=T), 12,-1))
  min_meta <- data.table::as.ITime(min(df$Sample_Time_PST))

  max_photo <- data.table::as.ITime(stringr::str_sub(max(meta$DateTimeOriginal, na.rm=T), 12,-1))
  max_meta <- data.table::as.ITime(max(df$Sample_Time_PST))

  time_dif_min <- as.numeric(difftime(min_meta,min_photo, units="hours"))
  time_dif_max <- as.numeric(difftime(max_meta,max_photo, units="hours"))

  #check if PST was used
  if(time_dif_min >= -0.9 & time_dif_max >= -0.9){
    message("Photos are already adjusted for PST")
    PST <- T #true if yes, false if it was in PDT

  }else{
    message("Photos are not adjusted for PST")
    PST <- F #true if yes, false if it was in PDT
  }

  #fix avenza metadata (doesn't pull in date/time taken)
  if(sum(avenza) > 0){
    avenza_time <- stringr::str_sub(meta$FileName[avenza], start=-10, end=-5)
    avenza_time <- paste(stringr::str_sub(avenza_time, start=1, end=2),
                         ":", stringr::str_sub(avenza_time, start=3, end=4),
                         ":", stringr::str_sub(avenza_time, start=5, end=6), sep="")
    avenza_date <- stringr::str_sub(meta$FileName[avenza], start=-19, end=-12)
    avenza_date <- paste(stringr::str_sub(avenza_date, start=1, end=4), ":",
                         stringr::str_sub(avenza_date, start=5, end=6), ":",
                         stringr::str_sub(avenza_date, start=7, end=8), sep="")
    meta$DateTimeOriginal[avenza] <- paste(avenza_date, avenza_time)
  }
  #remove any remaining photos without date
  before <- nrow(meta)
  remove <- meta[is.na(meta$DateTimeOriginal) == T,]
  meta <- meta[is.na(meta$DateTimeOriginal) == F,]
  after <- nrow(meta)
  warning(paste(paste(before-after, "samples were removed because they didn't have times:"),
                paste(" ", remove$FileName, collapse="\n"), sep="\n"))

  meta$time <- unlist(stringr::str_split(meta$DateTimeOriginal, pattern=" "))[seq(from=2, to=nrow(meta)*2, by=2)]
  meta$DateTimeOriginal <- as.POSIXct(meta$DateTimeOriginal, format="%Y:%m:%d %H:%M:%S")
  meta$time <- data.table::as.ITime(meta$time)

  if(PST == F ){
    #convert from PDT to PST on phone if not manually corrected and it's DST
    meta$time_PST <- meta$time - data.table::as.ITime("01:00")
  }else{meta$time_PST <- meta$time}

  meta$site <- NA

  #relabel photos
  sapply(1:nrow(meta), function(x){
    photo <- meta$FileName[x]
    time <- meta$time_PST[x]
    time_name <- stringr::str_replace_all(time, ":", "_")
    #site <- df$Site_Num[max(which(df$Sample_Time_PST <= time))
    site <- df$Site_Num[which.min(abs(df$Sample_Time_PST - time))]

    meta$site[x] <- site
    newfile <- paste(site, date, time_name, sep="_")
    newfile <- paste(newfile, ".jpg", sep="")

    file.rename(from=paste(wd, photo, sep="/"), to=paste(wd, newfile, sep="/"))
  })

  #check for inconsistencies
  check <- meta %>% group_by(site) %>% summarise(sd_time = sd(time))
  pot_issues <- subset(check, check$sd_time >= 40)
  warning(paste("sites", paste(pot_issues$site, collapse=", "), "need to be checked"))

  missing <- df$Site_Num[!(df$Site_Num %in% meta$site)]

  warning(paste("missing sites", paste(missing, collapse=", ")))

}
