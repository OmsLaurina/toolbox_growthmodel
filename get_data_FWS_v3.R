##################################################
######             get_model_FWS.R          ######
##################################################
# get FWS data for biovolume calculation    ######
# select data in study region               ######
##################################################
# input : data FLR6 or FLR25                ######
# output : files.txt                        ######
##################################################
# cruise: Protevsmed-SWOT 2018              ######
# author: Lloyd Izard                       ######
# modif by Roxane Tzortzis 09/2021          ######
# modif by Laurina Oms 01/2023              ######
##################################################


rm(list=ls())
graphics.off()
cat('\f')

# library
library(data.table)
library(stringr)
library(rowr)

# Set the working place
setwd("~/Bureau/PROTEVSWOT/")
inputs = "~/Bureau/PROTEVSWOT/"
outputs = "~/Bureau/PROTEVSWOT/"


#Uncomment the chosen FLR

#FLR6
path_cyto <- "~/Bureau/PROTEVSWOT/Processed_FLR6/"
wanted_group="cryptophyte" #"red_pico_2" #"red_pico_euk" #"synechococcus"
index_vol_FLR = 24
car_start = 25
car_end = 41
sep_files = ','


# # #FLR25
# path_cyto <- "~/Bureau/PROTEVSWOT/Processed_FLR25/"
# wanted_group= ""  #"microphyto" #"red_nano_euk" #"red-nano-euk-sws" #"red-pico-euk-3" #"pico_HFLR" #
# index_vol_FLR = 21
# car_start = 30
# car_end = 46
# sep_files = ';'

#Create empty vector to save data
list_res <- NULL
info_res <- NULL
mat_res <- NULL

# Select the good data for water 1 and 2
tmp_path <- paste(path_cyto, "for_model/", sep = "")
cat(tmp_path)
cat("\n")

#####
# GET LISTMODES
#####

my_list_per_day <- list.files(path = tmp_path, pattern = "*.csv", all.files = F, full.names = F) # All files
my_list_per_day <- my_list_per_day[my_list_per_day %like% "listmode"] # Listmodes
my_list_per_day <- as.data.frame(my_list_per_day[my_list_per_day %like% wanted_group]) # Select data with wanted group
my_list_per_day <- data.frame(column_name = my_list_per_day[5:105, ])

cat(nrow(my_list_per_day))
cat("\n")

#####
# Check length listmode
#####

list_res <- rbind(list_res, my_list_per_day)

#####
# GET INFO
#####

my_info_per_day <- list.files(path = tmp_path, pattern = "*.txt", all.files = F, full.names = F)
my_info_per_day <- my_info_per_day[my_info_per_day %like% "info"]
my_info_per_day <- as.data.frame(my_info_per_day[my_info_per_day %like% wanted_group])
my_info_per_day <- data.frame(column_name = my_info_per_day[5:105, ])

#####
# Check length info
#####

info_res <- rbind(info_res, my_info_per_day)

#####
# GET VOLUME
#####

for (i in 1:nrow(my_info_per_day)){
  tmp_file <- as.character(my_info_per_day[i,1])
  file <- paste(tmp_path, '/', tmp_file, sep = "")
  check <- paste("\t", file, sep = " ")
  cat(check)
  cat("\n")
  tmp_data <- read.table(file, header = T, sep = ",", fill=T)
  mat_res <- rbind(mat_res, as.character(tmp_data[index_vol_FLR,]))
}

#####
# GET CLEAN VOLUME
#####

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

vol <- as.data.frame(numextract(mat_res[,1]))

#####
# GET DATE
#####

list_date <- as.data.frame(sapply(list_res, substring, car_start, car_end))
colnames(list_date) <- 'Date'
list_date$Date <- as.POSIXct(list_date$Date, format = "%Y-%m-%d %Hh%M")

# date_debut <- as.POSIXct("2018-05-11 02:00:00", format = "%Y-%m-%d %H:%M:%S")
# date_fin <- as.POSIXct("2018-05-13 08:30:00", format = "%Y-%m-%d %H:%M:%S")
# 
# indices_dates_conformes <- which(list_date$Date >= date_debut & list_date$Date <= date_fin)
# list_date <- subset(list_date, Date >= date_debut & Date <= date_fin)

## arrondi a la 1/2 heure pres

my_time <- strptime(list_date$Date, "%Y-%m-%d %H:%M:%S")

for (i in 1:length(my_time)){
  
  my_time$sec[i] <- 0
  if((my_time$min[i] >= 15) && (my_time$min[i] <= 45)){
    my_time$min[i] <- 30
  }
  if ((my_time$min[i] <= 15) || (my_time$min[i] >= 45)){
    if(my_time$min[i] <= 15){my_time$min[i] <- 0}
    if (my_time$min[i] >= 45){
      my_time$min[i] <- 0
      my_time$hour[i] <- my_time$hour[i] + 1
    }
  }
}

list_date$Date <- as.POSIXct(my_time, format = "%Y-%m-%d %H:%M:%S")

#####
# GET PAR
#####
file_url <- file.path(inputs, "PAR.txt")
tmp_par <- read.table(file_url, header = T, sep = '\t')
tmp_par$Date <- as.POSIXct(tmp_par$Date, format = "%Y-%m-%d %H:%M:%S")

tmp <- which(tmp_par$Date %in% list_date$Date)
tmp_par <- tmp_par[tmp,] 
row.names(tmp_par) <- NULL

#####
# GET SELECTION SET
#####

sel_set <- rep(wanted_group, nrow(tmp_par))

#### FINAL DATA FRAME = FWS

final_mat <- cbind(list_date$Date, sel_set, list_res, vol, tmp_par$PAR)
colnames(final_mat) <- c("Sampling.Date", "Selection.Set", "File", "Volume", "PAR2")

file_url_mat <- file.path(inputs, paste("PAR_volume_files",wanted_group,".txt"))
write.table(final_mat, file_url_mat, col.names = T, row.names = F, sep = ';', quote = F)


### FWS
final_FWS <- NULL

for (i in 1:nrow(final_mat)){
  tmp <- paste(tmp_path, final_mat$File[i], sep = '')
  cat(tmp)
  cat('\n')
  tmp_data <- read.table(tmp, sep = sep_files, dec = '.', header = T,fill=T) #virgule for small
  cat(nrow(tmp_data))
  cat('\n')
  final_FWS <- cbind.fill(final_FWS, tmp_data$FWS.Total, fill = NA)
}

file_url_final <- file.path(inputs, paste("FWS_ALL",wanted_group,".txt"))
write.table(final_FWS, file_url_final, col.names = T, row.names = F, sep = ';', quote = F)
