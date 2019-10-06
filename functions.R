library(raster)
library(snowfall)
library(Hmisc)
library(MODIS)

# existing tiles
# Return a data.frame containing identifiers of all the tiles existing on the server.
# The data.frame contains two columns, tilesH and tilesV, with tilesH the horizontal
# of the tiles and tilesV the vertical identifier.
# Each line represents a tile.
existing_tiles <- function(){
    # Accessing a csv file in "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MOD11A2/2005/001"
    # containing the list of HDF files for the first date in 2005.
    # The date used does not really matter here as this csv file is used to extract a list of all tiles.
    hdf_files <- read.csv("https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MOD11A2/2005/001.csv")
    tile_names <- substr(hdf_files$name, 18, 23) # tile name is always in the same position in MODIS HDF names
    
    tiles <- data.frame(tilesH = as.numeric(substr(tile_names, 2,3)),
                        tilesV = as.numeric(substr(tile_names, 5,6)))
    return(tiles)
}

# mosaics_parameters
# Return a data.frame containing all the possible combinations of Product, Year, Month, Day_or_Night.
# The data.frame starts at the earliest month where HDF are available for this month 
# (March 2000 for MOD11A2, and August 2002 for MYD11A2) and stops one month before the 
# current month.
mosaics_parameters <- function(){
    current_date <- as.character(Sys.Date())
    current_month <- as.numeric(format(Sys.Date(), "%m"))
    current_year <- as.numeric(format(Sys.Date(), "%Y"))

    # Creates all the combinations of parameters
    combinations <- expand.grid(Product = c('MOD11A2', "MYD11A2"), Year = c(2000:current_year), 
                                Month = seq(12), Day_or_Night = c('Day', 'Night'))
    
    # Remove inexisting combinations (before 2000/03 and 2002/08 for MOD and MYD, 
    # and after the last month for both)
    combinations <- subset(combinations, (!(Year == current_year & Month >= current_month)))
    combinations <- subset(combinations, (!(Product == 'MOD11A2' & Year == 2000 & Month < 3)))
    combinations <- subset(combinations, (!(Product == 'MYD11A2' & Year < 2002)))
    combinations <- subset(combinations, (!(Product == 'MYD11A2' & Year == 2002 & Month < 7)))
    return(combinations)
}

# raster_high_qc
# Takes an LST and a QC raster and return a raster with LST values only for pixel
# with the highest quality (QC equal to "00 Pixel produced, good quality, not necessary
# to examine more detailed QA".).
# QC information is stored as an 8 bit integer, meaning that the first 2 bits of that integer should be 0.
# Note that since HDF files are encoded in big endian, the bits are read from right to left, 
# meaning that in reality, one should first flip the 8 bit integer, and then select the pixels based 
# upon the last two bits or flip the first two bits (this does not matter here since 00 fliped is
# still 00).
# Note that there is a very large percentage of data with subpar quality (i.e. not 00) which is
# not consistent with the metadata associated (for example, the metadata for MOD11A2.A2005177.h09v05.006
# indicates that 76% of the pixel have a good quality while the QC layer indicates that there is none).
# The percentage of cell with NA in the QC and a value in LST is roughly equivalent to this percentage in this case.
# Since I coulnd't explore further this discrepancy, I added a quick fix in the form of
# an option that considers the NAs in the QC layers as 00.
#
# Parameters:
# LST: raster containing LST information
# QC: raster containing the QC information for the LST raster
raster_high_qc <- function(LST, QC, goodQC_equals_NA = FALSE){
    # unique integer present in the QC raster
    unique_int <- na.omit(unique(QC[]))
    
    # position of the integer among unique_int that have the bits 1-0 == 00
    high_QC <- which((colSums(sapply(unique_int, 
                                     function(x) 
                                         as.integer(intToBits(x))[2:1]) == 0) - 1) == 1) 
    
    # vector of the integer with last two bits == 00
    int_high_QC <- unique_int[high_QC] 
    
    # creating a mask where pixel with high QC == 1
    QC_mask <- QC
    QC_mask[][QC[] %in% int_high_QC] <- 1
    QC_mask[][QC_mask[] != 1 & !is.na(QC_mask[])] <- 0
    
    if(goodQC_equals_NA == TRUE)
        QC_mask[][is.na(QC_mask[])] <- 1
    
    # Applying the mask to the LST raster to keep only the high QC pixels
    LST_high <- LST * QC_mask
    return(LST_high)
}

# monthly_average_tile
# For a given product and time of the day, at a given date (month within year) and a given location,
# download all the tiles for the month, select only the pixels with the highest QC,
# and creates a monthly composite.
# This function could be modified to implement a form of temporal gap filling by (1) 
# downloading HDF before and after the beginning of the month, (2) creating cloud masks,
# (3) for each raster within the month, averaging the LST value of the cloudy pixel 
# from the rasters located before and after.
#
# Parameters:
# product: charachter. Name of the MODIS product ("MOD11A2" or "MYD11A2")
# year: numeric. Year selected
# month: numeric. Month for which to  create a composite
# tileH: numerc. Horizontal number of the tile.
# tileV: numeric. Vertical numer of the tile.
# Day_or_Night: character. Time of the day for which the composite is computed ("Day" or "Night")
# tp: character. Path to the temporary directory where HDF and temporary GeoTiff are stored
monthly_average_tile <- function(product, year, month, tileH, tileV, Day_or_Night, tp){
    begin_date <- paste(year, month, "01", sep = ".") # first day of the month
    days_in_month <- monthDays(as.Date(paste(year, month, 1, sep="-"))) # max. number of days in the month
    end_date <- paste(year, month, days_in_month, sep = ".") # last day of the month
    
    # defining which layer should be downloaded (LST_Day + QC_day or LST_night + QC_night)
    if(Day_or_Night == "Day")
        sds <- "11"
    else
        sds <- "000011"
    
    # defining the directories in which files are temporary stored (files are actually stored in subdirectories)
    original_data_path <- file.path(tp, "MODIS_ARC")
    processed_data_path <- file.path(tp, "MODIS_ARC/PROCESSED")
    
    # changing the path whithin the options of the package MODIS
    MODISoptions(localArcPath = original_data_path)
    MODISoptions(outDirPath = processed_data_path)
    
    # download and store as GeoTiff (in a tempdir) the LST and QC layer for the product and for day or night.
    # the try function allows the function to keep going if a tile does not exist
    try(runGdal(product, begin = begin_date, end = end_date,
                tileH = tileH, tileV = tileV, SDSstring = sds,
                job = paste(product, year, month, "h", 
                            tileH, "v", tileV, sep = "."), stuborness = 10, wait = 10))
    
    
    # For each file in the tempdir, only the pixel with a QC equal to "00 Pixel produced, 
    # good quality, not necessary to examine more detailed QA" are kept.
    rasters_dir <- file.path(processed_data_path, paste(product, year, month, "h", 
                                                        tileH, "v", tileV, sep = "."))
    dates <- dir(rasters_dir)
    
    # stops the function if the tile does not exists
    if(length(dates) == 0){
        unlink(rasters_dir, recursive = T)
        return()
    }

    dates <- unique(substr(dates, 14, 16)) # defining the date identifier (julian day)
    
    # creates a list of every 8 days rasters in the month with only the pixel passing the QA criterion
    list_lst_high <- sapply(dates, function(r.date){
                        LST <- raster(file.path(rasters_dir, paste0(product, ".A", year, r.date, ".",
                                                        "LST_", Day_or_Night, "_1km.tif")))
                        QC <- raster(file.path(rasters_dir, paste0(product, ".A", year, r.date, ".",
                                                       "QC_", Day_or_Night, ".tif")))
            
                        raster_high_qc(LST, QC, goodQC_equals_NA = TRUE)
        })
    
    stacked_lst <- stack(list_lst_high)
    monthly_lst <- mean(stacked_lst)
    
    #emptying the directories original_data_path and processed_data_path
    unlink(rasters_dir, recursive = T)
    do.call(unlink, list(dir(original_data_path, full.names = T), recursive =T))
    
    # Returning the monthly composite for the tile
    return(monthly_lst)
}

# monthly_mosaics_lst
# Creates a mosaic of a product (either night or day) for a given month within a year.
# Acts as a wrapper for monthly_average_tile() and allows for multicore processing
#
# Parameters:
# product: charachter. Name of the MODIS product ("MOD11A2" or "MYD11A2")
# year: numeric. Year selected
# month: numeric. Month for which to  create a composite
# tilesH: numeric. A vector of the horizontal tiles to be processed. Should be the 
#           same length as tilesV. For example, if tilesH = c(1, 2) and tilesV = c(5,5),
#           only the tiles H1V5 and H2V5 will be processed.
# tileV: numeric. A vector of the vertical tiles to be downloaded. Should be the 
#           same length as tilesH
# Day_or_Night: character. Time of the day for which the composite is computed ("Day" or "Night")
# ncores: numeric. If ncores > 1, multicore threading will be used on a single machine
# cluster_type: character. Type of cluster ('SOCK', 'MPI', 'PVM' or 'NWS'.
#               Default to 'SOCK', but using 'MPI' should allow HPC (Note that this has not been tested yet).
# save_path: character. A path pointing to the directory where the composite geotiff should be saved.
#               If not specified, the raster is returned as an object in the global environment.
#               If specified, the raster is saved as a geotiff, follwing this nomenclature:
#                       Product.day_or_night.year.month.tiff
monthly_mosaic_lst <- function(product, month, year, Day_or_Night, tilesH, 
                               tilesV, ncores = NULL, cluster_type = 'SOCK',
                               save_path = NULL){
    if(length(tilesV) != length(tilesH))
        stop("tilesH and tilesV should be of equal length!")
    
    tp <- tempdir()
    shell.exec(tp)
    
    # Single core computation of a list of tile of different spatial extent. Each tile is a monthly composite.
    if(is.null(ncores)){
        list_lst <- sapply(seq(length(tilesV)), function(x){
            print(paste(tilesH[x], tilesV[x]))
            monthly_average_tile(product = product, year = year, month = month,
                                 tileH = tilesH[x], tileV = tilesV[x], 
                                 Day_or_Night = Day_or_Night, tp =tp)
        })
    }
    
    #Multicore computation of a list of tile of different spatial extent. Each tile is a monthly composite
    if(!is.null(ncores)){
        sfInit(parallel = T, cpus = ncores, type = cluster_type)
        sfExport(list = c("product", "year", "month", "tilesH", "tilesV",
                      "Day_or_Night", "monthly_average_tile",
                      "raster_high_qc", "tp"))
        sfLibrary(raster)
        sfLibrary(Hmisc)
        sfLibrary(MODIS)
        list_lst <- sfLapply(seq(length(tilesV)), function(x){
            #print(paste(tilesH[x], tilesV[x]))
            monthly_average_tile(product = product, year = year, month = month,
                                 tileH = tilesH[x], tileV = tilesV[x], 
                                 Day_or_Night = Day_or_Night, tp =tp)
        })
        
    }
    
    cat("\nMosaicing the tiles. Don't exit!")
    # Remove any NULL item in the list (created if non existing tile number were used)
    list_lst_non_null <- list_lst[which(sapply(list_lst, function(x) !is.null(x)))]
    
    # Creates a mosaic of monthly composite
    monthly_lst <- do.call(mosaic, c(list_lst_non_null, fun = mean))
    
    if(!is.null(save_path)){
        if(!(dir.exists(save_path))){
            cat("\nSaving directory does not exists, it will be created")
            dir.create(save_path, recursive = T)
        }
            
        writeRaster(monthly_lst, 
                    file.path(save_path, paste0(product, ".006.A", year, month, 
                                                ".", Day_or_Night, ".tif")),
                    format = "GTiff")
    }

    else
        return(monthly_lst)
    
}

# mosaics
# Wrap up function that creates multiple mosaics defined as combinations of product, year, month,
# time of the day
#
# Parameters:
# data.frame. DF containing the required combinations of parameters. Columns must be named
#               'Product', 'Year', 'Month' and 'Day_or_Night'
#               Can be produced by mosaic_parameters()
# tiles: data.frame. Contains all the tiles identifiers (one row per tile). Columns must be named
#               'tilesH' and 'tilesV.
#               Can be produced by existing_tiles()
# ncores: numeric. If ncores > 1, multicore threading will be used on a single machine
# cluster_type: character. Type of cluster ('SOCK', 'MPI', 'PVM' or 'NWS'.
#               Default to 'SOCK', but using 'MPI' should allow HPC (Note that this has not been tested yet).
# save_path: character. A path pointing to the directory where the mosaic geotiff should be saved.
#               If not specified, the raster is returned as an object in the global environment.
#               Not recommended as it will dramactically increase memory usage!
#               If specified, the raster is saved as a geotiff, follwing this nomenclature:
#                       Product.day_or_night.year.month.tiff
mosaics <- function(save_path, tiles, parameters, ncores = NULL, cluster_type = 'SOCK'){
    apply(parameters, 1, function(row){
        product <- row['Product']
        month <- row['Month']
        year <- row['Year']
        Day_or_Night <- row['Day_or_Night']
        
        monthly_mosaic_lst(product= product, year = year, month = month, tilesH = tiles$tilesH,
                           tilesV = tiles$tilesV, Day_or_Night = Day_or_Night, ncores = ncores,
                           cluster_type = cluster_type, save_path = save_path)
    })
}

