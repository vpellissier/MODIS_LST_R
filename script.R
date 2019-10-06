# Tested with i7-3740QM @ 2.70GHz / Win7 with 16 GB RAM, with R 3.5.1

# PRE-REQUISITES (Windows):
# This needs a version of GDAL , preferably from OSGeo4w.
# This GDAL is installed if you have GRASS or QGIS installed, otherwise it can be found
# here http://trac.osgeo.org/osgeo4w.
# The path might have to be specified with MODIS::MODISoptions(gdalPath = "path/to/gdal/bin")
# The script also requires to have curl and/or wget installed and their path added to the environment
# variables.

# Alogrithm:
# The functions are commented in details (in the functions.R file), but the algorithm 
# works with four main functions:
# * raster_high_qc() takes a LST and a QC raster as input, and return a raster with only the
#     highest quality pixels
# * monthly_average_tile() download all the HDFs for a tile for a unique combination of parameters
#     (product * year * month * day_or_night and saves the matching LST and QC layers as 
#     temporary GeoTiff.
#     The GeoTiff are then passed to raster_high_qc() to keep only the highest quality pixels
#     and averaged to create a monthly composite of the tile.
# * monthly_mosaic_lst() creates a monthly mosaic of several tiles for a unique combination of 
#     parameters (product * year * month * day_or_night) and saves it as a GeoTiff
# * mosaics() creates monthly mosaics of several tiles for several combinations of parameters and
#     saves them as GeoTiff
#
# NOTE: 
# The script can be use in multicore by building a tile monthly composite in each core.
# Only one tile is computed on each core to avoid using too much RAM (and too large overhead)
# HDF and temporary GeoTiff are removed from the hard drive after they have been used to avoid 
# excessive HD usage.
# Moreover, unique combinations of parameters and tiles identifiers can be generated automatically
# using functions described in the function.R file, meaning that this function can be run on a 
# cluster using a bash script. The only variables to pass to each node are (1) a path to save 
# the final raster and (2) an integer between 1 and the number of unique combination of 
# parameters (here 832).

# Installing the packages needed
pkgs <- c("MODIS", "raster", "snowfall", "Hmisc")
sapply(pkgs, function(pkg){
    if(!(pkg %in% pkgs))
       install.packages(pkg)
})

# Sourcing the functions (the functions are commented in the functions.R file)
project_path <- "path/to_the_project_directory"
source(file.path(project_path, "functions.R"))

# Creating a dataframe containing the identifiers of all possible tiles
tiles <- existing_tiles()

# Creating a data.frame containing all the possible combinations of parameters (product, year, month,
# day_or_night)
parameters <- mosaics_parameters()

# Settings the login informations for the NASA server (only done once per machine)
MODIS::EarthdataLogin(usr = 'username', pwd = "password")

# Computing and saving the mosaics for all parameters combinations found in parameters and for all
# the tiles found in tiles
mosaic_path <- "path/mosaic_dir"
mosaics(save_path = mosaic_path, tiles = tiles, parameters = parameters,
        ncores = NULL, cluster_type = 'SOCK')





