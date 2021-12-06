# srun -p test --mem 50g -t 0-04:00 -c 1 -N 1 --pty /bin/bash
library( sf)
library( raster)
library( data.table)
library( fst)
library( areal)

# select years to run. I suggest 5 years at a time
years.run <- 1999:2020

array_num <- as.numeric( Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_num <- ifelse( array_num == '' | is.na( array_num), 1, array_num)
years.run.sel <- years.run[array_num]

# get usa mask for masking
# download USA polygon from rnaturalearth
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

## ==================================================== ##
# read zcta shapefile and crosswalk
## ==================================================== ##
zip_sf_reader <- function( d = direct.dat){
  # zcta file downloaded from 'ftp://ftp2.census.gov/geo/tiger/GENZ2016/shp/cb_2016_us_zcta510_500k.zip'
  zcta_shapefile <- file.path( d, 'cb_2017_us_zcta510_500k.shp')
  # zcta_shapefile <- file.path( d, 'cb_2016_us_zcta510_500k.shp')

  # zcta-ZIP crosswalk file downloaded from 'http://mcdc2.missouri.edu/data/corrlst/'
  cw <- disperseR::crosswalk
  # cw <- fread( crosswalk_csv)
  # make sure ZCTA's are 5 digits to merge on zcta ID
  # cw$ZCTA <- formatC( cw$ZCTA, width = 5, format = "d", flag = "0")

  zips <- st_read(zcta_shapefile)
  zips <- st_transform( zips, p4s)
  setnames( zips, 'ZCTA5CE10', 'ZCTA')
  zips <- merge( zips, cw, by = "ZCTA", all = F, allow.cartesian = TRUE)
  # make sure ZIPs are 5 digits to merge on zcta ID
  zips$ZIP <- formatC( zips$ZIP, width = 5, format = "d", flag = "0")

  return( zips)
}

## ==================================================== ##
##          load counties - just first 30 in TX
## ==================================================== ##
direct.dat <- '/projects/HAQ_LAB/lhennem/disperseR/main/input/zcta_500k/'
# direct.dat <- '~/Dropbox/Harvard/Manuscripts/Energy_Transitions/Data_and_code/data/gis'
zips <- zip_sf_reader()[, 'ZIP']

## ==================================================== ##
##          Read in hysplit raster
## ==================================================== ##
grids_to_zips <- function( file.in,
			   path.out){
  print( file.in)
  # read in the file
  in.g <- read.fst( file.in, as.data.table = T)
  in.g[ is.na( in.g)] <- 0
  
  # remove units with all zeros
  in.unames <- names( in.g)[ !( names( in.g) %in% c( 'x', 'y'))]
  hyads.sums <- suppressWarnings( melt( in.g[, lapply( .SD, sum), .SDcols = ( in.unames)]))
  hyads.sums.use <- as.character( hyads.sums$variable)[hyads.sums$value > 0]
  
  # trim the dataset to include only units that are not 0
  in.g.trim <- in.g[, c( 'x', 'y', hyads.sums.use), with = F]
  
  # convert to raster
  in.r <- rasterFromXYZ( in.g.trim)
  crs( in.r) <- p4s
  names( in.r) <- names(in.g.trim)[!( names( in.g.trim) %in% c( 'x', 'y'))]
  
  # convert to sf object
  ncin_spatpoly <- rasterToPolygons( in.r)
  ncin_sf <- st_as_sf( ncin_spatpoly)
  ncin_sf <- st_transform( ncin_sf, crs( zips))
  ncin_sf$GID <- 1:nrow( ncin_sf)
  
  # convert to data table and melt
  ncin.dt <- data.table( ncin_sf)[, geometry := NULL]
  ncin.m <- melt( ncin.dt, id.vars = 'GID',
                  variable.name = 'uID', value.name = 'pm25')
  
  # define small in variable
  ncin.train <- ncin_sf[,c( 'GID')]
  
  # take areas of zips over grids
  weights <- aw_intersect( zips, source = ncin.train, areaVar = "area")
  weights.dt <- data.table( weights)
  
  # take total areas of zip fores
  zips.a <- data.table( ZIP = zips$ZIP,
                        ZIP.area = as.vector( st_area( zips)))
  
  # merge together for weighting dataset
  weights.m <- merge( weights.dt, zips.a, by = 'ZIP')
  weights.m[, areaWeight := area / ZIP.area]
  
  # intersect weights and grids
  ncin.m.intersect <- merge( weights.m, ncin.m, by = 'GID', allow.cartesian = TRUE)
  
  # multiple weights by grid pm25, sum by 
  ncin.m.intersect[, pm25a := areaWeight * pm25]
  zips.pm <- ncin.m.intersect[, .( pm25 = sum( pm25a)), by = .( ZIP, uID)]
  
  # cast for smaller file size
  zips.pm.c <- dcast( zips.pm, ZIP ~ uID, value.var = 'pm25')
  
  # define output file, write
  file.in.base <- basename( file.in)
  file.out.base <- gsub( 'grids_', 'zips_', file.in.base)
  file.out <- file.path( path.out, file.out.base)
  write.fst( zips.pm.c, path = file.out)
  
  return( file.out)
}

## ==================================================== ##
##          Run the function
## ==================================================== ##
# by unit, month
grid.files <- list.files( '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp25/grids',
                          pattern = 'grids_.*_\\d{2}\\.fst',
                          full.names = TRUE)
grids_pm.list <- lapply( grid.files, grids_to_zips, 
			 path.out = '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp25/zips')

# by unit, year
grid.files.yr <- list.files( '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp25/grids',
                             pattern = 'grids_.*\\d{4}\\.fst',
                             full.names = TRUE)
grids_pm.list <- lapply( grid.files.yr, grids_to_zips,
                         path.out = '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp25/zips')

# by total, month
grid.files.tot <- list.files( '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp25/grids',
                              pattern = 'grids_.*total_*\\d{2}\\.fst',
                              full.names = TRUE)
grids_pm.list <- lapply( grid.files.tot, grids_to_zips,
                         path.out = '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp25/zips')

# by total, year
grid.files.tot.yr <- list.files( '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/',
                              pattern = 'grids_.*total_\\d{4}\\.fst',
                              full.names = TRUE)
grids_pm.list <- lapply( grid.files.tot.yr, grids_to_zips,
                         path.out = '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp25/zips')




















