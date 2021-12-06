rm( list = ls())

source( "~/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R")

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

#======================================================================#
## Load meteorology as list of months
#======================================================================#
#define the layer names, do the actual downloading
Sys.setenv(TZ='UTC')
layer.names <- c( "air.2m.mon.mean.nc",
                  "apcp.mon.mean.nc",
                  "rhum.2m.mon.mean.nc",
                  "vwnd.10m.mon.mean.nc",
                  "uwnd.10m.mon.mean.nc")
names( layer.names) <- c( "temp", "apcp", "rhum", "vwnd", "uwnd")

# do the data downloading
# set destination parameter to where you want the data downloaded,
# for example, destination = '~/Desktop'
list.met <- lapply( layer.names,
                    downloader.fn, 
		    destination = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/met',
                    dataset = 'NARR')


# take over US
mets2005   <- usa.functioner( 2005, list.met, dataset = 'NARR', return.usa.sub = F)
mets2006   <- usa.functioner( 2006, list.met, dataset = 'NARR', return.usa.sub = F)
mets2011   <- usa.functioner( 2011, list.met, dataset = 'NARR', return.usa.sub = F)

#======================================================================#
## Load ddm as annual
#======================================================================#
ddm2005 <- ddm_to_zip( ddm_coal_file = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/COAL_impacts_2005_update.csv',
                       Year = 2005)
ddm2006 <- ddm_to_zip( ddm_coal_file = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/COAL_impacts_2006_update.csv',
                       Year = 2006)
names( ddm2005) <- 'cmaq.ddm'
names( ddm2006) <- 'cmaq.ddm'

#======================================================================#
## Load anuual hyads
#======================================================================#
hyads2005.dt <- read.fst( '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp/coal_hyads_grids/grids_exposures_total_2005.fst', as.data.table = TRUE)
hyads2006.dt <- read.fst( '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp/coal_hyads_grids/grids_exposures_total_2006.fst', as.data.table = TRUE)
hyads2011.dt <- read.fst( '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp/coal_hyads_grids/grids_exposures_total_2011.fst', as.data.table = TRUE)

hyads2005 <- rasterFromXYZ( hyads2005.dt[, .( x, y, hyads)], crs = p4s)
hyads2006 <- rasterFromXYZ( hyads2006.dt[, .( x, y, hyads)], crs = p4s)
hyads2011 <- rasterFromXYZ( hyads2011.dt[, .( x, y, hyads)], crs = p4s)

## ========================================================= ##
##                Read in emissions data
## ========================================================= ##
d_cems_cmaq.f <- "/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/2005_cemsum.txt"
d_nonegu.f <- "/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/ptinv_ptnonipm_xportfrac_cap2005v2_2005cs_orl_06jan2011_v4_orl_COAL.txt"

d_cmaq <- fread( d_cems_cmaq.f)
d_nonegu <- fread( d_nonegu.f, skip = "06029", header = F)[,1:63]
d_nonegu.names <- unlist( fread( d_nonegu.f, skip = 'FIPS,PLANTID,', header = F, nrows = 1))
names( d_nonegu) <- d_nonegu.names
d_nonegu.slim <- d_nonegu[ POLCODE == 'SO2', .( XLOC, YLOC, ANN_EMIS)]

## Convert to spatial object, take over CMAQ raster
d_nonegu.sp <- SpatialPointsDataFrame( d_nonegu.slim[, .( XLOC, YLOC)], 
                                       data.frame( d_nonegu.slim[, ANN_EMIS]),
                                       proj4string = CRS( "+proj=longlat +datum=WGS84 +no_defs"))
d_nonegu.sp <- spTransform( d_nonegu.sp, CRS( p4s))
d_nonegu.r <- rasterize( d_nonegu.sp, ddm.m.all)$d_nonegu.slim...ANN_EMIS.
d_nonegu.r[is.na(d_nonegu.r[])] <- 0

## ========================================================= ##
##                SOx inverse distance by year
## ========================================================= ##
idwe2005.dt <- fread( '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/ampd_dists_sox_weighted_2005_total.csv', drop = 'V1')
idwe2006.dt <- fread( '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/ampd_dists_sox_weighted_2006_total.csv', drop = 'V1')
idwe2011.dt <- fread( '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/ampd_dists_sox_weighted_2011_total.csv', drop = 'V1')
idwe2005 <- rasterFromXYZ( idwe2005.dt, crs = p4s)
idwe2006 <- rasterFromXYZ( idwe2006.dt, crs = p4s)
idwe2011 <- rasterFromXYZ( idwe2011.dt, crs = p4s)
names( idwe2005) <- 'idwe'
names( idwe2006) <- 'idwe'
names( idwe2011) <- 'idwe'

## ========================================================= ##
##                Plots
## ========================================================= ##
# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

#======================================================================#
# stack up and project annual data
#======================================================================#
dats2005.a <- project_and_stack( ddm2005, hyads2005, idwe2005, 
                                 mets2005, d_nonegu.r, mask.use = mask.usa)
dats2006.a <- project_and_stack( ddm2006, hyads2006, idwe2006, 
                                 mets2006, d_nonegu.r, mask.use = mask.usa)
dats2011.a <- project_and_stack( ddm2006, hyads2011, idwe2011, 
                                 mets2011, d_nonegu.r, mask.use = mask.usa)
dats2011.a$cmaq.ddm <- NA

#======================================================================#
## Combine into raster stack, train model
#======================================================================#
cov.names = c( "temp", "rhum", "vwnd", "uwnd", "wspd")

# predict annual 2006 using model trained in 2005
preds.ann.hyads06w05 <- lm.hyads.ddm.holdout( dat.stack = dats2005.a, dat.stack.pred = dats2006.a, 
					      name.idwe = 'idwe', x.name = 'hyads',
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)
preds.ann.idwe06w05  <- lm.hyads.ddm.holdout( dat.stack = dats2005.a, dat.stack.pred = dats2006.a, 
					      name.idwe = 'idwe', x.name = 'idwe',
                                              ho.frac = 0, covars.names = cov.names, return.mods = T)


#======================================================================#
## Save data
#======================================================================#
# annual stacks, 
# annual model
save( dats2005.a, dats2006.a, dats2011.a,
      d_nonegu.r,
      preds.ann.hyads06w05, #preds.ann.hyads05w06,
      preds.ann.idwe06w05,  #preds.ann.idwe05w06,
      file = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/hyads_to_cmaq_models20210715.RData')





















