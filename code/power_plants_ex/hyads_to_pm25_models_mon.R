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
mets2005.m <- usa.functioner( 2005, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)
mets2006.m <- usa.functioner( 2006, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)
mets2011.m <- usa.functioner( 2011, list.met, dataset = 'NARR', avg.period = 'month', return.usa.sub = F)

# combine monthly rasters into single list
mets.m.all <- append( append( mets2005.m, mets2006.m), mets2011.m)

#======================================================================#
## Load ddm as monthly
#======================================================================#
ddm2005.m <- ddm_to_zip( ddm_coal_file = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/COAL_impacts_2005_update.csv',
                         Year = 2005, avg.period = 'month')
ddm2006.m <- ddm_to_zip( ddm_coal_file = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/COAL_impacts_2006_update.csv',
                         Year = 2006, avg.period = 'month')

names( ddm2005.m) <- names( mets2005.m)
names( ddm2006.m) <- names( mets2006.m)

# combine into single list
ddm.m.all <- stack( ddm2005.m, ddm2006.m)

#======================================================================#
## Load monthly hyads
#======================================================================#
hyads_dir <- '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp/coal_hyads_grids'
hyads2005.m.f <- list.files( hyads_dir, pattern = 'grids_exposures_total_2005_', full.names = TRUE)
hyads2006.m.f <- list.files( hyads_dir, pattern = 'grids_exposures_total_2006_', full.names = TRUE)
hyads2011.m.f <- list.files( hyads_dir, pattern = 'grids_exposures_total_2011_', full.names = TRUE)

# create lists from monthly grid objects
hyads2005.m.l <- lapply( hyads2005.m.f, read.fst, as.data.table = TRUE)
hyads2006.m.l <- lapply( hyads2006.m.f, read.fst, as.data.table = TRUE)
hyads2011.m.l <- lapply( hyads2011.m.f, read.fst, as.data.table = TRUE)
names( hyads2005.m.l) <- names( mets2005.m)
names( hyads2006.m.l) <- names( mets2006.m)
names( hyads2011.m.l) <- names( mets2011.m)

# create lists of monthly rasters
HyADSrasterizer <- function( X){
  r <- rasterFromXYZ( X[, .( x, y, hyads)], crs = p4s)
  r[is.na( r)] <- 0
  return( r)
}

hyads2005.m <- lapply( hyads2005.m.l, HyADSrasterizer)
hyads2006.m <- lapply( hyads2006.m.l, HyADSrasterizer)
hyads2011.m <- lapply( hyads2011.m.l, HyADSrasterizer)

# combine into single list
hyads.m.all <- stack( stack( hyads2005.m), stack( hyads2006.m), stack( hyads2011.m))

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
idwe.m.dt <- fread( '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/ampd_dists_sox_weighted.csv', drop = 'V1')
idwe.m.l <- split( idwe.m.dt, by = 'yearmon')

# create lists of monthly rasters
IDWErasterizer <- function( X){
  r <- rasterFromXYZ( X[, .( x, y, tot.sum)], crs = p4s)
  r[is.na( r)] <- 0
  return( r)
}

idwe.m <- stack( lapply( idwe.m.l, IDWErasterizer))
names( idwe.m) <- names( hyads.m.all)

## ========================================================= ##
##                Plots
## ========================================================= ##
# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

#======================================================================#
## Combine into raster stack, train model
#======================================================================#
cov.names = c( "temp", "rhum", "vwnd", "uwnd", "wspd")

# predict each month in 2006 using model trained in 2005
preds.mon.hyads06w05 <- mapply( month.trainer, names( mets2005.m), names( mets2006.m),
                                MoreArgs = list( name.x = 'hyads', y.m = hyads.m.all,
                                                 ddm.m = ddm.m.all, mets.m = mets.m.all,
                                                 idwe.m = idwe.m, emiss.m = d_nonegu.r, 
                                                 .mask.use = mask.usa, cov.names = cov.names))
preds.mon.idwe06w05  <- mapply( month.trainer, names( mets2005.m), names( mets2006.m),
                                MoreArgs = list( name.x = 'idwe', y.m = idwe.m,
                                                 ddm.m = ddm.m.all, mets.m = mets.m.all,
                                                 idwe.m = idwe.m, emiss.m = d_nonegu.r, 
                                                 .mask.use = mask.usa, cov.names = cov.names))


#======================================================================#
## Save data
#======================================================================#
# annual stacks, 
# annual model
save( hyads.m.all, ddm.m.all, mets.m.all,
      idwe.m, d_nonegu.r,
      preds.mon.hyads06w05, #preds.mon.hyads05w06,
      preds.mon.idwe06w05,  #preds.mon.idwe05w06,
      file = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/hyads_to_cmaq_models_month20210715.RData')


#======================================================================#
## Plot the metrics
#======================================================================#
## extract evaluation statistics
## IDWE gets big change from bivariate spline, HyADS does not  
preds.metrics.hyads <- preds.mon.hyads06w05[ 'metrics',]
preds.metrics.idwe  <- preds.mon.idwe06w05[ 'metrics',]

metrics <- data.table( month = c( as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.hyads)))),
                                  as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.idwe))))),
                       model = c( rep( 'hyads', length( names( preds.metrics.hyads))),
                                  rep( 'idwe', length( names( preds.metrics.idwe)))),
                       class = c( rep( 'gam', 2 * length( names( preds.metrics.hyads))),
                                  rep( 'lm', 2 * length( names( preds.metrics.idwe)))),
                       'R^2' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$R),
                                  sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$R),
                                  sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$R),
                                  sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$R)),
                       NMB = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$NMB),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$NMB),
                                sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$NMB),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$NMB)),
                      MB = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$MB),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$MB),
                                sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$MB),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$MB)),
                        NME = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$NME),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$NME),
                                sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$NME),
                                sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$NME)),
                       RMSE = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$RMSE),
                                 sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$RMSE),
                                 sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$RMSE),
                                 sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$RMSE)))

metrics[model == 'hyads' & class == 'gam']
metrics.m <- melt( metrics, id.vars = c( 'model', 'month', 'class'), variable.name = 'metric')



















