# srun -p test --mem 50g -t 0-04:00 -c 1 -N 1 --pty /bin/bash
# devtools::install_github( 'lhenneman/disperseR@dev')

# select years to run. I suggest 5 years at a time
years.run <- 2019:2020

library(disperseR) # our package

# set directory structure
disperseR::create_dirs('/projects/HAQ_LAB/lhennem/disperseR')

# download data
disperseR::get_data(data = "all",
                    start.year = "1999",
                    start.month = "01",
                    end.year = "2020",
                    end.month = "12")

# read in emissions data
load( '/projects/HAQ_LAB/mrasel/R/ampd-raw-data-processing/data/units_coal_1997_2021.rda')
load( '/projects/HAQ_LAB/mrasel/R/ampd-raw-data-processing/data/units_monthly.rda')

# limit to specific unitsi
units_updated <- units_updated %>%
  mutate(uID=gsub("-", ".", ID))

unitsrun <- units_updated %>%
  dplyr::filter(year %in% years.run) %>% # only get data for called years
  dplyr::filter( SOx > 0) %>% # only get data for called years
  data.table()

# find unique combos of Latitude, Longitude, and Height
unitslatlonh <- unique( unitsrun[ ,.( Latitude, Longitude, Height, year)] )
unitslatlonh[, unit_run_ref:=1:nrow( unitslatlonh)]
unitsrun_notrim <- merge( unitsrun, unitslatlonh)
unitsrun_trim <- merge( unitsrun, unitslatlonh)[ !duplicated( unit_run_ref)]

# define yearmons for link to zips/grids/counties
yearmons <- disperseR::get_yearmon(start.year = as( years.run[1], 'character'),
                                   start.month = "01",
                                   end.year = as( years.run[length( years.run)], 'character'),
                                   end.month = '12')

# name the rdata file
rdata_file <- 'disperseR_grids_2019-2020.RData'

# need to:
#  read in combined_gridlinks
#  add in missing units (copy identical ones)
# save updated gridlinks
load( file.path( rdata_dir, rdata_file))
edits <- lapply( yearmons, 
  function( ym){
    print( ym)
    year.ym <- as( substr( ym, 1, 4), 'numeric')
    mont.ym <- substr( ym, 5, 6)
    name.ym <- paste0( "MAP", mont.ym, ".", year.ym)

    map.ym <- get( name.ym)

    unitsrun_use <- unitsrun_notrim[year == year.ym]
    units.map <- names( map.ym)[ !(names( map.ym) %in% c( 'ZIP', 'x', 'y'))]

    changes <- lapply( units.map, 
      function( u){
        u_tmp <- unitsrun_use[ID == u]
        if( nrow( u_tmp) == 0) {
          map.ym[, (u) := NULL]
          return( -1)
        }
        u_use <- unitsrun_use[ID != u & unit_run_ref == u_tmp$unit_run_ref]$ID
        if( length( u_use) == 0) return(0)
        map.ym[, (u_use) := lapply( u_use, function(u.) unlist( map.ym[, u, with = F]))]
        return( 1)
      })
    assign( name.ym, map.ym, envir = .GlobalEnv)
    return( changes)
})

# weight by emissions
month.locs.u <- lapply( years.run, function( y) 
  disperseR::calculate_exposure( year.E = y,
                                 year.D = y,
                                 link.to = 'grids',
                                 pollutant = 'SO2.tons',
                                 rda_file = "loaded",
                                 exp_dir = exp_dir,
                                 units.mo = units_monthly,
                                 source.agg = 'unit',
                                 time.agg = 'month',
                                 return.monthly.data = F)
)


month.locs <- lapply( years.run, function( y)
  disperseR::calculate_exposure( year.E = y,
                                 year.D = y,
                                 link.to = 'grids',
                                 pollutant = 'SO2.tons',
                                 rda_file = "loaded",
                                 exp_dir = exp_dir,
                                 units.mo = units_monthly,
                                 source.agg = 'total',
                                 time.agg = 'month',
                                 return.monthly.data = T)
)

year.locs <- lapply( years.run, function( y)
  disperseR::calculate_exposure( year.E = y,
                                 year.D = y,
                                 link.to = 'grids',
                                 pollutant = 'SO2.tons',
                                 rda_file = "loaded",
                                 exp_dir = exp_dir,
                                 units.mo = units_monthly,
                                 source.agg = 'total',
                                 time.agg = 'year',
                                 return.monthly.data = F)
)

year.locs.u <- lapply( years.run, function( y)
  disperseR::calculate_exposure( year.E = y,
                                 year.D = y,
                                 link.to = 'grids',
                                 pollutant = 'SO2.tons',
                                 rda_file = "loaded",
                                 exp_dir = exp_dir,
                                 units.mo = units_monthly,
                                 source.agg = 'unit',
                                 time.agg = 'year',
                                 return.monthly.data = F)
)








