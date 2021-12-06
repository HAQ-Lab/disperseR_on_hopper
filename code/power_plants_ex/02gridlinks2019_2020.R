# srun -p test --mem 15g -t 0-04:00 -c 1 -N 1 --pty /bin/bash
# devtools::install_github( 'lhenneman/disperseR@dev')

# select years to run. I suggest 5 years at a time
#years.run <- 2000:2005
#years.run <- 2005:2009
#years.run <- 2010:2015
years.run <-  2019:2020

library(disperseR) # our package

# set directory structure
disperseR::create_dirs('/projects/HAQ_LAB/lhennem/disperseR')

# hysplit output copied to scratch
#hysp_dir <- '/scratch/lhennem/disperseR/main/output/hysplit'

# download data
disperseR::get_data(data = "all",
                    start.year = "2000",
                    start.month = "01",
                    end.year = "2020",
                    end.month = "12")

# pick out units to run -- load updated coal dataset
load( '/projects/HAQ_LAB/mrasel/R/ampd-raw-data-processing/data/units_coal_1997_2021.rda')
#load( '/projects/HAQ_LAB/lhennem/data/disperseR/ampd_rasel/units_coal_1997_2021.rda')

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
unitsrun_trim <- merge( unitsrun, unitslatlonh)[ !duplicated( unit_run_ref)]

# define yearmons for link to zips/grids/counties
yearmons <- disperseR::get_yearmon(start.year = as( years.run[1], 'character'),
                                   start.month = "01",
                                   end.year = as( years.run[length( years.run)], 'character'),
                                   end.month = '12')

# select specific 
array_num <- as.numeric( Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_num <- ifelse( array_num == '' | is.na( array_num), 1, array_num)
set.seed( array_num)

# run linker
# 2721-5, 676-3
unitsuse <- unitsrun_trim[!duplicated(ID)]
linked_grid <- disperseR::link_all_units(units.run = unitsuse[c( array_num, sample( 1:nrow( unitsuse)))],
                                         link.to = 'grids',
                                         mc.cores = 1,
                                         year.mons = yearmons,
                                         pbl.height = pblheight,
                                         crosswalk. = crosswalk,
                                         duration.run.hours = 24*7,
                                         res.link = 36000,
                                         overwrite = FALSE,
                                         pbl.trim = FALSE,
					 crop.usa = TRUE,
                                         return.linked.data = FALSE)

# check for corrupt files

nums <- lapply( 1:nrow( unitsuse),
  function( num){
    rerun <- FALSE
    print( paste( 'testing input_refs row', num))
    subset <- unitsuse[!grep('\\*', ID)][num]
    
    unitfiles <- list.files( ziplink_dir, 
                             pattern = paste0( '_', subset$ID, '_'),
                             full.name = T)
    u.t <- rep( FALSE, length( unitfiles))
    for( uf in seq_along( unitfiles)){
      file <- unitfiles[uf]
      print( file)
      o <- system( paste0( "Rout=`Rscript -e 'fst::read.fst(\"",
                           file, "\", from = 1, to = 1)'`"))
    
    if( !file.exists( file) | o != 0){
      u.t[uf] <- TRUE
    } else
      in.fst <- tryCatch( read.fst( file))
        if( is.na( min( in.fst$y)) | max( in.fst$y) > 1e9 | 
              max( in.fst$x) > 1e9 | is.na( mean( in.fst$N)) | 
              T %in% (in.fst$y %% 1000 != 0)| T %in% (in.fst$x %% 1000 != 0))
          u.t[uf] <- TRUE
    }

    if( (TRUE %in% grepl( TRUE, u.t))){
      unlink( unitfiles[u.t])

      linked_grid <- disperseR::link_all_units(units.run = subset,
                                         link.to = 'grids',
                                         mc.cores = 1,
                                         year.mons = yearmons,
                                         pbl.height = pblheight,
                                         crosswalk. = crosswalk,
                                         duration.run.hours = 24*7,
                                         res.link = 36000,
                                         overwrite = FALSE,
                                         pbl.trim = FALSE,
                                         crop.usa = TRUE,
                                         return.linked.data = FALSE)
     }
  })



# combine monthly links
if( array_num == 1)
  combined_gridlinks <- disperseR::combine_monthly_links(month_YYYYMMs = yearmons,
                                                         link.to = 'grids',
                                                         filename = 'disperseR_grids_1999_2016-2018.RData')

# need to:
#  read in combined_gridlinks
#  add in missing units (copy identical ones)
# save updated gridlinks


















