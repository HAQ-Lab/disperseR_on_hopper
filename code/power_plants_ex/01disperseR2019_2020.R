# srun -p test --mem 15g -t 0-01:00 -c 1 -N 1 --pty /bin/bash
# devtools::install_github( 'lhenneman/disperseR@dev')

# select years to run. I suggest 5 years at a time
#years.run <- 2000:2004
#years.run <- 2005:2009
#years.run <- 2010:2015
years.run <- 2019:2020

library(disperseR) # our package

# set directory structure
#disperseR::create_dirs('/n/holyscratch01/zigler_lab/lhenneman/diseperseR')
disperseR::create_dirs('/projects/HAQ_LAB/lhennem/disperseR')

# we want proc_dir and hysp_dir to be in scratch
proc_dir <- '/scratch/lhennem/disperseR/main/process'
dir.create(proc_dir, showWarnings = FALSE)

# download data
disperseR::get_data(data = "all",
                    start.year = "1999",
                    start.month = "01",
                    end.year = "2020",
                    end.month = "12")

# pick out units to run -- load updated coal dataset
load( '/projects/HAQ_LAB/mrasel/R/ampd-raw-data-processing/data/units_coal_1997_2021.rda')
#load( '/projects/HAQ_LAB/lhennem/data/disperseR/ampd_rasel/units_coal_1997_2021.rda')

# limit to specific units
unitsrun <- units_updated %>%
  dplyr::filter(year %in% years.run) %>% # only get data for called years
  dplyr::filter( SOx > 0) %>% # only get data for called years
  data.table()

# find unique combos of Latitude, Longitude, and Height
unitslatlonh <- unique( unitsrun[ ,.( Latitude, Longitude, Height, year)] ) 
unitslatlonh[, unit_run_ref:=1:nrow( unitslatlonh)]
unitsrun_trim <- merge( unitsrun, unitslatlonh)[ !duplicated( unit_run_ref)]

# define data.table with all emissions events
input_refs <- disperseR::define_inputs(units = unitsrun_trim,
                                       startday = paste0( years.run[1], '-01-01'),
                                       endday = paste0( years.run[length( years.run)], '-12-31'),
                                       start.hours =  c(0, 6, 12, 18),
                                       duration = 24 * 7)

# select specific 
array_num <- as.numeric( Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_num <- ifelse( array_num == '' | is.na( array_num), 1, array_num)
set.seed( array_num^2 + 51)
refs.use <- sample( 1:nrow( input_refs))

# run disperser
hysp_raw <- disperseR::run_disperser_parallel(input.refs = input_refs[refs.use],
                                              pbl.height = pblheight,
                                              species = 'so2',
                                              proc_dir = proc_dir,
                                              mc.cores = 1)






















