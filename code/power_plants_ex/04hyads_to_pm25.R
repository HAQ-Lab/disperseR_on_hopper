# srun -p test --mem 50g -t 0-04:00 -c 1 -N 1 --pty /bin/bash
# devtools::install_github( 'lhenneman/disperseR@dev')

# select years to run. I suggest 5 years at a time
years.run <- 1999:2020

array_num <- as.numeric( Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_num <- ifelse( array_num == '' | is.na( array_num), 1, array_num)
years.run.sel <- years.run[array_num]

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

# load import functions/datasets
source( '~/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
# source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
load( '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/hyads_to_cmaq_models20210715.RData')
load( '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/hyads_to_cmaq_models_month20210715.RData')

# point to hyads exposures (raw, unitless)
hyads.dir <- '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp/coal_hyads_grids'
hyadsPM25.dir <- '/projects/HAQ_LAB/lhennem/disperseR/main/output/exp25/grids'
# hyads.dir <- '/Users/lhenneman/Dropbox/Harvard/RFMeval_Local/HyADS/disperseR_hyads'
# hyadsPM25.dir <- '~/Google Drive/Harvard/HyADS/'

# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

#======================================================================#
## Annual conversions
#======================================================================#
#mapply( hyads_to_pm25_unit, rep( 2013:2015, each = 12), 1:12,
lapply( years.run.sel,
        hyads_to_pm25_unit,
        fstart = file.path( hyads.dir, 'grids_exposures_byunit_'),
        fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_byunit_'),
        model.dataset = preds.ann.hyads06w05,
        model.name = 'model.cv', # 'model.gam'
        name.x = 'hyads',
        met.dest = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/met',
        mask.use = mask.usa
)

#======================================================================#
## Annual conversions -- total
#======================================================================#
#mapply( hyads_to_pm25_unit, rep( 2013:2015, each = 12), 1:12,
lapply( years.run.sel,
        hyads_to_pm25_unit,
        fstart = file.path( hyads.dir, 'grids_exposures_total_'),
        fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_total_'),
        model.dataset = preds.ann.hyads06w05,
        model.name = 'model.cv', # 'model.gam'
        name.x = 'hyads',
        mask.use = mask.usa,
        met.dest = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/met',
        total = T
)

#======================================================================#
## Monthly conversions
#======================================================================#
mapply( hyads_to_pm25_unit, rep( years.run.sel, each = 12), 1:12,
        MoreArgs = list(
                fstart = file.path( hyads.dir, 'grids_exposures_byunit_'),
                fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_byunit_'),
                model.dataset = preds.mon.hyads06w05,
                model.name = 'model.cv', # 'model.gam'
                name.x = 'hyads',
                mask.use = mask.usa,
		met.dest = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/met'
        )
)

#======================================================================#
## Monthly conversions -- total
#======================================================================#
#mapply( hyads_to_pm25_unit, rep( 2013:2015, each = 12), 1:12,
mapply( hyads_to_pm25_unit, rep( years.run.sel, each = 12), 1:12,
        MoreArgs = list(
                fstart = file.path( hyads.dir, 'grids_exposures_total_'),
                fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_total_'),
                model.dataset = preds.mon.hyads06w05,
                model.name = 'model.cv', # 'model.gam'
                name.x = 'hyads',
                mask.use = mask.usa,
                total = T,
                met.dest = '/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/met'
        )
)



