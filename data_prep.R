#### Preparing data for model input ####
######### From Phillips 2020 approach ##
# for oligotrophic lake metabolism ## 
rm(list=ls())
# first, an rstan refresher

library(tidyverse)
library(tidybayes)
library(lubridate)
library(rstan)
library(LakeMetabolizer)
`%!in%` <- Negate(`%in%`)
# options(mc.cores = parallel::detectCores())
setwd("/Users/josephvanderwall/Documents/Research/Data Depository/All data (for R)/")

# read in some data

# read in three masterfiles - oxygen 
ox.mas <- read_csv('Summerfiles/Masterfiles/doobs_master.csv')
# check sensor start times
ox.mas %>%
  group_by(lake) %>%
  summarise(startdate = datetime[1])

# write_csv(x = ox.mas %>%  # fixing three lakes with messed up datetimes (timezone issues)
#             mutate(datetime = ifelse(lake %in% c('HOL','UPH'),
#                                       strptime(datetime,format = '%Y-%m-%d %H:%M:%S') - 3*60*60,
#                                       strptime(datetime,format = '%Y-%m-%d %H:%M:%S') + 0 )) %>%
#             mutate(datetime = as.POSIXct(datetime, origin=origin) - 6*60*60),
#           file = 'Summerfiles/Masterfiles/doobs_master.csv')

# temperature
t.mas <- read_csv('Summerfiles/Masterfiles/wtr_master.csv')
# check sensor start dates
t.mas %>%
  group_by(lake) %>%
  summarise(startdate = datetime[1])
t.mas %>%
  filter(depth_m == 1) %>%
  group_by(lake) %>%
  summarise(startdate = datetime[1])

# write_csv(x = t.mas %>%  # fixing three lakes with messed up datetimes (timezone issues)
#             mutate(datetime = ifelse(lake %in% c('HOL','UPH','SAP'),
#                                       strptime(datetime,format = '%Y-%m-%d %H:%M:%S') - 6*60*60,
#                                       strptime(datetime,format = '%Y-%m-%d %H:%M:%S') + 0 )) %>%
#             mutate(datetime = as.POSIXct(datetime, origin=origin) - 6*60*60),
#           file = 'Summerfiles/Masterfiles/wtr_master.csv')

# light
i.mas <- read_csv('Summerfiles/Masterfiles/irr_master.csv')
i.mas <- i.mas %>%
  mutate(par = lux/42) %>% # rough approximation of lux to par conversion (from Thimijan and Heins 1983)
  dplyr::select(-lux)
# check sensor start dates
i.mas %>%
  group_by(lake) %>%
  summarise(startdate = datetime[1])

  
## save incase temp data is similarly fucked - can DELETE
# write_csv(x = i.mas %>%  # fixing four lakes with messed up datetimes (timezone issues)
#             mutate(datetime = ifelse(lake %in% c('ASH','CLW','SMT','RNY'),
#                                       strptime(datetime,format = '%Y-%m-%d %H:%M:%S') + 6*60*60,
#                                       strptime(datetime,format = '%Y-%m-%d %H:%M:%S') + 0 )) %>%
#             mutate(datetime = as.POSIXct(datetime, origin=origin) - 6*60*60),
#           file = 'Summerfiles/Masterfiles/irr_master.csv')

# metadata file with elevation etc.
meta.dat <- read_csv('Metadata/buoy_metadata.csv')
# wind speed distributions
wnd.dist <- read_csv('windspeed_simulation.csv') %>%
    pivot_wider(names_from = 'type',
                values_from = 'value')

## 'Clean' oxygen and temperature data by removing data 2 standard deviations away from daily means
ox.clean = ox.mas %>%
  mutate(year = year(datetime),
         yday = yday(datetime)) %>%
  group_by(lake,year,yday) %>%
  filter(abs(DO_mgL - mean(DO_mgL, na.rm=T)) < 2*sd(DO_mgL, na.rm=T)
         #abs(Temp - mean(Temp, na.rm=T)) < 2*sd(Temp, na.rm=T)
         ) %>%
  ungroup() %>%
  dplyr::select(-c(year,yday))

t.clean = t.mas %>%
  mutate(year = year(datetime),
         yday = yday(datetime)) %>%
  group_by(lake,year,yday,depth_m) %>%
  filter(abs(temp_C - mean(temp_C, na.rm=T)) < 2*sd(temp_C, na.rm=T)) %>%
  ungroup() %>%
  dplyr::select(-c(year,yday))

## first calculate extinction coefficients from light data in PAR
 # define function for calculating light extinction
 ext_fun = function(x){
   if (nrow(x) < 2) {
     z <- NA
   } else {
     m = lm(log(par) ~ depth_m, data = x)
     z1 = m$coef[['(Intercept)']]
     z2 = m$coef[['depth_m']]
   }
   data_frame(datetime = x$datetime[1],
              lake = x$lake[1],
              unique_ob = x$unique_ob[1],
              year = x$year[1], 
              yday = x$yday[1], 
              hour = x$hour[1],
              min = x$min[1],
              par0 = exp(z1), # intercept 
              ext = -z2) # slope
 } 
 # calculate extinction for all daytime profiles
 ext.dat <- i.mas %>%
   mutate(par =  ifelse(par>0,par,{i.mas %>% filter(par >0)}$par %>% min)) %>% # makes 0's into lowest observed light value
   mutate(year = year(datetime),
          yday = yday(datetime),
          # adds one to avoid a 0 hour (1:24 vs 0:23) (from Phillips 2020)
          hour = hour(datetime) + 1, 
          min = minute(datetime),
          sun_on = is.day(datetime, lat = 47.9)) %>% # adds a daylight index
   mutate(unique_ob = paste(lake,year,yday,hour-1,min, sep = '_')) %>% #unique id for each observation
   na.omit() %>%
   split(.$unique_ob) %>% # splits by unique ID, could take a while (~ 5 minutes; not very efficient) 
   map_df(~ ext_fun(.x)) %>% # calculates extinction using above function
   #filter(ext > 0) %>%
   mutate(par0 = ifelse(is.day(datetime, lat = 47.9) == FALSE,0,par0)) %>%
   left_join(i.mas %>%
               filter(depth_m == 0) %>%
               dplyr::select(datetime,lake,par)) %>%
   rename(par_surf = par) %>%
   left_join(i.mas %>%
               filter(depth_m == 2) %>%
               dplyr::select(datetime,lake,par)) %>%
   rename(par_2 = par)
 
ext.sum <- ext.dat %>%
  filter(ext > 0) %>%
  group_by(lake) %>%
  summarize(u_ext=mean(ext),
            sd_ext = sd(ext),
            n = length(ext)) %>%
  mutate(se_ext = sd_ext/sqrt(n))
 
  #write_csv(ext.dat %>% filter(ext>0),file = 'Summerfiles/Masterfiles/extinction_day_master.csv')
 
 ext.dat %>%
   ggplot(.,aes(x=yday,y=par0)) +
   geom_point(alpha = 0.1) +
   facet_wrap(~lake,
              scales = 'free') +
   tidybayes::theme_tidybayes()

ext.dat %>%
  filter(par0 > 1) %>% 
     ggplot(.,aes(x=yday,y=exp(ext))) +
     geom_point(alpha = 0.1) +
     facet_wrap(~lake,
                scales = 'free') +
  xlab('Day of Year') +
  ylab(expression(Extinction~Coefficient~(m^-1))) +
   tidybayes::theme_tidybayes()

ext.dat %>%
  group_by(lake) %>%
  summarise(u_ext = mean(ext),
            sd_ext = sd(ext),
            n = length(ext)) %>%
  left_join(meta.dat %>% 
              dplyr::select(lake,ele_m)) %>%
  ggplot(.,aes(x=ele_m,y=u_ext)) +
  geom_point()

ext.dat %>%
  left_join(i.mas %>%
              filter(depth_m == 2) %>%
              dplyr::select(datetime,lake,par)) %>%
  mutate(par_2 = par) %>%
  dplyr::select(-par) %>%
  ggplot(.,aes(x=log(par0),y=log(par_2))) +
  geom_point(alpha=0.08) +
  geom_abline(slope=1,intercept = 0) +
  facet_wrap(~lake, scales = 'free') +
  tidybayes::theme_tidybayes()
 

## Next get mixing depth for each timestep from temperature data

# check for temperature matching (can swap out for light easily)
# t.mas %>%
#   mutate(year = year(datetime),
#          yday = yday(datetime),
#          # adds one to avoid a 0 hour (1:24 vs 0:23) (from Phillips 2020)
#          hour = hour(datetime), 
#          min = minute(datetime)) %>%
#   mutate(yday_hour_min = yday + hour/24 + min/(24*60)) %>%
#   filter(yday %in% c(223:225)) %>%
#   filter(depth_m %in% c(0,1,2)) %>%
#   ggplot(.,aes(x=yday_hour_min,y=temp_C, color=as.factor(depth_m))) +
#   geom_point(alpha=0.5) +
#   facet_wrap(~lake, scales = 'free') +
#   tidybayes::theme_tidybayes()

z.dat <- t.clean %>%
  mutate(year = year(datetime),
         yday = yday(datetime),
         # adds one to avoid a 0 hour (1:24 vs 0:23) (from Phillips 2020)
         hour = hour(datetime), 
         min = minute(datetime)) %>%
  mutate(yday_hour_min = yday + hour/24 + min/(24*60),
         unique_ob = paste(lake,year,yday,hour,min, sep = '_')) %>%
  group_by(lake,unique_ob) %>%
  mutate(zmix = thermo.depth(wtr=temp_C, depths = depth_m)) %>%
  summarize(zmix = zmix[1]) %>%
  #replaces impossible values with NA
  mutate(zmix = ifelse(zmix>0.5 & zmix <15,zmix, NA)) %>%
  # fills NAs with previous zmix
  fill(zmix)

range(z.dat$zmix, na.rm = T)

## look at individual series of oxygen sat data - can delete 
# ox.clean %>%
#   mutate(year = year(datetime),
#          yday = yday(datetime),
#          # adds one to avoid a 0 hour (1:24 vs 0:23) (from Phillips 2020)
#          hour = hour(datetime),
#          min = minute(datetime)) %>%
#   mutate(yday_hour_min = yday + hour/24 + min/(24*60)) %>%
#   left_join(t.clean %>% filter(depth_m==1) %>%
#               mutate(year = year(datetime),
#                      yday = yday(datetime),
#                      # adds one to avoid a 0 hour (1:24 vs 0:23) (from Phillips 2020)
#                      hour = hour(datetime),
#                      min = minute(datetime)) %>%
#               mutate(yday_hour_min = yday + hour/24 + min/(24*60))) %>%
#   mutate(dosat = o2.at.sat.base(temp_C,altitude = with(meta.dat, ele_m[match(.$lake,lake)])),
#          isday = ifelse(is.day(datetime,lat = 47.5) == TRUE,1,0)) %>%
#    filter(yday %in% c(220:221)) %>%
#   # filter(depth_m %in% c(0,1,2)) %>%
#   ggplot(.,aes(x=yday_hour_min,y=100*(DO_mgL/dosat),color=as.factor(isday))) +
#   geom_point(alpha=0.5) +
#   facet_wrap(~lake, scales = 'free') +
#   tidybayes::theme_tidybayes()


### Wind speed simulation
# rnorm(1,1)
# rnorm(with())
# with(wnd.dist, value[match($ele_m,elevation)])

## merge data files by datetime
dat.prep.full <-  ox.clean %>%
  inner_join(.,t.clean %>% filter(depth_m == 1)) %>%
  dplyr::select(-depth_m) %>%
  inner_join(.,i.mas %>% filter(depth_m == 2)) %>%
  dplyr::select(-depth_m) %>%
  mutate(year = year(datetime),
         yday = yday(datetime),
         # adds one to avoid a 0 hour (1:24 vs 0:23) (from Phillips 2020)
         hour = hour(datetime) + 1, 
         min = minute(datetime),
         #par = lux/42,  # rough approximation of lux conversion (will use lit for real on)
         wspeed = 1) %>%
  rename(par_2 = par,
         temp_1 = temp_C) %>%
  # from Phillips 2020 make an identifier for unique series
  arrange(lake,year,yday,hour,min) %>%
  group_by(lake,year) %>%
  mutate(i = ifelse(is.na(DO_mgL)==T, 1, 0), 
         j = c(1,abs(diff(i)))) %>%
  filter(is.na(DO_mgL)==F, is.na(wspeed)==F) %>%
  mutate(series = cumsum(j)) %>% 
  ungroup() %>%
  # create unique index for each series
  # remove series with fewer than 24 observations
  mutate(unique_series = year + series/10) %>%
  group_by(lake, unique_series) %>%
  mutate(series_length = length(unique_series)) %>%
  filter(series_length > 47) %>% # removes any series with less than a day of observations
  ungroup() %>%
  mutate(yday_hour_min = yday + (hour-1)/24 + min/(24*60),
         unique_ob = as.character(paste(lake,year,yday,hour-1,min,sep = '_'))) %>%
  left_join(ext.dat %>% 
              mutate(unique_ob = as.character(paste(lake,year,yday,hour-1,min,sep = '_'))) %>%
                       dplyr::select(unique_ob,par0,ext,par_surf)) %>%
  # correct for negative ext, add average ext for that lake for all values 0 or less
  mutate(ext = ifelse(ext > 0,ext, with(ext.sum, u_ext[match(.$lake,lake)]))) %>%
  left_join(z.dat) %>%
  # calculate avg light in the mixing area as a funciton of modeled extinction and surface light
  mutate(par_int = (par0 - par0*exp(-abs(ext)*zmix))/(ext*zmix)) %>%
  # recreate series index and make unique index for days
  # create index for observations (for joining later)
  # replace 0 par_int with smallest non-zero value (~0.25)
  # also adds elevation (for wind simulation)
  mutate(unique_series = as.factor(unique_series) %>% as.numeric(),
         index = 1:length(DO_mgL),
         par_int = ifelse(par_int==0, 0.001, par_int),
         elevation = with(meta.dat, ele_m[match(.$lake,lake)]),
         hour = hour - 1) %>% # hour back for wind speed distribution
  # adds wind speed parameters
  merge(.,wnd.dist,
        by = c('elevation','hour')) %>%
  # # simulates wind speed from real data parameters for each hour of the day 
  # then scales it to 10 m above surface (our anemometer was 1.5 meters above ground)
  mutate(wspeed = wind.scale.base(1.5, abs(rnorm(mean,sd)))) %>%
  dplyr::select(-i, -j,-mean,-sd) %>%
  mutate(hour = hour + 1) %>% # switches hour back so no zeros 
  # dplyr::select(-i, -j) %>%
  mutate(
    # calculate from LakeAnayzer to compare later
    dosat = o2.at.sat.base(temp_1,altitude = 1000 + with(meta.dat, ele_m[match(.$lake,lake)])),
    # temperature in Kelvin
    temp_k = temp_1 + 273.15,
    # Schmidt number and associated conversion from CO2 to O2
    # based on Wanninkhoff 1992; Holtgrieve et al 2010; Staehr
    # error in original version; second term should be -120.10*temp, but was +120.10*temp (Lottig et al. 2022)
    sch_o2 = 1800.6 - 120.10*temp_1 + 3.7818*temp_1^2 - 0.047608*temp_1^3,
    sch_conv = (sch_o2/600)^(-0.5),
    # O2 solubility in mL/L
    # based on Weiss 1970
    do_sol = exp(-173.4292 + 249.6339*(100/temp_k) + 143.3483*log(temp_k/100) - 
                   21.8492*(temp_k/100)),
    # convert to mg/L of DO
    # use ideal gas law, solved for the number of moles n = P*V/(R*T)
    # use pressure in kPA corrected for Myvatn's elevation of ~300m = 98 kPA
    # use pressure from elevation conversion
    # Engineering ToolBox, (2003). Atmospheric Pressure vs. Elevation above Sea Level. 
    #[online] Available at: https://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html [Accessed Day Mo. Year].
    # equatiion: p = 101325 (1 - 2.25577 10-5 h)5.25588  give Pascals, convert to kpa  
    kpa = (101325 * (1-2.25577*10^-5 * with(meta.dat, ele_m[match(.$lake,lake)]))^5.25588)/1000,
    # use R = 8.3144598 L*kPa/(K*mol)  
    # use O2 molar mass of 2*32; multiply by 1000 to convert from g to mg 
    #do_eq = 32*(98*do_sol)/(8.3144598*temp_k),
    do_eq = 32*(kpa*do_sol)/(8.3144598*temp_k),
    # calculate equilibrium DO from back-calculation from sonde o2_sat for comparison
    do_eq2 = 100*DO_mgL/dosat,
    # calculate gas exchange coefficient using [wind speed] and lake area in meters squared
    k600 = k.vachon.base(wnd = wspeed,lake.area = 10000 * with(meta.dat, lakearea_hectare[match(.$lake,lake)]))) %>%
    mutate(kgas = k600.2.kGAS.base(k600 = k600,temperature = temp_1, gas = "O2")) %>%  #m/d
    mutate(k = (kgas/48)/zmix) %>% #convert gas to T^-1
    # finally, make an index for each day within each lake
    group_by(lake) %>%
    mutate(unique_day = paste(year, yday) %>% as.factor() %>% as.numeric()) %>% 
    ungroup() %>%
    arrange(lake,datetime)

# with constant wind
# before <- hist(dat.prep.full$k,breaks=50)
# after <- hist(dat.prep.full$k,breaks=50)
# with wind simulation corrected to 10 m
dat.prep.full <- read_csv(file = 'Summerfiles/dat.prep.full.csv') %>%
  arrange(lake,datetime)
write_csv(dat.prep.full, file = 'Summerfiles/dat.prep.full.csv')

## Troubleshooting plots

dat.prep.full %>%
  filter(lake =='UPH') %>% 
  filter(yday %in% c(225:230)) %>% 
  ggplot(.,aes(x=datetime,y=DO_mgL/dosat)) +
  geom_point() +
  geom_line() +
  theme_bw()
  




######
###### Package data 
###### (from Phillips 2020) loop through each lake added by JV
for (i in 1:length(unique(dat.prep.full$lake))) {
  lakename = unique(dat.prep.full$lake)[i]
  dat.prep = dat.prep.full %>% filter(lake == lakename) 
  
  o2_obs = 1000*dat.prep$DO_mgL # DO in mg/m^-3
  o2_eq = 1000*dat.prep$dosat # DO sat in mg/m^-3
  light = dat.prep$par_int # light in PAR/
  temp = dat.prep$temp_1 # temp in degrees C
  wspeed = dat.prep$wspeed # windspeed in m/s at 10 m above surface (corrected with LakeMetabolizer)
  sch_conv = dat.prep$sch_conv # conversion
  k = dat.prep$k # gas exchange 
  map_days = dat.prep$unique_day # unique day index for each lake
  days_per_year = c({dat.prep %>%  # how many unique days per year and lake
      group_by(year) %>%
      summarize(value = length(unique(unique_day)))}$value)
  obs_per_series = c({dat.prep %>% # observations per year and lake
      group_by(unique_series) %>%
      summarize(value = length(unique_series))}$value) 
  obs_per_day = c({dat.prep %>% # observations per day
      group_by(unique_day) %>%
      summarize(value = length(unique_day))}$value) 
  z = dat.prep$zmix  # 3.3 # mixing depth (use LakeMetabolizer function)
  k0 = 2.07
  k1 = 0.215
  k2 = 1.7
  n_obs = length(o2_obs) # number of observations
  n_series = length(obs_per_series) # number of continuous series
  n_days = sum(days_per_year) # days with full data coverage
  n_years = length(days_per_year) # number of years
  
  
  # export as .R file to be read by stan model
  stan_rdump(c("o2_obs","o2_eq","light","temp","wspeed","sch_conv","map_days","obs_per_series","days_per_year",
               "obs_per_day", "z","k0","k1","k2","n_obs","n_series","n_days","n_years","k"),
             file=paste(lakename,"stan_prep.R", sep = '_'))
}

  
