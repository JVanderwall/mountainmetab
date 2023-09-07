#### Working with Phillips et al. 2020 inhibition model
### modified by J. Vanderwall for elevation lake metabolism project
### this script loops through all lakes and runs and saves model output (and model itself)

rm(list=ls())

# load packages
library(tidyverse)
library(rstan)
library(loo)
library(patchwork)
library(lubridate)
library(tidybayes)
library(gridExtra)
library(shinystan)
library(bayesplot)
setwd('/Users/josephvanderwall/Documents/Research/Data Depository/All data (for R)/')
source("Rfiles/stan_ultility.R")

# windows specific
# memory.limit(size = 100000)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(javascript = FALSE)

# list of all the lakes in this study 
lake.list = c("ASH", "BON", "CLW", "HOL", "LNG", "LOR", "MCD", 
              "MOR", "MSO", "RNY", "SAP", "SMT", "SUM","UPH")
             ## "MSOL","LORL","MORL")

# data masterfile
full.dat <- read_csv(paste0('Summerfiles/dat.prep.full.csv'))
unique(full.dat$lake)

model <- "o2_model_inhibition_JVedit.stan" #Steele 2 param inhibition
model_path <- paste0("Rfiles/",model, sep = '')

sm <- stan_model(file = model_path,verbose = F)

output_path <- paste0("PhillipsModel/output/")

# model setup
chains <- 2
iter <- 1000
warmup <- 400
adapt_delta <- 0.99
max_treedepth <- 15
thin <- 1

i = 4
# for (i in 1:length(lake.list)) {
   rm(data,modelfit,fit_summary,fit_clean,out,out2,out3,p2)

  lakey = lake.list[i]
  
  data <- read_rdump(paste("PhillipsModel/inputs/",lakey,"_stan_prep.R",sep=""))
  
  data$sig_b0 <- 0.5 #pmax smoothing parameter
  data$sig_r <- 0.5  #respiration smoothing parameter
  data$sig_i0 <- 0.2  #light saturation smoothing parameter
  data$sig_obs <- 0.1
  data$days_per_year <- as.array(data$days_per_year)
  data$obs_per_series <- as.array(data$obs_per_series)
  data$o2_freq <- 48
  # data$o2_eq <- data$o2_eq
  # data$o2_obs <- data$o2_obs
  #data$i0 = 6 # makes inhibiton function a straight line
  data$gamma_1 = 1 # fixes gpp temp scaling to one (aka none)
  #data$gamma_2 = 1 # fixes er temp scaling to one (aka none)
  data$temp_ref <- data$temp %>% mean() 
  data$k <- data$k  ## testing windspeed
  modelfit <- sampling(sm, data = data, chains = chains, cores = chains, iter = iter, warmup = min(iter*0.5,warmup),
                      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth), 
                      seed=19990909,thin = thin,save_warmup=FALSE, verbose = T)
  
  fit_summary <- summary(modelfit, probs=c(0.025,0.5,0.975))$summary %>% 
    {as_tibble(.) %>%
        mutate(var = rownames(summary(modelfit)$summary))}
  
   # fit_summary2 <- summary(modelfit, probs=c(0.16,0.5,0.84))$summary %>% 
  #   {as_tibble(.) %>%
  #       mutate(var = rownames(summary(modelfit)$summary))}
  
  fit_clean <- fit_summary %>%
    rename(lower = '2.5%', middle = '50%',upper = '97.5%')  %>%
    mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
           index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           day = ifelse(name %in% c("GPP","ER","NEP","AIR","Flux","GPP_m2","ER_m2","NEP_m2","b0","r","i0","b"), 
                        index, 
                        data$map_days[index])) %>%
    dplyr::select(name, index, day, middle,lower,upper)
  
  # fit_clean2 <- fit_summary2 %>%
  #   rename(lower = '16%', middle = '50%',upper = '84%')  %>%
  #   mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
  #          index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
  #          day = ifelse(name %in% c("GPP","ER","NEP","AIR","Flux","GPP_m2","ER_m2","NEP_m2","b0","r","i0","b"), 
  #                       index, 
  #                       data$map_days[index])) %>%
  #   dplyr::select(name, index, day, middle,lower,upper)
  
  out <- fit_clean %>%
    filter(name %in% c("GPP","ER","NEP")) %>%
    rename(unique_day = day) %>% 
    left_join(.,full.dat %>% filter(lake == lakey) %>%
                dplyr::select(unique_day,yday,year) %>% distinct()) %>% 
    left_join(.,full.dat %>% filter(lake == lakey) %>%
                expand(year,yday,name=c("GPP","ER","NEP"))) %>% 
    mutate(middle = ifelse(name=="ER",-middle,middle),
           lower = ifelse(name=="ER",-lower,lower),
           upper = ifelse(name=="ER",-upper,upper),
           name = factor(name, levels=c("GPP","NEP","ER")),
           lake = lakey)
  out2 <- fit_clean %>%
    filter(name %in% c("GPP_m2","ER_m2","NEP_m2"))%>%
    rename(unique_day = day) %>% 
    left_join(full.dat %>% filter(lake == lakey) %>%
                dplyr::select(unique_day,yday,year) %>% distinct()) %>% 
    full_join(full.dat %>% filter(lake == lakey) %>%
                expand(year,yday,name=c("GPP_m2","ER_m2","NEP_m2"))) %>% 
    mutate(middle = ifelse(name=="ER_m2",-middle,middle),
           lower = ifelse(name=="ER_m2",-lower,lower),
           upper = ifelse(name=="ER_m2",-upper,upper),
           name = factor(name, levels=c("GPP_m2","NEP_m2","ER_m2")),
           lake = lakey)
  
  out3 <- rbind(out,out2)
  
  p2 <- ggplot(data = out2 %>% drop_na(year),aes(yday, middle, color = name))+
    geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = name),
                linetype = 0, alpha = 0.2)+
    geom_line() + geom_point(size = 0.5) +
    ggtitle(lakey) +
    # geom_point(data = out %>% left_join(c14),aes(x=yday,y=(p80/12.011),color="C14")) +
    scale_color_manual(values = c("dodgerblue","black","firebrick")) +
    scale_fill_manual(values = c("dodgerblue","black", "firebrick")) +
    theme_bw() +
    # tidybayes::theme_tidybayes() +
    labs(y=expression(mg~O[2]~m^-2~d^-1)) 
    # facet_wrap(vars(year))
  p2
  
  ## check distributions
ggplot(data = out2 %>% drop_na(year),aes(x=middle,y=lake,
                 fill = name)) +
      stat_slab(color='black',size = 0.7, scale = 0.6) + #xlim(-80,80) +
      # stat_dotsinterval(side = 'bottom',scale = 0.66) +
      scale_fill_manual(breaks = c('GPP_m2','ER_m2','NEP_m2'),
                        values = c('darkgreen','darkred','grey25'),
                        name = '',
                        labels = c('GPP','ER','NEP')) +
    # geom_hline(yintercept = 0, size = 0.3, color = "gray50")+
    # geom_ribbon(aes(ymin = lower, ymax = upper, fill = name),
    #             linetype = 0, alpha = 0.2)+
   # geom_line() + geom_point(size = 0.5) +
    ggtitle(lakey) +
    # geom_point(data = out %>% left_join(c14),aes(x=yday,y=(p80/12.011),color="C14")) +
    # scale_color_manual(values = c("dodgerblue","black","firebrick")) +
    # scale_fill_manual(values = c("dodgerblue","black", "firebrick")) +
    theme_bw() +
    # tidybayes::theme_tidybayes() +
    labs(x=expression(mg~O[2]~m^-2~d^-1)) 
  # facet_wrap(vars(year))
  
  ggsave(plot = p2,filename = paste0(output_path,lakey,'_timeseries_',
                                     month(Sys.Date()),'_',day(Sys.Date()),'.png'),
         width=11,height=8.5,dpi=300)
  
  # export
  write_csv(out3, paste0(output_path,lakey,"_","_daily_full_",
                         month(Sys.Date()),"_",day(Sys.Date()),".csv"))
  write_csv(fit_clean, paste0(output_path,lakey,"_summary_clean_",
                              month(Sys.Date()),"_",day(Sys.Date()),".csv"))
  # save model full output
  saveRDS(modelfit, paste0(output_path,lakey,"_fit_",
                          month(Sys.Date()),"_",day(Sys.Date()),".rds"))
  
# }
  
# functions for troubleshooting
  
  check_n_eff(modelfit)
  check_rhat(modelfit)
  check_div(modelfit)
  check_treedepth(modelfit,max_treedepth)
  check_energy(modelfit)
  
  launch_shinystan(modelfit)
  
 full.dat %>%
    group_by(lake) %>%
    summarise(start = datetime[1])
 
 full.dat %>%
   filter(lake == 'UPH')
 