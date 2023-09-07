#fixes columns (within load.data)
coln <- function(df) {
  colnames(df) <- c("sens_add", "date", "time", "windSpeed", 
                    "windDirection", "gustWindSpeed", "airTemp",
                    "xOrient", "yOrient", "null", "northWindSpeed", 
                    "eastWindSpeed")
  colnames(df)
}

# fixes times (within load.data)
# must externally define sensorStartDateTime
timefix <- function(df,ssdt) {
  sensorOriginDate = as.Date(x = "2000-01-01 00:00:00")
  df$date <- as.Date(df$date, format = "%d/%m/%y")
  df$senstime <- parse_date_time(paste(df$date, df$time), 
                                 "ymd_HMS")
  difftime(df$senstime, sensorOriginDate) + 
    ssdt
  
}

#loads data and cleans it up
#requires above packages
load.data <- function(file){
  wnd <- read.csv(file = file, header = T)
  
  colnames(wnd) <- coln(wnd)
  wnd$datetime <- timefix(wnd, ssdt = sensorStartDateTime)
  
  wnd
}

#plots windspeed with another variable
#must externally define location (Hol, Uph, or Sapp)
serplot <- function(df, location, var1 = "windSpeed", var2) {
  par(mar = c(2.2,4,2,4))
  v1 <- which(colnames(df) == var1)
  v2 <- which(colnames(df) == var2)
  yeet <- pretty(c(0,ceiling(range(df[[v1]])[2]))) #defines common axis tics
  yeet2 <- pretty(c(floor(range(df[[v2]])[1]), 
                    ceiling(range(df[[v2]])[2])))
  dateyeet <- pretty(range(df$datetime))
  
  plot(x=df$datetime, y = df[[v1]], col = "blue", ylim = range(yeet),
       xlab = "date and time", type = 'l', lwd = 0.75,
       main = location, xaxt = 'n', yaxt = 'n',
       ylab = "")
  axis(2, col = "blue", at = yeet, labels = yeet)
  mtext(2, text = "wind speed (m/s)", line = 2)
  axis(1, at = dateyeet, labels = dateyeet)
  grid(ny = length(dateyeet))
  
  par(new=T)
  plot(x=df$datetime, y=df[[v2]], col ="red",ylim = range(yeet2),
       axes=F, xlab ="", type = 'l', lwd = 0.5, main = "",
       xaxt = 'n', yaxt = 'n', ylab ="")
  axis(4, col = "red", at = yeet2, labels = yeet2, line = 0)
  mtext(side = 4, text = var2, 
        line = 2.5)
}

#writes just date/time and windspeed to output csv 
#calls external floc
wsfile <- function(wnd) {
  wnds <- dplyr::select(wnd, datetime, windSpeed)
  write.csv(wnds, file = paste("Research/Summer 2019/WindSpeed/",floc, "/MetabFiles/", 
                               sensorStartDate, ".txt", sep = ""))
}

##taken from github - loads minidot data hella
read_minidot <- function(fname, skip = 3, ...)
{
  files <- list.files(fname, full.names=TRUE)
  if(requireNamespace('data.table', quietly=TRUE))
  {
    dat <- data.table::rbindlist(lapply(files, data.table::fread, skip = skip, ...))
  } else {
    dat <- do.call(rbind, lapply(files, read.csv, skip = skip, ...))
  }
  
  # drop battery column
  if(ncol(dat) == 5) 
    dat <- dat[,-2]
  
  colnames(dat) <- c('time_sec', 'temperature', 'DO', 'q')
  dat$timestamp <- as.POSIXct(dat[[1]], origin="1970-01-01")
  #dat$timestamp <- strptime(df$timestamp)
  dat
}

##loads hobo pendant data from a folder with multiple depths and puts them in same dataframe
#make sure input csvs are labelled as depth.csv (single number depth)
load.hobo <- function(pat, do) {
  
  setwd(pat)
  df <- ""
  file.names <- dir(pat, pattern =".csv")
  sort(file.names)
  
  for(i in 1:length(file.names)){
    file <- read.csv(file.names[i],header=T, skip = 1)
    if (ncol(file) == 3) {
      file$irr = NA
      colnames(file) <- c(paste("index", "_", str_remove(file.names[i], ".csv"),sep =""),
                          "datetime",
                          paste("wtr", "_", str_remove(file.names[i], ".csv"),sep =""),
                          paste("irr", "_", str_remove(file.names[i], ".csv"),sep =""))
    }
    else {
      colnames(file) <- c(paste("index", "_", str_remove(file.names[i], ".csv"),sep =""),
                          "datetime",
                          paste("wtr", "_", str_remove(file.names[i], ".csv"),sep =""),
                          paste("irr", "_", str_remove(file.names[i], ".csv"),sep =""))
    }
    
    if (i == 1) {
      df <- cbind(df, file[2:4])
    } else {
      df <- merge(df, file[2:4], by = 'datetime')
    }
  }
  
  df$datetime <- strptime(df$datetime, format = '%m/%d/%y %I:%M:%S %p') # fixes time before final merge
  
  dooby <- do[,c(6,2)]  ## adds temp data from Minidot at 1 m
  names(dooby) <- c("datetime", "wtr_1")
  
  df$datetime <- as.POSIXct(df$datetime)
  
  df <- merge(df, dooby, by = 'datetime')
  
  # time is in GMT format, so we need to fix it to MDT timezone
  # df$datetime <- with_tz(ymd_hms(df$datetime),'America/Denver')
  #
  tdf <- df[,c(1,grep("wtr", names(df)))]
  irrdf <- df[,c(1,grep("irr", names(df)))]
  
  
  write.csv(tdf, file = paste(pat,"/wtr_summer.txt", sep = ""))
  write.csv(irrdf, file = paste(pat,"/irr_summer.txt", sep = ""))
}


## from LakeMetabolizer (now defunct?)
##
o2.at.sat.base <- function(temp, baro, altitude=0, salinity=rep(0,length(temp)), model='garcia-benson'){
  
  # Conversion from mL/L (the usual output of the garcia, weiss, etc. equations)
  # to mg/L per USGS memo 2011.03
  mgL.mlL <- 1.42905
  
  # Correction for air pressure; incorportes effects of altitude & vapor pressure of water
  mmHg.mb <- 0.750061683 # conversion from mm Hg to millibars
  if(missing(baro)){
    mmHg.inHg <- 25.3970886 # conversion from inches Hg to mm Hg
    standard.pressure.sea.level <- 29.92126 # Pb, inches Hg
    standard.temperature.sea.level <- 15 + 273.15 # Tb, 15 C = 288.15 K
    gravitational.acceleration <- 9.80665 # g0, m/s^2
    air.molar.mass <- 0.0289644 # M, molar mass of Earth's air (kg/mol)
    universal.gas.constant <- 8.31447 #8.31432 # R*, N*m/(mol*K)
    
    # estimate pressure by the barometric formula
    baro <- (1/mmHg.mb) * mmHg.inHg * standard.pressure.sea.level * 
      exp( (-gravitational.acceleration * air.molar.mass * altitude) / (universal.gas.constant * standard.temperature.sea.level) )
  }
  # pressure correction per USGS memos 81.11 and 81.15. calculate u by Antoine equation.
  u <- 10 ^ (8.10765 - 1750.286 / (235 + temp)) # u is vapor pressure of water; water temp is used as an approximation for water & air temp at the air-water boundary
  press.corr <- (baro*mmHg.mb - u) / (760 - u) # pressure correction is ratio of current to standard pressure after correcting for vapor pressure of water. 0.750061683 mmHg/mb
  
  # Estimate O2 at saturation in mL/L by several models
  if(tolower(model) == 'garcia'){
    
    Ts <- log((298.15 - temp)/(273.15 + temp))
    
    lnC <- 2.00856 + 3.22400 *Ts + 3.99063*Ts^2 + 4.80299*Ts^3 + 9.78188e-1*Ts^4 + 
      1.71069*Ts^5 - salinity*(6.24097e-3 + 6.93498e-3*Ts + 6.90358e-3*Ts^2 + 4.29155e-3*Ts^3) - 3.1168e-7*salinity^2
    
    o2.sat <- exp(lnC)
    
  } else if(tolower(model) == 'garcia-benson'){
    
    Ts <- log((298.15 - temp)/(273.15 + temp))
    
    lnC <- 2.00907 + 3.22014*Ts + 4.05010*Ts^2 + 4.94457*Ts^3 + -2.56847e-1*Ts^4 + 
      3.88767*Ts^5 - salinity*(6.24523e-3 + 7.37614e-3*Ts + 1.03410e-2*Ts^2 + 8.17083e-3*Ts^3) - 4.88682e-7*salinity^2
    
    o2.sat <- exp(lnC)
    
  } else if(tolower(model) == 'weiss'){
    tempk <- temp + 273.15
    
    lnC <- -173.4292 + 249.6339 * (100 / tempk) + 143.3483 *
      log(tempk / 100) - 21.8492 * (tempk / 100) + 
      salinity * (-0.033096 + 0.014259 * (tempk / 100) - 0.0017000 * (tempk / 100)^2)
    
    o2.sat <- exp(lnC)
    
  } else if(tolower(model) == 'benson'){
    ## TODO: Fix this to include salinity
    if(!all(salinity==0)){
      warning('Benson model does not currently include salinity')
    }
    
    o2.sat <- (-0.00006 * (temp)^3) + (0.00725 * (temp)^2) - (0.39571 * (temp)) + 14.59030
    
    o2.sat <- o2.sat / mgL.mlL # undo the conversion (below) from ml/L to mg/L; Benson model appears to already predict in mg/L
    
  } else {
    stop(paste0('unrecognized model: ', model))
  }
  o2.sat <- o2.sat * mgL.mlL * press.corr
  
  return(o2.sat)
}

k.vachon.base <- function(wnd, lake.area, params=c(2.51,1.48,0.39)){
  U10 <- wnd  #This function uses just the wind speed it is supplied
  k600 <- params[1] + params[2]*U10 + params[3]*U10*log10(lake.area/1000000) # units in cm h-1
  k600 <- k600*24/100 #units in m d-1
  return(k600)
}

getSchmidt	<-	function(temperature, gas){
  
  range.t	<-	c(4,35) # supported temperature range
  
  Schmidt	<-	data.frame(
    "He"=c(368,-16.75,0.374,-0.0036),
    "O2"=c(1568,-86.04,2.142,-0.0216),
    "CO2"=c(1742,-91.24,2.208,-0.0219),
    "CH4"=c(1824,-98.12,2.413,-0.0241),
    "SF6"=c(3255,-217.13,6.837,-0.0861),
    "N2O"=c(2105,-130.08,3.486,-0.0365),
    "Ar"=c(1799,-106.96,2.797,-0.0289),
    "N2"=c(1615,-92.15,2.349,-0.0240)
  )
  
  obsT <- is.finite(temperature) # logical for observed (not NA or NaN [or Inf or -Inf]) -RDB
  
  if (!is.character(gas)){stop(paste('gas must be a character. was given as',gas))}
  if (length(gas)>1){stop("only one gas can be specified for this version")}
  if (!any(names(Schmidt)==gas)){stop(paste(gas,'not found in list of coded gasses'))}
  if (any(temperature[obsT] < range.t[1] | temperature[obsT] > range.t[2])){ # This logical threw an error if any temperature were NA (or NaN, etc.) -RDB
    warning("temperature out of range")
  }
  A	<-	unlist(Schmidt[gas])[1]
  B	<-	unlist(Schmidt[gas])[2]
  C	<-	unlist(Schmidt[gas])[3]
  D	<-	unlist(Schmidt[gas])[4]
  
  Sc = as.numeric(A+B*temperature+C*temperature^2+D*temperature^3)
  
  return(Sc)
}

k600.2.kGAS.base	<-	function(k600,temperature,gas="O2"){
  
  n	<-	0.5
  schmidt	<-	getSchmidt(temperature,gas)
  Sc600	<-	schmidt/600
  
  kGAS	<-	k600*(Sc600^-n)
  return(kGAS)
}
