suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggpubr)
})


#### OBTAINING K WITH NIGHT TIME REGRESSIONS ####
### IT NEEDS A FILE WITH A  TIME, TEMPERATURE, DISSOLVED OXYGEN, and DISSOLVED OXYGEN AT 100% sat
### IMPORTANT, WITH THE HEADINGS "local.time", "temp.water", "DO.obs" , "DO.sat".
### It performs several regressions for each day, moving the window of data for the regression, and saving the best regression
### You need to add several parameters: 
###     site = the file with the data
###     reg_start = the time to start the regressions (as  HH:MM:SS; '18:00:00')
###     reg_size = thelength of the period for the regressions (in hours)
###     num_reg = the number of regressions you want to do
###     timesteps = the timesteps of your data, in minutes
### It should be run as: data.out <- night_fun(site=data, reg_start= '17:00:00', reg_size=2, num_reg = 5, timesteps=10)
### This example would do the night time regression of the periods: 17-19;18-20;19-21;20-22;21-23 and save the best one
### The output consists in a dataframe with the day, slope,  slope.se, intercept, r2, and k600 as (d-1)

night_fun <- function(site, reg_start, reg_size, num_reg , timesteps){
  
  #first step is to make a new file with the daily mean values. This used only to be able to select each day, and take the temperature for the k600 calcs
  daily = site %>% group_by(Daily = format(local.time, "%Y-%m-%d")) %>%
    summarise_all(funs(mean=mean(., na.rm=T)))
  
  #initialising variables
  daily$slope <- NA
  daily$r2    <- NA
  daily$intercept <- NA
  daily$slope.se <- NA
  daily$k600 <- NA
  

  
  #the main loop to go through each day/night
  for(i in 2:length(daily$Daily)-1){
    #print(paste("day", i))
    #extract the day and the day after to isolate the night period
    day1 <- as.Date(daily$local.time_mean[i])
    day2 <- as.Date(daily$local.time_mean[i+1])
    hour1 <- '15:00:00'
    
    #and make a POSIX.ct object to be able to subset the file later
    Date1 <- as.POSIXct( paste(day1 , hour1))
    Date2 <- as.POSIXct( paste(day2 , hour1))
    
    #this is the data for each night
    test<-subset(site, local.time> Date1 & local.time < Date2)
    
    
    #Smooth O2 data and calculate dC/dt
    oxyf1<-stats::filter(test$DO.obs,rep(1/9,9), sides=2)
    oxyf2<-  oxyf1[c(-1,-length(oxyf1))]
    #Convert the values to  mgO2 L-1 day-1. Careful, if your data is not with 10min frequency you need to change this
    test$deltaO2<-c(0,0,0,((oxyf2[-1]-oxyf2[-length(oxyf2)])/timesteps)*1440)
    
    #calculate the oxygen deficit in mgO2/L
    test$O2def <-test$DO.sat-test$DO.obs
    
    #a sequence to do several regressions. It can be changed, now it makes 6 regressions
    twindow <- seq(from=1, to=num_reg*3600, by=3600 )
    
    #new data frame to store the "varying regressions" for each dat
    varyreg <- data.frame(slope=double(length(twindow)),intercept=double(length(twindow)),slope.se=double(length(twindow)), r2=double(length(twindow)), reg_start=double(length(twindow)))
    
    #second loop to do the multiple regressions
    for( j in 1:length(twindow)){
       
      
        #first it defines the period to make the regression, using the twindow sequence to move the window
        reg_per1 <- as.POSIXct(as.POSIXct(paste(day1 , reg_start))+twindow[j])
        reg_per2 <-  reg_per1 + reg_size*3600
        
        #for each window, subset the data to make the regression
        nightr_d <- subset(test, local.time> reg_per1 & local.time < reg_per2)
        
        #make a linear regression and store all variables of interest
        fit1 <- lm( nightr_d$deltaO2~nightr_d$O2def)
        varyreg$slope[j] <- summary(fit1)$coefficients[2, 1]
        varyreg$r2[j]    <- summary(fit1)$r.squared
        varyreg$slope.se[j] <- summary(fit1)$coefficients[2, 2]
        varyreg$intercept[j] <-  summary(fit1)$coefficients[1, 1]
        varyreg$reg_start[j]    <- as.POSIXct(reg_per1)
    }
    
    #now we select only the values that have a positive slope
    varyreg2 <- subset(varyreg, slope>0)
    #in case of no values that are positive, in order to not get any error, I just use the full dataset, but that day will be useless
    if(nrow(varyreg2)==0){ varyreg2 =varyreg }
    
    #find which of the regressions was the best
    bestie <- which.max(varyreg2$r2)
    
    #and save the relevant variables
    daily$slope[i] <- varyreg2$slope[bestie]
    daily$intercept[i] <- varyreg2$intercept[bestie]
    daily$slope.se[i] <- varyreg2$slope.se[bestie]
    daily$r2[i]    <- varyreg2$r2[bestie]
    
    #finally convert the slope to k600 with the temperature correction
    daily$k600[i]<- ((600/(1800.6-(120.1*daily$temp.water_mean[i])+(3.7818*daily$temp.water_mean[i]^2)-(0.047608*daily$temp.water_mean[i]^3)))^-0.5)*varyreg2$slope[bestie]
        
    ##### THIS IS A BONUS, TO MAKE A COMPOSITE PLOT OF [O2], dO2, LIGHT AND THE REGRESSION ####
    ####It is disconnected by default, if you want to use it make sure you have light data also

       #recover the best regression and make it again
       reg_per1 <- as.POSIXct(as.POSIXct(paste(day1 , reg_start))+twindow[bestie])
       reg_per2 <-  reg_per1 + reg_size*3600
       nightr_d <- subset(test, local.time> reg_per1 & local.time < reg_per2)

       #make a plot of the regression
       night_regfit<- ggplot(data= nightr_d, aes(x=O2def, y=deltaO2))+
         geom_point()+
         geom_smooth(method=lm)+
         geom_text(x = quantile(nightr_d$O2def, 0.5, na.rm=T), y = quantile(nightr_d$deltaO2, 0.95, na.rm=T), label =lm_eqn(nightr_d$O2def,nightr_d$deltaO2), parse = TRUE)+
         labs(x= expression(O[2]~deficit~(mg~O[2]~L^-1) ),y=expression(Delta~O[2]~(mg~O[2]~L^-1~d^-1)))+
         theme(axis.text = element_text(size=14), axis.title = element_text(size=14))

       # a plot for the O2 deficit
       o2def <-ggplot(data=test)+
         geom_line(aes(x=as.POSIXct(local.time), y=DO.sat-DO.obs), size=2, color="blue3")+
         scale_x_datetime(date_breaks = "6 hour", labels=date_format("%H-%d"))+
         geom_vline(xintercept = as.numeric(reg_per1), linetype=4)+
         geom_vline(xintercept = as.numeric(reg_per2), linetype=4)+
         labs(x="hour - day" ,y=expression(O[2]~deficit~(mg~O[2]~L^-1) ))+
         theme(axis.text = element_text(size=14), axis.title = element_text(size=14))

        #and the delta O2
       o2_delta <- ggplot(data=test)+
         geom_line(aes(x=as.POSIXct(local.time), y=deltaO2), size=2, color="darkred")+
         scale_x_datetime(date_breaks = "6 hour", labels=date_format("%H-%d"))+
         geom_vline(xintercept = as.numeric(reg_per1), linetype=4)+
         geom_vline(xintercept = as.numeric(reg_per2), linetype=4)+
         labs(x="hour - day" ,y=expression(Delta~O[2]~(mg~O[2]~L^-1~d^-1)))+
         theme(axis.text = element_text(size=14), axis.title = element_text(size=14)) 
       #the light
       lightp <- ggplot(data=test)+
         geom_line(aes(x=as.POSIXct(local.time), y=light), size=2, color="goldenrod2")+
         scale_x_datetime(date_breaks = "6 hour", labels=date_format("%H-%d"))+
         labs(x="hour - day" ,y="light")+
         theme(axis.text = element_text(size=14), axis.title = element_text(size=14)) 

       #make a pannel with them
       pan1<-  ggarrange(o2def, o2_delta, lightp, night_regfit , ncol = 2, nrow = 2)

       name <- deparse(substitute(site))

       #and export it with the name of the site and the date as title
       pan1 %>%
         ggexport(filename = paste(name,'_', day1,'.png', sep=""),res = 170, width = 1000, height = 1000)
       ####
          
    }
  
  
  #and return the daily file
  daily
}

#this is a custum function to plot the equation of the regression in the plot, just a detail
lm_eqn <- function(x,y){
  m <- lm(y ~ (x));
  eq <- substitute(italic(y) == a +  b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
} 
