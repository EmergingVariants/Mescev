###############################################################################
#Simulate a renwal process (infection) with fixed reproduction number and serial
#interval, with multiple introductions at different times
###############################################################################

library(tidyverse)
library(rstanarm)

#rm(list=ls()) # WARNING: delete all objects in your current R session


args = commandArgs(trailingOnly=TRUE)

#Arguments
  if(length(args)<6){
    print("Need at least 6 arguments: 1 - file for parameters, 2 - file for introductions, 3- output_folder 4-first line to read 5-last line to read 
          6 output_file_prefix 7[optional] - generation time file. \n")
    exit(1)
  }
Parameters_file <- args[1]
Introduction_file <- args[2]
output_path <- args[3]
first_line <- as.numeric(args[4])
last_line <-as.numeric(args[5])
output_prefix<- args[6]
# 
# 
# Introduction_file='/home/fra/Projects/EmergingVariants/Omicron/Results/Introductions/BA5-Omicron/BA5-Omicron'
# Parameters_file = '/home/fra/Projects/EmergingVariants/Omicron/Params/Renewal_aug_params_Lit/Parameters_BA5-Omicron.csv'
# output_path = '/home/fra/Projects/EmergingVariants/Omicron/Results/Projections/Renewal_base/BA5-Omicron'
# first_line = 1
# last_line = 0
# output_prefix = 'BA5-Omicron_'


params<-read.table(Parameters_file, header=1,sep=',')
params <- na.omit(params) #remove nan
if(last_line==0){
  last_line=nrow(params)
}

for (line in seq(first_line,last_line)){
  R <- params$Rt[line]
  Ve <- params$Ve[line]     #Vaccine efficacy
  Ne <- params$Ne[line]     #Natural immunity efficacy

  source <- params$source_country[line]
  target <- params$target_country[line]
  multiplier_in <- params$undercounting_multiplier_source[line]
  multiplier_out <- params$undercounting_multiplier_target[line]

  if (multiplier_out%%1 == 0){
    #this means that multiplier out gets read as integer but I need a float
    multiplier_out<-paste(multiplier_out,".0",sep='')
  }
  if (multiplier_in%%1 == 0){
    #this means that multiplier out gets read as integer but I need a float
    multiplier_in<-paste(multiplier_in,".0",sep='')
  }


  date_start <- params$date_start[line]
  date_end <- params$date_end[line]
  if (as.Date(date_end) < as.Date(date_start)+21){
    date_end <- as.Date(date_end)+21
  }

  #date_end <- "2020-12-30"
  N <- as.numeric(as.Date(date_end) - as.Date(date_start))

  variant <- params$variant[line]
  introductions <- paste(Introduction_file,source,target,date_start, params$date_end[line],multiplier_in,multiplier_out, sep = "_")
  skip_to_next <- FALSE

  introductions <- paste(introductions,".csv",sep="")

  # Note that print(b) fails since b doesn't exist

  tryCatch(read.csv(introductions,comment.char = "#",header = 1), error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

  introductions<-read.csv(introductions,comment.char = "#",header = 1)

  I_0 <- introductions$cumulative_incidence_SRC[seq(1,N)]
  recovered_country <- introductions$cumulative_incidence_TRG[seq(1,N)]
  vaxxed_country <- introductions$cumulative_incidence_vax[seq(1,N)]
  #
  # I_0 <- c(1:N)*0
  # recovered_country <- c(1:N)*0
  # vaxxed_country <- c(1:N)*0
  #

  # for(day in seq(1,N)){
  #   if(day %in% introductions$X){
  #     index = match(day,introductions$X)
  #     I_0[day] = introductions$cumulative_incidence[index]
  #     recovered_country[day]=introductions$cumulative_incidence_TRG[index]
  #     vaxxed_country[day]=introductions$cumulative_incidence_vax[index]
  #   }
  # }
  #We assume to seed the epidemic at day 1 to N with n infected cases
  #Does need not to be integer
  #In this example I use I~I_0*exp((R_0-1)/T t), where
  #R_0~1.3
  #I_0 <- 5*exp((1.3-1)/5*seq(from=1,to=N))
  #We assume Rt is fixed at a constant value, but can vary on a daily basis
  #To mimic daily variations of Rt
  #R <- 1.4

  # The generation time distribution: the nth element is the fraction of
  # transmissions that occur n days after infection.

  if(length(args)>6){
    si<-scan(args[6], sep = ',')
  } else {
  # This is a discretised gamma distribution with mean 5.5 and std dev 2.14, as
  # befits covid (Ferretti, Ledda et al 2020)
  si <- c(0.0028993756, 0.0424730927, 0.1240492512, 0.1872169056, 0.1967815246,
          0.1645307574, 0.1174711149, 0.0747170618, 0.0435081787, 0.0236309301,
          0.0121317746, 0.0059451849, 0.0028018312, 0.0012772359, 0.0005657808)
  }
  max_si <- length(si)

  # R_tilde(t), which is the reproduction number. It has a weekly random change:
  # random draw from a normal distribution (that gets subsequently inverse logit
  # transformed, like all covariates predicting R, in order to keep R in the range
  # of interest). Specify the standard deviation
  # of that normal here. (Zero means no walk.)

  if (!is.null(params$Rt_stddev)) {
    R_daily_step_stddev <- 0
  } else {
    R_daily_step_stddev <- 0
  }

  R_t <- rep(R,N)
  if (R_daily_step_stddev != 0){
    for (day in seq(2,N)) {
      R_t[day] <-invlogit(rnorm(n = 1, mean = 0, sd = R_daily_step_stddev))-0.5+R_t[day-1]
    }
  }
#This simulates one renewal process
Renewal_process <- function(N,I_0,R_t){


  group_name <- 'country'
  pops <- data.frame(group = group_name)


  # A data frame with one row per day
  df <- tibble(group = group_name,
               day = seq(from = 1, to =  N),
               delta_I = I_0,
               cum_I = NA_real_)


  df$cum_I <- cumsum(df$delta_I)

  for (day in seq(from = 2, to =  N)) {
    # The time window in the recent past for which the incidence then contributes
    # to incidence now:
    infectiousness_range <- seq(from = max(day - max_si, 1), to = day - 1)
    # The incidence for those days:
    contributing_incidences <- df$delta_I[infectiousness_range]
    # The weights for their contribution to incidence now:
    weights <- si[day - infectiousness_range]
    # Put together for incidence now:
    incidence_today <- R_t[[day]] * sum(contributing_incidences * weights)
    #As of now The number of cases averted is just Ve*Vu, natural immunity efficacy
    incidence_today = incidence_today*(1-Ne*recovered_country[[day]])*(1-Ve*vaxxed_country[[day]])
    df$delta_I[[day]] <-  df$delta_I[[day]]+incidence_today
    df$cum_I[[day]] <- df$cum_I[[day - 1]] + incidence_today
  }
  return(df)
}

  contribution_day <- Renewal_process(N,I_0,R_t)



  # Infected=rep(0,N)
  # delta_I = rep(0,N)
  # for (day in seq(from=1,to=N-1)){
  #   contribution_day <- Renewal_process(N-day+1,I_0[day],R_t[day:N])
  #   Infected[day:N] <- Infected[day:N] + contribution_day$I
  #   delta_I[day:N] <- delta_I[day:N] + contribution_day$delta_I
  # }
  #Last day arrivals do not infect anybody
  #Infected[N] <- Infected[N]+I_0[N]

  group_name <- target
  epidemic_day_1 <- date_start
  epi <- tibble(group = group_name,
         date = seq(from = as.Date(epidemic_day_1),
                    by = 1, length.out = N ),
         cum_I = contribution_day$cum_I,
         Incidence=contribution_day$delta_I,
         imported = I_0,
         Rt = R_t[seq(N)])

  output_name = paste(output_path,"/",output_prefix,sep="")
  output_name = paste(output_name,source,target,date_start, date_end,variant,"undercounting_SRC",multiplier_in, "undercounting_TRG", multiplier_out, "R",R, "ne", Ne, "ve",Ve, sep = "_")
  output_name = paste(output_name,'.csv',sep="")
  print(output_name)
  write.csv(epi, output_name, row.names=F)

  #
  # '''
  #   output_name = paste(source,target,date_start, date_end,variant,"undercounting",multiplier,"R",R, sep = "_")
  #   output_name = paste(output_name,'.pdf',sep="")
  #
  #   ggplot() +
  #     geom_point(data = epi, aes(x = date, y = I),size=2) +
  #     labs(x = 'time (days)', y = "Cumulative Incidence", size=14)+
  #     theme_bw()+
  #     theme(text=element_text(size=12), axis.line = element_line(colour = 'black', size = 2, linetype = "solid", ),  panel.border = element_blank(),
  #           panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #   ggsave(output_name,path=paste(output_path, '/Images/',sep=""), device=cairo_pdf,width=5, height=5, units='in', dpi=1200)
  #
  #   ggplot(epi) +
  #     geom_point(aes(x=date,y=Incidence,colour="blue")) +
  #     geom_point(aes(x=date,y = imported,colour="orange")) +
  #     labs(x = 'time (days)', y = "Incidence")+
  #     theme_bw()+
  #     scale_colour_manual("",breaks=c("blue","orange"), values=c("#0C7EF5","#F5A318"),labels=c("Total", "Imported"), name="")+
  #     theme(text=element_text(size=16), axis.text.x = element_text(angle = 90), axis.line = element_line(colour = 'black',  size = 1, linetype = "solid", ),  panel.border = element_blank(),
  #           panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #
  #   ggsave(paste( output_prefix,'_epidemic_line_', line,'_Incidence.pdf',sep=""), path=paste(output_path, '/Images/',sep=""), device=cairo_pdf,width=5, height=5, units='in', dpi=1200)

}



#This simulates one renewal process
# Renewal_process <- function(N,I_0,R_t){
#   
#   
#   group_name <- 'single_process'
#   pops <- data.frame(group = group_name)
#   
#   
#   # A data frame with one row per day
#   df <- tibble(group = group_name,
#                day = seq(from = 1, to =  N),
#                delta_I = c(I_0,rep(0,N-1)),
#                I = NA_real_)
#   
#   
#   df$I[[1]]=I_0
#   df$I <- cumsum(df$delta_I)
#   
#   for (day in seq(from = 2, to =  N)) {
#     # The time window in the recent past for which the incidence then contributes
#     # to incidence now:
#     infectiousness_range <- seq(from = max(day - max_si, 1), to = day - 1)
#     # The incidence for those days:
#     contributing_incidences <- df$delta_I[infectiousness_range]
#     # The weights for their contribution to incidence now:
#     weights <- si[day - infectiousness_range]
#     # Put together for incidence now:
#     incidence_today <- R_t[[day]] * sum(contributing_incidences * weights)
#     #As of now The number of cases averted is just Ve*Vu
#     incidence_today = incidence_today*(1-Ve*Vu)
#     df$delta_I[[day]] <- incidence_today
#     df$I[[day]] <- df$I[[day - 1]] + incidence_today
#   }
#   return(df)
# }
