
## Custom EM Algorithm function:
custom_emalgo <- function(
  Nc, # number of cups till refill 
  n,  # number of sampled data points 
  n1, # number of refills observed
  n0 = n - n1, # number of non-refills observed
  pi,      # probability of me sampling a random cup duing the day
  m = 100, # larger number for numerical approximation
  p_new = 0.5,  # initial value
  eps_criteria = 0.00000001 # convergence criteria
){
  p_old <- p_new + 1
  z     <-  0:m 
  iter <- 0
  while( abs(p_new - p_old) > eps_criteria){
    # update p_old
    p_old <- p_new
    pstar <- p_old*(1 - pi)
    
    # get conditional of (Z = z | X = x):
    marg_z <- (1 - pstar)*(pstar^z)
    Cz_x0 <- ((Nc - 1) / (Nc + z)) * marg_z  # x = 0
    Cz_x1 <- ((1 + z) / (Nc + z)) * marg_z   # x = 1
    Cz_x0 <- Cz_x0 / sum(Cz_x0)
    Cz_x1 <- Cz_x1 / sum(Cz_x1)
    
    # get moment for (Z = z | X = x):
    m0 <- sum(z*Cz_x0)         # x = 0
    m1 <- sum(z*Cz_x1)         # x = 1
    
    # calculate update for p_new:
    p_new <- ((m0*n0 + m1*n1) / (n +  m0*n0 + m1*n1) ) * (1 / (1 - pi))
    iter <- iter + 1
  }
  return(list(percent = 100*p_new, iterations = iter))
}


## load libraries and data
library(ggplot2)
library(dplyr)
library(xlsx)
mydf <- read.xlsx(file = "C:/Users/Chris/Desktop/coffeeStudy/myData.xlsx", sheetIndex = 1)
## create time of day variable
tod <- lapply(strsplit(as.character(mydf$Time), " "), function(x) unlist(strsplit(x[2], ":")))
tod <- unlist(lapply(tod, function(x) as.numeric(x[1]) + as.numeric(x[2])/60))
mydf$tod <- tod
## refill indicator
refill <- mydf$Refilled
## values for later
n  <- nrow(mydf)
n1 <- sum(refill)
n0 <- n - n1
multiplot_colors <- c("#6C3483", "#FFC300")
singleplot_color <- "#0E6655"

## Overall
cat(paste0("Trials:   ", nrow(mydf),
           "\nRefills:   ",sum(mydf$Refilled), " (", round(100*mean(mydf$Refilled),digits=2),"%)",
           "\nCup Size:   8.0oz\nPot Size:  74.4oz\nExpected:  10.75%\n"))

## interesting time facts:
timediff <- split(mydf$tod, mydf$Date)
timediff <- timediff[unlist(lapply(timediff, length)) > 1]
timediff <- unlist(lapply(timediff, diff))
# average time between two cups on the same date is ~3hrs:
summary(timediff)
# minimun time difference is 48 minutes:
summary(60*timediff)
# times when all 3 pots were empty:
mydf[mydf$Strikeout == 1, "Time"]

## overall em estimates:
em_overall <- data.frame(approx_cup_size = c("8oz","9oz","10oz","12oz","16oz","20oz"),
                      cups_in_pot  = c(9,8,7,6,5,4),
                      delinquency_pct_pi05 = rep(NA, 6),
                      delinquency_pct_pi10 = rep(NA, 6))
for(cc in 1:6){
  em_overall$delinquency_pct_pi05[cc] <- round(custom_emalgo(Nc = em_overall$cups_in_pot[cc], 
                                  n = n,  n1 = n1, pi = 0.05)$percent, digits = 5)
  em_overall$delinquency_pct_pi10[cc] <- round(custom_emalgo(Nc = em_overall$cups_in_pot[cc], 
                                  n = n,  n1 = n1, pi = 0.10)$percent, digits = 5)
}
em_overall
## note increasing pi will also increase delinquency_pct

## Day of week
mydf$DayOfWeek <- factor(mydf$DayOfWeek, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
week_sum <- mydf %>% group_by(DayOfWeek) %>%
  summarise(Trials = length(Refilled), 
            Refills = sum(Refilled),
            RefillPct = round(100*mean(Refilled),digits=2))
week_sum$delinquency_pct <- NA
for(cc in 1:5){
  week_sum$delinquency_pct[cc] <- custom_emalgo(Nc = 8, 
                                   n = week_sum$Trials[cc],  
                                   n1 = week_sum$Refills[cc], 
                                   pi = 0.05)$percent
}
week_sum

plot_data <- data.frame(Day   = rep(week_sum$DayOfWeek,2),
                        Value = c(week_sum$RefillPct, week_sum$delinquency_pct),
                        Bars  = c(rep(c("My Refill %", "Coworker Delinquency %"), each = 5)))
p <- ggplot(plot_data, aes(x=Day, y=Value, fill = Bars ))+
  geom_bar(stat="identity",position="dodge") +
  scale_fill_manual(values = multiplot_colors) +
  ggtitle("Refill Risk by Day of Week") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Percent") + ylim(0,65)
p

## daily frequency (2-dim)
two_dim_sum = function(x, cases){
  lc <- length(cases)
  out_v <- rep(NA, lc)
  for(cc in 1:lc){
    out_v[cc] <- sum(x == cases[cc])
  }
  return(out_v)
}
TAB2  <- table(mydf$Date, mydf$DayOfWeek)
cases <- sort(unique(TAB2))
cases <- cases[cases != 0] 
out_two_sum <- as.data.frame(matrix(NA, ncol = ncol(TAB2), nrow = length(cases)))
names(out_two_sum) <- attributes(TAB2)$dimnames[[2]]
for(cc in 1:ncol(out_two_sum)){
  out_two_sum[,cc]  <- two_dim_sum(x = TAB2[,cc], cases = cases)
}
out_two_sum

## Kernel density estimate
kde <- density(tod, bw = "sj", kernel = "epanechnikov", adjust = 1)
plot_df <- data.frame(x = kde$x, y = kde$y)
plot_df <- plot_df[plot_df$x >= 7 & plot_df$x <= 16, ]
p <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_line(size = 2, colour = singleplot_color) +
  ylab("Density") +
  xlab("Time of Day") +
  scale_x_continuous(breaks = 7:16,  limits=c(7, 16)) +
  ggtitle("My Arrival Distribution") + 
  theme(plot.title = element_text(hjust = 0.5))
p

## Nadaraya-Watson kernel regression
reg_model <- ksmooth(tod, refill, "normal", bandwidth = 1.5)
## create plot over time
plot_df <- data.frame(Hour = reg_model$x, 
                      Refill = 100*reg_model$y)
plot_df <- plot_df[plot_df$Hour > 7 & plot_df$Hour < 15.5,]
## em-algo, using a hypothetical sample of 100 to create n, n1
plot_df$Delinquency <- NA
n1s <- plot_df$Refill
plot_df$Delinquency <- sapply(n1s, 
                              function(x){
                                custom_emalgo(Nc = 8, 
                                              n = 100,  
                                              n1 = x, 
                                              pi = 0.05)$percent
                              })
plot_df <- data.frame(Hour = rep(plot_df$Hour, 2),
                 Percent = c(plot_df$Refill, plot_df$Delinquency),
                 Lines   = c(rep("My Refill %", nrow(plot_df)),
                             rep("Coworker Delinquency %", nrow(plot_df))))
p <- ggplot(plot_df, aes(x = Hour, y = Percent, colour = Lines)) +
  geom_line(size = 2) +
  ylab("Percent") +
  ylim(0,100) +
  scale_colour_manual(values = multiplot_colors) +
  ggtitle("Refill Risk by Time of Day") + 
  theme(plot.title = element_text(hjust = 0.5))
p

## create a plot of [-10min, +10min] from the closest half hour
half_hrs <- (0:47)/2
from_half_hour <- sapply(tod, 
       function(x){
          tmp <- 60*(x - half_hrs)
          return(tmp[which.min(abs(tmp))])
       }
      )
halfhr_model <- ksmooth(from_half_hour, refill, "normal", bandwidth = 5)
## em-algo, using a hypothetical sample of 100 to create n0, n1
n1s <- 100*halfhr_model$y
halfhr_Delinquency  <- sapply(n1s, 
                              function(x){
                                custom_emalgo(Nc = 8, 
                                              n = 100,  
                                              n1 = x, 
                                              pi = 0.05)$percent
                              })
plot_df <- data.frame(Mins = rep(halfhr_model$x,2), 
                      Percent = c(100*halfhr_model$y, halfhr_Delinquency),
                      Lines   =  c(rep("My Refill %", length(halfhr_Delinquency)),
                                   rep("Coworker Delinquency %", length(halfhr_Delinquency))))
time_lag <- 7
plot_df <- plot_df[abs(plot_df$Mins) <= time_lag, ]
p <- ggplot(plot_df, aes(x = Mins, y = Percent, colour = Lines)) +
  geom_line(size = 2) +
  ylab("Percent") +
  xlab("Minutes from Half Hour") +
  scale_colour_manual(values = multiplot_colors) +
  ylim(0,100) +
  scale_x_continuous(breaks = 2*((-time_lag/2):(time_lag/2)), limits = c(-time_lag, time_lag)) +
  ggtitle("Refill Risk vs Time from Half Hour") + 
  theme(plot.title = element_text(hjust = 0.5))
p

## some bootstraped simulated CIs:
B = 100000
sim_result <- rep(NA, B)
set.seed(123)
sample_n1s <- rbinom(n = B, size = n, prob = n1/n)
tab_n1s    <- table(sample_n1s)
tab_values <- as.numeric(names(tab_n1s)) 
tab_pcts   <- rep(NA, length(tab_values))
for(bb in 1:length(tab_values)){
  tab_pcts[bb] <- custom_emalgo(Nc = 4, n = n,  n1 = tab_values[bb], pi = 0.05)$percent
}
boot_percentiles <- cumsum(tab_n1s)/B
# 2.5 percentile:
tab_pcts[tail(which(boot_percentiles <= 0.025),1)]
# 97.5 percentile:
tab_pcts[which(boot_percentiles >= 0.975)[1]]
# 95% Bootstrap CIs: 
# NC = 9: (36, 76)
# Nc = 8: (25, 71)
# Nc = 7: (10, 65)
# Nc = 6: (0, 56)
# Nc = 5: (0, 41)
# Nc = 4: (0, 14)


