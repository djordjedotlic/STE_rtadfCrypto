if(!require("RCurl")) install.packages("RCurl"); library("RCurl")
if(!require("jsonlite")) install.packages("jsonlite"); library("jsonlite")
if(!require("rtadfr")) install.packages("rtadfr"); library("rtadfr")
if(!require("dplyr")) install.packages("dplyr"); library("dplyr")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
source("dataStamp.R")


# Set up data----

set.seed(1203)  # for replicability
data <- fromJSON(getURL(paste0('http://thecrix.de/data/crix.json')))
#########
data$date <- as.Date(data$date)
data <- data[1:which(data$date == "2018-09-11"),]
#########

T   <- nrow(data)  # Sample size
r0  <- round(T*(0.01 + 1.8 / sqrt(T)))  # Minimal window size

#test <- rtadf(data[,2], r0 = r0, test = "adf")  # estimate test statistic and date-stamping sequence
#cvs  <- rtadfSimPar(T, nrep = 1000, r0 = r0, test = "adf")  # simulate critical values and date-stamping threshold

#testDf <- list("test statistic" = test$testStat, "critical values" = cvs)  # test results

#print(testDf)  

#if( testDf$`test statistic`[1] >= testDf$`critical values`[[1]]){
#  print("We can reject the null hypothesis and conclude that the time series express traits of a rational bubble")
#} else{ print("We cannot reject the null hypothesis")}



#Rolling Window Right tail  ADF Full sample----

windowSize <- r0  
    
test_rwADF_full <- list()
cvs_rwADF_full <- list()


for(rwadf in 0:(nrow(data) - windowSize)){
      
      windowStart <- rwadf + 1 
      windowEnd   <- windowStart + windowSize - 1
        test_rwADF_full[rwadf + 1] <- rtadf(data[windowStart:windowEnd,2], r0 = windowSize, test = "adf")  # estimate test statistic and date-stamping sequence
        #rwADF_results_test <- c(rwADF_results_test, test_rwADF)
        cvs_rwADF_full[[rwadf + 1]]  <- rtadfSimPar(windowSize, nrep = 1000, r0 = windowSize, test = "adf")
        #rwADF_results_cvs <- c(rwADF_results_cvs, rez_cvs)
        }

rwADF_results_full <- data.frame(row.names = (1:length(cvs_rwADF_full)))
rwADF_results_full$test <- as.numeric(test_rwADF_full)
for(n in 1:nrow(rwADF_results_full)){
rwADF_results_full$cvs_90[n] <- as.numeric(cvs_rwADF_full[[n]][[1]])
rwADF_results_full$cvs_95[n] <- as.numeric(cvs_rwADF_full[[n]][[2]])
rwADF_results_full$cvs_99[n] <- as.numeric(cvs_rwADF_full[[n]][[3]])
}
#row.names(rwADF_results) <- windowSize:nrow(data)
rwADF_results_full$date <- as.Date(data$date[windowSize:nrow(data)])

save(file = "./results/rwADF_results_full.RData", rwADF_results_full)

stamps_rwADF_full <- date.stamp(stat = rwADF_results_full$test, cv = rwADF_results_full$cvs_95, n = nrow(data))

X <-  data.frame(data[windowSize:nrow(data),], "test" = rwADF_results_full[,1], "cvs" = rwADF_results_full[,3])
X$date <- as.Date(X$date)

p_rwadf_full <- ggplot(data = X, aes(x = date)) +
  geom_line(aes(y = test, colour = "The rolling window ADF statistic sequence (left axis)")) +
  geom_line(aes(y = cvs, colour = "The 95% critical value sequence (left axis)")) +
  geom_line(aes(y = price/1000, colour = "CriX (right axis)")) +
  scale_y_continuous(limits = c(-5, 60), 
                   sec.axis = sec_axis(~.*1000, name = NULL)) +
  #scale_x_date(breaks = pt, labels = c("2014", "2015", "2016", "2017", "2018") ) +
  labs(y = NULL, x = "Date") + 
  ggtitle("Rolling-window ADF- Full sample")

p_rwadf_full <- p_rwadf_full + theme(
  panel.background =  element_blank(),
  panel.border = element_rect(linetype = 1, colour = "black", fill = NA),
  panel.grid.major = element_line(linetype = 2, color = "grey90"),
  legend.title=element_blank(),
  legend.position = c(0.2, 0.9),
  legend.background = element_rect(fill = "transparent", colour = "transparent"),
  legend.key = element_rect(fill = "transparent", colour = "transparent"),
  plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
)

for(i in 1:nrow(stamps_rwADF_full)){
  p_rwadf_full <- p_rwadf_full + annotate("rect", xmin = X$date[stamps_rwADF_full[i,1]], xmax = X$date[stamps_rwADF_full[i,2]], 
                    ymin= -5, ymax = 60,alpha=.35, fill="yellow")
}
pdf("./graphs/p_rwadf_full.pdf")
p_rwadf_full
dev.off() 


#Rolling Window Right tail  ADF Reduced sample----
data_reduced <- data[365:nrow(data),]
windowSize_reduced <- round(nrow(data_reduced)*(0.01 + 1.8 / sqrt(nrow(data_reduced))))

test_rwADF_reduced <- list()
cvs_rwADF_reduced <- list()
remove(rwadf)

for(rwadf in 0:(nrow(data_reduced) - windowSize_reduced)){
  
  windowStart <- rwadf + 1 
  windowEnd   <- windowStart + windowSize_reduced - 1
  test_rwADF_reduced[rwadf + 1] <- rtadf(data_reduced[windowStart:windowEnd,2], r0 = windowSize_reduced, test = "adf")  # estimate test statistic and date-stamping sequence
  #rwADF_results_test <- c(rwADF_results_test, test_rwADF)
  cvs_rwADF_reduced[[rwadf + 1]]  <- rtadfSimPar(windowSize_reduced, nrep = 1000, r0 = windowSize_reduced, test = "adf")
  #rwADF_results_cvs <- c(rwADF_results_cvs, rez_cvs)
}

rwADF_results_reduced <- data.frame(row.names = (1:length(cvs_rwADF_reduced)))
rwADF_results_reduced$test <- as.numeric(test_rwADF_reduced)
for(n in 1:nrow(rwADF_results_reduced)){
  rwADF_results_reduced$cvs_90[n] <- as.numeric(cvs_rwADF_reduced[[n]][[1]])
  rwADF_results_reduced$cvs_95[n] <- as.numeric(cvs_rwADF_reduced[[n]][[2]])
  rwADF_results_reduced$cvs_99[n] <- as.numeric(cvs_rwADF_reduced[[n]][[3]])
}
rwADF_results_reduced$date <- as.Date(data_reduced$date[windowSize_reduced:nrow(data_reduced)])

save(file = "./results/rwADF_results_reduced.RData", rwADF_results_reduced)

stamps_rwADF_reduced <- date.stamp(stat = rwADF_results_reduced$test, cv = rwADF_results_reduced$cvs_95, n = length(rwADF_results_reduced$test))


X <-  data.frame(data_reduced[windowSize_reduced:nrow(data_reduced),], "test" =rwADF_results_reduced$test, "cvs" = rwADF_results_reduced$cvs_95)
X$date <- as.Date(X$date)

p_rwadf_reduced <- ggplot(data = X, aes(x = date)) +
  geom_line(aes(y = test, colour = "The rolling window ADF statistic sequence (left axis)")) +
  geom_line(aes(y = cvs, colour = "The 95% critical value sequence (left axis)")) +
  geom_line(aes(y = price/1000, colour = "CriX (right axis)")) +
  scale_y_continuous(limits = c(-5, 60), 
                     sec.axis = sec_axis(~.*1000, name = NULL)) +
  #scale_x_date(breaks = pt, labels = c("2014", "2015", "2016", "2017", "2018") ) +
  labs(y = NULL, x = "Date") + 
  ggtitle("Rolling-window ADF- Reduced sample")

p_rwadf_reduced <- p_rwadf_reduced + theme(
  panel.background =  element_blank(),
  panel.border = element_rect(linetype = 1, colour = "black", fill = NA),
  panel.grid.major = element_line(linetype = 2, color = "grey90"),
  legend.title=element_blank(),
  legend.position = c(0.2, 0.9),
  legend.background = element_rect(fill = "transparent", colour = "transparent"),
  legend.key = element_rect(fill = "transparent", colour = "transparent"),
  plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
)

for(i in 1:nrow(stamps_rwADF_reduced)){
  p_rwadf_reduced <- p_rwadf_reduced + annotate("rect", xmin = X$date[stamps_rwADF_reduced[i,1]], xmax = X$date[stamps_rwADF_reduced[i,2]], 
                                                ymin= -5, ymax = 60,alpha=.35, fill="yellow")
}

pdf("./graphs/p_rwadf_reduced_95.pdf")
p_rwadf_reduced
dev.off() 





#SADF Critical values----


m <- 250 # number of replications
n <- nrow(data)
r0 <- 0.01 + 1.8 / sqrt(n)
swindow0 <- floor(n*r0)
d <- n - swindow0 + 1

badf.stat <-  MultipleBubbles::badf(m, n)$values


save(badf.stat, file = "./results/badf_MonteCarlo.RData")
#load("./results/badf_MonteCarlo.RData")

qe <- c(0.9, 0.95, 0.99)
quantile.badfs <- matrix(NA, nrow = length(qe), ncol = ncol(badf.stat))

for(i in 1:ncol(badf.stat)){
  quantile.badfs[, i] <- quantile(badf.stat[,i], qe, na.rm = T)
}

rownames(quantile.badfs) <- as.character(qe)

#GSADF CVs-----

m <- 100 # number of replications
n <- nrow(data)
r0 <- 0.01 + 1.8 / sqrt(n)
swindow0 <- floor(n*r0)
d <- n - swindow0 + 1

bsadf.stat <-  MultipleBubbles::bsadf(m, n)$values

save(bsadf.stat, file = "./results/bsadf_MonteCarlo.RData")
#load("./results/bsadf_MonteCarlo.RData")


qe <- c(0.9, 0.95, 0.99)
quantile.bsadfs <- matrix(NA, nrow = length(qe), ncol = ncol(bsadf.stat))

for(i in 1:ncol(bsadf.stat)){
  quantile.bsadfs[, i] <- quantile(bsadf.stat[,i], qe, na.rm = T)
}

rownames(quantile.bsadfs) <- as.character(qe)


#B(S)ADF----

var.name <- colnames(data)[-1]

sadf.gsadf.res <- matrix(0, 2, length(var.name))
colnames(sadf.gsadf.res) <- var.name

bsadf.stat.seq <- badf.stat.seq <- list()


  res.stat_full <- list()
  res.stat_reduced <- list()
  
  
  for(sample in c("full", "reduced")){
    if(sample == "full"){
      res.stat_full <- MultipleBubbles::sadf_gsadf(data$price, adflag = 7, mflag = 1, IC = 2, parallel = "TRUE")
} else {
        res.stat_reduced <- MultipleBubbles::sadf_gsadf(data$price[365:nrow(data)], adflag = 7, mflag = 1, IC = 2, parallel = "TRUE")
}}
  
  save(res.stat_full, file = "./results/b(s)adf_stat_full.RData")
  save(res.stat_reduced, file = "./results/b(s)adf_stat_reduced.Rdata")
  
  
#  sadf.gsadf.res[, ts] = c(res.stat$sadf, res.stat$gsadf)
#  bsadf.stat.seq[[ts]] = res.stat$bsadfs
#  badf.stat.seq[[ts]] = res.stat$badfs
  


stamps_badf_full  <- date.stamp(stat = res.stat_full[[1]]$badfs, cv = t(t(quantile.badfs[2,])))
stamps_bsadf_full <- date.stamp(stat = res.stat_full[[1]]$bsadfs, cv = t(t(quantile.bsadfs[2,])))

stamps_badf_reduced  <- date.stamp(stat = res.stat_reduced[[1]]$badfs, cv = t(t(quantile.badfs[2,1:length(res.stat_reduced[[1]]$badfs)])))
stamps_bsadf_reduced <- date.stamp(stat = res.stat_reduced[[1]]$bsadfs, cv = t(t(quantile.bsadfs[2,1:length(res.stat_reduced[[1]]$badfs)])))


### B(S)ADF Graphs----------

#P_badf_full----
X <-  data.frame(data[r0:nrow(data),] ,"test" = res.stat_full[[1]]$badfs, "cvs" = t(t(quantile.badfs[2,])))   #Generalna prica, posle cu da prilagodim
X$date <- as.Date(X$date)

p_badf_full <- ggplot(data = X, aes(x = date)) +
  geom_line(aes(y = test, colour = "The ADF statistic sequence (left axis)")) +
  geom_line(aes(y = cvs, colour = "The 95% critical value sequence (left axis)")) +
  geom_line(aes(y = price/1000, colour = "CriX (right axis)")) +
  scale_y_continuous(limits = c(-5, 60), 
                     sec.axis = sec_axis(~.*1000, name = NULL)) +
  #scale_x_date(breaks = pt, labels = c("2014", "2015", "2016", "2017", "2018") ) +
  labs(y = NULL, x = "Date") + 
  ggtitle("PWY- Full sample")

p_badf_full <- p_badf_full + theme(
  panel.background =  element_blank(),
  panel.border = element_rect(linetype = 1, colour = "black", fill = NA),
  panel.grid.major = element_line(linetype = 2, color = "grey90"),
  legend.title=element_blank(),
  legend.position = c(0.2, 0.9),
  legend.background = element_rect(fill = "transparent", colour = "transparent"),
  legend.key = element_rect(fill = "transparent", colour = "transparent"),
  plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
)

for(i in 1:nrow(stamps_badf_full)){
  p_badf_full <- p_badf_full + annotate("rect", xmin = X$date[stamps_badf_full[i,1]], xmax = X$date[stamps_badf_full[i,2]], 
                                        ymin= -5, ymax = 60,alpha=.35, fill="yellow")
}

pdf("./graphs/p_badf_full.pdf")
p_badf_full
dev.off() 


#p_badf_reduced----
data_reduced <- data[365:nrow(data),]
windowSize_reduced <- round(nrow(data_reduced)*(0.01 + 1.8 / sqrt(nrow(data_reduced))))

X <-  data.frame(data_reduced[windowSize_reduced:nrow(data_reduced),] ,"test" = res.stat_reduced[[1]]$badfs, "cvs" = t(t(quantile.badfs[2,1:length(res.stat_reduced[[1]]$badfs)])))   #Generalna prica, posle cu da prilagodim
X$date <- as.Date(X$date)

p_badf_reduced <- ggplot(data = X, aes(x = date)) +
  geom_line(aes(y = test, colour = "The ADF statistic sequence (left axis)")) +
  geom_line(aes(y = cvs, colour = "The 95% critical value sequence (left axis)")) +
  geom_line(aes(y = price/1000, colour = "CriX (right axis)")) +
  scale_y_continuous(limits = c(-5, 60), 
                     sec.axis = sec_axis(~.*1000, name = NULL)) +
  #scale_x_date(breaks = pt, labels = c("2014", "2015", "2016", "2017", "2018") ) +
  labs(y = NULL, x = "Date") + 
  ggtitle("PWY- Reduced sample")

p_badf_reduced <- p_badf_reduced + theme(
  panel.background =  element_blank(),
  panel.border = element_rect(linetype = 1, colour = "black", fill = NA),
  panel.grid.major = element_line(linetype = 2, color = "grey90"),
  legend.title=element_blank(),
  legend.position = c(0.2, 0.9),
  legend.background = element_rect(fill = "transparent", colour = "transparent"),
  legend.key = element_rect(fill = "transparent", colour = "transparent"),
  plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
)

for(i in 1:nrow(stamps_badf_reduced)){
  p_badf_reduced <- p_badf_reduced + annotate("rect", xmin = X$date[stamps_badf_reduced[i,1]], xmax = X$date[stamps_badf_reduced[i,2]], 
                                              ymin= -5, ymax = 60,alpha=.35, fill="yellow")
}

pdf("./graphs/p_badf_reduced.pdf")
p_badf_reduced
dev.off() 


#P_bsadf_full----
X <-  data.frame(data[r0:nrow(data),] ,"test" = res.stat_full[[1]]$bsadfs, "cvs" = t(t(quantile.bsadfs[2,])))   #Generalna prica, posle cu da prilagodim
X$date <- as.Date(X$date)

p_bsadf_full <- ggplot(data = X, aes(x = date)) +
  geom_line(aes(y = test, colour = "The BSADF statistic sequence (left axis)")) +
  geom_line(aes(y = cvs, colour = "The 95% critical value sequence (left axis)")) +
  geom_line(aes(y = price/1000, colour = "CriX (right axis)")) +
  scale_y_continuous(limits = c(-5, 60), 
                     sec.axis = sec_axis(~.*1000, name = NULL)) +
  #scale_x_date(breaks = pt, labels = c("2014", "2015", "2016", "2017", "2018") ) +
  labs(y = NULL, x = "Date") + 
  ggtitle("PSY- Full sample")

p_bsadf_full <- p_bsadf_full + theme(
  panel.background =  element_blank(),
  panel.border = element_rect(linetype = 1, colour = "black", fill = NA),
  panel.grid.major = element_line(linetype = 2, color = "grey90"),
  legend.title=element_blank(),
  legend.position = c(0.2, 0.9),
  legend.background = element_rect(fill = "transparent", colour = "transparent"),
  legend.key = element_rect(fill = "transparent", colour = "transparent"),
  plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
)

for(i in 1:nrow(stamps_bsadf_full)){
  p_bsadf_full <- p_bsadf_full + annotate("rect", xmin = X$date[stamps_bsadf_full[i,1]], xmax = X$date[stamps_bsadf_full[i,2]], 
                                          ymin= -5, ymax = 60,alpha=.35, fill="yellow")
}

pdf("./graphs/p_bsadf_full.pdf")
p_bsadf_full
dev.off() 


#p_bsadf_reduced----
data_reduced <- data[365:nrow(data),]
windowSize_reduced <- round(nrow(data_reduced)*(0.01 + 1.8 / sqrt(nrow(data_reduced))))

X <-  data.frame(data_reduced[windowSize_reduced:nrow(data_reduced),] ,"test" = res.stat_reduced[[1]]$bsadfs, "cvs" = t(t(quantile.bsadfs[2,1:length(res.stat_reduced[[1]]$bsadfs)])))   #Generalna prica, posle cu da prilagodim
X$date <- as.Date(X$date)

p_bsadf_reduced <- ggplot(data = X, aes(x = date)) +
  geom_line(aes(y = test, colour = "The BSADF statistic sequence (left axis)")) +
  geom_line(aes(y = cvs, colour = "The 95% critical value sequence (left axis)")) +
  geom_line(aes(y = price/1000, colour = "CriX (right axis)")) +
  scale_y_continuous(limits = c(-5, 60), 
                     sec.axis = sec_axis(~.*1000, name = NULL)) +
  #scale_x_date(breaks = pt, labels = c("2014", "2015", "2016", "2017", "2018") ) +
  labs(y = NULL, x = "Date") + 
  ggtitle("PSY- Reduced sample")

p_bsadf_reduced <- p_bsadf_reduced + theme(
  panel.background =  element_blank(),
  panel.border = element_rect(linetype = 1, colour = "black", fill = NA),
  panel.grid.major = element_line(linetype = 2, color = "grey90"),
  legend.title=element_blank(),
  legend.position = c(0.2, 0.9),
  legend.background = element_rect(fill = "transparent", colour = "transparent"),
  legend.key = element_rect(fill = "transparent", colour = "transparent"),
  plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
)

for(i in 1:nrow(stamps_bsadf_reduced)){
  p_bsadf_reduced <- p_bsadf_reduced + annotate("rect", xmin = X$date[stamps_bsadf_reduced[i,1]], xmax = X$date[stamps_bsadf_reduced[i,2]], 
                                                ymin= -5, ymax = 60,alpha=.35, fill="yellow")
}

pdf("./graphs/p_bsadf_reduced.pdf")
p_bsadf_reduced
dev.off() 





