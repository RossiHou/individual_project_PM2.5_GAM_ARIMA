library(ggplot2)
library(tidyr)
library(GGally)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(car)
library(PerformanceAnalytics)
library(mgcv)
library(nlme)
library(splines)
library(mgcViz)
library(ggplot2)
library(qgam)
library(tidyverse)
library(forecast)
library(tseries)


project <- read.csv("E:/BCNCS/BL5233/air_pollution1.csv")


time_interval <- data.frame(rep(1,2))
colnames(time_interval) <- c('time_interval')

project <- cbind(project,time_interval)



ggplot(project, aes(y=project$PM2.5,x=1)) + geom_violin(fill="#F8766D", width = 0.5) + 
  geom_boxplot(width = 0.1) + 
  stat_summary(fun=mean, geom="point", shape=20, size=4)+
  theme() + 
  labs(x="", y="pm2.5") +
  ggtitle("Boxplot and Density of PM2.5") + coord_flip() +
  annotate("text", x= 1.165, y=900, color="black",size = 4,
           label=paste("Min: ",fivenum(project$PM2.5)[1],"\n",
                       "Median: ",fivenum(project$PM2.5)[3],"\n",
                       "Mean:",round(mean(na.omit(project$PM2.5)),1),"\n",
                       "Max: ",fivenum(project$PM2.5)[5],"\n",
                       "NA:",sum(is.na(project$PM2.5)))) 

ks.test(project$PM2.5,"pnorm")
shapiro.test(project$PM2.5)
##非正态分布
dim(project)
##[1] 33648    13
project <- drop_na(project)
dim(project)
##[1] 31626    13
project <- project %>% mutate_at(vars(year, month, day), funs(factor))
str(project)


ggboxplot(project, x="month", y="PM2.5", color="month", legend="none", outlier.size=3, alpha=0.1) +
  stat_summary(fun=mean, geom="line", aes(group=1), color = "blue", size=1.1)+ ylim(0,500)+ 
  ggtitle("PM2.5 - month")

ggboxplot(project, x="year", y="PM2.5", color="year", legend="none", outlier.size=3, alpha=0.1) +
  stat_summary(fun=mean, geom="line", aes(group=1), color = "blue", size=1.1)+ ylim(0,500)+ 
  ggtitle("PM2.5 - year")



grid.arrange(
  ggplot(project, aes(x = PM2.5, fill = year)) + geom_density(alpha = 0.55) + ggtitle("PM2.5 density vs year") + theme(legend.position = "none"),
  ggplot(project, aes(x = PM2.5, fill = month)) + geom_density(alpha = 0.5) + ggtitle("PM2.5 density vs month")+ theme(legend.position = "none"),
  ggplot(project, aes(x = PM2.5, fill = CO)) + geom_density(alpha = 0.5) + ggtitle("PM2.5 density vs co")+ theme(legend.position = "none"),
  ggplot(project, aes(x = PM2.5, fill = PM10)) + geom_density(alpha = 0.5)+ ggtitle("PM2.5 density vs pm10")+ theme(legend.position = "none"),
  ggplot(project, aes(x = PM2.5, fill = TEMP)) + geom_density(alpha = 0.5)+ ggtitle("PM2.5 density vs temp")+ theme(legend.position = "none"),
  ggplot(project, aes(x = PM2.5, fill = PRES)) + geom_density(alpha = 0.5)+ ggtitle("PM2.5 density vs PRES")+ theme(legend.position = "none"),
  ggplot(project, aes(x = PM2.5, fill = DEWP)) + geom_density(alpha = 0.5)+ ggtitle("PM2.5 density vs DEWP")+ theme(legend.position = "none"),
  ggplot(project, aes(x = PM2.5, fill = WSPM)) + geom_density(alpha = 0.5)+ ggtitle("PM2.5 density vs WSPM")+ theme(legend.position = "none"),
  ncol=3
)

chart.Correlation(project %>% select(6:12), histogram=TRUE, pch=19)


ggplot(project) +
  geom_smooth(aes(x=TEMP, y=PRES))+
  facet_grid(.~month)+
  ggtitle("TEMP vs. PRES across month")+
  theme(legend.title=element_blank())





model1 <- gam(PM2.5~s(CO,bs="cs")+
                s(PM10,bs="cs")+
                s(TEMP,bs="cs")+
                s(PRES,bs="cs")+
                s(DEWP,bs="cs")+
                s(WSPM,bs="cs")+
                month+offset(log(time_interval)),
              method = "GCV.Cp",family=quasipoisson(link = "log"),data=projecttest_bind)

gam.check(model1)
summary(model1)

viz1 <- getViz(model1)

print(plot(viz1, allTerms = T), pages = 1)

model1.2 <-gam(PM2.5~s(CO,bs="cs")+
                 s(PM10,bs="cs")+
                 s(TEMP,bs="cs")+
                 s(PRES,bs="cs")+
                 s(DEWP,bs="cs")+
                 s(WSPM,bs="cs")+
                 s(CO,PM10)+
                 s(DEWP,PM10)+
                 s(CO,DEWP)+offset(log(time_interval))
                 month
               ,method = "GCV.Cp",family=quasipoisson(link = "log"),data=project)

summary(model1.2)

par(mfrow=c(2,2))


gam.check(model1.2)

vis.gam(model1.2, view = c("CO", "PM10"),
        theta = 50, n.grid = 50, lwd = 0.4,zlab = "",
        ticktype = "detailed", color = "topo")

vis.gam(model1.2, view = c("DEWP", "PM10"),
        theta = 50, n.grid = 50, lwd = 0.4,zlab = "",
        ticktype = "detailed", color = "topo")

vis.gam(model1.2, view = c("DEWP", "CO"),
        theta = 50, n.grid = 50, lwd = 0.4,zlab = "",
        ticktype = "detailed", color = "topo")

viz1.2 <- getViz(model1.2)

print(plot(viz1.2, allTerms = T), pages = 1)

viz1.3 <- getViz(model1.3)

print(plot(viz1.3, allTerms = T), pages = 1)



#########################################start:2013.3.1,end:2016.12.18
par(mfrow=c(1,1))
monthseries <- aggregate(project, by=list(project$day,project$month,project$year), FUN=mean)
myseries <- ts(monthseries$PM2.5,start=c(2013,59),end=c(2016,352),frequency=365)
plot(myseries)

###check stability
ndiffs(myseries)
dmyseries <- diff(myseries,differences = 1,lag=1)
adf.test(myseries)
##<0.05,no problem
##lag=0
##don't need to do that, just for the demonstration
######choose the model
Acf(myseries)
Pacf(myseries)
##0,0,0
arimafit <- arima(myseries,order=c(0,0,0))
arimafit

arimatest <- auto.arima(dmyseries)

qqPlot(arimafit$residuals)

Box.test(arimatest$residuals,type = "Ljung-Box")


#### forecast 2017.2.19
predictdata <- read.csv("E:/BCNCS/BL5233/individual/2007predictdata.csv")

predictdata <- predictdata %>% mutate_at(vars(year, month, day), funs(factor))
str(predictdata)
forecast(arimatest,63)###########126.68
aggregate(predictdata, by=list(predictdata$day,predictdata$month,predictdata$year), FUN=mean)########109.25
plot(forecast(arimatest,63))
####
predictdata$predict <- exp(predict.gam(model1.2,newdata = predictdata))
aggregate(predictdata, by=list(predictdata$day,predictdata$month,predictdata$year), FUN=mean)

##combine this 2 result. you will get a good prediction result