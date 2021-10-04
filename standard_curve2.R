#R code for creation of the standard curves found in Figure 2

setwd("loation/of/data")

library(ggplot2)

data = read.csv("file with creosote reads converted to relative abundance", header=T)

diet.data2 <- data[which(data$plant=="yes"),]

#standard curve from diet only samples, undigested 
diet.mod <- lm(Percent.Creo ~ 0 + Percent.Diet, data = diet.data2)
summary(diet.mod)

#prediction line
diet.predict <- ggpredict(diet.mod, "Percent.Diet", type="random" ) 
diet.predicted <- data.frame(Percent.Diet=diet.predict$x, Pred.y= diet.predict$predicted)

new.dat <- data.frame(Percent.Diet = seq(min(0),max(0.6), by = 0.01))

#generate prediction intervals on new data
diet.prediction.intervals <- predict(diet.mod, newdata = new.dat, interval = 'prediction')
diet.prediction.intervals <- as.data.frame(diet.prediction.intervals)

#generate confidence intervals on new data 
diet.confidence.intervals <- predict(diet.mod, newdata = new.dat, interval = 'confidence')
diet.confidence.intervals <- as.data.frame(diet.confidence.intervals )
names(diet.confidence.intervals)[1] <- "Percent.Creo"

#add to data for plotting
new.dat <- cbind(new.dat, diet.prediction.intervals)
names(new.dat )[2] <- "Percent.Creo.Predicted "
names(new.dat )[3] <- "lwr.pred "
names(new.dat )[4] <- "upr.pred "

#add to data for plotting
new.dat <- cbind(new.dat, diet.confidence.intervals)
names(new.dat )[5] <- "Percent.Creo.Predicted2 "
names(new.dat )[6] <- "lwr.conf"
names(new.dat )[7] <- "upr.conf "

#plot standard curve with prediction and confidence intervals
fig.3A <- ggplot()+
  geom_point(data=diet.data2, aes(y=Percent.Creo, x= Percent.Diet, shape=plant))+ scale_shape_manual(values=c(15)) +
  geom_line(data=diet.predicted,aes(y=Pred.y, x= Percent.Diet), color="black", size=1) +
  geom_ribbon(data=new.dat, aes(x=new.dat$Percent.Diet, ymin = new.dat$lwr.pred, ymax = new.dat$upr.pred), 
              fill = "orange",   alpha = 0.2) +
  xlab("Creosote in Diet") + ylab("Creosote Recovered in Diet Sample") +
  theme_classic() + 
  geom_ribbon(data=new.dat, aes(x=Percent.Diet, ymin=new.dat$lwr.conf, ymax=new.dat$upr.conf), fill=NA, linetype=2, color="black")
fig.3A

fig.3A2 <- fig.3A +  scale_x_continuous(labels = percent_format(accuracy=1), breaks = seq(0, .6, by = 0.10), limits=c(0,.6))
fig.3A2

fig.3A3 <- fig.3A2 + coord_cartesian(ylim = c(0, 1), xlim = c(0,.6)) 
fig.3A3

fig.3A4 <- fig.3A3 + scale_y_continuous(labels = percent_format(accuracy=1), breaks = seq(0, 1, by = 0.10))
fig.3A4

fig.3A5 <- fig.3A4 + theme(legend.position = "none")
fig.3A5


#standard curve produced from fecal samples
rat.mod <- lm(Percent.Creo ~ 0 + Percent.Diet, data = data)
summary(rat.mod)

#prediction line
rat.predict <- ggpredict(rat.mod, "Percent.Diet", type="random" ) 
rat.predicted <- data.frame(Percent.Diet=rat.predict$x, Pred.y= rat.predict$predicted)

new.dat2  <- data.frame(Percent.Diet = seq(min(0),max(0.6), by = 0.01))

#generate prediction intervals on new data
rat.prediction.intervals <- predict(rat.mod, newdata = new.dat2, interval = 'prediction')

#generate confidence intervals on new data
rat.confidence.intervals <- predict(rat.mod, newdata = new.dat2, interval = 'confidence')

#add to data 
new.dat2 <- cbind(new.dat2, rat.prediction.intervals)
names(new.dat2)[2] <- "Percent.Creo.Predicted "
names(new.dat2)[3] <- "lwr.pred "
names(new.dat2)[4] <- "upr.pred "

#add to data
new.dat2 <- cbind(new.dat2, rat.confidence.intervals)
names(new.dat2)[5] <- "Percent.Creo.Predicted2 "
names(new.dat2)[6] <- "lwr.conf"
names(new.dat2)[7] <- "upr.conf "

#plot standard curve with confidence and prediction intervals
fig.3B <- ggplot()+
  geom_point(data=data, aes(y=Percent.Creo, x= Percent.Diet, shape=trial)) +
  geom_line(data=rat.predicted,aes(y=Pred.y, x= Percent.Diet), color="black", size=1) +
  geom_ribbon(data=new.dat2, aes(x=Percent.Diet, ymin = new.dat2$lwr.pred, ymax = new.dat2$upr.pred), 
              fill = "blue",   alpha = 0.2) +
  xlab("Creosote in Diet") + ylab("Creosote Recovered in Feces") +
  theme_classic() + 
  geom_ribbon(data=new.dat2, aes(x=Percent.Diet, ymin=new.dat2$lwr.conf, ymax=new.dat2$upr.conf), fill=NA, linetype=2, color="black")
fig.3B

fig.3B2 <- fig.3B +  scale_x_continuous(labels = percent_format(accuracy=1), breaks = seq(0, .6, by = 0.10), limits=c(0,.6))
fig.3B2

fig.3B3 <- fig.3B2 + coord_cartesian(ylim = c(0, 1), xlim = c(0,.6)) 
fig.3B3

fig.3B4 <- fig.3B3 + scale_y_continuous(labels = percent_format(accuracy=1), breaks = seq(0, 1, by = 0.10))
fig.3B4

fig.3B5 <- fig.3B4 + theme(legend.position = "none")
fig.3B5

