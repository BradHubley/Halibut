fs3<-function()
{
#used in framework assessment.
#MUST LOAD library(MASS)
# fixed station CPUE analysis where stations have been sampled > 4 years
# DO NOT USE STRATA IN THIS ANALYSIS
#same as fs1 excpet no sql query and dropped all preliminary analyses
#see fs.bothwgt.sql
#Check the data out
#table(fs.dat$YEAR,fs.dat$STATION)
#temp<-table(fs.dat$YEAR,fs.dat$STATION)
#rowSums(temp)
#colSums(temp)
#min(colSums(temp))
#max(colSums(temp))

postscript("plot.fs3.ps",horizontal=FALSE)
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(2.5,3,4,.5),las=1,xpd=NA,cex=1)

####################################################################################################
# All stations in all nafo areas where stations sampled > 4 years
fs.dat<-read.table("fs.2010.txt",header=T)		#used for plotting 2010 estimate
x<-fs.dat[-dim(fs.dat)[1],]
names(x)<-c('year','nafo','station','stratum','duration','hooks','bothwgt')
x<-x[x$year<=2009,]

x$year<-as.factor(x$year)
x$station<-as.factor(x$station)
x$stratum<-as.factor(x$stratum)
# offset with hooks or log hooks.  Steve Smith suggested log hooks because link is log

#a1<-glm(bothwgt~year+station,offset=log(hooks),data=x,family=negative.binomial(theta = 1, link = "log"))
#a1<-glmmPQL(bothwgt~year+station+offset(hooks),random=~1|stratum/station,data=x,niter = 10,family=negative.binomial(theta = .37, link = "log"))
a1<-glm.nb(bothwgt~year+station+offset(log(hooks)),data=x)

xv<-seq(1998,2009,1)
p.s1<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(3,length(xv)))),type="response",se.fit=T)
p.s2<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(33,length(xv)))),type="response",se.fit=T)
p.s3<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(150,length(xv)))),type="response",se.fit=T)
#select out just year effects
#(a1$coef)[grep("year",names(a1$coef))]

#calculate CI for one predicted line
#say for p.s1
pr<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(33,length(xv)))),se.fit=T)
family <- family(a1) 
lower.all <- family$linkinv(pr$fit - qnorm(0.95) * pr$se.fit) 
upper.all <- family$linkinv(pr$fit + qnorm(0.95) * pr$se.fit) 
pr.all <- family$linkinv(pr$fit) 

temp1<-summary(a1)
temp2<-anova(a1,test='F')
temp3<-data.frame(x,predict(a1))

ylim1<-c(min(c(p.s1$fit,p.s2$fit,p.s3$fit),na.rm=T),max(c(p.s1$fit,p.s2$fit,p.s3$fit),na.rm=T))

plot(xv,p.s1$fit,type='l',lty=4,ylim=ylim1,xlab='',ylab='',lwd=3)
lines(xv,p.s2$fit,lty=2,lwd=3)
lines(xv,p.s3$fit,lty=1,lwd=3)
text(2009.2,p.s1$fit[10],"2",cex=1)
text(2009.2,p.s2$fit[10],"33",cex=1)
text(2009.17,p.s3$fit[10],"141",cex=1)
#mtext("Year",side=1,line=2.5,cex=1.1)
mtext("Standarized catch rate",side=2,line=3.5,cex=1.1,las=0)
mtext("3NOPs4VWX",side=3,line=.5,cex=1.3,las=0)

ylim1<-c(min(lower.all),max(upper.all))

plot(xv,pr.all,type='l',lty=1,ylim=ylim1,xlab='',ylab='',lwd=3)
lines(xv,upper.all,lty=2,lwd=3)
lines(xv,lower.all,lty=2,lwd=3)
#mtext("Year",side=1,line=2.5,cex=1.1)
mtext("Standarized catch rate",side=2,line=3.5,cex=1.1,las=0)
mtext("3NOPs4VWX",side=3,line=.5,cex=1.3,las=0)

#box plot station effects

write.table(x, file="fs.csv", sep=",", row.names=FALSE, col.names=TRUE)

print("OK")

####################################################################################################
# All stations in 4VWX only
x<-fs.dat
names(x)<-c('year','nafo','station','stratum','duration','hooks','bothwgt')
x<-x[x$nafo=='4V'|x$nafo=='4W'|x$nafo=='4X',]
x$year<-as.factor(x$year)
x$station<-as.factor(x$station)
x$stratum<-as.factor(x$stratum)
# offset with hooks or log hooks.  Steve Smith suggested log hooks because link is log

#a1<-glm(bothwgt~year+station,offset=log(hooks),data=x,family=negative.binomial(theta = 1, link = "log"))
#a1<-glmmPQL(bothwgt~year+station+offset(hooks),random=~1|stratum/station,data=x,niter = 10,family=negative.binomial(theta = .37, link = "log"))
a1<-glm.nb(bothwgt~year+station+offset(log(hooks)),data=x)

xv<-seq(1998,2009,1)
p.s4<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(2,length(xv)))),type="response",se.fit=T)
p.s5<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(33,length(xv)))),type="response",se.fit=T)
p.s6<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(141,length(xv)))),type="response",se.fit=T)
#select out just year effects
#(a1$coef)[grep("year",names(a1$coef))]

temp1<-summary(a1)
temp2<-anova(a1,test='F')
#print(temp1)
#print(temp2)
temp3<-data.frame(x,predict(a1))

ylim1<-c(min(c(p.s4$fit,p.s5$fit,p.s6$fit),na.rm=T),max(c(p.s4$fit,p.s5$fit,p.s6$fit),na.rm=T))
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(2.5,3,4,.5),las=1,xpd=NA,cex=1)

plot(xv,p.s4$fit,type='l',lty=4,ylim=ylim1,xlab='',ylab='',lwd=3)
lines(xv,p.s5$fit,lty=2,lwd=3)
lines(xv,p.s6$fit,lty=1,lwd=3)
text(2009.2,p.s4$fit[10],"2",cex=1)
text(2009.2,p.s5$fit[10],"33",cex=1)
text(2009.17,p.s6$fit[10],"141",cex=1)
mtext("Year",side=1,line=2.5,cex=1.1)
mtext("Standarized catch rate",side=2,line=3.5,cex=1.1,las=0)
mtext("4VWX",side=3,line=.5,cex=1.3,las=0)

write.table(x, file="fs3.csv", sep=",", row.names=FALSE, col.names=TRUE)

##test vessel effects
ves1<-read.table("ves1.txt")
names(ves1)<-c('year','nafo','vessel','station','stratum','bothwt')
table(ves1$station,ves1$year)

x<-ves1
x$year<-as.factor(x$year)
x$station<-as.numeric(x$station)
x$stratum<-as.numeric(x$stratum)
x$vessel<-as.numeric(x$vessel)

temp<-table(x$station,x$year)

#exclude stations done 4 years or less
# as.numeric(names(rowSums(temp)[rowSums(temp)<5]))
x1<-x[x$station!=68&x$station!=91&x$station!=95&x$station!=103&x$station!=117&x$station!=120&x$station!=124&x$station!=125&x$station!=126&x$station!=133&x$station!=134&x$station!=139&x$station!=145&x$station!=148&x$station!=152&x$station!=157&x$station!=168&x$station!=169&x$station!=171&x$station!=173&x$station!=174&x$station!=175&x$station!=176&x$station!=177&x$station!=178&x$station!=179&x$station!=180&x$station!=181&x$station!=184&x$station!=187
&x$station!=188&x$station!=190&x$station!=191&x$station!=192&x$station!=193&x$station!=194&x$station!=195&x$station!=196&x$station!=197&x$station!=199&x$station!=200&x$station!=201&x$station!=202&x$station!=203&x$station!=204&x$station!=205&x$station!=206&x$station!=207&x$station!=208&x$station!=210&x$station!=211&x$station!=212&x$station!=213&x$station!=214&x$station!=215&x$station!=216&x$station!=217&x$station!=218&x$station!=219&x$station!=221
&x$station!=222&x$station!=223&x$station!=224&x$station!=225&x$station!=226&x$station!=227&x$station!=228&x$station!=229&x$station!=230&x$station!=232&x$station!=233&x$station!=234&x$station!=235&x$station!=236&x$station!=237&x$station!=238&x$station!=239&x$station!=240&x$station!=241&x$station!=243&x$station!=244&x$station!=245&x$station!=246&x$station!=247&x$station!=248&x$station!=249&x$station!=250&x$station!=252&x$station!=253&x$station!=254
&x$station!=255&x$station!=256&x$station!=257&x$station!=258&x$station!=259&x$station!=261&x$station!=262&x$station!=263&x$station!=264&x$station!=265&x$station!=266&x$station!=267&x$station!=317&x$station!=335,]

temp<-table(x1$station,x1$year)

x2<-x1
x2$year<-as.factor(x2$year)
x2$station<-as.numeric(x2$station)
x2$stratum<-as.numeric(x2$stratum)
x2$vessel<-as.numeric(x2$vessel)

print("OK")

####################################################################################################
# All stations vessesl >=3 years
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(2.5,3,4,.5),las=1,xpd=NA,cex=1)

fs.v3.dat<-read.table("fs.v3.txt",header=T)		#used for plotting 2010 estimate
x<-fs.v3.dat[-dim(fs.v3.dat)[1],]
names(x)<-c('year','nafo','station','stratum','duration','hooks','bothwgt')
x<-x[x$year<=2009,]

x$year<-as.factor(x$year)
x$station<-as.factor(x$station)
x$stratum<-as.factor(x$stratum)
a1<-glm.nb(bothwgt~year+station+offset(log(hooks)),data=x)
xv<-seq(1998,2009,1)

 p.s1<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(3,length(xv)))),type="response",se.fit=T)
 p.s2<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(33,length(xv)))),type="response",se.fit=T)
 p.s3<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(150,length(xv)))),type="response",se.fit=T)
 temp1<-data.frame(1998:2009,p.s2$fit,p.s2$se)
 names(temp1)<-c("year","fit","se")

a7<-lm(fit~year,data=temp1,weights=1/(se^2/222))
anova(a7,test='F')

pr<-predict(a1,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(33,length(xv)))),se.fit=T)
family <- family(a1) 
lower.4yr <- family$linkinv(pr$fit - qnorm(0.95) * pr$se.fit) 
upper.4yr <- family$linkinv(pr$fit + qnorm(0.95) * pr$se.fit) 
pr2.4yr <- family$linkinv(pr$fit) 

ylim1<-c(0,200)
plot(xv,pr2.4yr,type='l',lty=1,xlim=c(1998,2010),ylim=ylim1,xlab='',ylab='',lwd=3)
lines(xv,upper.4yr,lty=2,lwd=3)
lines(xv,lower.4yr,lty=2,lwd=3)
points(2010,55.424,pch=3)	#2010 estimate
#mtext("Year",side=1,line=2.5,cex=1.1)
mtext("Standarized catch rate",side=2,line=3.5,cex=1.1,las=0)
mtext("Vessels >3 years",side=3,line=.5,cex=1.3,las=0)
print("OK")


####################################################################################################
# 50 stations sampled every year since 1999
s50.dat<-read.table("50.stations.txt",header=T)		#USED IN 2009 ASSESSMENT
x<-s50.dat[-dim(s50.dat)[1],]
names(x)<-c('year','nafo','station','stratum','duration','hooks','bothwgt')
x<-x[x$year<=2009,]

#s50.dat<-read.table("57.stations.2010.txt",header=T)		#used for plotting 2010 estimate
#x<-s50.dat[-dim(s50.dat)[1],]


x$year<-as.factor(x$year)
x$station<-as.factor(x$station)
x$stratum<-as.factor(x$stratum)
# offset with hooks or log hooks.  Steve Smith suggested log hooks because link is log

#a1<-glm(bothwgt~year+station,offset=log(hooks),data=x,family=negative.binomial(theta = 1, link = "log"))
#a1<-glmmPQL(bothwgt~year+station+offset(hooks),random=~1|stratum/station,data=x,niter = 10,family=negative.binomial(theta = .37, link = "log"))
a1.50<-glm.nb(bothwgt~year+station+offset(log(hooks)),data=x)

xv<-seq(1998,2009,1)
temp1<-summary(a1.50)
temp2<-anova(a1.50,test='F')
#print(temp1)
#print(temp2)
temp3<-data.frame(x,predict(a1.50))


#calculate CI for one predicted line
#say for p.s4
pr.50<-predict(a1.50,data.frame(year=as.factor(xv),hooks=(rep(1,length(xv))),station=as.factor(rep(69,length(xv)))),se.fit=T)
family <- family(a1.50) 
lower.50 <- family$linkinv(pr.50$fit - qnorm(0.95) * pr.50$se.fit) 
upper.50 <- family$linkinv(pr.50$fit + qnorm(0.95) * pr.50$se.fit) 
pr2.50 <- family$linkinv(pr.50$fit) 

par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(2.5,3,4,.5),las=1,xpd=NA,cex=1)
#ylim1<-c(min(c(p.s4$fit,p.s5$fit,p.s6$fit),na.rm=T),max(c(p.s4$fit,p.s5$fit,p.s6$fit),na.rm=T))
ylim1<-c(0,200)

####################################################################################################
# Halibut survey and commercial index unstandardized
#all stations
years<-1998:2010
#VDC output
HS.4vwx<-c(47.82,38.86,43.75,41.12,42.06,38.48,36.51,40.63,41.25,52.03,48.14,71.27,79.37)
HS.n.4vwx<-c()
HS.se.4vwx<-c(114.95,63.32,67.95,64.1,73.42,68.54,66.91,58.87,46.7,88.46,81.59,97.68,117.53)
HS<-c(45.7,43.05,48.12,46.22,44.28,40.27,42.67,45.59,41.21,50.89,50.42,78.53,78.88)
HS.n<-c(161,167,216,190,192,187,214,173,172,234,275,205,215)
HS.se<-c(111.54,74.19,74.14,72.94,74.03,70.07,77.38,65.57,46.07,88.33,85.6,116.89,114.4)/sqrt(HS.n)
upper.HS<-HS+2*HS.se
lower.HS<-HS-2*HS.se

CI<-c(108.96,110.7,93.09,117.39,102.69,99.81,80.63,87.18,103.82,96.51,116.03,158.42,130.66)
CI.n<-c(540,478,686,511,623,646,853,517,465,473,798,462,541)
CI.se<-c(111.43,104.23,88.56,112.88,89.05,85.43,82.9,81.18,112.38,128.82,150.86,165.96,131.48)/sqrt(CI.n)
upper<-CI+CI.se*2
lower<-CI-CI.se*2

avg.50<-c(72.91,46.46,56.23,46.97,56.18,43.4,47.81,44.7,52.47,74.49,57.4,71.55,96.89)
avg.50.n<-c(34,rep(50,12))
avg.50.se<-c(101.18,64.13,71.39,61.18,82.54,71.02,70.99,55.83,49.5,87.1,65.9,60.66,104.69)/sqrt(avg.50.n)
upper.avg.50<-avg.50+avg.50.se*2
lower.avg.50<-avg.50-avg.50.se*2



ylim1<-c(0,200)
plot(years,CI,type='l',lty=1,xlim=c(1998,2010),ylim=ylim1,xlab='',ylab='',lwd=3)
arrows(years,upper,years,lower,length=0.05,angle=90,code=3)  # change length to 0.05 if want horizontal lines on error bar and code=3
lines(years,HS,lty=1,lwd=3,col='dark green')
#points(2010,HS[13],pch=3,col='dark green')	#2010 estimate
arrows(years,upper.HS,years,lower.HS,length=0.05,angle=90,code=3,col='dark green')  # change length to 0.05 if want horizontal lines on error bar and code=3
lines(years,avg.50,lty=1,lwd=3,col='gold')
arrows(years,upper.avg.50,years,lower.avg.50,length=0.05,angle=90,code=3,col='gold')  # change length to 0.05 if want horizontal lines on error bar and code=3
mtext("Standarized catch rate",side=2,line=3.5,cex=1.1,las=0)
mtext("(kg / 1000 hooks / 10 hrs)",side=2,line=2.5,cex=1.1,las=0)
mtext("A) Halibut survey and commercial index",side=3,line=.5,cex=1.3,las=0,adj=0)
legend(1998,ylim1[2],legend=c("Commercial index","All stations","50 stations"),col=c("black","dark green","gold"),lwd=2,bty="n",cex=.9)
#mtext("Commercial index",side=3,line=.5,cex=1.3,las=0,adj=0)

####################################################################################################
# Halibut survey standardized
plot(xv,pr2.all,type='l',lty=1,xlim=c(1998,2010),ylim=ylim1,xlab='',ylab='',lwd=3)
points(2010,56.50796,pch=3)	#2010 estimate
arrows(xv,upper.all,xv,lower.all,length=0.05,angle=90,code=3)  # change length to 0.05 if want horizontal lines on error bar and code=3

lines(xv,pr2.4yr,lty=1,lwd=3,col='red')
points(2010,55.424,pch=3,col='red')	#2010 estimate
arrows(xv,upper.4yr,xv,lower.4yr,length=0.05,angle=90,code=3,col='red')  # change length to 0.05 if want horizontal lines on error bar and code=3

lines(xv,pr2.50,lty=1,lwd=3,col='gold')
points(2010,59.31692,pch=3,col='gold')	#2010 estimate
arrows(xv,upper.50,xv,lower.50,length=0.05,angle=90,code=3,col='gold')  # change length to 0.05 if want horizontal lines on error bar and code=3

#lines(xv,HS[-13],lty=1,lwd=3,col='dark green')
#points(2010,HS[13],pch=3,col='dark green')	#2010 estimate
#arrows(xv,upper.HS,xv,lower.HS,length=0.05,angle=90,code=3,col='dark green')  # change length to 0.05 if want horizontal lines on error bar and code=3
#mtext("Year",side=1,line=2.5,cex=1.1)
mtext("Standarized catch rate",side=2,line=3.5,cex=1.1,las=0)
mtext("(kg / 1000 hooks)",side=2,line=2.5,cex=1.1,las=0)
mtext("B) Halibut survey",side=3,line=.5,cex=1.3,las=0,adj=0)
legend(1998,ylim1[2],legend=c("GLM","Vessels > 3yrs","50 stations GLM"),col=c("black","red","gold"),lwd=2,bty="n",cex=.9)

####################################################################################################
# RV survey
rv.dat<-read.table("rv.4vwx.txt",header=T)
rv.dat<-rv.dat[order(rv.dat[,1]),]
years<-1970:2010
ylim1<-c(0,0.8)
upper<-rv.dat[,5]+rv.dat[,6]*2
lower<-rv.dat[,5]-rv.dat[,6]*2

plot(years,rv.dat[,5],type='l',lty=1,xlim=c(1970,2010),ylim=ylim1,xlab='',ylab='',lwd=3)
arrows(years,upper,years,lower,length=0.05,angle=90,code=3)  # change length to 0.05 if want horizontal lines on error bar and code=3
mtext("Mean numbers per tow",side=2,line=3.5,cex=1.1,las=0)
mtext("RV survey",side=3,line=.5,cex=1.3,las=0,adj=0)

dev.off()


return()

}
