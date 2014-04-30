plot.vpop4r<-function(update.x=T,update.std=F)
{

if(update.x==T){vpop4r<<-getadmout("~/halibut/adm/vpop4r.rep")}
x<-vpop4r

if(update.std==T)
{ 
#hal.std<<-read.table("c:\\kurtis\\halibut\\adm\\vpop4r.std",skip=1,row.names=NULL,
#col.names=c("number","name","value","std.dev"))
}
#z<-hal.std

minsize=.1 #used for bubble plots
maxsize=3 

postscript("plot.vpop4r.ps",horizontal=FALSE)
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(3,3,3,.5),las=1,xpd=NA,cex=.9)
 
#LANDINGS
lyears<-seq(1970,x$last.year,1)

plot(lyears,x$C.ll[11:45],type="l",xlab="",ylab="")
mtext(side=3,text="Longline landings",cex=1.5,line=.5,adj=.1)
mtext(outer=T,side=2,text="Tones (1000s)",cex=1.5,line=.5,las=0,adj=.75)
#plot(lyears,x$C.ot,type="l",xlab="",ylab="")
#mtext(side=3,text="OT",cex=0.8,line=0,adj=.1)

#print("OK")

#Proportion of LANDINGS
par(mfcol=c(3,2),omi=c(1,1,1,1),mar=c(3,3,3,.5),las=1,xpd=NA,cex=.9)

lyears<-seq(1960,x$last.year,1)
plot(lyears,x$Foreign/(x$Foreign+x$Canada),type="l",xlab="",ylab="")
mtext(side=2,text="Proportion Foreign landings",line=3,outer=F,cex=1.2,las=0)
#mtext(outer=F,side=1,text="Years",cex=1.5,line=3)

plot(lyears,x$C.ll/(x$C.ot+x$C.ll),type="l",xlab="",ylab="")
mtext(side=2,text="Proportion long line",line=3,outer=F,cex=1.2,las=0)
mtext(outer=F,side=1,text="Years",cex=1.5,line=3)



#SURVEY SELECTIVITIES
par(mfrow=c(4,2),omi=c(1,1,1,1),mar=c(3,3,3,.5),las=1,xpd=NA,cex=.9)
plot(x$ages,x$old.summer.rv.F.sel,type="n",xlab="",ylab="", ylim=c(0,1))
lines(x$ages,x$old.summer.rv.F.sel,lty=1)
title(main="Summer RV")

plot(x$ages,x$halibut.fixed.F.sel,type="n",xlab="",ylab="", ylim=c(0,1))
lines(x$ages,x$halibut.fixed.F.sel,lty=1)
lines(x$ages,x$halibut.fixed.M.sel,lty=2)
title(main="Halibut fixed station")

mtext(outer=T,side=1,text="Age",cex=1.5,line=-23)
mtext(outer=T,side=2,text="Selectivity",cex=1.5,las=0,adj=.8)
#mtext(outer=T,side=3,text="Survey Selectivities",cex=1.4,line=2)
 
#COMMERCIAL SELECTIVITIES
#par(mfrow=c(4,2),omi=c(1,1,1,1),mar=c(3,3,3,.5),las=1,xpd=NA)
 
plot(x$ages,x$ll.F.sel,type="n",xlab="",ylab="", ylim=c(0,1))
lines(x$ages,x$ll.F.sel,lty=1)
lines(x$ages,x$ll.M.sel,lty=2)
title(main="Longline")

#plot(x$ages,x$ot.F.sel,type="n",xlab="",ylab="", ylim=c(0,1))
#lines(x$ages,x$ot.F.sel,lty=1)
#title(main="Otter trawl")
legend(13,-0.3,legend=c("females","males"), lty=c(1,2),bty="n",cex=1.1)

############################################################################################
###########################################################################################
#Commercial length-at-age
par(las=1,mfrow=c(2,1),omi=c(1.5,1,1,1),mar=c(4,4,2,1),cex=.9)

#LL
plot(c(1970,x$last.year),c(min(x$lens.comm),max(x$lens.comm)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$obs.catlen.comm[x$obs.catlen.comm[,2]==1&x$obs.catlen.comm[,4]==1,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]

prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens.com))
yy<-sort(rep(x$lens.com,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  
print("OK symbols")
#mtext(side=3,text="Females",adj=.1,line=0)  
mtext(side=3,text="LL Proportions at Length",adj=.1,outer=F,line=.5)  
#mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9)  
#print("OK")

################################################################################################
#SURVEY PROPORTIONS:
#1=summer_rv
#2=spring_rv
#3=fall_rv
#6=4VWcod_rv

#########################################################
#summer_rv
plot(c(x$first.year,x$last.year),c(min(x$lens),max(x$lens)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$obs.catlen.rv[x$obs.catlen.rv[,2]==1&x$obs.catlen.rv[,4]==1,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]

prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens))
yy<-sort(rep(x$lens,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  
print("OK symbols")

#mtext(side=3,text="Females",adj=.1,line=0)  
mtext(side=3,text="Summer Survey Proportions at Length",adj=.1,outer=F,line=.5)  
#mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9)  
########

##############################################################################3
#halibut fixed station
plot(c(1997,x$last.year),c(min(x$lens),max(x$lens)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$obs.catlen.rv[x$obs.catlen.rv[,2]==7&x$obs.catlen.rv[,4]==1,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]

prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens))
yy<-sort(rep(x$lens,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  
print("OK symbols")

mtext(side=3,text="Females",adj=.1,line=0)  
mtext(side=3,text="Halibut Survey Proportions at Length",adj=.1,outer=T)  
#mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9)  
########
##############################################################################
# SPAWNER NUMBERS, BIOMASS AND NUMBER OF RECRUITS
true.ssb=c(14.4876, 11.6345, 9.08794, 6.95595, 5.24195, 3.92409, 2.94761, 2.2323, 1.71039, 1.52231, 1.47736, 1.51722, 1.61053, 1.73303, 1.86594, 1.99476, 2.1084, 2.19767, 2.25663, 2.2824, 2.27364, 2.45615, 2.61134, 2.73171, 2.81077, 2.82252, 2.75128, 2.61125, 2.4209, 2.19902, 1.9642, 1.73418, 1.52535, 1.35221, 1.22695)
true.ssn=c(780, 603.288, 458.316, 343.361, 254.39, 189.182, 143.412, 110.924, 87.5108, 91.4936, 99.6097, 110.06, 121.486, 132.841, 143.309, 152.253, 159.19, 163.753, 165.723, 165.022, 161.68, 172.029, 179.997, 185.123, 187.055, 183.268, 172.894, 158.177, 141.061, 123.233, 106.149, 91.041, 78.9224, 70.5818, 66.5798)
true.r=c(200, 429.95, 488.656, 543.828, 593.265, 634.998, 667.362, 689.067, 699.248, 697.499, 683.89, 658.963, 623.712, 579.542, 528.214, 471.775, 412.474, 352.676, 294.765, 241.049, 193.67, 154.517, 125.15, 106.741, 100.023, 105.264, 122.256, 150.32, 188.338, 234.794, 287.837, 345.351, 405.044, 464.536, 521.455)
true.f=c(rep(.2,20),rep(.1,15))

par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)
years<-seq(x$first.year,x$last.year,1)
lyears<-seq(x$first.year+1,x$last.year,1)
#R<-x$SSNrec[2,][-c(1,length(x$SSNrec[2,]))]	#recruits  drop 1970 and last.year+1
R<-x$SSNrec[2,][-1]	#recruits  drop 1970

ages<-seq(x$first.age,x$last.age,1)

plot(years,x$SSNrec[1,]/1000,type="n",xlab="",ylab="",ylim=c(0,max(x$SSNrec[2,])*1.5/1000),axes=F)
lines(years,x$SSNrec[1,]/1000,lty=1,lwd=2)	#spawner numbers in millions 
lines(years,true.ssn/1000,lty=1,lwd=2,col='red')	#true spawner numbers in millions 
lines(lyears,R/1000,lty=4,lwd=1)  		#recruits in millions
lines(lyears,true.r[-1]*2/1000,lty=4,lwd=1,col='blue')  		#true recruits in millions
axis(1)
axis(2)

par(new=T)
plot(years,x$SSBrec[1,],type="n",xlab="",ylab="",ylim=c(0,max(x$SSBrec[1,])),axes=F)
lines(years,x$SSBrec[1,],lty=2,lwd=2)		# biomass in 1000 tons
lines(years,true.ssb,lty=2,lwd=2,col='green',xpd=F)		# true biomass in 1000 tons
axis(4)
box()
#mtext(side=3,text="Spawner numbers and biomass (Mature females), and recruits",cex=1.5,line=1)
mtext(side=2,text="S and R (millions)",cex=1.5,line=3,las=0)
mtext(side=4,text="SSB (1000t)",cex=1.5,line=3,las=0)
#legend(1993,max(x$SSBrec[1,])*.9,legend=c("Spawner numbers","Recruits","SSB"),lty=c(1,4,2),lwd=2,bty="n",cex=1.1)
#legend(1993,-max(x$SSBrec[1,])/4,legend=c("Spawner numbers","Recruits","SSB"),lty=c(1,4,2),lwd=2,bty="n",cex=1.1)
legend(1993,-max(x$SSBrec[1,]/4),legend=rev(c("Numbers","Recruits","Biomass")),lty=rev(c(1,4,2)),lwd=2,bty="n",cex=1.1)

print("SSB")

# TOTAL NUMBERS AND BIOMASS
true.total.biomass=c(33.5135, 26.6279, 20.7738, 16.1094, 12.5873, 10.1555, 8.5745, 7.66713, 7.26534, 7.22672, 7.42206, 7.72339, 8.05914, 8.37482, 8.63008, 8.79633, 8.85504, 8.79455, 8.6125, 8.31511, 7.91358, 8.12586, 8.18632, 8.09098, 7.84656, 7.471, 6.99342, 6.45563, 5.90479, 5.38532, 4.93737, 4.59559, 4.38801, 4.3353, 4.45011)
true.total.numbers=c(4200, 3855, 3786.59, 3915.09, 4175.94, 4494.42, 4815.21, 5109.88, 5356.19, 5537.37, 5641.63, 5661.86, 5595.39, 5443.66, 5212.01, 4909.23, 4547.14, 4140, 3703.95, 3256.39, 2815.16, 2500.34, 2188.04, 1902.93, 1665.39, 1491.58, 1393.39, 1378.39, 1449.9, 1607.06, 1845.12, 2155.67, 2527.16, 2945.4, 3394.15)

plot(years,x$Total.N/1000,type="n",xlab="",ylab="",ylim=c(0,max(x$Total.N)/1000),axes=F)
lines(years,x$Total.N/1000,lty=1,lwd=2)	#total numbers in millions
lines(years,true.total.numbers/1000,lty=1,lwd=2,col='red')	#True total numbers in millions
axis(1)
axis(2)
par(new=T)
plot(years,x$Total.W,type="n",xlab="",ylab="",ylim=c(0,max((x$Total.W)[-1])),axes=F)
lines(years[-1],(x$Total.W)[-1],lty=2,lwd=2)		#total biomass in 1000 t
lines(years[-1],true.total.biomass[-1],lty=2,lwd=2,col='green')		#total biomass in 1000 t
axis(4)
box()
#mtext(side=3,text="Total numbers and biomass",cex=0.9,line=1)
mtext(side=2,text="N (millions)",cex=1.5,line=3,las=0)
mtext(side=4,text="Total biomass (1000t)",cex=1.5,line=3,las=0)

############################################################################################
###########################################################################
#mean weights
w.f<-(x$Females%*%x$w[1,])/rowSums(x$Females)
w.m<-(x$Males%*%x$w[2,])/rowSums(x$Males)
#drop last year
w.f<-w.f[-length(w.f)]
w.m<-w.m[-length(w.m)]

plot(years,w.f,type="n",xlab="",ylab="",ylim=c(0,max(w.f)),axes=F)
axis(1)
axis(2)
lines(years,w.f,lty=1,lwd=2)	
lines(years,w.m,lty=2,lwd=2)	
box()
mtext(side=1,text="Year",cex=0.9,line=2.5)
mtext(side=2,text="Mean weights (kg)",cex=0.9,line=3.5,las=0)
legend(1971,max(w.f)*.9,legend=c("Females","Males"),lty=c(1,2),lwd=2,bty="n",cex=1.1)

###########################################################################
#mean exploitation
par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)

u.f<-apply(x$Females[,5:19]*x$F.utotal[10:45,5:19],1,sum)/rowSums(x$Females[,5:19])
u.m<-apply(x$Males[,5:19]*x$M.utotal[10:45,5:19],1,sum)/rowSums(x$Males[,5:19])
#drop last year
u.f<-u.f[-length(u.f)]
u.m<-u.m[-length(u.m)]
#drop first year
years<-years[-1]
u.f<-u.f[-1]
u.m<-u.m[-1]

plot(years,u.f,type="n",xlab="",ylab="",ylim=c(0,max(u.f)),axes=F)
axis(1)
axis(2)
lines(years,u.f,lty=1,lwd=2)	
#lines(years,u.m,lty=2,lwd=2)	
lines(years,true.f[-1],lty=2,lwd=2,col='red')	
box()
mtext(side=1,text="Year",cex=0.9,line=2.5)
mtext(side=2,text="Fishing mortality",cex=0.9,line=3.5,las=0)
legend(1971,max(u.f)*.9,legend=c("Estimated","True"),lty=c(1,2),col=c("black","red"),lwd=2,bty="n",cex=1.1)

##############################################################################
##############################################################################
par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)

# SPAWNER-RECRUIT RELATIONSHIP

#S<-seq(min(x$SSB),max(x$SSB),length=100)
S<-seq(0,50,length=1000)
#R<-2.71828*x$p2/x$p1*S*exp(-S/x$p1);	# Ricker reparamertized   
R<-2.71828*x$p2/x$p1*(S-x$p3)*exp(-(S-x$p3)/x$p1);	# Ricker reparamertized   with depensation; x intercept offset estiamted
#B-H
R2<-(x$alpha*(S-x$p3)/(1+x$alpha*(S-x$p3)/x$R0))
R<-R/1000
x$SSNrec[2,]<-x$SSNrec[2,]/1000
lyears<-seq(x$first.year+1,x$last.year,1)

#plot 1970-2006 SSB vs 1971-2007 recruitment
plot(x$SSB[1:(length(x$SSB)-2)],x$SSNrec[2,][-c(1,length(x$SSNrec[2,]))],type='n',xlab="",ylab="", ylim=c(0,max(x$SSNrec[2,])),xlim=c(min(x$SSB),max(x$SSB)))
points(x$SSB[1:(length(x$SSB)-2)],x$SSNrec[2,][-c(1,length(x$SSNrec[2,]))],pch=as.character(substring(lyears,3,3)),cex=.7)
lines(S,R,lty=1,lwd=2,xpd=F)  #Ricker
#lines(S,R2,lty=1,lwd=2,xpd=F)	#B-H
mtext(side=3,text="Spawner-recruitment relationship",cex=0.9,line=1)
mtext(side=1,text="SSB (1000t)",cex=0.9,line=2.5)
mtext(side=2,text="Age-1 recruits (millions)",cex=0.9,line=3.5,las=0)

######################################################################
#####################################################################

# predicted vs observed 
#pred at len by year
#print("OK")

par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)

plot.rv<-function(x,series.code) 
	{
	
	obs.s<-x$obs.catlen.rv[x$obs.catlen.rv[,2]==series.code,]
	pred.s<-x$pred.catlen.rv[x$obs.catlen.rv[,2]==series.code,]
	
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
	lens <- x$lens
	len.mat <- function(lens,YEARS)
	{
	#
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	#a <- seq(1,ages)
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
#	summeryears.s <- unique(obs.s[,2])
	
	obs.females <- obs.s[obs.s[,4]==1,]
	summeryears.f <- obs.females[,3]
	obs.females <- obs.females[,seq(5,x$n.lens+4)]

	obs.males <- obs.s[obs.s[,4]==2,]
	summeryears.m <- obs.males[,3]
	obs.males <- obs.males[,seq(5,x$n.lens+4)]

	pred.females <- pred.s[obs.s[,4]==1,]
	pred.males <- pred.s[obs.s[,4]==2,]
	
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	lens.m <- len.mat(lens,summeryears.m)
	years.m <- cllyear.mat(lens,summeryears.m)

	cex <- 0.6
	ylim <- c(-0.3,0.3)

#temp<-as.vector(obs.females)
#print(temp)
junk<-data.frame(as.vector(years.f),as.vector(lens.f),as.vector(obs.females))
#junk<-data.frame((rep(lens,dim(obs.females)[1])),rep(summeryears.f,each=length(lens)),as.vector(obs.females))
#names(junk)<-c('len','year','p')
names(junk)<-c('year','len','p')
#print(junk)
junk<<-junk
  for(y in 1:length(summeryears.f))
    {
	plot(lens.f[y,],obs.females[y,],cex=cex,ylab="",xlab="",type='p')
	lines(lens.f[y,],pred.females[y,],lty=1)
	mtext(side=3,outer=T,"Females",line=0,cex=1.1,adj=0)
	mtext(side=3,outer=F,summeryears.f[y],line=.2,cex=.8,adj=0)
 	
	plot(lens.m[y,],obs.males[y,],cex=cex,ylab="",xlab="",type='p')
	lines(lens.m[y,],pred.males[y,],lty=1)
	mtext(side=3,outer=T,"Males",line=0,cex=1.1,adj=1)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Proportions at length in summer survey",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Proportions at length in spring survey",cex=1)}
  if(series.code==3){mtext(side=3,outer=T,line=-.5,text="Proportions at length in fall survey",cex=1)}
  if(series.code==6){mtext(side=3,outer=T,line=-.5,text="Proportions at length in 4VsW cod survey",cex=1)}
  if(series.code==7){mtext(side=3,outer=T,line=-.5,text="Proportions at length in fixed station survey",cex=1)}
  if(series.code==8){mtext(side=3,outer=T,line=-.5,text="Proportions at length in commercial index",cex=1)}
   }
	}
#1=summer_rv
#2=spring_rv
#3=fall_rv
#6=4VWcod_rv
#7=fixed station
#8=commercial index

##########################################################################
#####################################################################

# RESIDUALS
#print("OK")

# residuals of obs - expected  (proportions)

par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)

plot.res<-function(x,series.code) 
	{
	
	obs.s<-x$obs.catlen.rv[x$obs.catlen.rv[,2]==series.code,]
	pred.s<-x$pred.catlen.rv[x$obs.catlen.rv[,2]==series.code,]
	
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
	lens <- x$lens
	len.mat <- function(lens,YEARS)
	{
	#
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	#a <- seq(1,ages)
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
#	summeryears.s <- unique(obs.s[,2])
	
	obs.females <- obs.s[obs.s[,4]==1,]
	summeryears.f <- obs.females[,3]

	obs.males <- obs.s[obs.s[,4]==2,]
	summeryears.m <- obs.males[,3]

	pred.females <- pred.s[obs.s[,4]==1,]
	pred.males <- pred.s[obs.s[,4]==2,]
	
	res.females <- obs.females[,seq(5,x$n.lens+4)]-pred.females
	res.males <- obs.males[,seq(5,x$n.lens+4)]-pred.males
	
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	lens.m <- len.mat(lens,summeryears.m)
	years.m <- cllyear.mat(lens,summeryears.m)

	cex <- 0.6
	ylim <- c(-0.3,0.3)


#print(dim(res.females))
	plot(lens.f,res.females,cex=cex,ylim=ylim,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
	
	plot(years.f,res.females,cex=cex,ylim=ylim,xlim=c(1960,2003),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)

	plot(lens.m,res.males,cex=cex,ylim=ylim,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
	mtext(side=1,outer=F,text="Length (cm)",cex=0.95,line=2.5)
	
	plot(years.m,res.males,cex=cex,ylim=ylim,xlim=c(1960,2003),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=1,outer=F,text="Year",cex=0.95,line=2.5)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in summer survey",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in spring survey",cex=1)}
  if(series.code==3){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in fall survey",cex=1)}
  if(series.code==6){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in 4VsW cod survey",cex=1)}
  if(series.code==7){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in fixed station survey",cex=1)}
  if(series.code==8){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in commercial index",cex=1)}

	}
#1=summer_rv
#2=spring_rv
#3=fall_rv
#6=4VWcod_rv

##########################################################################
#CATCH@LENGTH PREDICTED VS OBSERVED
par(mfrow=c(4,4),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)

plot.comm<-function(x,series.code) 
	{

	obs.s<-x$obs.catlen.comm[x$obs.catlen.comm[,2]==series.code,]
	pred.s<-x$pred.catlen.comm[x$obs.catlen.comm[,2]==series.code,]
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
	lens <- x$lens.comm
#	lens <- seq(10,278,3)
	len.mat <- function(lens,YEARS)
	{
	#
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	#a <- seq(1,ages)
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
#	summeryears.s <- unique(obs.s[,3])
	
	obs.females <- obs.s[obs.s[,4]==1,]
	summeryears.f <- obs.females[,3]
	obs.females<-obs.females[,seq(5,x$n.lens.comm+4)]

	obs.males <- obs.s[obs.s[,4]==2,]
	summeryears.m <- obs.males[,3]
	obs.males<-obs.males[,seq(5,x$n.lens.comm+4)]

	pred.females <- pred.s[obs.s[,4]==1,]
	pred.males <- pred.s[obs.s[,4]==2,]
		
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	lens.m <- len.mat(lens,summeryears.m)
	years.m <- cllyear.mat(lens,summeryears.m)
	
	cex <- 0.6
	ylim <- c(-0.3,0.3)
	
  if(series.code==1)
    {
  for(y in 1:length(summeryears.f))
    {
	plot(lens.f[y,],obs.females[y,],cex=cex,ylab="",xlab="",type='p')
	lines(lens.f[y,],pred.females[y,],lty=1)
	mtext(side=3,outer=T,"Females",line=0,cex=1.1,adj=0)
	mtext(side=3,outer=F,summeryears.f[y],line=.2,cex=.8,adj=0)
 	
	plot(lens.m[y,],obs.males[y,],cex=cex,ylab="",xlab="",type='p')
	lines(lens.m[y,],pred.males[y,],lty=1)
	mtext(side=3,outer=T,"Males",line=0,cex=1.1,adj=1)
	mtext(side=3,outer=T,line=2,text="Proportions at length in LL catch",cex=1)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Proportions at length in LL",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Proportions at length in OT",cex=1)}
    }
    }
  if(series.code==2)
    {
  for(y in 1:length(summeryears.f))
    {
	plot(lens.f[y,],obs.females[y,],cex=cex,ylab="",xlab="",type='p')
	lines(lens.f[y,],pred.females[y,],lty=1)
	mtext(side=3,outer=T,"Females",line=0,cex=1.1,adj=0)
	mtext(side=3,outer=F,summeryears.f[y],line=.2,cex=.8,adj=0)
        mtext(side=3,outer=T,line=2,text="Proportions at length in OT catch",cex=1)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Proportions at length in LL",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Proportions at length in OT",cex=1)}
    }
    }
	}
#//1=LL
#//2=OT

##########################################################################################
#RESIDUALS CATCH@LENGTH
par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
#print("OK")

plot.res.comm<-function(x,series.code) 
	{
	
	obs.s<-x$obs.catlen.comm[x$obs.catlen.comm[,2]==series.code,]
	pred.s<-x$pred.catlen.comm[x$obs.catlen.comm[,2]==series.code,]
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
	lens <- x$lens.comm
#	lens <- seq(10,278,3)
	len.mat <- function(lens,YEARS)
	{
	#
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	#a <- seq(1,ages)
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
#	summeryears.s <- unique(obs.s[,3])
	
	obs.females <- obs.s[obs.s[,4]==1,]
	summeryears.f <- obs.females[,3]
	obs.males <- obs.s[obs.s[,4]==2,]
	summeryears.m <- obs.males[,3]
	pred.females <- pred.s[obs.s[,4]==1,]
	pred.males <- pred.s[obs.s[,4]==2,]
	
	res.females <- obs.females[,seq(5,x$n.lens.comm+4)]-pred.females
	res.males <- obs.males[,seq(5,x$n.lens.comm+4)]-pred.males
	
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	lens.m <- len.mat(lens,summeryears.m)
	years.m <- cllyear.mat(lens,summeryears.m)
	
	cex <- 0.6
	ylim <- c(-0.3,0.3)
	
	plot(lens.f,res.females,cex=cex,ylim=ylim,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
	
	plot(years.f,res.females,cex=cex,ylim=ylim,xlim=c(1960,2003),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)

	plot(lens.m,res.males,cex=cex,ylim=ylim,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
	mtext(side=1,outer=F,text="Length (cm)",cex=0.95,line=2.5)
	
	plot(years.m,res.males,cex=cex,ylim=ylim,xlim=c(1960,2003),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=1,outer=F,text="Year",cex=0.95,line=2.5)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in LL",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in OT",cex=1)}

	}

plot.res.comm2<-function(x,series.code) 
	{
	
	obs.s<-x$obs.catlen.comm[x$obs.catlen.comm[,2]==series.code,]
	pred.s<-x$pred.catlen.comm[x$obs.catlen.comm[,2]==series.code,]
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
	lens <- x$lens.comm
#	lens <- seq(10,278,3)
	len.mat <- function(lens,YEARS)
	{
	#
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	#a <- seq(1,ages)
	#years <- seq(YEARS[1],(YEARS[2])-1)
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
#	summeryears.s <- unique(obs.s[,3])
	
	obs.females <- obs.s[obs.s[,4]==1,]
	summeryears.f <- obs.females[,3]
	pred.females <- pred.s[obs.s[,4]==1,]
	
	res.females <- obs.females[,seq(5,x$n.lens.comm+4)]-pred.females
	
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	
	cex <- 0.6
	ylim <- c(-0.3,0.3)
	
print(lens.f)
print(dim(res.females))

	plot(lens.f,res.females,cex=cex,ylim=ylim,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
	
	plot(years.f,res.females,cex=cex,ylim=ylim,xlim=c(1960,2003),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in LL",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in OT",cex=1)}

	}


##########################################################################
#1=summer_rv
#2=spring_rv
#3=fall_rv
#6=4VWcod_rv
#7=fixed station
#8=commercial index

plot.rv(x,1)
#mtext(side=3,outer=T,line=0,text="Proportions at length in summer survey",cex=1)
plot.res(x,1)
#mtext(side=3,outer=T,line=0,text="Residuals (obs.-pred.) for the proportions of summer survey C@L",cex=1)
#	print("OK")

plot.rv(x,7)
#mtext(side=3,outer=T,line=0,text="Proportions at length in fixed station survey",cex=1)
#	print("OK")
plot.res(x,7)
#mtext(side=3,outer=T,line=0,text="Residuals (obs.-pred.) for the proportions of Halibut fixed station survey C@L",cex=1)
#	print("OK")

#//1=LL
#//2=OT

plot.comm(x,1)
plot.res.comm(x,1)
#mtext(side=3,outer=T,line=0,text="Residuals (obs.-pred.) for the proportions of LL C@L",cex=1)
print("OK")
#plot.res.comm2(x,2)
#mtext(side=3,outer=T,line=0,text="Residuals (obs.-pred.) for the proportions of OT C@L",cex=1)
#print("OK")

#######################################################################
#################################################################
#EXPLOITATION

lyears<-seq(1970,x$last.year,1)
#ps.options(horizontal=TRUE)
par(mfrow=c(4,2),omi=c(1,1,1,1),mar=c(3,3,3,.5),las=1,xpd=NA,cex=.9)
plot(lyears,x$u.Canada[11:45],type="l",xlab="",ylab="")
mtext(side=3,text="Canada",cex=0.8,line=0,adj=.1)
plot(lyears,x$u.Foreign[11:45],type="l",xlab="",ylab="")
mtext(side=3,text="Foreign",cex=0.8,line=0,adj=.1)
mtext(side=3,text="Exploitation rates",line=1,outer=T,cex=1.2,las=0)
mtext(outer=T,side=1,text="Years",cex=1.5,line=-23)
mtext(outer=T,side=2,text="Fully exploited U",cex=1.5,line=3,las=0,adj=.8)

plot(lyears,x$u.ll.F[11:45],type="l",xlab="",ylab="")
mtext(side=3,text="LL Females",cex=0.8,line=0,adj=.1)
plot(lyears,x$u.ll.M[11:45],type="l",xlab="",ylab="")
mtext(side=3,text="LL Males",cex=0.8,line=0,adj=.1)

plot(lyears,x$u.ot[11:45],type="l",xlab="",ylab="")
mtext(side=3,text="OT",cex=0.8,line=0,adj=.1)


plot(lyears,x$F.utotal[11:45,8],type="l",xlab="",ylab="")
mtext(side=3,text="LL Females",cex=0.8,line=0,adj=.1)
print("OK")

#########################################################################################
#########################################################################################
# CPUE FITS
par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)
lyears<-seq(1970,x$last.year,1)

plot(lyears,x$pred.summer.rv.total,type="l",xlab="",ylab="",ylim=c(0,max(c(x$summer.rv.total,x$pred.summer.rv.total))*1.2))
points(x$summer.rv.years,x$summer.rv.total,pch=16)
mtext(side=3,text="RV summer",cex=1.5,line=.5,adj=.1)

lyears<-seq(1998,x$last.year,1)
plot(lyears,x$pred.halibut.fixed.total,type="l",xlab="",ylab="",ylim=c(0,max(c(x$halibut.fixed.total,x$pred.halibut.fixed.total))*1.2))
points(x$halibut.fixed.years,x$halibut.fixed.total,pch=16)
mtext(side=3,text="Halibut survey",cex=1.5,line=.5,adj=.1)

#mtext(side=3,text="Survey Index",line=1,outer=T,cex=1.2)
mtext(outer=T,side=1,text="Years",cex=1,line=1)
mtext(outer=F,side=2,text="Abundance index",cex=1.5,line=3,las=0)
print("OK")

###########################################################################
#numbers @ age
par(mfcol=c(1,1),las=1,omi=c(.25,.25,.25,.25),mar=c(5,5,4,4),cex=.9)
plot(c(x$first.year,x$last.year+1),c(x$first.age,x$last.age),type="n",xlab="Year",ylab="Age") 
x2<-x$Females
pyears<-x$first.year:(x$last.year+1)
pages<-x$first.age:x$last.age
prob.vec<-as.vector(x2)
xx<-rep(pyears,length(pages))
yy<-sort(rep(pages,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.2)  
#wireframe(prob.vec ~ yy * xx)
print("OK symbols")

mtext(side=3,text="Numbers at Age: Females",adj=.1,line=1)  
mtext(side=3,text=paste("max point size = ",max(x$Females)),adj=.9,line=0,cex=.9)  

print("OK")

##################################################3
#Exploitation at age

par(mfcol=c(1,1),las=1,omi=c(.25,.25,.25,.25),mar=c(5,5,4,4),cex=.9)
plot(c(x$first.year,x$last.year),c(x$first.age,x$last.age),type="n",xlab="Year",ylab="Age") 
x2<-x$F.utotal[11:((x$last.year-1960)+1),]	# 1970 on
pyears<-x$first.year:(x$last.year)
pages<-x$first.age:(x$last.age)
prob.vec<-as.vector(x2)
xx<-rep(pyears,length(pages))
yy<-sort(rep(pages,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.2)  
mtext(side=3,text="Exploitation Rates: Females",adj=.1,line=1)  
mtext(side=3,text=paste("max point size = ",max(x$F.utotal)),adj=.9,line=0,cex=.9)  
print("OK symbols")


plot(c(x$first.year,x$last.year),c(x$first.age,x$last.age),type="n",xlab="Year",ylab="Age") 
x2<-x$M.utotal[11:((x$last.year-1960)+1),]
pyears<-x$first.year:x$last.year
pages<-x$first.age:(x$last.age)
prob.vec<-as.vector(x2)
xx<-rep(pyears,length(pages))
yy<-sort(rep(pages,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.2)  
print("OK symbols")

mtext(side=3,"Exploitation Rates: Males",adj=.1,outer=F,line=1)  
mtext(side=3,text=paste("max point size = ",max(x$M.utotal)),adj=.9,line=0,cex=.9)  
print("OK")

###########################


dev.off()

}

