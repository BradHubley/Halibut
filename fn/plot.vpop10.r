plot.vpop10<-function(update.x=T)
{

if(update.x==T)
 {
 source("fn/read.admb.r")	 
 vpop10<-read.admb("admb/vpop10")
 vpop10.std<-read.table("admb/vpop10.std",skip=1,row.names=NULL,
		 col.names=c("number","name","value","std.dev"))
 }
#vpop10<<-vpop10
#vpop10std<<-vpop10.std
x<-vpop10
z<-vpop10.std
names(x)<-gsub("_",".",names(x))

#browser()

minsize=.1 #used for bubble plots
maxsize=3 

postscript("plot.vpop10.ps",horizontal=FALSE)
par(mfrow=c(4,2),omi=c(1,1,1,1),mar=c(2,3,2,.5),las=1,xpd=NA,cex=.9)
 
#LANDINGS
lyears<-seq(1960,x$last.year,1)
plot(lyears,x$Canada,type="l",xlab="",ylab="")
mtext(side=3,text="Canada",cex=0.8,line=0,adj=.1)
#mtext(outer=F,side=2,text="Tones",cex=1,line=3.5,las=0)
plot(lyears,x$Foreign,type="l",xlab="",ylab="")
mtext(side=3,text="Foreign",cex=0.8,line=0,adj=.1)
mtext(side=3,text="Landings",line=0,outer=T,cex=1.5,las=0)
mtext(outer=T,side=1,text="Years",cex=1.5,line=-23)

plot(lyears,(x$C.ll.F+x$C.ll.M),type="l",xlab="",ylab="")
mtext(side=3,text="LL ",cex=0.8,line=0,adj=.1)
mtext(outer=T,side=2,text="Tones",cex=1.5,line=.5,las=0,adj=.75)
plot(lyears,x$C.ot,type="l",xlab="",ylab="")
mtext(side=3,text="OT",cex=0.8,line=0,adj=.1)
print("OK LANDINGS")

#Proportion of LANDINGS
par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)

lyears<-seq(1960,x$last.year,1)
plot(lyears,x$Foreign/(x$Foreign+x$Canada),type="l",xlab="",ylab="")
mtext(side=2,text="Proportion foreign landings",line=3,outer=F,cex=1.2,las=0)
#mtext(outer=F,side=1,text="Years",cex=1.5,line=3)

plot(lyears,(x$C.ll.M+x$C.ll.F)/(x$C.ot+x$C.ll.F+x$C.ll.M),type="l",xlab="",ylab="")
mtext(side=2,text="Proportion longline",line=3,outer=F,cex=1.2,las=0)
mtext(outer=F,side=1,text="Years",cex=1.5,line=3)



#SURVEY SELECTIVITIES
par(mfrow=c(4,2),omi=c(1,1,1,1),mar=c(3,3,3,.5),las=1,xpd=NA,cex=.9)

plot(x$ages,x$halibut.fixed.F.sel,type="n",xlab="",ylab="", ylim=c(0,1))
lines(x$ages,x$halibut.fixed.F.sel,lty=1)
lines(x$ages,x$halibut.fixed.M.sel,lty=2)
title(main="Halibut survey station")

plot(x$ages,x$old.summer.rv.F.sel,type="n",xlab="",ylab="", ylim=c(0,1))
lines(x$ages,x$old.summer.rv.F.sel,lty=1)
title(main="Summer RV")

mtext(outer=T,side=1,text="Age",cex=1.5,line=-23)
mtext(outer=T,side=2,text="Selectivity",cex=1.5,las=0,adj=.8)
#mtext(outer=T,side=3,text="Survey Selectivities",cex=1.4,line=2)
 
#COMMERCIAL SELECTIVITIES
#par(mfrow=c(4,2),omi=c(1,1,1,1),mar=c(3,3,3,.5),las=1,xpd=NA)
 
plot(x$ages,x$ll.F.sel,type="n",xlab="",ylab="", ylim=c(0,1))
lines(x$ages,x$ll.F.sel,lty=1)
lines(x$ages,x$ll.M.sel,lty=2)
title(main="Longline")

plot(x$ages,x$ot.F.sel,type="n",xlab="",ylab="", ylim=c(0,1))
lines(x$ages,x$ot.F.sel,lty=1)
title(main="Otter trawl")
legend(13,-0.3,legend=c("females","males"), lty=c(1,2),bty="n",cex=1.1)

############################################################################################
###########################################################################################
#Commercial length comp
par(las=1,mfrow=c(2,1),omi=c(1.5,1,1,1),mar=c(4,4,2,1),cex=.9)

#LL
plot(c(1988,x$last.year),c(min(x$lens.comm),max(x$lens.comm)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$out.obs.catlen.comm[x$out.obs.catlen.comm[,2]==1&x$out.obs.catlen.comm[,4]==1,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]

prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens.com))
yy<-sort(rep(x$lens.com,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  
print("OK symbols")
mtext(side=3,text="Females",adj=.1,line=0)  
mtext(side=3,text="LL Proportions at Length",adj=.1,outer=T)  
#mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9)  
#print("OK")

plot(c(1988,x$last.year),c(min(x$lens.comm),max(x$lens.comm)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$out.obs.catlen.comm[x$out.obs.catlen.comm[,2]==1&x$out.obs.catlen.comm[,4]==2,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]

prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens.com))
yy<-sort(rep(x$lens.com,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  

mtext(side=3,text="Males",adj=.1,line=0)  
#mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9)  
#print("OK")

#OT
plot(c(1977,x$last.year),c(min(x$lens.comm),max(x$lens.comm)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$out.obs.catlen.comm[x$out.obs.catlen.comm[,2]==2,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]

prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens.com))
yy<-sort(rep(x$lens.com,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  

mtext(side=3,text="Females & Males",adj=.1,line=0)  
mtext(side=3,text="OT Proportions at Length",adj=.1,outer=T)  
mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9)  
print("OK Commercial length comp")

################################################################################################
#SURVEY PROPORTIONS:
#1=summer_rv
#2=spring_rv
#3=fall_rv
#6=4VWcod_rv

par(las=1,mfrow=c(2,1),omi=c(1.5,1,1,1),mar=c(4,4,2,1),cex=.9)

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

mtext(side=3,text="Females",adj=.1,line=0)  
mtext(side=3,text="Summer Survey Proportions at Length",adj=.1,outer=T)  
mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9)  
########

plot(c(x$first.year,x$last.year),c(min(x$lens),max(x$lens)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$obs.catlen.rv[x$obs.catlen.rv[,2]==1&x$obs.catlen.rv[,4]==2,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]
prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens))
yy<-sort(rep(x$lens,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  
print("OK RV length comp")

mtext(side=3,text="Males",adj=.1,line=0)  
mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9) 

##############################################################################3
#halibut survey
plot(c(1997,x$last.year),c(min(x$lens),max(x$lens)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$obs.catlen.rv[x$obs.catlen.rv[,2]==7&x$obs.catlen.rv[,4]==1,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]

prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens))
yy<-sort(rep(x$lens,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  

mtext(side=3,text="Females",adj=.1,line=0)  
mtext(side=3,text="Halibut Survey Proportions at Length",adj=.1,outer=T)  
mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9)  
########
plot(c(1997,x$last.year),c(min(x$lens),max(x$lens)),type="n",xlab="Year",ylab="Length (cm)") 
x1<-x$obs.catlen.rv[x$obs.catlen.rv[,2]==7&x$obs.catlen.rv[,4]==2,]
x1[x1==0]<-NA
x2<-x1[,-c(1,2,3,4)]
pyears<-x1[,3]
prob.vec<-as.vector(x2)
xx<-rep(pyears,length(x$lens))
yy<-sort(rep(x$lens,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.15)  
print("OK HalSurv length comp")

mtext(side=3,text="Males",adj=.1,line=0)  
mtext(side=3,text=paste("max point size = ",max(x2,na.rm=T)),adj=.9,line=0,cex=.9) 


##############################################################################
#PLOT MEAN LENGTH IN CATCH
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(3,3,3,.5),las=1,xpd=NA,cex=.9)

  len1<-data.frame(x$out.obs.catlen.comm[,1:4],(x$obs.catlen.comm[,5:80]%*%x$lens.comm)/rowSums(x$obs.catlen.comm[,5:80]))
  names(len1)=c("n","method","year","sex","mean.len")
  
  len1.ll.f<-len1[len1$method==1&len1$sex==1,]
  len1.ll.m<-len1[len1$method==1&len1$sex==2,]
  len1.ot<-len1[len1$method==2,]
  
  llyears<-1988:x$last.year
  otyears<-1984:x$last.year
  ymax<-max(len1.ll.f[,5],len1.ll.m[,5],len1.ot[,5])
  
plot(llyears,len1.ll.f[,5],type="l",xlim=c(1980,x$last.year),ylim=c(50,150),xlab="",ylab="",col='red',lwd=2)
lines(llyears,len1.ll.m[,5],col='blue',lwd=2)
lines(otyears,len1.ot[,5],col='black',lwd=2,lty=2)
mtext(outer=F,side=2,text="Mean length (cm)",cex=1.5,line=3,las=0,adj=.5)
legend(1993,150,legend=c("LL female","LL male","OT"),lty=c(1,1,2),col=c('red','blue','black'),lwd=2,bty="n",cex=1.1)

##############################################################################
#PLOT MEAN WEIGHT IN CATCH
# a<-0.009112
# b<-3.062269
# x$obs.catlen.comm[,5:80]%*%(a*x$lens.comm^b)/rowSums(x$obs.catlen.comm[,5:80])

##############################################################################
# SPAWNER NUMBERS, BIOMASS AND NUMBER OF RECRUITS

# years 
years<-c(seq(x$first.year,x$last.year,1),NA)
lyears<-c(seq(x$first.year+1,x$last.year,1),NA)

par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)
R<-x$rec[-1]	#recruits  drop 1970

ages<-seq(x$first.age,x$last.age,1)

plot(years,x$SSN/1000,type="n",xlab="",ylab="",ylim=c(0,max(x$rec)*1.5/1000),axes=F)
lines(years,x$SSN/1000,lty=1,lwd=2)	#spawner numbers in 1000s 
lines(lyears,R/1000,lty=4,lwd=1)  		#recruits in 1000s
axis(1)
axis(2)

par(new=T)
plot(years,x$SSB/1000,type="n",xlab="",ylab="",ylim=c(0,max(x$SSB/1000)),axes=F)
lines(years,x$SSB/1000,lty=2,lwd=2)		# biomass in 1000 tons
axis(4)
box()
#mtext(side=3,text="Spawner numbers and biomass (Mature females), and recruits",cex=1.5,line=1)
mtext(side=2,text="S and R (1000s)",cex=1.5,line=3,las=0)
mtext(side=3,text="Females only",cex=1.5,line=.5,las=0)
mtext(side=4,text="SSB (1000 t)",cex=1.5,line=3,las=0)
#legend(1993,max(x$SSB[-length(x$SSB)]/1000)*.9,legend=c("Numbers","Recruits","Biomass"),lty=c(1,4,2),lwd=2,bty="n",cex=1.1)
legend(1993,-max(x$SSB[-length(x$SSB)]/1000)/4,legend=rev(c("Numbers","Recruits","Biomass")),lty=rev(c(1,4,2)),lwd=2,bty="n",cex=1.1)

print("SSB")

##############################################################################
# SPAWNER NUMBERS, BIOMASS AND NUMBER OF RECRUITS

par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)

R<-x$rec[-1]	#recruits  drop 1970
ages<-seq(x$first.age,x$last.age,1)
plot(years,x$SSN/1000,type="n",xlab="",ylab="",ylim=c(0,max(x$rec[-length(x$rec)])*1.5/1000),axes=F)
#lines(years,x$SSN[-length(x$SSN)]/1000,lty=1,lwd=2)	#spawner numbers in 1000s 
lines(lyears,R/1000,lty=4,lwd=1)  		#recruits in 1000s
axis(1)
axis(2)

par(new=T)
plot(years,x$SSB2/1000,type="n",xlab="",ylab="",ylim=c(0,max(x$SSB2/1000)),axes=F)
lines(years,x$SSB2/1000,lty=2,lwd=2)		# biomass in 1000 tons
axis(4)
box()
#mtext(side=3,text="Spawner numbers and biomass (Mature females), and recruits",cex=1.5,line=1)
mtext(side=2,text="Recruits (1000s)",cex=1.5,line=3,las=0)
mtext(side=4,text="SSB (1000 t)",cex=1.5,line=3,las=0)
mtext(side=3,text="A)",cex=1.5,line=.5,adj=0)
#legend(1993,max(x$SSB2/1000)*.9,legend=rev(c("Numbers","Recruits","Biomass")),lty=c(1,4,2),lwd=2,bty="n",cex=1.1)
#legend(1993,8.5,legend=rev(c("Recruits","Biomass")),lty=rev(c(4,2)),lwd=2,bty="n",cex=1.1)

print("SSB2")
print(x$SSB2/1000)

##############################################################################
# TOTAL NUMBERS AND BIOMASS
plot(years,x$Total.N/1000,type="n",xlab="",ylab="",ylim=c(0,max(x$Total.N)/1000),axes=F)
lines(years,x$Total.N/1000,lty=1,lwd=2)	#total numbers in millions
axis(1)
axis(2)
par(new=T)
plot(years,x$Total.W/1000,type="n",xlab="",ylab="",ylim=c(0,max((x$Total.W/1000)[-1])),axes=F)
lines(years[-1],(x$Total.W/1000)[-1],lty=2,lwd=2)		#total biomass in 1000 t
axis(4)
box()
mtext(side=3,text="B)",cex=1.5,line=.5,adj=0)
mtext(side=2,text="N (1000s)",cex=1.5,line=3,las=0)
mtext(side=4,text="Total biomass (1000 t)",cex=1.5,line=3,las=0)
legend(1993,max(x$Total.W/1000)+max(x$Total.W/1000)/3,legend=rev(c("Numbers","Recruits","Biomass")),lty=rev(c(1,4,2)),lwd=2,bty="n",cex=1.1)
print("Total N & B")

############################################################################################
###########################################################################
#mean weights

# years 
years<-seq(x$first.year,x$last.year,1)

par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)

#w.f<-apply(x$Females*x$w[1,],1,sum)/rowSums(x$Females)
#w.m<-apply(x$Males*x$w[2,],1,sum)/rowSums(x$Males)
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
mtext(side=1,text="Year",cex=1.5,line=2.5)
mtext(side=2,text="Mean weight (kg)",cex=1.5,line=3.5,las=0)
legend(1998,max(w.f)*.9,legend=c("Females","Males"),lty=c(1,2),lwd=2,bty="n",cex=1.1)
print("mean weights")

###########################################################################
#mean exploitation
u.f<-apply(x$Females[1:dim(x$Females)[1]-1,]*x$F.utotal[11:dim(x$F.utotal)[1],],1,sum)/rowSums(x$Males)[-dim(x$Females)[1]]
u.m<-apply(x$Males[1:dim(x$Males)[1]-1,]*x$M.utotal[11:dim(x$M.utotal)[1],],1,sum)/rowSums(x$Males)[-dim(x$Males)[1]]

u.f<-NULL
u.m<-NULL

for(y in 1:length(years))
  {
  u.f[y]<-x$Females[y,4:x$last.age]%*%x$F.utotal[y+10,4:x$last.age]/sum(x$Females[y,4:x$last.age])
  u.m[y]<-x$Males[y,4:x$last.age]%*%x$M.utotal[y+10,4:x$last.age]/sum(x$Males[y,4:x$last.age])
  }
  
#drop last year
#u.f<-u.f[-length(u.f)]
#u.m<-u.m[-length(u.m)]

plot(years,u.f,type="n",xlab="",ylab="",ylim=c(0,max(u.f)),axes=F)
axis(1)
axis(2)
lines(years,u.f,lty=1,lwd=2)	
lines(years,u.m,lty=2,lwd=2)	
box()
mtext(side=1,text="Year",cex=1.5,line=2.5)
mtext(side=2,text="Exploitation rate",cex=1.5,line=3.5,las=0)
legend(1971,max(u.f)*.9,legend=c("Females","Males"),lty=c(1,2),lwd=2,bty="n",cex=1.1)

print("Exploitation rate")
##############################################################################
# SPAWNER-RECRUIT RELATIONSHIP
## see plot.rk1.r

#S<-seq(min(x$SSB),max(x$SSB),length=100)
S<-seq(0,50,length=1000)
#R<-2.71828*x$p2/x$p1*S*exp(-S/x$p1);	# Ricker reparamertized   
R<-2.71828*x$p2/x$p1*(S-x$p3)*exp(-(S-x$p3)/x$p1);	# Ricker reparamertized   with depensation; x intercept offset estiamted
#B-H
R2<-(x$alpha*(S-x$p3)/(1+x$alpha*(S-x$p3)/x$R0))
R<-R/1000
x$rec<-x$rec/1000

#plot 1970-2006 SSB vs 1971-2007 recruitment
plot(x$SSB2[1:(length(x$SSB2)-2)]/1000,x$rec[-length(x$rec)][-c(1,length(x$rec[-length(x$rec)]))],type='n',xlab="",ylab="", ylim=c(0,max(x$rec[-length(x$rec)])),xlim=c(min(x$SSB/1000),max(x$SSB/1000)))
points(x$SSB2[1:(length(x$SSB2)-2)]/1000,x$rec[-length(x$rec)][-c(1,length(x$rec[-length(x$rec)]))],pch=as.character(substring(lyears,3,3)),cex=.7)
lines(S,R,lty=1,lwd=2,xpd=F)  #Ricker
#lines(S,R2,lty=1,lwd=2,xpd=F)	#B-H
mtext(side=3,text="Spawner-recruitment relationship",cex=1.5,line=1)
mtext(side=1,text="SSB (1000 t)",cex=1.5,line=2.5)
mtext(side=2,text="Age-1 recruits (1000s)",cex=1.5,line=3.5,las=0)


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
  if(series.code==7){mtext(side=3,outer=T,line=-.5,text="Proportions at length in halibut survey",cex=1)}
  if(series.code==8){mtext(side=3,outer=T,line=-.5,text="Proportions at length in commercial index",cex=1)}
   }
	}
#1=summer_rv
#2=spring_rv
#3=fall_rv
#6=4VWcod_rv
#7=halibut survey
#8=commercial index

##########################################################################
#####################################################################
#RV RAW residuals
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
	
	plot(years.f,res.females,cex=cex,ylim=ylim,xlim=c(1970,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)

	plot(lens.m,res.males,cex=cex,ylim=ylim,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
	mtext(side=1,outer=F,text="Length (cm)",cex=0.95,line=2.5)
	
	plot(years.m,res.males,cex=cex,ylim=ylim,xlim=c(1970,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=1,outer=F,text="Year",cex=0.95,line=2.5)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in summer survey",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in spring survey",cex=1)}
  if(series.code==3){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in fall survey",cex=1)}
  if(series.code==6){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in 4VsW cod survey",cex=1)}
  if(series.code==7){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in halibut survey",cex=1)}
  if(series.code==8){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in commercial index",cex=1)}

	}
#1=summer_rv
#2=spring_rv
#3=fall_rv
#6=4VWcod_rv

##########################################################################
#####################################################################
#RV standardized residuals
#print("OK")

plot.res.standardized<-function(x,series.code) 
	{
	
    resid.s<-x$catlen_rv_resid[x$catlen_rv_resid[,3]==series.code,]
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,1.5,.5),cex=.9)
	lens <- x$lens
	len.mat <- function(lens,YEARS)
	{

	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
	obs.females <- resid.s[resid.s[,4+1]==1,]
	summeryears.f <- obs.females[,3+1]

	obs.males <- resid.s[resid.s[,4+1]==2,]
	summeryears.m <- obs.males[,3+1]

	res.females <- obs.females[,seq(6,x$n.lens+5)]
	res.males <- obs.males[,seq(6,x$n.lens+5)]
	
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	lens.m <- len.mat(lens,summeryears.m)
	years.m <- cllyear.mat(lens,summeryears.m)

	cex <- 0.6
	ylim.f <- c(min(res.females),max(res.females))
	ylim.m <- c(min(res.males),max(res.males))

	plot(lens.f,res.females,cex=cex,ylim=ylim.f,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
	
	plot(years.f,res.females,cex=cex,ylim=ylim.f,xlim=c(1970,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)

	plot(lens.m,res.males,cex=cex,ylim=ylim.m,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
	mtext(side=1,outer=F,text="Length (cm)",cex=0.95,line=2.5)
	
	plot(years.m,res.males,cex=cex,ylim=ylim.m,xlim=c(1970,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=1,outer=F,text="Year",cex=0.95,line=2.5)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in summer survey",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in spring survey",cex=1)}
  if(series.code==3){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in fall survey",cex=1)}
  if(series.code==6){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in 4VsW cod survey",cex=1)}
  if(series.code==7){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in halibut survey",cex=1)}
  if(series.code==8){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in commercial index",cex=1)}

#  mtext(side=3,outer=T,line=1,text="Standardized",cex=1.2)
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

	obs.s<-x$out.obs.catlen.comm[x$out.obs.catlen.comm[,2]==series.code,]
	pred.s<-x$pred.catlen.comm[x$out.obs.catlen.comm[,2]==series.code,]
	
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

##########################################################################
plot.comm2<-function(x,series.code) 
	{

	obs.s<-x$out.obs.catlen.comm[x$out.obs.catlen.comm[,2]==series.code,]
	pred.s<-x$pred.catlen.comm[x$out.obs.catlen.comm[,2]==series.code,]
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
	lens <- x$lens.comm
	len.mat <- function(lens,YEARS)
	{
	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
	obs.females <- obs.s[obs.s[,4]==1,]
	summeryears.f <- obs.females[,3]
	obs.females<-obs.females[,seq(5,x$n.lens.comm+4)]

	pred.females <- pred.s[obs.s[,4]==1,]
		
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	lens.m <- len.mat(lens,summeryears.m)
	years.m <- cllyear.mat(lens,summeryears.m)
	
	cex <- 0.6
	ylim <- c(-0.3,0.3)
	
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
#RAW RESIDUALS CATCH@LENGTH
par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
#print("OK")

plot.res.comm<-function(x,series.code) 
	{
	
	obs.s<-x$out.obs.catlen.comm[x$out.obs.catlen.comm[,2]==series.code,]
	pred.s<-x$pred.catlen.comm[x$out.obs.catlen.comm[,2]==series.code,]
	
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
	
	plot(years.f,res.females,cex=cex,ylim=ylim,xlim=c(1985,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)

	plot(lens.m,res.males,cex=cex,ylim=ylim,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
	mtext(side=1,outer=F,text="Length (cm)",cex=0.95,line=2.5)
	
	plot(years.m,res.males,cex=cex,ylim=ylim,xlim=c(1985,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=1,outer=F,text="Year",cex=0.95,line=2.5)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in LL",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in OT",cex=1)}

	}

plot.res.comm2<-function(x,series.code) 
	{
	
	obs.s<-x$out.obs.catlen.comm[x$out.obs.catlen.comm[,2]==series.code,]
	pred.s<-x$pred.catlen.comm[x$out.obs.catlen.comm[,2]==series.code,]
	
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
	
	plot(years.f,res.females,cex=cex,ylim=ylim,xlim=c(1985,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in LL",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=-.5,text="Residual proportions at length in OT",cex=1)}

	}


##########################################################################################
#Standardized RESIDUALS CATCH@LENGTH
par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,.75,.5),cex=.9)
#print("OK")

plot.res.comm.standardized<-function(x,series.code) 
	{
	
	obs.s<-x$catlen_comm_resid[x$catlen_comm_resid[,3]==series.code,]
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,1.5,.5),cex=.9)
	lens <- x$lens.comm
	len.mat <- function(lens,YEARS)
	{
	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
	obs.females <- obs.s[obs.s[,4+1]==1,]
	summeryears.f <- obs.females[,3+1]
	obs.males <- obs.s[obs.s[,4+1]==2,]
	summeryears.m <- obs.males[,3+1]
	
	res.females <- obs.females[,seq(6,x$n.lens.comm+5)]
	res.males <- obs.males[,seq(6,x$n.lens.comm+5)]
	
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	lens.m <- len.mat(lens,summeryears.m)
	years.m <- cllyear.mat(lens,summeryears.m)
	
	cex <- 0.6
	ylim.f <- c(min(res.females),max(res.females))
	ylim.m <- c(min(res.males),max(res.males))

	plot(lens.f,res.females,cex=cex,ylim=ylim.f,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
	
	plot(years.f,res.females,cex=cex,ylim=ylim.f,xlim=c(1985,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)

	plot(lens.m,res.males,cex=cex,ylim=ylim.m,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
	mtext(side=1,outer=F,text="Length (cm)",cex=0.95,line=2.5)
	
	plot(years.m,res.males,cex=cex,ylim=ylim.m,xlim=c(1985,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=1,outer=F,text="Year",cex=0.95,line=2.5)
	mtext(side=3,outer=F,"Males",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in LL",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in OT",cex=1)}

	}

plot.res.comm2.standardized<-function(x,series.code) 
	{
	
	obs.s<-x$catlen_comm_resid[x$catlen_comm_resid[,3]==series.code,]
	
	par(mfrow=c(4,2),mar=c(3,4,1,1),omi=c(.75,1,1.5,.5),cex=.9)
	lens <- x$lens.comm
	len.mat <- function(lens,YEARS)
	{
	b <- rep(lens,length(YEARS))
	bb <- t(matrix(b,length(lens),length(YEARS)))
	return(bb)
	}
	
	cllyear.mat <- function(ages,YEARS)
	{
	b <- rep(YEARS,length(ages))
	bb <- matrix(b,length(YEARS),length(ages))
	return(bb)
	}
	
	obs.females <- obs.s[obs.s[,4+1]==1,]
	summeryears.f <- obs.females[,3+1]

	res.females <- obs.females[,seq(6,x$n.lens.comm+5)]
	
	lens.f <- len.mat(lens,summeryears.f)
	years.f <- cllyear.mat(lens,summeryears.f)
	
	cex <- 0.6
	ylim.f <- c(min(res.females),max(res.females))
	
	plot(lens.f,res.females,cex=cex,ylim=ylim.f,ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
	
	plot(years.f,res.females,cex=cex,ylim=ylim.f,xlim=c(1985,x$last.year),ylab="",xlab="")
	abline(h=0,lty=2,xpd=F)
	mtext(side=3,outer=F,"Females",line=0,cex=.8)
  if(series.code==1){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in LL",cex=1)}
  if(series.code==2){mtext(side=3,outer=T,line=0,text="Standardized residual proportions at length in OT",cex=1)}

	}


##########################################################################
#1=summer_rv
#2=spring_rv
#3=fall_rv
#6=4VWcod_rv
#7=halibut survey
#8=commercial index

plot.rv(x,1)
plot.res(x,1)
plot.res.standardized(x,1)
#	print("OK")

plot.rv(x,7)
plot.res(x,7)
plot.res.standardized(x,7)
#	print("OK")

#//1=LL
#//2=OT

plot.comm(x,1)
plot.res.comm(x,1)
plot.res.comm.standardized(x,1)
print("ok")

#not working
#plot.comm2(x,2)
#print("OK")
#plot.res.comm2(x,2)
#print("OK")
#plot.res.comm2.standardized(x,2)
#print("OK")

#######################################################################
#################################################################
#Fishing Mortality

lyears<-seq(1970,x$last.year,1)
y1<-length(1960:2009)
#ps.options(horizontal=TRUE)
par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)

#plot(lyears,x$F.utotal[11:y1,8],type="l",xlab="",ylab="",lty=1,lwd=2)
#lines(lyears,x$M.utotal[11:y1,8],lty=2,lwd=2)

plot(lyears,x$u.ll.F,type="l",xlab="",ylab="",ylim=c(0,1))
lines(lyears,x$u.ll.M,lty=2,lwd=2)
lines(lyears,x$u.ot[11:y1],type="l",lty=1,lwd=2,col='blue')
mtext(outer=F,side=2,text="Fully recruited fishing mortality",cex=1.5,line=3,las=0,adj=.5)
legend(1970,.9,legend=(c("LL female","LL male","OT")),col=(c("black","black","blue")),lwd=2,lty=c(1,2,1),bty="n",cex=.9)

print("F")

#########################################################################################
#################################################################
#Fishing Mortality
# compare with tagging
lyears<-seq(1970,x$last.year,1)
par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9,xpd=F)
tag.F<-c(0.169,0.249,0.197,0.096)
tag.F[4]<-tag.F[4]*2
tag.F.se<-0.04
upper.tag.F<-tag.F+2*tag.F.se
lower.tag.F<-tag.F-2*tag.F.se
tag.y<-2007:2010

u.ll.F.se<-z$std.dev[z$name=="u_ll_F"]

#don't do it this way
#cf=x$u.ll.F/(x$u.ll.F+.1)*x$Females[-41,8]*(1-exp(-(x$u.ll.F+.1)))
#cm=x$u.ll.M/(x$u.ll.M+.1)*x$Males[-41,8]*(1-exp(-(x$u.ll.M+.1)))
#u=(cf+cm)/(x$Females[-41,8]+x$Males[-41,8])
#F=(u*.3/(1-exp(-.3)))
F.bar=(rowSums(x$F.utotal[11:50,8:20]*x$Females[1:40,8:20])+rowSums(x$M.utotal[11:50,8:20]*x$Males[1:40,8:20]))/rowSums(x$Females[1:40,8:20]+x$Males[1:40,8:20])

plot(lyears,F.bar,type="l",xlab="",ylab="",ylim=c(0.05,.35),xlim=c(2005,2010))
lines(lyears,F.bar+u.ll.F.se*2,lty=8,lwd=2)
lines(lyears,F.bar-u.ll.F.se*2,lty=8,lwd=2)
points(tag.y,tag.F,pch=1,col='red')	#2010 estimate  0.04
arrows(tag.y,upper.tag.F,tag.y,lower.tag.F,length=0.05,angle=90,code=3,col='red')  # change length to 0.05 if want horizontal lines on error bar and code=3
mtext(outer=F,side=2,text="Fully exploited fishing mortality",cex=1.5,line=3,las=0,adj=.5)
#legend(1970,.9,legend=(c("LL female","LL male","OT")),col=(c("black","black","blue")),lwd=2,lty=c(1,2,1),bty="n",cex=.9)

print("F")

#########################################################################################
#sex ratio in the catch
s.rat<-rowSums(x$obs.catlen.comm[x$obs.catlen.comm[,2]==1&x$obs.catlen.comm[,4]==1,])/
(rowSums(x$obs.catlen.comm[x$obs.catlen.comm[,2]==1&x$obs.catlen.comm[,4]==1,])+
rowSums(x$obs.catlen.comm[x$obs.catlen.comm[,2]==1&x$obs.catlen.comm[,4]==2,]))


#########################################################################################
# CPUE FITS
par(mfcol=c(2,1),las=1,omi=c(1,.75,.5,.5),mar=c(5,5,4,4),cex=.9)
lyears<-seq(1970,x$last.year,1)

#ylim1<-c(0,max(c(x$summer.rv.total/10^6,x$pred.summer.rv.total/10^6))*1.2)
ylim1<-c(0,3)
plot(lyears,x$pred.summer.rv.total/10^6,type="l",axes=F,xlab="",ylab="",ylim=ylim1)
axis(1,seq(1970,2010,5))
axis(2)
points(x$summer.rv.years,x$summer.rv.total/10^6,pch=16)
points(2010,2941051/10^6,pch=3)	#2010 estimate  
mtext(side=3,text="Summer RV",cex=1.5,line=0,adj=.1)

lyears<-seq(1998,x$last.year,1)
plot(lyears,x$pred.halibut.fixed.total/10^6,type="l",axes=F,xlab="",ylab="",ylim=c(0,max(c(x$halibut.fixed.total/10^6,x$pred.halibut.fixed.total/10^6))*1.2))
axis(1,seq(1998,2010,2))
axis(2)
points(x$halibut.fixed.years,x$halibut.fixed.total/10^6,pch=16)
points(2010,55.424* 28701/10^6,pch=3)	#2010 estimate rescaled
mtext(side=3,text="Halibut survey",cex=1.5,line=0,adj=.1)
#mtext(outer=F,side=2,text="Standardized index of abundance (kg/1000 hooks/10 hrs)",cex=1.5,line=3,las=0)

mtext(outer=T,side=1,text="Years",cex=1.5,line=1)
mtext(outer=T,side=2,text="Index of abundance",cex=1.5,line=-1,las=0)
print("OK")

###########################################################################
#numbers @ age

par(mfcol=c(1,1),las=1,omi=c(1,1,1,1),mar=c(2,3,5,.5),cex=.9)
plot(c(x$first.year,x$last.year+1),c(x$first.age,20),type="n",xlab="Year",ylab="Age",cex=1.5) 
x2<-x$Females[,1:20]
pyears<-x$first.year:(x$last.year+1)
pages<-x$first.age:20
prob.vec<-as.vector(x2)
xx<-rep(pyears,length(pages))
yy<-sort(rep(pages,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.2)  
#wireframe(prob.vec ~ yy * xx)
print("OK symbols")

mtext(side=3,text="Numbers at Age: Females",adj=.1,line=1,cex=1.5)  
#mtext(side=3,text=paste("max point size = ",max(x$Females[,1:20])),adj=.9,line=0,cex=.9)  

print("OK")

###########################################################################
#RV residuals numbers @ length
#Females
par(mfcol=c(1,1),las=1,omi=c(1,1,1,.5),mar=c(2,3,5,2),cex=.9)

#age breaks at len
f.breaks<-(x$lenatage[1,1:19]+x$lenatage[1,2:20])/2
m.breaks<-(x$lenatage[2,1:19]+x$lenatage[2,2:20])/2
f.lab<-x$lenatage[1,1:15]
f.lab[3]<-f.lab[3]-3
f.lab[4]<-f.lab[4]+5
m.lab<-x$lenatage[2,1:15]
age.lab<-c("age 1","age 2","age 3","age 4","age 5","age 6","age 7","age 8","age 9","age 10","age 11","age 12","age 13","age 14","age 15")

plot(c(x$first.year,x$last.year+1),c(x$lens[1],x$lens[43]),type="n",xlab="Year",ylab="Length (cm)",cex=1.5) 
x2<-x$catlen_rv_resid[1:41,6:77]
x2<-x2[-37,]		# two estimates for 2005
pyears<-x$first.year:(x$last.year)
plens<-x$lens
prob.vec<-(as.vector(x2))
#prob.vec<-rnorm(length(prob.vec),0,1)
temp<-prob.vec
temp2<-prob.vec
temp[temp==0]<-NA
temp[temp>0]<-"red"
temp[temp<0]<-"black"
temp2[temp2>=0]<-temp2[temp2>=0]
temp2[temp2<0]<-temp2[temp2<0]
xx<-rep(pyears,length(plens))
yy<-sort(rep(plens,length(pyears)))
symbols(xx, yy, circles=abs(prob.vec), add = T, inches = 0.2,fg=temp,bg=temp)  
#wireframe(prob.vec ~ yy * xx)
for(i in 1:10)
  {
  abline(h=f.breaks[i],lwd=1,lty=8,col='dark green',xpd=T)
  text(2012,f.lab[i],age.lab[i],cex=.9,adj=0)
  }
print("OK symbols")

mtext(side=3,text="Standardized residual proportions at length: RV Females",adj=.1,line=1,cex=1.5)  

print("OK")

###########################################################################
#RV residuals numbers @ length
#Males
par(mfcol=c(1,1),las=1,omi=c(1,1,1,.5),mar=c(2,3,5,2),cex=.9)

plot(c(x$first.year,x$last.year+1),c(x$lens[1],x$lens[43]),type="n",xlab="Year",ylab="Length (cm)",cex=1.5) 
x2<-x$catlen_rv_resid[42:82,6:77]
x2<-x2[-37,]		# two estimates for 2005
pyears<-x$first.year:(x$last.year)
plens<-x$lens
prob.vec<-(as.vector(x2))
#prob.vec<-rnorm(length(prob.vec),0,1)
temp<-prob.vec
temp2<-prob.vec
temp[temp==0]<-NA
temp[temp>0]<-"red"
temp[temp<0]<-"black"
temp2[temp2>=0]<-temp2[temp2>=0]
temp2[temp2<0]<-temp2[temp2<0]
xx<-rep(pyears,length(plens))
yy<-sort(rep(plens,length(pyears)))
symbols(xx, yy, circles=abs(prob.vec), add = T, inches = 0.2,fg=temp,bg=temp)  
#wireframe(prob.vec ~ yy * xx)
for(i in 1:20)
  {
  abline(h=m.breaks[i],lwd=1,lty=8,col='dark green',xpd=T)
  }
for(j in 1:8)
  {
  text(2012,m.lab[j],age.lab[j],cex=.9,adj=0)
  }
print("OK symbols")

mtext(side=3,text="Standardized residual proportions at length: RV Males",adj=.1,line=1,cex=1.5)  

print("OK")

###########################################################################
#LL residuals numbers @ length
#Females
par(mfcol=c(1,1),las=1,omi=c(1,1,1,.5),mar=c(2,3,5,2),cex=.9)

plot(c(1988,x$last.year+1),c(x$lens[10],x$lens[62]),type="n",xlab="Year",ylab="Length (cm)",cex=1.5) 
x$catlen_comm_resid<-as.data.frame(x$catlen_comm_resid)
x2<-x$catlen_comm_resid[x$catlen_comm_resid[3]==1&x$catlen_comm_resid[5]==1,6:81]
pyears<-1988:(x$last.year)
plens<-x$lens.comm
prob.vec<-(as.vector(as.matrix(x2)))
temp<-prob.vec
temp2<-prob.vec
temp[temp==0]<-NA
temp[temp>=0]<-"red"
temp[temp<0]<-"black"
temp2[temp2>=0]<-temp2[temp2>=0]
temp2[temp2<0]<-temp2[temp2<0]
xx<-rep(pyears,length(plens))
yy<-sort(rep(plens,length(pyears)))
symbols(xx, yy, circles=abs(prob.vec), add = T, inches = 0.2,fg=temp,bg=temp)  
#wireframe(prob.vec ~ yy * xx)
for(i in 3:20)
  {
  abline(h=f.breaks[i],lwd=1,lty=8,col='dark green',xpd=T)
  }
for(j in 3:13)
  {
  text(2011,f.lab[j],age.lab[j],cex=.9,adj=0)
  }
print("OK symbols")

mtext(side=3,text="Standardized residual proportions at length: LL Females",adj=.1,line=1,cex=1.5)  

print("OK")

###########################################################################
#LL residuals numbers @ length
#Males
par(mfcol=c(1,1),las=1,omi=c(1,1,1,.5),mar=c(2,3,5,2),cex=.9)

plot(c(1988,x$last.year+1),c(x$lens[10],x$lens[53]),type="n",xlab="Year",ylab="Length (cm)",cex=1.5) 
x$catlen_comm_resid<-as.data.frame(x$catlen_comm_resid)
x2<-x$catlen_comm_resid[x$catlen_comm_resid[3]==1&x$catlen_comm_resid[5]==2,6:81]
pyears<-1988:(x$last.year)
plens<-x$lens.comm
prob.vec<-(as.vector(as.matrix(x2)))
temp<-prob.vec
temp2<-prob.vec
temp[temp==0]<-NA
temp[temp>=0]<-"red"
temp[temp<0]<-"black"
temp2[temp2>=0]<-temp2[temp2>=0]
temp2[temp2<0]<-temp2[temp2<0]
xx<-rep(pyears,length(plens))
yy<-sort(rep(plens,length(pyears)))
symbols(xx, yy, circles=abs(prob.vec), add = T, inches = 0.2,fg=temp,bg=temp)  
#wireframe(prob.vec ~ yy * xx)
for(i in 3:20)
  {
  abline(h=m.breaks[i],lwd=1,lty=8,col='dark green',xpd=T)
  }
for(j in 3:8)
  {
  text(2011,m.lab[j],age.lab[j],cex=.9,adj=0)
  }
print("OK symbols")

mtext(side=3,text="Standardized residual proportions at length: LL Males",adj=.1,line=1,cex=1.5)  

print("OK")

###########################################################################
#OT residuals numbers @ length
par(mfcol=c(1,1),las=1,omi=c(1,1,1,.5),mar=c(2,3,5,2),cex=.9)

plot(c(1984,x$last.year+1),c(x$lens[1],x$lens[43]),type="n",xlab="Year",ylab="Length (cm)",cex=1.5) 
x$catlen_comm_resid<-as.data.frame(x$catlen_comm_resid)
x2<-x$catlen_comm_resid[x$catlen_comm_resid[3]==2,6:81]
pyears<-1984:(x$last.year)
plens<-x$lens.comm
prob.vec<-(as.vector(as.matrix(x2)))
temp<-prob.vec
temp2<-prob.vec
temp[temp==0]<-NA
temp[temp>=0]<-"red"
temp[temp<0]<-"black"
temp2[temp2>=0]<-temp2[temp2>=0]
temp2[temp2<0]<-temp2[temp2<0]
xx<-rep(pyears,length(plens))
yy<-sort(rep(plens,length(pyears)))
symbols(xx, yy, circles=abs(prob.vec), add = T, inches = 0.2,fg=temp,bg=temp)  
#wireframe(prob.vec ~ yy * xx)
for(i in 1:10)
  {
  abline(h=f.breaks[i],lwd=1,lty=8,col='dark green',xpd=T)
  text(2011,f.lab[i],age.lab[i],cex=.9,adj=0)
  }
print("OK symbols")

mtext(side=3,text="Standardized residual proportions at length: OT",adj=.1,line=1,cex=1.5)  

print("OK")

##################################################
#Exploitation at age

par(mfcol=c(1,1),las=1,omi=c(1,1,1,1),mar=c(2,3,5,.5),cex=.9)
plot(c(x$first.year,x$last.year),c(x$first.age,x$last.age),type="n",xlab="Year",ylab="Age",cex=1.5) 
x2<-x$F.utotal[11:((x$last.year-1960)+1),]	# 1970 on
pyears<-x$first.year:(x$last.year)
pages<-x$first.age:(x$last.age)
prob.vec<-as.vector(x2)
xx<-rep(pyears,length(pages))
yy<-sort(rep(pages,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.2)  
mtext(side=3,text="Exploitation Rates: Females",adj=.1,line=1,cex=1.5)  
#mtext(side=3,text=paste("max point size = ",max(x$F.utotal)),adj=.9,line=0,cex=.9)  
print("OK symbols")


plot(c(x$first.year,x$last.year),c(x$first.age,x$last.age),type="n",xlab="Year",ylab="Age",cex=1.5) 
x2<-x$M.utotal[11:((x$last.year-1960)+1),]
pyears<-x$first.year:x$last.year
pages<-x$first.age:(x$last.age)
prob.vec<-as.vector(x2)
xx<-rep(pyears,length(pages))
yy<-sort(rep(pages,length(pyears)))
symbols(xx, yy, circles=prob.vec, add = T, inches = 0.2)  
print("OK symbols")

mtext(side=3,"Exploitation Rates: Males",adj=.1,outer=F,line=1,cex=1.5)  
#mtext(side=3,text=paste("max point size = ",max(x$M.utotal)),adj=.9,line=0,cex=.9)  
print("OK")

###########################


dev.off()

}

