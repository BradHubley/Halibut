plot.caa1<-function(update=T)
{ 
if(update==T)
{
hal.caa<-getadmout("admb/caa1.rep")	
hal.caa.std<-read.table("admb/caa1.std",skip=1,row.names=NULL,
			col.names=c("number","name","value","std.dev"))
}
if(update==T)
{
vpop10<-getadmout("admb/vpop10.rep")	
vpop10.std<-read.table("admb/vpop10.std",skip=1,row.names=NULL,
			col.names=c("number","name","value","std.dev"))
}

x<-hal.caa
x2<-hal.caa.std


first.year=x$first.year
last.year=x$last.year
first.age=x$first.age
last.age=x$last.age                
n.landings=x$n.landings
landing.year=x$landing.year
landings=x$landings
mat=x$mat                    
obs.comcaa=x$obs.comcaa
obs.rvcaa=x$obs.rvcaa
obs.rvwaa=x$obs.rvwaa
obs.rvspring=x$obs.rvspring          
Nage=x$Nage
M=x$M
Cage=x$Cage
RVage=x$RVage              
RVspring.age=x$RVspring.age
obs.rvspring.total.catch=x$obs.rvspring.total.catch
pred.rvspring.total.catch=x$pred.rvspring.total.catch
obs.Cpage=x$obs.Cpage                
pred.Cpage=x$pred.Cpage
obs.RVpage=x$obs.RVpage
pred.RVpage=x$pred.RVpage
Cpage.di=x$Cpage.di                
RVpage.di=x$RVpage.di
RVspring.page.di=x$RVspring.page.di
obs.RVspring.page=x$obs.RVspring.page
pred.RVspring.page=x$pred.RVspring.page       
F=x$F
F.full=x$F.full
ssb=x$ssb
C.Sfullold=x$C.Sfullold              
C.log.varLold=x$C.log.varLold
C.log.varRold=x$C.log.varRold
select=x$select
RV.Sfullold=x$RV.Sfullold             
RV.log.varLold=x$RV.log.varLold
RV.log.varRold=x$RV.log.varRold
RVspring.Sfull=x$RVspring.Sfull
RVspring.log.varL=x$RVspring.log.varL        
RVspring.log.varR=x$RVspring.log.varR
RV.log.qold=x$RV.log.qold
RV.select=x$RV.select
RVspring.select=x$RVspring.select          
obs.total.catch=x$obs.total.catch
pred.total.catch=x$pred.total.catch
year1=x$year1
mean.RVwt.age=x$mean.RVwt.age            
mean.rvwaa=x$mean.rvwaa
Wtage=x$Wtage
ctotal.like=x$ctotal.like
Cpage.like=x$Cpage.like               
rvtotal.like=x$rvtotal.like
RVpage.like=x$RVpage.like
rvspring_total.like=x$rvspring_total.like
RVspring_page.like=x$RVspring_page.like       
ctotal.sigma=x$ctotal.sigma
rvtotal.sigma=x$rvtotal.sigma
ofv=x$ofv                


postscript("plot.caa1.ps",horizontal=FALSE)
#win.metafile("Rplot%02d.wmf", pointsize = 10)
# PLOT MEAN WEIGHT AT AGE FOR COMMERCIAL CATCH 

par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(5,6,3,1),las=1)


# PLOT FITTED AND OBSERVED TOTAL CATCHES FOR COMMERCIAL AND RV DATA 

par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(5,6,3,1),las=1,xpd=NA)

# COMMERCIAL DATA 
year<-seq(1988,last.year,1)
age<-seq(first.age,last.age,1)

pred.total.catch<-pred.total.catch[19:40]
ymax<-max(c(obs.total.catch,pred.total.catch))/1000
plot(year,pred.total.catch/1000,pch=1,ylim=c(0,.25),xlim=c(1988,2010),xlab="",ylab="",type="l",axes=F,lwd=2)
axis(1,at=seq(1988,2010,2))
axis(2,seq(0,.25,.05))
points(year,obs.total.catch/1000,pch=16)  					# black dots are observed values
mtext(side=2,line=3,"Total catch (millions of fish)",cex=1,las=0)  
mtext(side=1,line=2,"Year",cex=1)  
mtext("Fig. 4",3,adj=1,line=5,cex=1,outer=T)
text(1988,ymax+(ymax/4),"(a)",cex=1.7)	

year<-seq(first.year,last.year,1)
obs.total.rv <- rowSums(obs.rvcaa)
pred.total.rv <- rowSums(RVage)
ymax<-max(c(obs.total.rv,pred.total.rv))/1000
plot(year,pred.total.rv/1000,pch=1,ylim=c(0,2.5),xlim=c(first.year,2010),xlab="",ylab="",type="l",axes=F,lwd=2)
lines(1998:last.year,pred.rvspring.total.catch/10^3,lty=2,lwd=2)
axis(1,at=seq(first.year,2010,10))
axis(2,seq(0,2.5,.5))
points(year,obs.total.rv/1000,pch=16)
points(1998:last.year,obs.rvspring.total.catch/10^3,pch=3)
mtext(side=2,line=4.5,"Index of abundance",cex=1,las=0)  
mtext(side=2,line=3,"(millions of fish)",cex=1,las=0)  
mtext(side=1,line=2,"Year",cex=1)  
text(first.year,ymax+(ymax/5),"(b)",cex=1.7)	
legend(1990,2.5, c("RV summer 1970-2009","Halibut survey 1998-2009"),pch=c(16,3),lty=c(1,2),bty='n',lwd=2,cex=.8) 


# PLOT YEAR 1 FOR ALL AGES 

#year1<-x2$value[x2$name=="year1"]/1000
#year1.se <- x2$std.dev[x2$name=="year1"]/1000
#age<-1:length(year1)  
#ymax<-max(year1+year1.se*2.1)
#ymin<-min(year1-year1.se*2.1)

#plot(age,year1,type="l",pch=1,xlab="",ylab="",ylim=c(ymin,ymax),lwd=2) 
#lines(age,year1+2*year1.se,lty=2,lwd=2)  
#lines(age,year1-2*year1.se,lty=2,lwd=2)
#mtext(side=2,line=4,"Millions of  fish",cex=0.9,las=0)  
#mtext(side=3,line=1,"Number of  fish in first year",cex=1) 
#mtext(side=1,line=3,"Age",cex=1)  

year1<-Nage[1,]
age<-1:last.age  
ymax<-max(year1*1.5)
ymin<-min(year1)

plot(age,year1,type="l",pch=1,xlab="",ylab="",ylim=c(ymin,ymax),lwd=2) 
mtext(side=2,line=4,"Millions of  fish",cex=0.9,las=0)  
mtext(side=3,line=1,"Number of fish in first year",cex=1) 
mtext(side=1,line=3,"Age",cex=1)  

# PLOT AGE 1 IN EACH YEAR 
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(5,6,3,1),las=1)
year2<-seq(first.year+1,last.year,1)

age1<-x2$value[x2$name=="age1"]
age1.se <-x2$std.dev[x2$name=="age1"]
se.l<-(age1-age1.se)/1000
se.u<-(age1+age1.se)/1000
age1<-(age1)/1000
ymax<-max(se.u*1.2)

plot(year2[1:37],age1[1:37],type="l",xlab="",ylab="",ylim=c(0,max(age1[1:37])*1.2),xlim=c(first.year,2010),axes=F,lwd=2)
axis(1,at=seq(first.year,last.year,5))
axis(2)
#first.year not plotted. Don't plot terminal year
lines(year2[1:10],se.u[1:10],lty=2,lwd=2)  
lines(year2[12:37],se.u[12:37],lty=2,lwd=2)  
lines(year2[1:10],se.l[1:10],lty=2,lwd=2)  
lines(year2[12:37],se.l[12:37],lty=2,lwd=2)  
mtext(side=1,line=3,"Year",cex=1)  
mtext(side=2,line=4.5,"Millions of hal recruits (age-1)",cex=1,las=0)  
mtext("Fig. 6",3,adj=1,line=5,cex=1,outer=T)
#text(first.year,ymax+(ymax/5),"(a)",cex=1.7)	


# PLOT SSB 
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(5,6,3,1),las=1)
ssb<-x2$value[x2$name=="ssb"]/1000
ssb.se <- x2$std.dev[x2$name=="ssb"]/1000
ymax<-max(ssb[1:length(year)-1]+ssb.se[1:length(year)-1]*2.1)
ymin<-min(ssb[1:length(year)-1]-ssb.se[1:length(year)-1]*2.1)

plot(year,ssb,type="l",xlab="",ylab="",ylim=c(0,max(ssb)*1.5),xlim=c(first.year,2010),axes=F)
axis(1,at=seq(first.year,2010,5))
axis(2)
lines(year[1:length(year)-1],ssb[1:length(year)-1]+2*ssb.se[1:length(year)-1],lty=2)  
lines(year[1:length(year)-1],ssb[1:length(year)-1]-2*ssb.se[1:length(year)-1],lty=2)
mtext(side=1,line=2,"Year",cex=1)  
mtext(side=2,line=4,"SSB (Thousands of tonnes)",cex=1,las=0)  
#mtext("Fig. ",3,adj=1,line=5,cex=1,outer=T)
text(first.year,ymax+(ymax/5),"(a)",cex=1.7)	

#Plot total biomass and total numbers
Total.W<-rowSums(x$Nage*x$obs.rvwaa)
Total.N<-rowSums(x$Nage)
years<-first.year:last.year
ymax<-max(Total.N)/1000

plot(years,Total.N/1000,type="n",xlab="",ylab="",ylim=c(0,ymax),axes=F)
lines(years,Total.N/1000,lty=1,lwd=2)	#total numbers in millions
axis(1,at=seq(first.year,2010,5))
axis(2,seq(0,2,.5))
legend(1993,2.2,legend=rev(c("Numbers","Biomass")),lty=rev(c(1,2)),lwd=2,bty="n",cex=1.1)
text(first.year,ymax+(ymax/4),"(b)",cex=1.7)	
par(new=T,xpd=NA)
plot(years,Total.W/1000,type="n",xlab="",ylab="",ylim=c(0,max((Total.W/1000)[-1])),axes=F)
axis(4,seq(0,20,5))
lines(years[-1],(Total.W/1000)[-1],lty=2,lwd=2)		#total biomass in 1000 t
mtext(side=2,text="N (millions)",cex=1.5,line=3,las=0)
mtext(side=4,text="Total biomass (1000 t)",cex=1.5,line=3,las=0)

# PLOT SSB from 3 models CAL, CAA and VPA
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(5,6,3,1),las=1)

years<-first.year:last.year
#CAL
ssb2<-vpop10$SSB2/1000			#males and females
ssb2<-vpop10.std$value[vpop10.std$name=="SSB"]/1000		#females only
ssb2.se <- vpop10.std$std.dev[vpop10.std$name=="SSB"]/1000
#CAA
ssb<-x2$value[x2$name=="ssb"]/1000
ssb.se <- x2$std.dev[x2$name=="ssb"]/1000
ymax<-max(ssb[1:length(year)-1]+ssb.se[1:length(year)-1]*2.1)
ymin<-min(ssb[1:length(year)-1]-ssb.se[1:length(year)-1]*2.1)
#VPA
#vpa<-read.table("vpa.txt",header=T)		#results from Bob's vpa model
#ssb3<-vpa$ssb/1000
browser()
plot(years,ssb,type="l",xlab="",ylab="",ylim=c(0,max(ssb)*1.5),xlim=c(first.year,2010),axes=F,col='red',lwd=2)  #CAA
axis(1,at=seq(first.year,2010,5))
axis(2)
lines(years,ssb2[1:length(years)],col='black',lwd=2)  				#CAL
#lines(years,ssb3[1:length(years)],col='green',lwd=2)  				#VPA
#lines(year[1:length(year)-1],ssb2[1:length(year)-1]+2*ssb2.se[1:length(year)-1],lty=2)  
#lines(year[1:length(year)-1],ssb2[1:length(year)-1]-2*ssb2.se[1:length(year)-1],lty=2)
mtext(side=1,line=2,"Year",cex=1)  
mtext(side=2,line=4,"SSB (Thousands of tonnes)",cex=1,las=0)  
text(first.year,ymax+ymax/10,"(a)",cex=1.7)
#legend(1995,max(c(ssb,ssb2)),legend=(c("CAL","CAA","VPA")),col=(c("black","red","green")),lwd=2,bty="n",cex=.9)
legend(1995,max(c(ssb,ssb2)),legend=(c("CAL","CAA")),col=(c("black","red")),lwd=2,bty="n",cex=.9)

# PLOT EXPLOITATION RATE BY YEAR
year<-first.year:last.year

#CAL
F.cal.f<-vpop10$F.utotal[11:50,10]
F.cal.m<-vpop10$M.utotal[11:50,7]
#CAA
F.full[1:18]<-NA
F.se <- x2$std.dev[x2$name=="F_full"][-(1:18)]
ymax<-max(c(F.full,F.cal.f,F.cal.m),na.rm=T)*1.2
ymin<-min(c(F.full,F.cal.f,F.cal.m),na.rm=T)*1.2
#VPA
#Favg<-vpa$Favg

plot(year[-(1:18)],F.full[-(1:18)],ylim=c(0,ymax*1.1),xlim=c(first.year,2010),type="l",xlab="",ylab="",pch=1,axes=F,col='red',lwd=2)
axis(1,at=seq(first.year,2010,5))
axis(2)
mtext(side=1,line=2,"Year",cex=1) 
lines(year,F.cal.f,lty=1,col='black',lwd=2)	#	CAL Females 
lines(year,F.cal.m,lty=2,col='black',lwd=2)  	# 	CAL Males
#lines(year,Favg,lty=2,col='green',lwd=2)  	# 	VPA
#lines(year,vpa$Fmax,lty=2,col='green')  	# 	VPA
mtext(side=2,line=3,"Fishing mortality rate",cex=1,las=0)   
text(first.year,ymax+ymax/5,"(b)",cex=1.7)
#legend(1995,ymax,legend=(c("CAL females","CAL males","CAA","VPA")),col=(c("black","black","red","green")),lwd=2,lty=c(1,2,1,1),bty="n",cex=.9)
legend(1995,ymax,legend=(c("CAL females","CAL males","CAA")),col=(c("black","black","red")),lwd=2,lty=c(1,2,1),bty="n",cex=.9)
print("OK")

# PLOT SELECTIVITY CURVES
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(5,6,3,1),las=1)
plot(first.age:last.age,select[1,],ylim=c(0,1.01),xlim=c(0,12),type="n",xlab="",ylab="",pch=1,axes=F)
axis(1,at=seq(0,12,2))
axis(2)
lines(first.age:last.age,select[1,],lty=8,lwd=2)
mtext(side=1,line=2,"Age",cex=1.2)
mtext(side=2,line=3,"Selectivity at age",cex=.9,las=0)   
mtext(side=3,line=1,"Commercial Catch",cex=1)   

plot(first.age:last.age,select[1,],ylim=c(0,1.01),xlim=c(0,12),type="n",xlab="",ylab="",pch=1,axes=F)
axis(1,at=seq(0,12,2))
axis(2)
lines(first.age:last.age,RV.select[1,],lty=1,lwd=2) 						# RV summer 1970-present
lines(first.age:last.age,RVspring.select[1,],lty=3,lwd=2) 					# RV halibut survey 1998-present
mtext(side=1,line=2,"Age",cex=1.2)
mtext(side=2,line=3,"Selectivity at age",cex=.9,las=0)   
mtext(side=3,line=1,"RV",cex=1)   
#legend(8,.5, c("first.year-1982",">1982"),lty=c(1,8)) 
legend(7,.75,c("summer 1970-2009","halibut survey 1998-2009"),lty=c(1,8),bty='n',lwd=c(2,2)) 



# PLOT FITTED AND OBSERVED VALUES BY AGE CLASS FOR EACH YEAR 

par(mfrow=c(4,3),las=1,omi=c(1,1,1.5,.5),mar=c(1,1,1,1)*3)


# COMMERCIAL CATCH DATA 

for (i in 1:length(1988:last.year)){ 
	ylim<-c(0,max(c(Cage[17+i,],obs.comcaa[i,])))  
	plot(1:length(Cage[1,]),Cage[18+i,],ylim=ylim,xlab="",
	ylab="",main=paste(1987+i),cex=0.5,axes=F,type="l",pch=1 )  
	points(1:length(obs.comcaa[1,]),obs.comcaa[i,],pch=16,cex=0.5)
	 axis(1,cex=.6)
         axis(2,cex=.6)
         box()

 	if(i==1|i==13|i==25|i==37) mtext(outer=T,side=2,line=2,"Number of fish (000's)",cex=1.2,las=0) 
	if(i==1|i==13|i==25|i==37) mtext(outer=T,side=3,line=0,"Commercial Catch-at-Age Data - ages within year",cex=1.2)   
	if(i==1|i==13|i==25|i==37) mtext(outer=T,side=1,line=2,"Year",cex=1.2)  

	} 

mtext(outer=T,side=2,line=2,"Number of fish (000's)",cex=1.2,las=0)  
mtext(outer=T,side=3,line=0,"Commercial Catch-at-Age Data - ages within year",cex=1.2)  
mtext(outer=T,side=1,line=2,"Year",cex=1.2)


# RV DATA summer 1970-2003

par(mfrow=c(4,3),las=1,omi=c(1,1,1.5,.5),mar=c(1,1,1,1)*3)

for (i in 1:length(RVage[,1])){ 
	ylim<-c(0,max(c(RVage[i,],obs.rvcaa[i,])))  
	plot(1:length(RVage[1,]),RVage[i,],type="l",ylim=ylim,xlab="",
         	ylab="",main=paste(1969+i),pch=1,cex=0.5,axes=0) 
	points(1:length(obs.rvcaa[1,]),obs.rvcaa[i,],pch=16,cex=0.5)
	 axis(1,cex=.6)
         axis(2,cex=.6)
	box()
	if(i==1|i==13|i==25|i==37)mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0) 
	if(i==1|i==13|i==25|i==37) mtext(outer=T,side=3,line=2,"RV Catch-at-Age Data - ages within year",cex=1.2)   
	if(i==1|i==13|i==25|i==37) mtext(outer=T,side=3,line=0,"Summer",cex=1.2)   
	if(i==1|i==13|i==25|i==37) mtext(outer=T,side=1,line=2,"Year",cex=1.2)  
	} 

mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0)  
mtext(outer=T,side=3,line=2,"RV Catch-at-Age Data - ages within year",cex=1.2)  
mtext(outer=T,side=3,line=0,"Summer",cex=1.2)  
mtext(outer=T,side=1,line=2,"Year",cex=1.2)


# RV DATA halibut survey 1998-2009
RVspring.age<-RVspring.age/10^3
obs.rvspring<-obs.rvspring/10^3

par(mfrow=c(4,3),las=1,omi=c(1,1,1.5,.5),mar=c(1,1,1,1)*3)

for (i in 1:length(RVspring.age[,1])){ 
	ylim<-c(0,max(c(RVspring.age[i,],obs.rvspring[i,])))  
	plot(1:length(RVspring.age[1,]),RVspring.age[i,],type="l",ylim=ylim,xlab="",
         	ylab="",main=paste(1997+i),pch=1,cex=0.5,axes=0) 
	points(1:length(obs.rvspring[1,]),obs.rvspring[i,],pch=16,cex=0.5)
	 axis(1,cex=.6)
         axis(2,cex=.6)
	box()
	if(i==1)mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0) 
	if(i==1) mtext(outer=T,side=3,line=2,"Halibut survey Catch-at-Age Data - ages within year",cex=1.2)  
	if(i==1) mtext(outer=T,side=3,line=0,"Spring",cex=1.2)  
	if(i==1) mtext(outer=T,side=1,line=2,"Year",cex=1.2)  
	} 

mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0)  
mtext(outer=T,side=3,line=2,"Halibut survey Catch-at-Age Data - ages within year",cex=1.2)  
mtext(outer=T,side=3,line=0,"Spring",cex=1.2)  
mtext(outer=T,side=1,line=2,"Year",cex=1.2)


# PLOT OBSERVED AND PREDICTED CATCH AT AGE FOR EACH COHORT 
par(mfrow=c(4,3),las=1,omi=c(1,1,1.5,.5),mar=c(1,1,1,1)*3)

# COMMERCIAL DATA 

year<-rep((first.year:last.year),length(first.age:last.age))
age1<-first.age:last.age
age<-rep(age1,each=length(first.year:last.year))
cohort<-year-age
pred.catch<-as.vector(Cage)
temp<-matrix(0,nrow=40,ncol=20)
temp[19:40,]<-obs.comcaa
#obs.catch<-as.vector(obs.comcaa) 
obs.catch<-as.vector(temp) 

pred.rv<-as.vector(RVage)
obs.rv<-as.vector(obs.rvcaa)
pred.rvspring<-as.vector(RVspring.age)
obs.rvspring<-as.vector(obs.rvspring)
pred.catch2<-as.data.frame(cbind(year,age,cohort,pred.catch,obs.catch,obs.rv,pred.rv,obs.rvspring,pred.rvspring))
# mtext(outer=T,side=2,line=0,"Number of fish ('000s)",cex=1.2,las=0)  
#pred.catch2<<-pred.catch2

for(i in 1987:2005)

  {

    pyear<-pred.catch2$year[pred.catch2$cohort==i]      
    ppred.fish<-pred.catch2$pred.catch[pred.catch2$cohort==i]
    pobs.fish<-pred.catch2$obs.catch[pred.catch2$cohort==i]
    ymax<-max(c(ppred.fish,pobs.fish))
         plot(pyear,ppred.fish,xlim=c(i+1,i+last.age),ylim=c(0,ymax),
         ylab="",xlab="",type="l",pch=1,lwd=1,axes=0,col=3,cex=0.5,
	 main=paste(i)  )
         axis(1,cex=.6)
         axis(2,cex=.6)
         box()

         points(pyear,pobs.fish,cex=0.5,pch=16,lwd=1)
	 if(i==1969|1985|1993|2005) mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0) 
	 if(i==1969|1985|1993|2005) mtext(outer=T,side=3,line=0,"Commercial Catch-at-Age Data - by cohort",cex=1.2)   
	 if(i==1969|1985|1993|2005) mtext(outer=T,side=1,line=2,"Year",cex=1.2)   
  }

mtext(outer=T,side=2,line=2,"Number of fish caught",cex=1.2,las=0)  
mtext(outer=T,side=3,line=0,"Commercial Catch-at-Age Data - by cohort",cex=1.2)  
mtext(outer=T,side=1,line=2,"Year",cex=1.2)

# RV DATA summer 1970-2003
par(mfrow=c(4,3),las=1,omi=c(1,1,1.5,.5),mar=c(1,1,1,1)*3)

for(i in 1969:2005)

  {
    pyear<-pred.catch2$year[pred.catch2$cohort==i]      
    ppred.fish<-pred.catch2$pred.rv[pred.catch2$cohort==i]
    pobs.fish<-pred.catch2$obs.rv[pred.catch2$cohort==i]
    ymax<-max(c(ppred.fish,pobs.fish))
         plot(pyear,ppred.fish,xlim=c(i+1,i+last.age),ylim=c(0,ymax),
         ylab="",xlab="",type="l",pch=1,cex=0.5,lwd=1,axes=0,col=3,
	 main=paste(i)  )
         axis(1,cex=.6)
         axis(2,cex=.6)
         box()
         points(pyear,pobs.fish,cex=.5,pch=16,lwd=1)  

 	 if(i==1969|i==1985|i==1993|i==2005) mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0) 
	 if(i==1969|i==1985|i==1993|i==2005) mtext(outer=T,side=3,line=2,"RV Catch-at-Age Data - by cohort",cex=1.2)  
	 if(i==1969|i==1985|i==1993|i==2005) mtext(outer=T,side=3,line=0,"Summer",cex=1.2)  
	 if(i==1969|i==1985|i==1993|i==2005) mtext(outer=T,side=1,line=2,"Year",cex=1.2)   

  }

mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0) 
mtext(outer=T,side=3,line=2,"RV Catch-at-Age Data - by cohort",cex=1.2)  
mtext(outer=T,side=3,line=0,"Summer",cex=1.2)  
mtext(outer=T,side=1,line=2,"Year",cex=1.2) 

# RV DATA halibut survey 1998-2009
par(mfrow=c(4,3),las=1,omi=c(1,1,1.5,.5),mar=c(1,1,1,1)*3)

for(i in 1997:2005)

  {
    pyear<-pred.catch2$year[pred.catch2$cohort==i]      
    ppred.fish<-pred.catch2$pred.rvspring[pred.catch2$cohort==i]
    pobs.fish<-pred.catch2$obs.rvspring[pred.catch2$cohort==i]
    ymax<-max(c(ppred.fish,pobs.fish))
         plot(pyear,ppred.fish,xlim=c(i+1,i+last.age),ylim=c(0,ymax),
         ylab="",xlab="",type="l",pch=1,cex=0.5,lwd=1,axes=0,col=3,
	 main=paste(i)  )
         axis(1,cex=.6)
         axis(2,cex=.6)
         box()
         points(pyear,pobs.fish,cex=.5,pch=16,lwd=1)  

	if(i==1997)mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0) 
	if(i==1997) mtext(outer=T,side=3,line=2,"Halibut survey Catch-at-Age Data - by cohort",cex=1.2)  
	if(i==1997) mtext(outer=T,side=3,line=0,"Spring",cex=1.2)  
	if(i==1997) mtext(outer=T,side=1,line=2,"Year",cex=1.2)  
	} 

mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0)  
mtext(outer=T,side=3,line=2,"Halibut survey Catch-at-Age Data - by cohort",cex=1.2)  
mtext(outer=T,side=3,line=0,"Spring",cex=1.2)  
mtext(outer=T,side=1,line=2,"Year",cex=1.2)


# PLOT OBSERVED AND PREDICTED RV CATCH OVER YEARS FOR EACH AGE 
#par(mfrow=c(4,3),las=1,omi=c(1,1,.5,.5))
par(mfrow=c(4,3),las=1,omi=c(1,1,.5,.5),mar=c(1,1,1,1)*3)

for(i in first.age:last.age)

  {
    pyear<-pred.catch2$year[pred.catch2$age==i]      
    ppred.fish<-pred.catch2$pred.rv[pred.catch2$age==i]
    pobs.fish<-pred.catch2$obs.rv[pred.catch2$age==i]

#    print(pyear)
#    print(ppred.fish)
#    print(pobs.fish)

    ymax<-max(c(ppred.fish,pobs.fish))
         plot(pyear,ppred.fish,ylim=c(0,ymax),xlim=c(first.year,2005),
         ylab="",xlab="",type="l",pch=1,cex=0.5,lwd=1,axes=F,col=3)
         axis(1,at=seq(first.year,last.year,10),labels=T,cex=.6)
         axis(2,cex=.6)
	 mtext(paste("age-",i),side=3,outer=F,line=.6,cex=.75)
         box()
         points(pyear,pobs.fish,cex=.5,pch=16,lwd=1)  

  }

mtext(outer=T,side=2,line=2,"Number of fish",cex=1.2,las=0) 
mtext(outer=T,side=3,line=.2,"RV Catch over Years - by age",cex=1.2)  
mtext(outer=T,side=1,line=0,"Year",cex=1.2) 


# RESIDUAL PLOTS COMMERCIAL DATA 
#NOTE: currently only uses 1970-2007
par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(5,6,3,1),las=1,xpd=NA)

year<-seq(1988,last.year,1)
age<-1:last.age
#res.comcaa<- log(obs.comcaa+0.001) - log(Cage) 

Cpage.rd<-matrix(data=NA,nrow=(last.year-1988+1),ncol=(last.age-first.age+1))

for(i in 1:length(year))
  {
  for(j in 1:length(age)) 
    { 
    Cpage.rd[i,j]<-sign(obs.Cpage[i,j]-pred.Cpage[i,j])*sqrt(Cpage.di[i,j])
    }
  }


maxres<-max(abs(Cpage.rd))  

plot(c(1988,last.year),c(1,last.age),type="n",xlab="",ylab="",xlim=c(first.year,2010),axes=F)
axis(1,at=seq(first.year,2010,5))
axis(2,seq(0,20,5))

for (i in 1:length(year))
  {
  for (j in 1:length(age)) 
    { 
    if(Cpage.rd[i,j]<0) points(year[i],age[j],pch=1,cex=(-Cpage.rd[i,j]/maxres)*1.75)
    if(Cpage.rd[i,j]>0) points(year[i],age[j],pch=16,cex=(Cpage.rd[i,j]/maxres)*1.75)
    }
  }

mtext("Age",side=2,line=3,cex=1,outer=F,las=0)
text(1965,last.age+(last.age/5),"(a) Commercial catch",cex=1.7,adj=0)

# RESIDUAL PLOTS RV summer

year<-seq(first.year,last.year,1)
RVpage.rd<-matrix(data=NA,nrow=(last.year-first.year+1),ncol=(last.age-first.age+1))

for(i in 1:length(year))
  {
  for(j in 1:length(age)) 
    { 
    RVpage.rd[i,j]<-sign(obs.RVpage[i,j]-pred.RVpage[i,j])*sqrt(RVpage.di[i,j])
    }
  }

maxres<-max(abs(RVpage.rd))  

plot(c(first.year,2010),c(1,last.age),type="n",xlab="",ylab="",xlim=c(first.year,2010),axes=F)
axis(1,at=seq(first.year,2010,5))
axis(2,seq(0,20,5))

for (i in 1:length(year))
  {
  for (j in 1:length(age)) 
    { 
    if(RVpage.rd[i,j]<0) points(year[i],age[j],pch=1,cex=(-RVpage.rd[i,j]/maxres)*1.75)
    if(RVpage.rd[i,j]>0) points(year[i],age[j],pch=16,cex=(RVpage.rd[i,j]/maxres)*1.75)
    }
  }

legend(2000,last.age+(last.age/3), c("obs>pred","pred>obs"),pch=c(16,1)) 
mtext("Age",side=2,line=3,cex=1,outer=F,las=0)
mtext("Year",side=1,line=3,cex=1,outer=F)
text(1965,last.age+(last.age/5),"(b) RV Summer",cex=1.7,adj=0)	
#text(2003,14,as.character(round(maxres,digits=2)),cex=1.2)



# RESIDUAL PLOTS RV halibut survey

par(mfrow=c(2,1),omi=c(1,1,1,1),mar=c(5,6,3,1),las=1,xpd=NA)

year<-1998:last.year
#res.rvspring<- log(obs.rvspring+0.001)-log(RVspring.age) 

RVspring.page.rd<-matrix(data=NA,nrow=(last.year-first.year+1),ncol=(last.age-first.age+1))

for(i in 1:length(year))
  {
  for(j in 1:length(age)) 
    { 
    RVspring.page.rd[i,j]<-sign(obs.RVspring.page[i,j]-pred.RVspring.page[i,j])*sqrt(RVspring.page.di[i,j])
    }
  }
#print(RVspring.page.rd)
max(abs(RVspring.page.rd[!is.na(RVspring.page.rd)]))


plot(c(1998,last.year),c(1,last.age),type="n",xlab="",ylab="",xlim=c(1998,2010),axes=F)
axis(1,at=seq(1998,2010,1))
axis(2,seq(0,20,5))

for (i in 1:length(year))
  {
  for (j in 1:length(age)) 
    { 
    if(RVspring.page.rd[i,j]<0) points(year[i],age[j],pch=1,cex=(-RVspring.page.rd[i,j]/maxres)*1.75)
    if(RVspring.page.rd[i,j]>0) points(year[i],age[j],pch=16,cex=(RVspring.page.rd[i,j]/maxres)*1.75)
    }
  }

mtext("Age",side=2,line=3,cex=1,outer=F,las=0)
text(1998.1,last.age+(last.age/5),"(c) Halibut survey",cex=1.7,adj=0)	

dev.off()


return()

}  


