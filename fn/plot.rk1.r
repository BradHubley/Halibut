plot.rk1<-function(update.x=T)
{
##Plot S-R curves 

if(update.x==T)
 {
 vpop10<-getadmout("~/halibut/adm/vpop10.rep")
 vpop10.std<-read.table("~/halibut/adm/vpop10.std",skip=1,row.names=NULL,
		 col.names=c("number","name","value","std.dev"))
 }
#vpop10<<-vpop10
#vpop10std<<-vpop10.std

x<-vpop10
z<-vpop10.std

postscript("plot.rk1.ps",horizontal=F)

############################################################################################################
##############################################################################
# SPAWNER-RECRUIT RELATIONSHIP
par(las=1,mfrow=c(2,1),omi=c(1.5,1,1,1),mar=c(4,4,2,1),cex=.9)

ssb<-x$SSB2[-40]/1000
r<-x$SSNrec[2,-1]/1000
years<-(1970:2008)
x1<-data.frame(years,ssb,r)
#drop 1981?? NO controlled with a random walk in .tpl
#x1<-x1[-11,]
#r<-r[-11]
#ssb<-ssb[-11]
a4<-nls(r~2.71828*p2/p1*(ssb)*exp(-(ssb)/p1),data=x1,start=c(p1=2,p2=200),na.action=na.omit)
a5<-ml.rklnl(ssb,r)



ssb2<-seq(min(ssb),max(ssb),length=100)
r5<-srfrk(a5[[1]][1],a5[[1]][2],ssb2)
r6<-srfrk(a5[[1]][1],a5[[1]][2],ssb)

r5<-srfrkpv(a5[[1]],ssb2)
r6<-srfrkpv(a5[[1]],ssb)

resid<-r-r6
temp<-cbind(years,resid)

#plot ricker curve
plot(ssb,r,xlab='',ylab='',type='n')
points(ssb,r,pch=substring(years+1,3,3),cex=.9,col='dark green')
points(ssb+.1,r,pch=substring(years+1,4,4),cex=.9,col='dark green')
lines(ssb2,r5,lty=1,lwd=2,xpd=F)
mtext(side=3,text="Atlantic halibut (NAFO: 3NOPs4VWX5Zc)",cex=1.3,line=.2,adj=.5,outer=F)
mtext(side=1,text="SSB (1000 t)",cex=1.5,line=2.5)
mtext(side=2,text="Age-1 recruits (1000s)",cex=1.5,line=3.5,las=0)

print("OK")



dev.off()
}
