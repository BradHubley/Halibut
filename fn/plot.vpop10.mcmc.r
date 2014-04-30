plot.vpop10.mcmc<-function(newdata=T)
{ 
if (newdata==T)
{
#mcdata<<-getmcout("c:\\kurtis\\halibut\\adm\\mcout.dat",burn=200)
mcdata<<-getmcout("~/halibut/adm/mcout.dat",burn=200)
halstd<<-read.table("~/halibut/adm/vpop10.std",skip=1,row.names=NULL,col.names=c("number","name","value","std.dev"))

} 

x<-mcdata
x.std<-halstd
#print(names(x$M.int1))

###################################################################
#graphsheet(orientation="portrait")
#windows(record = TRUE)
postscript("plot.vpop10.mcmc.ps",horizontal=TRUE)
par(mfrow=c(4,4),las=1,omi=c(0,.5,0.5,.5),mar=c(4,4,.5,1.5))

for(i in 1:length(names(x)))
{
 var1<-x[[i]]
 temp<-data.frame(var1,x$ofv)
# var1<-var1[!is.na(var1)]

if(is.matrix(var1)==T)
{
for(j in 1:dim(var1)[2])
{

if(!(min(var1[,j],na.rm=T)==max(var1[,j],na.rm=T)))
  {
  hist(var1[,j],breaks=seq(min(var1[,j],na.rm=T),max(var1[,j],na.rm=T),length=25),probability=T,cex=.7,xlab="",plot=T,main="")
  mtext(paste(names(x)[i],j,sep="-"),side=1,line=2,cex=.7)
  abline(v=x.std$value[x.std$name==names(mcdata)[i]][j],lty=8)
#  a<-hist(var1,plot=F,breaks=seq(min(var1),max(var1),length=25))  #play with nclass
#  plot(a$breaks,c(0,a$counts),type='l')
#  abline(v=x.std$value[x.std$name==names(mcdata)[i]],lty=8)
#  print(x.std$value[x.std$name==names(mcdata)[i]])
#junk<-hist(x$a.age0,breaks=seq(min(var1),max(var1),length=25),probability=T,plot=F)$counts
 
#  mtext("Probability Density",2,outer=T,cex=1.4,line=4)
#segments(z$a.age0,0,z$a.age0,max(junk),lty=8,lwd=3)

 plot(var1[,j],type="l",xlab="",ylab="")
 mtext(paste(names(x)[i],j,sep="-"),side=2,line=3.5,cex=.7,las=0)

 acf(var1[,j],plot=T,xlab="",ylab="") 
 mtext("Lag",side=1,line=2,cex=.7)
 mtext("ACF",side=2,line=3,cex=.7,las=0)

 plot(temp[,j],temp$x.ofv[!is.na(var1[,j])],xlab="",ylab="")
 mtext(paste(names(x)[i],j,sep="-"),side=1,line=2,cex=.7,las=0)
 mtext("ofv",side=2,line=3.5,cex=.7,las=0)

  }  #close if
}  #close j
}else{

if(!(min(var1,na.rm=T)==max(var1,na.rm=T)))
  {
  hist(var1,breaks=seq(min(var1,na.rm=T),max(var1,na.rm=T),length=25),probability=T,cex=.7,xlab="",plot=T,main="")
  mtext(names(x)[i],side=1,line=2,cex=.7)
  abline(v=x.std$value[x.std$name==names(mcdata)[i]],lty=8)
#  a<-hist(var1,plot=F,breaks=seq(min(var1),max(var1),length=25))  #play with nclass
#  plot(a$breaks,c(0,a$counts),type='l')
#  abline(v=x.std$value[x.std$name==names(mcdata)[i]],lty=8)
#  print(x.std$value[x.std$name==names(mcdata)[i]])
#junk<-hist(x$a.age0,breaks=seq(min(var1),max(var1),length=25),probability=T,plot=F)$counts
 
#  mtext("Probability Density",2,outer=T,cex=1.4,line=4)
#segments(z$a.age0,0,z$a.age0,max(junk),lty=8,lwd=3)

 plot(var1,type="l",xlab="",ylab="")
 mtext(names(x)[i],side=2,line=3.5,cex=.7,las=0)

 acf(var1,plot=T,xlab="",ylab="") 
 mtext("Lag",side=1,line=2,cex=.7)
 mtext("ACF",side=2,line=3,cex=.7,las=0)

# plot(temp$var1,temp$x.ofv[!is.na(var1)],xlab="",ylab="")
# mtext(names(x)[i],side=1,line=2,cex=.7,las=0)
# mtext("ofv",side=2,line=3.5,cex=.7,las=0)

breaks<-hist(var1,plot=F,breaks=10)[1]
breaks<-as.numeric(as.vector(unlist(breaks)))
group<-findInterval(var1,breaks)
group.index<-sort(unique(group))
temp<-data.frame(var1,temp$x.ofv[!is.na(var1)],group)
print(names(x)[i])
ofv<-NULL
bin<-NULL

for(k in 1:(length(group.index)))
  {
  ofv[k]<-min(split(temp,group)[[k]][2])
  bin[k]<-mean(split(temp,group)[[k]][1])
  }

plot(bin,ofv,type='l',lwd=2,col='red',xlab="",ylab="")
v=x.std$value[x.std$name==names(mcdata)[i]]
# print(v)
 abline(v=x.std$value[x.std$name==names(mcdata)[i]],lty=8)
 mtext(names(x)[i],side=1,line=2,cex=.7,las=0)
 mtext("ofv",side=2,line=3.5,cex=.7,las=0)


  } # close if
 }  #close else
 }  #close i
dev.off()
}








