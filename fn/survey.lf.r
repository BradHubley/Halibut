survey.lf<-function()
{

#commercial index
x2<-lf.10 
names(x2)<-c("year","gearcd.id","gear","hookck.id","hook.type","hooksize","tripcd.id","triptype","nafo","depth","sex","len","numatlen")
x2<-x2[!is.na(x2$year)&!is.na(x2$sex)&!is.na(x2$numatlen),]
#fixed staion survey
x3<-lf.4
names(x3)<-c("year","gearcd.id","gear","hookck.id","hook.type","hooksize","tripcd.id","triptype","nafo","depth","sex","len","numatlen")
x3<-x3[!is.na(x3$year)&!is.na(x3$sex)&!is.na(x3$numatlen),]

x2$nafo=as.character(x2$nafo)
x3$nafo=as.character(x3$nafo)

x2$nafo2<-substring(x2$nafo,1,2)	
x2$nafo2[x2$nafo=="4VS"]<-x2$nafo[x2$nafo=="4VS"]
x2$nafo2[x2$nafo=="4VN"]<-x2$nafo[x2$nafo=="4VN"]
x2$nafo2[x2$nafo=="4PS"]<-x2$nafo[x2$nafo=="4PS"]

x3$nafo2<-substring(x3$nafo,1,2)	
x3$nafo2[x3$nafo=="4VS"]<-x3$nafo[x3$nafo=="4VS"]
x3$nafo2[x3$nafo=="4VN"]<-x3$nafo[x3$nafo=="4VN"]
x3$nafo2[x3$nafo=="4PS"]<-x3$nafo[x3$nafo=="4PS"]

sex.index<-c(1,2)

postscript("plot.survey.lf.ps",horizontal=F)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Commercial index
par(mfrow=c(10,2),omi=c(1.5,1,2.5,.5)*.5,mar=c(5,3,3,.5)*.5,las=1)
 

year.index<-sort(unique(x2$year)) 
nafo.index<-sort(unique(x2$nafo2))
#"3N"  "3O"  "3P"  "4T"  "4V"  "4VN" "4VS" "4W"  "4X"  "5Z" 
#length categories
a1<-seq(10,223,3)	#for .dat file
#a1<-seq(72,223,3)	#for plotting
m10<-matrix(0,nrow=length(year.index),ncol=length(a1))
f10<-matrix(0,nrow=length(year.index),ncol=length(a1))

  for(y in 1:length(year.index))
    {
  for(i in 1:length(a1))
    {
    m10[y,i]<-sum(x2$numatlen[x2$sex==sex.index[1]&x2$year==year.index[y]&x2$len>a1[i]&x2$len<=a1[i+1]],na.rm=TRUE)
    f10[y,i]<-sum(x2$numatlen[x2$sex==sex.index[2]&x2$year==year.index[y]&x2$len>a1[i]&x2$len<=a1[i+1]],na.rm=TRUE)
    }
    if(sum(m10[y,],na.rm=TRUE)<=1)
     {
     plot(c(1,1),c(1,1),type='n',axes=FALSE)	
     title('NO DATA',line=-2)
     }
    if(sum(f10[y,],na.rm=TRUE)<=1)
     {
     plot(c(1,1),c(1,1),type='n',axes=FALSE)	
     title('NO DATA',line=-2)
     }
    if(sum(m10[y,],na.rm=TRUE)>1)
     {
      barplot(m10[y,],xlab='',axes=FALSE,names.arg=a1,xpd=FALSE,cex.axis=.7)
	axis(2,at=pretty(m10[y,],n=4))
  	mtext(year.index[y],side=3,line=0,adj=.5,outer=FALSE)
     }
    if(sum(f10[y,],na.rm=TRUE)>1)
     {
      barplot(f10[y,],xlab='',axes=FALSE,names.arg=a1,xpd=FALSE,cex.axis=.7)
	axis(2,at=pretty(f10[y,],n=4))
  	mtext(year.index[y],side=3,line=0,adj=.5,outer=FALSE)
     }
    mtext("Male",side=3,line=-1,adj=0,outer=TRUE)
    mtext("Female",side=3,line=-1,adj=1,outer=TRUE)
    mtext('Commercial index',side=3,outer=TRUE,line=1.5,cex=1.5)
    mtext('Length (cm)',side=1,outer=TRUE,line=1.5,cex=1.5)
    }
    

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Fixed station survey
par(mfrow=c(10,2),omi=c(1.5,1,2.5,.5)*.5,mar=c(5,3,3,.5)*.5,las=1)

year.index<-sort(unique(x3$year)) 
nafo.index<-sort(unique(x3$nafo2))
#"3N"  "3O"  "3P"  "4T"  "4V"  "4VN" "4VS" "4W"  "4X"  "5Z" 
#length categories
a1<-seq(10,223,3)	#for .dat file
#a1<-seq(72,223,3)	#for plotting
m4<-matrix(0,nrow=length(year.index),ncol=length(a1))
f4<-matrix(0,nrow=length(year.index),ncol=length(a1))

  for(y in 1:length(year.index))
    {
  for(i in 1:length(a1))
    {
    m4[y,i]<-sum(x3$numatlen[x3$sex==sex.index[1]&x3$year==year.index[y]&x3$len>a1[i]&x3$len<=a1[i+1]],na.rm=TRUE)
    f4[y,i]<-sum(x3$numatlen[x3$sex==sex.index[2]&x3$year==year.index[y]&x3$len>a1[i]&x3$len<=a1[i+1]],na.rm=TRUE)
    }
    if(sum(m4[y,],na.rm=TRUE)<=1)
     {
     plot(c(1,1),c(1,1),type='n',axes=FALSE)	
     title('NO DATA',line=-2)
     }
    if(sum(f4[y,],na.rm=TRUE)<=1)
     {
     plot(c(1,1),c(1,1),type='n',axes=FALSE)	
     title('NO DATA',line=-2)
     }
    if(sum(m4[y,],na.rm=TRUE)>1)
     {
      barplot(m4[y,],xlab='',axes=FALSE,names.arg=a1,xpd=FALSE,cex.axis=.7)
	axis(2,at=pretty(m4[y,],n=4))
  	mtext(year.index[y],side=3,line=0,adj=.5,outer=FALSE)
     }
    if(sum(f4[y,],na.rm=TRUE)>1)
     {
      barplot(f4[y,],xlab='',axes=FALSE,names.arg=a1,xpd=FALSE,cex.axis=.7)
	axis(2,at=pretty(f4[y,],n=4))
  	mtext(year.index[y],side=3,line=0,adj=.5,outer=FALSE)
     }
    mtext("Male",side=3,line=-1,adj=0,outer=TRUE)
    mtext("Female",side=3,line=-1,adj=1,outer=TRUE)
    mtext('Fixed station survey',side=3,outer=TRUE,line=1.5,cex=1.5)
    mtext('Length (cm)',side=1,outer=TRUE,line=1.5,cex=1.5)
    }

#output
f4<<-f4
m4<<-m4
f10<<-f10
m10<<-m10

#n	method	year	sex	10	13	16	19	22	25	28	31	34	37	40	43	46	49	52	55	58	61	64	67	70	73	76	79	82	85	88	91	94	97	100	103	106	109	112	115	118	121	124	127	130	133	136	139	142	145	148	151	154	157	160	163	166	169	172	175	178	181	184	187	190	193	196	199	202	205	208	211	214	217	220	223
f4b<-data.frame(rowSums(f4),rep(7,length(year.index)),year.index,rep(1,length(year.index)),f4)
m4b<-data.frame(rowSums(m4),rep(7,length(year.index)),year.index,rep(2,length(year.index)),m4)
f10b<-data.frame(rowSums(f10),rep(8,length(year.index)),year.index,rep(1,length(year.index)),f10)
m10b<-data.frame(rowSums(m10),rep(8,length(year.index)),year.index,rep(2,length(year.index)),m10)
names(f4b)<-c("n","method","year","sex",a1)
names(m4b)<-c("n","method","year","sex",a1)
names(f10b)<-c("n","method","year","sex",a1)
names(m10b)<-c("n","method","year","sex",a1)

survey.lf<-rbind(f4b,m4b,f10b,m10b)

write.table(survey.lf, file = "survey.lf.txt", sep = "\t", row.names = F, quote = F)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#summer survey LF data  
#see get.rv2.r


dev.off()

return()
}

