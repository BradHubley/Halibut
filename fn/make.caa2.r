make.caa2<-function()
{

x<-vpop10
#x<-hal22b
years<-1:length(1988:2009)

caa.f<-matrix(0,nrow=dim(x$obs.catlen.comm2[x$obs.catlen.comm2[,2]==1&x$obs.catlen.comm2[,4]==1,])[1],ncol=x$last.age)
caa.m<-matrix(0,nrow=dim(x$obs.catlen.comm2[x$obs.catlen.comm2[,2]==1&x$obs.catlen.comm2[,4]==2,])[1],ncol=x$last.age)
caa.ot<-matrix(0,nrow=dim(x$obs.catlen.comm2[x$obs.catlen.comm2[,2]==2,])[1],ncol=x$last.age)
caa.rv.f<-matrix(0,nrow=dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==1&x$obs.catlen.rv2[,4]==1,])[1],ncol=x$last.age)
caa.rv.m<-matrix(0,nrow=dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==1&x$obs.catlen.rv2[,4]==2,])[1],ncol=x$last.age)
caa.hal.survey.f<-matrix(0,nrow=dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==7&x$obs.catlen.rv2[,4]==1,])[1],ncol=x$last.age)
caa.hal.survey.m<-matrix(0,nrow=dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==7&x$obs.catlen.rv2[,4]==2,])[1],ncol=x$last.age)

#caa.f<-matrix(0,nrow=years,ncol=x$last.age)
#caa.m<-matrix(0,nrow=years,ncol=x$last.age)
#caa.ot<-matrix(0,nrow=years,ncol=x$last.age)

for(i in 1:(dim(x$obs.catlen.comm2)[1]))
{
if(x$obs.catlen.comm2[i,2]==1&x$obs.catlen.comm2[i,4]==1)
{
caa.f[i,]<-colSums(x$obs.catlen.comm2[i,(5:(dim(x$obs.catlen.comm2)[2]))]*t(x$F.lenage.comm))
}

if(x$obs.catlen.comm2[i,2]==1&x$obs.catlen.comm2[i,4]==2)
{
caa.m[(i-(dim(x$obs.catlen.comm2[x$obs.catlen.comm2[,2]==1&x$obs.catlen.comm2[,4]==1,])[1])),]<-colSums(x$obs.catlen.comm2[i,5:dim(x$obs.catlen.comm2)[2]]*t(x$M.lenage.comm))
}

if(x$obs.catlen.comm2[i,2]==2)
{
caa.ot[(i-(dim(x$obs.catlen.comm2[x$obs.catlen.comm2[,2]==1,])[1])),]<-colSums(x$obs.catlen.comm2[i,5:dim(x$obs.catlen.comm2)[2]]*t(x$F.lenage.comm))
}

}

for(i in 1:(dim(x$obs.catlen.rv2)[1]))
{
if(x$obs.catlen.rv2[i,2]==1&x$obs.catlen.rv2[i,4]==1)
{
caa.rv.f[i,]<-colSums(x$obs.catlen.rv2[i,5:dim(x$obs.catlen.rv2)[2]]*t(x$F.lenage))
}
if(x$obs.catlen.rv2[i,2]==1&x$obs.catlen.rv2[i,4]==1)
{
start<-dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==1&x$obs.catlen.rv2[,4]==1,])[1]
caa.rv.m[i,]<-colSums(x$obs.catlen.rv2[i,5:dim(x$obs.catlen.rv2)[2]]*t(x$M.lenage))
}
if(x$obs.catlen.rv2[i,2]==7&x$obs.catlen.rv2[i,4]==1)
{
start<-dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==1&x$obs.catlen.rv2[,4]==1,])[1] + 
		dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==1&x$obs.catlen.rv2[,4]==2,])[1] 
caa.hal.survey.f[i-start,]<-colSums(x$obs.catlen.rv2[i,5:dim(x$obs.catlen.rv2)[2]]*t(x$F.lenage))
}
if(x$obs.catlen.rv2[i,2]==7&x$obs.catlen.rv2[i,4]==2)
{
start<-dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==1&x$obs.catlen.rv2[,4]==1,])[1] + 
		dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==1&x$obs.catlen.rv2[,4]==2,])[1] +
		dim(x$obs.catlen.rv2[x$obs.catlen.rv2[,2]==7&x$obs.catlen.rv2[,4]==1,])[1]
caa.hal.survey.m[i-start,]<-colSums(x$obs.catlen.rv2[i,5:dim(x$obs.catlen.rv2)[2]]*t(x$M.lenage))
}

}

caa.rv.f<<-caa.rv.f
caa.rv.m<<-caa.rv.m
caa.rv.all<<-caa.rv.f+caa.rv.m
caa.hal.survey.f<<-caa.hal.survey.f
caa.hal.survey.m<<-caa.hal.survey.m
caa.hal.survey.all<<-caa.hal.survey.f+caa.hal.survey.m

#drop 1970-1987  
caa.rv.f2<<-caa.rv.f[-c(1:18),]
caa.rv.m2<<-caa.rv.m[-c(1:18),]
caa.rv.all2<<-caa.rv.all[-c(1:18),]

#drop 1984-1987  first 4 rows
caa.ot<-caa.ot[-c(1:4),]
caa.f<<-caa.f
caa.m<<-caa.m
caa.ot<<-caa.ot
caa.all<<-caa.f+caa.m+caa.ot


}
