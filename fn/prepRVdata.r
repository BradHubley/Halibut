prepRVdata<-function(RVdata,years,bins){
	
	RVdata[[1]]$year<-as.numeric(format(RVdata[[1]]$SDATE,"%Y"))
	RVdata[[2]]$year<-as.numeric(format(RVdata[[2]]$SDATE,"%Y"))
	if(missing(years))years<-sort(unique(RVdata[[2]]$year))
	if(missing(bins))bins<-seq(10,223,3)
	
	#RVdata[[1]]<-rbind(RVdata[[1]], subset(RVdata[[1]],CLEN==2), subset(RVdata[[1]],CLEN==3), subset(RVdata[[1]],CLEN==3))
	
	N<-c()
	B<-c()
	survF.lst<-list()
	stratF.lst<-list()
	F.lst<-list()
	survM.lst<-list()
	stratM.lst<-list()
	M.lst<-list()
	surv.lst<-list()
	strat.lst<-list()
	Total.lst<-list()
	
	print(paste("Year", "Nstr", "N(000s)", "B(t)"))
	for (i in 1:length(years)){
		
		# Index of abundance
		Wh<-with(subset(RVdata[[2]],year==years[i]),tapply(TUNITS,STRAT,unique))
		Nh<-with(subset(RVdata[[2]],year==years[i]),tapply(TOTNO,STRAT,mean))
		Bh<-with(subset(RVdata[[2]],year==years[i]),tapply(TOTWGT,STRAT,mean))
		N[i]<-sum(Nh*Wh)
		B[i]<-sum(Bh*Wh)
		print(paste(years[i],length(Wh),round(N[i]/1000),round(B[i]/1000)))
		
		survF.lst[[i]]<-subset(RVdata[[2]],year==years[i])
		survM.lst[[i]]<-subset(RVdata[[2]],year==years[i])
		surv.lst[[i]]<-subset(RVdata[[2]],year==years[i])
		# Length frequency
		for(j in 1:length(bins)){
			NbinhF<-with(subset(RVdata[[1]],year==years[i]&FLEN>=bins[j]&FLEN<bins[j+1]&FSEX==2),tapply(CLEN,SETNO,sum))
			NbinhF.dat<-data.frame(SETNO=as.numeric(names(NbinhF)),NbinhF)
			names(NbinhF.dat)[2]<-paste("bin",bins[j],sep='')
			survF.lst[[i]]<-merge(survF.lst[[i]],NbinhF.dat,all=T)
			
			NbinhM<-with(subset(RVdata[[1]],year==years[i]&FLEN>=bins[j]&FLEN<bins[j+1]&FSEX==1),tapply(CLEN,SETNO,sum))
			NbinhM.dat<-data.frame(SETNO=as.numeric(names(NbinhM)),NbinhM)
			names(NbinhM.dat)[2]<-paste("bin",bins[j],sep='')
			survM.lst[[i]]<-merge(survM.lst[[i]],NbinhM.dat,all=T)
			
			Nbinh<-with(subset(RVdata[[1]],year==years[i]&FLEN>=bins[j]&FLEN<bins[j+1]),tapply(CLEN,SETNO,sum))
			Nbinh.dat<-data.frame(SETNO=as.numeric(names(Nbinh)),Nbinh)
			names(Nbinh.dat)[2]<-paste("bin",bins[j],sep='')
			surv.lst[[i]]<-merge(surv.lst[[i]],Nbinh.dat,all=T)
		}
		
		# Females
		survF.lst[[i]][is.na(survF.lst[[i]])]<-0
		stratF.lst[[i]]<-aggregate(survF.lst[[i]], list(survF.lst[[i]]$STRAT), mean)
		F.lst[[i]]<-colSums(sweep(stratF.lst[[i]][paste("bin",bins,sep='')],1,FUN='*',stratF.lst[[i]]$TUNITS))
		# Males
		survM.lst[[i]][is.na(survM.lst[[i]])]<-0
		stratM.lst[[i]]<-aggregate(survM.lst[[i]], list(survM.lst[[i]]$STRAT), mean)
		M.lst[[i]]<-colSums(sweep(stratM.lst[[i]][paste("bin",bins,sep='')],1,FUN='*',stratM.lst[[i]]$TUNITS))
		# Combined
		surv.lst[[i]][is.na(surv.lst[[i]])]<-0
		strat.lst[[i]]<-aggregate(surv.lst[[i]], list(surv.lst[[i]]$STRAT), mean)
		Total.lst[[i]]<-colSums(sweep(strat.lst[[i]][paste("bin",bins,sep='')],1,FUN='*',strat.lst[[i]]$TUNITS))
	}
	RVcatlenF<-do.call("rbind",F.lst)
	RVcatlenM<-do.call("rbind",M.lst)
	RVcatlenC<-do.call("rbind",Total.lst)
	Index<-data.frame(Year=years,N=N,B=B)
		
	list(Index=Index,RVcatlenF=RVcatlenF,RVcatlenM=RVcatlenM,RVcatlenC=RVcatlenC)
}

