 setwd("C:/Halibut/Brad")

 source("fn/getRVdata.r")
 RVdata<-getRVdata()

 source("fn/prepRVdata.r")
 RVinput<-prepRVdata(RVdata)
 
 source("fn/write.dat.R")
 write.dat(RVinput[3],"RVcatlen.dat")
 write.table(RVinput$RVcatlenC, file = "RVcatlen.dat", sep = "\t",row.names = F, quote = F)

 
 
 source("fn/read.admb.r")
 output<-read.admb("admb/vpop10")
  
 
 setwd("C:/Halibut/halibut admb model")
 #system("vpop10")
 
 
	source("README.KT.r")
	source("fn/plot.vpop10.r")
	plot.vpop10()
	

	# test do you see me?
	
 setwd("C:/Halibut/Brad")
	
#system("caa1")
 source("fn/getadmout.r")
output<-getadmout("admb/caa1.rep")
 	source("fn/plot.caa1.r")
	plot.caa1()

	