getmcout<-function(infile="c:\\4x\\mcout.dat",burn=0)
{

# Gets the formatted output from the mcmc interations of Admodel Builder and brings
# it into Splus as a list of matrices and vectors
  x <- read.table(infile)
#print(x)
# remove burn-in samples
  afterburn<-(burn+1):length(x[,1])
  x<- x[afterburn,]
  x[,length(x[1,])]<-as.numeric(as.vector(x[,length(x[1,])]))
  # get locations of var names 

  names.locs<-names1<-NULL 

  for (i in 1:length(x[1,]))
        { 
        if(  !is.numeric(as.vector(x[1,i]))) 
                { 
                names.locs<-c(names.locs,i) 
                names1<-c(names1,as.character(x[1,i]))
                } 
        } 
 print(names1)
 print(names.locs) 
  diff<-NULL

  for (i in 1:(length(names.locs)))
        {
        if(i!=length(names.locs)) diff[i]<-names.locs[i+1] - names.locs[i] - 1  
        else diff[i] <- length(x[1,]) - names.locs[i] 
        } 
print(diff)



  data.list<-list() 
  #data.list[i]<-NULL
  for (i in 1:length(diff))
        { 
        if(diff[i]==1) data.list[[i]]<- as.vector(x[,(names.locs[i]+1):(names.locs[i]+diff[i])])
        else {data.list[[i]]<- as.matrix(x[(names.locs[i]+1):(names.locs[i]+diff[i])])}
#       names(data.list[i])<-paste(as.character(names1[i]))
        }  
#       names(data.list)<-paste(as.character(names1))
        names(data.list)<-names1
 
  return(data.list)

}
