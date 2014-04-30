getadmout<-function(infile="c:\\mydocu~1\\sjh99-~1\\stocka~1\\output\\walleye5.rep"){
#
# Gets the formatted output from Admodel Builder (ask SJH) and brings
# it into R as a list of matrices and vectors
  x <- scan(infile, what="", sep="")
  x.elements <- is.na(as.numeric(x))
  x.dims <- seq(1, sum(x.elements) - 1, 2)
  x.labs <- seq(2, sum(x.elements), 2)
  y <- vector("list", length(x.labs))
  names(y) <- x[x.elements][x.labs]
  x.omit <- rep(F, length(x))
  for (j in 2:length(x)){
        if(x.elements[j-1] == T & x.elements[j] == T){
           x.omit[j-1] <- T
        }
  }
  x1 <- x[!x.omit]
  x1.elements <- is.na(as.numeric(x1))
  x1.breaks   <- c((1:length(x1))[x1.elements], c(length(x1)+1))
  print(length(x1.breaks))
  x1.dims <-  x[x.elements][x.dims]
  print(length(x1.dims))
  for(i in 2:length(x1.breaks)){
        y[[i-1]] <- as.numeric(x1[(x1.breaks[i-1] + 1):(x1.breaks[i] - 1)])     
        if(substring(x1.dims[i-1],1,1) == "r"){
                by.row <- T
        } else {
                by.row <- F
        }
        n.row <- as.numeric(substring(x1.dims[i-1],2))
        if(n.row > 1){
                y[[i-1]] <- matrix(y[[i-1]], byrow=by.row, nrow=n.row)
        }       
  }
  return(y)
}
