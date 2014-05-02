build.admb.dat <- function() {
    # [build]s the [ADMB] [dat]a set
    lf$ll.m$sex <- 0
    lf$ll.f$sex <- 1
    lf$ot.unk$sex <- 1
    lf$ll.m$method <- 1
    lf$ll.f$method <- 1
    lf$ot.unk$method <- 2
    lf.combined <- rbind(lf$ll.m, lf$ll.f, lf$ot.unk)
    #sample <- ddply(ss.all, c("year", "sex", "gear.type", "nafo.group"), 
        #function(x) {
            #data.frame(sample.size = sum(x[, -c(1:4)]))
        #})

    #names(sample)[names(sample) == "gear.type"] <- "method"
    #sample$method <- as.character(sample$method)
    #sample[sample$method == "otter trawl", "method"] <- 2
    #sample[sample$method == "long line", "method"] <- 1
    #sample[sample$sex == 0, "sex"] <- 999
    #sample[sample$sex == 1, "sex"] <- 0
    #sample[sample$sex == 2, "sex"] <- 1
    #sample <- sample[-which(sample$sex == 999), ]
    #sample$method <- as.numeric(sample$method)
# but what do do about the unknown sample size?
# we have infered length-frequencies for male or female for it but do not know sample size
# for now just ignoring sample size
    
    #com.lf <- merge(sample, lf.combined)
    com.lf <- lf.combined

    columns.to.ignore <- (1:ncol(com.lf))[names(com.lf) %in% c("year", "nafo.group", "sex", "method")]
    com.lf.s1 <- aggregate(com.lf[, -columns.to.ignore], by = list(com.lf$year, 
        com.lf$method, com.lf$nafo.group), sum)
    names(com.lf.s1) <- c("year", "method", "nafo.group")

   # rearrange col.lf
   com.lf <- com.lf[,c(columns.to.ignore, (1:ncol(com.lf))[-columns.to.ignore])]

    write.table(com.lf, file = "output/com.lf.txt", sep = "\t", row.names = F, quote = F)
    write.table(com.lf.s1, file = "output/com.lf.s1.txt", sep = "\t", row.names = F, quote = F)
} 

