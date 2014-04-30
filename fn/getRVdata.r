getRVdata<-function(){
		
	require(RODBC) || stop("Package RODBC cannot be found")
	channel <- odbcConnect("bank", "hubleyb", "p35mghk")
	
	
##Brad's first cut
	a.flen = sqlQuery(channel,paste("select i.mission,i.setno,(dmin+dmax)/2 depth, strat, sdate, dist,flen,fwt,clen, fsex from groundfish.gsinf i, groundfish.gsdet c where i.mission=c.mission and i.setno=c.setno and spec=30 and to_char(sdate,'mm') in ('06','07','08') and strat between '440' and '495'  and type=1",sep=""))

	#total set info
	a = sqlQuery(channel,paste("select i.mission,i.setno,(dmin+dmax)/2 depth, strat, sdate, dist,totno,totwgt from groundfish.gsinf i, groundfish.gscat c where i.mission=c.mission and i.setno=c.setno and spec=30 and to_char(sdate,'mm') in ('06','07','08') and strat between '440' and '495'  and type=1",sep=""))
	b = sqlQuery(channel,paste("select i.mission,i.setno,(dmin+dmax)/2 depth, strat, sdate,dist, 0 totno, 0totwgt from groundfish.gsinf i where to_char(sdate,'mm') in ('06','07','08') and strat between '440' and '495' and type=1"))

	
###Nell's queries
	#NB missing 2004 and 2007 teleost surveys  ## this could be added right in by selecting 'SUMMER_TELEOST' in l.fk_series_id and only 2004 and 2007...I think the teleost and the needler dis a survey in 2005
#	a.flen = sqlQuery(channel,paste("select l.year, i.mission,i.setno,(dmin+dmax)/2 depth, i.strat, i.sdate, i.dist, c.flen, c.fwt, c.clen, c.fsex from groundfish.gsmission_list l, groundfish.gsinf i, groundfish.gsdet c,  groundfish.gsmgt m where i.mission=c.mission and i.mission=l.pk_mission and i.setno=c.setno and i.strat=m.strat and m.unit in ('4VWX') and l.fk_series_id='SUMMER' and c.spec=30 and i.type=1",sep=""))
	#dim(a.flen)

#	a = sqlQuery(channel,paste("select l.year, i.mission,i.setno,(dmin+dmax)/2 depth, i.strat, i.sdate, i.dist, c.totno, c.totwgt  from groundfish.gsmission_list l, groundfish.gsinf i, groundfish.gsmgt m, groundfish.gscat c where i.mission=c.mission and i.mission=l.pk_mission and i.setno=c.setno and i.strat=m.strat and m.unit in ('4VWX') and l.fk_series_id='SUMMER' and c.spec=30 and i.type=1",sep=""))
	#dim(a)

	##this is a work around an outer join on catches
#	b = sqlQuery(channel,paste("select l.year, i.mission,i.setno,(dmin+dmax)/2 depth, i.strat, i.sdate, i.dist, 0 totno, 0 totwgt  from groundfish.gsmission_list l, groundfish.gsinf i, groundfish.gsmgt m where  i.mission=l.pk_mission  and i.strat=m.strat and m.unit in ('4VWX') and l.fk_series_id='SUMMER' and i.type=1",sep="")) 
	
	a.flen$FSEX[is.na(a.flen$FSEX)]<-0
	
	d = sqlQuery(channel,paste("select strat,area from groundfish.gsstratum")) #strata areas area is in nm^2
	
	a$TOTWGT = a$TOTWGT*1.75/a$DIST
	a$TOTNO = a$TOTNO*1.75/a$DIST
	
	e = rbind(a,b)
	e = e[!duplicated(e[,c('MISSION','SETNO')]),]
	
	e$DEPTH = e$DEPTH*1.8288 #depth in m
	
	d1<-d
	d$TUNITS <- d$AREA/((35./6080.2)*1.75)
	d1$TUNITS <- d1$AREA/((41./6080.2)*1.75)
		
	e1<-merge(subset(e,as.Date(SDATE)<as.Date("1981-01-01")),d,by='STRAT')
	e2<-merge(subset(e,as.Date(SDATE)>=as.Date("1981-01-01")),d1,by='STRAT')
	e = rbind(e1,e2)

	odbcCloseAll()
	
	list(a.flen,e)

}
	
