channel <- odbcConnect("bank", username, password)

testQuery<-sqlQuery(channel,paste("SELECT year, SUM(units) units, FLEN, SUM(strclen) clen, SUM(strclen)/SUM(units) avgclen FROM (SELECT year, icld.strat, ROUND(floor(s.area/(1.75*(41/6080.2))+.5))  nits, flen, avg(clen) avgclen, variance(clen) varclen, s.area/((41./1000.0/6080.2)*1.75/1000.0)*AVG(CLEN)/1000000.0 strclen FROM gsstratum s, (SELECT year, icl.strat, icl.mission, icl.setno,  cl.flen, DECODE(NVL(sampwgt,0),0,0,totwgt/NVL(sampwgt,0)*nvl(clen*dncoef,0)*1.75/dist) clen FROM (SELECT mission,setno, 1+3*FLOOR(flen/3) flen, SUM(clen) clen, AVG(fwt) avg_fwt FROM gsdet WHERE  len IS NOT NULL AND spec=30 GROUP BY mission,setno, 1+3*FLOOR(flen/3) ) d, (SELECT year, mission, setno, strat, dist, dncoef, totwgt, sampwgt, flen FROM (SELECT class, flen FROM gs_lengths WHERE  lass=3 AND flen <= (SELECT max(flen) + 1 FROM gsdet WHERE spec=30 AND flen IS NOT NULL AND (mission, setno) IN (SELECT DISTINCT i.mission, i.setno FROM gsinf i, groundfish.gsmission_list l, gsmgt m WHERE i.mission=l.pk_mission AND i.strat=m.strat and ( m.unit in ('4VWX')) AND l.fk_series_id='SUMMER' AND l.year between 1970 and 2014 AND i.type=1))) l, (SELECT year, i.mission, i.setno, strat,  ist, dncoef, totwgt, sampwgt FROM (SELECT mission, setno, totwgt, sampwgt FROM gscat WHERE spec=30) c, (SELECT l.year, i.mission, i.setno, i.strat, dist, DECODE(CEIL((TO_NUMBER(TO_CHAR(i.sdate,' H24'))-(07-1))/(19-07)),1,1.0,1) dncoef FROM gsinf i, groundfish.gsmission_list l, gsmgt m WHERE i.mission=l.pk_mission AND i.strat=m.strat and ( m.unit in ('4VWX')) AND l.fk_series_id='SUMMER' AND  .year between 1970 and 2014 AND i.type=1) i WHERE i.mission=c.mission(+) AND i.setno=c.setno(+)) ic ) icl WHERE icl.mission=d.mission(+) AND icl.setno=d.setno(+) AND icl.flen=d.flen(+)) icld WHERE  .strat=icld.strat GROUP BY year, icld.strat, s.area, flen) iclds GROUP BY year, flen ORDER by year, flen"))
odbcCloseAll()
