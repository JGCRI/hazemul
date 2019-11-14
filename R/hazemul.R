library(fields)
library(RColorBrewer)

"%&%"<-function(x,y)paste0(x,y)
taper<-function(x,range=c(0.025,0.975)){
    xlims<-quantile(c(x),prob=range,na.rm=TRUE)
     x[x<xlims[1]]<-xlims[1]
     x[x>xlims[2]]<-xlims[2]
    return(x)}





load("Rdatasets/RData_GAT")

#allgat is a 3D array that contains yearly averages of global mean tas, for each year of historical/future simulations (in this case 1920-2100) (dimension 1), each ensemble member (in this case 10) (dimension 2) and each scenario (in this case 5) (dimension 3), note that one of the scenarios (position 3, "1pt5degC_OS" is not used


yearspan<-10  #10 year averages used to compute reference from which anomalies are defined

baseline<-seq(yearspan)
historical.offset<-75 #because we want our reference period to be 1995:2004

gatanom.rcp45<-allgat[,,"RCP45"]-mean(allgat[historical.offset+baseline,,"RCP45"]) #1.92
gatanom.rcp85<-allgat[,,"RCP85"]-mean(allgat[historical.offset+baseline,,"RCP85"]) #3.29

gatanom.1pt5<-allgat[,,"1pt5degC"]-mean(allgat[historical.offset+baseline,,"1pt5degC"])  #0.95
gatanom.2pt0<-allgat[,,"2pt0degC"]-mean(allgat[historical.offset+baseline,,"2pt0degC"])  #1.46


gatanom.rcp45<-apply(gatanom.rcp45,1,mean)  #averaging over ensemble members
gatanom.rcp85<-apply(gatanom.rcp85,1,mean)

gatanom.1pt5<-apply(gatanom.1pt5,1,mean)
gatanom.2pt0<-apply(gatanom.2pt0,1,mean)


#now compute target temperature at the end of the simulations in order to find corresponding years in other simulations used for time-shift
gatanom.rcp45.ave<-mean(rev(gatanom.rcp45)[-c(1:20)][1:yearspan])  #RCP45 ends at 2080, last 20 years are NAs, but that is just an oddity of this dataset
gatanom.1pt5.ave<-mean(rev(gatanom.1pt5)[1:yearspan])
gatanom.2pt0.ave<-mean(rev(gatanom.2pt0)[1:yearspan])
gatanom.rcp85.ave<-mean(rev(gatanom.rcp85)[1:yearspan])
gatanom.rcp85early.ave<-mean(rev(gatanom.rcp85)[-c(1:20)][1:yearspan])

#find years when origin temperature is closed to target temperature, for each pair of origin-target scenarios (note that the origin scenario is always "higher" than the target)
window<-which.min(abs(filter(gatanom.rcp45,filter=rep(1/yearspan,yearspan))-gatanom.1pt5.ave))
window.1pt5degC.RCP45<-(window-(yearspan/2-1)):(window+(yearspan/2))-historical.offset
window<-which.min(abs(filter(gatanom.rcp45,filter=rep(1/yearspan,yearspan))-gatanom.2pt0.ave))
window.2pt0degC.RCP45<-(window-(yearspan/2-1)):(window+(yearspan/2))-historical.offset

window<-which.min(abs(filter(gatanom.rcp85,filter=rep(1/yearspan,yearspan))-gatanom.1pt5.ave))   #these are referenced to 1920
window.1pt5degC.RCP85<-(window-(yearspan/2-1)):(window+(yearspan/2))-historical.offset   #these are referenced to 1995
window<-which.min(abs(filter(gatanom.rcp85,filter=rep(1/yearspan,yearspan))-gatanom.2pt0.ave))
window.2pt0degC.RCP85<-(window-(yearspan/2-1)):(window+(yearspan/2))-historical.offset

window<-which.min(abs(filter(gatanom.rcp85,filter=rep(1/yearspan,yearspan))-gatanom.rcp45.ave))   #these are referenced to 1920
window.RCP45.RCP85<-(window-(yearspan/2-1)):(window+(yearspan/2))-historical.offset   #these are referenced to 1995

window<-which.min(abs(filter(gatanom.2pt0,filter=rep(1/yearspan,yearspan))-gatanom.1pt5.ave))   #these are referenced to 1920
window.1pt5degC.2pt0degC<-(window-(yearspan/2-1)):(window+(yearspan/2))-historical.offset   #these are referenced to 1995


save(list=c("window.RCP45.RCP85","window.1pt5degC.2pt0degC","window.1pt5degC.RCP45","window.1pt5degC.RCP85","window.2pt0degC.RCP45","window.2pt0degC.RCP85"),file="Rdatasets/RData_10yrwindows_timeshift")  #these are all referenced to 1995 and the reason is that the arrays containing the extreme indices start at 1995.



#pattern scaling

for(index in c("txx","tnn","wsdi","gsl","fd","cdd","sdii","r10mm","r95ptot","rx5day")){
    attach("Rdatasets/RData_allscenarios2_"%&%index)  #these R archives contain arrays by index/scenario and two objects with coordinates (lat/lon). each array dimensions are length(lon) by length(lat) by number of years (1995:2100) for all scenarios but RCP45, which ends at 2080
    lat<-get("lat",pos=2)
    lon<-get("lon",pos=2)
    lon2<-lon
    lon2[lon>180]<-lon2[lon>180]-360
    oo<-order(lon2)  #for plotting
    lon2<-sort(lon2)



    for(originscenario in c("RCP45","RCP85")){
        for(targetscenario in c("1pt5degC","2pt0degC")){
            index.target<-get(paste(index,targetscenario,sep="."),pos=2)
            index.origin<-get(paste(index,originscenario,sep="."),pos=2)
            allgat.origin<-allgat[-seq(historical.offset),,originscenario][1:ifelse(originscenario=="RCP45",86,106),]
            allgat.target<-apply(allgat[-seq(historical.offset),,targetscenario][(ifelse(index=="gsl",105,106)-yearspan+1):ifelse(index=="gsl",105,106),],2,mean)  #this is because "gsl" is computed over a year that starts in July in the Southern hemisphere so the last half year in the arrays is full of NAs. for lat < 0

            index.baseline<-apply(index.target[,,seq(yearspan),],c(1,2,4),mean)


            pattern<-array(dim=c(2,288,192,10))
            for(em in 1:10)for(j in 1:192)for(i in 1:288)if(!all(is.na(index.origin[i,j,,em])))pattern[,i,j,em]<-lm(index.origin[i,j,,em]~allgat.origin[,em])$coef  #sadly this is done by gridpoint with a big loop
            index.emulated<-array(dim=c(288,192,10))
            for(em in 1:10)for(j in 1:192)for(i in 1:288)if(!all(is.na(index.origin[i,j,,em])))index.emulated[i,j,em]<-pattern[1,i,j,em]+pattern[2,i,j,em]*allgat.target[em]
            index.change.emulated<-index.emulated-index.baseline

            index.change.truth<-apply(index.target[,,(ifelse(index=="gsl",105,106)-yearspan+1):ifelse(index=="gsl",105,106),],c(1,2,4),mean)-index.baseline

#computing the three components of the error metric
            iv.truth<-apply(index.change.truth,c(1,2),var)
            iv.emulation<-apply(index.change.emulated, c(1,2),var)
            systemerror<-(apply(index.change.truth,c(1,2),mean)-apply(index.change.emulated,c(1,2),mean))^2

            totalmse<-iv.truth+iv.emulation+systemerror

            iv.truth.component<-100*(taper(iv.truth)/taper(totalmse))
            iv.emulation.component<-100*(taper(iv.emulation)/taper(totalmse))
            systemerror.component<-100*(taper(systemerror)/taper(totalmse))

            #compute the components at an aggregate level (global average)
            ww<-cos(lat/180*pi)
            ww<-t(matrix(ww,length(lat),length(lon)))
            ww[is.na(totalmse)]<-0

            totalmse.global<-sum((totalmse*ww)/sum(ww),na.rm=TRUE)
            global.iv.emulation.scaled<-sum((ww*iv.emulation/iv.truth)/sum(ww),na.rm=TRUE)
            global.systemerror.scaled<-sum((ww*systemerror/iv.truth)/sum(ww),na.rm=TRUE)


            save(list=c("lon2","lat","oo","index.change.emulated","index.change.truth","iv.truth","iv.emulation","systemerror","totalmse", "iv.truth.component","iv.emulation.component",
                        "systemerror.component","totalmse.global","global.iv.emulation.scaled","global.systemerror.scaled","pattern","allgat.target"),
                 file=paste("Rdatasets/RData",index,targetscenario,"SPS_from",originscenario,sep="_"))    #SPS = Simple Pattern Scaling :-)

        }
    }
    detach(2)
}




for(index in c("txx","tnn","wsdi","gsl","fd","cdd","sdii","r10mm","r95ptot","rx5day")){
    print(index)
    attach("Rdatasets/RData_allscenarios2_"%&%index)
    lat<-get("lat",pos=2)
    lon<-get("lon",pos=2)
    lon2<-lon
    lon2[lon>180]<-lon2[lon>180]-360
    oo<-order(lon2)
    lon2<-sort(lon2)


    originscenario<-"RCP85"
    targetscenario<-"RCP45"

    index.target<-get(paste(index,targetscenario,sep="."),pos=2)
    index.origin<-get(paste(index,originscenario,sep="."),pos=2)
    allgat.origin<-allgat[-seq(historical.offset),,originscenario][1:106,]
    allgat.target<-apply(allgat[-seq(historical.offset),,targetscenario][(ifelse(index=="gsl",85,86)-yearspan+1):ifelse(index=="gsl",85,86),],2,mean)

    index.baseline<-apply(index.target[,,seq(yearspan),],c(1,2,4),mean)


    pattern<-array(dim=c(2,288,192,10))
    for(em in 1:10)for(j in 1:192)for(i in 1:288)if(!all(is.na(index.origin[i,j,,em])))pattern[,i,j,em]<-lm(index.origin[i,j,,em]~allgat.origin[,em])$coef
    index.emulated<-array(dim=c(288,192,10))
    for(em in 1:10)for(j in 1:192)for(i in 1:288)if(!all(is.na(index.origin[i,j,,em])))index.emulated[i,j,em]<-pattern[1,i,j,em]+pattern[2,i,j,em]*allgat.target[em]
    index.change.emulated<-index.emulated-index.baseline

    index.change.truth<-apply(index.target[,,(ifelse(index=="gsl",85,86)-yearspan+1):ifelse(index=="gsl",85,86),],c(1,2,4),mean)-index.baseline


    iv.truth<-apply(index.change.truth,c(1,2),var)
    iv.emulation<-apply(index.change.emulated, c(1,2),var)
    systemerror<-(apply(index.change.truth,c(1,2),mean)-apply(index.change.emulated,c(1,2),mean))^2

    totalmse<-iv.truth+iv.emulation+systemerror

    iv.truth.component<-100*(taper(iv.truth)/taper(totalmse))
    iv.emulation.component<-100*(taper(iv.emulation)/taper(totalmse))
    systemerror.component<-100*(taper(systemerror)/taper(totalmse))

    ww<-cos(lat/180*pi)
    ww<-t(matrix(ww,length(lat),length(lon)))
    ww[is.na(totalmse)]<-0

    totalmse.global<-sum((totalmse*ww)/sum(ww),na.rm=TRUE)
    global.iv.emulation.scaled<-sum((ww*iv.emulation/iv.truth)/sum(ww),na.rm=TRUE)
    global.systemerror.scaled<-sum((ww*systemerror/iv.truth)/sum(ww),na.rm=TRUE)


    save(list=c("lon2","lat","oo","index.change.emulated","index.change.truth","iv.truth","iv.emulation","systemerror","totalmse", "iv.truth.component","iv.emulation.component",
                "systemerror.component","totalmse.global","global.iv.emulation.scaled","global.systemerror.scaled","pattern","allgat.target"),
         file=paste("Rdatasets/RData",index,targetscenario,"SPS_from",originscenario,sep="_"))
detach(2)}


for(index in c("txx","tnn","wsdi","gsl","fd","cdd","sdii","r10mm","r95ptot","rx5day")){
    print(index)
    attach("Rdatasets/RData_allscenarios2_"%&%index)
    lat<-get("lat",pos=2)
    lon<-get("lon",pos=2)
    lon2<-lon
    lon2[lon>180]<-lon2[lon>180]-360
    oo<-order(lon2)
    lon2<-sort(lon2)


    originscenario<-"2pt0degC"
    targetscenario<-"1pt5degC"

    index.target<-get(paste(index,targetscenario,sep="."),pos=2)
    index.origin<-get(paste(index,originscenario,sep="."),pos=2)
    allgat.origin<-allgat[-seq(historical.offset),,originscenario][1:106,]
    allgat.target<-apply(allgat[-seq(historical.offset),,targetscenario][(ifelse(index=="gsl",105,106)-yearspan+1):ifelse(index=="gsl",105,106),],2,mean)

    index.baseline<-apply(index.target[,,seq(yearspan),],c(1,2,4),mean)


    pattern<-array(dim=c(2,288,192,10))
    for(em in 1:10)for(j in 1:192)for(i in 1:288)if(!all(is.na(index.origin[i,j,,em])))pattern[,i,j,em]<-lm(index.origin[i,j,,em]~allgat.origin[,em])$coef
    index.emulated<-array(dim=c(288,192,10))
    for(em in 1:10)for(j in 1:192)for(i in 1:288)if(!all(is.na(index.origin[i,j,,em])))index.emulated[i,j,em]<-pattern[1,i,j,em]+pattern[2,i,j,em]*allgat.target[em]
    index.change.emulated<-index.emulated-index.baseline

    index.change.truth<-apply(index.target[,,(ifelse(index=="gsl",105,106)-yearspan+1):ifelse(index=="gsl",105,106),],c(1,2,4),mean)-index.baseline


    iv.truth<-apply(index.change.truth,c(1,2),var)
    iv.emulation<-apply(index.change.emulated, c(1,2),var)
    systemerror<-(apply(index.change.truth,c(1,2),mean)-apply(index.change.emulated,c(1,2),mean))^2

    totalmse<-iv.truth+iv.emulation+systemerror

    iv.truth.component<-100*(taper(iv.truth)/taper(totalmse))
    iv.emulation.component<-100*(taper(iv.emulation)/taper(totalmse))
    systemerror.component<-100*(taper(systemerror)/taper(totalmse))

    ww<-cos(lat/180*pi)
    ww<-t(matrix(ww,length(lat),length(lon)))
    ww[is.na(totalmse)]<-0

    totalmse.global<-sum((totalmse*ww)/sum(ww),na.rm=TRUE)
    global.iv.emulation.scaled<-sum((ww*iv.emulation/iv.truth)/sum(ww),na.rm=TRUE)
    global.systemerror.scaled<-sum((ww*systemerror/iv.truth)/sum(ww),na.rm=TRUE)


    save(list=c("lon2","lat","oo","index.change.emulated","index.change.truth","iv.truth","iv.emulation","systemerror","totalmse", "iv.truth.component","iv.emulation.component",
                "systemerror.component","totalmse.global","global.iv.emulation.scaled","global.systemerror.scaled","pattern","allgat.target"),
         file=paste("Rdatasets/RData",index,targetscenario,"SPS_from",originscenario,sep="_"))
detach(2)}






#Now time shift

####time shift

for(index in c("txx","tnn","wsdi","gsl","fd","cdd","sdii","r10mm","r95ptot","rx5day")){
    print(index)

    attach("Rdatasets/RData_allscenarios2_"%&%index)
    lat<-get("lat",pos=2)
    lon<-get("lon",pos=2)
    lon2<-lon
    lon2[lon>180]<-lon2[lon>180]-360
    oo<-order(lon2)
    lon2<-sort(lon2)



    for(originscenario in c("RCP45","RCP85")){
        print(originscenario)
        for(targetscenario in c("1pt5degC","2pt0degC")){
            print(targetscenario)
            index.target<-get(paste(index,targetscenario,sep="."),pos=2)
            index.origin<-get(paste(index,originscenario,sep="."),pos=2)

            index.baseline<-apply(index.target[,,seq(yearspan),],c(1,2,4),mean)

            index.emulated<-apply(index.origin[,,get(paste("window",targetscenario,originscenario,sep=".")),],c(1,2,4),mean)
            index.change.emulated<-index.emulated-index.baseline

            index.change.truth<-apply(index.target[,,(ifelse(index=="gsl",105,106)-yearspan+1):ifelse(index=="gsl",105,106),],c(1,2,4),mean)-index.baseline


            iv.truth<-apply(index.change.truth,c(1,2),var)
            iv.emulation<-apply(index.change.emulated, c(1,2),var)
            systemerror<-(apply(index.change.truth,c(1,2),mean)-apply(index.change.emulated,c(1,2),mean))^2

            totalmse<-iv.truth+iv.emulation+systemerror

            iv.truth.component<-100*(taper(iv.truth)/taper(totalmse))
            iv.emulation.component<-100*(taper(iv.emulation)/taper(totalmse))
            systemerror.component<-100*(taper(systemerror)/taper(totalmse))

            ww<-cos(lat/180*pi)
            ww<-t(matrix(ww,length(lat),length(lon)))
            ww[is.na(totalmse)]<-0

            totalmse.global<-sum((totalmse*ww)/sum(ww),na.rm=TRUE)
            global.iv.emulation.scaled<-sum((ww*iv.emulation/iv.truth)/sum(ww),na.rm=TRUE)
            global.systemerror.scaled<-sum((ww*systemerror/iv.truth)/sum(ww),na.rm=TRUE)


            save(list=c("lon2","lat","oo","index.change.emulated","index.change.truth","iv.truth","iv.emulation","systemerror","totalmse", "iv.truth.component","iv.emulation.component",
                        "systemerror.component","totalmse.global","global.iv.emulation.scaled","global.systemerror.scaled",paste("window",targetscenario,originscenario,sep=".")),
                 file=paste("Rdatasets/RData",index,targetscenario,"TS_from",originscenario,sep="_"))

        }
    }
    detach(2)
}






for(index in c("txx","tnn","wsdi","gsl","fd","cdd","sdii","r10mm","r95ptot","rx5day")){
    print(index)

    attach("Rdatasets/RData_allscenarios2_"%&%index)
    lat<-get("lat",pos=2)
    lon<-get("lon",pos=2)
    lon2<-lon
    lon2[lon>180]<-lon2[lon>180]-360
    oo<-order(lon2)
    lon2<-sort(lon2)



    for(originscenario in c("RCP85")){
        print(originscenario)
        for(targetscenario in c("RCP45")){
            print(targetscenario)
            index.target<-get(paste(index,targetscenario,sep="."),pos=2)
            index.origin<-get(paste(index,originscenario,sep="."),pos=2)

            index.baseline<-apply(index.target[,,seq(yearspan),],c(1,2,4),mean)

            index.emulated<-apply(index.origin[,,get(paste("window",targetscenario,originscenario,sep=".")),],c(1,2,4),mean)
            index.change.emulated<-index.emulated-index.baseline

            index.change.truth<-apply(index.target[,,(ifelse(index=="gsl",85,86)-yearspan+1):ifelse(index=="gsl",85,86),],c(1,2,4),mean)-index.baseline


            iv.truth<-apply(index.change.truth,c(1,2),var)
            iv.emulation<-apply(index.change.emulated, c(1,2),var)
            systemerror<-(apply(index.change.truth,c(1,2),mean)-apply(index.change.emulated,c(1,2),mean))^2

            totalmse<-iv.truth+iv.emulation+systemerror

            iv.truth.component<-100*(taper(iv.truth)/taper(totalmse))
            iv.emulation.component<-100*(taper(iv.emulation)/taper(totalmse))
            systemerror.component<-100*(taper(systemerror)/taper(totalmse))

            ww<-cos(lat/180*pi)
            ww<-t(matrix(ww,length(lat),length(lon)))
            ww[is.na(totalmse)]<-0

            totalmse.global<-sum((totalmse*ww)/sum(ww),na.rm=TRUE)
            global.iv.emulation.scaled<-sum((ww*iv.emulation/iv.truth)/sum(ww),na.rm=TRUE)
            global.systemerror.scaled<-sum((ww*systemerror/iv.truth)/sum(ww),na.rm=TRUE)


            save(list=c("lon2","lat","oo","index.change.emulated","index.change.truth","iv.truth","iv.emulation","systemerror","totalmse", "iv.truth.component","iv.emulation.component",
                        "systemerror.component","totalmse.global","global.iv.emulation.scaled","global.systemerror.scaled",paste("window",targetscenario,originscenario,sep=".")),
                 file=paste("Rdatasets/RData",index,targetscenario,"TS_from",originscenario,sep="_"))
}}
    detach(2)
}





for(index in c("txx","tnn","wsdi","gsl","fd","cdd","sdii","r10mm","r95ptot","rx5day")){
    print(index)

    attach("Rdatasets/RData_allscenarios2_"%&%index)
    lat<-get("lat",pos=2)
    lon<-get("lon",pos=2)
    lon2<-lon
    lon2[lon>180]<-lon2[lon>180]-360
    oo<-order(lon2)
    lon2<-sort(lon2)



    for(originscenario in c("2pt0degC")){
        print(originscenario)
        for(targetscenario in c("1pt5degC")){
            print(targetscenario)
            index.target<-get(paste(index,targetscenario,sep="."),pos=2)
            index.origin<-get(paste(index,originscenario,sep="."),pos=2)

            index.baseline<-apply(index.target[,,seq(yearspan),],c(1,2,4),mean)

            index.emulated<-apply(index.origin[,,get(paste("window",targetscenario,originscenario,sep=".")),],c(1,2,4),mean)
            index.change.emulated<-index.emulated-index.baseline

            index.change.truth<-apply(index.target[,,(ifelse(index=="gsl",105,106)-yearspan+1):ifelse(index=="gsl",105,106),],c(1,2,4),mean)-index.baseline


            iv.truth<-apply(index.change.truth,c(1,2),var)
            iv.emulation<-apply(index.change.emulated, c(1,2),var)
            systemerror<-(apply(index.change.truth,c(1,2),mean)-apply(index.change.emulated,c(1,2),mean))^2

            totalmse<-iv.truth+iv.emulation+systemerror

            iv.truth.component<-100*(taper(iv.truth)/taper(totalmse))
            iv.emulation.component<-100*(taper(iv.emulation)/taper(totalmse))
            systemerror.component<-100*(taper(systemerror)/taper(totalmse))

            ww<-cos(lat/180*pi)
            ww<-t(matrix(ww,length(lat),length(lon)))
            ww[is.na(totalmse)]<-0

            totalmse.global<-sum((totalmse*ww)/sum(ww),na.rm=TRUE)
            global.iv.emulation.scaled<-sum((ww*iv.emulation/iv.truth)/sum(ww),na.rm=TRUE)
            global.systemerror.scaled<-sum((ww*systemerror/iv.truth)/sum(ww),na.rm=TRUE)


            save(list=c("lon2","lat","oo","index.change.emulated","index.change.truth","iv.truth","iv.emulation","systemerror","totalmse", "iv.truth.component","iv.emulation.component",
                        "systemerror.component","totalmse.global","global.iv.emulation.scaled","global.systemerror.scaled",paste("window",targetscenario,originscenario,sep=".")),
                 file=paste("Rdatasets/RData",index,targetscenario,"TS_from",originscenario,sep="_"))

        }
    }
    detach(2)
}






