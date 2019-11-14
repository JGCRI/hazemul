#code to plot paper's figures


"%&%"<-function(x,y)paste0(x,y)

load("Rdatasets/RData_cooLE")
lon2<-lonLE
lon2[lonLE>180]<-lon2[lonLE>180]-360
oo<-order(lon2)
lon2<-sort(lon2)
lat<-latLE

library(RColorBrewer)
percentcol<-c(brewer.pal(9,"YlOrRd"),"darkorchid4")
varpalette<-c("yellowgreen",percentcol)

library(fields)

indices<-c("fd","wsdi","tnn","txx","gsl","cdd","sdii","r95ptot","r10mm","rx5day")
indices.percent<-c("fd","wsdi","gsl","cdd","sdii","r95ptot","r10mm","rx5day")

taper<-function(x,range=c(0.025,0.975)){
    xlims<-quantile(c(x),prob=range,na.rm=TRUE)
     x[x<xlims[1]]<-xlims[1]
     x[x>xlims[2]]<-xlims[2]
    return(x)}




                                        #only a couple of indices to show: SDII and WSDI

originscenario<-"RCP45"
for( targetscenario in c("1pt5degC")){

    for(index in indices){

        load(paste("Rdatasets/RData",index,targetscenario,"SPS_from",originscenario,sep="_"))

        truemean<-apply(index.change.truth,c(1,2),mean)
        truemean[,lat<(-60)]<-NA
        truevar<-iv.truth
        truevar[,lat<(-60)]<-NA


        spsmean<-apply(index.change.emulated,c(1,2),mean)
        spsmean[,lat<(-60)]<-NA
        spsvar<-iv.emulation
        spsvar[,lat<(-60)]<-NA

        load(paste("Rdatasets/RData",index,targetscenario,"TS_from",originscenario,sep="_"))
        tsmean<-apply(index.change.emulated,c(1,2),mean)
        tsmean[,lat<(-60)]<-NA
        tsvar<-iv.emulation
        tsvar[,lat<(-60)]<-NA

        tsmeandiff<-(truemean-tsmean)^2
        spsmeandiff<-(truemean-spsmean)^2


        spsmeandiff<-taper(spsmeandiff)
        truevar<-taper(truevar)
        spsvar<-taper(spsvar)
        spsvar[truevar==0]<-0
        spsmeandiff[truevar==0]<-0


        spstotalvar<-taper((spsmeandiff+spsvar)/truevar)

        tsmeandiff<-taper(tsmeandiff)
        tsvar<-taper(tsvar)
        tsmeandiff[truevar==0]<-0
        tsvar[truevar==0]<-0

        tstotalvar<-taper((tsmeandiff+tsvar)/truevar)




        col<-brewer.pal(11,"Spectral")
        varpalette<-rev(col[-c(1,2,11)])

        errpalette<-rev(col)[c(5:11,1)]

        jpeg(paste0("newpics3/",paste(index,"totalvariance_alexdec",targetscenario,"from",originscenario,sep="_"),".jpg"),quality=100,width=1800,height=1800)
        par(mfrow=c(2,2),mar=c(5.1,5.1, 4.1 ,15.1))

        temp<-taper(sqrt(spsmeandiff[oo,]/truevar[oo,]))
        temp[temp>4]<-4
        image(lon2,lat,temp,zlim=c(0,4),breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4),col=errpalette,
              main="Error of Approximation/Internal Variance",xlab="",ylab="",axes=F,cex.main=3,ylim=c(-60,90),axis.args=list(cex.axis=3),legend.width=4)
        map("world",add=TRUE,interior=F)
        temp<-taper(sqrt(tsmeandiff[oo,]/truevar[oo,]))
        temp[temp>4]<-4
        image.plot(lon2,lat,temp,zlim=c(0,4),breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4),col=errpalette,
              main="Error of Approximation/Internal Variance",xlab="",ylab="",axes=F,cex.main=3,ylim=c(-60,90),axis.args=list(cex.axis=3),legend.width=4)
        map("world",add=TRUE,interior=F)
        temp<-taper(sqrt(spsvar[oo,]/truevar[oo,]))
        temp[temp>2]<-2
        image(lon2,lat,temp,zlim=c(0,2),breaks=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2),col=varpalette,
                   main="Pattern Scaling Variance/Internal Variance",axes=F,xlab="",ylab="",cex.main=3,ylim=c(-60,90),axis.args=list(cex.axis=3),legend.width=4)
        map("world",add=TRUE,interior=F)
        temp<-taper(sqrt(tsvar[oo,]/truevar[oo,]))
        temp[temp>2]<-2
        image.plot(lon2,lat,temp,zlim=c(0,3),breaks=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2),col=varpalette,
                   main="Time-Shift Variance/Internal Variance",axes=F,xlab="",ylab="",cex.main=3,ylim=c(-60,90),
                   axis.args=list(cex.axis=3),legend.width=4)
        map("world",add=TRUE,interior=F)
        dev.off()
    }
}













allfractionalvariances<-array(dim=c(2,10,2,4))
dimnames(allfractionalvariances)<-list(c("1pt5degC","2pt0degC"),indices,c("ts","sps"),c("totalvar","meandiff","truevar","methodvar"))

originscenario<-"RCP45"

for( targetscenario in c("1pt5degC","2pt0degC")){

    for(index in indices){

        load(paste("Rdatasets/RData",index,targetscenario,"SPS_from",originscenario,sep="_"))

        truemean<-apply(index.change.truth,c(1,2),mean)
        truemean[,lat<(-60)]<-NA
        truevar<-iv.truth
        truevar[,lat<(-60)]<-NA


        spsmean<-apply(index.change.emulated,c(1,2),mean)
        spsmean[,lat<(-60)]<-NA
        spsvar<-iv.emulation
        spsvar[,lat<(-60)]<-NA

        load(paste("Rdatasets/RData",index,targetscenario,"TS_from",originscenario,sep="_"))
        tsmean<-apply(index.change.emulated,c(1,2),mean)
        tsmean[,lat<(-60)]<-NA
        tsvar<-iv.emulation
        tsvar[,lat<(-60)]<-NA

        tsmeandiff<-(truemean-tsmean)^2
        spsmeandiff<-(truemean-spsmean)^2


        spsmeandiff<-taper(spsmeandiff)
        truevar<-taper(truevar)
        spsvar<-taper(spsvar)

        spstotalvar<-spsmeandiff+spsvar+truevar

        tsmeandiff<-taper(tsmeandiff)
        tsvar<-taper(tsvar)

        tstotalvar<-tsmeandiff+tsvar+truevar

        ww<-wwLE
        ww[is.na(truevar)]<-0
        allfractionalvariances[targetscenario,index,,"truevar"]<-weighted.mean(sqrt(truevar),w=ww)

        ww<-wwLE
        temp<-tstotalvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","totalvar"]<-weighted.mean(sqrt(tstotalvar/truevar),w=ww)

        ww<-wwLE
        temp<-tsmeandiff/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","meandiff"]<-weighted.mean(sqrt(tsmeandiff/truevar),w=ww)

        ww<-wwLE
        temp<-tsvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","methodvar"]<-weighted.mean(sqrt(tsvar/truevar),w=ww)

        ww<-wwLE
        temp<-spstotalvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","totalvar"]<-weighted.mean(sqrt(spstotalvar/truevar),w=ww)

        ww<-wwLE
        temp<-spsmeandiff/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","meandiff"]<-weighted.mean(sqrt(spsmeandiff/truevar),w=ww)

        ww<-wwLE
        temp<-spsvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","methodvar"]<-weighted.mean(sqrt(spsvar/truevar),w=ww)
    }}

allfractionalvariances.RCP45<-allfractionalvariances





jpeg("newpics3/barplot_totalvariance_absolutevalues.jpg",quality=100,width=1200,height=1200)
par(mfrow=c(2,2),mar=c(7,7,7,2))

for(targetscen in c("1pt5degC","2pt0degC")){
    temp<-t(allfractionalvariances[targetscen,,"sps",-1])
    oo<-order(temp[1,])
    barplot(temp[c(3,1),oo],names.arg=indices[oo],legend.text=FALSE,beside=FALSE,col=brewer.pal(11,"Spectral")[10:11],main=targetscennames[targetscen]%&%"\n Simple Pattern Scaling",xlab="",cex.names=2,cex.main=2,xlim=c(0,12),ylim=c(0,4),las=2,args.legend=list(cex=2),cex.axis=2)
    segments(0,1,12.5,1,col=brewer.pal(11,"Spectral")[10],lwd=3,lty=2)

    temp<-t(allfractionalvariances[targetscen,,"ts",-1])
    oo<-order(temp[1,])
    barplot(temp[c(3,1),oo],names.arg=indices[oo],legend.text=FALSE,beside=FALSE,col=brewer.pal(11,"Spectral")[10:11],main=targetscennames[targetscen]%&%"\n Timeshift",xlab="",cex.names=2,cex.main=2,xlim=c(0,12),las=2,ylim=c(0,4),args.legend=list(cex=2),cex.axis=2)
    segments(0,1,12.5,1,col=brewer.pal(11,"Spectral")[10],lwd=3,lty=2)

}
dimnames(temp[c(2,3,1),])
dev.off()



allfractionalvariances<-array(dim=c(3,10,2,4))
dimnames(allfractionalvariances)<-list(c("1pt5degC","2pt0degC","RCP45"),indices,c("ts","sps"),c("totalvar","meandiff","truevar","methodvar"))


originscenario<-"RCP85"

for( targetscenario in c("1pt5degC","2pt0degC","RCP45")){

    for(index in indices){

        load(paste("Rdatasets/RData",index,targetscenario,"SPS_from",originscenario,sep="_"))

        truemean<-apply(index.change.truth,c(1,2),mean)
        truemean[,lat<(-60)]<-NA
        truevar<-iv.truth
        truevar[,lat<(-60)]<-NA


        spsmean<-apply(index.change.emulated,c(1,2),mean)
        spsmean[,lat<(-60)]<-NA
        spsvar<-iv.emulation
        spsvar[,lat<(-60)]<-NA

        load(paste("Rdatasets/RData",index,targetscenario,"TS_from",originscenario,sep="_"))
        tsmean<-apply(index.change.emulated,c(1,2),mean)
        tsmean[,lat<(-60)]<-NA
        tsvar<-iv.emulation
        tsvar[,lat<(-60)]<-NA

        tsmeandiff<-(truemean-tsmean)^2
        spsmeandiff<-(truemean-spsmean)^2


        spsmeandiff<-taper(spsmeandiff)
        truevar<-taper(truevar)
        spsvar<-taper(spsvar)

        spstotalvar<-spsmeandiff+spsvar+truevar

        tsmeandiff<-taper(tsmeandiff)
        tsvar<-taper(tsvar)

        tstotalvar<-tsmeandiff+tsvar+truevar



        ww<-wwLE
        ww[is.na(truevar)]<-0
        allfractionalvariances[targetscenario,index,,"truevar"]<-weighted.mean(sqrt(truevar),w=ww)

        ww<-wwLE
        temp<-tstotalvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","totalvar"]<-weighted.mean(sqrt(tstotalvar/truevar),w=ww)

        ww<-wwLE
        temp<-tsmeandiff/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","meandiff"]<-weighted.mean(sqrt(tsmeandiff/truevar),w=ww)

        ww<-wwLE
        temp<-tsvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","methodvar"]<-weighted.mean(sqrt(tsvar/truevar),w=ww)

        ww<-wwLE
        temp<-spstotalvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","totalvar"]<-weighted.mean(sqrt(spstotalvar/truevar),w=ww)

        ww<-wwLE
        temp<-spsmeandiff/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","meandiff"]<-weighted.mean(sqrt(spsmeandiff/truevar),w=ww)

        ww<-wwLE
        temp<-spsvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","methodvar"]<-weighted.mean(sqrt(spsvar/truevar),w=ww)

    }}

allfractionalvariances.RCP85<-allfractionalvariances




allfractionalvariances<-array(dim=c(1,10,2,4))
dimnames(allfractionalvariances)<-list(c("1pt5degC"),indices,c("ts","sps"),c("totalvar","meandiff","truevar","methodvar"))


originscenario<-"2pt0degC"

for( targetscenario in c("1pt5degC")){

    for(index in indices){

        load(paste("Rdatasets/RData",index,targetscenario,"SPS_from",originscenario,sep="_"))

        truemean<-apply(index.change.truth,c(1,2),mean)
        truemean[,lat<(-60)]<-NA
        truevar<-iv.truth
        truevar[,lat<(-60)]<-NA


        spsmean<-apply(index.change.emulated,c(1,2),mean)
        spsmean[,lat<(-60)]<-NA
        spsvar<-iv.emulation
        spsvar[,lat<(-60)]<-NA

        load(paste("Rdatasets/RData",index,targetscenario,"TS_from",originscenario,sep="_"))
        tsmean<-apply(index.change.emulated,c(1,2),mean)
        tsmean[,lat<(-60)]<-NA
        tsvar<-iv.emulation
        tsvar[,lat<(-60)]<-NA

        tsmeandiff<-(truemean-tsmean)^2
        spsmeandiff<-(truemean-spsmean)^2


        spsmeandiff<-taper(spsmeandiff)
        truevar<-taper(truevar)
        spsvar<-taper(spsvar)

        spstotalvar<-spsmeandiff+spsvar+truevar

        tsmeandiff<-taper(tsmeandiff)
        tsvar<-taper(tsvar)

        tstotalvar<-tsmeandiff+tsvar+truevar



        ww<-wwLE
        ww[is.na(truevar)]<-0
        allfractionalvariances[targetscenario,index,,"truevar"]<-weighted.mean(sqrt(truevar),w=ww)

        ww<-wwLE
        temp<-tstotalvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","totalvar"]<-weighted.mean(sqrt(tstotalvar/truevar),w=ww)

        ww<-wwLE
        temp<-tsmeandiff/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","meandiff"]<-weighted.mean(sqrt(tsmeandiff/truevar),w=ww)

        ww<-wwLE
        temp<-tsvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"ts","methodvar"]<-weighted.mean(sqrt(tsvar/truevar),w=ww)

        ww<-wwLE
        temp<-spstotalvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","totalvar"]<-weighted.mean(sqrt(spstotalvar/truevar),w=ww)

        ww<-wwLE
        temp<-spsmeandiff/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","meandiff"]<-weighted.mean(sqrt(spsmeandiff/truevar),w=ww)

        ww<-wwLE
        temp<-spsvar/truevar
        ww[is.na(temp)|temp==Inf]<-0
        allfractionalvariances[targetscenario,index,"sps","methodvar"]<-weighted.mean(sqrt(spsvar/truevar),w=ww)

    }}

allfractionalvariances.2pt0degC<-allfractionalvariances






biaserror.from2pt0degC.rescaled<-allfractionalvariances.2pt0degC[,,,"meandiff"]

biaserror.fromRCP45.rescaled<-allfractionalvariances.RCP45[,,,"meandiff"]

biaserror.fromRCP85.rescaled<-allfractionalvariances.RCP85[,,,"meandiff"]

iverror.from2pt0degC.rescaled<-allfractionalvariances.2pt0degC[,,,"methodvar"]

iverror.fromRCP45.rescaled<-allfractionalvariances.RCP45[,,,"methodvar"]

iverror.fromRCP85.rescaled<-allfractionalvariances.RCP85[,,,"methodvar"]




gatdifferences<-c(gatanom.2pt0.ave-gatanom.1pt5.ave,gatanom.rcp45.ave-gatanom.2pt0.ave,gatanom.rcp45.ave-gatanom.1pt5.ave,
                  gatanom.rcp85early.ave-gatanom.rcp45.ave, gatanom.rcp85.ave-gatanom.2pt0.ave,gatanom.rcp85.ave-gatanom.1pt5.ave)
names(gatdifferences)<-c("2pt0-1pt5","RCP45-2pt0","RCP45-1pt5","RCP85-RCP45","RCP85-2pt0","RCP85-1pt5")

library(RColorBrewer)

indices.colors<-brewer.pal(11,"Spectral")[c(5:1,7:11)]

lateralmove<-0.02


origin.scenarios<-c("2.0C","RCP4.5","RCP4.5","RCP8.5","RCP8.5","RCP8.5")
target.scenarios<-c("1.5C","2.0C","1.5C","RCP4.5","2.0C","1.5C")


jpeg("newpics3/biaserrors_asfunctionof_GATgap.jpg",quality=100,width=1800,height=1200)
par(mar=c(10,7,7,2))
plot(gatdifferences,c(0,0,1.5,1.5,3.5,3.5),type="n",xlab="GSAT (degrees C)",cex.lab=2, main="Emulation error vs.  difference in GSAT \n between origin and target scenario by the end of the century", ylab="",las=1,axes=F,cex.main=3,xlim=c(0,4))
axis(1, at=seq(0,4),las=1,cex.axis=2)
axis(2, at=seq(0,3,by=0.5),las=1,cex.axis=2)
abline(h=seq(0,3,by=0.5),col="grey")
abline(h=1)
points(rep(gatdifferences["2pt0-1pt5"],10)-lateralmove,biaserror.from2pt0degC.rescaled[,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["2pt0-1pt5"],10)+lateralmove,biaserror.from2pt0degC.rescaled[,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP45-2pt0"],10)-lateralmove,biaserror.fromRCP45.rescaled["2pt0degC",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP45-2pt0"],10)+lateralmove,biaserror.fromRCP45.rescaled["2pt0degC",,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP45-1pt5"],10)-lateralmove,biaserror.fromRCP45.rescaled["1pt5degC",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP45-1pt5"],10)+lateralmove,biaserror.fromRCP45.rescaled["1pt5degC",,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP85-RCP45"],10)-lateralmove,biaserror.fromRCP85.rescaled["RCP45",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP85-RCP45"],10)+lateralmove,biaserror.fromRCP85.rescaled["RCP45",,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP85-2pt0"],10)-lateralmove,biaserror.fromRCP85.rescaled["2pt0degC",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP85-2pt0"],10)+lateralmove,biaserror.fromRCP85.rescaled["2pt0degC",,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP85-1pt5"],10)-lateralmove,biaserror.fromRCP85.rescaled["1pt5degC",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP85-1pt5"],10)+lateralmove,biaserror.fromRCP85.rescaled["1pt5degC",,"ts"],col=indices.colors,pch=18, cex=4)

text(c(0,gatdifferences),rep(2.75, 7), label=c("from:\nto:\ndiff:",
                                         paste(origin.scenarios, target.scenarios,paste0(round(gatdifferences,dig=1),"C"), sep="\n")),cex=1.5,adj=1)
legend(0.1,3.5,bty="n",pch=19,col=indices.colors,cex=2.5,legend=indices,ncol=10)
dev.off()








jpeg("newpics3/iverrors_asfunctionof_GATgap.jpg",quality=100,width=1800,height=1200)
par(mar=c(10,7,7,2))
plot(gatdifferences,c(0,0,1.5,1.5,3.5,3.5),type="n",xlab="GSAT (degrees C)",cex.lab=2, main="Emulation variability vs.  difference in GSAT \n between origin and target scenario by the end of the century", ylab="",las=1,axes=F,cex.main=3,xlim=c(0,4))
axis(1, at=seq(0,4),las=1,cex.axis=2)
axis(2, at=seq(0,3,by=0.5),las=1,cex.axis=2)
abline(h=seq(0,3,by=0.5),col="grey")
abline(h=1)
points(rep(gatdifferences["2pt0-1pt5"],10)-lateralmove,iverror.from2pt0degC.rescaled[,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["2pt0-1pt5"],10)+lateralmove,iverror.from2pt0degC.rescaled[,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP45-2pt0"],10)-lateralmove,iverror.fromRCP45.rescaled["2pt0degC",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP45-2pt0"],10)+lateralmove,iverror.fromRCP45.rescaled["2pt0degC",,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP45-1pt5"],10)-lateralmove,iverror.fromRCP45.rescaled["1pt5degC",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP45-1pt5"],10)+lateralmove,iverror.fromRCP45.rescaled["1pt5degC",,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP85-RCP45"],10)-lateralmove,iverror.fromRCP85.rescaled["RCP45",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP85-RCP45"],10)+lateralmove,iverror.fromRCP85.rescaled["RCP45",,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP85-2pt0"],10)-lateralmove,iverror.fromRCP85.rescaled["2pt0degC",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP85-2pt0"],10)+lateralmove,iverror.fromRCP85.rescaled["2pt0degC",,"ts"],col=indices.colors,pch=18, cex=4)

points(rep(gatdifferences["RCP85-1pt5"],10)-lateralmove,iverror.fromRCP85.rescaled["1pt5degC",,"sps"],col=indices.colors,pch=20, cex=4)
points(rep(gatdifferences["RCP85-1pt5"],10)+lateralmove,iverror.fromRCP85.rescaled["1pt5degC",,"ts"],col=indices.colors,pch=18, cex=4)

text(c(0,gatdifferences),rep(2.75, 7), label=c("from:\nto:\ndiff:",
                                         paste(origin.scenarios, target.scenarios,paste0(round(gatdifferences,dig=1),"C"), sep="\n")),cex=1.5,adj=1)
legend(0.1,3.5,bty="n",pch=19,col=indices.colors,cex=2.5,legend=indices,ncol=10)
dev.off()

scenarios<-c("RCP85","RCP45","2pt0degC","1pt5degC")
scenarionames<-c("RCP8.5","RCP4.5","2.0C","1.5C")
for(i in 1:2){
    origin<-scenarios[i]
    oname<-scenarionames[i]
    for(j in (i+1):4){
        target<-scenarios[j]
        tname<-scenarionames[j]
        jpeg("newpics3/scatterplots_"%&%origin%&%"_"%&%target%&%".jpg",quality=100,width=600,height=600)
        iverrorsps<-get("iverror.from"%&%origin%&%".rescaled")[target,,"sps"]
        biaserrorsps<-get("biaserror.from"%&%origin%&%".rescaled")[target,,"sps"]
        iverrorts<-get("iverror.from"%&%origin%&%".rescaled")[target,,"ts"]
        biaserrorts<-get("biaserror.from"%&%origin%&%".rescaled")[target,,"ts"]
        plot(iverrorsps,biaserrorsps,col=indices.colors,cex=3,pch=20,xlim=range(0,1.5),ylim=range(-0.1,3),
             xlab="IV of emulation", ylab="Emulation Error",main="From "%&%oname%&%" to "%&%tname,las=1,cex.main=2,cex.axis=2,cex.lab=1.5)
        segments(rep(-0.1,100),seq(0,1,length=100),rep(1.7,100),seq(0,1,length=100),col="grey"%&%seq(0,99))
     points(iverrorts,biaserrorts,col=indices.colors,cex=3,pch=18)
     points(iverrorsps,biaserrorsps,col=indices.colors,cex=3,pch=20)
        abline(h=c(0,1),v=1,lwd=2)
   legend(0,3.0,bty="n",cex=1.5,pch=19,col=indices.colors[c(1,6,2,7,3,8,4,9,5,10)],legend=indices[c(1,6,2,7,3,8,4,9,5,10)],ncol=5)
             dev.off()}}



for(i in 3){
    origin<-scenarios[i]
    oname<-scenarionames[i]
    for(j in 4){
        target<-scenarios[j]
        tname<-scenarionames[j]
        jpeg("newpics3/scatterplots_"%&%origin%&%"_"%&%target%&%".jpg",quality=100,width=600,height=600)
        iverrorsps<-get("iverror.from"%&%origin%&%".rescaled")[,"sps"]
        biaserrorsps<-get("biaserror.from"%&%origin%&%".rescaled")[,"sps"]
        iverrorts<-get("iverror.from"%&%origin%&%".rescaled")[,"ts"]
        biaserrorts<-get("biaserror.from"%&%origin%&%".rescaled")[,"ts"]
        plot(iverrorsps,biaserrorsps,col=indices.colors,cex=3,pch=20,xlim=range(0,1.5),ylim=range(-0.1,3),
             xlab="IV of emulation", ylab="Emulation Error",main="From "%&%oname%&%" to "%&%tname,las=1,cex.main=2,cex.axis=2,cex.lab=1.5)
        segments(rep(-0.1,100),seq(0,1,length=100),rep(1.7,100),seq(0,1,length=100),col="grey"%&%seq(0,99))
     points(iverrorts,biaserrorts,col=indices.colors,cex=3,pch=18)
     points(iverrorsps,biaserrorsps,col=indices.colors,cex=3,pch=20)
        abline(h=c(0,1),v=1,lwd=2)
   legend(0,3.0,bty="n",cex=1.5,pch=19,col=indices.colors[c(1,6,2,7,3,8,4,9,5,10)],legend=indices[c(1,6,2,7,3,8,4,9,5,10)],ncol=5)
             dev.off()}}









jpeg("newpics3/barplot_totalvariance_absolutevalues_RCP85.jpg",quality=100,width=1200,height=1200)
par(mfrow=c(2,2),mar=c(7,7,7,2))

for(targetscen in c("1pt5degC","2pt0degC")){
    temp<-t(allfractionalvariances.RCP85[targetscen,,"sps",-1])
    oo<-order(temp[1,])
    barplot(temp[c(3,1),oo],names.arg=indices[oo],legend.text=FALSE,beside=FALSE,col=brewer.pal(11,"Spectral")[10:11],main=targetscennames[targetscen]%&%"\n Simple Pattern Scaling",xlab="",cex.names=2,cex.main=2,xlim=c(0,12),ylim=c(0,4),las=2,args.legend=list(cex=2),cex.axis=2)
    segments(0,1,12.5,1,col=brewer.pal(11,"Spectral")[9],lwd=3,lty=2)

    temp<-t(allfractionalvariances.RCP85[targetscen,,"ts",-1])
    oo<-order(temp[1,])
    barplot(temp[c(3,1),oo],names.arg=indices[oo],legend.text=FALSE,beside=FALSE,col=brewer.pal(11,"Spectral")[10:11],main=targetscennames[targetscen]%&%"\n Timeshift",xlab="",cex.names=2,cex.main=2,xlim=c(0,12),las=2,ylim=c(0,4),args.legend=list(cex=2),cex.axis=2)
    segments(0,1,12.5,1,col=brewer.pal(11,"Spectral")[9],lwd=3,lty=2)
}

dev.off()

jpeg("newpics3/barplot_totalvariance_absolutevalues_RCP45.jpg",quality=100,width=1200,height=1200)
par(mfrow=c(2,2),mar=c(7,7,7,2))

for(targetscen in c("1pt5degC","2pt0degC")){
    temp<-t(allfractionalvariances.RCP45[targetscen,,"sps",-1])
    oo<-order(temp[1,])
    barplot(temp[c(3,1),oo],names.arg=indices[oo],legend.text=FALSE,beside=FALSE,col=brewer.pal(11,"Spectral")[10:11],main=targetscennames[targetscen]%&%"\n Simple Pattern Scaling",xlab="",cex.names=2,cex.main=2,xlim=c(0,12),ylim=c(0,4),las=2,args.legend=list(cex=2),cex.axis=2)
    segments(0,1,12.5,1,col=brewer.pal(11,"Spectral")[9],lwd=3,lty=2)

    temp<-t(allfractionalvariances.RCP45[targetscen,,"ts",-1])
    oo<-order(temp[1,])
    barplot(temp[c(3,1),oo],names.arg=indices[oo],legend.text=FALSE,beside=FALSE,col=brewer.pal(11,"Spectral")[10:11],main=targetscennames[targetscen]%&%"\n Timeshift",xlab="",cex.names=2,cex.main=2,xlim=c(0,12),las=2,ylim=c(0,4),args.legend=list(cex=2),cex.axis=2)
    segments(0,1,12.5,1,col=brewer.pal(11,"Spectral")[9],lwd=3,lty=2)
}

dev.off()


