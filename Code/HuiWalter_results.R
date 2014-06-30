library(rjags)
library(batch)
library(stringr)
library(ggplot2)
library(scales)
#Global parameters
got.mdrive<-length(dir("M:"))>0
is.win<-grepl("w32",R.Version()$platform)
is.girion<-sum(grep("girion",system("hostname",intern=TRUE)))
dropbox.sac<-"C:/Documents and Settings/GRydevik/My Documents/Dropbox"

dropbox.bioss<-"D:\\Dropbox"
dropbox.osx<-"/Users/gustafrydevik/Dropbox"
girion.path<-"/home/gustaf"
#Giles.path<-..... ##Add your directory here.
dropbox.path<-c(dropbox.osx,dropbox.sac,dropbox.bioss)[got.mdrive+is.win+1]

main.path<-if(is.girion){girion.path}else{file.path(dropbox.path,"PhD folder")}

##Project specific parameters
project.path<-file.path(main.path,"btb-sample-size")

data.dir<-file.path(project.path,"Data")
script.dir<-file.path(project.path,"Code")
output.dir<-file.path(project.path,"Output")
sim.dir<-file.path(project.path,"simulation-output")
session.id<-"2014062023"
hw.results.list<-lapply(grep(session.id,dir(sim.dir,full.names=T),value=T),read.table,header=T)
hw.results.df<-do.call("rbind",hw.results.list)


##adding run variables
hw.results.df$samplesize<-rep(as.numeric(str_replace(str_extract(grep(session.id,dir(sim.dir),value=T),"samplesize([[:digit:]]+)"),"samplesize","")),
                              times=sapply(hw.results.list,nrow))
hw.results.df[hw.results.df$samplesize==1,"samplesize"]<-100000
hw.results.df$scenario<-as.factor(rep((str_replace(str_extract(grep(session.id,dir(sim.dir),value=T),"scenario_([[:alpha:]]+)"),"scenario_","")),
                              times=sapply(hw.results.list,nrow)))
hw.results.df$vaccine.efficacy<-as.factor(rep((str_replace(str_extract(grep(session.id,dir(sim.dir),value=T),"vaccEff0\\.([[:alnum:]]+)"),"vaccEff","")),
                                      times=sapply(hw.results.list,nrow)))
hw.results.df$gold.standard<-as.factor(paste(rep(str_extract(grep(session.id,dir(sim.dir),value=T),"nposGold([[:digit:]]+)"),
                                              times=sapply(hw.results.list,nrow)),
                                              rep(str_extract(grep(session.id,dir(sim.dir),value=T),"nnegGold([[:digit:]]+)"),
                                                  times=sapply(hw.results.list,nrow)),sep=";"))
levels(hw.results.df$scenario)<-c("70%se, 99.90%sp","85%se, 99.90%sp","70%se, 99.99%sp","50%se, 99.90%sp","70%se, 99.87%sp","50%se, 99.87%sp","85%se, 99.99%sp")
hw.results.df$scenario<-factor(hw.results.df$scenario,levels(hw.results.df$scenario)[c(1,4,5,6,2,3,7)])
levels(hw.results.df$gold.standard)<-c("300+/1000- GS tests","600+/2000- GS tests","900+/3000- GS tests")


hw.results.q20<-aggregate(hw.results.df[,2],by=hw.results.df[,c("scenario","vaccine.efficacy","gold.standard","samplesize","variable")],quantile,0.2)
hw.results.q20.SPDIVA<-subset(hw.results.q20,variable=="SpDIVA vaccine pop")
hw.results.q20.SPDIVA$success<-with(hw.results.q20.SPDIVA,ave(x,scenario,gold.standard,vaccine.efficacy,FUN=function(x)any(x>0.9985)))
hw.results.q20.SPDIVA$samplesize<-factor(str_replace(sprintf("%d",hw.results.q20.SPDIVA$samplesize),"([[:digit:]]{3}$)","k"))
hw.results.q20.SPDIVA$samplesize<-factor(hw.results.q20.SPDIVA$samplesize,levels(hw.results.q20.SPDIVA$samplesize)[c(2,4,1,3)])

#hw.results.q20.SPDIVA$samplesize<-str_replace(sprintf("%d",hw.results.q20.SPDIVA$samplesize),"([[:digit:]]{3}$)"," \\1")
#hw.results.q20.SPDIVA$samplesize<-factor(as.character(hw.results.q20.SPDIVA$samplesize),levels=unique(hw.results.q20.SPDIVA$samplesize))


library(plyr)
hw.results.df$sp<-str_extract(as.character(hw.results.df$scenario),"[[:digit:]]+\\.[[:digit:]]+\\%sp")
hw.results.df$se<-str_extract(as.character(hw.results.df$scenario),"[[:digit:]]+\\.?[[:digit:]]*\\%se")

hw.results.plyr<-ddply(hw.results.df,.(se,sp,vaccine.efficacy,gold.standard,samplesize,variable),.fun=function(x)c(lower.ci.ave=mean(x$lower),upper.ci.ave=mean(x$upper)))
hw.results.plyr.70se<-subset(hw.results.plyr,se=="70%se"&variable=="SpDIVA vaccine pop"&sp!="99.87%sp")
hw.table.ci<-hw.results.plyr.70se[,-c(1,6)]
hw.table.ci$lower.ci.ave<-round(hw.table.ci$lower.ci.ave*100,3)
hw.table.ci$upper.ci.ave<-pmin(round(hw.table.ci$upper.ci.ave*100,3),100)
hw.table.ci<-hw.table.ci[order(hw.table.ci$samplesize,hw.table.ci$sp,hw.table.ci$gold.standard,hw.table.ci$vaccine.efficacy),]
hw.table.ci$ci.width<-hw.table.ci$upper.ci.ave-hw.table.ci$lower.ci.ave
hw.table.ci$Diff.99.85<-ifelse(hw.table.ci$lower.ci.ave>=99.85,"*","")
hw.table.ci$samplesize<-str_replace(sprintf("%d",hw.table.ci$samplesize),"000$"," 000")
hw.table.ci$sp<-str_replace(hw.table.ci$sp,"sp","")
names(hw.table.ci)[1]<-"True SP"
names(hw.table.ci)[8]<-"Significantly>99.85%?"

print(xtable(hw.table.ci[,c("samplesize",
               "True SP",
               "gold.standard",
               "vaccine.efficacy",
               "lower.ci.ave","upper.ci.ave","ci.width","Significantly>99.85%?")],digits=3),
      ,type="html",
      file=file.path(output.dir,"SampleSize_table.html"),
      include.rownames=FALSE)


my.theme<-theme_minimal()+theme(text=element_text(size=20))
ggplot(subset(hw.results.q20.SPDIVA),aes(x=samplesize,y=x,group=vaccine.efficacy,col=vaccine.efficacy))+
  geom_smooth(size=1.2)+facet_grid(gold.standard~scenario)+
  geom_hline(yintercept=0.9985)+
  scale_y_continuous(labels = percent_format())+my.theme



#### Sensitivity run
session.id<-"2014062613"
file.names.full<-grep(session.id,dir(sim.dir,full.names=T),value=T)
sp.results.list<-lapply(grep(session.id,dir(sim.dir,full.names=T),value=T),read.table,header=T)

par.samplesize<-sprintf("%d",
        as.numeric(str_replace(str_extract(grep(session.id,dir(sim.dir),value=T),
                                           "samplesize([[:digit:]]+)"),"samplesize",""
                               )
        ))

par.sp<-as.numeric(
  str_replace(str_extract(
    grep(session.id,dir(sim.dir),value=T),"sp0\\.([[:digit:]]+)"),"sp","")
)


par.gs<-as.factor(paste(str_extract(grep(session.id,dir(sim.dir),value=T),"nposGold([[:digit:]]+)"),
                str_extract(grep(session.id,dir(sim.dir),value=T),"nnegGold([[:digit:]]+)"),
                    sep=";")
          )

levels(par.gs)<-c("1500+/5000- GS tests","30+/100- GS tests","300+/1000- GS tests")

sp.results.list<-lapply(1:length(sp.results.list),function(x){
  data.frame(sp.results.list[[x]],
             samplesize=par.samplesize[x],
             gold.standard=par.gs[x],
             sp=par.sp[x])
})

sp.results.df<-do.call("rbind",sp.results.list)

##adding run variables

sp.results.q20<-aggregate(sp.results.df[,2],by=sp.results.df[,c("sp","gold.standard","samplesize","variable")],quantile,0.2)
sp.results.q20$samplesize<-as.numeric(as.character(sp.results.q20$samplesize))
#sp.results.q20$samplesize<-factor(str_replace(as.character(sp.results.q20$samplesize),"([[:digit:]]{3}$)","k"))
#sp.results.q20$samplesize<-factor(sp.results.q20$samplesize,levels(sp.results.q20$samplesize)[c(2:8,1)])


my.theme<-theme_minimal()+theme(text=element_text(size=20))
ggplot(subset(sp.results.q20,variable=="SpDIVA vaccine pop"),aes(x=samplesize,y=x,group=vaccine.efficacy))+
  geom_smooth(size=1.5)+facet_grid(gold.standard~sp)+
  geom_hline(yintercept=0.9985)+
  scale_y_continuous(labels = percent_format())+my.theme


