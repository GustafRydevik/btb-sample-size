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

levels(hw.results.df$gold.standard)<-c("300+/1000- GS tests","600+/2000- GS tests","900+/3000- GS tests")


hw.results.q20<-aggregate(hw.results.df[,2],by=hw.results.df[,c("scenario","vaccine.efficacy","gold.standard","samplesize","variable")],quantile,0.2)
hw.results.q20.SPDIVA<-subset(hw.results.q20,variable=="SpDIVA vaccine pop")
hw.results.q20.SPDIVA$success<-with(hw.results.q20.SPDIVA,ave(x,scenario,gold.standard,vaccine.efficacy,FUN=function(x)any(x>0.9985)))
hw.results.q20.SPDIVA$samplesize<-str_replace(sprintf("%d",hw.results.q20.SPDIVA$samplesize),"([[:digit:]]{3}$)"," \\1")
hw.results.q20.SPDIVA$samplesize<-factor(as.character(hw.results.q20.SPDIVA$samplesize),levels=unique(hw.results.q20.SPDIVA$samplesize))
ggplot(subset(hw.results.q20.SPDIVA),aes(x=samplesize,y=x,group=vaccine.efficacy,col=vaccine.efficacy))+
  geom_smooth()+facet_grid(scenario~gold.standard)+
  geom_hline(yintercept=0.9985)+
  scale_y_continuous(labels = percent_format())

