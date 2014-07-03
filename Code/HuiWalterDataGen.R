
#Code for testing HuiWalter estimation on generated data.

library(batch)
library(rjags)
library(xtable)
library(stringr)
#Global parameters
got.mdrive<-length(dir("M:"))>0
is.win<-grepl("w32",R.Version()$platform)
is.girion<-sum(grep("girion",system("hostname",intern=TRUE)))
dropbox.sac<-"C:/Documents and Settings/GRydevik/My Documents/Dropbox"

dropbox.bioss<-"D:\\Dropbox"
dropbox.osx<-"/Users/gustafrydevik/Dropbox"
girion.path<-"/home/gustaf"
#Giles.path<-..... ##Add your 
dropbox.path<-c(dropbox.osx,dropbox.sac,dropbox.bioss)[got.mdrive+is.win+1]

main.path<-if(is.girion){girion.path}else{file.path(dropbox.path,"PhD folder")}

##Project specific parameters
project.path<-file.path(main.path,"btb-sample-size")

data.dir<-file.path(project.path,"Data")
script.dir<-file.path(project.path,"Code")
output.dir<-file.path(project.path,"Output")
sim.dir<-file.path(project.path,"simulation-output")

###Parameters for the run
#Pars for jags
nreps<-1
rep.prefix=0
nchains<-5
niter<-1000
n.mcmc.samples<-1000
save.samples=T

##Pars for the scenario
scenario.name<-"baseline"
SeStd<-c(Vacc=0.7,nonVacc=0.7)
SeDIVA<-c(Vacc=0.7,nonVacc=0.7)
SpStd<-c(Vacc=0.999,nonVacc=0.999)  ## a specificity of 0.5 mirrors that it reacts to vaccinated animals
SpDIVA<-c(Vacc=0.999,nonVacc=0.999)
vaccine.efficacy<-0.6
Prevalence<-c(High=6/100,Low=1/100)
samplesize<-40000
props<-c(1/4,1/4,1/4,1/4)# balance between populations High/vacc,Low/vacc,High/nonvacc,low/nonvacc
gold.scale<-1/10
npos.gold<-300*c(Vacc=1/2,nonVacc=1/2)*gold.scale
nneg.gold<-1000*c(Vacc=1/2,nonVacc=1/2)*gold.scale

time<-as.integer(as.numeric(format(Sys.time(),"%Y%m%d%H")))
seed<-1000

##Or read in parameters from a batch call
parseCommandArgs()

names(SeStd)<-c("Vacc","nonVacc")
names(SeDIVA)<-c("Vacc","nonVacc")
names(SpStd)<-c("Vacc","nonVacc")
names(SpDIVA)<-c("Vacc","nonVacc")
names(Prevalence)<-c("High","Low")
names(npos.gold)<-c("Vacc","nonVacc")
names(nneg.gold)<-c("Vacc","nonVacc")

set.seed(seed)

#Set up the data objects
nCattle<-matrix(c(nHigh.vacc=round(samplesize*props[1]),
                  nLow.vacc=round(samplesize*props[2]),
                  nHigh.nonvacc=round(samplesize*props[3]),
                  nLow.nonvacc=round(samplesize*props[4])),ncol=2)
btb.testdata<-data.frame(Pop=c(rep("High",sum(nCattle[1:2])),rep("Low",sum(nCattle[1:2]))))

btb.testdata$Vaccinated<-c(rep("Vacc",nCattle[1]),
                           rep("nonVacc",nCattle[2]),    
                           rep("Vacc",nCattle[3]),
                           rep("nonVacc",nCattle[4])
)
btb.testdata$Vaccinated<-factor(btb.testdata$Vaccinated,levels=c("Vacc","nonVacc"))


#Now run the following $reps times to generate multiple estimates

set.seed(seed)
sp.est<-vector(length=nreps)
posterior.intervals.all<-data.frame()
for(i in 1:nreps){
  
  #simulate some data
  btb.testdata$Case[btb.testdata$Pop=="High"&btb.testdata$Vaccinated=="Vacc"]<-rbinom(nCattle[1],1,p=Prevalence[1]*vaccine.efficacy)
  btb.testdata$Case[btb.testdata$Pop=="Low"&btb.testdata$Vaccinated=="Vacc"]<-rbinom(nCattle[2],1,p=Prevalence[2]*vaccine.efficacy)
  btb.testdata$Case[btb.testdata$Pop=="High"&btb.testdata$Vaccinated=="nonVacc"]<-rbinom(nCattle[3],1,p=Prevalence[1])
  btb.testdata$Case[btb.testdata$Pop=="Low"&btb.testdata$Vaccinated=="nonVacc"]<-rbinom(nCattle[4],1,p=Prevalence[2])
  
  for(vacc.status in c("Vacc","nonVacc")){
    SubsetPos<-with(btb.testdata,Case==1&Vaccinated==vacc.status)
    SubsetNeg<-with(btb.testdata,Case==0&Vaccinated==vacc.status)
    
    btb.testdata$Std[SubsetPos]<-rbinom(sum(SubsetPos),1,SeStd[vacc.status])
    btb.testdata$Std[SubsetNeg]<-rbinom(sum(SubsetNeg),1,1-SpStd[vacc.status])
    btb.testdata$DIVA[SubsetPos]<-rbinom(sum(SubsetPos),1,SeDIVA[vacc.status])
    btb.testdata$DIVA[SubsetNeg]<-rbinom(sum(SubsetNeg),1,1-SpDIVA[vacc.status])
    
  }
  
  
  
  btb.testdata.agg<-as.array(table(btb.testdata[,c("Pop","Vaccinated","Std","DIVA")]))
  btb.dimnames<-c(dimnames(btb.testdata.agg)[1:2],Std.vs.DIVA=list(c("--","+-","-+","++")))
  dim(btb.testdata.agg)<-c(2,2,4)
  dimnames(btb.testdata.agg)<-btb.dimnames
  class(btb.testdata.agg)<-"array"
  
  ###Gold standard data 
  SeDIVA.gold<-c(Vacc=rbinom(1,size=npos.gold["Vacc"],prob=SeDIVA["Vacc"]),
                 nonVacc=rbinom(1,size=npos.gold["nonVacc"],prob=SeDIVA["nonVacc"]))
  SeStd.gold<-c(Vacc=rbinom(1,size=npos.gold["Vacc"],prob=SeStd["Vacc"]),
                nonVacc=rbinom(1,size=npos.gold["nonVacc"],prob=SeStd["nonVacc"]))
  SpDIVA.gold<-c(Vacc=rbinom(1,size=nneg.gold["Vacc"],prob=SpDIVA["Vacc"]),
                 nonVacc=rbinom(1,size=nneg.gold["nonVacc"],prob=SpDIVA["nonVacc"]))
  SpStd.gold<-c(Vacc=rbinom(1,size=nneg.gold["Vacc"],prob=SpStd["Vacc"]),
                nonVacc=rbinom(1,size=nneg.gold["nonVacc"],prob=SpStd["nonVacc"]))
  
  huiwalter.jags<-jags.model(file.path(script.dir,"HuiWalterMultinomial_vaccinePop.txt"),
                             data=list(
                               counts=btb.testdata.agg,
                               nCattle=nCattle,
                               SeDIVA.gold=SeDIVA.gold,
                               SeStd.gold=SeStd.gold,
                               SpDIVA.gold=SpDIVA.gold,
                               SpStd.gold=SpStd.gold,
                               npos.gold=npos.gold,
                               nneg.gold=nneg.gold),
                             n.chains=nchains)
  
  update(huiwalter.jags,niter)
  huiwalter.samples<-coda.samples(huiwalter.jags,c("SeStd","SeDIVA","SpStd","SpDIVA","Prevalence","probResult","v.efficacy"),n.mcmc.samples)
  huiwalter.allchains<-do.call("rbind",huiwalter.samples)
  class(huiwalter.allchains)<-"mcmc"
  
  
  ###Generate posterior estimates
  posterior.intervals<-HPDinterval(as.mcmc(huiwalter.allchains[1:nrow(huiwalter.allchains),
                                                               which(str_detect(colnames(huiwalter.allchains),
                                                                                c("SeStd|SeDIVA|SpStd|SpDIVA|Prevalence|v.efficacy")))]))
  posterior.intervals<-data.frame(posterior.intervals)
  posterior.intervals$median<-sapply(data.frame(huiwalter.allchains[1:nrow(huiwalter.allchains),
                                                                    which(str_detect(colnames(huiwalter.allchains),
                                                                                     c("SeStd|SeDIVA|SpStd|SpDIVA|Prevalence|v.efficacy")))]
  ),
  median)
  posterior.intervals$mean<-sapply(data.frame(huiwalter.allchains[1:nrow(huiwalter.allchains),
                                                                  which(str_detect(colnames(huiwalter.allchains),
                                                                                   c("SeStd|SeDIVA|SpStd|SpDIVA|Prevalence|v.efficacy")))]
  ),
  mean)
  rownames(posterior.intervals)<-gsub("\\[1\\]"," vaccine pop",rownames(posterior.intervals))
  rownames(posterior.intervals)<-gsub("\\[2\\]"," nonvaccine pop",rownames(posterior.intervals))
  rownames(posterior.intervals)[1:2]<-c("High prevalence pop","Low prevalence pop")
  posterior.intervals<-data.frame(variable=rownames(posterior.intervals),posterior.intervals)
  posterior.intervals[,"True.value"]<-c(
    Prevalence,
    SeStd,
    SeDIVA,
    SpStd,  
    SpDIVA,
    vaccine.efficacy)
  posterior.intervals$samplesize<-samplesize
  posterior.intervals$spDIVA.vacc<-SpDIVA[1]
  posterior.intervals$gold.standard<-paste(sum(npos.gold),"+/",sum(nneg.gold),"-",collapse="",sep="")
  posterior.intervals$iteration<-i
  posterior.intervals.all<-rbind(posterior.intervals.all,posterior.intervals)
  #print(xtable(posterior.intervals,digits=4),type="html",file= "PosteriorIntervals.html")
}

  write.table(posterior.intervals.all,file=paste(
    sim.dir,"/HuiWalterCI",
    "_time_",time,
    "_sp",SpDIVA[1],
    "_se",SeDIVA[1],
    "_samplesize",sprintf("%d",samplesize),
    "_nposGold",sum(npos.gold),
    "_nnegGold",sum(nneg.gold),
    "_vaccEff",vaccine.efficacy,
    ".csv",sep=""))



