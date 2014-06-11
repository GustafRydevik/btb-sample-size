library(INLA)


#generate some data for a hui walter paradigm

SeStd<-c(Vacc=0.7,nonVacc=0.7)
SeViva<-c(Vacc=0.7,nonVacc=0.7)
SpStd<-c(Vacc=0.55,nonVacc=0.997)  ## a specificity of 0.5 mirrors that it reacts to vaccinated animals
SpViva<-c(Vacc=0.999,nonVacc=0.999)

Prevalence<-c(High=5/100,Low=2/100)
nCattle<-matrix(c(nHigh.vacc=7500,nLow.vacc=7500,nHigh.nonvacc=7500,nLow.nonvacc=7500),ncol=2)
propVaccinated<-1/2
btb.testdata<-data.frame(Pop=c(rep("High",sum(nCattle[1:2])),rep("Low",sum(nCattle[1:2]))))
btb.testdata$Case[btb.testdata$Pop=="High"]<-rbinom(nCattle[1],1,p=Prevalence[1])
btb.testdata$Case[btb.testdata$Pop=="Low"]<-rbinom(nCattle[2],1,p=Prevalence[2])
btb.testdata$Vaccinated<-c(rep("Vacc",nCattle[1]),
                rep("nonVacc",nCattle[2]),    
                rep("Vacc",nCattle[3]),
                rep("nonVacc",nCattle[4])
)
btb.testdata$Vaccinated<-factor(btb.testdata$Vaccinated,levels=c("Vacc","nonVacc"))            

  for(vacc.status in c("Vacc","nonVacc")){
    SubsetPos<-with(btb.testdata,Case==1&Vaccinated==vacc.status)
    SubsetNeg<-with(btb.testdata,Case==0&Vaccinated==vacc.status)

    btb.testdata$Std[SubsetPos]<-rbinom(sum(SubsetPos),1,SeStd[vacc.status])
    btb.testdata$Std[SubsetNeg]<-rbinom(sum(SubsetNeg),1,1-SpStd[vacc.status])
    btb.testdata$Viva[SubsetPos]<-rbinom(sum(SubsetPos),1,SeViva[vacc.status])
    btb.testdata$Viva[SubsetNeg]<-rbinom(sum(SubsetNeg),1,1-SpViva[vacc.status])
    
   }



btb.testdata.agg<-as.array(table(btb.testdata[,c("Pop","Vaccinated","Std","Viva")]))
btb.dimnames<-c(dimnames(btb.testdata.agg)[1:2],Std.vs.Viva=list(c("--","+-","-+","++")))
dim(btb.testdata.agg)<-c(2,2,4)
dimnames(btb.testdata.agg)<-btb.dimnames
class(btb.testdata.agg)<-"array"
library(rjags)

temp<-jags.model("./Code/HuiWalterMultinomial_vaccinePop.txt",data=list(counts=btb.testdata.agg,
                                     nCattle=nCattle),
                                     n.chains=5)

update(temp,10000)
tmp.samples<-coda.samples(temp,c("SeStd","SeViva","SpStd","SpViva","Prevalence","probResult"),5000)
tmp2<-do.call("rbind",tmp.samples)
class(tmp2)<-"mcmc"
posterior.intervals<-HPDinterval(tmp2)[1:10,]
rownames(posterior.intervals)<-gsub("\\[1\\]"," vaccine pop",rownames(posterior.intervals))
rownames(posterior.intervals)<-gsub("\\[2\\]"," nonvaccine pop",rownames(posterior.intervals))
rownames(posterior.intervals)[1:2]<-c("High prevalence pop","Low prevalence pop")
posterior.intervals<-data.frame(posterior.intervals)
posterior.intervals[,"True.value"]<-c(
  Prevalence<-c(High=5/100,Low=2/100),
  SeStd<-c(Vacc=0.7,nonVacc=0.7),
SeViva<-c(Vacc=0.7,nonVacc=0.7),
SpStd<-c(Vacc=0.55,nonVacc=0.997),  ## a specificity of 0.5 mirrors that it reacts to vaccinated animals
SpViva<-c(Vacc=0.999,nonVacc=0.999))


print(xtable(posterior.intervals,digits=4),type="html",file= "PosteriorIntervals.html")
