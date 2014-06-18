library(rjags)
library(batch)
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



Prevalence<-c(High=6/100,Low=1/100)
sim.scenarios<-list(
  baseline=list(
    SeStd=c(Vacc=0.7,nonVacc=0.7),
    SeViva=c(Vacc=0.7,nonVacc=0.7),
    SpStd=c(Vacc=0.999,nonVacc=0.999),  ## a specificity of 0.5 mirrors that it reacts to vaccinated animals
    SpViva=c(Vacc=0.999,nonVacc=0.999),
    vaccine.efficacy=0.6,
    Prevalence=c(High=6/100,Low=1/100),
    props=c(1/4,1/4,1/4,1/4)# balance between populations High/vacc,Low/vacc,High/nonvacc,low/nonvacc
    
  ))

sample.size.range<-seq(5000,100000,by=5000)
current.seed<-as.integer(as.numeric(format(Sys.time(),"%Y%m%d%H")))

for(rep in 0:1){
  for(simvars in sim.scenarios){
    for(Sample.size in sample.size.range){
               rbatch(rfile=shQuote(file.path(script.dir,"HuiWalterDataGen.R")),
               SeStd=simvars$SeStd,
               SeViva=simvars$SeViva,
               SpStd=simvars$SpStd,
               SpViva=simvars$SpViva,
               vaccine.efficacy=simvars$vaccine.efficacy,
               Prevalence=simvars$Prevalence,
               props=simvars$props,
               samplesize=Sample.size,
               ##Controlling parameters here
               nreps=9,
               rep.prefix=rep,
               nchains=5,
               niter=1000,
               n.mcmc.samples=1000,
               save.samples=F,
               seed=current.seed
        )
      }
    }
  }

rbatch.local.run(2)
