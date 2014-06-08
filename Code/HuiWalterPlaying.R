library(INLA)


#generate some data for a hui walter paradigm

Se1<-0.4
Se2<-0.7
Sp1<-0.95
Sp2<-0.90
PrevalenceA<-10/100
PrevalenceB<-30/100
SizeA<-10000
SizeB<-5000
Y<-data.frame(Pop=c(rep("A",SizeA),rep("B",SizeB)))
Y$Case[Y$Pop=="A"]<-rbinom(SizeA,1,p=PrevalenceA)
Y$Case[Y$Pop=="B"]<-rbinom(SizeB,1,p=PrevalenceB)


Y$t1[Y$Case==1]<-rbinom(sum(Y$Case==1),1,Se1)
Y$t1[Y$Case==0]<-rbinom(sum(Y$Case==0),1,1-Sp1)
Y$t2[Y$Case==1]<-rbinom(sum(Y$Case==1),1,Se2)
Y$t2[Y$Case==0]<-rbinom(sum(Y$Case==0),1,1-Sp2)

y.agg<-rbind(data.frame(Pop=c("A","B"),Pop.size=c(SizeA,SizeB),test="t1",result=table(Y$Pop,Y$t1)[,2],
                        Cases=table(Y$Pop,Y$Case)[,2]
),
             data.frame(Pop=c("A","B"),Pop.size=c(SizeA,SizeB),test="t2",result=table(Y$Pop,Y$t2)[,2],
                        Cases=table(Y$Pop,Y$Case)[,2]))
y.agg


library(rjags)

temp<-jags.model("HuiWalter.txt",data=list(T1A.pos=y.agg[1,4],
                                     T2A.pos=y.agg[3,4],
                                     T1B.pos=y.agg[2,4],
                                     T2B.pos=y.agg[4,4],
                                     T1A.neg=SizeA-y.agg[1,4],
                                     T2A.neg=SizeA-y.agg[3,4],
                                     T1B.neg=SizeB-y.agg[2,4],
                                     T2B.neg=SizeB-y.agg[4,4],
                                     N.A=SizeA,
                                     N.B=SizeB),n.chains=5)

update(temp,10000)
tmp.samples<-coda.samples(temp,c("Se1","Se2","Sp1","Sp2","PrA","PrB"),5000)
summary(tmp.samples)
gelman.diag(tmp.samples)
plot(tmp.samples)
#something wrong in the model assumption here...


temp2<-jags.model("HuiWalterMultiNomial.txt",data=list(
  OA=c(with(subset(Y,Pop=="A"),table(t1,t2))),
                                            OB=c(with(subset(Y,Pop=="B"),table(t1,t2))),
                                            N.A=SizeA,
                                            N.B=SizeB),n.chains=5)
  update(temp2,10000)
  tmp2.samples<-coda.samples(temp2,c("Se1","Se2","Sp1","Sp2","PrA","PrB"),2000)



formula = result~Pop+Pop:test
result=inla(formula,data=y.agg,
            family="binomial",
            Ntrials=Pop.size,
            control.predictor=list(compute=T),
            control.compute=list(dic=T))




formula