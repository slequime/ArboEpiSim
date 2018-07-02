#Variables to change####
EIP50 <- 7
slopeInfec <-2.5
Pmax <- 0.8
df = expand.grid(EIP50=EIP50,slopeInfec=slopeInfec,Pmax=Pmax)

for(condition in 1:nrow(df)){
  
ongoing.condition=df[condition,]

#Run in parallel 

library(foreach)
library(doParallel)

cores=detectCores()
cores=2
cl = makeCluster(cores[1]-1)
registerDoParallel(cl)

#Number of replicates

replicates=100

Outbreak_indices = foreach(j=1:replicates, .combine = rbind) %dopar% {

library(plyr)
library(dplyr)
library(reshape2)
library(data.table)
library(grid)
library(stringr)

####New Functions####

#viral.transmission.mos = function(x, Pmax, EIP, slope) {Pmax*1/2*(1+tanh((x-EIP)/slope))}
viral.transmission.mos = function(x, Pmax, EIP, slope) {Pmax / (1 + exp(-slope*(x-EIP))) }

viremia = function(x, Pmax, IIP, slopeIIP, Recov, slopeRecov) {
  (Pmax / (1 + exp(-slopeIIP*(x-IIP)))) - (Pmax/(1 + exp(-slopeRecov*(x-Recov))))
}

viremia2 = function(x, Pmax, IIP, slopeIIP, Recov, slopeRecov, dist) {
  ((Pmax / (1 + exp(-slopeIIP*(x-IIP)))) - (Pmax/(1 + exp(-slopeRecov*(x-Recov)))))-dist
}

rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

#####Initial parameters####

#Variable parameters Mosquito dynamics
EIPmean= ongoing.condition[,4]#M parameter

slopeMos=ongoing.condition[,3]#B parameter

Pmax.mos.putative=ongoing.condition[,2]#K parameter
    
#Same for all sim
Nvert=10000 #population (number of individuals) of vertebrates
Dvector=3 #relative density of vectors

#Dynamics of viremia in vertebrate (see viremia validation for justification)
IIPmean=4.844
IIPsd=0.5
slopeIIP=2.536
slopeIIP.sd=0.2
Pmax.vert.putative=1
Pmax.vert.var=0.002
Pmax.vert=Pmax.vert.putative-Pmax.vert.var
  #Beta law parameters
    alpha.vert=(Pmax.vert*(1-Pmax.vert)/Pmax.vert.var)*Pmax.vert
    beta.vert=(alpha.vert/Pmax.vert)-alpha.vert
Recov=8.392
Recov.sd=0.5
slopeRecov=4.038
slopeRecov.sd=0.2

# IIPmean=7 #mean length (in days) of the intrasinc incubation period
# IIPsd=1 #standard deviation (in days) of the intrinsinc incubation period
# 
# mean.length.dis=4 #mean length (in days) of the total infectious period (after IIP; this date is Recov50)
# sd.length.dis=2 #standard deviation (in days) of the total infectious period (after IIP)

surv.rate.vector=0.85 #probability of a vector to die each day
biting.rate.vector=0.3 #daily biting rate of a vector

# inf.biting.vector=0.8 #proportion of bite that are infectious
# inf.biting.vertebrate=0.85 #proportion of bite that gets the vector infected if host is infected

sim.length=400 #simulation length in days
N.vectors.infected.ini=0 #number of mosquitoes initially infected
N.hosts.infected.ini=1 #number of hosts initially infected

###Tables 

ID.vert=c(1:Nvert)
ID.vect=c(1:(Nvert*Dvector))
inf.by = as.numeric(rep(NA,Nvert))

attempt <-10
for(i in 1:attempt){
table1.vert = data.table(ID.vert,inf.by)
table1.vert$infected = as.numeric(rep(NA,Nvert))
table1.vert$IIP50 = rnorm(Nvert,IIPmean,IIPsd)
table1.vert$slopeInfec = rnorm(Nvert,slopeIIP,slopeIIP.sd)
table1.vert$Pmax = rbeta(Nvert,alpha.vert,beta.vert)
table1.vert$slopeRecov = rnorm(Nvert,slopeRecov,slopeRecov.sd)
table1.vert$Infec.pop = as.numeric(rep(NA,Nvert))

table1.vert = ddply(table1.vert, .(ID.vert), mutate, Recov50=rtnorm(1,Recov,Recov.sd, a=IIP50+1, b=Inf))

table1.vert = tryCatch(ddply(table1.vert, .(ID.vert), mutate, End.date=uniroot(viremia2, Pmax=Pmax, IIP=IIP50, slopeIIP=slopeInfec, Recov=Recov50, slopeRecov=slopeRecov, dist=0.01, upper=Recov50*2, lower=Recov50, extendInt="yes", maxiter=100000)$root), error=function(e) "skip")
if(is.data.frame(table1.vert)==TRUE){break}
}

table1.vert=as.data.table(table1.vert)
setkey(table1.vert, ID.vert)

#Vector table

Total.mos.in.sim=(Nvert*Dvector)+(Nvert*Dvector)*(1-surv.rate.vector)*sim.length

final=NULL
for (i in 1:sim.length){
  interm = rep(i, (Nvert*Dvector)*(1-surv.rate.vector))
  final = c(final, interm)}

ID.vect=c(1:Total.mos.in.sim)
birth.date = c(rep(0,(Nvert*Dvector)), final)

table1.vect = data.frame(ID.vect,birth.date)
table1.vect$dead.date = as.numeric(rep(NA,(Total.mos.in.sim)))
table1.vect$inf.by = as.numeric(rep(NA,(Total.mos.in.sim)))
table1.vect$infected = as.numeric(rep(NA,(Total.mos.in.sim)))
table1.vect$Infec.pop = as.numeric(rep(NA,(Total.mos.in.sim)))

result = c()
for(i in 1:Total.mos.in.sim){
  u = rbinom(1,1, Pmax.mos.putative)
  if(u == 1) {
    result[i] = EIPmean
  }
  else if(u == 0) {
    result[i] = Inf
  }
}

result[result < 0] = 0

table1.vect$EIP50 = result

table1.vect = as.data.table(table1.vect)

setkey(table1.vect, ID.vect)

###Initial infections###
infected.vert.ID = NULL
if(N.hosts.infected.ini !=0){
table1.vert[1:N.hosts.infected.ini,]$infected = 0
table1.vert[1:N.hosts.infected.ini,]$inf.by = 0
infected.vert.ID = c(1:N.hosts.infected.ini)
}

infected.vect.ID = NULL
if(N.vectors.infected.ini !=0){
table1.vect[1:N.vectors.infected.ini,]$infected = 0
table1.vect[1:N.vectors.infected.ini,]$inf.by = 0
infected.vect.ID = c(1:N.vectors.infected.ini)
}

ID.generator=(Nvert*Dvector)
ID.vect = 1:(Nvert*Dvector)

####First layer, which mosquitoes bite####

for (t in 1:sim.length) {
#First step, surviving mosquitoes that day
dead.ID = sample(ID.vect, (Nvert*Dvector)*(1-surv.rate.vector), replace=FALSE)
surviving.ID = ID.vect[!(ID.vect %in% dead.ID)]

new.mos.ID = as.numeric((ID.generator+1):((ID.generator+length(dead.ID))))
ID.generator=ID.generator+length(dead.ID)

active.mosquitoes.ID=c(surviving.ID,new.mos.ID)
table1.vect[dead.ID,dead.date:=as.numeric(t)]
ID.vect = active.mosquitoes.ID

active.infected.vert = table1.vert[t < (infected+End.date),ID.vert]
recovered.infected.vert = table1.vert[t > (infected+End.date),ID.vert]

#Second step, biting mosquitoes that day

stinging.ID =  sample(active.mosquitoes.ID, (Nvert*Dvector)*biting.rate.vector, replace=FALSE)
stinged.ID = sample(ID.vert, length(stinging.ID), replace=TRUE)
sting.matrix = data.table(stinging.ID,stinged.ID)

stinging.infected = stinging.ID[stinging.ID %in% infected.vect.ID] #mosquitoes biting with infection
stinged.infected = unique(stinged.ID[stinged.ID %in% active.infected.vert]) #vertebrates biten with active infection
stinged.recovered = unique(stinged.ID[stinged.ID %in% recovered.infected.vert]) #vertebrates biten with active infection

events=sting.matrix
events$Infected.vert = ifelse(events$stinged.ID %in% stinged.infected, TRUE, FALSE)
events$Infected.vect = ifelse(events$stinging.ID %in% stinging.infected, TRUE, FALSE)
events$Recovered.vert = ifelse(events$stinged.ID %in% stinged.recovered, TRUE, FALSE)

events = subset(events, Recovered.vert==FALSE & (Infected.vert==TRUE | Infected.vect==TRUE))

####Second layer, what gives the bite####

if(nrow(events) != 0){

for (st in 1:nrow(events)){
  
  Event.ongoing = events[st]

  ##If mosquito is infected & Host not infected
  
  Minfected = Event.ongoing[,Infected.vect]
  Vinfected = Event.ongoing[,Infected.vert]

  if(Minfected == TRUE && Vinfected == FALSE){
    
    #take data from mosquito.ID table
    
    b.ongoing = Event.ongoing[,stinging.ID]
    x = t-table1.vect[b.ongoing, infected]
    b.EIP50 = table1.vect[b.ongoing, EIP50]

    Ptransmit.t.vect = ifelse(x >= b.EIP50, 1, 0)
    transmitVect2Vert=sample(c(TRUE,FALSE),1,prob=c(Ptransmit.t.vect, 1-Ptransmit.t.vect))

    if(transmitVect2Vert==TRUE){
      table1.vert[Event.ongoing[,stinged.ID],inf.by:=as.numeric(b.ongoing)]
      table1.vert[Event.ongoing[,stinged.ID],infected:=as.numeric(t)]
      table1.vert[Event.ongoing[,stinged.ID],Infec.pop:=as.numeric(Ptransmit.t.vect)]
      infected.vert.ID = unique(c(infected.vert.ID, Event.ongoing[,stinged.ID]))
    }
    }

  #If mosquito is not infected & Host is infected
    
  if(Minfected == FALSE && Vinfected == TRUE){
    a.ongoing = Event.ongoing[,stinged.ID]
    x = t-table1.vert[a.ongoing, infected]
    a.Pmax = table1.vert[a.ongoing, Pmax]
    a.IIP50 =table1.vert[a.ongoing, IIP50]
    a.slopeInfec =table1.vert[a.ongoing, slopeInfec]
    a.Recov50 =table1.vert[a.ongoing, Recov50]
    a.slopeRecov =table1.vert[a.ongoing, slopeRecov]
    
Ptransmit.t.vert2 = viremia(x,
                            a.Pmax,
                            a.IIP50,
                            a.slopeInfec ,
                            a.Recov50,
                            a.slopeRecov)
    
Ptransmit.t.vert = ifelse(Ptransmit.t.vert2 > 0, Ptransmit.t.vert2, 0)
transmitVert2Vect=sample(c(TRUE,FALSE),1,prob=c(Ptransmit.t.vert, 1-Ptransmit.t.vert))
  
    if( transmitVert2Vect==TRUE){ 
      table1.vect[Event.ongoing[,stinging.ID],inf.by:=as.numeric(a.ongoing)]
      table1.vect[Event.ongoing[,stinging.ID],infected:=as.numeric(t)]
      table1.vect[Event.ongoing[,stinging.ID],Infec.pop:=as.numeric(Ptransmit.t.vert)]
      infected.vect.ID =unique(c(infected.vect.ID, Event.ongoing[,stinging.ID]))
    }
} 
}
}
}

test=Sys.time() - it.Start
#}

#})

data.summary = NULL
for(t in 1:sim.length){
  sum1=table1.vect[!is.na(table1.vect$infected),]
  sum2=table1.vert[!is.na(table1.vert$infected),]
  
  data.t=data.frame(t=t, 
                    vect.infected.cumulative=nrow(sum1[sum1$infected <= t,]),
                    vert.infected.cumulative = nrow(sum2[sum2$infected <= t,]))
  data.summary = rbind.fill(data.summary,data.frame(data.t))
}

#Modelisation of the cumulative outbreak dynamic:

outbreak.dynamic <- data.summary
outbreak.dynamic2 <- melt(outbreak.dynamic,id="t", variable_name="host")
outbreak.dynamic2$variable <- str_replace(outbreak.dynamic2$variable, "vect.infected.cumulative", "vector")
outbreak.dynamic2$variable <- str_replace(outbreak.dynamic2$variable, "vert.infected.cumulative", "vertebrate")

outbreak.dynamic2$Simulation_nb = j
outbreak.dynamic2$EIPmean=ongoing.condition[,1]
outbreak.dynamic2$slopeMos=ongoing.condition[,2]
outbreak.dynamic2$Pmax.mos=ongoing.condition[,3]
outbreak.dynamic2
}
stopCluster(cl)
write.table(Outbreak_indices, file=paste("EIP",ongoing.condition[,1],"_slope",ongoing.condition[,2],"_Pmax",ongoing.condition[,3],"_fixed_threshold_simulation.txt", sep=""), quote = F, sep = "\t", row.names = T)
}