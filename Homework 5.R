#####################################################################
#Homework: estimating MSY, Bmsy, umsy from a 
#general stock assessment model
#FISH 458 School of Aquatic and Fishery Sciences
#University of Washington
#Michael Paschal, paschalm@uw.edu
#May 1, 2018 - May 7, 2018
#####################################################################

rm(list=ls())     #get rid of existing global variables!!!

#list of global model parameter values
alpha <- 0.0001   #length-weight alpha
beta <- 3.05      #length-weight beta
Linfinity <- 80   #von Bertalanffy Linfinity
k <- 0.2          #von Bertalanffy K
t0 <- -0.2        #von Bertalanffy t0
h <- 0.7          #Beverton-Holt steepness is R/R0 when S = 0.2SSB0. 
R0 <- 1000000     #number of recruits in unfished population
Lv50 <- 40        #length at which 50% of fish are vulnerable to gear
Lv95 <- 45        #length at which 95% of fish are vulnerable to gear
pfemale <- 0.5    #proportion of the population that is female
Lm50 <- 30        #length at which 50% of females are mature
Lm95 <- 35        #length at which 95% of females are mature
surv <- 0.8       #natural survival, assumed constant for all ages
n <- 10           #plus group age 

#=========PART ONE================================================================
#Calculates unfished equilibrium numbers at age, including the plus group.
#Inputs: unfished rec, natural survival, and
#plus group age. Assumes R0 is for age 1 individuals. 
#=================================================================================
unfished.num <- function(R0, surv, n) {
   N.vec <- vector(length=n)
   N.vec[1] = R0
   for(yr in 2:(n-1)) {
     N.vec[yr] = N.vec[yr-1] * surv
   }
   N.vec[n] = (surv/(1-surv))*N.vec[n-1]
   return(N.vec)
}
x <- unfished.num(R0=R0, surv=surv, n=n)
print(x)

#=========PART TWO================================================================
#given a vector of ages, and parameters, returns a vector of weights at those 
#ages. 
#Inputs: age.vec a vector of ages
# alpha, beta of length-weight; Linfinity, k, t0 of Von Bertalanffy
#=================================================================================
weight.at.age <- function(age.vec, alpha, beta, Linfinity, k, t0) {
   N.ages <- length(age.vec)
   length.vec = vector(length = N.ages)
   for(t in 1:N.ages) {
     length.vec[t] = Linfinity*(1-exp(-k*(t-t0)))
   }
   weight.vec = alpha*(length.vec^beta)
   
   return(weight.vec)
}
x <- weight.at.age(age.vec=1:n,alpha=alpha, beta=beta, Linfinity=Linfinity, k=k,t0=t0)
print(x)

#=========PART THREE==============================================================
#Proportion of population at age that is spawning females
#Given a vector of ages, returns the proportion at each that is female and spawning.
#Assumes a logistic curve with Lm50 the length at which 50% are mature and 
#Lm95 the length at which 95% are mature
#Inputs: age.vec a vector of ages, Linfinity, k, t0, Lm50, Lm95, pfemale
#=================================================================================
prop.fem.spawners <- function(age.vec, Linfinity, k, t0, Lm50, Lm95, pfemale) {
   N.ages <- length(age.vec)
   length.vec = vector(length = N.ages)
   for(t in 1:N.ages) {
     length.vec[t] = Linfinity*(1-exp(-k*(t-t0)))
   }
   prop.mature = pfemale / (1+exp(-log(19)*(length.vec-Lm50)/(Lm95-Lm50)))
   
   return(prop.mature)
}
x <- prop.fem.spawners(age.vec=1:n,Linfinity=Linfinity, k=k,t0=t0, Lm50=Lm50, 
                       Lm95=Lm95, pfemale=pfemale)
print(x)

#=========PART FOUR===============================================================
#Calculate SSB0, which comes from the sum of the spawning biomass per recruit
#multiplied by the number of recruits. In this model this is also the unfished 
#numbers-at-age multiplied by the weight-at-age 
#multiplied by the proportion of the fish that are spawning females.
#Inputs: n, alpha, beta, Linfinity, k, t0, Lm50, Lm95, pfemale, surv, R0
#=================================================================================
calc.SSB0 <- function(n, alpha, beta, Linfinity, k, t0, Lm50, Lm95, pfemale, surv, R0) {
   N.vec = unfished.num(R0=1, surv = surv, n = n)
   prop.mature = prop.fem.spawners(age.vec=1:n, Linfinity=Linfinity, k=k, t0=t0,
                                      Lm50=Lm50, Lm95=Lm95, pfemale=pfemale)
   weight.vec = weight.at.age(age.vec=1:n, alpha=alpha, beta=beta, k=k, t0=t0, Linfinity=Linfinity)
   SBPR0 = 0
   for(i in 1:length(N.vec)) {
     SBPR0 = SBPR0 + (N.vec[i]*weight.vec[i]*prop.mature[i])
   }
   SSB0 = SBPR0 * R0
   return(SSB0)
}
x <- calc.SSB0(n=n, alpha=alpha, beta=beta, Linfinity=Linfinity, k=k, t0=t0, 
          Lm50=Lm50, Lm95=Lm95, pfemale=pfemale, surv=surv, R0=R0)
print(x)

#==========PART FIVE===============================================================
#Vulnerability to fishing, as a function of length, converts an age-vector into 
#vulnerability values by age. First converts age to length, then converts length
#to vulnerability assuming a logistic selectivity curve 
#=================================================================================
prop.vulnerable <- function(age.vec, Linfinity, k, t0, Lv50, Lv95) {
   N.ages <- length(age.vec)
   length.vec = vector(length = N.ages)
   for(t in 1:N.ages) {
     length.vec[t] = Linfinity*(1-exp(-k*(t-t0)))
   }
   
   prop.vuln = 1 / (1 + exp(-log(19) * (length.vec-Lv50)/(Lv95-Lv50)))
   return(prop.vuln)
}
x <- prop.vulnerable(age.vec=1:n,Linfinity=Linfinity, k=k,t0=t0, Lv50=Lv50, Lv95=Lv95)
print(x)

#==========PART SIX===============================================================
#Numbers of fish in the population as function of harvest rate u, recruits R, plus 
#age a, the vulnerability parameters, and von Bertalanffy parameters.
#Returns a vector of numbers at age
#===================================================================================
calc.pop.size <- function(R, n, u, surv, Linfinity, k, t0, Lv50, Lv95) {
   N.vec = vector(length = n)
   N.vec[1] = R
   v = prop.vulnerable(age.vec=1:n, Linfinity = Linfinity, k=k, t0=t0,
                       Lv50=Lv50, Lv95=Lv95)
   for(i in 2:(n-1)) {
     N.vec[i] = N.vec[i-1]*surv*(1-v[i-1]*u)
   }
   N.vec[n] = (surv*(1-v[n]*u))/(1-surv*(1-v[n]*u))*N.vec[n-1]
   return(N.vec)
}
x <- calc.pop.size(R=1, n=n, u=0.1, surv=surv, Linfinity=Linfinity, k=k, 
              t0=t0, Lv50=Lv50, Lv95=Lv95)
print(x)

#==========PART SEVEN==============================================================
#Calculate equilibrium yield and spawning biomass. 
#Requires all the parameters of the model
#Calls all the little bits and pieces calculated so far. 
#==================================================================================
eqm.Y.SSB <- function(alpha, beta, Linfinity, k, t0, h, R0, Lv50, Lv95,
                     pfemale, Lm50, Lm95, surv, n, u) {
   
   #numbers at age from 1 recruit as a function of R0
   N.vec = calc.pop.size(R=R0, n=n, u=u, surv=surv, Linfinity=Linfinity, k=k, 
                         t0=t0, Lv50=Lv50, Lv95=Lv95)
   
   #proportion of spawners, and weight at age
   prop.spawners = prop.fem.spawners(age.vec=1:n,Linfinity=Linfinity, k=k,t0=t0, Lm50=Lm50, 
                                     Lm95=Lm95, pfemale=pfemale)
   
   #SBPR(u) spawning biomass per recruit as a function of u
   weight.vec = weight.at.age(age.vec=1:n,alpha=alpha, beta=beta, Linfinity=Linfinity, k=k,t0=t0)
   mat.fem = prop.fem.spawners(age.vec=1:n,Linfinity=Linfinity, k=k,t0=t0, Lm50=Lm50, 
                               Lm95=Lm95, pfemale=pfemale)
   SBPR = 0
   for(i in 1:n) {
     SBPR = SBPR + N.vec[i]*weight.vec[i]*mat.fem[i]
   }
   
   #YPR(u) yield per recruit at harvest rate u
   v = prop.vulnerable(age.vec=1:n, Linfinity = Linfinity, k=k, t0=t0,
                       Lv50=Lv50, Lv95=Lv95)
   YPR = 0
   for(i in 1:n) {
     YPR = YPR + N.vec[i]*weight.vec[i]*v[i]*u
   }
   
   #caculate SSB0 (for Beverton-Holt)
   SSB0 = calc.SSB0(n=n, alpha=alpha, beta=beta, Linfinity=Linfinity, k=k, t0=t0, 
                    Lm50=Lm50, Lm95=Lm95, pfemale=pfemale, surv=surv, R0=R0)
   
   a = ((1-h) / (4*h*R0)) * SSB0
   b = ((5*h)-1) / (4*h*R0)
   
   #calculate equilibrium number of recruits from Beverton-Holt and SBPR(u)
   eqm.recruits = (SBPR - a) / (b*SBPR)
   
   #calculate total equilibrium yield and spawning biomass at this exploitation rate
   yield = YPR * eqm.recruits
   SSB = SBPR * eqm.recruits
   
   #returns a list, to access elements use x[[1]] and x[[2]]
   return(list(yield=yield, SSB=SSB, YPR=YPR, SBPR=SBPR, eqm.recruits=eqm.recruits))
   
}
#calculate equilibrium yield and spawning biomass for a given harvest rate
x <- eqm.Y.SSB(alpha=alpha, beta=beta, Linfinity=Linfinity, k=k, t0=t0,
                    h=h, R0=R0, Lv50=Lv50, Lv95=Lv95,
                    pfemale=pfemale, Lm50=Lm50, Lm95=Lm95, surv=surv, n=n, u=0.1)
print(x)

#==========PART EIGHT==============================================================
#Calculate yield, SSB for a range of exploitation rates from 0 to 1
#Requires all the parameters of the model
#Calls all the little bits and pieces calculated so far. 
#==================================================================================
calc.ref.points <- function(alpha, beta, Linfinity, k, t0, h, R0, Lv50, Lv95,
                             pfemale, Lm50, Lm95, surv, n) {
   harvest.rates = seq(from = 0, to = 1, by = 0.001)
   ncalcs = length(harvest.rates)
   results = matrix(nrow=ncalcs, ncol=3)
   results[,1] = harvest.rates
   for(i in 1:ncalcs) {
     temp = eqm.Y.SSB(alpha=alpha, beta=beta, Linfinity=Linfinity, k=k, t0=t0,
                      h=h, R0=R0, Lv50=Lv50, Lv95=Lv95,
                      pfemale=pfemale, Lm50=Lm50, Lm95=Lm95, surv=surv, n=n, u=harvest.rates[i])
     results[i,2] = temp$yield
     results[i,3] = temp$SSB
   }
   return(results)
}
x <- calc.ref.points(alpha=alpha, beta=beta, Linfinity=Linfinity, k=k, t0=t0,
                           h=h, R0=R0, Lv50=Lv50, Lv95=Lv95,
                           pfemale=pfemale, Lm50=Lm50, Lm95=Lm95, surv=surv, n=n)
head(x)   #print out the first few rows

#write the results to a file for later inspection to find the MSY etc. 
write.csv(x=x, file="MSY results.csv")   

#create plot for the writeup
plot(x[,1], x[,2], xlab = "Harvest Rate", ylab = "Yield", main = "Yield given various harvest rates", type="l")
which.max(x[,2])
x[390,]
#maximum yield occurs at u=0.389, yield = 7.002e12, SSB=1.127e13



