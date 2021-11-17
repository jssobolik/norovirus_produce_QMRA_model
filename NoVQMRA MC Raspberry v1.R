#R code for QMRA model of human norovirus contamination of produce by harvester and packer hands in an agriculture setting;
# For additional information, refer to the published paper at doi: 10.1016/j.ijfoodmicro.2021.109365;
#license: GNU General Public License v3.0;

###################################
#Note: for the raspberry analysis there is assumed to be no packer contamination of the produce


rm(list=ls()) #clear all variables (this is good to have at the top of the script to clear any preexisting variables in your environment that might interfere)

#load packages
library(mc2d) 
library(mvtnorm)
#run 10000 trials of simulation
ndvar(10001)
#adding in uncertainty parameters
ndunc(101)
set.seed(12345)

#assign variable distributions
#Harvester variables
#Concentration of NoV in stool
conc<-mcstoc(rpert,type="VU",min=100,mode=1000000,max=100000000000,shape=10)
#mass of feces per both hands log (g feces/both hands) Jacxsens et al 2017
fh <- mcstoc(rbetagen, type="VU", shape1=4.57, shape2=2.55, min=-8.00, max=-1.00)
#total surface area of one hand that touches produce
harea<-245
#hand wash efficacy log removal 0, 0.5, 6
hweff<-mcstoc(rtriang,type="V", min=6,mode=6,max=6)
#Probability of harvester hand wash 0.27, 0.73
Phhw<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(1,0))
#Probability of harvester using gloves (Jaykus et al 2009)-updated following discussion with Kira to empirical distribution
#0.25 0.75
Phwg<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(0,1))
#Norovirus prevalence among people without gastroenteritis; (0.04,0.96) 
prevh<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(1,0))

#Packer variables
#concentration of NoV in stool
concp<-mcstoc(rpert,type="VU",min=0,mode=0,max=0,shape=0)
#mass of feces on both hands log (g feces/both hands) Jacxsens et al 2017
fhp <- mcstoc(rbetagen, type="VU", shape1=4.57, shape2=2.55, min=-8.00, max=-1.00)
#total surface area of one hand that touches produce
hareap<-245
#Hand wash efficacy log removal
hweffp<-mcstoc(rtriang,type="V", min=0,mode=0,max=0)
#Probability of packer hand wash
Pphw<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(0,1))
#Probability of packer using gloves (Jaykus et al 2009)-updated following discussion with Kira to empirical distribution
Ppwg<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(0,1))
#Norovirus prevalence among people without gastroenteritis 
prevp<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(1,0))

#transfer rates

#Harvest hands to produce
Trhhp<-mcstoc(rlnorm,type="VU",meanlog=-2.32,sdlog=0.15)
#Harvest gloved-hands to produce (Jaykus 2009)
Trhgp<-mcstoc(rbeta,type="VU",shape1=1.165,shape2=2.630)
#Packing produce to hands
Trpph<-mcstoc(rbeta,type="VU",shape1=7.42,shape2=39.51)
#Packing hands to produce
Trphp<-mcstoc(rlnorm,type="VU",meanlog=-2.32,sdlog=0.15)
#Packing gloved-hands to produce (Jaykus 2009)
Trpgp<-mcstoc(rbeta,type="VU",shape1=1.165,shape2=2.630)

#norovirus virion average daily reduction log-10 units per day (5C 0.011 log-10 units) versus
#0.151 for 20C; assumed produce are eaten within 7 days of harvesting
norodecaylow<-0.011

P <-mcstoc(rtriang,type="V", min=0.63,mode=0.722,max=0.8)
mu <-mcstoc(rtriang,type="V", min=399,mode=1106,max=2428)

#Harvest produce to bin
Trhpbi<-mcstoc(rtriang,type="VU",min=0.009,mode=0.028,max=0.038)
#Packing produce to belt
Trppbe<-mcstoc(rtriang,type="VU",min=0.009,mode=0.028,max=0.038)
#Harvester NoV transfer from hand/restroom envir. to glove
Trhhg<-mcstoc(runif,type="VU",min=0,max=0.444)
#Packer NoV transfer from hand/restroom envir. to glove
Trphg<-mcstoc(runif,type="VU",min=0,max=0.444)
#Produce to produce NoV transfer rates in bin (Jaykus 2009) Note these are produce to hand
Trhppbi<-mcstoc(rbeta, type="VU",shape1=2.257,shape2=15.502)

#Number of contacts
#Harvest hand produce (produce touch from harvest to bin)
nhPh<-mcstoc(rempiricalD,type="V",values=c(1),prob=c(1))
#Packing hand produce (from belt to final packed box)
npPh<-mcstoc(rempiricalD,type="V",values=c(1),prob=c(1))
#Produce bin
nPbi<-mcstoc(rempiricalD,type="VU",values=c(1,2,3,4,5),prob=c(0.2,0.2,0.2,0.2,0.2))
#Produce belt
nPbe<-mcstoc(rempiricalD,type="VU",values=c(1,2,3,4,5),prob=c(0.2,0.2,0.2,0.2,0.2))
#Produce to produce contacts in bin
nPPbi<-mcstoc(rempiricalD,type="VU",values=c(1,2,3,4,5,6,7,8,9,10),prob=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))

#calculated variables based on distributions 
#NoV shedding (virus/gram feces) for the harvester - with global prevalence factor
shedH<-conc*prevh
#g feces/harvester's one hand
mfh<-10^(fh/2)
#NoV on harvester's hand - autoinnoculation
nvhA<-(mfh*shedH)
#NoV on harvester's hand following conditional of hand washing (assuming restroom environment is not a source of NoV)
nvhH<-nvhA/(10^(hweff*Phhw))
#NoV on gloved harvester hands 
nvhG<-(nvhH - nvhH*((1-Trhhg)*Phwg))
#NoV shedding for the packer - with global prevalence factor
shedP<-concp*prevp
#g feces on packer's hand 
mfhp<-10^(fhp/2)
#NoV on packer's hand - autoinnoculation
nvpA<-(mfhp*shedP)
#removed environmental contamination; virus on hand leaving restroom after conditional hand washing
nvpH<-nvpA/(10^(hweffp*Pphw))
#NoV on gloved packer hands 
nvpG<-(nvpH - nvpH*((1-Trphg)*Ppwg))
#NoV on produce after harvest (with transfer switch based on gloves y/n) updated per conversation with Kira
#accounting for small raspberry size here
ifelse(Phwg>0, nvhP1<-(2.1/245)*(nvhG- nvhG*((1-Trhgp)^nhPh)), nvhP1<-(2.1/245)*(nvhG- nvhG*((1-Trhhp)^nhPh)))
#NoV on bin after contact with produce (assumed to be clamshell)
nvBi<- nvhP1-nvhP1*((1-Trhpbi)^nPbi)
#NoV on produce after contacting bin (assumed to be clamshell)
nvhP2<-nvhP1-nvBi
#Produce-produce dilution: NoV on produce in clamshell after produce contact
nvPPBi<- nvhP2-nvhP2*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in clamshell
nvhP3<-nvhP2-nvPPBi
#NoV on belt after contacting produce
#nvBe<- nvhP3-nvhP3*(1-Trppbe)^nPbe
#NoV on produce after contacting belt
#nvpP1<-nvhP3-nvBe
#NoV on produce after packing by hand, nvpP2 is variable for NoV on produce after packing; update ifelse for glove switch per conversation with Kira
#ifelse(Ppwg>0, 
  #     ifelse(nvpP1>nvpG, nvpP2<-nvpP1- (nvpP1-nvpP1*((1-Trpgp)^npPh)), nvpP2<-nvpP1+( (2.1/245)*(nvpG-nvpG*((1-Trpgp)^npPh)))),
 #      ifelse(nvpP1>nvpG, nvpP2<-nvpP1- (nvpP1-nvpP1*((1-Trpph)^npPh)), nvpP2<-nvpP1+((2.1/245)*(nvpG-nvpG*((1-Trphp)^npPh))))
#)
#NoV on packer hand after packing produce (corrected by deleting extra -nvpH);update ifelse for glove switch per conversation with Kira
#ifelse(Ppwg>0, 
 #      ifelse(nvpP1>nvpG, nvpH1<-(nvpP1-nvpP1*((1-Trpgp)^npPh))+nvpG, nvpH1<- nvpG*((1-Trpgp)^npPh)), 
  #     ifelse(nvpP1>nvpG, nvpH1<-(nvpP1-nvpP1*((1-Trpph)^npPh))+nvpG, nvpH1<- nvpG*((1-Trphp)^npPh))
#)


#calculating total viral load per produce items factoring in viral decay at 5C/20C, assumed 7 days total of viral decay
nvpP2decayl <- nvhP3/(10^(norodecaylow*7))


#body weight pulled for adults from exposure handbook, kg)
bodyweight<-80

#USDA reported weights in grams per produce item
#https://hannaone.com/Recipe/weightcantaloupe.html
#https://hannaone.com/Recipe/weightlettuce.html
#https://hannaone.com/Recipe/weighttomato.html
weightlettucehead<-mcstoc(rtriang,type="V", min=309,mode=360,max=626)
#weightmelon<-mcstoc(rtriang,type="V", min=441,mode=552,max=814)
weighttomato<-mcstoc(rtriang,type="V", min=91,mode=123,max=182)
weightberry<-mcstoc(rtriang,type="V", min=3,mode=4,max=5)

#Exposure handbook intake of specific fruits and veggies (daily consumption g/kg-day)
lettuceconsumpt<-mcstoc(rnorm,type="V",mean=0.44, sd=0.01)
#melonconsumpt<-mcstoc(rnorm,type="V",mean=0.70, sd=0.05)
tomatoconsumpt<-mcstoc(rnorm,type="V",mean=0.83, sd=0.02)
berryconsumpt<-mcstoc(rnorm,type="V",mean=0.45, sd=0.02)

#dose calculation from RASPBERRY ending with total virions/day
dose1rasp <- (bodyweight*berryconsumpt*nvpP2decayl)/weightberry


#calculating risk from RASPBERRY with fractional Poisson model aggregated viral decay at 20C

riskfracpoissR1 <- P*(1-exp(-dose1rasp/mu))

###putting everything through the mc function
#sensitivity for risk raspberry
QMRAberry<-mc(conc,fh,hweff, Phhw,Phwg, prevh, concp, fhp, hweffp, Pphw, Ppwg, prevp, Trhhp, Trhgp,Trpph,Trphp,Trpgp, P, mu,
              Trhpbi,Trppbe, Trhhg,Trphg,Trhppbi,nhPh,npPh,nPbi,nPbe,nPPbi,shedH, mfh, nvhA, nvhH, nvhG, shedP, mfhp, nvpA, nvpH, nvpG, nvhP1,
              nvBi,nvhP2,nvPPBi,nvhP3,nvpP2decayl, weightlettucehead, weighttomato, weightberry, lettuceconsumpt, tomatoconsumpt, berryconsumpt, 
              dose1rasp,riskfracpoissR1)

berrydata<-summary(QMRAberry, digits=2)
dfberry<-as.data.frame(print(berrydata))
write.csv(dfberry, "C:\\Users\\jsoboli\\Desktop\\fig3Aberry6log100.csv", row.names=TRUE)


####################################################################################################################################################################################################################################################################


