#R code for QMRA model of human norovirus contamination of produce by harvester and packer hands in an agriculture setting;
# For additional information, refer to the published paper at doi: 10.1016/j.ijfoodmicro.2021.109365;
#license: GNU General Public License v3.0;

###################################
#load packages
library(mc2d) 
library(mvtnorm)
#run 10000 trials of simulation
ndvar(10001)
#adding in uncertainty parameters
ndunc(101)

#assign variable distributions
#Harvester variables
#Concentration of NoV in stool
conc<-mcstoc(rpert,type="VU",min=100,mode=1000000,max=100000000000,shape=10)
#mass of feces per both hands log (g feces/both hands) Jacxsens et al 2017
fh <- mcstoc(rbetagen, type="VU", shape1=4.57, shape2=2.55, min=-8.00, max=-1.00)
#total surface area of one hand that touches produce
harea<-245
#hand wash efficacy log removal
hweff<-mcstoc(rtriang,type="V", min=0,mode=0.5,max=6)
#Probability of harvester hand wash
Phhw<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(0.27,0.73))
#Probability of harvester using gloves (Jaykus et al 2009)-updated following discussion with Kira to empirical distribution
Phwg<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(0.25,0.75))
#Norovirus prevalence among people without gastroenteritis 
prevh<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(1,0))

#Packer variables
#concentration of NoV in stool
concp<-mcstoc(rpert,type="VU",min=100,mode=1000000,max=100000000000,shape=10)
#mass of feces on both hands log (g feces/both hands) Jacxsens et al 2017
fhp <- mcstoc(rbetagen, type="VU", shape1=4.57, shape2=2.55, min=-8.00, max=-1.00)
#total surface area of one hand that touches produce
hareap<-245
#Hand wash efficacy log removal
hweffp<-mcstoc(rtriang,type="V", min=0,mode=0.5,max=6)
#Probability of packer hand wash
Pphw<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(0.27,0.73))
#Probability of packer using gloves (Jaykus et al 2009)-updated following discussion with Kira to empirical distribution
Ppwg<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(0.25,0.75))
#Norovirus prevalence among people without gastroenteritis 
prevp<-mcstoc(rempiricalD,type="V",values=c(1,0),prob=c(1,0))

#transfer rates

#Harvest hands to produce - Bouwknet et al 2015 
Trhhp<-mcstoc(rlnorm,type="VU",meanlog=-2.22,sdlog=0.17)
#Harvest gloved-hands to produce (Jaykus 2009)
Trhgp<-mcstoc(rbeta,type="VU",shape1=1.165,shape2=2.630)
#Packing produce to hands
Trpph<-mcstoc(runif,type="VU",min=0.105,max=0.175)
#Packing hands to produce Bouwknet et al 2015
Trphp<-mcstoc(rlnorm,type="VU",meanlog=-2.22,sdlog=0.17)
#Packing gloved-hands to produce (Jaykus 2009)
Trpgp<-mcstoc(rbeta,type="VU",shape1=1.165,shape2=2.630)

#norovirus virion average daily reduction log-10 units per day (5C 0.011 log-10 units) versus
#0.151 for 20C; assumed produce are eaten within 7 days of harvesting
norodecaylow<-0.151

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
nvhG<-nvhH - nvhH*((1-Trhhg)*Phwg)
#NoV shedding for the packer - with global prevalence factor
shedP<-concp*prevp
#g feces on packer's hand 
mfhp<-10^(fhp/2)
#NoV on packer's hand - autoinnoculation
nvpA<-(mfhp*shedP)
#removed environmental contamination; virus on hand leaving restroom after conditional hand washing
nvpH<-nvpA/(10^(hweffp*Pphw))
#NoV on gloved packer hands
nvpG<-nvpH - nvpH*((1-Trphg)*Ppwg)
#NoV on produce after harvest (with transfer switch based on gloves y/n) updated per conversation with Kira
ifelse(Phwg>0, nvhP1<-nvhG- nvhG*((1-Trhgp)^nhPh), nvhP1<-nvhG- nvhG*((1-Trhhp)^nhPh))
#NoV on bin after contact with produce
nvBi<- nvhP1-nvhP1*((1-Trhpbi)^nPbi)
#NoV on produce after contacting bin
nvhP2<-nvhP1-nvBi
#Produce-produce dilution: NoV on produce in bin after produce contact
nvPPBi<- nvhP2-nvhP2*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP3<-nvhP2-nvPPBi
#NoV on belt after contacting produce
nvBe<- nvhP3-nvhP3*(1-Trppbe)^nPbe
#NoV on produce after contacting belt
nvpP1<-nvhP3-nvBe
#NoV on produce after packing by hand, nvpP2 is variable for NoV on produce after packing; update ifelse for glove switch per conversation with Kira
ifelse(Ppwg>0, 
       ifelse(nvpP1>nvpG, nvpP2<-nvpP1- (nvpP1-nvpP1*((1-Trpgp)^npPh)), nvpP2<-nvpP1+( nvpG-nvpG*((1-Trpgp)^npPh))),
       ifelse(nvpP1>nvpG, nvpP2<-nvpP1- (nvpP1-nvpP1*((1-Trpph)^npPh)), nvpP2<-nvpP1+( nvpG-nvpG*((1-Trphp)^npPh)))
)
#NoV on packer hand after packing produce (corrected by deleting extra -nvpH);update ifelse for glove switch per conversation with Kira
ifelse(Ppwg>0, 
       ifelse(nvpP1>nvpG, nvpH1<-(nvpP1-nvpP1*((1-Trpgp)^npPh))+nvpG, nvpH1<- nvpG*((1-Trpgp)^npPh)), 
       ifelse(nvpP1>nvpG, nvpH1<-(nvpP1-nvpP1*((1-Trpph)^npPh))+nvpG, nvpH1<- nvpG*((1-Trphp)^npPh))
)

#additional produce items (expanded)
#produce item #2
#NoV on harvesters hand after produce #1
nvhH1<-(nvhG-nvhP1)
#NoV on produce #2
ifelse(Phwg>0, nvhP12<-nvhH1-(nvhH1*((1-Trhgp)^nhPh)), nvhP12<-nvhH1- (nvhH1*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin after contact with produce
ifelse(nvhP12>nvBi, nvBi2<-(nvhP12-nvhP12*((1-Trhpbi)^nPbi))+nvBi, nvBi2<-nvBi*((1-Trhpbi)^nPbi))
#norovirus on produce after contacting bin
ifelse(nvhP12>nvBi, nvhP22<-nvhP12-(nvhP12-nvhP12*((1-Trhpbi)^nPbi)), nvhP22<-nvhP12+(nvBi-nvBi*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin (norovirus on produce in bin after produce contact)
nvPPBi2<-nvhP22-nvhP22*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP23<-nvhP22-nvPPBi2
#NoV on belt after contacting produce
ifelse(nvhP23>nvBe, nvBe2<-(nvhP23-nvhP23*((1-Trppbe)^nPbe))+nvBe, nvBe2<-nvBe*((1-Trppbe)^nPbe))
#NoV on produce after contacting belt
ifelse(nvhP23>nvBe, nvpP12<-nvhP23-(nvhP23-nvhP23*((1-Trppbe)^nPbe)), nvpP12<-nvhP23+(nvBe-nvBe*((1-Trppbe)^nPbe)))
#revised nvpH2 formula; final NoV on packer's hands and produce item
ifelse(Ppwg>0, 
       ifelse(nvpP12>nvpH1, nvpH2<-(nvpP12-nvpP12*((1-Trpgp)^npPh))+nvpH1, nvpH2<- nvpH1*((1-Trpgp)^npPh)),
       ifelse(nvpP12>nvpH1, nvpH2<-(nvpP12-nvpP12*((1-Trpph)^npPh))+nvpH1, nvpH2<- nvpH1*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 2
ifelse(Ppwg>0, 
       ifelse(nvpP12>nvpH1, nvpP22<-nvpP12-(nvpP12-nvpP12*((1-Trpgp)^npPh)), nvpP22<-nvpP12+( nvpH1-nvpH1*((1-Trpgp)^npPh))),
       ifelse(nvpP12>nvpH1, nvpP22<-nvpP12-(nvpP12-nvpP12*((1-Trpph)^npPh)), nvpP22<-nvpP12+( nvpH1-nvpH1*((1-Trphp)^npPh)))
)
#produce item #3
#NoV on harvesters hand
nvhH2<-nvhH1-nvhP12
#NoV on produce #3
ifelse(Phwg>0,nvhP13<-nvhH2- (nvhH2*((1-Trhgp)^nhPh)), nvhP13<-nvhH2- (nvhH2*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin after contact with produce
ifelse(nvhP13>nvBi2, nvBi3<-(nvhP13-nvhP13*((1-Trhpbi)^nPbi))+nvBi2, nvBi3<-nvBi2*((1-Trhpbi)^nPbi))
#bin interactions - norovirus on produce after contact with bin
ifelse(nvhP13>nvBi2, nvhP24<-nvhP13-(nvhP13-nvhP13*((1-Trhpbi)^nPbi)), nvhP24<-nvhP13+(nvBi2-nvBi2*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin - norovirus on produce in bin after produce contact
nvPPBi3<-nvhP24-nvhP24*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP25<-nvhP24-nvPPBi3
#belt interactions - norovirus on belt after produce contact
ifelse(nvhP25>nvBe2, nvBe3<-(nvhP25-nvhP25*((1-Trppbe)^nPbe))+nvBe2, nvBe3<-nvBe2*((1-Trppbe)^nPbe))
#belt interactions - norovirus on produce after belt contact
ifelse(nvhP25>nvBe2, nvpP13<-nvhP25-(nvhP25-nvhP25*((1-Trppbe)^nPbe)), nvpP13<-nvhP25+(nvBe2-nvBe2*((1-Trppbe)^nPbe)))
#revised nvpH3 formula, final norovirus on packer's hand
ifelse(Ppwg>0,
       ifelse(nvpP13>nvpH2, nvpH3<-(nvpP13-nvpP13*((1-Trpgp)^npPh))+nvpH2, nvpH3<- nvpH2*((1-Trpgp)^npPh)),
       ifelse(nvpP13>nvpH2, nvpH3<-(nvpP13-nvpP13*((1-Trpph)^npPh))+nvpH2, nvpH3<- nvpH2*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 3
ifelse(Ppwg>0,
       ifelse(nvpP13>nvpH2, nvpP23<-nvpP13- (nvpP13-nvpP13*((1-Trpgp)^npPh)), nvpP23<-nvpP13+( nvpH2-nvpH2*((1-Trpgp)^npPh))),
       ifelse(nvpP13>nvpH2, nvpP23<-nvpP13- (nvpP13-nvpP13*((1-Trpph)^npPh)), nvpP23<-nvpP13+( nvpH2-nvpH2*((1-Trphp)^npPh)))
)
#produce item #4
#NoV on harvesters hand
nvhH3<-nvhH2-nvhP13
#NoV on produce #4
ifelse(Phwg>0, nvhP14<-nvhH3- (nvhH3*((1-Trhgp)^nhPh)), nvhP14<-nvhH3- (nvhH3*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin after contact with produce
ifelse(nvhP14>nvBi3, nvBi4<-(nvhP14-nvhP14*((1-Trhpbi)^nPbi))+nvBi3, nvBi4<-nvBi3*((1-Trhpbi)^nPbi))
#bin interactions - norovirus on produce after contact with bin
ifelse(nvhP14>nvBi3, nvhP26<-nvhP14-(nvhP14-nvhP14*((1-Trhpbi)^nPbi)), nvhP26<-nvhP14+(nvBi3-nvBi3*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin - norovirus on produce in bin after produce contact
nvPPBi4<-nvhP26-nvhP26*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP27<-nvhP26-nvPPBi4
#belt interactions - norovirus on belt after produce contact
ifelse(nvhP27>nvBe3, nvBe4<-(nvhP27-nvhP27*((1-Trppbe)^nPbe))+nvBe3, nvBe4<-nvBe3*((1-Trppbe)^nPbe))
#belt interactions - norovirus on produce after belt contact
ifelse(nvhP27>nvBe3, nvpP14<-nvhP27-(nvhP27-nvhP27*((1-Trppbe)^nPbe)), nvpP14<-nvhP27+(nvBe3-nvBe3*((1-Trppbe)^nPbe)))
#revised nvpH4 formula, final norovirus on packer's hand 
ifelse(Ppwg>0,
       ifelse(nvpP14>nvpH3, nvpH4<-(nvpP14-nvpP14*((1-Trpgp)^npPh))+nvpH3, nvpH4<- nvpH3*((1-Trpgp)^npPh)),
       ifelse(nvpP14>nvpH3, nvpH4<-(nvpP14-nvpP14*((1-Trpph)^npPh))+nvpH3, nvpH4<- nvpH3*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 4
ifelse(Ppwg>0,
       ifelse(nvpP14>nvpH3, nvpP24<-nvpP14- (nvpP14-nvpP14*((1-Trpgp)^npPh)), nvpP24<-nvpP14+( nvpH3-nvpH3*((1-Trpgp)^npPh))),
       ifelse(nvpP14>nvpH3, nvpP24<-nvpP14- (nvpP14-nvpP14*((1-Trpph)^npPh)), nvpP24<-nvpP14+( nvpH3-nvpH3*((1-Trphp)^npPh)))
)
#produce item #5
#NoV on harvesters hand
nvhH4<-nvhH3-nvhP14
#NoV on produce #5
ifelse(Phwg>0, nvhP15<-nvhH4- (nvhH4*((1-Trhgp)^nhPh)), nvhP15<-nvhH4- (nvhH4*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin after contact with produce
ifelse(nvhP15>nvBi4, nvBi5<-(nvhP15-nvhP15*((1-Trhpbi)^nPbi))+nvBi4, nvBi5<-nvBi4*((1-Trhpbi)^nPbi))
#bin interactions - norovirus on produce after contact with bin
ifelse(nvhP15>nvBi4, nvhP28<-nvhP15-(nvhP15-nvhP15*((1-Trhpbi)^nPbi)), nvhP28<-nvhP15+(nvBi4-nvBi4*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin - norovirus on produce in bin after produce contact
nvPPBi5<-nvhP28-nvhP28*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP29<-nvhP28-nvPPBi5
#belt interactions - norovirus on belt after produce contact
ifelse(nvhP29>nvBe4, nvBe5<-(nvhP29-nvhP29*((1-Trppbe)^nPbe))+nvBe4, nvBe5<-nvBe4*((1-Trppbe)^nPbe))
#belt interactions - norovirus on produce after belt contact
ifelse(nvhP29>nvBe4, nvpP15<-nvhP29-(nvhP29-nvhP29*((1-Trppbe)^nPbe)), nvpP15<-nvhP29+(nvBe4-nvBe4*((1-Trppbe)^nPbe)))
#revised nvpH5 formula, final norovirus on packer's hand
ifelse(Ppwg>0,
       ifelse(nvpP15>nvpH4, nvpH5<-(nvpP15-nvpP15*((1-Trpgp)^npPh))+nvpH4, nvpH5<- nvpH4*((1-Trpgp)^npPh)),
       ifelse(nvpP15>nvpH4, nvpH5<-(nvpP15-nvpP15*((1-Trpph)^npPh))+nvpH4, nvpH5<- nvpH4*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 5
ifelse(Ppwg>0,
       ifelse(nvpP15>nvpH4, nvpP25<-nvpP15- (nvpP15-nvpP15*((1-Trpgp)^npPh)), nvpP25<-nvpP15+( nvpH4-nvpH4*((1-Trpgp)^npPh))),
       ifelse(nvpP15>nvpH4, nvpP25<-nvpP15- (nvpP15-nvpP15*((1-Trpph)^npPh)), nvpP25<-nvpP15+( nvpH4-nvpH4*((1-Trphp)^npPh)))
)
#produce item #6
#NoV on harvesters hand
nvhH5<-nvhH4-nvhP15
#NoV on produce #6
ifelse(Phwg>0, nvhP16<-nvhH5- (nvhH5*((1-Trhgp)^nhPh)), nvhP16<-nvhH5- (nvhH5*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin after contact with produce
ifelse(nvhP16>nvBi5, nvBi6<-(nvhP16-nvhP16*((1-Trhpbi)^nPbi))+nvBi5, nvBi6<-nvBi5*((1-Trhpbi)^nPbi))
#bin interactions - norovirus on produce after contact with bin
ifelse(nvhP16>nvBi5, nvhP30<-nvhP16-(nvhP16-nvhP16*((1-Trhpbi)^nPbi)), nvhP30<-nvhP16+(nvBi5-nvBi5*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin - norovirus on produce in bin after produce contact
nvPPBi6<-nvhP30-nvhP30*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP31<-nvhP30-nvPPBi6
#belt interactions - norovirus on belt after produce contact
ifelse(nvhP31>nvBe5, nvBe6<-(nvhP31-nvhP31*((1-Trppbe)^nPbe))+nvBe5, nvBe6<-nvBe5*((1-Trppbe)^nPbe))
#belt interactions - norovirus on produce after belt contact
ifelse(nvhP31>nvBe5, nvpP16<-nvhP31-(nvhP31-nvhP31*((1-Trppbe)^nPbe)), nvpP16<-nvhP31+(nvBe5-nvBe5*((1-Trppbe)^nPbe)))
#revised nvpH6 formula, final norovirus on packer's hand
ifelse(Ppwg>0,
       ifelse(nvpP16>nvpH5, nvpH6<-(nvpP16-nvpP16*((1-Trpgp)^npPh))+nvpH5, nvpH6<- nvpH5*((1-Trpgp)^npPh)),
       ifelse(nvpP16>nvpH5, nvpH6<-(nvpP16-nvpP16*((1-Trpph)^npPh))+nvpH5, nvpH6<- nvpH5*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 6
ifelse(Ppwg>0,
       ifelse(nvpP16>nvpH5, nvpP26<-nvpP16- (nvpP16-nvpP16*((1-Trpgp)^npPh)), nvpP26<-nvpP16+( nvpH5-nvpH5*((1-Trpgp)^npPh))),
       ifelse(nvpP16>nvpH5, nvpP26<-nvpP16- (nvpP16-nvpP16*((1-Trpph)^npPh)), nvpP26<-nvpP16+( nvpH5-nvpH5*((1-Trphp)^npPh)))
)
#produce item #7
#NoV on harvesters hand
nvhH6<-nvhH5-nvhP16
#NoV on produce #7
ifelse(Phwg>0,nvhP17<-nvhH6- (nvhH6*((1-Trhgp)^nhPh)), nvhP17<-nvhH6- (nvhH6*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin after contact with produce
ifelse(nvhP17>nvBi6, nvBi7<-(nvhP17-nvhP17*((1-Trhpbi)^nPbi))+nvBi6, nvBi7<-nvBi6*((1-Trhpbi)^nPbi))
#bin interactions - norovirus on produce after contact with bin
ifelse(nvhP17>nvBi6, nvhP32<-nvhP17-(nvhP17-nvhP17*((1-Trhpbi)^nPbi)), nvhP32<-nvhP17+(nvBi6-nvBi6*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin - norovirus on produce in bin after produce contact
nvPPBi7<-nvhP32-nvhP32*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP33<-nvhP32-nvPPBi7
#belt interactions - norovirus on belt after produce contact
ifelse(nvhP33>nvBe6, nvBe7<-(nvhP33-nvhP33*((1-Trppbe)^nPbe))+nvBe6, nvBe7<-nvBe6*((1-Trppbe)^nPbe))
#belt interactions - norovirus on produce after belt contact
ifelse(nvhP33>nvBe6, nvpP17<-nvhP33-(nvhP33-nvhP33*((1-Trppbe)^nPbe)), nvpP17<-nvhP33+(nvBe6-nvBe6*((1-Trppbe)^nPbe)))
#revised nvpH7 formula, final norovirus on packer's hand
ifelse(Ppwg>0,
       ifelse(nvpP17>nvpH6, nvpH7<-(nvpP17-nvpP17*((1-Trpgp)^npPh))+nvpH6, nvpH7<- nvpH6*((1-Trpgp)^npPh)),
       ifelse(nvpP17>nvpH6, nvpH7<-(nvpP17-nvpP17*((1-Trpph)^npPh))+nvpH6, nvpH7<- nvpH6*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 7
ifelse(Ppwg>0,
       ifelse(nvpP17>nvpH6, nvpP27<-nvpP17- (nvpP17-nvpP17*((1-Trpgp)^npPh)), nvpP27<-nvpP17+( nvpH6-nvpH6*((1-Trpgp)^npPh))),
       ifelse(nvpP17>nvpH6, nvpP27<-nvpP17- (nvpP17-nvpP17*((1-Trpph)^npPh)), nvpP27<-nvpP17+( nvpH6-nvpH6*((1-Trphp)^npPh)))
)
#produce item #8
#NoV on harvesters hand
nvhH7<-nvhH6-nvhP17
#NoV on produce #8
ifelse(Phwg>0,nvhP18<-nvhH7- (nvhH7*((1-Trhgp)^nhPh)),nvhP18<-nvhH7- (nvhH7*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin after contact with produce
ifelse(nvhP18>nvBi7, nvBi8<-(nvhP18-nvhP18*((1-Trhpbi)^nPbi))+nvBi7, nvBi8<-nvBi7*((1-Trhpbi)^nPbi))
#bin interactions - norovirus on produce after contact with bin
ifelse(nvhP18>nvBi7, nvhP34<-nvhP18-(nvhP18-nvhP18*((1-Trhpbi)^nPbi)), nvhP34<-nvhP18+(nvBi7-nvBi7*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin - norovirus on produce in bin after produce contact
nvPPBi8<-nvhP34-nvhP34*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP35<-nvhP34-nvPPBi8
#belt interactions - norovirus on belt after produce contact
ifelse(nvhP35>nvBe7, nvBe8<-(nvhP35-nvhP35*((1-Trppbe)^nPbe))+nvBe7, nvBe8<-nvBe7*((1-Trppbe)^nPbe))
#belt interactions - norovirus on produce after belt contact
ifelse(nvhP35>nvBe7, nvpP18<-nvhP35-(nvhP35-nvhP35*((1-Trppbe)^nPbe)), nvpP18<-nvhP35+(nvBe7-nvBe7*((1-Trppbe)^nPbe)))
#revised nvpH8 formula, final norovirus on packer's hand
ifelse(Ppwg>0,
       ifelse(nvpP18>nvpH7, nvpH8<-(nvpP18-nvpP18*((1-Trpgp)^npPh))+nvpH7, nvpH8<- nvpH7*((1-Trpgp)^npPh)),
       ifelse(nvpP18>nvpH7, nvpH8<-(nvpP18-nvpP18*((1-Trpph)^npPh))+nvpH7, nvpH8<- nvpH7*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 8
ifelse(Ppwg>0,
       ifelse(nvpP18>nvpH7, nvpP28<-nvpP18- (nvpP18-nvpP18*((1-Trpgp)^npPh)), nvpP28<-nvpP18+( nvpH7-nvpH7*((1-Trpgp)^npPh))),
       ifelse(nvpP18>nvpH7, nvpP28<-nvpP18- (nvpP18-nvpP18*((1-Trpph)^npPh)), nvpP28<-nvpP18+( nvpH7-nvpH7*((1-Trphp)^npPh)))
)
#produce item #9
#NoV on harvesters hand
nvhH8<-nvhH7-nvhP18
#NoV on produce #9
ifelse(Phwg>0,nvhP19<-nvhH8- (nvhH8*((1-Trhgp)^nhPh)),nvhP19<-nvhH8- (nvhH8*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin after contact with produce
ifelse(nvhP19>nvBi8, nvBi9<-(nvhP19-nvhP19*((1-Trhpbi)^nPbi))+nvBi8, nvBi9<-nvBi8*((1-Trhpbi)^nPbi))
#bin interactions - norovirus on produce after contact with bin
ifelse(nvhP19>nvBi8, nvhP36<-nvhP19-(nvhP19-nvhP19*((1-Trhpbi)^nPbi)), nvhP36<-nvhP19+(nvBi8-nvBi8*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin - norovirus on produce in bin after produce contact
nvPPBi9<-nvhP36-nvhP36*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP37<-nvhP36-nvPPBi9
#belt interactions - norovirus on belt after produce contact
ifelse(nvhP37>nvBe8, nvBe9<-(nvhP37-nvhP37*((1-Trppbe)^nPbe))+nvBe8, nvBe9<-nvBe8*((1-Trppbe)^nPbe))
#belt interactions - norovirus on produce after belt contact
ifelse(nvhP37>nvBe8, nvpP19<-nvhP37-(nvhP37-nvhP37*((1-Trppbe)^nPbe)), nvpP19<-nvhP37+(nvBe8-nvBe8*((1-Trppbe)^nPbe)))
#revised nvpH9 formula, final norovirus on packer's hand
ifelse(Ppwg>0,
       ifelse(nvpP19>nvpH8, nvpH9<-(nvpP19-nvpP19*((1-Trpgp)^npPh))+nvpH8, nvpH9<- nvpH8*((1-Trpgp)^npPh)),
       ifelse(nvpP19>nvpH8, nvpH9<-(nvpP19-nvpP19*((1-Trpph)^npPh))+nvpH8, nvpH9<- nvpH8*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 9
ifelse(Ppwg>0,
       ifelse(nvpP19>nvpH8, nvpP29<-nvpP19- (nvpP19-nvpP19*((1-Trpgp)^npPh)), nvpP29<-nvpP19+( nvpH8-nvpH8*((1-Trpgp)^npPh))),
       ifelse(nvpP19>nvpH8, nvpP29<-nvpP19- (nvpP19-nvpP19*((1-Trpph)^npPh)), nvpP29<-nvpP19+( nvpH8-nvpH8*((1-Trphp)^npPh)))
)
#produce item #10
#NoV on harvesters hand
nvhH9<-nvhH8-nvhP19
#NoV on produce #10
ifelse(Phwg>0,nvhP110<-nvhH9- (nvhH9*((1-Trhgp)^nhPh)), nvhP110<-nvhH9- (nvhH9*((1-Trhhp)^nhPh)))
#bin interactions - norovirus on bin afetr contact with produce
ifelse(nvhP110>nvBi9, nvBi10<-(nvhP110-nvhP110*((1-Trhpbi)^nPbi))+nvBi9, nvBi10<-nvBi9*((1-Trhpbi)^nPbi))
#bin interactions - norovirus on produce after contact with bin
ifelse(nvhP110>nvBi9, nvhP38<-nvhP110-(nvhP110-nvhP110*((1-Trhpbi)^nPbi)), nvhP38<-nvhP110+(nvBi9-nvBi9*((1-Trhpbi)^nPbi)))
#produce-produce dilution in bin
nvPPBi10<-nvhP38-nvhP38*((1-Trhppbi)^nPPbi)
#NoV on produce item after dilution to produce in bin
nvhP39<-nvhP38-nvPPBi10
#belt interactions - norovirus on belt after produce contact
ifelse(nvhP39>nvBe9, nvBe10<-(nvhP39-nvhP39*((1-Trppbe)^nPbe))+nvBe9, nvBe10<-nvBe9*((1-Trppbe)^nPbe))
#belt interactions - norovirus on produce after belt contact
ifelse(nvhP39>nvBe9, nvpP110<-nvhP39-(nvhP39-nvhP39*((1-Trppbe)^nPbe)), nvpP110<-nvhP39+(nvBe9-nvBe9*((1-Trppbe)^nPbe)))
#revised nvpH10 formula, final norovirus on packer's hand
ifelse(Ppwg>0,
       ifelse(nvpP110>nvpH9, nvpH10<-(nvpP110-nvpP110*((1-Trpgp)^npPh))+nvpH9, nvpH10<- nvpH9*((1-Trpgp)^npPh)),
       ifelse(nvpP110>nvpH9, nvpH10<-(nvpP110-nvpP110*((1-Trpph)^npPh))+nvpH9, nvpH10<- nvpH9*((1-Trphp)^npPh))
)
#NoV on final packed produce item number 10
ifelse(Ppwg>0,
       ifelse(nvpP110>nvpH9, nvpP210<-nvpP110- (nvpP110-nvpP110*((1-Trpgp)^npPh)), nvpP210<-nvpP110+( nvpH9-nvpH9*((1-Trpgp)^npPh))),
       ifelse(nvpP110>nvpH9, nvpP210<-nvpP110- (nvpP110-nvpP110*((1-Trpph)^npPh)), nvpP210<-nvpP110+( nvpH9-nvpH9*((1-Trphp)^npPh)))
)

#calculating total viral load per produce items factoring in viral decay at 5C/20C, assumed 7 days total of viral decay
nvpP2decayl <- nvpP2/(10^(norodecaylow*7))
nvpP22decayl <- nvpP22/(10^(norodecaylow)*7)
nvpP23decayl <- nvpP23/(10^(norodecaylow)*7)
nvpP24decayl <- nvpP24/(10^(norodecaylow)*7)
nvpP25decayl <- nvpP25/(10^(norodecaylow)*7)
nvpP26decayl <- nvpP26/(10^(norodecaylow)*7)
nvpP27decayl <- nvpP27/(10^(norodecaylow)*7)
nvpP28decayl <- nvpP28/(10^(norodecaylow)*7)
nvpP29decayl <- nvpP29/(10^(norodecaylow)*7)
nvpP210decayl <-nvpP210/(10^(norodecaylow)*7)

#body weight pulled for adults from exposure handbook, kg)
bodyweight<-80

#USDA reported weights in grams per produce item
#https://hannaone.com/Recipe/weightcantaloupe.html
#https://hannaone.com/Recipe/weightlettuce.html
#https://hannaone.com/Recipe/weighttomato.html
weightlettucehead<-mcstoc(rtriang,type="V", min=309,mode=360,max=626)
weightmelon<-mcstoc(rtriang,type="V", min=441,mode=552,max=814)
weighttomato<-mcstoc(rtriang,type="V", min=91,mode=123,max=182)
weightberry<-mcstoc(rtriang,type="V", min=3,mode=4,max=5)

#Exposure handbook intake of specific fruits and veggies (daily consumption g/kg-day)
lettuceconsumpt<-mcstoc(rnorm,type="V",mean=0.44, sd=0.01)
melonconsumpt<-mcstoc(rnorm,type="V",mean=0.70, sd=0.05)
tomatoconsumpt<-mcstoc(rnorm,type="V",mean=0.83, sd=0.02)
berryconsumpt<-mcstoc(rnorm,type="V",mean=0.45, sd=0.02)

#dose calculation from TOMATO ending with total virions/day
dose1tomato <- (bodyweight*tomatoconsumpt*nvpP2decayl)/weighttomato
dose2tomato <- (bodyweight*tomatoconsumpt*nvpP22decayl)/weighttomato
dose3tomato <- (bodyweight*tomatoconsumpt*nvpP23decayl)/weighttomato
dose4tomato <- (bodyweight*tomatoconsumpt*nvpP24decayl)/weighttomato
dose5tomato <- (bodyweight*tomatoconsumpt*nvpP25decayl)/weighttomato
dose6tomato <- (bodyweight*tomatoconsumpt*nvpP26decayl)/weighttomato
dose7tomato <- (bodyweight*tomatoconsumpt*nvpP27decayl)/weighttomato
dose8tomato <- (bodyweight*tomatoconsumpt*nvpP28decayl)/weighttomato
dose9tomato <- (bodyweight*tomatoconsumpt*nvpP29decayl)/weighttomato
dose10tomato <- (bodyweight*tomatoconsumpt*nvpP210decayl)/weighttomato

#calculating risk from TOMATO with fractional Poisson model aggregated viral decay at 5C

riskfracpoissT1 <- P*(1-exp(-dose1tomato/mu))
riskfracpoissT2 <- P*(1-exp(-dose2tomato/mu))
riskfracpoissT3 <- P*(1-exp(-dose3tomato/mu))
riskfracpoissT4 <- P*(1-exp(-dose4tomato/mu))
riskfracpoissT5 <- P*(1-exp(-dose5tomato/mu))
riskfracpoissT6 <- P*(1-exp(-dose6tomato/mu))
riskfracpoissT7 <- P*(1-exp(-dose7tomato/mu))
riskfracpoissT8 <- P*(1-exp(-dose8tomato/mu))
riskfracpoissT9 <- P*(1-exp(-dose9tomato/mu))
riskfracpoissT10 <- P*(1-exp(-dose10tomato/mu))

###putting everything through the mc function
#sensitivity for risk tomato
QMRAtomato<-mc(conc,fh,hweff, Phhw,Phwg, prevh, concp, fhp, hweffp, Pphw, Ppwg, prevp, Trhhp, Trhgp,Trpph,Trphp,Trpgp, P, mu,
              Trhpbi,Trppbe, Trhhg,Trphg,Trhppbi,nhPh,npPh,nPbi,nPbe,nPPbi,shedH, mfh, nvhA, nvhH, nvhG, shedP, mfhp, nvpA, nvpH, nvpG, nvhP1,
              nvBi,nvhP2,nvPPBi,nvhP3,nvBe,nvpP1,nvpH1,nvpP2, nvpP22,nvpP23,nvpP24,nvpP25,nvpP26,nvpP27,nvpP28,nvpP29,nvpP210,nvpP2decayl,nvpP22decayl,nvpP23decayl,nvpP24decayl,nvpP25decayl,
              nvpP26decayl, nvpP27decayl, nvpP28decayl,nvpP29decayl,nvpP210decayl,weightlettucehead, weightmelon, weighttomato, weightberry, lettuceconsumpt,melonconsumpt, tomatoconsumpt, berryconsumpt, 
              dose1tomato,dose2tomato, dose3tomato, dose4tomato, dose5tomato, dose6tomato, dose7tomato, dose8tomato, dose9tomato, dose10tomato,
              riskfracpoissT2, riskfracpoissT3, riskfracpoissT4, riskfracpoissT5, riskfracpoissT6, riskfracpoissT7, riskfracpoissT8, riskfracpoissT9, riskfracpoissT1)
summary(QMRAtomato)
summary(QMRAtomato$riskfracpoissT1)
summary(QMRAtomato$nvpP2)


tortomato<-tornado(QMRAtomato)
plot(tortomato)
print(tortomato)
mcratio(QMRAtomato)
unctomato=tornadounc(QMRAtomato)
unctomato
####################################################################################################################################################################################################################################################################


