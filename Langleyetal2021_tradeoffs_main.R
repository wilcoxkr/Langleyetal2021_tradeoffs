### Script for analyses reported in Langley et al. 2021 Ecology
### Authors: Adam Langley (adam.langley@villanova.edu) and Emily Grman (egrman@emich.edu)

### DOI to the dataset used in the script below is https://doi.org/10.5061/dryad.rfj6q57c1

library(tidyverse)
library(reshape2)
library(ggpubr)
library(lme4)
library(grid)

theme_set(theme_bw())
theme_eg=theme_update(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), strip.background=element_rect(color="white", fill="white"))

#read in data
effect.sizes=read.csv("Langley et al 2022 Ecology data.csv")

#making nice labels for treatment names
trt.labels=data.frame(global.change.treatment=unique(effect.sizes$global.change.treatment))
trt.labels$trt.labels=as.factor(trt.labels$global.change.treatment)
levels(trt.labels$trt.labels)=c("Disturbance", "CO2", "Drought", "Irrigation", "Multiple Nutrients", "Nitrogen", "Phosphorus", "Temperature")

#summary statistics on effect size E:
E.real=effect.sizes[effect.sizes$permutation=="real" & !is.na(effect.sizes$E),]
E.real=merge(E.real, trt.labels, by="global.change.treatment")
E.real$respE=as.factor(ifelse(E.real$E>0, "pos", ifelse(E.real$E<0, "neg", "null")))
tableS2=dcast(E.real, trt.labels~respE, value.var="E", fun=length)
tableS2$number=tableS2$neg+tableS2$null+tableS2$pos
tableS2$prop.pos=round(tableS2$pos/tableS2$number, 2)
tableS2$prop.neg=round(tableS2$neg/tableS2$number, 2); tableS2
dcast(E.real, .~respE, value.var="E", fun=length)
1830+61+1900
1830/3791 #proportion negative responses
1900/3791 #proportion positive responses


#Figure S1:
alltrt=qplot(E, data=E.real, ylab="Number of species", xlab="Effect of any treatment (E)", xlim=c(-1.1, 1.1), main="a)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
BMCT=qplot(E, data=E.real[E.real$trt.labels=="Disturbance",], ylab="Number of species", xlab="Effect of disturbance (E)", xlim=c(-1.1, 1.1), main="b)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
CO2=qplot(E, data=E.real[E.real$trt.labels=="CO2",], ylab="Number of species", xlab="Effect of CO2 (E)", xlim=c(-1.1, 1.1), main="c)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
drought=qplot(E, data=E.real[E.real$trt.labels=="Drought",], ylab="Number of species", xlab="Effect of drought (E)", xlim=c(-1.1, 1.1), main="d)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
irr=qplot(E, data=E.real[E.real$trt.labels=="Irrigation",], ylab="Number of species", xlab="Effect of irrigation (E)", xlim=c(-1.1, 1.1), main="e)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
multnut=qplot(E, data=E.real[E.real$trt.labels=="Multiple Nutrients",], ylab="Number of species", xlab="Effect of multiple nutrients (E)", xlim=c(-1.1, 1.1), main="f)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
nitrogen=qplot(E, data=E.real[E.real$trt.labels=="Nitrogen",], ylab="Number of species", xlab="Effect of nitrogen (E)", xlim=c(-1.1, 1.1), main="g)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
phos=qplot(E, data=E.real[E.real$trt.labels=="Phosphorus",], ylab="Number of species", xlab="Effect of phosphorus (E)", xlim=c(-1.1, 1.1), main="h)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
temp=qplot(E, data=E.real[E.real$trt.labels=="Temperature",], ylab="Number of species", xlab="Effect of temperature (E)", xlim=c(-1.1, 1.1), main="i)", fill=code) + scale_fill_manual(values=c("gray50", "black", "gray35")) + theme(legend.position="none")
ggarrange(alltrt, BMCT, CO2, drought, irr, multnut, nitrogen, phos, temp, ncol=3, nrow=3)
ggsave("figs/Fig S1.pdf", width=8, height=6)


#----------WHICH TRADEOFFS DO WE HAVE THE DATA TO TEST?

trtsinsite=dcast(effect.sizes, site_project_comm~global.change.treatment, value.var="E", fun=length)

trtsinsite$BMCTxCO2=ifelse(trtsinsite$BMCT>0 & trtsinsite$CO2>0, 1, 0)
trtsinsite$BMCTxdrought=ifelse(trtsinsite$BMCT>0 & trtsinsite$drought>0, 1, 0)
trtsinsite$BMCTxirr=ifelse(trtsinsite$BMCT>0 & trtsinsite$irr>0, 1, 0)
trtsinsite$BMCTxmult_nutrient=ifelse(trtsinsite$BMCT>0 & trtsinsite$mult_nutrient>0, 1, 0)
trtsinsite$BMCTxN=ifelse(trtsinsite$BMCT>0 & trtsinsite$N>0, 1, 0)
trtsinsite$BMCTxP=ifelse(trtsinsite$BMCT>0 & trtsinsite$P>0, 1, 0)
trtsinsite$BMCTxtemp=ifelse(trtsinsite$BMCT>0 & trtsinsite$temp>0, 1, 0)
trtsinsite$CO2xdrought=ifelse(trtsinsite$CO2>0 & trtsinsite$drought>0, 1, 0)
trtsinsite$CO2xirr=ifelse(trtsinsite$CO2>0 & trtsinsite$irr>0, 1, 0)
trtsinsite$CO2xmult_nutrient=ifelse(trtsinsite$CO2>0 & trtsinsite$mult_nutrient>0, 1, 0)
trtsinsite$CO2xN=ifelse(trtsinsite$CO2>0 & trtsinsite$N>0, 1, 0)
trtsinsite$CO2xP=ifelse(trtsinsite$CO2>0 & trtsinsite$P>0, 1, 0)
trtsinsite$CO2xtemp=ifelse(trtsinsite$CO2>0 & trtsinsite$temp>0, 1, 0)
trtsinsite$droughtxirr=ifelse(trtsinsite$drought>0 & trtsinsite$irr>0, 1, 0)
trtsinsite$droughtxmult_nutrient=ifelse(trtsinsite$drought>0 & trtsinsite$mult_nutrient>0, 1, 0)
trtsinsite$droughtxN=ifelse(trtsinsite$drought>0 & trtsinsite$N>0, 1, 0)
trtsinsite$droughtxP=ifelse(trtsinsite$drought>0 & trtsinsite$P>0, 1, 0)
trtsinsite$droughtxtemp=ifelse(trtsinsite$drought>0 & trtsinsite$temp>0, 1, 0)
trtsinsite$irrxmult_nutrient=ifelse(trtsinsite$irr>0 & trtsinsite$mult_nutrient>0, 1, 0)
trtsinsite$irrxN=ifelse(trtsinsite$irr>0 & trtsinsite$N>0, 1, 0)
trtsinsite$irrxP=ifelse(trtsinsite$irr>0 & trtsinsite$P>0, 1, 0)
trtsinsite$irrxtemp=ifelse(trtsinsite$irr>0 & trtsinsite$temp>0, 1, 0)
trtsinsite$mult_nutrientxN=ifelse(trtsinsite$mult_nutrient>0 & trtsinsite$N>0, 1, 0)
trtsinsite$mult_nutrientxP=ifelse(trtsinsite$mult_nutrient>0 & trtsinsite$P>0, 1, 0)
trtsinsite$mult_nutrientxtemp=ifelse(trtsinsite$mult_nutrient>0 & trtsinsite$temp>0, 1, 0)
trtsinsite$NxP=ifelse(trtsinsite$N>0 & trtsinsite$P>0, 1, 0)
trtsinsite$Nxtemp=ifelse(trtsinsite$N>0 & trtsinsite$temp>0, 1, 0)
trtsinsite$Pxtemp=ifelse(trtsinsite$P>0 & trtsinsite$temp>0, 1, 0)

#how many studies for each tradeoff? formatting for use later in loops
temp=subset(trtsinsite, select=-c(site_project_comm, BMCT, CO2, drought, irr, mult_nutrient, N, P, temp))
temp2=colSums(temp)
temp2=data.frame(studies=temp2)
temp2$tradeoffs=row.names(temp2)
trtcombos.avail=temp2[temp2$studies>0,]
temp=matrix(unlist(strsplit(as.character(unique(trtcombos.avail$tradeoffs)), "x")), ncol=2, byrow=T)
trtcombos.avail$var1=temp[,1]
trtcombos.avail$var1pretty=as.factor(trtcombos.avail$var1)
levels(trtcombos.avail$var1pretty)=c("Disturbance", "CO2", "Drought", "Irrigation", "Multiple Nutrients", "Nitrogen")
trtcombos.avail$var2=temp[,2]
trtcombos.avail$var2pretty=as.factor(trtcombos.avail$var2)
levels(trtcombos.avail$var2pretty)=c("Drought", "Irrigation", "Multiple Nutrients", "Nitrogen", "Phosphorus", "Temperature")
trtcombos.avail.m=melt(trtcombos.avail, id="tradeoffs", measure=c("var1", "var2"))
names(trtcombos.avail.m)=c("tradeoff", "variable", "trt")
ntradeoffs.all=unique(trtcombos.avail.m$tradeoff)

#identifying which tradeoffs can be tested in which site(s)
trtcombos.m=melt(trtsinsite, id="site_project_comm", measure=c(row.names(trtcombos.avail)))
trtcombos.m=trtcombos.m[trtcombos.m$value>0,]


#---------COMPARING SIMULATED DATA TO REAL DATA-----ALL SPECIES


#looping through all tradeoffs we want
perm.list=unique(effect.sizes$permutation)
output.lm.E.allsp=numeric(0)
datallreal.E=numeric(0)
datallsim.E=numeric(0)

for (h in 1:length(ntradeoffs.all)) {
  trdf=as.data.frame(trtcombos.m[trtcombos.m$variable==as.character(ntradeoffs.all[h]), c("site_project_comm")])
  names(trdf)="site_project_comm"
  dath=merge(effect.sizes, trdf, by="site_project_comm")
  names(dath)[names(dath)=="global.change.treatment"]="trt"
  names(dath)[names(dath)=="E"]="eff"
  trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff==as.character(ntradeoffs.all[h]),]
  dath=merge(dath, trdf.vars.m, by="trt")
  want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
  wantsim=want[!want$permutation=="real",]
  datallsim.E=rbind(datallsim.E, wantsim)
  wantreal=want[want$permutation=="real",]
  datallreal.E=rbind(datallreal.E, wantreal)
  nsites=unique(want$site_project_comm)
  temp.i=numeric(0)
  
  for(i in 1:length(perm.list)) {
    dati=want[want$permutation==as.character(perm.list[i]),]
    temp.j=numeric(0)
    
    for(j in 1:length(nsites)) {
      datj=dati[dati$site_project_comm==as.character(nsites[j]),]
      datj=datj[complete.cases(datj[,c("var1", "var2")]),]
      datj$QI=ifelse(datj$var1>0 & datj$var2>0, 1, 0)
      datj$QII=ifelse(datj$var1<0 & datj$var2>0, 1, 0)
      datj$QIII=ifelse(datj$var1<0 & datj$var2<0, 1, 0)
      datj$QIV=ifelse(datj$var1>0 & datj$var2<0, 1, 0)
      nspecies=length(datj$species)
      temp.j2=data.frame(site_project_comm=as.character(nsites[j]), prop.QI=sum(datj$QI)/nspecies, prop.QII=sum(datj$QII)/nspecies, prop.QIII=sum(datj$QIII)/nspecies, prop.QIV=sum(datj$QIV)/nspecies)
      temp.j=rbind(temp.j, temp.j2)
    }
    
    temp.i2=data.frame(temp.j, permutation=as.character(perm.list[i]))
    temp.i=rbind(temp.i, temp.i2)
    #collect all that for each permutation and label with the permutation
  }
  
  #loop through sites to rank observed proportion in quadrants against all permuted proportions
  temp.h=numeric(0)
  for(g in 1:length(nsites)) {
    dat=temp.i[temp.i$site_project_comm==as.character(nsites[g]),]
    dat$prop.correspondence=dat$prop.QI + dat$prop.QIII
    dat$prop.tradeoffs=dat$prop.QII + dat$prop.QIV
    dat$QI.ranks=rank(dat$prop.QI)
    dat$QII.ranks=rank(dat$prop.QII)
    dat$QIII.ranks=rank(dat$prop.QIII)
    dat$QIV.ranks=rank(dat$prop.QIV)	
    dat$correspondence.ranks=rank(dat$prop.correspondence)
    dat$tradeoffs.ranks=rank(dat$prop.tradeoffs)
    temp.h=rbind(temp.h, dat)
    
  }
  
  temp.h2=data.frame(temp.h, tradeoff=as.character(ntradeoffs.all[h]))
  output.lm.E.allsp=rbind(output.lm.E.allsp, temp.h2)
  #collect all that for each tradeoff and label with the tradeoff
}

output.E=melt(output.lm.E.allsp, id=c("site_project_comm", "permutation", "tradeoff"))

#getting nice labels for plotting later:
tradeoff.labels=data.frame(tradeoff=unique(output.E$tradeoff))
tradeoff.labels$tradeoff2=as.factor(tradeoff.labels$tradeoff)
levels(tradeoff.labels$tradeoff2)=c("Disturbance x Drought", "Disturbance x Irrigation", "Disturbance x Nitrogen", "Disturbance x Phosphorus", "Disturbance x Temperature", "CO2 x Irrigation", "CO2 x Nitrogen", "CO2 x Temperature", "Drought x Irrigation", "Drought x Nitrogen", "Drought x Temperature", "Irrigation x Multiple Nutrients", "Irrigation x Nitrogen", "Irrigation x Phosphorus", "Irrigation x Temperature", "Multiple Nutrients x Temperature", "Nitrogen x Phosphorus", "Nitrogen x Temperature")
#getting sample sizes (number of studies) for each tradeoff:
tradeoff.n=dcast(output.E[output.E$permutation=="real",], tradeoff~variable, fun=length)
tradeoff.n=data.frame(tradeoff=tradeoff.n$tradeoff, numberstudies=tradeoff.n$prop.QI)
tradeoff.labels=merge(tradeoff.labels, tradeoff.n, by="tradeoff")


# FIGURE 2: species responses to treatment pairs

#first cleaning up long list of species responses to pairs of treatments for each tradeoff:
datallreal.E=datallreal.E[complete.cases(datallreal.E[,c("var1", "var2")]),]
datallreal.E=merge(datallreal.E, tradeoff.labels, by="tradeoff")

#adding in mean abundance in control plots
control.abund=dcast(effect.sizes[effect.sizes$permutation=="real",], site_project_comm+species~., value.var="relabund.control", fun=mean)
names(control.abund)[names(control.abund)=="."]<-"relabund.control"
datallreal.E=merge(datallreal.E, control.abund, by=c("site_project_comm", "species"))

#now cleaning up data for simulated communities:
datallsim.E=datallsim.E[complete.cases(datallsim.E[,c("var1", "var2")]),]
datallsim.E=merge(datallsim.E, tradeoff.labels, by="tradeoff")

#Figure 2: 
ggplot(data=datallreal.E, aes(var1, var2)) + geom_point(aes(size=relabund.control, alpha=I(0.2))) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + xlab("Effect of treatment 1 (E)") + ylab("Effect of treatment 2 (E)") + theme(legend.position="none") 
ggsave("figs/Fig 2.pdf", height=6, width=6)

#calculating proportion of species in 4 quadrants (text in fig 2):
datallreal.E$quad=as.factor(ifelse(datallreal.E$var1<0 & datallreal.E$var2<0, "QIII", ifelse(datallreal.E$var1>0 & datallreal.E$var2>0, "QI", ifelse(datallreal.E$var1>0 & datallreal.E$var2<0, "QIV", "QII"))))
nosp.quad.real=dcast(datallreal.E, quad~., value.var="permutation", fun=length) #number of species in each quadrant
datallsim.E$quad=as.factor(ifelse(datallsim.E$var1<0 & datallsim.E$var2<0, "QIII", ifelse(datallsim.E$var1>0 & datallsim.E$var2>0, "QI", ifelse(datallsim.E$var1>0 & datallsim.E$var2<0, "QIV", "QII"))))
nosp.quad.simt=dcast(datallsim.E, permutation+quad~., value.var="permutation", fun=length) #number of species in each quadrant
nosp.quad.sim=dcast(nosp.quad.simt, quad~., value.var=".", fun=mean)
nosp.quad=merge(nosp.quad.real, nosp.quad.sim, by="quad")
names(nosp.quad)=c("quad", "real", "mean.sim")
nosp.quad$real.lab=paste("O:", nosp.quad$real, sep="")
nosp.quad$sim.lab=paste("E:", round(nosp.quad$mean.sim, 0), sep="")
nosp.quad

#not shown: real data by treatment combination
ggplot(data=datallreal.E, aes(var1, var2)) + geom_point(aes(size=relabund.control, alpha=I(0.2))) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + xlab("Effect of treatment 1 (E)") + ylab("Effect of treatment 2 (E)") + theme(legend.position="none") + facet_wrap(~tradeoff2, ncol=3) 
ggsave("figs/fig not shown E.pdf", height=12, width=6.5)

#calculating proportion of species in 4 quadrants (text in fig not shown):
output.realprop=dcast(output.E[output.E$permutation=="real" & output.E$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs"),], site_project_comm+tradeoff+variable~.)
names(output.realprop)[names(output.realprop)=="."]<-"obsprop"
output.simprop=dcast(output.E[output.E$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs") & !output.E$permutation=="real",], site_project_comm+tradeoff+variable~., mean)
names(output.simprop)[names(output.simprop)=="."]<-"simulatedprop"
temp=dcast(output.simprop, tradeoff~variable, value.var="simulatedprop", fun=mean); temp 
output.proportions=merge(output.realprop, output.simprop, by=c("site_project_comm", "tradeoff", "variable"))
names(output.proportions)[names(output.proportions)=="variable"]<-"quadrant"
output.proportions$dif.in.prop=output.proportions$obsprop-output.proportions$simulatedprop; output.proportions


#----------------COMPARING OBSERVED TO EXPECTED PROPORTIONS

# FIGURE 3: histogram of observed-expected proportions

ggplot(prop.bystudy, aes(x=dif.in.prop)) + geom_histogram(binwidth=0.04) + geom_vline(xintercept=mean(prop.bystudy$dif.in.prop)) + geom_vline(xintercept=0, lty=2) + xlab("Observed - Expected proportion\nin Quadrants II and IV") + ylab("Number of studies") + geom_text(data=ttestlab, aes(x=-Inf, y=Inf, label="p<0.001", hjust=-0.5, vjust=1.5), size=3)
ggsave("figs/fig 3.pdf", width=4, height=3)

#t-test omitting SERC (strongest tradeoff)
t.test(prop.bystudy$obsprop[!prop.bystudy$site_project_comm=="SERC_CXN_0"], prop.bystudy$simulatedprop[!prop.bystudy$site_project_comm=="SERC_CXN_0"], paired=T) #p<0.001


# FIGURE 4: quilt of quadrants

#---ttests for individual treatment combinations where we have at least 3 studies:

ttestlab=tradeoff.labels[tradeoff.labels$numberstudies>2,c("tradeoff", "numberstudies")]
tradeoff.list=ttestlab$tradeoff
ttestpval=numeric(0)
for(i in 1:length(tradeoff.list)) {
  dat=prop.bystudy[prop.bystudy$tradeoff==as.character(tradeoff.list[i]),]
  ttest=t.test(dat$obsprop, dat$simulatedprop, paired=T)
  ttestpval=rbind(ttestpval, ttest$p.value)
}
ttestlab$ttestpval=ttestpval
ttestlab$numberstudies=NULL

#fig 4: average difference between observed and expected across studies with a treatment combination
forquilt=dcast(prop.bystudy, tradeoff~., value.var="dif.in.prop", fun=mean)
names(forquilt)[names(forquilt)=="."]<-"dif"
forquilt=merge(forquilt, ttestlab, by="tradeoff", all=T)
forquilt=merge(forquilt, tradeoff.labels, by="tradeoff")
forquilt$trtlabels1=matrix(unlist(strsplit(levels(forquilt$tradeoff2), " x ")), ncol=2, byrow=T)[,1]
forquilt$trtlabels2=matrix(unlist(strsplit(levels(forquilt$tradeoff2), " x ")), ncol=2, byrow=T)[,2]
forquilt$nlab=paste("n=", forquilt$numberstudies, sep="")
forquilt$plab=ifelse(is.na(forquilt$ttestpval), "", ifelse(forquilt$ttestpval>0.1, "NS", paste("p=", round(forquilt$ttestpval, 2), sep="")))
forquilt$pnlab=paste(forquilt$plab, " (", forquilt$numberstudies, ")", sep="")
forquilt$plab=ifelse(is.na(forquilt$ttestpval), "NA", ifelse(forquilt$ttestpval>0.1, "NS", paste("p=", round(forquilt$ttestpval, 2), sep="")))
forquilt$pnlab=paste(forquilt$plab, " (", forquilt$numberstudies, ")", sep="")
ggplot(aes(trtlabels1, trtlabels2), data=forquilt[!forquilt$tradeoff %in% c("irrxmult_nutrient", "mult_nutrientxtemp"),]) + geom_tile(aes(fill=dif)) + scale_fill_gradient2(low="firebrick4", mid="gray90", high="deepskyblue4", midpoint=0, na.value="white") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(fill="Observed - expected proportion") + geom_text(aes(label=pnlab), size=3)
ggsave("figs/fig 4.pdf", height=4, width=8)


# FIGURE S6: histogram of observed-expected proportions in each study separately

#getting p-values on our proportions 
nperms=length(perm.list)
output.realpropranks=dcast(output.E[output.E$permutation=="real" & output.E$variable %in% "tradeoffs.ranks",], site_project_comm+tradeoff~variable)
temp=output.proportions[output.proportions$quadrant=="prop.tradeoffs",]
output.realpropranks=merge(temp, output.realpropranks, by=c("site_project_comm", "tradeoff"))
output.realpropranks$p=ifelse(output.realpropranks$tradeoffs.ranks<=nperms/2, output.realpropranks$tradeoffs.ranks/nperms, (nperms-output.realpropranks$tradeoffs.ranks)/nperms)*2

#formatting data to plot
quadhist.bystudy=dcast(output.E[output.E$variable=="prop.tradeoffs",], site_project_comm+permutation+tradeoff~variable, fun=mean)
quadhist.bystudy$lab=as.factor(paste(quadhist.bystudy$tradeoff, quadhist.bystudy$site_project_comm, sep="::"))
quadhist.bystudy=merge(quadhist.bystudy, tradeoff.labels, by="tradeoff")
quadhist.bystudy.simmeans=dcast(quadhist.bystudy[!quadhist.bystudy$permutation=="real",], site_project_comm+tradeoff+tradeoff2+lab~., fun=mean, value.var="prop.tradeoffs")
names(quadhist.bystudy.simmeans)[5]="sim.mean"

#formatting p-values to put on figs
plabs=output.realpropranks[,c("site_project_comm", "tradeoff", "p")]
plabs=merge(plabs, quadhist.bystudy[quadhist.bystudy$permutation=="real",], by=c("site_project_comm", "tradeoff"))
plabs$lab=as.factor(paste(plabs$tradeoff, plabs$site_project_comm, sep="::"))
plabs$plab=paste("p=", plabs$p, sep="")

#all in one huge plot:
ggplot(data=quadhist.bystudy[!quadhist.bystudy$permutation=="real",], aes(prop.tradeoffs, fill=I("gray"))) + geom_histogram() + facet_wrap(~lab, ncol=6) + geom_vline(data=quadhist.bystudy[quadhist.bystudy$permutation=="real",], size=1, aes(xintercept=prop.tradeoffs)) + geom_vline(data=quadhist.bystudy.simmeans, size=1, lty=3, aes(xintercept=sim.mean)) + xlab("Proportion of species in QII and QIV") + ylab("Number of communities") + geom_text(data=plabs, aes(x=Inf, y=Inf, label=plab, hjust=1, vjust=1.5), size=3)
ggsave("figs/fig S6.pdf", width=12, height=15)



#---------COMPARING SIMULATED DATA TO REAL DATA-----OMITTING SPECIES ABSENT FROM EITHER ALL TREATMENT OR ALL CONTROL PLOTS


#looping through all tradeoffs we want
perm.list=unique(effect.sizes$permutation)
output.lm.E.d=numeric(0)
datallreal.E.d=numeric(0)
datallsim.E.d=numeric(0)

for (h in 1:length(ntradeoffs.all)) {
  trdf=as.data.frame(trtcombos.m[trtcombos.m$variable==as.character(ntradeoffs.all[h]), c("site_project_comm")])
  names(trdf)="site_project_comm"
  dath=merge(effect.sizes[effect.sizes$code=="present",], trdf, by="site_project_comm")
  names(dath)[names(dath)=="global.change.treatment"]="trt"
  names(dath)[names(dath)=="E"]="eff"
  trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff==as.character(ntradeoffs.all[h]),]
  dath=merge(dath, trdf.vars.m, by="trt")
  want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
  wantsim=want[!want$permutation=="real",]
  datallsim.E.d=rbind(datallsim.E.d, wantsim)
  wantreal=want[want$permutation=="real",]
  datallreal.E.d=rbind(datallreal.E.d, wantreal)
  nsites=unique(want$site_project_comm)
  temp.i=numeric(0)
  
  for(i in 1:length(perm.list)) {
    dati=want[want$permutation==as.character(perm.list[i]),]
    temp.j=numeric(0)
    
    for(j in 1:length(nsites)) {
      datj=dati[dati$site_project_comm==as.character(nsites[j]),]
      datj=datj[complete.cases(datj[,c("var1", "var2")]),]
      datj$QI=ifelse(datj$var1>0 & datj$var2>0, 1, 0)
      datj$QII=ifelse(datj$var1<0 & datj$var2>0, 1, 0)
      datj$QIII=ifelse(datj$var1<0 & datj$var2<0, 1, 0)
      datj$QIV=ifelse(datj$var1>0 & datj$var2<0, 1, 0)
      nspecies=length(datj$species)
      temp.j2=data.frame(site_project_comm=as.character(nsites[j]), prop.QI=sum(datj$QI)/nspecies, prop.QII=sum(datj$QII)/nspecies, prop.QIII=sum(datj$QIII)/nspecies, prop.QIV=sum(datj$QIV)/nspecies)
      temp.j=rbind(temp.j, temp.j2)
    }
    
    temp.i2=data.frame(temp.j, permutation=as.character(perm.list[i]))
    temp.i=rbind(temp.i, temp.i2)
    #collect all that for each permutation and label with the permutation
  }
  
  #loop through sites to rank observed proportion in quadrants against all permuted proportions
  temp.h=numeric(0)
  for(g in 1:length(nsites)) {
    dat=temp.i[temp.i$site_project_comm==as.character(nsites[g]),]
    dat$prop.correspondence=dat$prop.QI + dat$prop.QIII
    dat$prop.tradeoffs=dat$prop.QII + dat$prop.QIV
    dat$QI.ranks=rank(dat$prop.QI)
    dat$QII.ranks=rank(dat$prop.QII)
    dat$QIII.ranks=rank(dat$prop.QIII)
    dat$QIV.ranks=rank(dat$prop.QIV)	
    dat$correspondence.ranks=rank(dat$prop.correspondence)
    dat$tradeoffs.ranks=rank(dat$prop.tradeoffs)
    temp.h=rbind(temp.h, dat)
    
  }
  
  temp.h2=data.frame(temp.h, tradeoff=as.character(ntradeoffs.all[h]))
  output.lm.E.d=rbind(output.lm.E.d, temp.h2)
  #collect all that for each tradeoff and label with the tradeoff
}

output.E.d=melt(output.lm.E.d, id=c("site_project_comm", "permutation", "tradeoff"))



# FIGURE S3: species responses to treatment pairs

#first cleaning up long list of species responses to pairs of treatments for each tradeoff:
datallreal.E.d=datallreal.E.d[complete.cases(datallreal.E.d[,c("var1", "var2")]),]
datallreal.E.d=merge(datallreal.E.d, tradeoff.labels, by="tradeoff")

#adding in mean abundance in control plots
control.abund=dcast(effect.sizes[effect.sizes$permutation=="real",], site_project_comm+species~., value.var="relabund.control", fun=mean)
names(control.abund)[names(control.abund)=="."]<-"relabund.control"
datallreal.E.d=merge(datallreal.E.d, control.abund, by=c("site_project_comm", "species"))

#now cleaning up data for simulated communities:
datallsim.E.d=datallsim.E.d[complete.cases(datallsim.E.d[,c("var1", "var2")]),]
datallsim.E.d=merge(datallsim.E.d, tradeoff.labels, by="tradeoff")

#Figure S3: 
ggplot(data=datallreal.E.d, aes(var1, var2)) + geom_point(aes(size=relabund.control, alpha=I(0.2))) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + xlab("Effect of treatment 1 (E)") + ylab("Effect of treatment 2 (E)") + theme(legend.position="none") 
ggsave("figs/Fig S3a.pdf", height=6, width=6)

#calculating proportion of species in 4 quadrants (text in fig S3a):
datallreal.E.d$quad=as.factor(ifelse(datallreal.E.d$var1<0 & datallreal.E.d$var2<0, "QIII", ifelse(datallreal.E.d$var1>0 & datallreal.E.d$var2>0, "QI", ifelse(datallreal.E.d$var1>0 & datallreal.E.d$var2<0, "QIV", "QII"))))
nosp.quad.real=dcast(datallreal.E.d, quad~., value.var="permutation", fun=length) #number of species in each quadrant
datallsim.E.d$quad=as.factor(ifelse(datallsim.E.d$var1<0 & datallsim.E.d$var2<0, "QIII", ifelse(datallsim.E.d$var1>0 & datallsim.E.d$var2>0, "QI", ifelse(datallsim.E.d$var1>0 & datallsim.E.d$var2<0, "QIV", "QII"))))
nosp.quad.simt=dcast(datallsim.E.d, permutation+quad~., value.var="permutation", fun=length) #number of species in each quadrant
nosp.quad.sim=dcast(nosp.quad.simt, quad~., value.var=".", fun=mean)
nosp.quad=merge(nosp.quad.real, nosp.quad.sim, by="quad")
names(nosp.quad)=c("quad", "real", "mean.sim")
nosp.quad$real.lab=paste("O:", nosp.quad$real, sep="")
nosp.quad$sim.lab=paste("E:", round(nosp.quad$mean.sim, 0), sep="")
nosp.quad


#----------------COMPARING OBSERVED TO EXPECTED PROPORTIONS

#calculating proportion of species in 4 quadrants:
output.realprop=dcast(output.E.d[output.E.d$permutation=="real" & output.E.d$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs"),], site_project_comm+tradeoff+variable~.)
names(output.realprop)[names(output.realprop)=="."]<-"obsprop"
output.simprop=dcast(output.E.d[output.E.d$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs") & !output.E.d$permutation=="real",], site_project_comm+tradeoff+variable~., mean)
names(output.simprop)[names(output.simprop)=="."]<-"simulatedprop"
temp=dcast(output.simprop, tradeoff~variable, value.var="simulatedprop", fun=mean); temp 
output.proportions=merge(output.realprop, output.simprop, by=c("site_project_comm", "tradeoff", "variable"))
names(output.proportions)[names(output.proportions)=="variable"]<-"quadrant"
output.proportions$dif.in.prop=output.proportions$obsprop-output.proportions$simulatedprop; output.proportions


# FIGURE S4: histogram of observed-expected proportions

#is the proportion of species in QII and QIV (tradeoffs) greater than expected? 
prop.bystudy=output.proportions[output.proportions$quadrant=="prop.tradeoffs",]
ttest=t.test(prop.bystudy$obsprop, prop.bystudy$simulatedprop, paired=T); ttest #p<0.0001
hist(prop.bystudy$dif.in.prop) #looks pretty normal?
ttestlab=data.frame(p=ttest$p.value)

#plotting histogram of differences:
figS4a=ggplot(prop.bystudy, aes(x=dif.in.prop)) + geom_histogram(binwidth=0.04) + geom_vline(xintercept=mean(prop.bystudy$dif.in.prop)) + geom_vline(xintercept=0, lty=2) + xlab("Observed - Expected proportion\nin Quadrants II and IV") + ylab("Number of studies") + geom_text(data=ttestlab, aes(x=-Inf, y=Inf, label="p<0.001", hjust=-0.5, vjust=1.5), size=3); figS4a
ggsave("figs/fig S4a.pdf", width=4, height=3)

#omitting SERC (strongest tradeoff)
t.test(prop.bystudy$obsprop[!prop.bystudy$site_project_comm=="SERC_CXN_0"], prop.bystudy$simulatedprop[!prop.bystudy$site_project_comm=="SERC_CXN_0"], paired=T) #p<0.001


# FIGURE S5: quilt of quadrants

#---ttests for individual treatment combinations where we have at least 3 studies:

ttestlab=tradeoff.labels[tradeoff.labels$numberstudies>2,c("tradeoff", "numberstudies")]
tradeoff.list=ttestlab$tradeoff
ttestpval=numeric(0)
for(i in 1:length(tradeoff.list)) {
  dat=prop.bystudy[prop.bystudy$tradeoff==as.character(tradeoff.list[i]),]
  ttest=t.test(dat$obsprop, dat$simulatedprop, paired=T)
  ttestpval=rbind(ttestpval, ttest$p.value)
}
ttestlab$ttestpval=ttestpval
ttestlab$numberstudies=NULL

#fig S5: average difference between observed and expected across studies with a treatment combination
forquilt=dcast(prop.bystudy, tradeoff~., value.var="dif.in.prop", fun=mean)
names(forquilt)[names(forquilt)=="."]<-"dif"
forquilt=merge(forquilt, ttestlab, by="tradeoff", all=T)
forquilt=merge(forquilt, tradeoff.labels, by="tradeoff")
forquilt$trtlabels1=matrix(unlist(strsplit(levels(forquilt$tradeoff2), " x ")), ncol=2, byrow=T)[,1]
forquilt$trtlabels2=matrix(unlist(strsplit(levels(forquilt$tradeoff2), " x ")), ncol=2, byrow=T)[,2]
forquilt$nlab=paste("n=", forquilt$numberstudies, sep="")
forquilt$plab=ifelse(is.na(forquilt$ttestpval), "", ifelse(forquilt$ttestpval>0.1, "NS", paste("p=", round(forquilt$ttestpval, 2), sep="")))
forquilt$pnlab=paste(forquilt$plab, " (", forquilt$numberstudies, ")", sep="")
forquilt$plab=ifelse(is.na(forquilt$ttestpval), "NA", ifelse(forquilt$ttestpval>0.1, "NS", paste("p=", round(forquilt$ttestpval, 2), sep="")))
forquilt$pnlab=paste(forquilt$plab, " (", forquilt$numberstudies, ")", sep="")
ggplot(aes(trtlabels1, trtlabels2), data=forquilt[!forquilt$tradeoff %in% c("irrxmult_nutrient", "mult_nutrientxtemp"),]) + geom_tile(aes(fill=dif)) + scale_fill_gradient2(low="firebrick4", mid="gray90", high="deepskyblue4", midpoint=0, na.value="white") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(fill="Observed - expected proportion") + geom_text(aes(label=pnlab), size=3)
ggsave("figs/fig S5.pdf", height=4, width=8)


# FIGURE S7: histogram of observed-expected proportions in each study separately

#getting p-values on our proportions 
nperms=length(perm.list)
output.realpropranks=dcast(output.E.d[output.E.d$permutation=="real" & output.E.d$variable %in% "tradeoffs.ranks",], site_project_comm+tradeoff~variable)
temp=output.proportions[output.proportions$quadrant=="prop.tradeoffs",]
output.realpropranks=merge(temp, output.realpropranks, by=c("site_project_comm", "tradeoff"))
output.realpropranks$p=ifelse(output.realpropranks$tradeoffs.ranks<=nperms/2, output.realpropranks$tradeoffs.ranks/nperms, (nperms-output.realpropranks$tradeoffs.ranks)/nperms)*2

#formatting data to plot
quadhist.bystudy=dcast(output.E.d[output.E.d$variable=="prop.tradeoffs",], site_project_comm+permutation+tradeoff~variable, fun=mean)
quadhist.bystudy$lab=as.factor(paste(quadhist.bystudy$tradeoff, quadhist.bystudy$site_project_comm, sep="::"))
quadhist.bystudy=merge(quadhist.bystudy, tradeoff.labels, by="tradeoff")
quadhist.bystudy.simmeans=dcast(quadhist.bystudy[!quadhist.bystudy$permutation=="real",], site_project_comm+tradeoff+tradeoff2+lab~., fun=mean, value.var="prop.tradeoffs")
names(quadhist.bystudy.simmeans)[5]="sim.mean"

#formatting p-values to put on figs
plabs=output.realpropranks[,c("site_project_comm", "tradeoff", "p")]
plabs=merge(plabs, quadhist.bystudy[quadhist.bystudy$permutation=="real",], by=c("site_project_comm", "tradeoff"))
plabs$lab=as.factor(paste(plabs$tradeoff, plabs$site_project_comm, sep="::"))
plabs$plab=paste("p=", plabs$p, sep="")

#all in one huge plot:
ggplot(data=quadhist.bystudy[!quadhist.bystudy$permutation=="real",], aes(prop.tradeoffs, fill=I("gray"))) + geom_histogram() + facet_wrap(~lab, ncol=6) + geom_vline(data=quadhist.bystudy[quadhist.bystudy$permutation=="real",], size=1, aes(xintercept=prop.tradeoffs)) + geom_vline(data=quadhist.bystudy.simmeans, size=1, lty=3, aes(xintercept=sim.mean)) + xlab("Proportion of species in QII and QIV") + ylab("Number of communities") + geom_text(data=plabs, aes(x=Inf, y=Inf, label=plab, hjust=1, vjust=1.5), size=3)
ggsave("figs/fig S7.pdf", width=12, height=15)




#---------COMPARING SIMULATED DATA TO REAL DATA-----LRR (OMITS SPECIES ABSENT FROM EITHER ALL TREATMENT OR ALL CONTROL PLOTS)


#looping through all tradeoffs we want
perm.list=unique(effect.sizes$permutation)
output.lm.LRR=numeric(0)
datallreal.LRR=numeric(0)
datallsim.LRR=numeric(0)

for (h in 1:length(ntradeoffs.all)) {
  trdf=as.data.frame(trtcombos.m[trtcombos.m$variable==as.character(ntradeoffs.all[h]), c("site_project_comm")])
  names(trdf)="site_project_comm"
  dath=merge(effect.sizes, trdf, by="site_project_comm")
  names(dath)[names(dath)=="global.change.treatment"]="trt"
  names(dath)[names(dath)=="LRR"]="eff"
  trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff==as.character(ntradeoffs.all[h]),]
  dath=merge(dath, trdf.vars.m, by="trt")
  want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
  wantsim=want[!want$permutation=="real",]
  datallsim.LRR=rbind(datallsim.LRR, wantsim)
  wantreal=want[want$permutation=="real",]
  datallreal.LRR=rbind(datallreal.LRR, wantreal)
  nsites=unique(want$site_project_comm)
  temp.i=numeric(0)
  
  for(i in 1:length(perm.list)) {
    dati=want[want$permutation==as.character(perm.list[i]),]
    temp.j=numeric(0)
    
    for(j in 1:length(nsites)) {
      datj=dati[dati$site_project_comm==as.character(nsites[j]),]
      datj=datj[complete.cases(datj[,c("var1", "var2")]),]
      datj$QI=ifelse(datj$var1>0 & datj$var2>0, 1, 0)
      datj$QII=ifelse(datj$var1<0 & datj$var2>0, 1, 0)
      datj$QIII=ifelse(datj$var1<0 & datj$var2<0, 1, 0)
      datj$QIV=ifelse(datj$var1>0 & datj$var2<0, 1, 0)
      nspecies=length(datj$species)
      temp.j2=data.frame(site_project_comm=as.character(nsites[j]), prop.QI=sum(datj$QI)/nspecies, prop.QII=sum(datj$QII)/nspecies, prop.QIII=sum(datj$QIII)/nspecies, prop.QIV=sum(datj$QIV)/nspecies)
      temp.j=rbind(temp.j, temp.j2)
    }
    
    temp.i2=data.frame(temp.j, permutation=as.character(perm.list[i]))
    temp.i=rbind(temp.i, temp.i2)
    #collect all that for each permutation and label with the permutation
  }
  
  #loop through sites to rank observed proportion in quadrants against all permuted proportions
  temp.h=numeric(0)
  for(g in 1:length(nsites)) {
    dat=temp.i[temp.i$site_project_comm==as.character(nsites[g]),]
    dat$prop.correspondence=dat$prop.QI + dat$prop.QIII
    dat$prop.tradeoffs=dat$prop.QII + dat$prop.QIV
    dat$QI.ranks=rank(dat$prop.QI)
    dat$QII.ranks=rank(dat$prop.QII)
    dat$QIII.ranks=rank(dat$prop.QIII)
    dat$QIV.ranks=rank(dat$prop.QIV)	
    dat$correspondence.ranks=rank(dat$prop.correspondence)
    dat$tradeoffs.ranks=rank(dat$prop.tradeoffs)
    temp.h=rbind(temp.h, dat)
    
  }
  
  temp.h2=data.frame(temp.h, tradeoff=as.character(ntradeoffs.all[h]))
  output.lm.LRR=rbind(output.lm.LRR, temp.h2)
  #collect all that for each tradeoff and label with the tradeoff
}

output.LRR=melt(output.lm.LRR, id=c("site_project_comm", "permutation", "tradeoff"))



# FIGURE S3: species responses to treatment pairs

#first cleaning up long list of species responses to pairs of treatments for each tradeoff:
datallreal.LRR=datallreal.LRR[complete.cases(datallreal.LRR[,c("var1", "var2")]),]
datallreal.LRR=merge(datallreal.LRR, tradeoff.labels, by="tradeoff")

#adding in mean abundance in control plots
control.abund=dcast(effect.sizes[effect.sizes$permutation=="real",], site_project_comm+species~., value.var="relabund.control", fun=mean)
names(control.abund)[names(control.abund)=="."]<-"relabund.control"
datallreal.LRR=merge(datallreal.LRR, control.abund, by=c("site_project_comm", "species"))

#now cleaning up data for simulated communities:
datallsim.LRR=datallsim.LRR[complete.cases(datallsim.LRR[,c("var1", "var2")]),]
datallsim.LRR=merge(datallsim.LRR, tradeoff.labels, by="tradeoff")

#Figure S3: 
ggplot(data=datallreal.LRR, aes(var1, var2)) + geom_point(aes(size=relabund.control, alpha=I(0.2))) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + xlab("Effect of treatment 1 (E)") + ylab("Effect of treatment 2 (E)") + theme(legend.position="none") 
ggsave("figs/Fig S3b.pdf", height=6, width=6)

#calculating proportion of species in 4 quadrants (text in fig S3b):
datallreal.LRR$quad=as.factor(ifelse(datallreal.LRR$var1<0 & datallreal.LRR$var2<0, "QIII", ifelse(datallreal.LRR$var1>0 & datallreal.LRR$var2>0, "QI", ifelse(datallreal.LRR$var1>0 & datallreal.LRR$var2<0, "QIV", "QII"))))
nosp.quad.real=dcast(datallreal.LRR, quad~., value.var="permutation", fun=length) #number of species in each quadrant
datallsim.LRR$quad=as.factor(ifelse(datallsim.LRR$var1<0 & datallsim.LRR$var2<0, "QIII", ifelse(datallsim.LRR$var1>0 & datallsim.LRR$var2>0, "QI", ifelse(datallsim.LRR$var1>0 & datallsim.LRR$var2<0, "QIV", "QII"))))
nosp.quad.simt=dcast(datallsim.LRR, permutation+quad~., value.var="permutation", fun=length) #number of species in each quadrant
nosp.quad.sim=dcast(nosp.quad.simt, quad~., value.var=".", fun=mean)
nosp.quad=merge(nosp.quad.real, nosp.quad.sim, by="quad")
names(nosp.quad)=c("quad", "real", "mean.sim")
nosp.quad$real.lab=paste("O:", nosp.quad$real, sep="")
nosp.quad$sim.lab=paste("E:", round(nosp.quad$mean.sim, 0), sep="")
nosp.quad


#----------------COMPARING OBSERVED TO EXPECTED PROPORTIONS

#calculating proportion of species in 4 quadrants:
output.realprop=dcast(output.LRR[output.LRR$permutation=="real" & output.LRR$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs"),], site_project_comm+tradeoff+variable~.)
names(output.realprop)[names(output.realprop)=="."]<-"obsprop"
output.simprop=dcast(output.LRR[output.LRR$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs") & !output.LRR$permutation=="real",], site_project_comm+tradeoff+variable~., mean)
names(output.simprop)[names(output.simprop)=="."]<-"simulatedprop"
temp=dcast(output.simprop, tradeoff~variable, value.var="simulatedprop", fun=mean); temp 
output.proportions=merge(output.realprop, output.simprop, by=c("site_project_comm", "tradeoff", "variable"))
names(output.proportions)[names(output.proportions)=="variable"]<-"quadrant"
output.proportions$dif.in.prop=output.proportions$obsprop-output.proportions$simulatedprop; output.proportions


# FIGURE S4: histogram of observed-expected proportions. note that this is identical to the fig created by E.d

#is the proportion of species in QII and QIV (tradeoffs) greater than expected? 
prop.bystudy=output.proportions[output.proportions$quadrant=="prop.tradeoffs",]
ttest=t.test(prop.bystudy$obsprop, prop.bystudy$simulatedprop, paired=T); ttest #p<0.0001
hist(prop.bystudy$dif.in.prop) #looks pretty normal?
ttestlab=data.frame(p=ttest$p.value)

#plotting histogram of differences:
figS4a=ggplot(prop.bystudy, aes(x=dif.in.prop)) + geom_histogram(binwidth=0.04) + geom_vline(xintercept=mean(prop.bystudy$dif.in.prop)) + geom_vline(xintercept=0, lty=2) + xlab("Observed - Expected proportion\nin Quadrants II and IV") + ylab("Number of studies") + geom_text(data=ttestlab, aes(x=-Inf, y=Inf, label="p<0.001", hjust=-0.5, vjust=1.5), size=3); figS4a

#omitting SERC (strongest tradeoff)
t.test(prop.bystudy$obsprop[!prop.bystudy$site_project_comm=="SERC_CXN_0"], prop.bystudy$simulatedprop[!prop.bystudy$site_project_comm=="SERC_CXN_0"], paired=T) #p<0.001


# FIGURE S5: quilt of quadrants

#---ttests for individual treatment combinations where we have at least 3 studies:

ttestlab=tradeoff.labels[tradeoff.labels$numberstudies>2,c("tradeoff", "numberstudies")]
tradeoff.list=ttestlab$tradeoff
ttestpval=numeric(0)
for(i in 1:length(tradeoff.list)) {
  dat=prop.bystudy[prop.bystudy$tradeoff==as.character(tradeoff.list[i]),]
  ttest=t.test(dat$obsprop, dat$simulatedprop, paired=T)
  ttestpval=rbind(ttestpval, ttest$p.value)
}
ttestlab$ttestpval=ttestpval
ttestlab$numberstudies=NULL

#fig S5: average difference between observed and expected across studies with a treatment combination
forquilt=dcast(prop.bystudy, tradeoff~., value.var="dif.in.prop", fun=mean)
names(forquilt)[names(forquilt)=="."]<-"dif"
forquilt=merge(forquilt, ttestlab, by="tradeoff", all=T)
forquilt=merge(forquilt, tradeoff.labels, by="tradeoff")
forquilt$trtlabels1=matrix(unlist(strsplit(levels(forquilt$tradeoff2), " x ")), ncol=2, byrow=T)[,1]
forquilt$trtlabels2=matrix(unlist(strsplit(levels(forquilt$tradeoff2), " x ")), ncol=2, byrow=T)[,2]
forquilt$nlab=paste("n=", forquilt$numberstudies, sep="")
forquilt$plab=ifelse(is.na(forquilt$ttestpval), "", ifelse(forquilt$ttestpval>0.1, "NS", paste("p=", round(forquilt$ttestpval, 2), sep="")))
forquilt$pnlab=paste(forquilt$plab, " (", forquilt$numberstudies, ")", sep="")
forquilt$plab=ifelse(is.na(forquilt$ttestpval), "NA", ifelse(forquilt$ttestpval>0.1, "NS", paste("p=", round(forquilt$ttestpval, 2), sep="")))
forquilt$pnlab=paste(forquilt$plab, " (", forquilt$numberstudies, ")", sep="")
ggplot(aes(trtlabels1, trtlabels2), data=forquilt[!forquilt$tradeoff %in% c("irrxmult_nutrient", "mult_nutrientxtemp"),]) + geom_tile(aes(fill=dif)) + scale_fill_gradient2(low="firebrick4", mid="gray90", high="deepskyblue4", midpoint=0, na.value="white") + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + labs(fill="Observed - expected proportion") + geom_text(aes(label=pnlab), size=3)


# FIGURE S7: histogram of observed-expected proportions in each study separately

#getting p-values on our proportions 
nperms=length(perm.list)
output.realpropranks=dcast(output.LRR[output.LRR$permutation=="real" & output.LRR$variable %in% "tradeoffs.ranks",], site_project_comm+tradeoff~variable)
temp=output.proportions[output.proportions$quadrant=="prop.tradeoffs",]
output.realpropranks=merge(temp, output.realpropranks, by=c("site_project_comm", "tradeoff"))
output.realpropranks$p=ifelse(output.realpropranks$tradeoffs.ranks<=nperms/2, output.realpropranks$tradeoffs.ranks/nperms, (nperms-output.realpropranks$tradeoffs.ranks)/nperms)*2

#formatting data to plot
quadhist.bystudy=dcast(output.LRR[output.LRR$variable=="prop.tradeoffs",], site_project_comm+permutation+tradeoff~variable, fun=mean)
quadhist.bystudy$lab=as.factor(paste(quadhist.bystudy$tradeoff, quadhist.bystudy$site_project_comm, sep="::"))
quadhist.bystudy=merge(quadhist.bystudy, tradeoff.labels, by="tradeoff")
quadhist.bystudy.simmeans=dcast(quadhist.bystudy[!quadhist.bystudy$permutation=="real",], site_project_comm+tradeoff+tradeoff2+lab~., fun=mean, value.var="prop.tradeoffs")
names(quadhist.bystudy.simmeans)[5]="sim.mean"

#formatting p-values to put on figs
plabs=output.realpropranks[,c("site_project_comm", "tradeoff", "p")]
plabs=merge(plabs, quadhist.bystudy[quadhist.bystudy$permutation=="real",], by=c("site_project_comm", "tradeoff"))
plabs$lab=as.factor(paste(plabs$tradeoff, plabs$site_project_comm, sep="::"))
plabs$plab=paste("p=", plabs$p, sep="")

#all in one huge plot:
ggplot(data=quadhist.bystudy[!quadhist.bystudy$permutation=="real",], aes(prop.tradeoffs, fill=I("gray"))) + geom_histogram() + facet_wrap(~lab, ncol=6) + geom_vline(data=quadhist.bystudy[quadhist.bystudy$permutation=="real",], size=1, aes(xintercept=prop.tradeoffs)) + geom_vline(data=quadhist.bystudy.simmeans, size=1, lty=3, aes(xintercept=sim.mean)) + xlab("Proportion of species in QII and QIV") + ylab("Number of communities") + geom_text(data=plabs, aes(x=Inf, y=Inf, label=plab, hjust=1, vjust=1.5), size=3)






# COMPARING E TO E.d, LRR, AND PERCENT STIMULATION


# subsetting to real data only:

real.effect.sizes=effect.sizes[effect.sizes$permutation=="real",]

#getting table of number of species dropped with various methods:
res.m=melt(real.effect.sizes, id=c("site_project_comm", "species", "global.change.treatment"), measure=c("E", "LRR", "PS"))
res.m=res.m[!is.na(res.m$value),]
sample.sizes=dcast(res.m, site_project_comm~variable, fun=length)
res.m2=melt(real.effect.sizes, id=c("site_project_comm", "species", "global.change.treatment"), measure="code")
sample.sizes2=dcast(res.m2, site_project_comm~value, fun=length)
sample.sizes=merge(sample.sizes, sample.sizes2, by="site_project_comm"); sample.sizes


#plotting comparisons:

EvsLRR=qplot(E, LRR, data=real.effect.sizes, alpha=I(0.1))
cor.test(real.effect.sizes$E, real.effect.sizes$LRR, method="spearman")

EvsPS=qplot(E, PS, data=real.effect.sizes, ylab="Percent stimulation", alpha=I(0.1), ylim=c(-10, 200))
cor.test(real.effect.sizes$E, real.effect.sizes$PS, method="spearman")

ggarrange(EvsLRR, EvsPS, ncol=1)
ggsave("figs/fig S2.pdf", width=3, height=5)

#results with E.d are identical, because E.d drops plots that must already be dropped from LRR and PS
E.dvsLRR=qplot(E, LRR, data=real.effect.sizes[real.effect.sizes$code=="present",], xlab="E.d", alpha=I(0.1))
cor.test(real.effect.sizes[real.effect.sizes$code=="present",]$E, real.effect.sizes[real.effect.sizes $code=="present",]$LRR, method="spearman")
E.dvsPS=qplot(E, PS, data=real.effect.sizes[real.effect.sizes$code=="present",], xlab="E.d", ylab="Percent stimulation", alpha=I(0.1), ylim=c(-10, 200))
cor.test(real.effect.sizes[real.effect.sizes$code=="present",]$E, real.effect.sizes[real.effect.sizes$code=="present",]$PS, method="spearman")
ggarrange(E.dvsLRR, E.dvsPS, ncol=1)


# DEMONSTRATING IMPORTANCE OF NULL MODEL: CORRELATION COEFFICIENTS FOR LRR AND E

#FIG S8

Esim1=ggplot(data=datallsim.E[datallsim.E$permutation==1,], aes(var1, var2)) + geom_point(aes(alpha=I(0.2))) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + geom_smooth(method="lm", color="black") + xlab("Effect of treatment 1 (E)") + ylab("Effect of treatment 2 (E)") + theme(legend.position="none") + facet_wrap(~tradeoff2, ncol=3) + ggtitle("a) Effect size E")
LRRsim1=ggplot(data=datallsim.LRR[datallsim.LRR$permutation==1,], aes(var1, var2)) + geom_point(aes(alpha=I(0.2))) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + geom_smooth(method="lm", color="black") + xlab("Effect of treatment 1 (LRR)") + ylab("Effect of treatment 2 (LRR)") + theme(legend.position="none") + facet_wrap(~tradeoff2, ncol=3) + ggtitle("b) Log response ration (LRR)")
ggarrange(Esim1, LRRsim1, ncol=2)
ggsave("figs/fig S8.pdf", height=10, width=14)

#TABLE S3: correlation coefficients for all 999 communities

tableS3=numeric(0)

temp=cor(datallsim.E$var1, datallsim.E$var2); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="BMCTxdrought"], datallsim.E$var2[datallsim.E$tradeoff=="BMCTxdrought"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="BMCTxirr"], datallsim.E$var2[datallsim.E$tradeoff=="BMCTxirr"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="BMCTxN"], datallsim.E$var2[datallsim.E$tradeoff=="BMCTxN"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="BMCTxP"], datallsim.E$var2[datallsim.E$tradeoff=="BMCTxP"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="BMCTxtemp"], datallsim.E$var2[datallsim.E$tradeoff=="BMCTxtemp"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="CO2xirr"], datallsim.E$var2[datallsim.E$tradeoff=="CO2xirr"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="CO2xN"], datallsim.E$var2[datallsim.E$tradeoff=="CO2xN"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="CO2xtemp"], datallsim.E$var2[datallsim.E$tradeoff=="CO2xtemp"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="droughtxirr"], datallsim.E$var2[datallsim.E$tradeoff=="droughtxirr"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="droughtxN"], datallsim.E$var2[datallsim.E$tradeoff=="droughtxN"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="droughtxtemp"], datallsim.E$var2[datallsim.E$tradeoff=="droughtxtemp"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="irrxmult_nutrient"], datallsim.E$var2[datallsim.E$tradeoff=="irrxmult_nutrient"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="irrxN"], datallsim.E$var2[datallsim.E$tradeoff=="irrxN"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="irrxP"], datallsim.E$var2[datallsim.E$tradeoff=="irrxP"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="irrxtemp"], datallsim.E$var2[datallsim.E$tradeoff=="irrxtemp"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="mult_nutrientxtemp"], datallsim.E$var2[datallsim.E$tradeoff=="mult_nutrientxtemp"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="NxP"], datallsim.E$var2[datallsim.E$tradeoff=="NxP"]); tableS3=rbind(tableS3, temp)
temp=cor(datallsim.E$var1[datallsim.E$tradeoff=="Nxtemp"], datallsim.E$var2[datallsim.E$tradeoff=="Nxtemp"]); tableS3=rbind(tableS3, temp)

temp2=numeric(0)
temp=cor(datallsim.E.d$var1, datallsim.E.d$var2); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="BMCTxdrought"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="BMCTxdrought"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="BMCTxirr"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="BMCTxirr"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="BMCTxN"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="BMCTxN"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="BMCTxP"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="BMCTxP"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="BMCTxtemp"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="BMCTxtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="CO2xirr"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="CO2xirr"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="CO2xN"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="CO2xN"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="CO2xtemp"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="CO2xtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="droughtxirr"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="droughtxirr"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="droughtxN"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="droughtxN"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="droughtxtemp"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="droughtxtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="irrxmult_nutrient"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="irrxmult_nutrient"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="irrxN"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="irrxN"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="irrxP"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="irrxP"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="irrxtemp"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="irrxtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="mult_nutrientxtemp"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="mult_nutrientxtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="NxP"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="NxP"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.E.d$var1[datallsim.E.d$tradeoff=="Nxtemp"], datallsim.E.d$var2[datallsim.E.d$tradeoff=="Nxtemp"]); temp2=rbind(temp2, temp)

tableS3=cbind(tableS3, temp2)


temp2=numeric(0)
temp=cor(datallsim.LRR$var1, datallsim.LRR$var2); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="BMCTxdrought"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="BMCTxdrought"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="BMCTxirr"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="BMCTxirr"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="BMCTxN"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="BMCTxN"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="BMCTxP"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="BMCTxP"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="BMCTxtemp"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="BMCTxtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="CO2xirr"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="CO2xirr"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="CO2xN"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="CO2xN"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="CO2xtemp"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="CO2xtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="droughtxirr"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="droughtxirr"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="droughtxN"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="droughtxN"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="droughtxtemp"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="droughtxtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="irrxmult_nutrient"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="irrxmult_nutrient"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="irrxN"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="irrxN"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="irrxP"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="irrxP"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="irrxtemp"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="irrxtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="mult_nutrientxtemp"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="mult_nutrientxtemp"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="NxP"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="NxP"]); temp2=rbind(temp2, temp)
temp=cor(datallsim.LRR$var1[datallsim.LRR$tradeoff=="Nxtemp"], datallsim.LRR$var2[datallsim.LRR$tradeoff=="Nxtemp"]); temp2=rbind(temp2, temp)

tableS3=cbind(tableS3, temp2)
write.csv(tableS3, "figs/table S3.csv", row.names=F)
