library(plyr)
library(sjPlot)

## Setup colours
col.XX_M <- rgb(245/255, 140/255, 2/255, 3/4)
col.XX_F <- rgb(249/255, 187/255, 235/255, 3/4)
col.YA1a <- rgb(219/255, 35/255, 33/255, 3/4)
col.YB10 <- rgb(50/255, 163/255, 10/255, 3/4)
col.YB1a <- rgb(35/255, 101/255, 164/255, 3/4)
col.YB20 <- rgb(158/255, 206/255, 141/255, 3/4)
col.YB2a <- rgb(62/255, 172/255, 65/255, 3/4)
col.YB2b <- rgb(62/255, 172/255, 65/255, 3/4) # same as a
col.YB30 <- rgb(199/255, 216/255, 22/255, 3/4) # same as 0
col.YB40 <- rgb(199/255, 216/255, 22/255, 3/4) # same as 0
col.YB50 <- rgb(199/255, 216/255, 22/255, 3/4) # same as 0

cbblue <- rgb(0/255, 114/255, 178/255)
cbred <- rgb(219/255, 35/255, 33/255)
cbpink <- rgb(240/255, 150/255, 150/255)
cbgreen <- rgb(0/255, 158/255, 115/255)
cborange <- rgb(230/255, 159/255, 0/255)
cbpurple <- rgb(204/255, 121/255, 167/255)
cbbluegreen <- rgb(0/255, 125/255, 125/255)

## test the colours
par(mfrow=c(2,1)) 
plot(seq(1,10), rep(1,10), pch=16, cex=3, col=c(col.XX_M,col.YA1a,col.YB10,col.YB1a,col.YB20,col.YB2a,col.YB2b,col.YB30,col.YB40,col.YB50))
plot(seq(1,7), rep(1,7), pch=16, cex=3, col=c(cbblue, cbred, cbpurple, cborange, cbgreen, cbpink, cbbluegreen))

# Analysis of phenotypic data from Meitreile
datapath <- '~/git/amplexus_success/input'
outpath <- '~/git/amplexus_success/output'
frogdata <- read.table(file.path(outpath, 'pheno_data.txt'), header=T)

str(frogdata)
frogdata$year <- as.factor(frogdata$year)

### All data summarised by genotype
summary(frogdata$genotype)


# frogdata <- subset(frogdata_All, frogdata_All$year==2015) 

sjt.xtab(frogdata$genotype, frogdata$year, title='Continency table of genotype by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genotype_year.html'))

### Define by recombination between microsatellites and Dmrt region
frogdata$genodiff <- as.character(frogdata$genotype)
frogdata$genodiff[frogdata$genotype=="YB1a" | frogdata$genotype=="YB2a" | frogdata$genotype=="YB2b" | frogdata$genotype=="YA1a"] <- "diff"
frogdata$genodiff[frogdata$genotype=="YB10" | frogdata$genotype=="YB20" | frogdata$genotype=="YB30" | frogdata$genotype=="YB40" | frogdata$genotype=="YB50"] <- "semidiff"
frogdata$genodiff[frogdata$genotype=="XX"] <- "undiff"
frogdata$genodiff <- as.factor(frogdata$genodiff)
summary(frogdata$genodiff)

### Define by major dmrt regions
frogdata$genodmrt <- as.character(frogdata$genotype)
frogdata$genodmrt[frogdata$genotype=="YB1a" | frogdata$genotype=="YB10"] <- "YB1"
frogdata$genodmrt[frogdata$genotype=="YB2a" | frogdata$genotype=="YB2b" | frogdata$genotype=="YB20"] <- "YB2"
frogdata$genodmrt[frogdata$genotype=="YB30" | frogdata$genotype=="YB40" | frogdata$genotype=="YB50"] <- "YBother"
frogdata$genodmrt <- as.factor(frogdata$genodmrt)
summary(frogdata$genodmrt)

# Plot summarising types of male
par(mfrow=c(1,2)) 
par(las=2)
par(mar=c(5,5,4,3))
barplot(table(frogdata$genodiff), horiz=T, main='Males by Y differentiation', cex.main=1.8, cex.axis=1.5, cex.names=1.35,col=c(cbbluegreen, cbpurple, cbred))
barplot(table(frogdata$genodmrt), horiz=T, main='Males by dmrt region', cex.main=1.8, cex.axis=1.5, cex.names=1.35, col=c(col.XX_M, col.YB1a, col.YB2a, col.YB40))
dev.copy(pdf,file.path(outpath, 'pop_summary.pdf'), width=16, height=8)
dev.off()


## data without YA1a
dmrtPlot <- subset(frogdata, frogdata$genodmrt!='YA1a') 
dmrtPlot$genodmrt <- factor(dmrtPlot$genodmrt)
summary(dmrtPlot$genodmrt)

sjt.xtab(dmrtPlot$genodiff, dmrtPlot$year, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodiff_year.html'))
sjt.xtab(dmrtPlot$genodmrt, dmrtPlot$year, title='Continency table of dmrt regions by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodmrt_year.html'))

## data with A N amplexus only (not with dead or males)
normalPlexusData <- subset(dmrtPlot, dmrtPlot$amplexus=='A' | dmrtPlot$amplexus=='N') 
normalPlexusData$amplexus <- factor(normalPlexusData$amplexus)
sjt.xtab(normalPlexusData$genodiff, normalPlexusData$amplexus, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodiff_amplexus.html'))
sjt.xtab(normalPlexusData$genodmrt, normalPlexusData$amplexus, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodmrt_amplexus.html'))

### Data without X males
YPlexusData <- subset(normalPlexusData, normalPlexusData$genotype!='XX')
YPlexusData$genotype <- factor(YPlexusData$genotype)
YPlexusData$genodiff <- factor(YPlexusData$genodiff)
YPlexusData$genodmrt <- factor(YPlexusData$genodmrt)

sjt.xtab(YPlexusData$genodiff, YPlexusData$amplexus, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodiff_amplexus_Yonly.html'))
sjt.xtab(YPlexusData$genodmrt, YPlexusData$amplexus, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodmrt_amplexus_Yonly.html'))

## data where amplexus with dead and male male interactions are included
dmrtPlot$amplexusAll <- dmrtPlot$amplexus
dmrtPlot$amplexusAll[dmrtPlot$amplexusAll=='AG'] <- 'A'
dmrtPlot$amplexusAll[dmrtPlot$amplexusAll=='D'] <- 'A'
dmrtPlot$amplexusAll <- factor(dmrtPlot$amplexusAll)
sjt.xtab(dmrtPlot$genodiff, dmrtPlot$amplexusAll, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodiff_amplexusAll.html'))
sjt.xtab(dmrtPlot$genodmrt, dmrtPlot$amplexusAll, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodmrt_amplexusAll.html'))

## Just extreme best and worst genotypes
XXXYB1data <- subset(dmrtPlot, dmrtPlot$genodmrt=='XX' | dmrtPlot$genodmrt=='YB1')
XXXYB1data$genodmrt <- factor(XXXYB1data$genodmrt)
sjt.xtab(XXXYB1data$genodiff, XXXYB1data$amplexusAll, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodiff_amplexusXXXYB1.html'))
sjt.xtab(XXXYB1data$genodmrt, XXXYB1data$amplexusAll, title='Continency table of Y differentiation by year', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodmrt_amplexusXXXYB1.html'))

sjt.xtab(dmrtPlot$genodmrt, dmrtPlot$genodiff, title='Continency table of genodmrt by genodiff', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genodmrt_genodiff.html'))
# ?sjt.xtab


## GLM

# normalPlexusData
amplexus_PA <- cbind(normalPlexusData$amplexus=='A',normalPlexusData$amplexus!='A')

AN_genomod <- glm(amplexus_PA ~ genodiff, family='binomial', data=normalPlexusData)
summary(AN_genomod)
anova(AN_genomod, test='Chisq') # Y differentiation does not influence amplexus

AN_dmrtmod <- glm(amplexus_PA ~ genodmrt, family='binomial', data=normalPlexusData)
summary(AN_dmrtmod)
anova(AN_dmrtmod, test='Chisq') # Dmrt alleles do not influence amplexus

AN_typemod <- glm(amplexus_PA ~ genotype, family='binomial', data=normalPlexusData)
summary(AN_typemod)
anova(AN_typemod, test='Chisq') # YB2a stands out with p=0.0857, or 0.0656 if all amplexus (including dead) are used

tab_model(AN_genomod, AN_dmrtmod, AN_typemod)

# Phenotypic analysis
pheno2014 <- subset(normalPlexusData, normalPlexusData$year=='2014')
# pheno2014 <- subset(normalPlexusData, normalPlexusData$year=='2014' & normalPlexusData$amplexus!='A')
phenoNot2014 <- subset(normalPlexusData, normalPlexusData$year!='2014')
phenoData <- rbind(pheno2014, phenoNot2014)
phenoData$genodiff <- relevel(phenoData$genodiff, ref=3) # level 1 is now undiff

boxplot(pheno2014$w[pheno2014$day=='A'], pheno2014$w[pheno2014$day=='B'], pheno2014$w[pheno2014$day=='C'], pheno2014$w[pheno2014$day=='D'], pheno2014$w[pheno2014$day=='E'], pheno2014$w[pheno2014$day=='F'], pheno2014$w[pheno2014$day=='G'], pheno2014$w[pheno2014$day=='H'], pheno2014$w[pheno2014$day=='I'], pheno2014$w[pheno2014$day=='J'], main="2014 weight by collection day (last 3 are amplexus)", ylab="weight", xlab="day")

AN_data <- subset(normalPlexusData, normalPlexusData$year=='2015' & normalPlexusData$l1!='NA')
amplexus_PA <- cbind(AN_data$amplexus=='A', AN_data$amplexus!='A')
AN_data$wl1 <- AN_data$w/AN_data$l1
AN_data$wl2 <- AN_data$w/AN_data$l2
AN_data$l1l2 <- AN_data$l1/AN_data$l2

AN_w_mod <- glm(amplexus_PA ~ w, family='binomial', data=AN_data)
summary(AN_w_mod)
AN_l1_mod <- glm(amplexus_PA ~ l1, family='binomial', data=AN_data)
summary(AN_l1_mod)
AN_l2_mod <- glm(amplexus_PA ~ l2, family='binomial', data=AN_data)
summary(AN_l2_mod)
AN_wl1_mod <- glm(amplexus_PA ~ wl1, family='binomial', data=AN_data)
summary(AN_wl1_mod)
AN_wl2_mod <- glm(amplexus_PA ~ wl2, family='binomial', data=AN_data)
summary(AN_wl2_mod)
AN_l1l2_mod <- glm(amplexus_PA ~ l1l2, family='binomial', data=AN_data)
summary(AN_l1l2_mod)
# str(anova(AN_wl2_mod, test='Chisq')$'Pr(>Chi)')

tapply(AN_data$wl1, list(AN_data$genotype, AN_data$year), mean, na.rm=T) #summarise cloud mean by 2 factors,
 
tab_model(AN_w_mod, AN_l1_mod, AN_l2_mod, AN_wl1_mod, AN_wl2_mod, AN_l1l2_mod, emph.p=T, show.aic=T, file=file.path(outpath, 'AN_pheno.html'))


diff_mod1 <- lm(w ~ genodiff + year, data=phenoData)
anova(diff_mod1)
summary(diff_mod1)

diff_mod2 <- lm(l1 ~ genodiff + year, data=phenoData)
anova(diff_mod2)

diff_mod3 <- lm(l2 ~ genodiff + year, data=phenoData)
anova(diff_mod3)

diff_mod4 <- lm(w/l1 ~ genodiff + year, data=phenoData)
anova(diff_mod4)

diff_mod5 <- lm(w/l2 ~ genodiff + year, data=phenoData)
anova(diff_mod5)

diff_mod6 <- lm(l1/l2 ~ genodiff + year, data=phenoData)
anova(diff_mod6)

tab_model(diff_mod1, diff_mod2, diff_mod3, diff_mod4, diff_mod5, diff_mod6, emph.p=T, show.aic=T, file=file.path(outpath, 'lm_genodiff.html'))
# Some 2015 individuals have lengths but not weight.
# All 2016 individuals have weight only

dmrt_mod1 <- lm(w ~ genodmrt + year, data=phenoData)
anova(dmrt_mod1)

dmrt_mod2 <- lm(l1 ~ genodmrt + year, data=phenoData)
anova(dmrt_mod2)

dmrt_mod3 <- lm(l2 ~ genodmrt + year, data=phenoData)
anova(dmrt_mod3)

dmrt_mod4 <- lm(w/l1 ~ genodmrt + year, data=phenoData)
anova(dmrt_mod4)

dmrt_mod5 <- lm(w/l2 ~ genodmrt + year, data=phenoData)
anova(dmrt_mod5)

dmrt_mod6 <- lm(l1/l2 ~ genodmrt + year, data=phenoData)
anova(dmrt_mod6)

tab_model(dmrt_mod1, dmrt_mod2, dmrt_mod3, dmrt_mod4, dmrt_mod5, dmrt_mod6, emph.p=T, show.aic=T, file=file.path(outpath, 'lm_genodmrt.html'))



## Morphometrics figures

str(frogdata) 

cyears <- c('2014', '2015')

for(m in 1:length(cyears)) {

morphodata <- subset(frogdata, frogdata$year==cyears[m]) 

## look for size difference of xx and xy individuals

par(mfrow=c(2,2))  
par(mar=c(5,5,4,3), oma=c(0,0,2,0))
plot(log(morphodata$w),log(morphodata$l2), type='n', main='Log w vs log l2', cex.main=1.8, cex.lab=1.3, xlab='log(w)', ylab='log(l2)')
points(log(morphodata$w[morphodata$genodiff=='undiff']), log(morphodata$l2[morphodata$genodiff=='undiff']), pch=19, col=cbred)
points(log(morphodata$w[morphodata$genodiff=='diff']), log(morphodata$l2[morphodata$genodiff=='diff']), pch=19, col=cbbluegreen)
points(log(morphodata$w[morphodata$genodiff=='semidiff']), log(morphodata$l2[morphodata$genodiff=='semidiff']), pch=19, col=cbpurple)

plot(log(morphodata$w),log(morphodata$l1), type='n', main='Log w vs log l1', cex.main=1.8, cex.lab=1.3, xlab='log(w)', ylab='log(l1)')
points(log(morphodata$w[morphodata$genodiff=='undiff']), log(morphodata$l1[morphodata$genodiff=='undiff']), pch=19, col=cbred)
points(log(morphodata$w[morphodata$genodiff=='diff']), log(morphodata$l1[morphodata$genodiff=='diff']), pch=19, col=cbbluegreen)
points(log(morphodata$w[morphodata$genodiff=='semidiff']), log(morphodata$l1[morphodata$genodiff=='semidiff']), pch=19, col=cbpurple)

plot(log(morphodata$l1),log(morphodata$l2), type='n', main='Log w vs log l1', cex.main=1.8, cex.lab=1.3,xlab='log(l1)', ylab='log(l1)')
points(log(morphodata$l1[morphodata$genodiff=='undiff']), log(morphodata$l2[morphodata$genodiff=='undiff']), pch=19, col=cbred)
points(log(morphodata$l1[morphodata$genodiff=='diff']), log(morphodata$l2[morphodata$genodiff=='diff']), pch=19, col=cbbluegreen)
points(log(morphodata$l1[morphodata$genodiff=='semidiff']), log(morphodata$l2[morphodata$genodiff=='semidiff']), pch=19, col=cbpurple)

frame()
legend('topleft', inset=0.05, legend=c('undiff', 'semidiff', 'diff'), pch =c(19,19,19), col=c(cbred,cbpurple,cbbluegreen))

mtext(paste('Linear by Y differentiation',cyears[m],""), outer = TRUE, cex = 1.5)

dev.copy(pdf,file.path(outpath, paste('genotype_linear_', cyears[m],'.pdf', sep='')), width=12, height=9)

dev.off()



par(mfrow=c(2,3)) 
par(mar=c(5,5,4,3), oma=c(0,0,2,0))
boxplot( morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='undiff'],
	morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='undiff'],
	morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff'],
	morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff'],
	morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='diff'],
	morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='diff'],
	names=list("A-undiff","N-undiff","A-semidiff","N-semidiff","A-diff","N-diff"), col=c(cbred,cbred,cbpurple,cbpurple,cbbluegreen,cbbluegreen), main ="w", xlab="amplexus", ylab="w", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

boxplot( morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='undiff'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='undiff'],
	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff'],
	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='diff'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='diff'],
	 names=list("A-undiff","N-undiff","A-semidiff","N-semidiff","A-diff","N-diff"), col=c(cbred,cbred,cbpurple,cbpurple,cbbluegreen,cbbluegreen), main ='l1', xlab="amplexus", ylab="l1", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)


boxplot( morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='undiff'],
	 morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='undiff'],
	 morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff'],
	 morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff'],
	 morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='diff'],
	 morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='diff'],
	names=list("A-undiff","N-undiff","A-semidiff","N-semidiff","A-diff","N-diff"), col=c(cbred,cbred,cbpurple,cbpurple,cbbluegreen,cbbluegreen), main ='l2', xlab="amplexus", ylab="l2", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

boxplot( morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='undiff']/morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='undiff'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='undiff']/morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='undiff'],
	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff']/morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff']/morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff'],
		 morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='diff']/morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='diff'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='diff']/morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='diff'],
 names=list("A-undiff","N-undiff","A-semidiff","N-semidiff","A-diff","N-diff"), col=c(cbred,cbred,cbpurple,cbpurple,cbbluegreen,cbbluegreen), main ='w/l1', xlab="amplexus", ylab="w/l1", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

boxplot( morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='undiff']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='undiff'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='undiff']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='undiff'],
	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff'],
	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodiff=='diff']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='diff'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodiff=='diff']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='diff'],
	 names=list("A-undiff","N-undiff","A-semidiff","N-semidiff","A-diff","N-diff"), col=c(cbred,cbred,cbpurple,cbpurple,cbbluegreen,cbbluegreen), main ='w/l2', xlab="amplexus", ylab="w/l2", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

boxplot( morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='undiff']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='undiff'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='undiff']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='undiff'],
	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='semidiff'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='semidiff'],
	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodiff=='diff']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodiff=='diff'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodiff=='diff']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodiff=='diff'],
	 names=list("A-undiff","N-undiff","A-semidiff","N-semidiff","A-diff","N-diff"), col=c(cbred,cbred,cbpurple,cbpurple,cbbluegreen,cbbluegreen), main ='l1/l2', xlab="amplexus", ylab="l1/l2", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

mtext(paste('Boxplots by Y differentiation',cyears[m],""), outer = TRUE, cex = 1.5)

dev.copy(pdf,file.path(outpath, paste('genotype_boxplots_', cyears[m],'.pdf', sep='')), width=14, height=8)

dev.off()


par(mfrow=c(1,2))	
par(mar=c(5,5,4,3), oma=c(0,0,2,0))

boxplot(morphodata$w[morphodata$amplexus=="A"], morphodata$w[morphodata$amplexus=="N"],
	morphodata$l1[morphodata$amplexus=="A"], morphodata$l1[morphodata$amplexus=="N"],
	morphodata$l2[morphodata$amplexus=="A"], morphodata$l2[morphodata$amplexus=="N"],
	names=list("A-w","N-w","A-l1","N-l1","A-l2","N-l2"), col=c(3,5), main ='A vs N', xlab="Phenotype", ylab="A vs N", cex.main=1.8, cex.lab=1.3,  las=2, cex.axis=0.65)
	
boxplot(morphodata$w[morphodata$amplexus=="A"]/morphodata$l1[morphodata$amplexus=="A"], morphodata$w[morphodata$amplexus=="N"]/morphodata$l1[morphodata$amplexus=="N"],
	morphodata$w[morphodata$amplexus=="A"]/morphodata$l2[morphodata$amplexus=="A"], morphodata$w[morphodata$amplexus=="N"]/morphodata$l2[morphodata$amplexus=="N"],
	morphodata$l1[morphodata$amplexus=="A"]/morphodata$l2[morphodata$amplexus=="A"], morphodata$l1[morphodata$amplexus=="N"]/morphodata$l2[morphodata$amplexus=="N"],
	names=list("A-w/l1","N-w/l1","A-w/l2","N-w/l2","A-l1/l2","N-l1/l2"), col=c(3,5), main ='A vs N', xlab="Phenotype (ratio)", ylab="A vs N", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

mtext(paste('A vs N phenotypes by Y differentiation',cyears[m],""), outer = TRUE, cex = 1.5)

dev.copy(pdf,file.path(outpath, paste('AvsN_', cyears[m],'.pdf', sep='')), width=9, height=8)
dev.off()


par(mfrow=c(1,1))	
par(mar=c(5,5,4,3))
boxplot(morphodata$w[morphodata$genodiff=="undiff"]/morphodata$l1[morphodata$genodiff=="undiff"], 
		morphodata$w[morphodata$genodiff=="semidiff"]/morphodata$l1[morphodata$genodiff=="semidiff"],
		morphodata$w[morphodata$genodiff=="diff"]/morphodata$l1[morphodata$genodiff=="diff"], 
		morphodata$w[morphodata$genodiff=="undiff"]/morphodata$l2[morphodata$genodiff=="undiff"],
		morphodata$w[morphodata$genodiff=="semidiff"]/morphodata$l2[morphodata$genodiff=="semidiff"],
		morphodata$w[morphodata$genodiff=="diff"]/morphodata$l2[morphodata$genodiff=="diff"],
		morphodata$l1[morphodata$genodiff=="undiff"]/morphodata$l2[morphodata$genodiff=="undiff"],
		morphodata$l1[morphodata$genodiff=="semidiff"]/morphodata$l2[morphodata$genodiff=="semidiff"],
		morphodata$l1[morphodata$genodiff=="diff"]/morphodata$l2[morphodata$genodiff=="diff"],
		names=list("Undiff\nw/l1","semidiff\nw/l1","Diff\nw/l1","Undiff\nw/l2","Semidiff\nw/l2","Diff\nw/l2","Undiff\nl1/l2","Semidiff\nl1/l2","Diff\nl1/l2"), col=c(cbred,cbpurple,cbbluegreen), main =paste('Differentiation effect on morphometric ratios', cyears[m]), xlab="Subset", ylab="Phenotype ratio", cex.main=1.8, cex.lab=1.3,las=2, cex.axis=0.65)

dev.copy(pdf,file.path(outpath, paste('genotype_Ratio_', cyears[m],'.pdf', sep='')), width=9, height=9)
dev.off()

}

cyears <- c('2014', '2015')

for(m in 1:length(cyears)) {

morphodata <- subset(frogdata, frogdata$year==cyears[m]) 

## look for size YB2erence of xx and xy individuals

par(mfrow=c(2,2))  
par(mar=c(5,5,4,3), oma=c(0,0,2,0))
plot(log(morphodata$w),log(morphodata$l2), type='n', main='Log w vs log l2', cex.main=1.8, cex.lab=1.3, xlab='log(w)', ylab='log(l2)')
points(log(morphodata$w[morphodata$genodmrt=='XX']), log(morphodata$l2[morphodata$genodmrt=='XX']), pch=19, col=col.XX_M)
points(log(morphodata$w[morphodata$genodmrt=='YB1']), log(morphodata$l2[morphodata$genodmrt=='YB1']), pch=19, col=col.YB1a)
points(log(morphodata$w[morphodata$genodmrt=='YB2']), log(morphodata$l2[morphodata$genodmrt=='YB2']), pch=19, col=col.YB2a)
points(log(morphodata$w[morphodata$genodmrt=='YBother']), log(morphodata$l2[morphodata$genodmrt=='YBother']), pch=19, col=col.YB40)

plot(log(morphodata$w),log(morphodata$l1), type='n', main='Log w vs log l1', cex.main=1.8, cex.lab=1.3, xlab='log(w)', ylab='log(l1)')
points(log(morphodata$w[morphodata$genodmrt=='XX']), log(morphodata$l1[morphodata$genodmrt=='XX']), pch=19, col=col.XX_M)
points(log(morphodata$w[morphodata$genodmrt=='YB1']), log(morphodata$l1[morphodata$genodmrt=='YB1']), pch=19, col=col.YB1a)
points(log(morphodata$w[morphodata$genodmrt=='YB2']), log(morphodata$l1[morphodata$genodmrt=='YB2']), pch=19, col=col.YB2a)
points(log(morphodata$w[morphodata$genodmrt=='YBother']), log(morphodata$l1[morphodata$genodmrt=='YBother']), pch=19, col=col.YB40)

plot(log(morphodata$l1),log(morphodata$l2), type='n', main='Log w vs log l1', cex.main=1.8, cex.lab=1.3,xlab='log(l1)', ylab='log(l1)')
points(log(morphodata$l1[morphodata$genodmrt=='XX']), log(morphodata$l2[morphodata$genodmrt=='XX']), pch=19, col=col.XX_M)
points(log(morphodata$l1[morphodata$genodmrt=='YB1']), log(morphodata$l2[morphodata$genodmrt=='YB1']), pch=19, col=col.YB1a)
points(log(morphodata$l1[morphodata$genodmrt=='YB2']), log(morphodata$l2[morphodata$genodmrt=='YB2']), pch=19, col=col.YB2a)
points(log(morphodata$l1[morphodata$genodmrt=='YBother']), log(morphodata$l2[morphodata$genodmrt=='YBother']), pch=19, col=col.YB40)

frame()
legend('topleft', inset=0.05, legend=c('XX', 'YB1', 'YB2', 'YBother'), pch =c(19,19,19), col=c(col.XX_M,col.YB1a,col.YB2a, col.YB40))

mtext(paste('Linear by Y dmrt type',cyears[m],""), outer = TRUE, cex = 1.5)

dev.copy(pdf,file.path(outpath, paste('Dmrtype_linear_', cyears[m],'.pdf', sep='')), width=12, height=9)

dev.off()



par(mfrow=c(2,3)) 
par(mar=c(5,5,4,3), oma=c(0,0,2,0))

boxplot( morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='XX'],
	morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='XX'],
	morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1'],
	morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1'],
	morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2'],
	morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2'],
	morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother'],
	morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother'],
	names=list("A-XX","N-XX","A-YB1","N-YB1","A-YB2","N-YB2","A-YBother","N-YBother"), col=c(col.XX_M,col.XX_M, col.YB1a,col.YB1a, col.YB2a,col.YB2a, col.YB40,col.YB40), main ="w", xlab="amplexus", ylab="w", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

boxplot( morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='XX'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='XX'],
	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1'],
	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2'],
 	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother'],
	 names=list("A-XX","N-XX","A-YB1","N-YB1","A-YB2", "N-YB2","A-YBother", "N-YBother"), col=c(col.XX_M,col.XX_M, col.YB1a,col.YB1a, col.YB2a,col.YB2a, col.YB40,col.YB40), main ='l1', xlab="amplexus", ylab="l1", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)


boxplot( morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='XX'],
	 morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='XX'],
	 morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1'],
	 morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1'],
	 morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2'],
	 morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2'],
 	 morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother'],
	 morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother'],
	 names=list("A-XX","N-XX","A-YB2","N-YB2","A-YB2","N-YB2", "A-YBother","N-YBother"), col=c(col.XX_M,col.XX_M, col.YB1a,col.YB1a, col.YB2a,col.YB2a, col.YB40,col.YB40), main ='l2', xlab="amplexus", ylab="l2", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

boxplot( morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='XX']/morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='XX'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='XX']/morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='XX'],
	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1']/morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1']/morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1'],
	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2']/morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2']/morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2'],
 	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother']/morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother']/morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother'],
	 names=list("A-XX","N-XX","A-YB2","N-YB2","A-YB2","N-YB2","A-YBother","N-YBother"), col=c(col.XX_M,col.XX_M, col.YB1a,col.YB1a, col.YB2a,col.YB2a, col.YB40,col.YB40), main ='w/l1', xlab="amplexus", ylab="w/l1", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

boxplot( morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='XX']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='XX'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='XX']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='XX'],
	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1'],
	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2'],
	 	 morphodata$w[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother'],
	 morphodata$w[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother'],
	 names=list("A-XX","N-XX","A-YB2","N-YB2","A-YB2","N-YB2", "A-YBother","N-YBother"), col=c(col.XX_M,col.XX_M, col.YB1a,col.YB1a, col.YB2a,col.YB2a, col.YB40,col.YB40), main ='w/l2', xlab="amplexus", ylab="w/l2", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

boxplot( morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='XX']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='XX'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='XX']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='XX'],
	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YB1'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YB1'],
	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YB2'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YB2'],
 	 morphodata$l1[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother']/morphodata$l2[morphodata$amplexus=="A" & morphodata$genodmrt=='YBother'],
	 morphodata$l1[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother']/morphodata$l2[morphodata$amplexus=="N" & morphodata$genodmrt=='YBother'],
	 names=list("A-XX","N-XX","A-YB2","N-YB2","A-YB2","N-YB2","A-YBother","N-YBother"), col=c(col.XX_M,col.XX_M, col.YB1a,col.YB1a, col.YB2a,col.YB2a, col.YB40, col.YB40), main ='l1/l2', xlab="amplexus", ylab="l1/l2", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

mtext(paste('Boxplots by Y dmrt type', cyears[m],""), outer = TRUE, cex = 1.5)

dev.copy(pdf,file.path(outpath, paste('dmrtype_boxplots_', cyears[m],'.pdf', sep='')), width=14, height=8)

dev.off()


par(mfrow=c(1,2))	
par(mar=c(5,5,4,3), oma=c(0,0,2,0))

boxplot(morphodata$w[morphodata$amplexus=="A"], morphodata$w[morphodata$amplexus=="N"],
	morphodata$l1[morphodata$amplexus=="A"], morphodata$l1[morphodata$amplexus=="N"],
	morphodata$l2[morphodata$amplexus=="A"], morphodata$l2[morphodata$amplexus=="N"],
	names=list("A-w","N-w","A-l1","N-l1","A-l2","N-l2"), col=c(3,5), main =paste('A vs N', cyears[m]), xlab="Phenotype", ylab="A vs N", cex.main=1.8, cex.lab=1.3,  las=2, cex.axis=0.65)
	
boxplot(morphodata$w[morphodata$amplexus=="A"]/morphodata$l1[morphodata$amplexus=="A"], morphodata$w[morphodata$amplexus=="N"]/morphodata$l1[morphodata$amplexus=="N"],
	morphodata$w[morphodata$amplexus=="A"]/morphodata$l2[morphodata$amplexus=="A"], morphodata$w[morphodata$amplexus=="N"]/morphodata$l2[morphodata$amplexus=="N"],
	morphodata$l1[morphodata$amplexus=="A"]/morphodata$l2[morphodata$amplexus=="A"], morphodata$l1[morphodata$amplexus=="N"]/morphodata$l2[morphodata$amplexus=="N"],
	names=list("A-w/l1","N-w/l1","A-w/l2","N-w/l2","A-l1/l2","N-l1/l2"), col=c(3,5), main = paste('A vs N', cyears[m]), xlab="Phenotype (ratio)", ylab="A vs N", cex.main=1.8, cex.lab=1.3, las=2, cex.axis=0.65)

mtext(paste('A vs N phenotypes by Y dmrt type',cyears[m],""), outer = TRUE, cex = 1.5)

dev.copy(pdf,file.path(outpath, paste('AvsN_', cyears[m],'.pdf', sep='')), width=9, height=8)
dev.off()


par(mfrow=c(1,1))	
par(mar=c(5,5,4,3))
boxplot(morphodata$w[morphodata$genodmrt=="XX"]/morphodata$l1[morphodata$genodmrt=="XX"], 
		morphodata$w[morphodata$genodmrt=="YB1"]/morphodata$l1[morphodata$genodmrt=="YB1"], 
		morphodata$w[morphodata$genodmrt=="YB2"]/morphodata$l1[morphodata$genodmrt=="YB2"],
		morphodata$w[morphodata$genodmrt=="YBother"]/morphodata$l1[morphodata$genodmrt=="YBother"],
		morphodata$w[morphodata$genodmrt=="XX"]/morphodata$l2[morphodata$genodmrt=="XX"],
		morphodata$w[morphodata$genodmrt=="YB1"]/morphodata$l2[morphodata$genodmrt=="YB1"],
		morphodata$w[morphodata$genodmrt=="YB2"]/morphodata$l2[morphodata$genodmrt=="YB2"],
		morphodata$w[morphodata$genodmrt=="YBother"]/morphodata$l2[morphodata$genodmrt=="YBother"],
		morphodata$l1[morphodata$genodmrt=="XX"]/morphodata$l2[morphodata$genodmrt=="XX"],
		morphodata$l1[morphodata$genodmrt=="YB1"]/morphodata$l2[morphodata$genodmrt=="YB1"],
		morphodata$l1[morphodata$genodmrt=="YB2"]/morphodata$l2[morphodata$genodmrt=="YB2"],
		morphodata$l1[morphodata$genodmrt=="YBother"]/morphodata$l2[morphodata$genodmrt=="YBother"],
		names=list("XX-w/l1","YB1-w/l1","YB2-w/l1","Yother-w/l1", "XX-w/l2","YB1-w/l2","YB2-w/l2","Yother-w/l2", "XX-l1/l2","YB1-l1/l2","YB2-l1/l2","Yother-l1/l2"), col=c(col.XX_M, col.YB1a, col.YB2a, col.YB40), main =paste('dmrt type effect on morphometric ratios', cyears[m]), xlab="Subset", ylab="Phenotype ratio", cex.main=1.8, cex.lab=1.3,las=2, cex.axis=0.65)

dev.copy(pdf,file.path(outpath, paste('dmrtype_Ratio_', cyears[m],'.pdf', sep='')), width=9, height=9)
dev.off()

}


cyears <- c('2014', '2015')

for(m in 1:length(cyears)) {

morphodata <- subset(frogdata, frogdata$year==cyears[m]) 

morphoplexus_data <- subset(morphodata, morphodata$amplexus!='D' & morphodata$amplexus!='AG' & morphodata$genodmrt!='YA1a')
morphoplexus_data$day <- factor(morphoplexus_data$day)
morphoplexus_data$genodmrt <- factor(morphoplexus_data$genodmrt)


par(mfrow=c(1,4)) 
par(mar=c(5,5,4,3)) 
undiff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='undiff']))
semidiff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='semidiff']))
diff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='diff']))
genodiff_day <- data.frame(undiff_day, semidiff_day, diff_day)
barplot(t(as.matrix(genodiff_day)), beside=T, names=levels(morphoplexus_data$day), main="All", cex.main=1.8, cex.lab=1.3, col=c(cbred,cbpurple,cbbluegreen))

undiff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='undiff' & morphoplexus_data$amplexus=='A']))
semidiff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='semidiff' & morphoplexus_data$amplexus=='A']))
diff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='diff' & morphoplexus_data$amplexus=='A']))
genodiff_day <- data.frame(undiff_day, semidiff_day, diff_day)
barplot(t(as.matrix(genodiff_day)), beside=T, names=levels(morphoplexus_data$day), main="Amplexus", cex.main=1.8, cex.lab=1.3, col=c(cbred,cbpurple,cbbluegreen))

undiff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='undiff' & morphoplexus_data$amplexus=='N']))
semidiff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='semidiff' & morphoplexus_data$amplexus=='N']))
diff_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodiff=='diff' & morphoplexus_data$amplexus=='N']))
genodiff_day <- data.frame(undiff_day, semidiff_day, diff_day)
barplot(t(as.matrix(genodiff_day)), beside=T, names=levels(droplevels(morphoplexus_data$day)), main="Non-amplexus", cex.main=1.8, cex.lab=1.3, col=c(cbred,cbpurple,cbbluegreen))

frame()
legend('topleft', inset=0.05, legend=c('undiff', 'semidiff', 'diff'), pch =c(15), col=c(cbred,cbpurple,cbbluegreen) )

mtext(paste('Differentiation by day',cyears[m],""), outer = TRUE, cex = 1.5)

dev.copy(pdf,file.path(outpath,paste('Differentiation_byday_', cyears[m], '.pdf', sep='')), width=16, height=5)
dev.off()

par(mfrow=c(1,4)) 
par(mar=c(5,5,4,3)) 
XX_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='XX']))
YB1_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YB1']))
YB2_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YB2']))
YBother_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YBother']))
genodmrt_day <- data.frame(XX_day, YB1_day, YB2_day, YBother_day)
barplot(t(as.matrix(genodmrt_day)), beside=T, names=levels(morphoplexus_data$day), main="All", cex.main=1.8, cex.lab=1.3, col=c(col.XX_M, col.YB1a, col.YB2a, col.YB40))

XX_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='XX' & morphoplexus_data$amplexus=='A']))
YB1_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YB1' & morphoplexus_data$amplexus=='A']))
YB2_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YB2' & morphoplexus_data$amplexus=='A']))
YBother_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YBother' & morphoplexus_data$amplexus=='A']))
genodmrt_day <- data.frame(XX_day, YB1_day, YB2_day, YBother_day)
barplot(t(as.matrix(genodmrt_day)), beside=T, names=levels(morphoplexus_data$day), main="Amplexus", cex.main=1.8, cex.lab=1.3, col=c(col.XX_M, col.YB1a, col.YB2a, col.YB40))

XX_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='XX' & morphoplexus_data$amplexus=='N']))
YB1_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YB1' & morphoplexus_data$amplexus=='N']))
YB2_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YB2' & morphoplexus_data$amplexus=='N']))
YBother_day <- as.integer(summary(morphoplexus_data$day[morphoplexus_data$genodmrt=='YBother' & morphoplexus_data$amplexus=='N']))
genodmrt_day <- data.frame(XX_day, YB1_day, YB2_day, YBother_day)
barplot(t(as.matrix(genodmrt_day)), beside=T, names=levels(morphoplexus_data$day), main="non-Amplexus", cex.main=1.8, cex.lab=1.3, col=c(col.XX_M, col.YB1a, col.YB2a, col.YB40))

frame()
legend('topleft', inset=0.05, legend=c('XX', 'XYB1', 'XYB2', 'XYBother'), pch =c(15), col=c(col.XX_M, col.YB1a, col.YB2a, col.YB40) )

mtext(paste('Dmrt by day',cyears[m],""), outer = TRUE, cex = 1.5)

dev.copy(pdf,file.path(outpath,paste('Dmrt_byday_', cyears[m], '.pdf', sep='')), width=16, height=5)
dev.off()


}




