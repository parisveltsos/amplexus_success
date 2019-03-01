library(adegenet)
datapath <- '~/git/amplexus_success/input'
outpath <- '~/git/amplexus_success/output/pca'

dir.create(outpath)
# Manually made from dropbox all genotype data, *, 260-2, 360-2, ? are some of the issues in that file that needs manual attention

metData <- read.table(file.path(datapath, 'pca_data.txt'), header=T)
str(metData)

# remove very rate YA1a genotypes for good resolution in the graph
metData_subset <- subset(metData, metData$genotype!='YA1a') 
genodata <- as.matrix(metData_subset[5:13], dimnames=list(paste(colnames(metData_subset[5:13])))) # all markers

obj <- df2genind(genodata, ploidy=1, sep=",")
obj@tab
obj@pop <- metData_subset$genotype
# check conversion worked by backconverting
genind2df(obj, sep="")

# setup colours
col.XX   <- rgb(206/255, 34/255, 43/255)
col.XX_F  <- rgb(249/255, 187/255, 235/255)
col.YA1a <- rgb(209/255, 127/255, 21/255)
col.YB10 <- rgb(100/255, 179/255, 223/255)
col.YB1a <- rgb(23/255, 87/255, 120/255)
col.YB20 <- rgb(156/255, 225/255, 89/255)
col.YB2a <- rgb(88/255, 135/255, 37/255)
col.YB2b <- rgb(122/255, 174/255, 61/255)
col.YB30 <- rgb(0/255, 222/255, 217/255)
col.YB40 <- rgb(113/255, 250/255, 241/255)
col.YB50 <- rgb(205/255, 243/255, 244/255)
col.scheme <- c(col.XX, col.XX_F, col.YA1a, col.YB10, col.YB1a, col.YB20, col.YB2a, col.YB2b, col.YB30, col.YB40, col.YB50)



X <- scaleGen(obj, NA.method='mean')
class(X)
dim(X)

pca1 <- dudi.pca(X, cent=FALSE, scale=FALSE, scannf=FALSE, nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
pca1
# s.label(pca1$li)
# title("PCA axes 1-2") 
# add.scatter.eig(pca1$eig[1:20], 3,1,2)

dev.copy(pdf,file.path(outpath,'pca_3pages.pdf'), width=12, height=12)
par(mfrow=c(2,1)) 
par(mar=c(5,5,4,3), oma=c(0,0,2,0))
s.class(pca1$li, pop(obj), col=transp(col.scheme,.6))
title("PCA\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li,pop(obj), xax=1, yax=3, csub=2, col=transp(col.scheme,.6)) 
title("PCA\naxes 1-3")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)

s.class(pca1$li, pop(obj), xax=1, yax=2, col=transp(col.scheme,.6), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
title("PCA\naxes 1-2")
s.class(pca1$li, pop(obj), xax=1, yax=3, col=transp(col.scheme,.6), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
title("PCA\naxes 1-3")
mtext("Meitreile", outer = TRUE, cex = 1.5)
dev.off()

legend('topleft', inset=0.05, legend=c('opa','eipa'), pch =15, col=2 )

legend('topleft', inset=0.05, legend=levels(pop(obj)), pch =15, col=transp(col.scheme,.6) )

# clean plot
par(mar=c(5,5,4,3))
plot(pca1$li$Axis1, pca1$li$Axis2, type='n', xlab="PC1", ylab="PC2", main="Meitreile", cex.main=1.8, cex.lab=1.3)
points(pca1$li$Axis1[metData_subset$genotype=='XX'], pca1$li$Axis2[metData_subset$genotype=='XX'], pch=16, col=transp(col.XX, .6))
points(pca1$li$Axis1[metData_subset$genotype=='XX_F'], pca1$li$Axis2[metData_subset$genotype=='XX_F'], pch=16, col=transp(col.XX_F, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YA1a'], pca1$li$Axis2[metData_subset$genotype=='YA1a'], pch=16, col=transp(col.YA1a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB10'], pca1$li$Axis2[metData_subset$genotype=='YB10'], pch=16, col=transp(col.YB10, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB1a'], pca1$li$Axis2[metData_subset$genotype=='YB1a'], pch=16, col=transp(col.YB1a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB20'], pca1$li$Axis2[metData_subset$genotype=='YB20'], pch=16, col=transp(col.YB20, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB2a'], pca1$li$Axis2[metData_subset$genotype=='YB2a'], pch=16, col=transp(col.YB2a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB2b'], pca1$li$Axis2[metData_subset$genotype=='YB2b'], pch=16, col=transp(col.YB2b, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB30'], pca1$li$Axis2[metData_subset$genotype=='YB30'], pch=16, col=transp(col.YB30, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB40'], pca1$li$Axis2[metData_subset$genotype=='YB40'], pch=16, col=transp(col.YB40, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB50'], pca1$li$Axis2[metData_subset$genotype=='YB50'], pch=16, col=transp(col.YB50, .6))

legend('topright', inset=0.05, legend=levels(pop(obj)), pch =16, col=transp(col.scheme,.6), cex=0.5 )

dev.copy(pdf,file.path(outpath,'pca_clean.pdf'), width=9, height=9)
dev.off()

# plot including year as shape (illustrates no difference between years)
par(mar=c(5,5,4,3))
plot(pca1$li$Axis1, pca1$li$Axis2, type='n', xlab="PC1", ylab="PC2", main="Meitreile", cex.main=1.8, cex.lab=1.3)
points(pca1$li$Axis1[metData_subset$genotype=='XX' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='XX' & metData_subset$year=='2016'], pch=1, col=transp(col.XX, .6))
points(pca1$li$Axis1[metData_subset$genotype=='XX' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='XX' & metData_subset$year=='2015'], pch=2, col=transp(col.XX, .6))
points(pca1$li$Axis1[metData_subset$genotype=='XX' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='XX' & metData_subset$year=='2014'], pch=3, col=transp(col.XX, .6))
points(pca1$li$Axis1[metData_subset$genotype=='XX_F' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='XX_F' & metData_subset$year=='2016'], pch=1, col=transp(col.XX_F, 1))
points(pca1$li$Axis1[metData_subset$genotype=='XX_F' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='XX_F' & metData_subset$year=='2015'], pch=2, col=transp(col.XX_F, 1))
points(pca1$li$Axis1[metData_subset$genotype=='XX_F' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='XX_F' & metData_subset$year=='2014'], pch=3, col=transp(col.XX_F,1))
points(pca1$li$Axis1[metData_subset$genotype=='YA1a' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YA1a' & metData_subset$year=='2016'], pch=1, col=transp(col.YA1a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YA1a' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YA1a' & metData_subset$year=='2015'], pch=2, col=transp(col.YA1a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YA1a' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YA1a' & metData_subset$year=='2014'], pch=3, col=transp(col.YA1a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB10' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YB10' & metData_subset$year=='2016'], pch=1, col=transp(col.YB10, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB10' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YB10' & metData_subset$year=='2015'], pch=2, col=transp(col.YB10, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB10' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YB10' & metData_subset$year=='2014'], pch=3, col=transp(col.YB10, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB1a' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YB1a' & metData_subset$year=='2016'], pch=1, col=transp(col.YB1a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB1a' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YB1a' & metData_subset$year=='2015'], pch=2, col=transp(col.YB1a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB1a' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YB1a' & metData_subset$year=='2014'], pch=3, col=transp(col.YB1a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB20'& metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YB20'& metData_subset$year=='2016'], pch=1, col=transp(col.YB20, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB20'& metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YB20'& metData_subset$year=='2015'], pch=2, col=transp(col.YB20, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB20' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YB20'& metData_subset$year=='2014'], pch=3, col=transp(col.YB20, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB2a' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YB2a' & metData_subset$year=='2016'], pch=1, col=transp(col.YB2a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB2a' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YB2a' & metData_subset$year=='2015'], pch=2, col=transp(col.YB2a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB2a' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YB2a' & metData_subset$year=='2014'], pch=3, col=transp(col.YB2a, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB2b' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YB2b' & metData_subset$year=='2016'], pch=1, col=transp(col.YB2b, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB2b' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YB2b' & metData_subset$year=='2015'], pch=2, col=transp(col.YB2b, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB2b' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YB2b' & metData_subset$year=='2014'], pch=3, col=transp(col.YB2b, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB30' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YB30' & metData_subset$year=='2016'], pch=1, col=transp(col.YB30, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB30' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YB30' & metData_subset$year=='2015'], pch=2, col=transp(col.YB30, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB30' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YB30' & metData_subset$year=='2014'], pch=3, col=transp(col.YB30, .6))
points(pca1$li$Axis1[metData_subset$genotype=='YB40' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YB40' & metData_subset$year=='2016'], pch=1, col=transp(col.YB40, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB40' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YB40' & metData_subset$year=='2015'], pch=2, col=transp(col.YB40, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB40' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YB40' & metData_subset$year=='2014'], pch=3, col=transp(col.YB40, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB50' & metData_subset$year=='2016'], pca1$li$Axis2[metData_subset$genotype=='YB50' & metData_subset$year=='2016'], pch=1, col=transp(col.YB50, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB50' & metData_subset$year=='2015'], pca1$li$Axis2[metData_subset$genotype=='YB50' & metData_subset$year=='2015'], pch=2, col=transp(col.YB50, 1))
points(pca1$li$Axis1[metData_subset$genotype=='YB50' & metData_subset$year=='2014'], pca1$li$Axis2[metData_subset$genotype=='YB50' & metData_subset$year=='2014'], pch=3, col=transp(col.YB50, 1))
legend('topright', inset=0.05, legend=levels(pop(obj)), pch =16, col=transp(col.scheme,.6), cex=0.5 )
legend('bottomright', inset=0.05, legend=c('2014','2015','2016'), pch =c(1,2,3), col=1, cex=0.5 )

dev.copy(pdf,file.path(outpath,'pca_by_year.pdf'), width=9, height=9)
dev.off()


# colour shows axis 3
# colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2") 
# title("PCA\naxes 1-2")
# abline(v=0,h=0,col="grey", lty=2)
# colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
# title("PCA\naxes 1-3")
# abline(v=0,h=0,col="grey", lty=2)


