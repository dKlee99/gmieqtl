#install.packages("reshape")
library(reshape)
a=read.table("C:/Users/dakyung/Desktop/Pathway/0627/GNO/TMM_log_anoikis.Counts.table.txt")
x=a[,c(1:12)] #LRO
#x=a[,c(1:5,11:15)] #TGFb_LRO
#x=a[,c(1:5,11:15,29,30)] #TGFb_LRO_043_included
#x=a[,c(1:10)] #Wnt_LRO
#x=a[,c(1:10,28,29)] #Wnt_LRO
#x=a[,c(31,32,34,35,37,38,40,41,43,44,52,53)] #GNO_Wnt
#x=a[,c(1:10)] #LRO
#x=a[,c(1:10)] #GNO
#x=a[,c(31,33,34,36,37,39,40,42,43,45,52,54)] #GNO_TGF
pca=prcomp(t(x),scale=F)
eigs=pca$sdev^2
proportion=eigs/sum(eigs)
proportion[1:5]
pca$x[,1:5]

#pca=prcomp(t(log10(x+1)))
#eigs=pca$sdev^2
#proportion=eigs/sum(eigs)
#proportion[1:5]
#pca$x[,1:5]

#project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
#barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

#par(cex=1.0, cex.axis=0.8, cex.main=0.8)
#pairs(project.pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
#pairs(project.pca$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
#par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)

#Plots scatter plot for PC 1 and 2
#plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
#points(project.pca$x, col="black", pch=16, cex=1)

#Plots scatter plot for PC 1 and 3
#plot(project.pca$x[,1], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
#points(project.pca$x[,1], project.pca$x[,3], col="black", pch=16, cex=1)

#Plots scatter plot for PC 2 and 3
#plot(project.pca$x[,2], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
#points(project.pca$x[,2], project.pca$x[,3], col="black", pch=16, cex=1)

t <- data.frame(sort(abs(pca$rotation[,"PC2"]), decreasing=TRUE)[1:50])

hm = cbind(gene = rownames(t), t)
head(hm, n=20)
hmm = melt(hm, id = c("gene"))
head(hmm)

write.table(hm,"C:/Users/dakyung/Desktop/Pathway/0627/LRO/PC1_TGF.txt")
