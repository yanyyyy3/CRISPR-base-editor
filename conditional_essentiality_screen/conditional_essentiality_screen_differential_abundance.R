library("ggplot2")
# library("RUVSeq")
library("edgeR")
library("plyr")
library(scales)
library(stringr)
library(reshape2)
library(RColorBrewer)
# library("forcats")
tab <- read.table("sample_counts.csv", sep="\t", header=T, stringsAsFactors = F, quote="\"")
rownames(tab) <- tab$guide_name

gene <- vapply(strsplit(unlist(strsplit(tab$guide_name,",")),"_"),'[',1,FUN.VALUE=character(1))
gene <- str_replace(gene,"random.*","NC")
gene <- unique(gene)
length(gene)

pdf(file = "sample_check.pdf")
cormat <- cor(log(tab[,4:length(tab[1,])] + 1))
m <- cormat[hclust(dist(cormat))$order ,hclust(dist(t(cormat)))$order]
rows <- dim(m)[1]
cols <- dim(m)[2]
melt.m <- cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows),
                reshape2::melt(m))
g <- ggplot2::ggplot(data=melt.m)

g <- g + ggplot2::geom_rect(ggplot2::aes(xmin=colInd-1,xmax=colInd,
                                         ymin=rowInd-1,ymax=rowInd, fill=value)
                            ,colour='white')

g <- g + ggplot2::scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
g <- g + ggplot2::scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))

g <- g + ggplot2::theme(panel.grid.minor=ggplot2::element_line(colour=NA),
                        panel.grid.major=ggplot2::element_line(colour=NA),
                        panel.background=ggplot2::element_rect(fill=NA,
                                                               colour=NA),axis.text.x = element_text(angle = 90))

print(g)
colors <-c('#588cc0', '#a3d3e6', '#e9f6e8', '#fee99d', '#fca55d')
plot_input <- function(guide_names, counts, classes,title){
    s1t1 <- data.frame(cbind(X1=guide_names, X2=counts, X3=classes))
    s1t1$X3 <- factor(classes, levels=c('essential_in_m9','essential_in_mops', 'essential_in_both', 'non-essential','NT'))
    s1t1$X2 <- as.numeric(counts)
    r <- ggplot(s1t1, aes(x=reorder(X1, -X2), y=log(X2+1, base=10), colour=X3, fill=X3)) + geom_bar(stat="identity",width = 0.001) + theme_classic() + theme(legend.title = element_blank(),axis.title.x=element_blank(),
                                                                                                                                                                     axis.text.x=element_blank(),
                                                                                                                                                                     axis.line.x =element_blank(), axis.line.y =element_blank() ,axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5)) + ylab("log10 counts") +ggtitle(title)+ scale_color_manual(values = colors) + scale_fill_manual(values = colors)
    print(r)
}

for (sample in c(  "initial_library",'LB_1', 'LB_2', 'BE_1', 'BE_2', 'killing_1', 'killing_2', 'M9_1', 'M9_2', 'MOPS_1', 'MOPS_2')){
  plot_input(tab$guide_name, unlist(tab[sample]),tab$character,sample)
}
ess <- tab[grep("essential_in",tab$character),]
for (sample in c(  "initial_library",'LB_1', 'LB_2', 'BE_1', 'BE_2', 'killing_1', 'killing_2', 'M9_1', 'M9_2', 'MOPS_1', 'MOPS_2')){
  plot_input(ess$guide_name, unlist(ess[sample]),ess$character,sample)
}

dev.off()


plot_logFC <- function(guide_names, logFC, classes,title){
    s1t1 <- data.frame(cbind(X1=guide_names, X2=logFC, X3=classes))
    s1t1$X3 <- factor(classes, levels=c('essential_in_m9','essential_in_mops', 'essential_in_both', 'non-essential','NT'))
    s1t1$X2 <- as.numeric(logFC)
    r <- ggplot(s1t1, aes(x=reorder(X1, -X2), y=X2, colour=X3, fill=X3)) + geom_bar(stat="identity") + theme_classic() + theme(legend.title = element_blank(), axis.title.x=element_blank(),
                                                                                                                                         axis.text.x=element_blank(), 
                                                                                                                                         axis.line.x =element_blank(), axis.line.y =element_blank() ,axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = colors) + scale_fill_manual(values =colors)+ ylab("logFC")+ggtitle(title)                      
    print(r)
}

plot_FDR <- function(guide_names, FDR, classes,title){
    s1t1 <- data.frame(cbind(X1=guide_names, X2=FDR, X3=classes))
    s1t1$X3 <- factor(classes, levels=c('essential_in_m9','essential_in_mops', 'essential_in_both', 'non-essential','NT'))
    s1t1$X2 <- as.numeric(FDR)
    r <- ggplot(s1t1, aes(x=reorder(X1, -X2), y=-log(X2, base=10), colour=s1t1$X3, fill=s1t1$X3)) + geom_bar(stat="identity") + theme_classic() + theme(legend.title = element_blank(), axis.title.x=element_blank(),
                                                                                                                                                        axis.text.x=element_blank(), 
                                                                                                                                                        axis.line.x =element_blank(), axis.line.y =element_blank() ,axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5)) + scale_color_manual(values =colors) + scale_fill_manual(values =colors) + xlab(length(guide_names))+ ylab("-log10 FDR") + ggtitle(title)                        
    print(r)
}
timepoint_bar <- function(value,title){
    timepoints<-c(rep("LB-initial",5),rep("BE-LB",5),rep("killing-BE",5),rep("M9-killing",5),rep("MOPS-killing",5))
    classes<-rep(c("essential_in_both","essential_in_m9","essential_in_mops","non-essential","NT"),5)
    value <-value
    data <- data.frame(timepoints,classes,value)
    data$timepoints <- as.character(data$timepoints)
    data$timepoints <- factor(data$timepoints, levels=c("LB-initial","BE-LB","killing-BE","M9-killing","MOPS-killing"))
    r <- ggplot(data,aes(fill=classes,y=value,x=timepoints))+ geom_bar(position="dodge", stat="identity")+ggtitle(paste(title," logFC between timepoints",sep=""))+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90))+ylab("logFC")
    print(r)
    write.table(as.data.frame(cbind(rownames(data),data)),file = paste("./",title,"_logFC_between_timepoints.csv",sep=""),sep = '\t',row.names = FALSE)
}
lcpm_group <- function(y,title) {
    nt <- cpm(y[y$genes$type=="NT",], log=TRUE, normalized.lib.sizes=TRUE)
    essential_in_both<- cpm(y[y$genes$type=="essential_in_both",], log=TRUE, normalized.lib.sizes=TRUE)
    essential_in_m9<- cpm(y[y$genes$type=="essential_in_m9",], log=TRUE, normalized.lib.sizes=TRUE)
    essential_in_mops<- cpm(y[y$genes$type=="essential_in_mops",], log=TRUE, normalized.lib.sizes=TRUE)
    noness <- cpm(y[y$genes$type=="non-essential",], log=TRUE, normalized.lib.sizes=TRUE)
    lcpm.m <- rbind(melt(essential_in_m9),melt(essential_in_mops),melt(essential_in_both),melt(noness),melt(nt))
    classes<-c(rep("essential_in_m9",dim(essential_in_m9)[1]*11),rep("essential_in_mops",dim(essential_in_mops)[1]*11),rep("essential_in_both",dim(essential_in_both)[1]*11),rep("non-essential",dim(noness)[1]*11),rep("NT",dim(nt)[1]*11))
    value <-lcpm.m$value
    samples <- lcpm.m$Var2
    data <- data.frame(samples,classes,value)
    r <- ggplot(data,aes(fill=samples,y=value,x=classes))+ geom_boxplot()+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90))+xlab("classes")+ylab("lcpm")
    print(r)
}
fraction_group <- function(y,title){
    nt <- cpm(y[y$genes$type=="NT",], log=TRUE, normalized.lib.sizes=TRUE)
    essential_in_both<- cpm(y[y$genes$type=="essential_in_both",], log=TRUE, normalized.lib.sizes=TRUE)
    essential_in_m9<- cpm(y[y$genes$type=="essential_in_m9",], log=TRUE, normalized.lib.sizes=TRUE)
    essential_in_mops<- cpm(y[y$genes$type=="essential_in_mops",], log=TRUE, normalized.lib.sizes=TRUE)
    noness <- cpm(y[y$genes$type=="non-essential",], log=TRUE, normalized.lib.sizes=TRUE)
    cs <- rbind(colSums(essential_in_m9),colSums(essential_in_mops),colSums(essential_in_both),colSums(noness),colSums(nt))
    cs <- cbind(library=c(dim(essential_in_m9)[1],dim(essential_in_mops)[1],dim(essential_in_both)[1],dim(noness)[1],dim(nt)[1]),cs) 
                #
    cs <- prop.table(cs,margin = 2) *100
    rownames(cs) <- c('essential_in_m9','essential_in_mops', 'essential_in_both', 'non-essential','NT')
    cs.m <- melt(cs)
    classes <- cs.m$Var1
    data <- data.frame(classes,cs.m$Var2,cs.m$value)
    r <- ggplot(data,aes(fill=classes,y=cs.m$value,x=cs.m$Var2))+ geom_bar(position="fill", stat="identity")+ggtitle(title)+theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))+xlab("samples")+ylab("percentage of cpm")
    print(r)
    # sink(file='./fraction_percentage.txt',append = TRUE)
    print(cs)
    # sink()
    write.table(as.data.frame(cbind(rownames(cs),cs)),file = "./group_fraction.csv",sep = '\t',row.names = FALSE)
}
cpm_logFC <- function(y,title){
    type <- y$genes$type
    y <-cpm(y$counts, log=TRUE, normalized.lib.sizes=TRUE, prior.count=2)
    cpm_cols <- c("LB_1-initial_libary","LB_2-initial_libary","BE_1-LB_1","BE_2-LB_2","killing_1-BE_1","killing_2-BE_2","M9_1-killing_1","M9_2-killing_2","MOSP_1-killing_1","MOSP_2-killing_2")   
    print(cpm_cols)
    logFC <- cbind(y[,2]-y[,1],y[,3]-y[,1],y[,4]-y[,2],y[,5]-y[,3],y[,6]-y[,4],y[,7]-y[,5],y[,8]-y[,6],y[,9]-y[,7],y[,10]-y[,6],y[,11]-y[,7])
    colnames(logFC) <-cpm_cols
    print(dim(logFC))
    nt <- logFC[which(type=="NT"),]
    nt <-nt[,1:length(cpm_cols)]
    ess<- logFC[which(type=="essential_in_both"),]
    ess <- ess[,1:length(cpm_cols)]
    essential_in_m9<- logFC[which(type=="essential_in_m9"),]
    essential_in_m9 <- essential_in_m9[,1:length(cpm_cols)]
    essential_in_mops<- logFC[which(type=="essential_in_mops"),]
    essential_in_mops <- essential_in_mops[,1:length(cpm_cols)]
    noness <- logFC[which(type=="non-essential"),]
    noness <- noness[,1:length(cpm_cols)]
    print(dim(melt(nt)))
    cpm.log2FC <- rbind(melt(ess),melt(essential_in_m9),melt(essential_in_mops),melt(noness),melt(nt))
    classes<-c(rep("essential_in_both",dim(ess)[1]*length(cpm_cols)),rep("essential_in_m9",dim(essential_in_m9)[1]*length(cpm_cols)),rep("essential_in_mops",dim(essential_in_mops)[1]*length(cpm_cols)),rep("non-essential",dim(noness)[1]*length(cpm_cols)),rep("NT",dim(nt)[1]*length(cpm_cols)))
    value <-cpm.log2FC$value
    samples <- cpm.log2FC$Var2
    data <- data.frame(samples,classes,value)
    r <- ggplot(data,aes(fill=samples,y=value,x=classes))+ geom_boxplot()+ggtitle(paste(title,"cpm log2FC between time points",sep=" "))+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90))+xlab("classes")+ylab("cpm logFC")
    print(r)
    d <- ggplot_build(r)
    log2FC_mean <- data.frame(aggregate(value,list("samples" = cpm.log2FC$Var2,'classes' = data$classes) ,  mean))
    r <- ggplot(log2FC_mean,aes(fill=samples,y=x,x=classes))+ geom_bar(position="dodge", stat="identity")+ggtitle(paste(title,"mean of cpm logFC",sep=" "))+theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90))+xlab("classes")+ylab("cpm logFC")
    print(r)

    log2FC_median <- data.frame(aggregate(value,list("samples" = cpm.log2FC$Var2,'classes' = data$classes) ,  median))
    # log2FC_median <- log(2**log2FC_median-1,base=2)
    # classes<-c(rep("NC",6),rep("essential",6),rep("non-essential",6))
    # data <- data.frame(log2FC_median[,1],classes,rownames(log2FC_median))
    r <- ggplot(log2FC_median,aes(fill=samples,y=x,x=classes))+ geom_bar(position="dodge", stat="identity")+ggtitle(paste(title,"median of cpm logFC",sep=" "))+theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90))+xlab("classes")+ylab("cpm logFC")
    print(r)
}

#build count mat s enzyme
count_mat <- tab
groups <- c('initial', 'LB', 'LB', 'BE', 'BE', 'killing', 'killing', 'M9', 'M9', 'MOPS', 'MOPS')
y <- DGEList(count_mat[,4:14], group=groups, genes=count_mat[,2])
y$genes$type<-count_mat[,1]

#filter by cpm
dim(y)
keep <- rowSums(cpm(y$counts) > 5) >= 4
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)


# try norm on non-targeting
y_n <- y
nt <- y[y$genes$type=="NT",]
dim(nt)
nt <- calcNormFactors(nt, method="TMM")
y_n$samples$norm.factors <- nt$samples$norm.factors
print(y_n$samples)

pdf(file = "cpm_plots.pdf")
#plots
fraction_group(y,"Percentage of counts of different classes")
# unnormalized barplot
lcpm <- cpm(y, log=TRUE, normalized.lib.sizes=TRUE)
boxplot(lcpm, las=2, main="")
title(main="Unnormalised data", ylab="Log-cpm")
#normalized bar plot
sink("norm.factors.txt",append = FALSE)
print(y_n$samples)
sink()
lcpm <- cpm(y_n, log=TRUE,normalized.lib.sizes=TRUE)
boxplot(lcpm, las=2, main="",pars=list(par(mar=c(8,4,4,2))))
title(main="Normalised data", ylab="Log-cpm")
write.table(lcpm,file="lcpm.csv",sep = "\t",col.names=NA)

lcpm_group(y,"unnormalized")
lcpm_group(y_n,"normalized")
cpm_logFC(y,"unnormalized")
cpm_logFC(y_n,"normalized")
dev.off()


#MDS libraries
pdf(file = "./dispersion.pdf")
# plotMDS(y)
plotMDS(y_n)
#design matrix
groups <- factor(groups)
design <- model.matrix(~0+groups)
colnames(design) <- sub("groups","",colnames(design))
# dispersions
y_n <- estimateDisp(y_n, design, robust=TRUE)
plotBCV(y_n)
fit <- glmQLFit(y_n, design, robust=TRUE)
plotQLDisp(fit)
dev.off()
# fit <- glmFit(y_n, design)

# differential expression analysis
mean_value<-vector()
median_value<-vector()
values_essential_in_both <- vector()
values_essential_in_m9 <- vector()
values_essential_in_mops <- vector()
values_noness <- vector()
values_nc <- vector()
pdf(file = "DE.pdf")
for (i in c("LB-initial","BE-LB","killing-BE","M9-killing","MOPS-killing")){
    # str <- "MOPS-killing"
    str <- i
    cont <- makeContrasts(str, levels=design)
    res <- glmQLFTest(fit, contrast=cont)
    tt <- topTags(res, n=Inf)
    
    tt_mod <- cbind(guides = tt$table[,1], type=tt$table[,2], tt$table[,3:7],stringsAsFactors=FALSE)
    write.table(tt_mod,file = paste(str,"_QLFTest.csv",sep=""),sep = "\t",row.names = FALSE)
    
    hist(tt_mod$PValue,breaks=100,main = paste("p value distribution (",str,")",sep = ""),xlab = "p value")
    plot(tt_mod[which(tt_mod$type=="non-essential"),"logCPM"], tt_mod[which(tt_mod$type=="non-essential"),"logFC"], pch=20, main=str, xlab="logCPM", ylab="logFC", col="grey")
    points(tt_mod[which(tt_mod$type=="NT"),"logCPM"], tt_mod[which(tt_mod$type=="NT"),"logFC"], pch=20, col=alpha("black", 0.4))
    points(tt_mod[which(tt_mod$type=="essential_in_both"),"logCPM"], tt_mod[which(tt_mod$type=="essential_in_both"),"logFC"], pch=20, col=alpha("red", 0.4))
    points(tt_mod[which(tt_mod$type=="essential_in_m9"),"logCPM"], tt_mod[which(tt_mod$type=="essential_in_m9"),"logFC"], pch=20, col=alpha("yellow", 0.4))
    points(tt_mod[which(tt_mod$type=="essential_in_mops"),"logCPM"], tt_mod[which(tt_mod$type=="essential_in_mops"),"logFC"], pch=20, col=alpha("green", 0.4))
    
    legend("bottomright", legend=c("essential_in_both","essential_in_M9","essential_in_MOPS","non-essential","NT"),col=c("red", "yellow","green","grey","black"), pch=c(20),cex=0.6)#     #points(sel_tab$dTA_TE, sel_tab$TA_TE, pch=21, col="red")
    plot_logFC(tt_mod$guides,tt_mod$logFC,tt_mod$type,str)
    
    for (j in c("essential_in_both","essential_in_m9","essential_in_mops","non-essential","NT")){
        print(dim(tt_mod[tt_mod$type==j,]))
        print(mean(tt_mod[tt_mod$type==j,]$logFC))
        print(median(tt_mod[tt_mod$type==j,]$logFC))
        mean_value<-append(mean_value,mean(tt_mod[tt_mod$type==j,]$logFC),after = length(mean_value))
        median_value<-append(median_value,median(tt_mod[tt_mod$type==j,]$logFC),after = length(median_value))
    }
    
    values_essential_in_both <- cbind(values_essential_in_both,tt_mod[which(tt_mod$type=="essential_in_both"),"logFC"])
    colnames(values_essential_in_both)[dim(values_essential_in_both)[2]] <- str
    values_essential_in_m9 <- cbind(values_essential_in_m9,tt_mod[which(tt_mod$type=="essential_in_m9"),"logFC"])
    colnames(values_essential_in_m9)[dim(values_essential_in_m9)[2]] <- str
    values_essential_in_mops <- cbind(values_essential_in_mops,tt_mod[which(tt_mod$type=="essential_in_mops"),"logFC"])
    colnames(values_essential_in_mops)[dim(values_essential_in_mops)[2]] <- str
    values_noness <- cbind(values_noness,tt_mod[which(tt_mod$type=="non-essential"),"logFC"])
    colnames(values_noness)[dim(values_noness)[2]] <- str
    values_nc <- cbind(values_nc,tt_mod[which(tt_mod$type=="NT"),"logFC"])
    colnames(values_nc)[dim(values_nc)[2]] <- str
    
}
timepoint_bar(mean_value,"mean")
timepoint_bar(median_value,"median")

tt_logFC <- rbind(melt(values_essential_in_both),melt(values_essential_in_m9),melt(values_essential_in_mops),melt(values_noness),melt(values_nc))
classes<-c(rep("essential_in_both",dim(values_essential_in_both)[1]*dim(values_essential_in_both)[2]),rep("essential_in_m9",dim(values_essential_in_m9)[1]*dim(values_essential_in_m9)[2]),rep("essential_in_mops",dim(values_essential_in_mops)[1]*dim(values_essential_in_mops)[2]),rep("non-essential",dim(values_noness)[1]*dim(values_noness)[2]),rep("NT",dim(values_nc)[1]*dim(values_nc)[2]))
value <-tt_logFC$value
data <- data.frame(tt_logFC$Var2,classes,value)
data$classes <- as.character(data$classes)
data$classes <- factor(data$classes, levels=c("essential_in_both","essential_in_m9","essential_in_mops","non-essential","NT"))
r <- ggplot(data,aes(fill=tt_logFC$Var2,y=value,x=classes))+ geom_boxplot()+ggtitle("logFC between timepoints")+theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90))+xlab("classes")+ylab("logFC")+ylim(-10,10)
print(r)
dev.off()



