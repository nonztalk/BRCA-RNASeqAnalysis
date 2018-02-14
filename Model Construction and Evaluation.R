# Load the data
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("DOSE")
#biocLite("clusterProfiler")
library(limma)
library(DOSE)
library(clusterProfiler)
# load("C:\\Users\\hengshi\\Box Sync\\2017 Fall\\BIOINF 527\\Group Project\\BRCA\\Matrix_for_Model.Rdata")
## Data Analysis
library("ggplot2")
library("GGally")
library("reshape2")
library("gplots")
library("glmnet") #
library("nnet")
library("MASS")
library("caret") #
library("gbm") #
library("tree") #
library("randomForest") #
library("Hmisc")
library("gridExtra")

# ind<-rep(NA,dim(matrixModel)[2] )
# 
# for(i in 1:dim(matrixModel)[2]){
#   ind[i]<-ifelse(class(matrixModel[, i])!="numeric",1,0)
#   }
# names(matrixModel[, ind==1])
# # "er_status_by_ihc"   "pr_status_by_ihc"   "her2_status_by_ihc"  "subtype"
# 
# ##Dimension-reduced unsupervised methods
# usedata <- matrixModel

# write.csv(usedata, row.names = FALSE, file="C:\\Users\\hengshi\\Box Sync\\2017 Fall\\BIOINF 527\\Group Project\\BRCA\\usedata.csv")
usedata <- read.csv("./data/matrixModel_6Group.csv", header=TRUE)
unique(usedata$er_status_by_ihc)
unique(usedata$pr_status_by_ihc)
unique(usedata$her2_status_by_ihc)
unique(usedata$subtype)

usedata$er_status_by_ihc<-factor(usedata$er_status_by_ihc)
usedata$pr_status_by_ihc<- factor(usedata$pr_status_by_ihc)
usedata$her2_status_by_ihc <- factor(usedata$her2_status_by_ihc)
usedata$subtype <- factor(usedata$subtype)





# Unsupervised Learning


undata<-usedata[, c(-1358, -1359, -1360, -1361)]


# plot scatter matrix, "ggpairs": a function in the package "GGally"
pm <- ggpairs(usedata)
pm 
# pdf(file="E:\\2017 Fall\\BIOINF 527\\Group Project\\BRCA\\BRCA\\pairs.pdf", width=10, height=10)
# dev.off()
 

# PCA
pca <- prcomp(undata, center = TRUE, scale. = TRUE) 


# standard deviation of PCs
print(pca)
# PCA summary
summary(pca)
pca$rotation[, 1:3]

# pcloading <-  melt(pca$rotation[,1:3])
# pcloading
# barplot <- ggplot(data=pcloading) +
#   geom_bar(aes(x=Var1, y=value), stat="identity", position = "identity") + 
#   facet_wrap(~Var2, nrow=1)
# barplot


pcpoints <- data.frame(usedata[, c(1358, 1359, 1360, 1361)],pca$x)


p1 <-qplot(x=PC1, y=PC2, data=pcpoints, colour=er_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "ER status") + 
  ggtitle("a")
p2 <- qplot(x=PC1, y=PC3, data=pcpoints, colour=er_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "ER status") + 
  ggtitle("b")
p3<-qplot(x=PC2, y=PC3, data=pcpoints, colour=er_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "ER status") + 
  ggtitle("c")

# grid.arrange(p1, p2, p3, nrow=1, ncol=3)


p4 <-qplot(x=PC1, y=PC2, data=pcpoints, colour=pr_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "PR status") + 
  ggtitle("d")
p5 <- qplot(x=PC1, y=PC3, data=pcpoints, colour=pr_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "PR status") + 
  ggtitle("e")
p6<-qplot(x=PC2, y=PC3, data=pcpoints, colour=pr_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "PR status") + 
  ggtitle("f")

# grid.arrange(p4, p5, p6, nrow=1, ncol=3)


p7 <-qplot(x=PC1, y=PC2, data=pcpoints, colour=her2_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "HER2 status") + 
  ggtitle("g")
p8 <- qplot(x=PC1, y=PC3, data=pcpoints, colour=her2_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "HER2 status") + 
  ggtitle("h")
p9<-qplot(x=PC2, y=PC3, data=pcpoints, colour=her2_status_by_ihc) +
  theme(legend.position="right") + theme_bw() + scale_color_discrete(name = "HER2 status") + 
  ggtitle("i")

# grid.arrange(p7, p8, p9, nrow=1, ncol=3)



p10 <-qplot(x=PC1, y=PC2, data=pcpoints, colour=subtype) +
  theme(legend.position="right") + theme_bw() + 
  ggtitle("j") + scale_color_discrete(name = "Combinations")
p11 <- qplot(x=PC1, y=PC3, data=pcpoints, colour=subtype) +
  theme(legend.position="right") + theme_bw() + 
  ggtitle("k") + scale_color_discrete(name = "Combinations")
p12 <-qplot(x=PC2, y=PC3, data=pcpoints, colour=subtype) +
  theme(legend.position="right") + theme_bw() + 
  ggtitle("l") + scale_color_discrete(name = "Combinations")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow=4, ncol=3)



# Hierarchical clustering
# cl.single <- hclust(dist(undata), method="single")
# plot(cl.single, xlab = "Genes")
# 
# cl.average <- hclust(dist(undata), method="average")
# plot(cl.average, xlab = "Genes")
# 
# cl.complete <- hclust(dist(undata))
# plot(cl.complete, xlab = "Genes")
# 
# heatmap.2(as.matrix(undata), trace="none", col=greenred(10))

# generate train dataset.
set.seed(12345)
trainn<-round(dim(usedata)[1]*0.8)
train.idex<-sample(1:dim(usedata)[1], trainn, replace=FALSE)

## Logistic Ridge Regression
# subtype:
ridgepath <- glmnet(x= as.matrix(undata[train.idex, ]), y=usedata$subtype[train.idex], family = "multinomial",alpha=0)
# summary(ridgepath)
# plot(ridgepath)
ridgecv <- cv.glmnet(x= as.matrix(undata[train.idex, ]), y=usedata$subtype[train.idex],  family = "multinomial",alpha=0)
lambda <- ridgecv$lambda.min
# plot(ridgecv)
# predict(ridgecv,s=lambda,type="coefficient")
classhat <- factor(predict(ridgecv, new = as.matrix(undata[-train.idex, ]), s=lambda, type = "class"))
# table(classhat, usedata$subtype[-train.idex])
A <- confusionMatrix(classhat, usedata$subtype[-train.idex])
# xtable(A$byClass[, c(1, 2, 3, 4, 5, 11)])
dat <- data.frame(class= usedata$subtype[-train.idex],classhat=classhat)
# p_LR_subtype <- ggplot(data=dat) + 
#   geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill") +
#   labs(title="Ridge regression fitting -- Subtype", x = "Subtype", y = "Proportion") + theme_bw() + 
#   scale_fill_discrete(name = "Subtype") + scale_color_discrete(name = "Subtype")
  
# separate:
# ER
ridgecv <- cv.glmnet(x= as.matrix(undata[train.idex, ]), y=usedata$er_status_by_ihc[train.idex],  family = "multinomial",alpha=0)
lambda <- ridgecv$lambda.min
classhat_er <- factor(predict(ridgecv, new = as.matrix(undata[-train.idex, ]), s=lambda, type = "class"))
B_er <- confusionMatrix(classhat_er, usedata$er_status_by_ihc[-train.idex])
dat_er <- data.frame(class= usedata$er_status_by_ihc[-train.idex],classhat=classhat_er)
p_LR_er <- ggplot(data=dat_er) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill",width = 0.5) +
  labs(title="a\nRidge regression fitting -- ER prediction", x = "ER status", y = "Proportion") + theme_bw() + 
  scale_fill_discrete(name = "ER status") + scale_color_discrete(name = "ER status")

# PR
ridgecv <- cv.glmnet(x= as.matrix(undata[train.idex, ]), y=usedata$pr_status_by_ihc[train.idex],  family = "multinomial",alpha=0)
lambda <- ridgecv$lambda.min
classhat_pr <- factor(predict(ridgecv, new = as.matrix(undata[-train.idex, ]), s=lambda, type = "class"))
B_pr <- confusionMatrix(classhat_pr, usedata$pr_status_by_ihc[-train.idex])
dat_pr <- data.frame(class= usedata$pr_status_by_ihc[-train.idex],classhat=classhat_pr)
p_LR_pr <- ggplot(data=dat_pr) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill",width = 0.5) +
  labs(title="\nRidge regression fitting -- PR prediction", x = "PR status", y = "Proportion") + theme_bw() + 
  scale_fill_discrete(name = "PR status") + scale_color_discrete(name = "PR status")

# her2
ridgecv <- cv.glmnet(x= as.matrix(undata[train.idex, ]), y=usedata$her2_status_by_ihc[train.idex],  family = "multinomial",alpha=0)
lambda <- ridgecv$lambda.min
classhat_her2 <- factor(predict(ridgecv, new = as.matrix(undata[-train.idex, ]), s=lambda, type = "class"))
B_her2 <- confusionMatrix(classhat_her2, usedata$her2_status_by_ihc[-train.idex])
dat_her2 <- data.frame(class= usedata$her2_status_by_ihc[-train.idex],classhat=classhat_her2)
p_LR_her2 <- ggplot(data=dat_her2) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill",width = 0.5) +
  labs(title="\nRidge regression fitting -- HER2 prediction", x = "HER2 status", y = "Proportion") + theme_bw() + 
  scale_fill_discrete(name = "HER2 status") + scale_color_discrete(name = "HER2 status")

dat_combine <- data.frame(class = usedata$subtype[-train.idex], 
                          ER = dat_er$classhat, PR = dat_pr$classhat, HER2 = dat_her2$classhat)
dat_combine$prediction <- apply(dat_combine[, 2:4], 1, subtypeName)

p_RL_comb <- ggplot(data = dat_combine) +
  geom_bar(aes(x = class, color = prediction, fill = prediction), alpha=0.5, position = "fill") + 
  theme_bw() + labs(title = "\nRidge Regression Fitting -- Combination", x = "Combination Status", y = "Proportion") +
  scale_fill_discrete(name = "Combinations") + scale_color_discrete(name = "Combinations")

grid.arrange(p_LR_er, p_LR_pr, p_LR_her2, p_RL_comb, nrow = 2, ncol = 2)
confusionMatrix(dat_combine$prediction, dat_combine$class)

## Random Forest
p <- ncol(usedata)
N <- nrow(usedata)

## subtype:
set.seed(123)
colnames(usedata)[1:1357] <- paste("Gene", 1:1357, sep = "")
rf.ff <- randomForest(formula=subtype~.-er_status_by_ihc-pr_status_by_ihc-her2_status_by_ihc, data=usedata[train.idex, ], importance=TRUE)
yhat.rf_subtype <- predict(rf.ff, newdata=usedata[-train.idex, ],  ntree=500,  type="class")
# table(yhat.rf, usedata$subtype[-train.idex])
C_subtype <- confusionMatrix(yhat.rf_subtype, usedata$subtype[-train.idex])
# xtable(B$byClass[, c(1, 2, 3, 4, 5, 11)])
datrf_subtype <- data.frame(class= usedata$subtype[-train.idex],classhat=yhat.rf_subtype)
# p_RF_subtype <- ggplot(data=datrf_subtype) + 
#   geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill") +
#   labs(title="A: Random Forest fitting -- subtype", x = "Subtype", y = "Proportion") + theme_bw() +
#   scale_color_discrete(name = "Subtype") + scale_fill_discrete(name = "Subtype")

## ER:
set.seed(234)
rf.ff <- randomForest(formula=er_status_by_ihc~.-subtype-pr_status_by_ihc-her2_status_by_ihc, data=usedata[train.idex, ], importance=TRUE)
yhat.rf_er <- predict(rf.ff, newdata=usedata[-train.idex, ],  ntree=500,  type="class")
C_er <- confusionMatrix(yhat.rf_er, usedata$er_status_by_ihc[-train.idex])
datrf_er <- data.frame(class= usedata$er_status_by_ihc[-train.idex],classhat=yhat.rf_er)
p_RF_er <- ggplot(data=datrf_er) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill",width=0.5) +
  labs(title="b\nRandom Forest fitting -- ER prediction", x = "ER status", y = "Proportion") + theme_bw() +
  scale_color_discrete(name = "ER status") + scale_fill_discrete(name = "ER status")

# PR:
set.seed(345)
rf.ff <- randomForest(formula=pr_status_by_ihc~.-subtype-er_status_by_ihc-her2_status_by_ihc, data=usedata[train.idex, ], importance=TRUE)
yhat.rf_pr <- predict(rf.ff, newdata=usedata[-train.idex, ],  ntree=500,  type="class")
C_pr <- confusionMatrix(yhat.rf_pr, usedata$pr_status_by_ihc[-train.idex])
datrf_pr <- data.frame(class= usedata$pr_status_by_ihc[-train.idex],classhat=yhat.rf_pr)
p_RF_pr <- ggplot(data=datrf_pr) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill", width = 0.5) +
  labs(title="\nRandom Forest fitting -- PR prediction", x = "PR status", y = "Proportion") + theme_bw() +
  scale_color_discrete(name = "PR status") + scale_fill_discrete(name = "PR status")

# HER2
set.seed(456)
rf.ff <- randomForest(formula=her2_status_by_ihc~.-subtype-er_status_by_ihc-pr_status_by_ihc, data=usedata[train.idex, ], importance=TRUE)
yhat.rf_her2 <- predict(rf.ff, newdata=usedata[-train.idex, ],  ntree=500,  type="class")
C_her2 <- confusionMatrix(yhat.rf_her2, usedata$her2_status_by_ihc[-train.idex])
datrf_her2 <- data.frame(class= usedata$her2_status_by_ihc[-train.idex],classhat=yhat.rf_her2)
p_RF_her2 <- ggplot(data=datrf_her2) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill", width = 0.5) +
  labs(title="\nRandom Forest fitting -- HER2 prediction", x = "HER2 status", y = "Proportion") + theme_bw() +
  scale_color_discrete(name = "HER2 status") + scale_fill_discrete(name = "HER2 status")

dat_combine <- data.frame(class = usedata$subtype[-train.idex], 
                          ER = datrf_er$classhat, PR = datrf_pr$classhat, HER2 = datrf_her2$classhat)
dat_combine$prediction <- apply(dat_combine[, 2:4], 1, subtypeName)

p_RF_comb <- ggplot(data = dat_combine) +
  geom_bar(aes(x = class, color = prediction, fill = prediction), alpha=0.5, position = "fill") + 
  theme_bw() + labs(title = "\nRandom Forest fitting -- Combination", x = "Combination Status", y = "Proportion") +
  scale_fill_discrete(name = "Combinations") + scale_color_discrete(name = "Combinations")

grid.arrange(p_RF_er, p_RF_pr, p_RF_her2, p_RF_comb, nrow = 2, ncol = 2)
confusionMatrix(dat_combine$prediction, dat_combine$class)

## boosting
# subtype:
set.seed(1234)
boost.ff <- gbm(formula= subtype~.-er_status_by_ihc-pr_status_by_ihc-her2_status_by_ihc, data=usedata[train.idex, ], distribution="multinomial",n.trees=500, interaction.depth=1)
# summary(boost.ff)
yhat.boost <- predict(boost.ff, newdata=usedata[-train.idex, ],  n.trees=500, type="response")
yhat.boost <- factor(colnames(yhat.boost)[apply(yhat.boost, 1, which.max)])
# table(yhat.boost, usedata$subtype[-train.idex])
D_subtype <-confusionMatrix(yhat.boost, usedata$subtype[-train.idex])
# xtable(C$byClass[, c(1, 2, 3, 4, 5, 11)])
datbt_subtype <- data.frame(class= usedata$subtype[-train.idex],classhat=yhat.boost)
# p_BT_subtype <- ggplot(data=datbt_subtype) + 
#   geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill") +
#   labs(title="A: Boosting fitting -- Subtype", x = "Subtype", y = "Proportion") + theme_bw() +
#   scale_color_discrete(name = "Subtype") + scale_fill_discrete(name = "Subtype")

# ER:
set.seed(2345)
boost.ff <- gbm(formula=er_status_by_ihc~.-subtype-pr_status_by_ihc-her2_status_by_ihc, data=usedata[train.idex, ], distribution="multinomial",n.trees=500, interaction.depth=1)
yhat.boost <- predict(boost.ff, newdata=usedata[-train.idex, ],  n.trees=500, type="response")
yhat.boost <- factor(colnames(yhat.boost)[apply(yhat.boost, 1, which.max)])
D_er <- confusionMatrix(yhat.boost, usedata$er_status_by_ihc[-train.idex])
datbt_er <- data.frame(class= usedata$er_status_by_ihc[-train.idex],classhat=yhat.boost)
p_BT_er <- ggplot(data=datbt_er) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill",width = 0.5) +
  labs(title="c\nBoosting fitting -- ER prediction", x = "ER status", y = "Proportion") + theme_bw() +
  scale_color_discrete(name = "ER status") + scale_fill_discrete(name = "ER status")

# PR:
set.seed(3456)
boost.ff <- gbm(formula=pr_status_by_ihc~.-subtype-er_status_by_ihc-her2_status_by_ihc, data=usedata[train.idex, ], distribution="multinomial",n.trees=500, interaction.depth=1)
yhat.boost <- predict(boost.ff, newdata=usedata[-train.idex, ],  n.trees=500, type="response")
yhat.boost <- factor(colnames(yhat.boost)[apply(yhat.boost, 1, which.max)])
D_pr <- confusionMatrix(yhat.boost, usedata$pr_status_by_ihc[-train.idex])
datbt_pr <- data.frame(class= usedata$pr_status_by_ihc[-train.idex],classhat=yhat.boost)
p_BT_pr <- ggplot(data=datbt_pr) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill",width = 0.5) +
  labs(title="\nBoosting fitting -- PR prediction", x = "PR status", y = "Proportion") + theme_bw() +
  scale_color_discrete(name = "PR status") + scale_fill_discrete(name = "PR status")

# HER2
set.seed(4567)
boost.ff <- gbm(formula=her2_status_by_ihc~.-subtype-er_status_by_ihc-pr_status_by_ihc, data=usedata[train.idex, ], distribution="multinomial",n.trees=500, interaction.depth=1)
yhat.boost <- predict(boost.ff, newdata=usedata[-train.idex, ],  n.trees=500, type="response")
yhat.boost <- factor(colnames(yhat.boost)[apply(yhat.boost, 1, which.max)])
D_her2 <- confusionMatrix(yhat.boost, usedata$her2_status_by_ihc[-train.idex])
datbt_her2 <- data.frame(class= usedata$her2_status_by_ihc[-train.idex],classhat=yhat.boost)
p_BT_her2 <- ggplot(data=datbt_her2) + 
  geom_bar(aes(x=class,color=classhat,fill=classhat),alpha=0.3,position = "fill",width = 0.5) +
  labs(title="\nBoosting fitting -- HER2 prediction", x = "HER2 status", y = "Proportion") + theme_bw() +
  scale_color_discrete(name = "HER2 status") + scale_fill_discrete(name = "HER2 status")

dat_combine <- data.frame(class = usedata$subtype[-train.idex], 
                          ER = datbt_er$classhat, PR = datbt_pr$classhat, HER2 = datbt_her2$classhat)
dat_combine$prediction <- apply(dat_combine[, 2:4], 1, subtypeName)

p_BT_comb <- ggplot(data = dat_combine) +
  geom_bar(aes(x = class, color = prediction, fill = prediction), alpha=0.5, position = "fill") + 
  theme_bw() + labs(title = "\nBoosting fitting -- Combination", x = "Combination Status", y = "Proportion") +
  scale_fill_discrete(name = "Combinations") + scale_color_discrete(name = "Combinations")
  
grid.arrange(p_BT_er, p_BT_pr, p_BT_her2, p_BT_comb, nrow = 2, ncol = 2)
confusionMatrix(dat_combine$prediction, dat_combine$class)
