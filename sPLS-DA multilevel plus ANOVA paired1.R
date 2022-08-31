
##=====Follicular fluid=====Metaphase II samples
##http://mixomics.org/graphics/sample-plots/plotindiv/

#setwd("F:/PhD/work/01 TP STUDY/TP metabolomica") 

#====INSTALL PACKAGES================================================================================

# Lista de paquetes de funciones a instalar
.packages = c("BiocManager","devtools","ggplot2","pcadapt","outliers","igraph",
              "rgl","graphics","reshape2","dplyr","ggpubr","remotes",
              "FactoMineR", "factoextra","corrplot","ggpubr","fpc", "NbClust","mixOmics", 
              "phyloseq","lme4","nlme","car","plotly","RadaFDR","GOplot",
              "GOsummaries","randomcoloR","parallel","doParallel","afex","ggbeeswarm",
              "emmeans","psych","caret","pheatmap","circlize","ComplexHeatmap")
#"RFunrichWebService",
# Instala los paquetes sin? los tienes instalados
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  install_github("vqv/ggbiplot")
  install_github("fxia22/RadaFDR")
  BiocManager::install(c("mixOmics","phyloseq","RFunrichWebService","GOsummaries"))
  BiocManager::install("RFunrichWebService",version = "3.8")
  BiocManager::install("RFmarkerDetector")
  BiocManager::install("ComplexHeatmap")
}

# Carga los paquetes sin? los tienes cargados
lapply(.packages, require, character.only=TRUE)

#======END INSTALL PACKAGES============================================================================

#=====Loading dataset================================================================================== 
Big.table.ipp <- as.data.frame(readxl::read_excel("./input data/Testosterone_metabolites_Preliminary_data (4ipp2).xlsx")) # Open the data
rownames(Big.table.ipp) <- Big.table.ipp[,1]
class(Big.table.ipp)
main_table.ipp <- as.data.frame(Big.table.ipp[,2:ncol(Big.table.ipp)])

Big.table <- as.data.frame(readxl::read_excel("./input data/normalized_data_testosterone_220221_ipp.xlsx")) # Open the data
rownames(Big.table) <- Big.table[,1]
class(Big.table)


main_table <- as.data.frame(Big.table[,2:ncol(Big.table)])

ClinicalData <- as.data.frame(readxl::read_excel("./input data/Clinical Data_Barb.xlsx")) # Open the data
row.names(ClinicalData) <- ClinicalData$Patient
head(ClinicalData)
ClinicalData1 <- ClinicalData[,2:ncol(ClinicalData)]

Annotations <- as.data.frame(readxl::read_excel("./Annotations.xlsx")) # Open the data
row.names(Annotations) <- Annotations$sample
head(Annotations)

#=== Checking normalization of the data =================================
colors = c(rep("green",29),rep("lightblue",29),rep("red",29))
boxplot(main_table, las=2)
boxplot(main_table,las = 2, col=colors)+
  legend("topright", legend = c("Basal","Low","restored"), 
         col = c("red", "green", "blue"), 
         pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)


gg<- melt(main_table)
colnames(gg)[1] <- "sample"
gg1 <- plyr::join_all(list(gg,Annotations),by="sample")
gg1 <- na.omit(gg1)

ggplot(gg1, aes(x = Patient, y = value, fill = Group)) +
  geom_boxplot(alpha=0.7)+theme_bw()

colors1 = c(rep("darkgreen",29),rep("blue",29),rep("red",29))


p0 <- ggplot(gg1, aes(x = Patient, y = value, fill = Group)) +
  geom_boxplot(alpha=0.7)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title="Normalization: Original normalization",
       x = "patients", y= "Intensity (Log2 scaled)")
p0

p1 <- ggplot(gg1, aes(x=sample, y=value, fill = Group)) +
  geom_boxplot(alpha=0.7)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title="Normalization: Original normalization",
       x = "samples", y= "Intensity (Log2 scaled)")
p1

p2 <- ggplot(gg1, aes(x=value, color=sample)) +
  geom_density()+theme_bw()+theme(legend.position="none")+
  scale_x_continuous(limits = c(-12, 12))+
  labs(title="Normalization: Original normalization",
       x = "Intensity (Log2 scaled)")
p2

gg0<- melt(main_table.ipp)
colnames(gg0)[1] <- "sample"
gg01 <- plyr::join_all(list(gg0,Annotations),by="sample")
gg01 <- na.omit(gg01)

p01 <- ggplot(gg01, aes(x = Patient, y = value, fill = Group)) +
  geom_boxplot(alpha=0.7)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title="Normalization: Log2 + sample median substraction",
       x = "patients", y= "Intensity (Log2 scaled)")
p01

# Multiple plot function http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/#:~:text=multiplot%20function%20This%20is%20the%20definition%20of%20multiplot.,a%20list%20of%20plot%20objects%20passed%20to%20plotlist.
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(p0,p1, p2,p01, cols=2)


#====sPLS-DA multilevel and PCA ============MULTIVARIATE ANALYSIS=======mixomics package==


X <- as.data.frame(t(main_table))
Y <- as.factor(Annotations$Group)

# Patient indicates the repeated measurements
# setup the design matrix by indicating the repeated measurements
design <- data.frame(Patient = Annotations$Patient)


#----PCA tuning----------------------------------------

tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
result <- mixOmics::pca(X, ncomp = 3, center = TRUE, scale = FALSE)
result

#--PCA--classic
pca.classic.TP <- mixOmics::pca(X, ncomp = 3, scale = F, center = T)
pca.classic.TP

mixOmics::plotIndiv(pca.classic.TP, ind.names = F, legend = T,ellipse = T,cex = 1,
                    group = Y,title = 'PCA (classic)')
#----PCA--multilevel-----------------------------------
pca.multilevel.TP <- mixOmics::pca(X, ncomp = 3, scale = F, center = T,multilevel = design)
pca.multilevel.TP

mixOmics::plotIndiv(pca.multilevel.TP, ind.names = F, legend = T,ellipse = T,cex = 1,
                    group = Y,title = 'PCA (multilevel)')


#----PLS-DA multilevel model----------------

TP.plsda.multilevel <- mixOmics::plsda(X, Y=Y,
                             multilevel = design, 
                             ncomp = 3)


mixOmics::plotIndiv(TP.plsda.multilevel, ind.names = F, ellipse =T, ellipse.level=.95,star = F, legend=T,
          style = 'ggplot2',cex=1.5,title = "PLS-DA multilevel")

# plotIndiv(TP.plsda.multilevel, ind.names = Y, 
#           style = '3d')
# plotIndiv(TP.plsda.multilevel, ind.names = F, ellipse =T,ellipse.level=.95,star = F, legend=T,
#           style = '3d',cex=0.5,title = "PLS-DA")

col.ID <- randomcoloR::distinctColorPalette(k=29)[as.integer(Annotations$Patient)]

col.groups <- color.mixo(1:3)[as.integer(Annotations$Groups.num)]

mixOmics::cim(TP.plsda.multilevel,comp=c(1:2),
    row.sideColors = cbind(col.groups, col.ID),row.cex = 0.8,clust.method=c("complete","complete"),
    row.names = paste(Annotations$Group, Annotations$Patient, sep = "_"),
    col.names = FALSE, legend=list(legend = c(levels(as.factor(Annotations$Groups.num))), 
                                   col = color.mixo(1:3),
                                   title = "Time point", cex = 0.8))
mixOmics::cim(TP.plsda.multilevel,comp=c(1:2),transpose = T,margins = c(2,23),zoom=F,keysize = c(1,1),
    row.sideColors = col.groups,row.cex = 0.7,col.cex = 0.7,clust.method=c("complete","complete"),
    row.names = Annotations$Group,
    #title = "sPLS-DA multilevel (Comp 1-2)",
    #cut.tree = c(1,1),
    #row.names=F,
    col.names = T, legend=list(legend = c(levels(as.factor(Annotations$Groups.num))), 
                               col = color.mixo(1:3),
                               title = "Time point", cex = 0.8))


# The plot Loading function displays the loading weights, 
mixOmics::plotLoadings(TP.plsda.multilevel, comp = 1, title = 'Loadings on comp 1', size.name = 0.6,
             contrib = 'max', method = 'mean',ndisplay=30)#,ndisplay=50)

# saving matrixs

loading.matrix.X <- TP.plsda.multilevel$loadings$X
#write.csv(loading.matrix.X,"loading.matrix.X.plsda.csv")


#------sPLS-DA multilevel model -------before Tuning----------------------------

TP.splsda.multilevel <- mixOmics::splsda(X, Y = Y,
                               multilevel = design, 
                               ncomp = 5)



plotIndiv(TP.splsda.multilevel, ind.names = F, ellipse =T,ellipse.level=.95,star = F, legend=T,
          style = 'ggplot2',cex=1,title = "sPLS-DA multilevel")

# plotIndiv(TP.splsda.multilevel, ind.names = Y, style = '3d')


#----Tuning parameters and numerical outputs___sPLS-DA multilevel

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat

## 1 - The number of components to retain

cl <- makePSOCKcluster(6)                                    #  nCores=6 ??? For parallelism
registerDoParallel(cl)
MyPerf.plsda1 <- perf(TP.splsda.multilevel, validation = "Mfold", folds = 5, dist="all",         # folds :: https://machinelearningmastery.com/k-fold-cross-validation/#:~:text=Cross%2Dvalidation%20is%20a%20resampling,k%2Dfold%20cross%2Dvalidation. 
                      progressBar = T, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda1, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda1$error.rate
MyPerf.plsda1$error.rate.class
MyPerf.plsda1$choice.ncomp
#MyPerf.plsda1$predict         # Prediction values for each component

## 2 - The number of variables keepX

list.keepX <- seq(1, 717,10)
head(list.keepX)


cl <- makePSOCKcluster(6)                                    #  nCores=6 ??? For parallelism
registerDoParallel(cl)

tune.splsda.srbct <- tune.splsda(X=X, Y=Y, ncomp = 3, # ncomp = 3 as suggested in the previous step
                                 multilevel = design,
                                 validation = 'Mfold',
                                 folds = 5, dist = 'max.dist',
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50,
                                 signif.threshold=0.01,
                                 progressBar = T)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct)


##------ sPLS-DA after tunning-----------------
#keepX = c(1, 21)
TP.splsda.multilevel <- mixOmics::splsda(X, Y = Y,
                               multilevel = design, 
                               ncomp = 2, 
                               keepX = c(1, 100))




plotIndiv(TP.splsda.multilevel, ind.names = F, ellipse =T,ellipse.level=.95,star = F, legend=T,
          style = 'ggplot2',cex=1,title = "sPLS-DA multilevel")

# ---performance of our final sPLS-DA model-------

perf.obj<-perf(TP.splsda.multilevel, 
               dist = "max.dist",
               validation = "Mfold",
               folds = 5, nrepeat =50, auc = T, progressBar = TRUE)

perf.obj$error.rate
perf.obj$error.rate.class
plot(perf.obj, dist = "max.dist",measure = "BER",
     xlab = NULL, ylab = NULL, overlay="dist",legend.position="vertical",sd = TRUE)

#---ploting heatmap-----------------------------

col.ID <- randomcoloR::distinctColorPalette(k=29)[as.integer(Annotations$Patient)]

col.groups <- color.mixo(1:3)[as.integer(Annotations$Groups.num)]

cim(TP.splsda.multilevel,comp=c(1:2),
    row.sideColors = cbind(col.groups, col.ID),row.cex = 0.8,
    clust.method=c("complete","complete"),
    row.names = paste(Annotations$Group, Annotations$Patient, sep = "_"),
    col.names = FALSE, legend=list(legend = c(levels(as.factor(Annotations$Groups.num))), 
                                   col = color.mixo(1:3),
                                   title = "Time point", cex = 0.8))
cim.results <- cim(TP.splsda.multilevel,comp=c(1:2),transpose = T,
    margins = c(2,23),zoom=F,keysize = c(1,1),
    row.sideColors = col.groups,row.cex = 0.9,
    col.cex = 0.7,cluster = "column",
    dist.method = c("euclidean","euclidean"),
    clust.method= c("average","average"),
    row.names = Annotations$Group,
    title = "sPLS-DA multilevel (Comp 1-2)",
    #cut.tree = c(1,1),
    #row.names=F,
    col.names = T, legend=list(legend =  c("Basal","Low","Restored"), 
                               col = color.mixo(1:3),
                               title = "Time point", cex = 0.8))
cim.results

cim.matrix <- cim.results$mat
#write.csv(cim.matrix,"cim.matrix_sPLS-DA_C1.1_C2.100.csv")

cim.matrix1 <- t(cim.matrix)

col3<- colorRampPalette(c("dodgerblue4","dodgerblue3", "white","coral2","brown3"))(100)

heatmap.func <- function(m,plot.name="Title", cluster_cols=T){
  
  fontsize = 7
  hm<-pheatmap(m, main=plot.name,cutree_cols = 1,
               cutree_rows = 4, color = col3,show_rownames = F, show_colnames = F,
               cluster_cols = F, cluster_rows = T,border_color="white",
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
               clustering_distance_cols = "euclidean", 
               gaps_col = c(29,58),fontsize = fontsize, legend = T, 
               annotation_legend = T,fontsize_number = fontsize,
               #filename = "sPLS-DA.pdf", cellwidth = 10, cellheight = 12,
               treeheight_row = 0, breaks = round(seq(-2,2,length=101 ),2)) #, filename = "test.pdf", cellwidth = 10, cellheight = 12
  #dev.off()
  print(hm)
  return(hm)
}

heatmap.obj <- heatmap.func(cim.matrix1,
             plot.name = "cim-heatmap (101 metabolite from sPLS-DA multilevel)",
             cluster_cols=T)

#dev.off()

dendog <- as.dendrogram(heatmap.obj$tree_row)
plot(dendog)
dendog1 <- as.hclust(heatmap.obj$tree_row)
clusters <-as.data.frame(cutree(dendog1, h=27))
colnames(clusters) <- "cluster"
clusters$compound <- row.names(clusters)
as.factor(clusters$cluster)


#----Ploting Loading function displays the loading weights, ---
 # The plotLoading function displays the loading weights, where colors indicate 
 # the class for which the selected variable has a maximal mean value

plotLoadings(TP.splsda.multilevel, comp = 2, title = 'Loadings on comp 2', size.name = 0.6,
             contrib = 'max', method = 'mean',ndisplay=21)#,ndisplay=50)

# saving matrixs

loading.matrix.X <- as.data.frame(TP.splsda.multilevel$loadings$X)
#write.csv(loading.matrix.X,"loading.matrix.X_sPLS-DA_C1.1_C2.100.csv")

vip.feature <- as.data.frame(vip(TP.splsda.multilevel))
colnames(vip.feature) <- paste("VIP_", colnames(vip.feature),sep = "")
#write.csv(vip.feature,"VIP.matrix.X_sPLS-DA_C1.1_C2.100.csv")

vip.feature$compound <- rownames(vip.feature)
loading.matrix.X$compound <- rownames(loading.matrix.X)

sPLS.DA_result <- plyr::join_all(list(loading.matrix.X, vip.feature,clusters),by="compound")
rownames(sPLS.DA_result) <- sPLS.DA_result$compound

#write.csv(sPLS.DA_result,"sPLS.DA_result_C1.1_C2.100.csv")

#=======Integrative analysis============================================

design=design

X1 = X
Y1 = ClinicalData1

list.keepY <- seq(1, 38,1)
head(list.keepX)

#---Tuning clinical data--------------------

tune.spls.CD<-tune.splslevel(X1, 
                             Y1,
                             multilevel=design, 
                             ncomp = 3, 
                             test.keepY = list.keepY,
                             test.keepX = list.keepX,
                             already.tested.X = c(1,100),
                             already.tested.Y = c(1,38),
                             mode = "canonical")

class(tune.spls.CD$cor.value)
colnames(tune.spls.CD$cor.value)
max(tune.spls.CD$cor.value)
which(tune.spls.CD$cor.value == max(tune.spls.CD$cor.value), arr.ind=TRUE)

error <- tune.spls.CD$error.rate
error
ncomp <- tune.spls.CD$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.spls.CD$choice.keepX  # optimal number of variables to select
select.keepX

plot(tune.spls.CD)


#-------spls correlation between metabolomics and clinical data------

clinic.spls.multilevel <- spls(X = X1,
                                 Y = Y1,
                                 multilevel = design,
                                 ncomp = 2,
                                 keepX = c(1, 50), 
                                 keepY = c(38, 38),
                                 mode = 'canonical')
plotVar(clinic.spls.multilevel, comp = 1:2, var.names = TRUE, cex = c(2,2)) 

tune.spls <- perf(clinic.spls.multilevel, 
                  validation = "Mfold", 
                  folds = 10, 
                  progressBar = T, 
                  nrepeat = 50)
plot(tune.spls$Q2.total)
abline(h = 0.0975)

stim.col <- c("darkblue", "green4","red3")

# showing the XY variables.

cim.spls.multilevel <- cim(clinic.spls.multilevel, mapping="XY", 
                      comp = c(1,2), margins = c(6,23))

#dev.off()
cim.spls.matrix <- cim.spls.multilevel$mat.cor

#write.csv(cim.spls.matrix,"cim.matrix_spls.multilevel_C1.1_C2.100.csv")

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(cim.spls.matrix,cutree_cols = 4,
         cutree_rows = 4,show_rownames = T, 
         show_colnames = T,
         cluster_cols = T, cluster_rows = T,
         border_color="white",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         clustering_distance_cols = "euclidean", 
         fontsize = 7, #filename = "test.pdf", cellwidth = 10, cellheight = 12,
         fontsize_number = 7,#clustering_callback = callback,
         #filename = "sPLS-DA.pdf", cellwidth = 10, cellheight = 12,
         treeheight_row = 10,treeheight_col = 10, 
         breaks = round(seq(-0.8,0.8,length=101 ),2)) #, filename = "test.pdf", cellwidth = 10, cellheight = 12
dev.off()
##====ANOVA paired=======================================================
library(afex)
library(ggbeeswarm)
library(emmeans)



table.aov <- cbind(t(main_table.ipp[1,]),Annotations)

# Within.aov.2<-aov_car(DV ~ Within_Cond*Within_Time + Error(Subject/Within_Cond*Within_Time), data=Within_Data,
#                       observed = "Within_Time") #this corrects GES if a variable is observed and not manipulat
#aov_car(dv ~ Group + Error(Patient/Group), data=tab) 

#===ANOVA paired

# m = matrix with samples as columns and proteins as rows
# annotation.table = Annotation matrix
# groups = name of the column that contain the groups that should be compared (factor)  

# z-score normalization by rows to plot a heatmap

z.score.func <- function(x){
  (x - mean(x))/sd(x)
}


m=main_table
annotation.table = Annotations
groups = "Group"
id_ind = "Patient"
i=1

DEP_prot <- function(m, annotation.table, groups, id_ind, plot.var=T) {
  
  groups1 <- annotation.table %>% dplyr::select(all_of(groups))
  groups1 <- as.factor(groups1[,1])
  n = length(levels(groups1))  # number of groups en la variable "groups" = 3
  k = 2
  
  comb <- factorial(n)/(factorial(k)*factorial(n-k))  # number of combinations
  comb
  
  #===
  
  p.value <- vector(mode = "numeric")
  pv <- matrix(nrow=nrow(m),ncol = comb, dimnames = list(row.names(m),c("pv(Basal-Low)","pv(Basal-Restor.)","pv(Low-Restor.)")))  ######
  diff <- matrix(nrow=nrow(m),ncol = comb,dimnames = list(row.names(m),c("Log2FC(Basal-Low)","Log2FC(Basal-Restor.)","Log2FC(Low-Restor.)")))
  means <- matrix(nrow=nrow(m),ncol = comb, dimnames = list(row.names(m),c("mean(Basal)","mean(Low)","mean(Restored)")))
  # lwr <- matrix(nrow=nrow(m),ncol = comb,dimnames = list(row.names(m),c("pv(Basal-Low)","pv(Basal-Restor.)","pv(Low-Restor.)")))
  # upr <- matrix(nrow=nrow(m),ncol = comb,dimnames = list(row.names(m),c("pv(Basal-Low)","pv(Basal-Restor.)","pv(Low-Restor.)")))
  
  dir.create("Boxplots")
  m1 = t(m)
  groups2 <- annotation.table %>% dplyr::select(all_of(c(groups,id_ind,"sample")))
  
  for (i in 1:ncol(m1)) {
    
    prot.data <- as.data.frame(m1[,i])
    
    prot.data$sample <- rownames(prot.data)
    
    data <- plyr::join_all(list(prot.data,groups2),by="sample")
    row.names(data) <- data$sample
    colnames(data)[1] <- "analyte"
    data <- data[,-2]
    data$Group <- as.factor(data$Group)
    #===ANOVA=======
    
    aov <- afex::aov_ez(id = id_ind, dv = 'analyte',data = data,
                        fun_aggregate = mean, within = groups,
                        return = afex_options("return_aov"))
    
    p.value[i] <- aov$anova_table$`Pr(>F)`
    
    
    
    #post hoc
    
    posthoc <- emmeans(aov,~ Group)
  
    
    posthoc1 <- as.data.frame(posthoc)
    means[i,] <- posthoc1$emmean
    
    pv.posthoc <- as.data.frame(pairs(posthoc, adjust="tukey"))
    
    pv[i,] <- pv.posthoc$p.value
    diff[i,] <- pv.posthoc$estimate
    
    posthoc1 <- as.data.frame(posthoc)
    
    if (plot.var == TRUE){
      
    
    #plots
    
    tuk.sig <- cut(pv[i,], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("p<.001", "p<.01", "p<.05","N/S"))

    png(filename = paste("Boxplots/","plot_",i,"_",".png",sep=""),width = 500, height = 480, units = "px", pointsize = 16,
        bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))


    #https://datavizpyr.com/how-to-make-boxplot-with-a-line-connecting-mean-values-in-r/#google_vignette
    df1_median <- data %>% 
      group_by(Group) %>% 
      summarize(average = median(analyte)) %>%
      ungroup()
    
     plot1 <- data %>% 
      ggplot(mapping = aes(x = Group, y = analyte, fill = Group)) + 
      geom_boxplot(alpha=0.5) +
      geom_beeswarm(priority='density',cex=1.5, alpha=0.5,show.legend = FALSE)+
      theme_bw(base_size = 16)+
      geom_point(data = df1_median, show.legend = FALSE,
                 mapping = aes(x = Group, y = average),
                 color="red", size=3) +
      geom_line(data = df1_median, colour="red",
                mapping = aes(x = Group, y = average, group=1))+
      labs(subtitle=rownames(m)[i],
           caption = paste("Tukey p.val: ",
                           "Basal-Low = ",tuk.sig[1],
                           "; Basal-Restor. = ",tuk.sig[2],
                           "; Low-Restor. = ",tuk.sig[3],sep=""),
           x = "Time point", y= "Intensity (Log2 scaled)")
    
    print(plot1)
    
    
    dev.off()
    
    
    
    }
  }  
  
  final.data <- cbind(means,diff,p.value,pv)
  row.names(final.data) <- row.names(m)
  
  final.data <- as.data.frame(final.data)
  final.data$q.value <- p.adjust(final.data$p.value, method="fdr")
  
  for (j in 7:(7+comb+1)){
    
    final.data$sig <- cut(final.data[,j], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
    
    colnames(final.data) <- c(colnames(final.data)[-c(length(colnames(final.data)))],paste("Sig_",colnames(final.data)[j]))
  }
  
  #final.data$Sig.Q.value <- cut(final.data$q.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  return(final.data)
}


DEP <- as.data.frame(DEP_prot(m=main_table,
                              annotation.table=Annotations,
                              groups="Group",
                              id_ind="Patient",
                              plot.var = F))


DEP1 <- cbind(Big.table,DEP)

#write.csv(DEP1, "DEP_fdr_ANOVA paired.csv")


DEP$compound <- rownames(DEP)

sPLS.DA_ANOVA <- plyr::join_all(list(DEP,sPLS.DA_result),by="compound")
#write.csv(sPLS.DA_ANOVA,"sPLS-DA plus ANOVAez.csv")

#=======Selecting the significant proteins from the sPLS-DA analysis=====

sPLS.DA_ANOVA.Sig <- subset(sPLS.DA_ANOVA,sPLS.DA_ANOVA$cluster!="NA")
sPLS.DA_ANOVA.Sig_mean <- sPLS.DA_ANOVA.Sig %>% dplyr::select(contains("mean"))  

sPLS.DA_ANOVA.Sig_Zmean <- as.matrix(sPLS.DA_ANOVA.Sig_mean[,1:3])

sPLS.DA_ANOVA.Sig_Zmean <- as.data.frame(t(apply(sPLS.DA_ANOVA.Sig_Zmean,1,z.score.func)))

sPLS.DA_ANOVA.Sig_Zmean$compound <- sPLS.DA_ANOVA.Sig$compound


# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#clustering
library(ComplexHeatmap)
set.seed(123)

cim.matrix1_A <- as.data.frame(cim.matrix1)
cim.matrix1_A$compound <- rownames(cim.matrix1)

cim.matrix1_B <- plyr::join_all(list(cim.matrix1_A,
                                     sPLS.DA_ANOVA.Sig_Zmean),by="compound")
mat1 <- as.matrix(cim.matrix1_B %>% select(contains("TP")))
mat2 <- as.matrix(cim.matrix1_B %>% select(contains("mean")))


library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "orange"))
col_fun(seq(-2, 2))

top_annotation_Hm1 = HeatmapAnnotation(`Testosterone level` = c(rep("Basal",29),
                                                            rep("Low",29),
                                                            rep("High",29)),
                                   col = list(`Testosterone level` = c("Basal" = "deepskyblue3",
                                                                       "Low" = "darkorange2",
                                                                       "High" = "grey")))

top_annotation_Hm2 = HeatmapAnnotation(`Testosterone level` = c(rep("Basal",1),
                                                                rep("Low",1),
                                                                rep("High",1)),
                                       col = list(`Testosterone level` = c("Basal" = "deepskyblue3",
                                                                           "Low" = "darkorange2",
                                                                           "High" = "grey")))


Hm1 <- Heatmap(mat1, name = "sPLS-DA",
               cluster_columns = F,
               row_km = 4, #column_split = rep(c("A", "B","C"), 29),
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "ward.D2",
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "ward.D2",
               width = unit(9, "cm"), height = unit(12, "cm"),
               top_annotation = top_annotation_Hm1,border = T,
               show_column_names = F)
Hm1

Hm2 <- Heatmap(mat2, name = "Log2 Intensity(scaled)",col = col_fun,
               cluster_columns = F,cluster_rows = F,border = T,
               width = unit(2, "cm"), height = unit(12, "cm"),
               top_annotation = top_annotation_Hm2,
               show_column_names = F)

Hm2


ht_list = Hm1 + Hm2

ht_opt(legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE)

draw(ht_list, ht_gap = unit(4, "mm"))

row_order_list <- row_order(ht_list)

write.infile(row_order_list,"row_order_list.txt")

###======Significant metabilites from ANOVA paired test========
DEP1$abs_Log2FC_AB <- abs(DEP1$`Log2FC(Basal-Low)`)
DEP1$abs_Log2FC_BC <- abs(DEP1$`Log2FC(Low-Restor.)`)


Sig.DEP1 <- subset(DEP1,DEP1$q.value<=0.05)
Sig.DEP1 <- subset(Sig.DEP1,Sig.DEP1$abs_Log2FC_AB>1)
Sig.DEP1 <- subset(Sig.DEP1,Sig.DEP1$abs_Log2FC_BC>1)

Sig.DEP2 <- Sig.DEP1 %>% dplyr::select(contains("TP"))

#write.csv(Sig.DEP2,"Sig.DEP2.csv")

# z-score normalization by rows to plot a heatmap

z.score.func <- function(x){
  (x - mean(x))/sd(x)
}


Sig.DEP2zscore <- as.data.frame(t(apply(Sig.DEP2,1, z.score.func)))

# heatmap

library(pheatmap)
pheatmap(Sig.DEP2zscore,show_rownames = F,cluster_cols = F)

annotation_rows <- sPLS.DA_ANOVA %>% dplyr::select(contains(c("q.v","VIP")))
rownames(annotation_rows) <- sPLS.DA_ANOVA$compound
annotation_rows <- subset(annotation_rows,annotation_rows$q.value<=0.05)

annotation_col <- as.data.frame(annotation.table[,c("Group")])
rownames(annotation_col)<-rownames(annotation.table)
colnames(annotation_col) <- "Time point"

ann_colors = list(`Time point`=c(Basal="deepskyblue3", Low="darkorange2", Restored="grey"))
ann_colorss = list(Sig.C1 = c(Non.Sig= "darkgrey"),
                   Sig.C2 = c(Down = "darkgreen", Non.Sig ="darkgrey",Up = "darkorange2"),
                   Sig.C3 = c(Down = "darkgreen", Non.Sig = "darkgrey"),
                   `Cluster PLS-Cox`= c(Clust.1 = "darkgoldenrod1", Clust.2 = "red", Clust.3 = "blue"),
                   `Disease stage`= c(stage.3 = "mediumorchid1",stage.4 = "darkolivegreen1"),
                   Lund = c(High.imm = "steelblue2", Normal = "lightsalmon", Pigmentation = "lightgreen", Proliferative ="lightgoldenrod4", nan="lightgrey"),
                   TCGA = c(immune ="violet", keratin ="turquoise1", MITF.low="yellow3", nan = "lightgrey"),
                   `Melanoma type` = c(ALM = "purple4", LMM = "seagreen3", NM = "sienna1", SSM = "yellow2", Unknownprimary ="violetred1", nan ="lightgrey"),
                   `BRAF status` = c(V600A = "orange", V600E = "maroon1", V600K = "paleturquoise2", WT = "limegreen", nan = "lightgrey"))
#____________

break1 <- round(seq(-2,2,length=101 ),2)
break2 <- round(seq(-6,6,length=101 ),2)
#==== OTRA FORMA DE GRAFICAR

col<- colorRampPalette(c("red4", "coral3","coral3","lightgoldenrodyellow","steelblue","steelblue", "royalblue4"))(100)
col01<- colorRampPalette(c("red4", "coral3","lightgoldenrodyellow","steelblue", "royalblue4"))(100)
col1 <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)
col2 <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
col3<- colorRampPalette(c("dodgerblue4","dodgerblue3", "lemonchiffon","coral2","brown3"))(100)

#=

heatmap.func <- function(m,plot.name="Title", cluster_cols=T){
  
  fontsize = 7
  hm<-pheatmap(m, main=plot.name,cutree_cols = 1,
               cutree_rows = 1, color = col3,show_rownames = T, show_colnames = T,
               cluster_cols = T, cluster_rows = T,
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
               gaps_col = c(29,58),fontsize = fontsize, legend = T, 
               annotation_legend = T,fontsize_number = fontsize,
               treeheight_row = 0, breaks = break1) #, filename = "test.pdf", cellwidth = 10, cellheight = 12
  print(hm)
}

heatmap.func(Sig.DEP2zscore,plot.name = "DEP (ANOVA paired), z-score",cluster_cols=T)

heatmap.func(Sig.DEP2,plot.name = "DEP (ANOVA paired)",cluster_cols=F)
heatmap.func(Sig.DEP2zscore,plot.name = "DEP (ANOVA paired), z-score",cluster_cols=F)


# only proteins thta change A-B and B-C

Sig.DEP1 <- subset(DEP1,DEP1$q.value<=0.05)
Sig.DEP1 <- subset(Sig.DEP1,Sig.DEP1$`pv(A-B)` <= 0.05)
Sig.DEP1 <- subset(Sig.DEP1,Sig.DEP1$`pv(B-C)` <= 0.05)

Sig.DEP2.a <- Sig.DEP1 %>% dplyr::select(contains("TP"))
Sig.DEP2.a_zscore <- as.data.frame(t(apply(Sig.DEP2.a,1, z.score.func)))

heatmap.func(Sig.DEP2.a,plot.name = "DEP A-B & B-C (posthoc ANOVA paired)",cluster_cols=T)
heatmap.func(Sig.DEP2.a,plot.name = "DEP A-B & B-C (posthoc ANOVA paired), z-score",cluster_cols=T)

heatmap.func(Sig.DEP2,plot.name = "DEP (ANOVA paired)",cluster_cols=F)
heatmap.func(Sig.DEP2.a_zscore,plot.name = "DEP (ANOVA paired), z-score",cluster_cols=F)

#== only metabolites that changed in sPLS-DA analysis===

sPLS.DA_ANOVA.1 <- subset(sPLS.DA_ANOVA,sPLS.DA_ANOVA$cluster %in% c(1:4))
metab.splsDA <- sPLS.DA_ANOVA.1$compound

Sig.DEP.splsDA <- subset(DEP1, row.names(DEP1) %in% metab.splsDA)
Sig.DEP.splsDA <- Sig.DEP.splsDA %>% dplyr::select(contains("TP"))
Sig.DEP.splsDA.zscore <- as.data.frame(t(apply(Sig.DEP.splsDA,1, z.score.func)))

heatmap.func(Sig.DEP.splsDA,plot.name = "DEP A-B & B-C (posthoc ANOVA paired)",cluster_cols=T)
heatmap.func(Sig.DEP.splsDA.zscore,plot.name = "DEP A-B & B-C (posthoc ANOVA paired), z-score",cluster_cols=T)

heatmap.func(Sig.DEP.splsDA,plot.name = "DEP (ANOVA paired)",cluster_cols=F)
heatmap.func(Sig.DEP.splsDA.zscore,plot.name = "DEP (ANOVA paired), z-score",cluster_cols=F)
