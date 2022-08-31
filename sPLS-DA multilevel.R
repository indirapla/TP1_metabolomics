#========01 TP STUDY(proteomics) =============

#   Paper: Plasma metabolome study reveals metabolic changes induced by pharmacological castration and testosterone supplementation in healthy young men.
# Authors: Jéssica de Siqueira Guedes; Indira Pla (indira.pla_parada@gmail.com); K. Barbara Sahlin; Gustavo Monnerat; Roger Appelqvist; György Marko-Varga; 
#          Aleksander Giwercman; Gilberto Barbosa Domont; Aniel Sanchez; Fábio César Sousa Nogueira; Johan Malm.



#====INSTALL PACKAGES================================================================================

# List of packages to install
.packages = c("BiocManager","devtools","mixOmics", "pheatmap","parallel","doParallel")


# Install packages if not installed
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  BiocManager::install(c("mixOmics"))
}

# Loading packages
lapply(.packages, require, character.only=TRUE)

#======END INSTALL PACKAGES============================================================================

#=====Loading dataset================================================================================== 
Big.table <- as.data.frame(readxl::read_excel("./data/normalized_data_testosterone_220221_ipp.xlsx")) # Open the data
rownames(Big.table) <- Big.table[,1]

main_table <- dplyr::select_if(Big.table,is.numeric)

Annotations <- as.data.frame(readxl::read_excel("./data/Annotations.xlsx")) # Open the data
row.names(Annotations) <- Annotations$sample
head(Annotations)


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
  print(hm)
  return(hm)
}

heatmap.obj <- heatmap.func(cim.matrix1,
                            plot.name = "cim-heatmap (101 metabolite from sPLS-DA multilevel)",
                            cluster_cols=T)


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

#=======Selecting the significant proteins from the sPLS-DA analysis=====

sPLS.DA_ANOVA.Sig <- subset(sPLS.DA_ANOVA,sPLS.DA_ANOVA$cluster!="NA")
sPLS.DA_ANOVA.Sig_mean <- sPLS.DA_ANOVA.Sig %>% dplyr::select(contains("mean"))  

sPLS.DA_ANOVA.Sig_Zmean <- as.matrix(sPLS.DA_ANOVA.Sig_mean[,1:3])

sPLS.DA_ANOVA.Sig_Zmean <- as.data.frame(t(apply(sPLS.DA_ANOVA.Sig_Zmean,1,z.score.func)))

sPLS.DA_ANOVA.Sig_Zmean$compound <- sPLS.DA_ANOVA.Sig$compound