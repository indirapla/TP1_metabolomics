#========01 TP STUDY(proteomics) =============

#   Paper: Plasma metabolome study reveals metabolic changes induced by pharmacological castration and testosterone supplementation in healthy young men.
# Authors: Jéssica de Siqueira Guedes; Indira Pla (indira.pla_parada@gmail.com); K. Barbara Sahlin; Gustavo Monnerat; Roger Appelqvist; György Marko-Varga; 
#          Aleksander Giwercman; Gilberto Barbosa Domont; Aniel Sanchez; Fábio César Sousa Nogueira; Johan Malm.



#====INSTALL PACKAGES================================================================================

# List of packages to install
.packages = c("BiocManager","devtools","ggplot2","reshape2","dplyr","ggpubr","lme4","nlme","car","plotly","RadaFDR","GOplot",
              "GOsummaries","randomcoloR","afex","ggbeeswarm",
              "emmeans","psych","caret")


# Install packages if not installed
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  install_github("vqv/ggbiplot")
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


# Function: z-score normalization by rows to plot a heatmap

z.score.func <- function(x){
  (x - mean(x))/sd(x)
}

##====Function ANOVA paired=======================================================

# This function performe one-way repeated measurements ANOVA per protein and then 
# performs a post-hoc test based on pairwise 'emmeans' function. 

#                  m: dataframe contating protein expresion per sample. Matrix with samples as columns and proteins as rows (Ex: m = main_table)
#   annotation.table: Annotation matrix. Matrix should contain a column called 'sample' containing samples IDs. (Ex: annotation.table = Annotations)
#             groups: name of the column that contain the groups that should be compared (factor) (Ex: groups = "Group")
#             id_ind: name of the column that containg the patients ID (Ex: id_ind="Patient")
#           plot.var: whether to plot or not boxplot per protein. If TRUE, plots will be saved in a folder called 'Boxplots' in dir() (Ex: FALSE)


library(afex)
library(ggbeeswarm)
library(emmeans)

DEP_prot <- function(m, annotation.table, groups, id_ind, plot.var=T) {
  
  groups1 <- annotation.table %>% dplyr::select(all_of(groups))
  groups1 <- as.factor(groups1[,1])
  n = length(levels(groups1))  # number of groups en la variable "groups" = 3
  k = 2                        # binary comparison
  
  comb <- factorial(n)/(factorial(k)*factorial(n-k))  # number of combinations
  
  
  #===create matrix where the results will be stored.
  
  p.value <- vector(mode = "numeric")
  pv <- matrix(nrow=nrow(m),ncol = comb, dimnames = list(row.names(m),c("pv(Basal-Low)","pv(Basal-Restor.)","pv(Low-Restor.)")))  ######
  diff <- matrix(nrow=nrow(m),ncol = comb,dimnames = list(row.names(m),c("Log2FC(Basal-Low)","Log2FC(Basal-Restor.)","Log2FC(Low-Restor.)")))
  means <- matrix(nrow=nrow(m),ncol = comb, dimnames = list(row.names(m),c("mean(Basal)","mean(Low)","mean(Restored)")))
  # lwr <- matrix(nrow=nrow(m),ncol = comb,dimnames = list(row.names(m),c("pv(Basal-Low)","pv(Basal-Restor.)","pv(Low-Restor.)")))
  # upr <- matrix(nrow=nrow(m),ncol = comb,dimnames = list(row.names(m),c("pv(Basal-Low)","pv(Basal-Restor.)","pv(Low-Restor.)")))
  
  #===create local folder where boxplot will be saved.
  dir.create("Boxplots")
  
  # ANOvA
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

# EXAMPLE
DEP <- as.data.frame(DEP_prot(m=main_table,
                              annotation.table=Annotations,
                              groups="Group",
                              id_ind="Patient",
                              plot.var = F))

#write.csv(DEP1, "DEP_fdr_ANOVA paired.csv")


###======Significant metabilites from ANOVA paired test========
DEP1$abs_Log2FC_AB <- abs(DEP1$`Log2FC(Basal-Low)`)
DEP1$abs_Log2FC_BC <- abs(DEP1$`Log2FC(Low-Restor.)`)


Sig.DEP1 <- subset(DEP1,DEP1$q.value<=0.05)
Sig.DEP1 <- subset(Sig.DEP1,Sig.DEP1$abs_Log2FC_AB>1)
Sig.DEP1 <- subset(Sig.DEP1,Sig.DEP1$abs_Log2FC_BC>1)

Sig.DEP2 <- Sig.DEP1 %>% dplyr::select(contains("TP"))

#write.csv(Sig.DEP2,"Sig.DEP2.csv")



