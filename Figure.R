library(ggplot2)
library(readr)
library(data.table)
library(dplyr)
library(gplots)
library(pheatmap)

setwd("/Users/khoahuynh/Library/Mobile Documents/com~apple~CloudDocs/Documents/Whole slide segmentation/Spatial_annotation/")

#--------------------------------------------------
# Load the cell annotations
#--------------------------------------------------

S1P1 = as.data.frame(fread("data_full_S1P1.csv"))
S2P1=as.data.frame(fread("data_full_S2P1.csv"))
S3P1=as.data.frame(fread("data_full_S3P1.csv"))
S3P2=as.data.frame(fread("data_full_S3P2.csv"))
S4P1=as.data.frame(fread("data_full_S4P1.csv"))
S4P2=as.data.frame(fread("data_full_S4P2.csv"))
S5P1=as.data.frame(fread("data_full_S5P1.csv"))
S5P2=as.data.frame(fread("data_full_S5P2.csv"))
S6P1=as.data.frame(fread("data_full_S6P1.csv"))
S6P1b=as.data.frame(fread("data_full_S6P1b.csv"))

#--------------------------------------------------
# Load the functions
#--------------------------------------------------

heatmap_anb=function(data_pos_neg,data_anb,group,name_anb){
  
  dat100=data.frame(data_pos_neg,group)
  colnames(dat100)[-ncol(dat100)]=colnames(data_pos_neg)
  
  SLIDE_1_value_anb=data.frame(data_anb,group)
  mean_value=count=NULL
  
  for(i in 1:length(names(table(group)))){
    
    dat9=dat100[which(dat100$group==names(table(group))[i]),]
    per_cent_count=as.numeric(sapply(dat9[,c(1:ncol(data_pos_neg))], sum))/nrow(dat9)
    count=append(count,per_cent_count)
    
    data_final9=SLIDE_1_value_anb[which(group==names(table(group))[i]),-ncol(SLIDE_1_value_anb)]
    per_cent_mean=as.numeric(sapply(as.data.frame(data_final9), median))
    mean_value=append(mean_value,per_cent_mean)
  }
  
  antibody_name=rep(colnames(dat9[,c(1:ncol(data_pos_neg))]),length(names(table(group))))
  Cluster=rep(names(table(group)),each=ncol(dat9[,c(1:ncol(data_pos_neg))]))
  plot_dat=data.frame(mean_value,count,antibody_name,Cluster)
  
  markers <- plot_dat$antibody_name %>% unique()
  
  plot_dat2=plot_dat[,c(1,3,4)]
  
  Cluster_order=names(table(group))
  
  
  
  
  ############ NOTICE THE t() above)
  
  
  
  plot_dat2=data.frame(plot_dat2,count)
  
  Anb_order=name_anb
  plot_dat2=plot_dat2[which(plot_dat2$antibody_name%in%Anb_order),]
  
  
  my_colors <- c("black", "grey70", "firebrick")
  
  Anb_order=name_anb
  plot_dat2[is.na(plot_dat2)==T]=0
  p1=plot_dat2 %>% filter(antibody_name %in% markers) %>% 
    mutate(`% Expressing` = count * 100,
           antibody_name = factor(antibody_name, levels = Anb_order),
           Cluster = factor(Cluster, levels = Cluster_order)) %>% 
    ggplot(aes(x=antibody_name, y = Cluster, color = mean_value, size = `% Expressing`)) + 
    geom_point() + 
    cowplot::theme_cowplot() + 
    #  theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 40, hjust=1)) +
    ylab('') +
    #  theme(axis.ticks = element_blank()) +
    scale_color_gradientn(colours = my_colors, oob = scales::squish, name = 'Mean value', limits = c(-1,1))  +
    scale_y_discrete(position = "right")+xlab("")+ggtitle("")
  
  return(p1)
}


#--------------------------------------------------
# Code to reproduce figure SF4A
#--------------------------------------------------


name_anb=c("EPCAM"     ,  "Keratin8_18", "PanCK"     ,  "CD66", "Keratin14"  , "SMA"     ,    "VIM"      ,   "CD31"     ,   "CD34"  ,     
           "Podoplanin" , "CD45"   ,     "HLA_DR"  ,    "CD11c"       ,      "CD141"   ,    "CD56"    ,    "CD3e"   ,     "FOXP3"  ,    
           "CD4"    ,     "CD8"     ,    "CD20" ,"HHP3"  ,     "GATA3"    ,  "GZMB" ,      "Ki67"    ,   "Galectin3" , "HLA_A"   ,   "CD45RO"   ,  "CD21"    ,   "PD_L1"  ,   
           "PD1"   ,     "IDO1"   ,    "CD107a"  ,   "ICOS"   ,    "Bcl_2"  ,          "IFNG" )

Healthy_data=rbind(cbind(predicted=S2P1$TACIT, S2P1[,paste0("P_",name_anb,sep="")]),
                   cbind(predicted=S4P1$TACIT, S4P1[,paste0("P_",name_anb,sep="")]),
                   cbind(predicted=S4P2$TACIT, S4P2[,paste0("P_",name_anb,sep="")]),
                   cbind(predicted=S5P1$TACIT, S5P1[,paste0("P_",name_anb,sep="")]),
                   cbind(predicted=S5P2$TACIT, S5P2[,paste0("P_",name_anb,sep="")]),
                   cbind(predicted=S6P1$TACIT, S6P1[,paste0("P_",name_anb,sep="")]),
                   cbind(predicted=S6P1b$TACIT, S6P1b[,paste0("P_",name_anb,sep="")])
)



Healthy_data_expression=rbind(cbind(scale(S2P1[,name_anb])),
                              cbind(scale(S4P1[,name_anb])),
                              cbind(scale(S4P2[,name_anb])),
                              cbind(scale(S5P1[,name_anb])),
                              cbind(scale(S5P2[,name_anb])),
                              cbind(scale(S6P1[,name_anb])),
                              cbind(scale(S6P1b[,name_anb])))


data_pos_dotplot=Healthy_data[,-1]


data_expression_dotplot=Healthy_data_expression

colnames(data_pos_dotplot)=name_anb

Signature <- readxl::read_excel("Signature_blake.xlsx")
Healthy_data$predicted=factor(Healthy_data$predicted,levels = c(Signature$cell_type[-14],"Others"))
data_expression_dotplot=as.data.frame(data_expression_dotplot)
data_pos_dotplot=as.data.frame(data_pos_dotplot)
heatmap_anb(data_anb = scale(data_expression_dotplot),data_pos_neg = data_pos_dotplot,group = Healthy_data$predicted,name_anb = name_anb)















#--------------------------------------------------
# Code to reproduce figure SF4B
#--------------------------------------------------


#####Group tissues by healthy and disease
Group_Healthy=c(rep("Healthy",nrow(S2P1)),rep("Healthy",nrow(S4P1)),rep("Healthy",nrow(S4P2)),rep("Healthy",nrow(S5P1)),rep("Healthy",nrow(S5P2)),rep("Healthy",nrow(S6P1)),rep("Healthy",nrow(S6P1b)))
Group_Disease=c(rep("Disease",nrow(S1P1)),rep("Disease",nrow(S3P1)),rep("Disease",nrow(S3P2)))
Group=c(Group_Healthy,Group_Disease)


#####Group tissues by cell type
CT_healthy=c(S2P1$TACIT,S4P1$TACIT,S4P2$TACIT,S5P1$TACIT,S5P2$TACIT,S6P1$TACIT,S6P1b$TACIT)
CT_disease=c(S1P1$TACIT,S3P1$TACIT,S3P2$TACIT)
CT=c(CT_healthy,CT_disease)



data_plot_proportion=data.frame(Group,CT)
data_plot_proportion$Group=factor(data_plot_proportion$Group,levels = c("Healthy","Disease"))
# Compute the proportions
prop_data <- data_plot_proportion %>%
  group_by(Group, CT) %>%
  tally() %>%
  group_by(Group) %>%
  mutate(proportion = n/sum(n))


# Plot using ggplot2 
ggplot(prop_data, aes(x=CT, y=proportion, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  labs(y="Proportion of cells") +
  theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label=round(proportion,2)), position=position_dodge(width=0.9), vjust=-0.5)+xlab("")+
  scale_fill_manual(values = c("Healthy"="#F47D20","Disease"="#4372C4")) 





#--------------------------------------------------
# Code to reproduce figure SF4C
#--------------------------------------------------




name_anb=c("EPCAM"     ,  "Keratin8_18", "PanCK"     ,  "CD66", "Keratin14"  , "SMA"     ,    "VIM"      ,   "CD31"     ,   "CD34"  ,     
           "Podoplanin" , "CD45"   ,     "HLA_DR"  ,    "CD11c"       ,      "CD141"   ,    "CD56"    ,    "CD3e"   ,     "FOXP3"  ,    
           "CD4"    ,     "CD8"     ,    "CD20" ,"HHP3"  ,     "GATA3"    ,  "GZMB" ,      "Ki67"    ,   "Galectin3" , "HLA_A"   ,   "CD45RO"   ,  "CD21"    ,   "PD_L1"  ,   
           "PD1"   ,     "IDO1"   ,    "CD107a"  ,   "ICOS"   ,    "Bcl_2"  ,          "IFNG" )


Healthy_data=rbind(cbind(predicted=S2P1$TACIT,S2P1[,name_anb]),
                   cbind(predicted=S4P1$TACIT,S4P1[,name_anb]),
                   cbind(predicted=S4P2$TACIT,S4P2[,name_anb]),
                   cbind(predicted=S5P1$TACIT,S5P1[,name_anb]),
                   cbind(predicted=S5P2$TACIT,S5P2[,name_anb]),
                   cbind(predicted=S6P1$TACIT,S6P1[,name_anb]),
                   cbind(predicted=S6P1b$TACIT,S6P1b[,name_anb])
)




Disease_data=rbind(cbind(predicted=S1P1$TACIT,S1P1[,name_anb]),
                   cbind(predicted=S3P1$TACIT,S3P1[,name_anb]),
                   cbind(predicted=S3P2$TACIT,S3P2[,name_anb]))



data_expression_healthy=data.frame(Group="Healthy",ct=Healthy_data$predicted,Healthy_data[,-1])
data_expression_disease=data.frame(Group="Disease",ct=Disease_data$predicted,Disease_data[,-1])
data_expression_plot=rbind(data_expression_healthy,data_expression_disease)
colnames(data_expression_plot)=c("Group","ct",name_anb)
data_expression_plot=data_expression_plot
data_expression_plot=data_expression_plot[which(data_expression_plot$ct!="Others"),]
healthy_data <- subset(data_expression_plot, Group == "Healthy")
disease_data <- subset(data_expression_plot, Group == "Disease")
# Define the markers you are interested in
markers <- colnames(healthy_data)[grep("EPCAM|Keratin8_18|PanCK|CD66|Keratin14|SMA|VIM|CD31|CD34|Podoplanin|CD45|HLA_DR|CD11c|CD141|CD56|CD3e|FOXP3|CD4|CD8|CD20|HHP3|GATA3|GZMB|Ki67|Galectin3|HLA_A|CD45RO|CD21|PD_L1|PD1|IDO1|CD107a|ICOS|Bcl_2|IFNG", 
                                       colnames(healthy_data))]
# Create an empty data frame to store the results
results <- data.frame(ct = character(0), marker = character(0), log2FC = numeric(0))

# Loop through each cell type and marker
for (ct in unique(healthy_data$ct)) {
  for (marker in markers) {
    # Extract the expression values for the current cell type and marker
    healthy_expr <- healthy_data[healthy_data$ct == ct, marker]
    disease_expr <- disease_data[disease_data$ct == ct, marker]
    
    # Calculate the log2 fold change
    log2FC <- log2(mean(disease_expr) / mean(healthy_expr))
    
    # Append the result to the data frame
    results <- rbind(results, data.frame(ct = ct, marker = marker, log2FC = log2FC))
  }
}

# Print the results
print(results)

# Load the required library for creating heatmaps
library(gplots)

# Pivot the results data frame for heatmap
heatmap_data <- reshape(results, idvar = "marker", timevar = "ct", direction = "wide")

# Set row names to marker and remove the marker column
rownames(heatmap_data) <- heatmap_data$marker
heatmap_data <- heatmap_data[, -1]
colnames(heatmap_data) <- sub("log2FC\\.", "", colnames(heatmap_data))
my.breaks <- c(seq(-2, 0, by=0.1),seq(0.1, 2, by=0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "red"))(length(my.breaks)/2))
pheatmap(scale(heatmap_data), cluster_rows=T, show_rownames=T,
            cluster_cols=T,  color=my.colors,breaks = my.breaks,main = "")




#--------------------------------------------------
# Code to reproduce figure F4G
#--------------------------------------------------



user_colors <-  c(
  "Acini" = "orange", "Ducts" = "green", "Myoepithelial" = "#8C33FF",
  "Fibroblast" = "blue", "Myofibroblast" = "yellow4",
  "VEC" = "black", "VEC Progen" = "#e7298a", "LECs" = "darkcyan",
  "B cells" = "red", "DC cells" = "pink", "DN" = "violet", "NK cells" = "skyblue", "Tc" = "#E76F51", "Treg" = "darkblue",
  "Macrophage" = "tan", "Th" = "cyan", "DP" = "royalblue", "Others" = "grey"
 )

plot_colors <- user_colors
plot_colors[!(names(user_colors) %in% input$selectedTACIT)] <- "grey"

data_plot=data.frame(X=S5P2$X[which(S5P2$Group.x=="S5P1_1")],Y=S5P2$Y[which(S5P2$Group.x=="S5P1_1")],TACIT=S5P2$TACIT[which(S5P2$Group.x=="S5P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = TACIT)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values = plot_colors) +
  theme_classic(base_size = 15)+ guides(color = guide_legend(override.aes = list(size = 3),title = "Cell type"))+ggtitle("Spatial plot with annotation")





data_plot=data.frame(X=S5P2$X[which(S5P2$Group.x=="S5P1_1")],Y=S5P2$Y[which(S5P2$Group.x=="S5P1_1")],CD3e=S5P2$CD3e[which(S5P2$Group.x=="S5P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = CD3e)) +
  geom_point(size = 0.1) +ggtitle("CD3e expressing")+
  scale_color_gradient(low = "grey", high = "red", limits = c(0, 25), na.value = "red")+theme_classic(base_size = 15)



data_plot=data.frame(X=S5P2$X[which(S5P2$Group.x=="S5P1_1")],Y=S5P2$Y[which(S5P2$Group.x=="S5P1_1")],CD8=S5P2$CD8[which(S5P2$Group.x=="S5P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = CD8)) +
  geom_point(size = 0.1) +ggtitle("CD8 expressing")+
  scale_color_gradient(low = "grey", high = "red", limits = c(0, 25), na.value = "red")+theme_classic(base_size = 15)




data_plot=data.frame(X=S5P2$X[which(S5P2$Group.x=="S5P1_1")],Y=S5P2$Y[which(S5P2$Group.x=="S5P1_1")],HLA_A=S5P2$HLA_A[which(S5P2$Group.x=="S5P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = HLA_A)) +
  geom_point(size = 0.1) +ggtitle("HLA_A expressing")+
  scale_color_gradient(low = "grey", high = "red", limits = c(0, 60), na.value = "red")+theme_classic(base_size = 15)



data_plot=data.frame(X=S5P2$X[which(S5P2$Group.x=="S5P1_1")],Y=S5P2$Y[which(S5P2$Group.x=="S5P1_1")],Keratin8_18=S5P2$Keratin8_18[which(S5P2$Group.x=="S5P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = Keratin8_18)) +
  geom_point(size = 0.1) +ggtitle("Keratin8_18 expressing")+
  scale_color_gradient(low = "grey", high = "red", limits = c(0, 25), na.value = "red")+theme_classic(base_size = 15)










data_plot=data.frame(X=S1P1$X[which(S1P1$Group.x=="S1P1_1")],Y=S1P1$Y[which(S1P1$Group.x=="S1P1_1")],TACIT=S1P1$TACIT[which(S1P1$Group.x=="S1P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = TACIT)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values = plot_colors) +
  theme_classic(base_size = 15)+ guides(color = guide_legend(override.aes = list(size = 3),title = "Cell type"))+ggtitle("Spatial plot with annotation")





data_plot=data.frame(X=S1P1$X[which(S1P1$Group.x=="S1P1_1")],Y=S1P1$Y[which(S1P1$Group.x=="S1P1_1")],CD3e=S1P1$CD3e[which(S1P1$Group.x=="S1P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = CD3e)) +
  geom_point(size = 0.1) +ggtitle("CD3e expressing")+
  scale_color_gradient(low = "grey", high = "red", limits = c(0, 25), na.value = "red")+theme_classic(base_size = 15)



data_plot=data.frame(X=S1P1$X[which(S1P1$Group.x=="S1P1_1")],Y=S1P1$Y[which(S1P1$Group.x=="S1P1_1")],CD8=S1P1$CD8[which(S5P2$Group.x=="S1P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = CD8)) +
  geom_point(size = 0.1) +ggtitle("CD8 expressing")+
  scale_color_gradient(low = "grey", high = "red", limits = c(0, 25), na.value = "red")+theme_classic(base_size = 15)




data_plot=data.frame(X=S1P1$X[which(S1P1$Group.x=="S1P1_1")],Y=S1P1$Y[which(S1P1$Group.x=="S1P1_1")],HLA_A=S1P1$HLA_A[which(S1P1$Group.x=="S1P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = HLA_A)) +
  geom_point(size = 0.1) +ggtitle("HLA_A expressing")+
  scale_color_gradient(low = "grey", high = "red", limits = c(0, 60), na.value = "red")+theme_classic(base_size = 15)



data_plot=data.frame(X=S1P1$X[which(S1P1$Group.x=="S1P1_1")],Y=S1P1$Y[which(S1P1$Group.x=="S1P1_1")],Keratin8_18=S1P1$Keratin8_18[which(S1P1$Group.x=="S1P1_1")])

ggplot(data_plot, aes(x = X, y = Y, color = Keratin8_18)) +
  geom_point(size = 0.1) +ggtitle("Keratin8_18 expressing")+
  scale_color_gradient(low = "grey", high = "red", limits = c(0, 25), na.value = "red")+theme_classic(base_size = 15)












