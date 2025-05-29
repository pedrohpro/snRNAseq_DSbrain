
# setwd("/mnt/sandbox-SSD-1/poliveira/bruna/scDown")
# getwd() # checking

rm(list = ls())     # clears Environment 

set.seed(123)


##################
###  Functions ###
##################

markersFunction <- function(){
  
  celltypes_markers <- markers_df %>% 
    filter(marker %in% rownames(scdown)) %>% 
    pull(marker)
  
  markers_df <- markers_df %>% 
    filter(marker %in% rownames(scdown))
  
  auxmark <- markers_df %>% 
    group_by(celltype) %>% 
    summarise(n_markers = n()) %>% 
    mutate(celltype = factor(celltype, levels = unique(markers_df$celltype))) %>% 
    arrange(celltype) %>%  
    mutate(
      y = length(unique(scdown@meta.data$seurat_clusters))*1.05,
      yend = y,
      x = 1, # (just for the loop to work)
      xend = n_markers, # provisory (just for the loop to work)
      colour = "#000000",
      celltype = case_when(
        celltype == "Ast" ~ "Astrocyte",
        celltype == "End" ~ "Endothelial",
        celltype == "Exc" ~ "Excitatory neuron",
        celltype == "Inh" ~ "Inhibitory neuron",
        celltype == "Mic" ~ "Microglia",
        celltype == "Oli" ~ "Oligodendrocyte",
        celltype == "OPC" ~ "Oligodendrocyte Progenitor",
        celltype == "Per" ~ "Pericyte",
        celltype == "NPC" ~ "Neural Progenitor"
      )
    ) %>% 
    as.data.frame()
  for (i in 2:nrow(auxmark)) {
    auxmark[i,"x"] <- auxmark[i-1,"xend"] + 1
    auxmark[i,"xend"] <- auxmark[i,"x"] + auxmark[i,"n_markers"] -1
  }
  
  # canonical markers
  dot_plot <- DotPlot(
    object = scdown, 
    features = celltypes_markers, 
    cols = c("blue","red"),
    cluster.idents = T,
    scale = T 
  )+
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),  # Remove panel borders
      axis.line = element_line(color = "black"),  # Set axis line color
      legend.background = element_rect(fill = "white"),  # Set legend background color to white
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
      axis.text.y = element_text(size = 17, face = "bold"),
      #legend.position = "none" # no legend
    )
  for (i in 1:nrow(auxmark)) {
    dot_plot <- dot_plot + 
      annotate(
        "segment", # type of annotation to be added
        x = auxmark[i,"x"],
        y = auxmark[i,"y"],
        xend = auxmark[i,"xend"],
        yend = auxmark[i,"yend"],
        #colour = auxmark[i,"colour"],
        
      ) +
      annotate(
        "text", # type of annotation to be added
        x = (auxmark[i,"x"] + auxmark[i,"xend"])/2, 
        y = auxmark[i,"y"] + 0.15, 
        label = auxmark[i,"celltype"], 
        #hjust = 1.1, 
        #vjust = -0.5, 
        color = "black",
        fontface = "bold",
        size = 4
      ) +
      coord_cartesian(clip = "off")
  }
  
  return(dot_plot)
}

barfunc <- function(
    df,
    x_axis,
    y_axis,
    fill_groups = F,
    annot_loc = "center", # "center" or "repel" or "top" or F
    fill_col = F,
    title_name = NULL,
    annot_text = y_axis,
    size_annot = 4,
    size_xaxis_text = 15,
    size_xaxis_title = 15,
    size_yaxis_text = 15,
    size_yaxis_title = 15,
    rotate_x = 0,
    x_title = F,
    y_title = F,
    extra_annot = NULL
){
  # Base plot
  bar_plot <- ggplot(df, aes(x = !!sym(x_axis), y = !!sym(y_axis))) +
    geom_bar(stat = "identity") +
    labs(
      title = title_name,
      x = x_axis,
      y = y_axis
    ) +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) + # remove extra padding between bars and x axis
    theme(
      plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
      panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
      panel.border = element_blank(),  # Remove panel borders
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
      axis.title.x = element_text(face = "bold", size = size_xaxis_title),  # Make axis titles bold
      axis.title.y = element_text(face = "bold", size = size_yaxis_title),  # Make axis titles bold
      axis.text.y = element_text(face = "bold", size = size_yaxis_text),  # Make axis text bold
      axis.text.x = element_text(face = "bold", angle = rotate_x, hjust = ifelse(rotate_x != 0, 1, 0.5), size = size_xaxis_text), # rotate x axis test 45 degrees
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.ticks.y = element_line(color = "black", linewidth = 0.5)  # Make y-axis ticks bold
      #legend.position = "none" # No legend
    )
  # y axis ticks
  if(size_yaxis_text == 0){
    bar_plot <- bar_plot + theme(axis.ticks.y = element_blank())
  }
  
  # Conditionally add fill aesthetic and scale_fill_manual
  if (fill_groups != F) {
    bar_plot <- bar_plot +
      aes(fill = !!sym(fill_groups)) +
      labs(fill = fill_groups)
  }
  
  # coloring bars
  if(is.list(fill_col)){
    bar_plot <- bar_plot +
      scale_fill_manual(values = fill_col)
  } else if (fill_col != F) {
    bar_plot <- bar_plot +
      geom_bar(stat = "identity", fill = fill_col) +
      theme(legend.position = "none")
  }
  
  # Annotations location
  if(!is.null(extra_annot)){
    bar_plot <- bar_plot +
      annotate(
        geom = "text",
        x = extra_annot[[x_axis]],
        y = extra_annot$y,
        label = extra_annot$annotation_text,
        color = extra_annot$color,
        size = 3,
        angle = 90,
        fontface = "bold"
      )
  }else if (annot_loc == "center") {
    bar_plot <- bar_plot +
      geom_text(
        aes(label = !!sym(annot_text)),
        position = position_stack(vjust = 0.5),
        size = size_annot
      )
  } else if (annot_loc == "repel") {
    bar_plot <- bar_plot +
      geom_text_repel(
        aes(label = !!sym(annot_text)),
        direction = "y",
        size = size_annot
      )
  } else if(annot_loc == "top"){
    bar_plot <- bar_plot +
      scale_y_continuous(expand = c(0,0), limits = c(0, 1.1*max(ifelse(fill_groups != F, df %>% group_by(!!sym(fill_groups)) %>% pull(!!sym(y_axis)), df[[y_axis]])))) +
      geom_text(
        aes(label = !!sym(annot_text)),
        vjust = -0.2,
        size = size_annot
      )
  }
  # x label title
  if(x_title != F){
    bar_plot <- bar_plot +
      xlab(x_title)
  }
  # y label title
  if(y_title != F){
    bar_plot <- bar_plot +
      ylab(y_title)
  }
  
  return(bar_plot)
}


########################
### loading library  ###
########################

library(devtools)
library(tidyverse) #do it all
library(Seurat) # single cell analysis
library(data.table)
library(org.Hs.eg.db) #Homo sapiens OrgDb
library(biomaRt) #gene annotations
library(RColorBrewer) # for pretty colors
library(ggrepel)
library(umap)
library(VennDiagram) #venn digrams
library(clusterProfiler) #functional analysis
library(pathview)
library(cowplot)
library(GOplot)
library(enrichplot)
library(gage)
library(vsn)
library(RColorBrewer) # colors pallets
library(viridis) # color pallets
library(scales) # colors pallets
library(metap) # for pvalues
library(rsvd) # for ALRA
library(ALRA) # zero imputation
library(SeuratWrappers) # third party seurat functions (eg.: RunALRA())
library(DESeq2)
library(patchwork)
library(scDblFinder) # doublet finder
library(SingleCellExperiment)
library(AUCell) # get set score
library(GSEABase) # build gene set function
library(sessioninfo) # package versions
library(chromo)


############################
### creating directories ###
############################

if (!file.exists("results_noALRA")){   
  dir.create("results_noALRA")
  dir.create("results_noALRA/exploratory_analysis")
}


############################
###         counts       ###
############################

# Runs names lists (maybe better to get this from meta data)
#runs <- sapply(43:71, function(i) paste0("EGAN000033605", i), simplify = T) # for all samples
runs <- c("EGAN00003360559", "EGAN00003360560", "EGAN00003360561", "EGAN00003360562", "EGAN00003360563", "EGAN00003360564", # young samples
          "EGAN00003360565", "EGAN00003360566", "EGAN00003360567", "EGAN00003360568", "EGAN00003360569", "EGAN00003360570")

# Raw counts
counts <- list()
for (run in runs){ 
  counts[[run]] <- Read10X(paste0("data/counts/",run,"/"))
}

# Separate objects
seurats <- list()
for (run in names(counts)) {
  seurats[[run]] <- CreateSeuratObject(counts = counts[[run]], project = run)
}

# merging seurat objects (without batch correction)
scdown <- merge(
  x = seurats$EGAN00003360559, # change the starting sample according to samples used
  y = seurats[names(seurats) != "EGAN00003360559"], # change the starting sample according to samples used
  add.cell.ids = runs,
  project = "allRuns"
)

# joining layers
scdown <- JoinLayers(scdown)


############################
###      meta data       ###
############################

meta <- read.delim("data/meta.csv", sep = ",") %>% 
  select(accession_id, biological_sex, subject_id, phenotype) %>%  # filtering columns
  rename( # renaming columns
    condition = phenotype,
    sex = biological_sex,
    donor = subject_id
  ) %>% 
  mutate(
    condition = case_when( # renaming elements
      condition == "non-diseased" ~ "CT",
      condition == "Down Syndrome" ~ "DS"
    ),
    age = substr(donor, start = nchar(donor) - 1, stop = nchar(donor) - 1) # adding age column
  ) %>% 
  filter( # getting only young samples to avoid the impact of Alzheimer's Disease neuroinflammation
    age == "Y"
  ) %>% 
  mutate_all(factor) # factoring all columns

# Including the sample meta data into the merged cell meta data
scdown@meta.data <- scdown@meta.data %>% 
  rownames_to_column("cell") %>% 
  left_join(meta, by = c("orig.ident" = "accession_id")) %>% 
  column_to_rownames("cell")

# Creating mitochondrial reads percentage
scdown[["percent.mt"]] <- PercentageFeatureSet(scdown, pattern = "^MT-") # high percentage could indicate either cell death or high metabolism cell type


#############################
###   genes of interest   ###
#############################

# # genes of interest
int_genes_df <- read.delim("data/genes_down.csv", sep = ";")
int_genes <- int_genes_df$genesofinterest
comp_genes <- int_genes_df$complement[int_genes_df$complement != ""]

# cell type markers
markers_df <- read.delim("data/celltypes_markers.csv", sep = ";")
celltypes_markers <- markers_df$marker

genes_fig1 <- c('C1R', 'C1S', 'C4A', 'C4B', 'C5', 'C8G', 'CD46', 'CD59', 'CFI', 'FCN2', 'SLC1A2', 'SLC1A3')


############################
###   filtering cells    ###
############################

scdown <- subset(scdown, subset = nCount_RNA > 200 & nCount_RNA < 60000) # Number of reads per cell

scdown <- subset(scdown, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000) # Number of features per cell

scdown <- subset(scdown, subset = percent.mt < 10)# Percentage of reads mapped to mitochondrial genes 

# doublets detection
sce <- as.SingleCellExperiment(scdown)
sce <- scDblFinder(sce)
scdown$doublet <- sce$scDblFinder.class
scdown <- subset(scdown, subset = doublet == "singlet")


###############################
###   filtering features    ###
###############################

# Removing mitochondrial genes
scdown <- scdown[!grepl("MT-", rownames(scdown)),]

# features information
chrObj <- chromoInitiate( # filters for chromosome and MT transcripts
    data.frame(Symbol = rownames(scdown), log2fc = 0, pval = 1),
    gene_col = "Symbol",
    fc_col = "log2fc",
    p_col = "pval"
  )
all_features <- chrObj@data %>% 
  dplyr::select(-log2fc, -pval, -DEG, -entrezgene_id)
scdown <- scdown[rownames(scdown) %in% all_features$Symbol,]

# Removing genes with zero expression
nonzero_genes <- rowSums(scdown[["RNA"]]@layers[["counts"]]) > 0
scdown <- scdown[nonzero_genes, ]
all_features <- all_features[all_features$Symbol %in% rownames(scdown),]


######################
###   Processing   ###
######################

scdown <- NormalizeData(scdown)

gc() # clears garbage (next step requires a lot of memory)

scdown <- FindVariableFeatures(scdown) # default 2000 features

top10 <- head(VariableFeatures(scdown), 10) # Highlight the 10 most highly variable genes

scdown <- ScaleData(scdown, features = rownames(scdown))

scdown <- RunPCA(scdown, features = VariableFeatures(object = scdown)) 

scdown <- RunUMAP(scdown, dims = 1:10)

scdown <- FindNeighbors(scdown, dims = 1:10, k.param = 30)

scdown <- FindClusters(scdown, resolution = 0.2)


#########################################
###      Batch effect verification    ###
#########################################

# UMAP
for (i in c("seurat_clusters","donor","sex","condition")) {
  clusters <- DimPlot(scdown, reduction = "umap", group.by = i, pt.size = 0.1, label = TRUE, raster=FALSE) +
    ggtitle("") +
    labs(color = i) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    coord_fixed() +
    theme(
      legend.text = element_text(size = 26),  # Adjust legend text size and style
      legend.title = element_text(size = 28),  # Adjust legend text size and style
      axis.text = element_text(size = 16),    # Adjust axis text size and style
      axis.title = element_text(size = 22)
    )
  
  ggsave(paste0("results_noALRA/exploratory_analysis/umap_",i,".png"), plot = clusters, width = 10, height = 10)
}


#################################
###      Bias verification    ###
#################################

# condition vs. sex
aux <- scdown@meta.data %>%
  select(condition, sex) %>%
  group_by(condition, sex) %>%
  summarise(num_cells = n())%>%
  group_by(condition) %>%
  mutate(pct = num_cells / sum(num_cells) * 100) %>%
  ungroup()

# pct value bar plot
bar_plot <- barfunc(aux, x_axis = "condition", y_axis = "pct", fill_groups = "sex", annot_loc = "center", annot_text = "num_cells")
ggsave("results_noALRA/exploratory_analysis/pct_condition_sex.png", plot = bar_plot, width = 3, height = 3)


####################################
###       Finding markers        ###
####################################

# Verifying assay type. Should be RNA
DefaultAssay(scdown)

# Changing identity to cluster column
Idents(scdown) <- "seurat_clusters"

dot_plot <- markersFunction()
ggsave("results_noALRA/exploratory_analysis/canonical_markers.png", plot = dot_plot, height = 5, width = 20)


###############################
###       Annotating        ###
###############################
Idents(scdown) <- "seurat_clusters"

scdown$before_annotating <- scdown$seurat_clusters

scdown <- RenameIdents(scdown, c(
  `0` = "Oli",
  `1` = "Ast",
  `2` = "Exc",
  `3` = "OPC",
  `4` = "Mic",
  `5` = "Exc",
  `6` = "Inh",
  `7` = "Inh",
  `8` = "Exc",
  `9` = "End",
  `10` = "End",
  `11` = "Ast",
  `12` = "OPC"
  )
)

scdown@meta.data <- scdown@meta.data %>%
  mutate(
    seurat_clusters = case_when(
      seurat_clusters == "0" ~ "Oli",
      seurat_clusters == "1" ~ "Ast",
      seurat_clusters == "2" ~ "Exc",
      seurat_clusters == "3" ~ "OPC",
      seurat_clusters == "4" ~ "Mic",
      seurat_clusters == "5" ~ "Exc",
      seurat_clusters == "6" ~ "Inh",
      seurat_clusters == "7" ~ "Inh",
      seurat_clusters == "8" ~ "Exc",
      seurat_clusters == "9" ~ "End",
      seurat_clusters == "10" ~ "End",
      seurat_clusters == "11" ~ "Ast",
      seurat_clusters == "12" ~ "OPC"
    )
  )

# ordering in decreasing number of cells
scdown@meta.data$seurat_clusters <- factor(scdown@meta.data$seurat_clusters, levels = c("Oli","Exc","Ast","OPC","Inh","Mic","End"))


########################
###       RDS        ###
########################

#saveRDS(scdown, file = "data/RDS/noALRA.rds")
saveRDS(scdown, file = "data/RDS/noALRA.rds")

#scdown <- readRDS("data/RDS/noALRA.rds")
scdown <- readRDS("data/RDS/noALRA.rds")

chrObj <- chromoInitiate( # filters for chromosome and MT transcripts
  data.frame(Symbol = rownames(scdown), log2fc = 0, pval = 1),
  gene_col = "Symbol",
  fc_col = "log2fc",
  p_col = "pval"
)
all_features <- chrObj@data %>% 
  dplyr::select(-log2fc, -pval, -DEG, -entrezgene_id)


###################################
###       Annotated UMAP        ###
###################################

Idents(scdown) <- "seurat_clusters"

# All cells
umap_plot <- DimPlot(
  scdown, 
  reduction = "umap", 
  label = F, 
  pt.size = 0.1, 
  raster = FALSE, 
  label.size = 12, 
  cols = c("Oli" = brewer.pal(7, "Set3")[1], "Exc" = brewer.pal(8, "Set3")[8], "OPC" = brewer.pal(7, "Set3")[3], 
           "Inh" = brewer.pal(7, "Set3")[4], "Ast" = brewer.pal(7, "Set3")[7], "Mic" = brewer.pal(7, "Set3")[6], 
           "End" = brewer.pal(7, "Set3")[5])
) +
  NoAxes() + # removes axis
  theme(
    legend.position = "none"
  )
ggsave("results_noALRA/exploratory_analysis/ANNOTATED_umap.png", plot = umap_plot, height = 12, width = 12, dpi = 600)

# CT cells
sub_scdown <- subset(scdown, condition == 'CT')

umap_plot <- DimPlot(
  sub_scdown, 
  reduction = "umap", 
  label = F, 
  pt.size = 0.1, 
  raster = FALSE, 
  label.size = 12,
  cols = c("Oli" = brewer.pal(7, "Set3")[1], "Exc" = brewer.pal(8, "Set3")[8], "OPC" = brewer.pal(7, "Set3")[3], 
           "Inh" = brewer.pal(7, "Set3")[4], "Ast" = brewer.pal(7, "Set3")[7], "Mic" = brewer.pal(7, "Set3")[6], 
           "End" = brewer.pal(7, "Set3")[5])
) +
  NoAxes() + # removes axis
  theme(
    legend.position = "none"
  )
rm(sub_scdown)
ggsave("results_noALRA/exploratory_analysis/ANNOTATED_umap_CT.png", plot = umap_plot, height = 12, width = 12, dpi = 600)

# DS cells
sub_scdown <- subset(scdown, condition == 'DS')

umap_plot <- DimPlot(
  sub_scdown, 
  reduction = "umap", 
  label = F, 
  pt.size = 0.1, 
  raster = FALSE, 
  label.size = 12,
  cols = c("Oli" = brewer.pal(7, "Set3")[1], "Exc" = brewer.pal(8, "Set3")[8], "OPC" = brewer.pal(7, "Set3")[3], 
           "Inh" = brewer.pal(7, "Set3")[4], "Ast" = brewer.pal(7, "Set3")[7], "Mic" = brewer.pal(7, "Set3")[6], 
           "End" = brewer.pal(7, "Set3")[5])
  #cols = rep("#22222277", length(unique(scdown@meta.data$seurat_clusters)))
) +
  NoAxes() + # removes axis
  theme(
    legend.position = "none"
  )
rm(sub_scdown)
ggsave("results_noALRA/exploratory_analysis/ANNOTATED_umap_DS.png", plot = umap_plot, height = 12, width = 12, dpi = 600)


####################################################
###      Repeating Batch effect verification     ###
####################################################

# Bar plots
for (i in c("donor","sex","condition")) {
  aux <- scdown@meta.data %>%
    select(seurat_clusters, !!sym(i)) %>%
    group_by(seurat_clusters, !!sym(i)) %>%
    summarise(num_cells = n())%>%
    group_by(seurat_clusters) %>%
    mutate(pct = num_cells / sum(num_cells) * 100) %>%
    ungroup()

  # pct value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "pct", fill_groups = i, annot_loc = "center", annot_text = "num_cells")
  ggsave(paste0("results_noALRA/exploratory_analysis/ANNOTATED_pctbar_",i,".png"), plot = bar_plot, width = 5, height = 4)
}


##################################################
###          Repeating markers plots           ###
##################################################

dot_plot <- markersFunction()
ggsave("results_noALRA/exploratory_analysis/ANNOTATED_dotplot.png", plot = dot_plot, height = 4, width = 18)


###########################################
###        Cell type composition        ###
###########################################

# Absolute value
aux <- scdown@meta.data %>%
 select(orig.ident, seurat_clusters, condition) %>%
 group_by(orig.ident, condition, seurat_clusters) %>%
 summarise(num_cells = n()) %>%
 group_by(orig.ident) %>%
 mutate(pct = num_cells / sum(num_cells) * 100) %>%
 ungroup()

# separate CT and DS horizontally
rowOrder <- meta %>%
 arrange(condition) %>%
 pull(accession_id)
aux$orig.ident <- factor(aux$orig.ident, levels = rowOrder)

# Percentage value
pct_composition <- barfunc(
 aux, x_axis = "orig.ident", y_axis = "pct", fill_groups = "seurat_clusters", annot_loc = "center", annot_text = "num_cells",
 #fill_col = c(brewer.pal(n = 10, name = "Set3"), brewer.pal(n = nlevels(aux$seurat_clusters)-10, name = "Paired")),
 rotate_x = 45, size_xaxis_text = 8)
ggsave("results_noALRA/exploratory_analysis/ANNOTATED_pct_composition.png", plot = pct_composition, height = 5, width = 6)


##############################################
###    Differential expression analysis    ###
##############################################

# creating composite column
scdown@meta.data$type_cond <- paste0(scdown@meta.data$seurat_clusters, "_", scdown@meta.data$condition)

# Changing identity to composite column
Idents(scdown) <- "type_cond"

DElist <- list()

# Running differential expression (DS vs CT for each cell type)
for (i in levels(scdown$seurat_clusters)){
  DElist[[i]] <- FindMarkers(
    scdown,
    ident.1 = paste0(i,"_DS"), # find markers for this identity
    ident.2 = paste0(i,"_CT"), # against this identity
    logfc.threshold = 0.0, # to analyse every gene
    min.pct = 0.0, # to analyse every gene
    min.cells.feature = 0
  )
}

# adding DEG and chromosome information and cell type information
DElist <- imap(DElist, function(i, ct){ # i is a df and ct is list element name
  i <- i %>% 
    # adding DEG and alteration information and gene name column
    mutate( # creates new columns
      celltype = ct,
      gene = rownames(i),
      DEG = case_when(
        avg_log2FC >= 1 & p_val_adj < 0.05 ~ "UP",
        avg_log2FC <= -1 & p_val_adj < 0.05 ~ "DOWN",
        TRUE ~ "NO" # else none of the other conditions are true
      ),
      alterat = case_when(
        avg_log2FC > 0 & p_val_adj < 0.05 ~ "UP",
        avg_log2FC < 0 & p_val_adj < 0.05 ~ "DOWN",
        TRUE ~ "NO" # else none of the other conditions are true
      )
    ) %>% 
    # adding gene info from biomart
    left_join(all_features, by = c("gene" = "Symbol"))
  return(i)
})

scdown@meta.data$type_cond <- factor(scdown@meta.data$type_cond, levels = c("Oli_DS","Oli_CT","Exc_DS","Exc_CT",
                                                                            "Inh_DS","Inh_CT","Ast_DS","Ast_CT",
                                                                            "OPC_DS","OPC_CT","Mic_DS","Mic_CT",
                                                                            "End_DS","End_CT"))

features_to_plot <- c("C1R","C1S","C4A","C4B","C5","C8G","CD46","CD59","CFI","FCN2","SLC1A2","SLC1A3")
View(DElist$Ast %>% filter(gene %in% comp_genes))

all_DEdf <- do.call(rbind, DElist)
write.table(all_DEdf, file = "supplementary_tables/allDEdf_noALRA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#################################
###     Composition Plot      ###
#################################

# All celltypes
DEdf <- do.call(rbind, DElist) %>% 
  filter(gene != "XIST")

DEdf$celltype <- factor(DEdf$celltype, levels = c("Mic","Ast","Oli","End","OPC","Inh","Exc"))

chrObj <- chromoInitiate(
  DEdf %>% dplyr::select(1:7), 
  gene_col = 'gene', 
  fc_col = 'avg_log2FC', 
  p_col = 'p_val_adj',
  celltype_col = 'celltype'
)

chrObj <- chromoComposition(chrObj, only_expr_features = T, separate_by = "celltype", score_method = 'n_DEGs')

comp_plot <- chromoCompositionPlot(chrObj)

ggsave("results_noALRA/compositionAllCelltypes.png", plot = comp_plot, height = 4, width = 4)

# Just Astrocyte
chrObj <- chromoInitiate(
  DElist$Ast %>% dplyr::select(1:7), 
  gene_col = 'gene', 
  fc_col = 'avg_log2FC', 
  p_col = 'p_val_adj'
)

chrObj <- chromoComposition(chrObj)

comp_plot <- chromoCompositionPlot(chrObj)

ggsave("results_noALRA/compositionAst.png", plot = comp_plot, height = 4, width = 10)


##############################################
###           Functional analysis          ###
##############################################

# ORA
ORA_ast <- enrichGO(
  keyType = "SYMBOL", #there is no preference for using SYMBOL, ENSEMBL or ENTREZ!
  gene = DElist$Ast %>% filter(DEG != "NO") %>% arrange(-abs(avg_log2FC)) %>% pull(gene),
  universe = rownames(scdown), # sample genes. IMPORTANT!!!
  OrgDb = org.Hs.eg.db,
  ont = "all",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

aux <- ORA_ast@result %>% 
  filter(ID %in% c('GO:0006958', 'GO:0050808', 'GO:0006956', 'GO:0019955'))

write.table(aux, file = "supplementary_tables/astORA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# ontologies of interest
int_ont <- c("GO:0098883", "GO:0150062", "GO:1905805", "GO:1905806", "GO:1905807", "GO:1905808", "GO:1905810", "GO:1905811", 
             "GO:0016322", "GO:1904799", "GO:1904800", "GO:1904801", "GO:0006958", "GO:0006956", "GO:0004875", "GO:0002430", 
             "GO:0030449", "GO:0008066", "GO:0050808", "GO:0045088", "GO:0050766", "GO:0019955")

# USAR CHROMO ORA
chrObj <- chromoInitiate( # placeholder to use chromoORAPlot()
  data.frame(Symbol = rownames(scdown), log2fc = 0, pval = 1),
  gene_col = "Symbol",
  fc_col = "log2fc",
  p_col = "pval"
)

chrObj@ora[["UP_DOWN"]][[1]] <- aux %>% 
  mutate(
    score = -log10(p.adjust),
    ONTOLOGY = 'BP'     
  )

funcenrich_plot <- chromoORAPlot(chrObj)

ggsave("results_noALRA/ORA_ast.pdf", plot = funcenrich_plot, height = 3, width = 5)


#################################
###      Package versions     ###
#################################

session_info()
