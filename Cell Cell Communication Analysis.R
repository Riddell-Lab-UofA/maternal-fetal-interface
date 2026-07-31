#integrated multinichenet

garcia <- readRDS("D:/Calgary placenta collaboration/Full integrated/garcia_rivas.rds")
garcia_subset <- subset(garcia, garcia$tissue = 'Myometrium SC', invert = T)


library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)


organism = "human"


options(timeout = 120)

if(organism == "human"){
  
  lr_network_all = 
    readRDS(url("https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
    )) %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = organism), 
      receptor = convert_alias_to_symbols(receptor, organism = organism))
  
  lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 
  
  lr_network = lr_network_all %>% 
    distinct(ligand, receptor)
  
  ligand_target_matrix = readRDS(
    url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
  
} 


Integrated <- integrated_w_labels

#add the active ident as a column in the metadata
#Labored$clusters <- Idents(object = Labored)
#see our sample identities
Idents(object = Integrated) <- "orig.ident"

levels(Integrated)

#you can change the levels back to the cell types
Idents(object = Integrated) <- "clusters"

levels(Integrated)

#now make a new column to add the different cell type groupings
garcia_subset$column <- garcia_subset$orig.ident

Non_labored <- 'Non_labored'

# for all the Non_labored samples: "25_6", "25_7"
garcia_subset$column <- garcia_subset$orig.ident
garcia_subset$column[garcia_subset$orig.ident == "Non-labored 1"] <- Non_labored
garcia_subset$column[garcia_subset$orig.ident == "Non-labored 2"] <- Non_labored


#now the other samples
Labored <- 'Labored'
garcia_subset$column[garcia_subset$orig.ident == "Labored 1"] <- Labored
garcia_subset$column[garcia_subset$orig.ident == "Labored 2"] <- Labored
garcia_subset$column[garcia_subset$orig.ident == "Labored 3"] <- Labored
garcia_subset$column[garcia_subset$orig.ident == "Labored 4"] <- Labored




sce = as.SingleCellExperiment(garcia_subset, assay = "RNA")
make.names(c("cluster.ident", "condition", "orig.ident", "cca_clusters"), unique = TRUE)

#gotta make sure all the columns are syntactically valid:
SummarizedExperiment::colData(sce)$cluster.ident = SummarizedExperiment::colData(sce)$cluster.ident %>% make.names()
SummarizedExperiment::colData(sce)$column = SummarizedExperiment::colData(sce)$column %>% make.names()
SummarizedExperiment::colData(sce)$orig.ident = SummarizedExperiment::colData(sce)$orig.ident %>% make.names()

#make sure they have the correct labels

sample_id = "orig.ident"
group_id = "column"
celltype_id = "cluster.ident"
covariates = NA
batches = NA

library(magrittr) # for the pipe %>%

sce$orig.ident <- make.names(sce$orig.ident)
sce$column <- make.names(sce$column)
sce$cluster.ident <- make.names(sce$cluster.ident)

#state what comparisons we are interested in
contrasts_oi = c("'Non_labored-Labored','Labored-Non_labored'")

contrast_tbl = tibble(contrast =
                        c("Non_labored-Labored", "Labored-Non_labored"),
                      group = c("Non_labored", "Labored"))


#define our senders and recievers of interest. For now we will just label all cells as cells of interest. We can filter this more/choose specific cell types later on
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]


#make a minimum cell number for 10 per analysis. If you don't have a lot of cells go down to min=5 cells
min_cells = 10
min_sample_prop = 0.50
fraction_cutoff = 0.05
empirical_pval = FALSE

#setting thresholds
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 

top_n_target = 250

#computer number of cores
n.cores = 32

scenario = "regular"

#to only look at DOWN regulated genes in your CONDITION OF INTEREST
#ligand_activity_down = FALSE

#to look at both UP and DOWN
ligand_activity_down = TRUE


#to run mulitnichenet output
multinichenet_output = multi_nichenet_analysis(
  sce = sce, 
  celltype_id = celltype_id, sample_id = sample_id, group_id = group_id, 
  batches = batches, covariates = covariates, 
  lr_network = lr_network, ligand_target_matrix = ligand_target_matrix, 
  contrasts_oi = contrasts_oi, contrast_tbl = contrast_tbl, 
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, 
  min_sample_prop = min_sample_prop,
  scenario = scenario, 
  ligand_activity_down = ligand_activity_down,
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj, 
  empirical_pval = empirical_pval, 
  top_n_target = top_n_target, 
  n.cores = n.cores, 
  verbose = TRUE
)

#Normalized pseudobulk expression for each cell type - sample combination

multinichenet_output$celltype_info$pb_df %>% head()

multinichenet_output$celltype_info$pb_df_group %>% head()



#DE information for each cell type - contrast combination
multinichenet_output$celltype_de %>% head()


#Output of the NicheNet ligand activity analysis, and the NicheNet ligand-target inference
multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% head()


#Tables with the final prioritization scores (results per group and per sample)
multinichenet_output$prioritization_tables$group_prioritization_tbl %>% head()


### Visualization of results

#Summarizing ChordDiagram circos plots
prioritized_tbl_oi_all = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 50, 
  rank_per_group = TRUE
)

prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(senders_receivers)) %>% magrittr::set_names(senders_receivers)

colors_receiver <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(senders_receivers)) %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)


#color_palette <- c('#FFD7D7', '#FFA7A6','#E7DCF9', '#D3EDDB','#F9D7F6','#8ca9ff','#ECCAFF','#CDD0F8','#ffd09b', '#B8DAED', '#C4E4DF', '#D1E2D0','#baffc9',  '#feffa3','#ea7b7b', '#faac68', '#afe4d6', '#a1d6b2' , '#ffcdc9',  '#fee2ad', '#ffc7a7')

#Interpretable bubble plots

group_oi = "Non_labored"

prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 50, 
  groups_oi = group_oi)

plot_oi = make_sample_lr_prod_activity_plots(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_M_50)
plot_oi


#'[RUNNING THE LIANA WRAPPER Package]

library(tidyverse)
library(SingleCellExperiment)
library(liana)
#Turning the Seurat object into SCE object
big_boss.sce <- as.SingleCellExperiment(x = big_boss)

# 2. Run LIANA on the SCE object
# Note: You must specify the column in your metadata that contains your cell type
#Wrapping and testing function
#For Cell Chat
liana_test <- liana_wrap(big_boss.sce, idents_col = 'cluster.ident', resource = 'CellChatDB')
liana_test <- liana_test %>% liana_aggregate()
liana_trunc <- liana_test %>% filter(aggregate_rank <= 0.05)
write.csv (liana_trunc, file = 'D:/Calgary placenta collaboration/liana_CellNet.csv')

liana_test <- liana_wrap(big_boss.sce, idents_col = 'cluster.ident', resource = 'CellPhoneDB')
liana_test <- liana_wrap(big_boss.sce, idents_col = 'cluster.ident', resource = 'CellTalkDB')
liana_test <- liana_wrap(big_boss.sce, idents_col = 'cluster.ident', resource = 'ICELLNET')
liana_test <- liana_test %>% liana_aggregate()
liana_trunc <- liana_test %>% filter(aggregate_rank <= 0.05)
write.csv (liana_trunc, file = 'D:/Calgary placenta collaboration/liana_CellNet.csv')



heat_freq(liana_trunc)

