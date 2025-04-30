#----------------- GENERAL PREPARATIONS ---------------------------------------
# Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(igraph)
library(GGally)
library(ggpp)
library(tidyr)
library(tibble)

aggregate_mechanism <- function(classes) {
  mechanism <- rep(NA, times = length(classes))
  mechanism[grep("aac", classes)] <- "AAC"
  mechanism[grep("aph", classes)] <- "APH"
  mechanism[grep("class_", classes)] <- "Beta-lactamase"
  mechanism[grep("erm", classes)] <- "Erm"
  mechanism[grep("mph", classes)] <- "Mph"
  mechanism[grep("mcr", classes)] <- "MCR"
  mechanism[grep("qnr", classes)] <- "Qnr"
  mechanism[grep("tet_efflux", classes)] <- "Tet efflux"
  mechanism[grep("tet_enzyme", classes)] <- "Tet enzyme"
  mechanism[grep("tet_rpg", classes)] <- "Tet RPG"
  
  return(mechanism)
}

metadata <- data.frame(fread("arg_metadata.tsv", header=TRUE)) %>% 
  tibble::column_to_rownames('V1')

metadata[is.na(metadata)] <- 0
metadata["taxonomy_score", metadata["taxonomy_score",] == "0" | metadata["taxonomy_score",] == "1" | metadata["taxonomy_score",] == "2"] <- "≤ Genus"
metadata["taxonomy_score", metadata["taxonomy_score",] == "3"] <- "Family"
metadata["taxonomy_score", metadata["taxonomy_score",] == "4"] <- "Order"
metadata["taxonomy_score", metadata["taxonomy_score",] == "5"] <- "Class"
metadata["taxonomy_score", metadata["taxonomy_score",] == "6"] <- "Phylum"

data_ww <- data.frame(fread("rarefied_counts_wastewater.tsv", header=TRUE)) %>% 
  tibble::column_to_rownames('V1')
data_ww <- data_ww[,metadata["arg_class",] != "mcr"]

data_hm <- data.frame(fread("rarefied_counts_human.tsv", header=TRUE)) %>% 
  tibble::column_to_rownames('V1') 
data_hm <- data_hm[,metadata["arg_class",] != "mcr"]

metadata <- metadata[,metadata["arg_class",] != "mcr"]

gene_class <- unique(sort(as.character(metadata["arg_class",])))

rel_ab_hm <- colSums(data_hm >= 3)/nrow(data_hm)
rel_ab_ww <- colSums(data_ww >= 3)/nrow(data_ww)

cut_off <- 0.05

hg_prom <- colnames(data_hm)[rel_ab_hm >= cut_off & rel_ab_ww < cut_off]
ww_prom <- colnames(data_hm)[rel_ab_hm < cut_off & rel_ab_ww >= cut_off]
co_prom <- colnames(data_hm)[rel_ab_hm >= cut_off & rel_ab_ww >= cut_off]
non_prom <- colnames(data_hm)[(rel_ab_hm < cut_off & rel_ab_ww < cut_off) & (rel_ab_hm > 0 | rel_ab_ww > 0)]
missing <- colnames(data_hm)[rel_ab_hm == 0 & rel_ab_ww == 0]

# Define color palettes
arg_category_palette <- c(
  "Co-promoted"='#000000',
  "HG-promoted" ="#D55E00",
  "WW-promoted" = "#0072B2",
  "Non-promoted"="#F0E442"
)

mechanism_pallete <- c(
  "AAC"='#337538',
  "APH" ="#5da899",
  "Beta-lactamase"='#dccd7d',
  "Erm"="#2e2585",
  "Mph" = "#94caec",
  "Qnr" = "#dddddd",
  "Tet efflux" = "#c26a77",
  "Tet enzyme" = "#9f4a96",
  "Tet RPG"='#7e2954'
)

phylum_palette <- c(
  "Actinomycetota"='#6A4C93',
  "Bacillota" ="#457B9D",
  "Bacteroidota" = "#F4A261",
  "Campylobacterota"="#A8DADC",
  "Pseudomonadota"="#E63946",
  "Other" = "#2A2A2A"
)

#----------------- FIGURE 1 ---------------------------------------------------

plot_data <- data.frame(human = colSums(data_hm >= 3)/nrow(data_hm), 
                        wastewater = colSums(data_ww >= 3)/nrow(data_ww), 
                        name = as.character(metadata["label",]), 
                        status = as.character(metadata["category",]),
                        class = as.character(metadata["arg_class",])) %>%
  .[.$human + .$wastewater > 0,]

plot_data$mechanism <- aggregate_mechanism(plot_data$class)

plot_data$lab <- ""
for (i in 1:nrow(plot_data)) { 
  if (plot_data$status[i] == "Established" && (plot_data$human[i] >= 0.25 | plot_data$wastewater[i] >= 0.25)) { 
    plot_data$lab[i] <- plot_data$name[i] 
  } 
}

fig_1 <- ggplot(data=plot_data, aes(x = wastewater, y = human, color = mechanism)) +
  geom_point(aes(shape = factor(status)), size = 2) +
  geom_text(aes(label=lab), color="black", size = 3, fontface = "italic") +
  theme_minimal() +
  ylab("Human gut") +
  xlab("Wastewater") +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  scale_x_continuous(labels = percent, limits = c(0,1)) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  scale_color_manual(values=mechanism_pallete, name = "Resistance mechanism") +
  geom_line(aes(x = 0.05), color="black", linewidth = 0.5) +
  geom_line(aes(y = 0.05), color="black", linewidth = 0.5) +
  theme(axis.text.x=element_text(size = rel(1.5)),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1)),
        plot.title = element_text(size = rel(3), , face = "bold"))

pdf("figure_1.pdf", height = 9, width = 10.7)
plot(fig_1)
dev.off()

#----------------- FIGURE 2a --------------------------------------------------

plot_data <- data.frame(Type = rep(c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted"), 2),
                        Status = c(rep("Established", 4), rep("Latent", 4)),
                        n = c(41, 42, 20, 229, 3, 70, 28, 1159),
                        prop = c(41/length(co_prom), 42/length(hg_prom), 20/length(ww_prom), 229/length(non_prom), 3/length(co_prom), 70/length(hg_prom), 28/length(ww_prom), 1159/length(non_prom)))

fig_2a <- ggplot(plot_data, aes(x = factor(Type, levels =c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")), y = prop, fill = Status)) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25, size=5) +
  theme_minimal() +
  ylab("Proportion") +
  xlab(NULL) +
  ggtitle("a") +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  scale_fill_manual(values=c(
    "Established"='#595959',
    "Latent" ="#ababab"
  ), name = "Status") +
  theme(axis.text.x=element_text(size = rel(1.5), angle = 45, hjust = 1, vjust = 1),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1)),
        plot.title = element_text(size = rel(3), face = "bold"),
        plot.title.position = "plot")

#----------------- FIGURE 2b --------------------------------------------------

plot_data <- data.frame(class = rep(gene_class, 8),
                        category = rep(c(rep("Co-promoted", length(gene_class)),
                                         rep("HG-promoted", length(gene_class)), 
                                         rep("WW-promoted", length(gene_class)),
                                         rep("Non-promoted", length(gene_class))), 2),
                        status = c(rep("Established", 88), rep("Latent", 88)),
                        value = NA)


co_prom_metadata <- metadata[,co_prom]
hg_prom_metadata <- metadata[,hg_prom]
ww_prom_metadata <- metadata[,ww_prom]
non_prom_metadata <- metadata[,non_prom]

for (class in gene_class){
  for (stat in c("Established", "Latent")) {
    plot_data$value[plot_data$class == class & plot_data$category == "Co-promoted" & plot_data$status == stat] <- sum(co_prom_metadata["arg_class",] == class & co_prom_metadata["category",] == stat)/ncol(co_prom_metadata[,co_prom_metadata["category",] == stat])
    plot_data$value[plot_data$class == class & plot_data$category == "HG-promoted" & plot_data$status == stat] <- sum(hg_prom_metadata["arg_class",] == class & hg_prom_metadata["category",] == stat)/ncol(hg_prom_metadata[,hg_prom_metadata["category",] == stat])
    plot_data$value[plot_data$class == class & plot_data$category == "WW-promoted" & plot_data$status == stat] <- sum(ww_prom_metadata["arg_class",] == class & ww_prom_metadata["category",] == stat)/ncol(ww_prom_metadata[,ww_prom_metadata["category",] == stat])
    plot_data$value[plot_data$class == class & plot_data$category == "Non-promoted" & plot_data$status == stat] <- sum(non_prom_metadata["arg_class",] == class & non_prom_metadata["category",] == stat)/ncol(non_prom_metadata[,non_prom_metadata["category",] == stat])
  }
}

plot_data$mechanism <- aggregate_mechanism(gene_class)

fig_2b <- ggplot(plot_data, aes(fill=mechanism, y=value, x=factor(category, levels = c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")))) + 
  geom_bar(position="fill", stat="identity") +
  theme_minimal() +
  xlab(NULL) +
  ylab("Proportion") +
  ggtitle("b") +
  scale_fill_manual(values=mechanism_pallete, name = "Resistance mechanism") +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x=element_text(size = rel(1.5), angle = 45, hjust = 1, vjust = 1),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1)),
        plot.title = element_text(size = rel(3), , face = "bold"),
        plot.title.position = "plot",
        strip.text = element_text(size = rel(1.5))) +
  facet_grid(~ status)

#----------------- FIGURE 2c --------------------------------------------------

ph_metadata <- data.frame(fread("arg_phylum_distr.tsv", header=TRUE)) %>% na.omit()
rownames(ph_metadata) <- ph_metadata$V1
ph_metadata <- ph_metadata[,-1]

major_phyla <- c("Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Pseudomonadota")
phylum_data <- ph_metadata[major_phyla,]
phylum_data["Other", ] <- as.numeric(colSums(ph_metadata[!(rownames(ph_metadata) %in% major_phyla),]) > 1)

plot_phyla <- function(gene_ids, metadata, phylum_data, title, tag) {
  
  subset_est <- phylum_data[,gene_ids[metadata["category" ,gene_ids] == "Established"]]
  subset_lat <- phylum_data[,gene_ids[metadata["category", gene_ids] == "Latent"]]
  
  plot_data <- data.frame(prop=c(rowSums(subset_est)/ncol(subset_est), rowSums(subset_lat)/ncol(subset_lat)),
                          n=c(rowSums(subset_est), rowSums(subset_lat)))
  plot_data$Phylum <- factor(rep(c("Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Pseudomonadota", "Other"), 2), 
                             levels = c("Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Pseudomonadota", "Other"))
  plot_data$Status <- factor(c(rep(paste(c("Established (n=", ncol(subset_est), ")"), collapse=""), 6), rep(paste(c("Latent (n=", ncol(subset_lat), ")"), collapse=""), 6))) 
  
  p <- ggplot(data=plot_data, aes(x = Phylum, y = prop, fill = Phylum, alpha = Status)) +
    geom_bar(position="dodge", stat="identity") +
    geom_text(aes(label=n), position=position_dodge(width=0.9), hjust=-0.5, size=5) +
    theme_minimal() +
    ylab("Proportion") +
    xlab(NULL) +
    ggtitle(title) +
    labs(tag = tag) +
    scale_fill_manual(values=phylum_palette, name = "Phylum") +
    coord_flip() +
    scale_alpha_discrete(range=c(1, 0.5)) +
    scale_y_continuous(labels = percent, limits = c(0,1)) +
    scale_x_discrete(limits = rev(levels(plot_data$Phylum))) +
    theme(axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
          axis.text.y=element_text(size = rel(1.3)),
          axis.title=element_text(size=rel(1.3)),
          legend.title=element_text(size=rel(1.2)), 
          legend.text=element_text(size=rel(1)),
          plot.tag = element_text(size=rel(2.5), face = "bold"),
          plot.title = element_text(size = rel(1.5)),
          plot.margin = margin(t = 10,
                               r = 30,
                               b = 10,
                               l = 20))
  
  return(p)
}

fig_2c <- plot_phyla(co_prom, metadata, phylum_data, "Co-promoted", "c")
fig_2d<- plot_phyla(hg_prom, metadata, phylum_data, "HG-promoted", "d")
fig_2e <- plot_phyla(ww_prom, metadata, phylum_data, "WW-promoted", "e")
fig_2f <- plot_phyla(non_prom, metadata, phylum_data, "Non-promoted", "f")


pdf("figure_2.pdf", height = 15, width = 15)
(fig_2a | fig_2b) /
  (fig_2c | fig_2d) /
  (fig_2e | fig_2f) 
dev.off()

#----------------- FIGURE 3a-b ------------------------------------------------

perform_permutation <- function(df, n_iter, n_subsample) {
  # Pre-allocate list to collect results
  results_list <- vector("list", length = n_iter * length(unique(df$taxonomy_score)))
  
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  taxonomy_scores <- unique(df$taxonomy_score)
  group_levels <- c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")
  
  result_idx <- 1
  
  for (i in seq_len(n_iter)) {
    df$grp <- sample(df$grp)  # scramble groups once per iteration
    
    for (score in taxonomy_scores) {
      tmp_mat <- matrix(0, nrow = n_iter, ncol = 4)  # preallocate temporary results
      
      for (j in seq_len(n_iter)) {
        # Sample once per group
        subsample_indices <- unlist(
          lapply(group_levels, function(g) {
            sample(which(df$grp == g), n_subsample, replace = FALSE)
          })
        )
        subsample <- df[subsample_indices, ]
        
        n_args <- sum(subsample$taxonomy_score == score)
        
        if (n_args == 0) {
          next
        }
        
        # Calculate fraction of each group matching score
        group_counts <- sapply(group_levels, function(g) {
          sum(subsample$taxonomy_score[subsample$grp == g] == score)
        })
        
        tmp_mat[j, ] <- group_counts / n_args
      }
      
      # Take column means and store result
      tmp_means <- colMeans(tmp_mat, na.rm = TRUE)
      results_list[[result_idx]] <- c(tmp_means, score)
      result_idx <- result_idx + 1
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Combine list into data.frame
  output <- do.call(rbind, results_list)
  colnames(output) <- c("co_prom", "hg_prom", "ww_prom", "non_prom", "score")
  
  output <- as.data.frame(output)
  return(output)
}


plot_tax_score <- function(metadata, status, n_iter, tag) {
  co_prom_metadata <- metadata[,co_prom] %>% select(which(.["category", ] == status))
  hg_prom_metadata <- metadata[,hg_prom] %>% select(which(.["category", ] == status))
  ww_prom_metadata <- metadata[,ww_prom] %>% select(which(.["category", ] == status))
  non_prom_metadata <- metadata[,non_prom] %>% select(which(.["category", ] == status))
  
  if (status == "Established") {
    plot_data <- data.frame(taxonomy_score = rep(unique(as.character(metadata["taxonomy_score",])), 4),
                            category = c(rep(paste(c("Co-promoted (n=", ncol(co_prom_metadata), ")"), collapse = ""), 5), 
                                         rep(paste(c("HG-promoted (n=", ncol(hg_prom_metadata), ")"), collapse = ""), 5), 
                                         rep(paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = ""), 5), 
                                         rep(paste(c("Non-promoted (n=", ncol(non_prom_metadata), ")"), collapse = ""), 5)),
                            value = NA)
    
    status_metadata <- metadata[,c(co_prom, hg_prom, ww_prom, non_prom)] %>% select(which(.["category", ] == status))
    
    n_subsample <- min(ncol(co_prom_metadata), ncol(hg_prom_metadata), ncol(ww_prom_metadata), ncol(non_prom_metadata))
  }
  
  else if (status == "Latent") {
    plot_data <- data.frame(taxonomy_score = rep(unique(as.character(metadata["taxonomy_score",])), 3),
                            category = c(rep(paste(c("HG-promoted (n=", ncol(hg_prom_metadata), ")"), collapse = ""), 5), 
                                         rep(paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = ""), 5), 
                                         rep(paste(c("Non-promoted (n=", ncol(non_prom_metadata), ")"), collapse = ""), 5)),
                            value = NA)
    
    status_metadata <- metadata[,c(hg_prom, ww_prom, non_prom)] %>% select(which(.["category", ] == status))
    
    n_subsample <- min(ncol(hg_prom_metadata), ncol(ww_prom_metadata), ncol(non_prom_metadata))
  }
  
  for (score in unique(plot_data$taxonomy_score)) {
    prop_co_prom <- c()
    prop_hg_prom <- c()
    prop_ww_prom <- c()
    prop_non_prom <- c()
    
    for (i in 1:n_iter) {
      if (status == "Established") {
        subsample_co_prom <- co_prom_metadata[,sample(colnames(co_prom_metadata), n_subsample, replace=FALSE)]
      }
      
      subsample_hg_prom <- hg_prom_metadata[,sample(colnames(hg_prom_metadata), n_subsample, replace=FALSE)]
      subsample_ww_prom <- ww_prom_metadata[,sample(colnames(ww_prom_metadata), n_subsample, replace=FALSE)]
      subsample_non_prom <- non_prom_metadata[,sample(colnames(non_prom_metadata), n_subsample, replace=FALSE)]
      
      if (status == "Established") {
        subsample_comb <- cbind(subsample_co_prom, subsample_hg_prom, subsample_ww_prom, subsample_non_prom)
      }
      else if (status == "Latent") {
        subsample_comb <- cbind(subsample_hg_prom, subsample_ww_prom, subsample_non_prom)
      }
      
      n_args <- sum(subsample_comb["taxonomy_score",] == score)
      
      if (n_args == 0) {
        next
      }
      else {
        if (status == "Established") {
          prop_co_prom <- c(prop_co_prom, sum(as.character(subsample_co_prom["taxonomy_score",]) == score)/n_args)
        }
        
        prop_hg_prom <- c(prop_hg_prom, sum(as.character(subsample_hg_prom["taxonomy_score",]) == score)/n_args)
        prop_ww_prom <- c(prop_ww_prom, sum(as.character(subsample_ww_prom["taxonomy_score",]) == score)/n_args)
        prop_non_prom <- c(prop_non_prom, sum(as.character(subsample_non_prom["taxonomy_score",]) == score)/n_args)
      }
      
    }
    
    
    if (status == "Established") {
      plot_data$value[plot_data$taxonomy_score == score & plot_data$category == paste(c("Co-promoted (n=", ncol(co_prom_metadata), ")"), collapse = "")] <- mean(prop_co_prom)
    }
    
    plot_data$value[plot_data$taxonomy_score == score & plot_data$category == paste(c("HG-promoted (n=", ncol(hg_prom_metadata), ")"), collapse = "")] <- mean(prop_hg_prom)
    plot_data$value[plot_data$taxonomy_score == score & plot_data$category == paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = "")] <- mean(prop_ww_prom)
    plot_data$value[plot_data$taxonomy_score == score & plot_data$category == paste(c("Non-promoted (n=", ncol(non_prom_metadata), ")"), collapse = "")] <- mean(prop_non_prom)
  }
  
  plot_data$category <- factor(plot_data$category, levels = c(paste(c("Co-promoted (n=", ncol(co_prom_metadata), ")"), collapse = ""), 
                                                              paste(c("HG-promoted (n=", ncol(hg_prom_metadata), ")"), collapse = ""), 
                                                              paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = ""), 
                                                              paste(c("Non-promoted (n=", ncol(non_prom_metadata), ")"), collapse = "")))
  
  p <- ggplot(data=plot_data, aes(x = factor(taxonomy_score, levels=c("Phylum", "Class", "Order", "Family", expression("≤ Genus"))), y = value, fill = category)) +
    geom_bar(position="dodge", stat="identity") +
    theme_minimal() +
    ylab("Proportion") +
    xlab(NULL) +
    ggtitle(status) +
    labs(tag = tag) +
    scale_x_discrete(labels=c(expression("Phylum", "Class", "Order", "Family", ""<="Genus"))) +
    scale_y_continuous(labels = percent, limits = c(0,1)) +
    scale_fill_manual(
      values = setNames(
        c("#000000", "#D55E00", "#0072B2", "#F0E442"),
        c(paste(c("Co-promoted (n=", ncol(co_prom_metadata), ")"), collapse = ""),
          paste(c("HG-promoted (n=", ncol(hg_prom_metadata), ")"), collapse = ""),
          paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = ""),
          paste(c("Non-promoted (n=", ncol(non_prom_metadata), ")"), collapse = ""))
      ), name = ""
    ) +
    theme(axis.text.x=element_text(size = rel(1.3), angle = 45, vjust = 1, hjust = 1),
          axis.text.y=element_text(size = rel(1.3)),
          axis.title=element_text(size=rel(1.3)),
          legend.title=element_text(size=rel(1.2)), 
          legend.text=element_text(size=rel(1)),
          plot.tag = element_text(size=rel(2.5), face = "bold"),
          plot.title = element_text(size = rel(1.5)),
          plot.margin = margin(t = 10,
                               r = 40,
                               b = 10,
                               l = 40))
  
  if (status == "Established") {
    print("Performing permutation test")
    permutation_df <- cbind(co_prom_metadata, 
                            hg_prom_metadata, 
                            ww_prom_metadata, 
                            non_prom_metadata) %>% 
      t() %>% 
      as.data.frame()
    
    permutation_df$grp <- c(rep("Co-promoted", ncol(co_prom_metadata)), 
                            rep("HG-promoted", ncol(hg_prom_metadata)), 
                            rep("WW-promoted", ncol(ww_prom_metadata)), 
                            rep("Non-promoted", ncol(non_prom_metadata)))
    
    permutation_results <- perform_permutation(permutation_df, n_iter, n_subsample)
    
    p_val_co_prom_phylum <- sum(as.numeric(permutation_results$co_prom[permutation_results$score == "Phylum"]) > plot_data$value[plot_data$category == paste(c("Co-promoted (n=", ncol(co_prom_metadata), ")"), collapse = "") & plot_data$taxonomy_score == "Phylum"]) / (n_iter+1)
    p_val_ww_prom_class <- sum(as.numeric(permutation_results$ww_prom[permutation_results$score == "Class"]) > plot_data$value[plot_data$category == paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = "") & plot_data$taxonomy_score == "Class"]) / (n_iter+1)
    
    print(paste(c("Co-promoted phylum-level permutation p-value:", p_val_co_prom_phylum), collapse = " "))
    print(paste(c("WW-promoted phylum-level permutation p-value:", p_val_ww_prom_class), collapse = " "))
  }
  return(p)
}

fig_3a <- plot_tax_score(metadata, "Established", 1000, "a")
fig_3b <- plot_tax_score(metadata, "Latent", 1000, "b")

#----------------- FIGURE 3c-f ------------------------------------------------

extract_phylum_comb <- function(phyla, gene, df) {
  count <- as.numeric(df[phyla[1], gene] > 0 & df[phyla[2], gene] > 0)
  
  return(count)
}

plot_phylum_combinations <- function(gene_ids, df, title, tag) {
  phylum_comb <- tidyr::crossing(x = rownames(df), y = rownames(df)) %>% data.frame()
  phylum_comb$count <- 0
  
  for (gene in gene_ids) {
    hits <- apply(phylum_comb, 1, extract_phylum_comb, gene=gene, df = phylum_data)
    phylum_comb$count <- phylum_comb$count + hits
  }
  
  phylum_comb_filtered <- phylum_comb %>%
    filter(as.numeric(factor(x, levels = c("Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Pseudomonadota", "Other"))) <=
             as.numeric(factor(y, levels = c("Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Pseudomonadota", "Other")))) 
  
  nsize <- phylum_comb_filtered$count[phylum_comb_filtered$x == phylum_comb_filtered$y]
  names(nsize) <- phylum_comb_filtered$x[phylum_comb_filtered$x == phylum_comb_filtered$y]
  
  phylum_comb_filtered <- phylum_comb_filtered[phylum_comb_filtered$x != phylum_comb_filtered$y,]
  phylum_comb_filtered$count[phylum_comb_filtered$count == 0] <- 0.001
  phylum_comb_filtered$prop <- phylum_comb_filtered$count/length(gene_ids)
  
  g <- graph_from_data_frame(d = phylum_comb_filtered, directed = FALSE)
  nsize <- nsize[match(V(g)$name, names(nsize))]
  
  E(g)$lty <- rep(1, nrow(phylum_comb_filtered))
  E(g)$lty[phylum_comb_filtered$prop < 0.10] <- 0
  
  node_color <- phylum_palette
  node_color[names(nsize)[nsize == 0]] <- "#ffffff"
  
  set.seed(123)
  
  p <- ggnet2(g, mode = "kamadakawai", color=V(g)$name, palette = node_color, alpha = 1,
              size = nsize, legend.size = 9, layout.exp = 0.1, edge.lty = E(g)$lty, 
              edge.size = 4*phylum_comb_filtered$prop, edge.alpha = 0.6) +
    guides(size = "none") +
    ggtitle(title) +
    labs(tag = tag) +
    theme(legend.title=element_text(size=rel(1.2)), 
          legend.text=element_text(size=rel(1)),
          plot.tag = element_text(size=rel(2.5), face = "bold"),
          plot.title = element_text(size = rel(1.5)),
          plot.margin = margin(t = 10,
                               r = 30,
                               b = 10,
                               l = 30))
  
  return(p)
}

ph_metadata <- data.frame(fread("arg_phylum_distr.tsv", header=TRUE)) %>% na.omit()
rownames(ph_metadata) <- ph_metadata$V1
ph_metadata <- ph_metadata[,-1]

major_phyla <- c("Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Pseudomonadota")
phylum_data <- ph_metadata[major_phyla,]
phylum_data["Other", ] <- as.numeric(colSums(ph_metadata[!(rownames(ph_metadata) %in% major_phyla),]) > 1)

inter_phylum_cop <- colnames(metadata[,colnames(metadata) %in% co_prom & metadata["taxonomy_score",] == "Phylum"])
inter_phylum_hgp <- colnames(metadata[,colnames(metadata) %in% hg_prom & metadata["taxonomy_score",] == "Phylum"])
inter_phylum_wwp <- colnames(metadata[,colnames(metadata) %in% ww_prom & metadata["taxonomy_score",] == "Phylum"])
inter_phylum_nop <- colnames(metadata[,colnames(metadata) %in% non_prom & metadata["taxonomy_score",] == "Phylum"])

fig_3c <- plot_phylum_combinations(inter_phylum_cop, phylum_data, "Co-promoted (n=20)", "c")
fig_3d <- plot_phylum_combinations(inter_phylum_hgp, phylum_data, "HG-promoted (n=17)", "d")
fig_3e <- plot_phylum_combinations(inter_phylum_wwp, phylum_data, "WW-promoted (n=4)", "e")
fig_3f <- plot_phylum_combinations(inter_phylum_nop, phylum_data, "Non-promoted (n=28)", "f")


pdf("figure_3.pdf", height = 14, width = 14)
(fig_3a | fig_3b) /
  (fig_3c | fig_3d) /
  (fig_3e | fig_3f)
dev.off()

#----------------- FIGURE 4 ---------------------------------------------------

read_mges <- function(filename, gene_ids, metadata, mge_type, arg_type, status) {
  df <- data.frame(fread(filename)) %>% 
    tibble::column_to_rownames('V1') %>%
    select(subset = -ncol(.)) 
  
  metadata_subset <- metadata[,gene_ids]
  df_subset <- df[,metadata_subset["category",] == status]
  
  output <- data.frame(prop = rowSums(df_subset[mge_type$subset != "arg",]) / ncol(df_subset), type = mge_type$subset[mge_type$subset != "arg"])
  
  args <- df_subset[mge_type$subset == "arg",]
  arg_summary <- data.frame(prop = rep(NA, 6), type = rep("arg", 6), row.names=unique(arg_type$V2))
  
  for (i in 1:nrow(arg_summary)) {
    arg_summary$prop[i] <- sum(colSums(df_subset[arg_type$V1[arg_type$V2 == rownames(arg_summary)[i]],]) > 0) / ncol(df_subset)
  }
  
  output <- rbind(output, arg_summary) %>% 
    filter(.$type != "is")
  
  output$status <- status
  output$n <- ncol(df_subset)
  
  return(output)
}

mge_type <- data.frame(fread("co_prom_mges.tsv")) %>% 
  tibble::column_to_rownames('V1') %>%
  select(subset = ncol(.))

arg_type <- data.frame(fread("arg_types.tsv", header=FALSE))

est_co_prom_mge <- read_mges("co_prom_mges.tsv", co_prom, metadata, mge_type, arg_type, "Established")
est_hg_prom_mge <- read_mges("hg_prom_mges.tsv", hg_prom, metadata, mge_type, arg_type, "Established")
est_ww_prom_mge <- read_mges("ww_prom_mges.tsv", ww_prom, metadata, mge_type, arg_type, "Established")
est_non_prom_mge <- read_mges("non_prom_mges.tsv", non_prom, metadata, mge_type, arg_type, "Established")

mge_type <- mge_type[mge_type$subset != "is",]

plot_data <- data.frame(element=rep(rownames(est_co_prom_mge),4),
                        type=rep(est_co_prom_mge$type,4),
                        value=c(est_co_prom_mge$prop, 
                                est_hg_prom_mge$prop, 
                                est_ww_prom_mge$prop, 
                                est_non_prom_mge$prop),
                        category=c(rep("Co-promoted", length(est_co_prom_mge$prop)),
                                   rep("HG-promoted", length(est_co_prom_mge$prop)), 
                                   rep("WW-promoted", length(est_co_prom_mge$prop)),
                                   rep("Non-promoted", length(est_co_prom_mge$prop))),
                        status=c(est_co_prom_mge$status, 
                                 est_hg_prom_mge$status, 
                                 est_ww_prom_mge$status, 
                                 est_non_prom_mge$status),
                        n=c(est_co_prom_mge$n, 
                            est_hg_prom_mge$n, 
                            est_ww_prom_mge$n, 
                            est_non_prom_mge$n)) %>%
  filter(type == "mpf")

plot_data$category <- paste0(plot_data$category, "\n(n=", plot_data$n, ")")

fig_4 <- ggplot(
  plot_data, 
  aes(
    y=value, 
    x=factor(category, levels=c("Non-promoted\n(n=229)", "WW-promoted\n(n=20)", "HG-promoted\n(n=42)", "Co-promoted\n(n=41)")), 
    fill=factor(element, levels = rev(unique(element))))
) + 
  geom_bar(position="dodge", stat="identity")  +
  theme_minimal() +
  ylab("Proportion")  +
  xlab(NULL) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  scale_fill_manual(
    values=c(
      "B"='#337538',
      "F" ="#dccd7d",
      "FA" = "#2e2585",
      "FATA"="#94caec",
      "G" = "#dddddd",
      "I"="#c26a77",
      "T" = "#9f4a96"
    ), 
    name = "MPF type",
    breaks = unique(plot_data$element)
  ) +
  theme(
    axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
    axis.text.y=element_text(size = rel(1.3)),
    axis.title=element_text(size=rel(1.3)),
    legend.title=element_text(size=rel(1.2)), 
    legend.text=element_text(size=rel(1)),
    plot.tag = element_text(size=rel(2.5), face = "bold"),
    plot.title = element_text(size = rel(1.5))
  ) +
  coord_flip()

pdf("figure_4.pdf", height = 5, width = 7)
plot(fig_4)
dev.off()

plot_data$value <- plot_data$value*plot_data$n
plot_data$category <- gsub("\\\n.*", "", plot_data$category)

contingency <- plot_data %>%
  select(element, category, value) %>%
  pivot_wider(names_from = category, values_from = value, values_fill = 0) %>%
  column_to_rownames("element")

#----------------- FIGURE 5a-b ------------------------------------------------

plot_n_pathogens <- function(metadata, status, n_iter, tag) {
  co_prom_metadata <- metadata[,co_prom] %>% select(which(.["category", ] == status))
  hg_prom_metadata <- metadata[,hg_prom] %>% select(which(.["category", ] == status))
  ww_prom_metadata <- metadata[,ww_prom] %>% select(which(.["category", ] == status))
  non_prom_metadata <- metadata[,non_prom] %>% select(which(.["category", ] == status))
  
  if (status == "Established") {
    plot_data <- data.frame(n_pathogens = rep(0:max(as.numeric(metadata["n_pathogens", metadata["category",] == status])), 4),
                            category = c(rep(paste(c("Co-promoted (n=", ncol(co_prom_metadata), ")"), collapse = ""), max(as.numeric(metadata["n_pathogens", metadata["category",] == status]))+1), 
                                         rep(paste(c("HG-promoted (n=", ncol(hg_prom_metadata), ")"), collapse = ""), max(as.numeric(metadata["n_pathogens", metadata["category",] == status]))+1), 
                                         rep(paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = ""), max(as.numeric(metadata["n_pathogens", metadata["category",] == status]))+1), 
                                         rep(paste(c("Non-promoted (n=", ncol(non_prom_metadata), ")"), collapse = ""), max(as.numeric(metadata["n_pathogens", metadata["category",] == status]))+1)),
                            value = NA)
  }
  
  else if (status == "Latent") {
    plot_data <- data.frame(n_pathogens = rep(0:max(as.numeric(metadata["n_pathogens", metadata["category",] == status])), 3),
                            category = c(rep(paste(c("HG-promoted (n=", ncol(hg_prom_metadata), ")"), collapse = ""), max(as.numeric(metadata["n_pathogens", metadata["category",] == status]))+1), 
                                         rep(paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = ""), max(as.numeric(metadata["n_pathogens", metadata["category",] == status]))+1), 
                                         rep(paste(c("Non-promoted (n=", ncol(non_prom_metadata), ")"), collapse = ""), max(as.numeric(metadata["n_pathogens", metadata["category",] == status]))+1)),
                            value = NA)
  }
  
  for (n in unique(plot_data$n_pathogens)) {
    if (status == "Established") {
      n_subsample <- min(ncol(co_prom_metadata), ncol(hg_prom_metadata), ncol(ww_prom_metadata), ncol(non_prom_metadata))
    }
    else if (status == "Latent") {
      n_subsample <- min(ncol(hg_prom_metadata), ncol(ww_prom_metadata), ncol(non_prom_metadata))
    }
    
    prop_co_prom <- c()
    prop_hg_prom <- c()
    prop_ww_prom <- c()
    prop_non_prom <- c()
    
    for (i in 1:n_iter) {
      if (status == "Established") {
        subsample_co_prom <- co_prom_metadata[,sample(colnames(co_prom_metadata), n_subsample, replace=FALSE)]
      }
      
      subsample_hg_prom <- hg_prom_metadata[,sample(colnames(hg_prom_metadata), n_subsample, replace=FALSE)]
      subsample_ww_prom <- ww_prom_metadata[,sample(colnames(ww_prom_metadata), n_subsample, replace=FALSE)]
      subsample_non_prom <- non_prom_metadata[,sample(colnames(non_prom_metadata), n_subsample, replace=FALSE)]
      
      if (status == "Established") {
        subsample_comb <- cbind(subsample_co_prom, subsample_hg_prom, subsample_ww_prom, subsample_non_prom)
      }
      else if (status == "Latent") {
        subsample_comb <- cbind(subsample_hg_prom, subsample_ww_prom, subsample_non_prom)
      }
      
      n_args <- sum(as.numeric(subsample_comb["n_pathogens",]) >= n)
      
      if (n_args == 0) {
        next
      }
      
      else {
        if (status == "Established") {
          prop_co_prom <- c(prop_co_prom, sum(as.numeric(subsample_co_prom["n_pathogens",]) >= n)/n_args)
        }
        
        prop_hg_prom <- c(prop_hg_prom, sum(as.numeric(subsample_hg_prom["n_pathogens",]) >= n)/n_args)
        prop_ww_prom <- c(prop_ww_prom, sum(as.numeric(subsample_ww_prom["n_pathogens",]) >= n)/n_args)
        prop_non_prom <- c(prop_non_prom, sum(as.numeric(subsample_non_prom["n_pathogens",]) >= n)/n_args)
      }
    }
    
    if (status == "Established") {
      plot_data$value[plot_data$n_pathogens == n & plot_data$category == paste(c("Co-promoted (n=", ncol(co_prom_metadata), ")"), collapse = "")] <- mean(prop_co_prom)
    }
    
    plot_data$value[plot_data$n_pathogens == n & plot_data$category == paste(c("HG-promoted (n=", ncol(hg_prom_metadata), ")"), collapse = "")] <- mean(prop_hg_prom)
    plot_data$value[plot_data$n_pathogens == n & plot_data$category == paste(c("WW-promoted (n=", ncol(ww_prom_metadata), ")"), collapse = "")] <- mean(prop_ww_prom)
    plot_data$value[plot_data$n_pathogens == n & plot_data$category == paste(c("Non-promoted (n=", ncol(non_prom_metadata), ")"), collapse = "")] <- mean(prop_non_prom)
  }
  
  levels_labels <- c(
    paste0("Co-promoted (n=", ncol(co_prom_metadata), ")"),
    paste0("HG-promoted (n=", ncol(hg_prom_metadata), ")"),
    paste0("WW-promoted (n=", ncol(ww_prom_metadata), ")"),
    paste0("Non-promoted (n=", ncol(non_prom_metadata), ")")
  )
  
  p <- ggplot(data=plot_data, aes(x = n_pathogens, y = value, color = factor(category, levels =  levels_labels))) +
    geom_line(linewidth=1) +
    theme_minimal() +
    scale_x_continuous(breaks = 0:15) +
    scale_y_continuous(labels = percent) +
    ylab("Proportion") +
    xlab(expression("Number of Pathogens (" >= ")")) +
    ggtitle(status) +
    labs(tag = tag) +
    scale_color_manual(
      values = setNames(
        c("#000000", "#D55E00", "#0072B2", "#F0E442"),
        levels_labels
      ),
      name = ""
    ) +
    theme(axis.text.x=element_text(size = rel(1.5)),
          axis.text.y=element_text(size = rel(1.3)),
          axis.title=element_text(size=rel(1.3)),
          legend.title=element_text(size=rel(1.2)), 
          legend.text=element_text(size=rel(1)),
          plot.tag = element_text(size=rel(2.5), face = "bold"),
          plot.title = element_text(size = rel(1.5)))
  
  return(p)
}

fig_5a <- plot_n_pathogens(metadata, "Established", 1000, "a")
fig_5b <- plot_n_pathogens(metadata, "Latent", 1000, "b")

#----------------- FIGURE 5c-f ------------------------------------------------

pathogen_data <- data.frame(fread("arg_pathogen_distr.tsv", header=TRUE)) %>% 
  tibble::column_to_rownames('V1')

taxonomy <- data.frame(fread("../../HGT inference project/R files/assembly_id_w_full_taxonomy.txt", header=FALSE)) %>%
  select(-1) %>% unique() %>% na.omit()
rownames(taxonomy) <- taxonomy$V8
taxonomy <- taxonomy[rownames(pathogen_data),]

plot_pathogens <- function(gene_ids, pathogen_data, metadata, taxonomy, title, tag) {
  subset_est <- pathogen_data[,gene_ids[metadata["category" ,gene_ids] == "Established"]]
  subset_lat <- pathogen_data[,gene_ids[metadata["category", gene_ids] == "Latent"]]
  
  plot_data <- data.frame(n=c(rowSums(subset_est)/ncol(subset_est), rowSums(subset_lat)/ncol(subset_lat)),
                          Status = factor(c(rep(paste(c("Established (n=", ncol(subset_est), ")"), collapse=""), 15), rep(paste(c("Latent (n=", ncol(subset_lat), ")"), collapse=""), 15))),
                          Pathogen = rep(rownames(pathogen_data),2))
  
  plot_data$Phylum <- taxonomy[plot_data$Pathogen,2]
  
  p <- ggplot(data=plot_data, aes(x = Pathogen, y = n, fill = Phylum, alpha = Status)) +
    geom_bar(position="dodge", stat="identity") +
    theme_minimal() +
    ylab("Proportion") +
    xlab(NULL) +
    ggtitle(title) +
    labs(tag = tag) +
    scale_fill_manual(values=phylum_palette, name = "Phylum") +
    scale_alpha_discrete(range=c(1, 0.5)) +
    scale_y_continuous(labels = percent, limits = c(0,0.75)) +
    theme(axis.text.x=element_text(size = rel(1.3), angle = 45, vjust = 1, hjust = 1, face = "italic"),
          axis.text.y=element_text(size = rel(1.3)),
          axis.title=element_text(size=rel(1.3)),
          legend.title=element_text(size=rel(1.2)), 
          legend.text=element_text(size=rel(1)),
          plot.tag = element_text(size=rel(2.5), face = "bold"),
          plot.title = element_text(size = rel(1.5)),
          plot.margin = ggplot2::margin(t = 20,
                                        r = 10,
                                        b = 20,
                                        l = 60))
  
  return(p)
}

fig_5c <- plot_pathogens(co_prom, pathogen_data, metadata, taxonomy, "Co-promoted", "c")
fig_5d <- plot_pathogens(hg_prom, pathogen_data, metadata, taxonomy, "HG-promoted", "d")
fig_5e <- plot_pathogens(ww_prom, pathogen_data, metadata, taxonomy, "WW-promoted", "e")
fig_5f <- plot_pathogens(non_prom, pathogen_data, metadata, taxonomy, "Non-promoted", "f")

pdf("figure_5.pdf", height = 15, width = 16)
(fig_5a | fig_5b) /
  (fig_5c | fig_5d) /
  (fig_5e | fig_5f)
dev.off()

#----------------- FIGURE 6 ---------------------------------------------------

read_genetic_compatibility <- function(filename, status, category) {
  df <- data.frame(fread(filename)) %>% 
    tibble::column_to_rownames('X0') 
  
  status_ids <- metadata[, colnames(df)] %>% 
    select(which(.["category", ] == status)) %>%
    colnames()
  
  df <- df[,status_ids]
  
  value <- c()
  species <- c()
  for (sp in rownames(df)) {
    value <- c(value, as.numeric(df[sp,]))
    species <- c(species, rep(sp, ncol(df)))
  }
  
  output <- data.frame(Value = value, 
                       Species = species, 
                       Category = category)
  
  return(output)
}

compile_pathogen_comp <- function(filename, gene_ids, pathogen_data, category, status) {
  subset <- pathogen_data[,gene_ids[metadata["category" ,gene_ids] == status]]
  gcomp <- read_genetic_compatibility(filename, status, category)
  
  df <- data.frame(prop=rowSums(subset)/ncol(subset),
                   Category = category,
                   Pathogen = rownames(pathogen_data),
                   median_distance = NA)
  
  for (i in 1:nrow(pathogen_data)) {
    df$median_distance[i] <- median(gcomp$Value[gcomp$Species == rownames(pathogen_data)[i]])
  }
  
  return(df)
  
}

pathogen_data <- data.frame(fread("arg_pathogen_distr.tsv", header=TRUE)) %>% 
  tibble::column_to_rownames('V1')

df_co_prom <- compile_pathogen_comp("gcomp_co_prom.tsv", co_prom, pathogen_data, "Co-promoted", "Established")
df_hg_prom <- compile_pathogen_comp("gcomp_hg_prom.tsv", hg_prom, pathogen_data, "HG-promoted", "Established")
df_ww_prom <- compile_pathogen_comp("gcomp_ww_prom.tsv", ww_prom, pathogen_data, "WW-promoted", "Established")
df_non_prom <- compile_pathogen_comp("gcomp_non_prom.tsv", non_prom, pathogen_data, "Non-promoted", "Established")

plot_data <- rbind(df_co_prom, df_hg_prom, df_ww_prom, df_non_prom)

fig_6b <- ggplot(data=plot_data, aes(x=median_distance, y = prop, color=factor(Category, levels = c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")))) +
  annotate("segment", x=min(plot_data$median_distance)*0.95, xend=max(plot_data$median_distance)*1.05, y=min(plot_data$prop), yend=min(plot_data$prop),
           arrow=arrow(length=unit(0.3, "cm"), type="closed")) + 
  annotate("segment", x=min(plot_data$median_distance)*0.95, xend=min(plot_data$median_distance)*0.95, y=min(plot_data$prop), yend=max(plot_data$prop)*1.05,
          arrow=arrow(length=unit(0.3, "cm"), type="closed")) +
  geom_point(size=3) +
  ylab("Prevalence in pathogen") +
  xlab("Genetic incompatibility with pathogen\n(distance between 5mer distributions)") +
  theme_minimal() +
  labs(tag="b") +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values=arg_category_palette, name = "") +
  theme(axis.text.x=element_text(size = rel(1.3)),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1.2)),
        plot.tag = element_text(size=rel(2.5), face = "bold"),
        plot.title = element_text(size = rel(1.5))) 

fig_6a <- ggplot(plot_data, aes(x = factor(Category, levels = c("Co-promoted", "WW-promoted", "HG-promoted", "Non-promoted")),
                                y = median_distance,
                                fill = factor(Category, levels = c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")))) +
  geom_boxplot() +
  scale_fill_manual(values = arg_category_palette,
                    name = "", guide="none") +
  theme_minimal() +
  labs(x = "", y = "Genetic incompatibility\n(distance between 5mer distributions)", tag = "a") +
  theme(axis.text.x=element_text(size = rel(1.3)),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1)),
        plot.tag = element_text(size=rel(2.5), face = "bold"),
        plot.title = element_text(size = rel(1.5)))

pdf("figure_6.pdf", height = 6, width = 14)
(fig_6a | fig_6b) 
dev.off()

#----------------- FIGURE S1a-b -----------------------------------------------

norm_data_hm <- data.frame(fread("normalized_counts_human.tsv", header=TRUE)) %>% 
  tibble::column_to_rownames('V1')
norm_data_hm[norm_data_hm == 0] <- NA

norm_data_ww <- data.frame(fread("normalized_counts_wastewater.tsv", header=TRUE)) %>% 
  tibble::column_to_rownames('V1')
norm_data_ww[norm_data_ww == 0] <- NA

hg_prom_hg <- norm_data_hm[,hg_prom]
hg_prom_ww <- norm_data_ww[,hg_prom]
hg_prom_metadata <- metadata[,hg_prom]

ww_prom_hg <- norm_data_hm[,ww_prom]
ww_prom_ww <- norm_data_ww[,ww_prom]
ww_prom_metadata <- metadata[,ww_prom]

manhattan <- function(counts, metadata, title) {
  counts_reshape <- c()
  label <- c()
  class <- c()
  gname <- c()
  
  for (i in 1:ncol(counts)) {
    counts_reshape <- c(counts_reshape, counts[,i])
    label <- c(label, rep(colnames(counts)[i], nrow(counts)))
    class <- c(class, rep(metadata["arg_class",i], nrow(counts)))
    
    if (metadata["category",i] == "Established") {
      gname <- c(gname, rep(metadata["label",i], nrow(counts)))
    }
    else {
      gname <- c(gname, rep("", nrow(counts)))
    }
  }
  
  plot_data <- data.frame(counts = counts_reshape, label = label, class = class, gname = gname)
  plot_data$mechanism <- aggregate_mechanism(plot_data$class)
  class_order <- plot_data$label[order(plot_data$class)]
  gname_order <- plot_data$gname[order(plot_data$class)]
  
  p <- ggplot(data=plot_data, aes(x = label, y = counts, color = mechanism)) +
    geom_point() +
    theme_minimal() +
    scale_x_discrete(limits = class_order, labels = gname_order) +
    scale_color_manual(values=mechanism_pallete, name = "Resistance mechanism") +
    xlab(NULL) +
    ylab(expression("logCPM"["B"])) +
    ylim(c(-1,9)) +
    labs(title = title) +
    theme(axis.text.x = element_blank(),
          axis.text.y=element_text(size = rel(1.95)),
          axis.title=element_text(size=rel(1.95)),
          legend.title=element_text(size=rel(1.8)), 
          legend.text=element_text(size=rel(1.5)),
          plot.title = element_text(size = rel(4.5), face = "bold"),
          plot.title.position = "plot",
          plot.margin = margin(t = 10,
                               r = 10,
                               b = 10,
                               l = 10))
  
  return(p)
}

fig_s1a <- manhattan(hg_prom_hg, hg_prom_metadata, "a")
fig_s1b <- manhattan(ww_prom_hg, ww_prom_metadata, "b")
fig_s1c <- manhattan(hg_prom_ww, hg_prom_metadata, "c")
fig_s1d <- manhattan(ww_prom_ww, ww_prom_metadata, "d")


pdf("figure_s1.pdf", width = 35, height = 17)
fig_s1a + fig_s1b + fig_s1c + fig_s1d
dev.off()

#----------------- FIGURE S2a -------------------------------------------------

do_fisher <- function(gene_ids, missing, mechanism, metadata) {
  m <- matrix(c(0,0,0,0), nrow=2, ncol=2)
  m[1,1] <- sum(aggregate_mechanism(as.character(metadata["arg_class",gene_ids])) == mechanism) 
  m[1,2] <- sum(aggregate_mechanism(as.character(metadata["arg_class",!(colnames(metadata) %in% gene_ids) & !(colnames(metadata) %in% missing)])) == mechanism)
  m[2,1] <- sum(aggregate_mechanism(as.character(metadata["arg_class",gene_ids])) != mechanism) 
  m[2,2] <- sum(aggregate_mechanism(as.character(metadata["arg_class",!(colnames(metadata) %in% gene_ids) & !(colnames(metadata) %in% missing)])) != mechanism)
  
  f <- fisher.test(m)
  return(f)
}

co_prom_metadata <- metadata[,co_prom]
hg_prom_metadata <- metadata[,hg_prom]
ww_prom_metadata <- metadata[,ww_prom]
non_prom_metadata <- metadata[,non_prom]

plot_data <- data.frame(mechanism = rep(unique(aggregate_mechanism(gene_class)), 4),
                        category=c(rep("Co-promoted", 9), rep("HG-promoted", 9), 
                                   rep("WW-promoted", 9), rep("Non-promoted", 9)),
                        ratio = NA,
                        p_val = NA)

missing_est <- colnames(metadata[,missing] %>% select(which(.["category", ] == "Established")))
metadata_est <- metadata[,metadata["category", ] == "Established"]

for(mech in unique(plot_data$mechanism)) {
  
  plot_data$ratio[plot_data$mechanism == mech & plot_data$category == "Co-promoted"] <- as.numeric(do_fisher(colnames(co_prom_metadata[, co_prom_metadata["category",] == "Established"]), missing_est, mech, metadata_est)$estimate)
  plot_data$p_val[plot_data$mechanism == mech & plot_data$category == "Co-promoted"] <- as.numeric(do_fisher(colnames(co_prom_metadata[, co_prom_metadata["category",] == "Established"]), missing_est, mech, metadata_est)$p.value) 
  
  plot_data$ratio[plot_data$mechanism == mech & plot_data$category == "HG-promoted"] <- as.numeric(do_fisher(colnames(hg_prom_metadata[, hg_prom_metadata["category",] == "Established"]), missing_est, mech, metadata_est)$estimate)
  plot_data$p_val[plot_data$mechanism == mech & plot_data$category == "HG-promoted"] <- as.numeric(do_fisher(colnames(hg_prom_metadata[, hg_prom_metadata["category",] == "Established"]), missing_est, mech, metadata_est)$p.value)
  
  plot_data$ratio[plot_data$mechanism == mech & plot_data$category == "WW-promoted"] <- as.numeric(do_fisher(colnames(ww_prom_metadata[, ww_prom_metadata["category",] == "Established"]), missing_est, mech, metadata_est)$estimate)
  plot_data$p_val[plot_data$mechanism == mech & plot_data$category == "WW-promoted"] <- as.numeric(do_fisher(colnames(ww_prom_metadata[, ww_prom_metadata["category",] == "Established"]), missing_est, mech, metadata_est)$p.value)
  
  plot_data$ratio[plot_data$mechanism == mech & plot_data$category == "Non-promoted"] <- as.numeric(do_fisher(colnames(non_prom_metadata[, non_prom_metadata["category",] == "Established"]), missing_est, mech, metadata_est)$estimate)
  plot_data$p_val[plot_data$mechanism == mech & plot_data$category == "Non-promoted"] <- as.numeric(do_fisher(colnames(non_prom_metadata[, non_prom_metadata["category",] == "Established"]), missing_est, mech, metadata_est)$p.value)
  
}

plot_data$significance <- ""
plot_data$significance[plot_data$p_val < 0.01] <- "*"

fig_s2a <- ggplot(plot_data, aes(y=log(ratio), 
                               x=factor(category, levels = c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")),
                               fill = mechanism)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=significance), position=position_dodge(width=0.9), vjust=-0.25, size=7) +
  xlab(NULL) +
  ylab("log(Odds ratio)") +
  labs(tag = "a", title = "Established") +
  theme_minimal() +
  scale_fill_manual(values=mechanism_pallete, name = "Resistance mechanism") +
  theme(axis.text.x=element_text(size = rel(1.7)),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.5)), 
        legend.text=element_text(size=rel(1.3)),
        plot.tag = element_text(size=rel(3), face = "bold"),
        plot.title = element_text(size = rel(2)),
        plot.margin = margin(t = 10,
                             r = 10,
                             b = 10,
                             l = 20))

#----------------- FIGURE S2b -------------------------------------------------

plot_data <- data.frame(mechanism = rep(unique(aggregate_mechanism(gene_class)), 4),
                        category=c(rep("Co-promoted", 9), rep("HG-promoted", 9), 
                                   rep("WW-promoted", 9), rep("Non-promoted", 9)),
                        ratio = NA,
                        p_val = NA)

missing_lat <- colnames(metadata[,missing] %>% select(which(.["category", ] == "Latent")))
metadata_lat <- metadata[,metadata["category", ] == "Latent"]

for(mech in unique(plot_data$mechanism)) {
  
  plot_data$ratio[plot_data$mechanism == mech & plot_data$category == "Co-promoted"] <- as.numeric(do_fisher(colnames(co_prom_metadata[, co_prom_metadata["category",] == "Latent"]), missing_lat, mech, metadata_lat)$estimate)
  plot_data$p_val[plot_data$mechanism == mech & plot_data$category == "Co-promoted"] <- as.numeric(do_fisher(colnames(co_prom_metadata[, co_prom_metadata["category",] == "Latent"]), missing_lat, mech, metadata_lat)$p.value) 
  
  plot_data$ratio[plot_data$mechanism == mech & plot_data$category == "HG-promoted"] <- as.numeric(do_fisher(colnames(hg_prom_metadata[, hg_prom_metadata["category",] == "Latent"]), missing_lat, mech, metadata_lat)$estimate)
  plot_data$p_val[plot_data$mechanism == mech & plot_data$category == "HG-promoted"] <- as.numeric(do_fisher(colnames(hg_prom_metadata[, hg_prom_metadata["category",] == "Latent"]), missing_lat, mech, metadata_lat)$p.value)
  
  plot_data$ratio[plot_data$mechanism == mech & plot_data$category == "WW-promoted"] <- as.numeric(do_fisher(colnames(ww_prom_metadata[, ww_prom_metadata["category",] == "Latent"]), missing_lat, mech, metadata_lat)$estimate)
  plot_data$p_val[plot_data$mechanism == mech & plot_data$category == "WW-promoted"] <- as.numeric(do_fisher(colnames(ww_prom_metadata[, ww_prom_metadata["category",] == "Latent"]), missing_lat, mech, metadata_lat)$p.value)
  
  plot_data$ratio[plot_data$mechanism == mech & plot_data$category == "Non-promoted"] <- as.numeric(do_fisher(colnames(non_prom_metadata[, non_prom_metadata["category",] == "Latent"]), missing_lat, mech, metadata_lat)$estimate)
  plot_data$p_val[plot_data$mechanism == mech & plot_data$category == "Non-promoted"] <- as.numeric(do_fisher(colnames(non_prom_metadata[, non_prom_metadata["category",] == "Latent"]), missing_lat, mech, metadata_lat)$p.value)
  
}

plot_data$significance <- ""
plot_data$significance[plot_data$p_val < 0.01] <- "*"

fig_s2b <- ggplot(plot_data, aes(y=log(ratio), 
                               x=factor(category, levels = c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")),
                               fill = mechanism)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=significance), position=position_dodge(width=0.9), vjust=-0.25, size=7) +
  xlab(NULL) +
  ylab("log(Odds ratio)") +
  labs(tag = "b", title = "Latent") +
  theme_minimal() +
  scale_fill_manual(values=mechanism_pallete, name = "Resistance mechanism") +
  theme(axis.text.x=element_text(size = rel(1.7)),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.5)), 
        legend.text=element_text(size=rel(1.3)),
        plot.tag = element_text(size=rel(3), face = "bold"),
        plot.title = element_text(size = rel(2)),
        plot.margin = margin(t = 10,
                             r = 10,
                             b = 10,
                             l = 20))


pdf("figure_s2.pdf", height = 15, width = 17)
fig_s2a /
  fig_s2b 
dev.off()

#----------------- FIGURE S3a-b -----------------------------------------------

inter_phyla_histogram <- function(gene_ids, metadata, title, tag) {
  plot_data <- data.frame(n = 2:9,
                          value = NA)
  
  for (i in 1:nrow(plot_data)) {
    plot_data$value[i] <- sum(as.numeric(metadata["n_phyla", gene_ids]) == plot_data$n[i])/length(gene_ids)
  }
  
  p <- ggplot(data = plot_data, aes(x = n, y = value)) +
    geom_bar(position="dodge", stat="identity") +
    theme_minimal() +
    ylab("Proportion")  +
    xlab("Number of host phyla") +
    scale_y_continuous(labels = percent, limits = c(0,0.8)) +
    scale_x_continuous(breaks = 1:11) +
    labs(tag = tag, title = title) +
    theme(axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
          axis.text.y=element_text(size = rel(1.3)),
          axis.title=element_text(size=rel(1.3)),
          legend.title=element_text(size=rel(1.2)), 
          legend.text=element_text(size=rel(1)),
          plot.tag = element_text(size=rel(2.5), face = "bold"),
          plot.title = element_text(size = rel(1.5)))
}

inter_phylum_cop <- colnames(metadata[,colnames(metadata) %in% co_prom & metadata["taxonomy_score",] == "Phylum"])
inter_phylum_hgp <- colnames(metadata[,colnames(metadata) %in% hg_prom & metadata["taxonomy_score",] == "Phylum"])
inter_phylum_wwp <- colnames(metadata[,colnames(metadata) %in% ww_prom & metadata["taxonomy_score",] == "Phylum"])
inter_phylum_nop <- colnames(metadata[,colnames(metadata) %in% non_prom & metadata["taxonomy_score",] == "Phylum"])

fig_s3a <- inter_phyla_histogram(inter_phylum_cop, metadata, "Co-promoted (n=20)", "a")
fig_s3b <- inter_phyla_histogram(inter_phylum_hgp, metadata, "HG-promoted (n=17)", "b")
fig_s3c <- inter_phyla_histogram(inter_phylum_wwp, metadata, "WW-promoted (n=4)", "c")
fig_s3d <- inter_phyla_histogram(inter_phylum_nop, metadata, "Non-promoted (n=28)", "d")


pdf("figure_s3.pdf", height = 10, width = 12)
fig_s3a + fig_s3b + fig_s3c + fig_s3d
dev.off()

#----------------- FIGURE S4 prep ---------------------------------------------

read_mges <- function(filename, gene_ids, metadata, mge_type, arg_type, status) {
  df <- data.frame(fread(filename)) %>% 
    tibble::column_to_rownames('V1') %>%
    select(subset = -ncol(.)) 
  
  metadata_subset <- metadata[,gene_ids]
  df_subset <- df[,metadata_subset["category",] == status]
  
  output <- data.frame(prop = rowSums(df_subset[mge_type$subset != "arg",]) / ncol(df_subset), type = mge_type$subset[mge_type$subset != "arg"])
  
  args <- df_subset[mge_type$subset == "arg",]
  arg_summary <- data.frame(prop = rep(NA, 6), type = rep("arg", 6), row.names=unique(arg_type$V2))
  
  for (i in 1:nrow(arg_summary)) {
    arg_summary$prop[i] <- sum(colSums(df_subset[arg_type$V1[arg_type$V2 == rownames(arg_summary)[i]],]) > 0) / ncol(df_subset)
  }
  
  output <- rbind(output, arg_summary) %>% 
    filter(.$type != "is")
  
  output$status <- status
  output$n <- ncol(df_subset)
  
  return(output)
}

do_fisher_test <- function(df) {
  df$category <- gsub("\\\n.*", "", df$category)
  
  contingency <- df %>%
    select(element, category, value) %>%
    pivot_wider(names_from = category, values_from = value, values_fill = 0) %>%
    column_to_rownames("element")
  
  categories <- colnames(contingency)
  
  # Safe pairwise Fisher test with proper error handling
  pairwise_results <- combn(categories, 2, function(pair) {
    sub_table <- contingency[, pair, drop = FALSE]
    sub_table <- sub_table[rowSums(sub_table) > 0, , drop = FALSE]
    
    # Must be at least 2 rows and 2 columns
    if (nrow(sub_table) < 2 || ncol(sub_table) != 2) {
      return(data.frame(cat1 = pair[1], cat2 = pair[2], p_value = NA, method = "skipped"))
    }
    
    mat <- as.matrix(sub_table)
    if (any(is.na(mat)) || any(mat < 0)) {
      return(data.frame(cat1 = pair[1], cat2 = pair[2], p_value = NA, method = "invalid"))
    }
    
    # Try Fisher test
    test <- tryCatch({
      fisher.test(mat)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(test)) {
      data.frame(cat1 = pair[1], cat2 = pair[2], p_value = test$p.value, method = test$method)
    } else {
      data.frame(cat1 = pair[1], cat2 = pair[2], p_value = NA, method = "error")
    }
  }, simplify = FALSE)
  
  results_df <- do.call(rbind, pairwise_results)
  print(results_df)
}

mge_type <- data.frame(fread("co_prom_mges.tsv")) %>% 
  tibble::column_to_rownames('V1') %>%
  select(subset = ncol(.))

arg_type <- data.frame(fread("arg_types.tsv", header=FALSE))

est_co_prom_mge <- read_mges("co_prom_mges.tsv", co_prom, metadata, mge_type, arg_type, "Established")
lat_co_prom_mge <- read_mges("co_prom_mges.tsv", co_prom, metadata, mge_type, arg_type, "Latent")

est_hg_prom_mge <- read_mges("hg_prom_mges.tsv", hg_prom, metadata, mge_type, arg_type, "Established")
lat_hg_prom_mge <- read_mges("hg_prom_mges.tsv", hg_prom, metadata, mge_type, arg_type, "Latent")

est_ww_prom_mge <- read_mges("ww_prom_mges.tsv", ww_prom, metadata, mge_type, arg_type, "Established")
lat_ww_prom_mge <- read_mges("ww_prom_mges.tsv", ww_prom, metadata, mge_type, arg_type, "Latent")

est_non_prom_mge <- read_mges("non_prom_mges.tsv", non_prom, metadata, mge_type, arg_type, "Established")
lat_non_prom_mge <- read_mges("non_prom_mges.tsv", non_prom, metadata, mge_type, arg_type, "Latent")

mge_type <- mge_type[mge_type$subset != "is",]

plot_data <- data.frame(element=rep(rownames(est_co_prom_mge),8),
                        type=rep(est_co_prom_mge$type,8),
                        value=c(est_co_prom_mge$prop, 
                                est_hg_prom_mge$prop, 
                                est_ww_prom_mge$prop, 
                                est_non_prom_mge$prop,
                                lat_co_prom_mge$prop, 
                                lat_hg_prom_mge$prop, 
                                lat_ww_prom_mge$prop, 
                                lat_non_prom_mge$prop),
                        category=c(rep(c(rep("Co-promoted", length(est_co_prom_mge$prop)),
                                         rep("HG-promoted", length(est_co_prom_mge$prop)), 
                                         rep("WW-promoted", length(est_co_prom_mge$prop)), 
                                         rep("Non-promoted", length(est_co_prom_mge$prop))),2)),
                        status=c(est_co_prom_mge$status, 
                                 est_hg_prom_mge$status, 
                                 est_ww_prom_mge$status, 
                                 est_non_prom_mge$status,
                                 lat_co_prom_mge$status, 
                                 lat_hg_prom_mge$status, 
                                 lat_ww_prom_mge$status, 
                                 lat_non_prom_mge$status),
                        n=c(est_co_prom_mge$n, 
                            est_hg_prom_mge$n, 
                            est_ww_prom_mge$n, 
                            est_non_prom_mge$n,
                            lat_co_prom_mge$n, 
                            lat_hg_prom_mge$n, 
                            lat_ww_prom_mge$n, 
                            lat_non_prom_mge$n))

#----------------- FIGURE S4a -------------------------------------------------

plot_data_subset <- plot_data[plot_data$type == "mpf" & plot_data$status == "Established",]
plot_data_subset$category <- paste0(plot_data_subset$category, "\n(n=", plot_data_subset$n, ")")

fig_a <- ggplot(
  plot_data_subset, 
  aes(
    y=value, 
    x=factor(category, levels=c("Non-promoted\n(n=229)", "WW-promoted\n(n=20)", "HG-promoted\n(n=42)", "Co-promoted\n(n=41)")), 
    fill=factor(element, levels = rev(unique(element))))
  ) + 
  geom_bar(position="dodge", stat="identity")  +
  theme_minimal() +
  ylab("Proportion")  +
  xlab(NULL) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(tag = "a", title = "Established") +
  scale_fill_manual(
    values=c(
      "B"='#337538',
      "F" ="#dccd7d",
      "FA" = "#2e2585",
      "FATA"="#94caec",
      "G" = "#dddddd",
      "I"="#c26a77",
      "T" = "#9f4a96"
      ), 
    name = "MPF type",
    breaks = unique(plot_data_subset$element)
  ) +
  theme(
    axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
    axis.text.y=element_text(size = rel(1.3)),
    axis.title=element_text(size=rel(1.3)),
    legend.title=element_text(size=rel(1.2)), 
    legend.text=element_text(size=rel(1)),
    plot.tag = element_text(size=rel(2.5), face = "bold"),
    plot.title = element_text(size = rel(1.5))
  ) +
  coord_flip()

print("Established MPF")
plot_data_subset$value <- plot_data_subset$value*plot_data_subset$n
do_fisher_test(plot_data_subset)

#----------------- FIGURE S4b -------------------------------------------------

plot_data_subset <- plot_data[plot_data$type == "mpf" & plot_data$status == "Latent",]
plot_data_subset$category <- paste0(plot_data_subset$category, "\n(n=", plot_data_subset$n, ")")

fig_b <- ggplot(
  plot_data_subset, 
  aes(
    y=value, 
    x=factor(category, levels=c("Non-promoted\n(n=1159)", "WW-promoted\n(n=28)", "HG-promoted\n(n=70)", "Co-promoted\n(n=3)")), 
    fill=factor(element, levels = rev(unique(element))))
  ) + 
  geom_bar(position="dodge", stat="identity")  +
  theme_minimal() +
  ylab("Proportion")  +
  xlab(NULL) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(tag = "b", title = "Latent") +
  scale_fill_manual(
    values=c(
      "B"='#337538',
      "F" ="#dccd7d",
      "FA" = "#2e2585",
      "FATA"="#94caec",
      "G" = "#dddddd",
      "I"="#c26a77",
      "T" = "#9f4a96"
      ), 
    name = "MPF type",
    breaks = unique(plot_data_subset$element)
  ) +
  theme(
    axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
    axis.text.y=element_text(size = rel(1.3)),
    axis.title=element_text(size=rel(1.3)),
    legend.title=element_text(size=rel(1.2)), 
    legend.text=element_text(size=rel(1)),
    plot.tag = element_text(size=rel(2.5), face = "bold"),
    plot.title = element_text(size = rel(1.5))
  ) +
  coord_flip()

print("Latent MPF")
plot_data_subset$value <- plot_data_subset$value*plot_data_subset$n
do_fisher_test(plot_data_subset)

#----------------- FIGURE S4c -------------------------------------------------

plot_data_subset <- plot_data[plot_data$type == "mob" & plot_data$status == "Established",]
plot_data_subset$category <- paste0(plot_data_subset$category, "\n(n=", plot_data_subset$n, ")")

fig_c <- ggplot(
  plot_data_subset, 
  aes(
    y=value, 
    x=factor(category, levels=c("Non-promoted\n(n=229)", "WW-promoted\n(n=20)", "HG-promoted\n(n=42)", "Co-promoted\n(n=41)")), 
    fill=factor(element, levels = rev(unique(element))))
  ) + 
  geom_bar(position="dodge", stat="identity")  +
  theme_minimal() +
  ylab("Proportion") +
  xlab(NULL) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(tag = "c", title = "Established") +
  scale_fill_manual(
    values=c(
      "MOBC"='#337538',
      "MOBF" ="#dccd7d",
      "MOBH" = "#2e2585",
      "MOBP"="#94caec",
      "MOBQ" = "#dddddd",
      "MOBT"="#c26a77",
      "MOBV" = "#9f4a96"
      ), 
    name = "Relaxase type",
    breaks = unique(plot_data_subset$element)
  ) +
  theme(
    axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
    axis.text.y=element_text(size = rel(1.3)),
    axis.title=element_text(size=rel(1.3)),
    legend.title=element_text(size=rel(1.2)), 
    legend.text=element_text(size=rel(1)),
    plot.tag = element_text(size=rel(2.5), face = "bold"),
    plot.title = element_text(size = rel(1.5))
  ) +
  coord_flip()

print("Established MOB")
plot_data_subset$value <- plot_data_subset$value*plot_data_subset$n
do_fisher_test(plot_data_subset)

#----------------- FIGURE S4d -------------------------------------------------

plot_data_subset <- plot_data[plot_data$type == "mob" & plot_data$status == "Latent",]
plot_data_subset$category <- paste0(plot_data_subset$category, "\n(n=", plot_data_subset$n, ")")

fig_d <- ggplot(
  plot_data_subset, 
  aes(
    y=value, 
    x=factor(category, levels=c("Non-promoted\n(n=1159)", "WW-promoted\n(n=28)", "HG-promoted\n(n=70)", "Co-promoted\n(n=3)")), 
    fill=factor(element, levels = rev(unique(element))))
  ) + 
  geom_bar(position="dodge", stat="identity")  +
  theme_minimal() +
  ylab("Proportion") +
  xlab(NULL) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(tag = "d", title = "Latent") +
  scale_fill_manual(
    values=c(
      "MOBC"='#337538',
      "MOBF" ="#dccd7d",
      "MOBH" = "#2e2585",
      "MOBP"="#94caec",
      "MOBQ" = "#dddddd",
      "MOBT"="#c26a77",
      "MOBV" = "#9f4a96"
    ), 
    name = "Relaxase type",
    breaks = unique(plot_data_subset$element)
  ) +
  theme(
    axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
    axis.text.y=element_text(size = rel(1.3)),
    axis.title=element_text(size=rel(1.3)),
    legend.title=element_text(size=rel(1.2)), 
    legend.text=element_text(size=rel(1)),
    plot.tag = element_text(size=rel(2.5), face = "bold"),
    plot.title = element_text(size = rel(1.5))
  ) +
  coord_flip()

print("Latent MOB")
plot_data_subset$value <- plot_data_subset$value*plot_data_subset$n
do_fisher_test(plot_data_subset)

#----------------- FIGURE S4e -------------------------------------------------

plot_data_subset <- plot_data[plot_data$type == "arg" & plot_data$status == "Established",]
plot_data_subset$category <- paste0(plot_data_subset$category, "\n(n=", plot_data_subset$n, ")")
plot_data_subset$element <- factor(
  plot_data_subset$element,
  levels = rev(c("Aminoglycoside", "Beta-lactam", "Macrolide", "Tetracycline", "Quinolone", "Other"))
)

fig_e <- ggplot(
  plot_data_subset, 
  aes(
    y=value, 
    x=factor(category, levels=c("Non-promoted\n(n=229)", "WW-promoted\n(n=20)", "HG-promoted\n(n=42)", "Co-promoted\n(n=41)")),
    fill=factor(element, levels=rev(c("Aminoglycoside", "Beta-lactam", "Macrolide", "Tetracycline", "Quinolone", "Other"))))
  ) + 
  geom_bar(position="dodge", stat="identity")  +
  theme_minimal() +
  ylab("Proportion")  +
  xlab(NULL) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(tag = "e", title = "Established") +
  scale_fill_manual(
    values=c(
      "Aminoglycoside"='#337538',
      "Beta-lactam" ="#dccd7d",
      "Macrolide" = "#2e2585",
      "Tetracycline"="#94caec",
      "Quinolone" = "#dddddd",
      "Other" = "#c26a77"
    ), 
    name = "Antibiotic class",
    breaks = c("Aminoglycoside", "Beta-lactam", "Macrolide", "Tetracycline", "Quinolone", "Other")
  ) +
  theme(
    axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
    axis.text.y=element_text(size = rel(1.3)),
    axis.title=element_text(size=rel(1.3)),
    legend.title=element_text(size=rel(1.2)), 
    legend.text=element_text(size=rel(1)),
    plot.tag = element_text(size=rel(2.5), face = "bold"),
    plot.title = element_text(size = rel(1.5))
  ) +
  coord_flip()

print("Established co-localized ARGs")
plot_data_subset$value <- plot_data_subset$value*plot_data_subset$n
do_fisher_test(plot_data_subset)

#----------------- FIGURE S4f -------------------------------------------------

plot_data_subset <- plot_data[plot_data$type == "arg" & plot_data$status == "Latent",]
plot_data_subset$category <- paste0(plot_data_subset$category, "\n(n=", plot_data_subset$n, ")")
plot_data_subset$element <- factor(
  plot_data_subset$element,
  levels = rev(c("Aminoglycoside", "Beta-lactam", "Macrolide", "Tetracycline", "Quinolone", "Other"))
)

fig_f <- ggplot(
  plot_data_subset, 
  aes(
    y=value, 
    x=factor(category, levels=c("Non-promoted\n(n=1159)", "WW-promoted\n(n=28)", "HG-promoted\n(n=70)", "Co-promoted\n(n=3)")),
    fill=element)
  ) + 
  geom_bar(position="dodge", stat="identity")  +
  theme_minimal() +
  ylab("Proportion")  +
  xlab(NULL) +
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  labs(tag = "f", title = "Latent") +
  scale_fill_manual(
    values=c(
      "Aminoglycoside"='#337538',
      "Beta-lactam" ="#dccd7d",
      "Macrolide" = "#2e2585",
      "Tetracycline"="#94caec",
      "Quinolone" = "#dddddd",
      "Other" = "#c26a77"
      ), 
    name = "Antibiotic class",
    breaks = c("Aminoglycoside", "Beta-lactam", "Macrolide", "Tetracycline", "Quinolone", "Other")
  ) +
  theme(
    axis.text.x=element_text(size = rel(1.3), angle = 0, vjust = 1, hjust = 1),
    axis.text.y=element_text(size = rel(1.3)),
    axis.title=element_text(size=rel(1.3)),
    legend.title=element_text(size=rel(1.2)), 
    legend.text=element_text(size=rel(1)),
    plot.tag = element_text(size=rel(2.5), face = "bold"),
    plot.title = element_text(size = rel(1.5))
  ) +
  coord_flip()

print("Latent co-localized ARGs")
plot_data_subset$value <- plot_data_subset$value*plot_data_subset$n
do_fisher_test(plot_data_subset)

pdf("figure_s4.pdf", height = 14, width = 14)
(fig_a | fig_b) /
  (fig_c | fig_d) /
  (fig_e | fig_f)
dev.off()

#----------------- FIGURE S5 --------------------------------------------------

perform_permutation <- function(df, n_iter, n_subsample) {
  # Pre-allocate list to collect results
  results_list <- vector("list", length = n_iter * length(unique(df$n_pathogens)))
  df$n_pathogens <- as.numeric(df$n_pathogens)
  
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  n_pathogens <- sort(unique(df$n_pathogens))
  group_levels <- c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")
  
  result_idx <- 1
  
  for (i in seq_len(n_iter)) {
    df$grp <- sample(df$grp)  # scramble groups once per iteration
    
    for (n in n_pathogens) {
      tmp_mat <- matrix(0, nrow = n_iter, ncol = 4)  # preallocate temporary results
      
      for (j in seq_len(n_iter)) {
        # Sample once per group
        subsample_indices <- unlist(
          lapply(group_levels, function(g) {
            sample(which(df$grp == g), n_subsample, replace = FALSE)
          })
        )
        subsample <- df[subsample_indices, ]
        
        n_args <- sum(subsample$n_pathogens >= n)
        
        if (n_args == 0) {
          next
        }
        
        # Calculate fraction of each group matching score
        group_counts <- sapply(group_levels, function(g) {
          sum(subsample$n_pathogens[subsample$grp == g] >= n)
        })
        
        tmp_mat[j, ] <- group_counts / n_args
      }
      
      # Take column means and store result
      tmp_means <- colMeans(tmp_mat, na.rm = TRUE)
      tmp_comb <- data.frame(mean = tmp_means,category = names(group_counts), n_pathogens = n)
      
      results_list[[result_idx]] <- tmp_comb
      result_idx <- result_idx + 1
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Combine list into data.frame
  output <- do.call(rbind, results_list)
  output <- as.data.frame(output)
  
  return(output)
}

plot_n_pathogens <- function(metadata, status, n_iter) {
  co_prom_metadata <- metadata[,co_prom] %>% select(which(.["category", ] == status))
  hg_prom_metadata <- metadata[,hg_prom] %>% select(which(.["category", ] == status))
  ww_prom_metadata <- metadata[,ww_prom] %>% select(which(.["category", ] == status))
  non_prom_metadata <- metadata[,non_prom] %>% select(which(.["category", ] == status))
  
  plot_data <- data.frame(n_pathogens = rep(0:max(as.numeric(metadata["n_pathogens", metadata["category",] == status])), 4),
                          value = NA)
  
  n_subsample <- min(ncol(co_prom_metadata), ncol(hg_prom_metadata), ncol(ww_prom_metadata), ncol(non_prom_metadata))
  
  
  for (n in unique(plot_data$n_pathogens)) {
    prop_co_prom <- c()
    
    for (i in 1:n_iter) {
      subsample_co_prom <- co_prom_metadata[,sample(colnames(co_prom_metadata), n_subsample, replace=FALSE)]
      subsample_hg_prom <- hg_prom_metadata[,sample(colnames(hg_prom_metadata), n_subsample, replace=FALSE)]
      subsample_ww_prom <- ww_prom_metadata[,sample(colnames(ww_prom_metadata), n_subsample, replace=FALSE)]
      subsample_non_prom <- non_prom_metadata[,sample(colnames(non_prom_metadata), n_subsample, replace=FALSE)]
      
      if (status == "Established") {
        subsample_comb <- cbind(subsample_co_prom, subsample_hg_prom, subsample_ww_prom, subsample_non_prom)
      }
      else if (status == "Latent") {
        subsample_comb <- cbind(subsample_hg_prom, subsample_ww_prom, subsample_non_prom)
      }
      
      n_args <- sum(as.numeric(subsample_comb["n_pathogens",]) >= n)
      
      if (n_args == 0) {
        next
      }
      
      else {
        prop_co_prom <- c(prop_co_prom, sum(as.numeric(subsample_co_prom["n_pathogens",]) >= n)/n_args)
      }
    }
    
    plot_data$value[plot_data$n_pathogens == n] <- mean(prop_co_prom)
  }
  
  print("Performing permutation analysis")
  permutation_df <- cbind(co_prom_metadata, 
                          hg_prom_metadata, 
                          ww_prom_metadata, 
                          non_prom_metadata) %>% 
    t() %>% 
    as.data.frame()
  
  permutation_df$grp <- c(rep("Co-promoted", ncol(co_prom_metadata)), 
                          rep("HG-promoted", ncol(hg_prom_metadata)), 
                          rep("WW-promoted", ncol(ww_prom_metadata)), 
                          rep("Non-promoted", ncol(non_prom_metadata)))
  
  permutation_results <- perform_permutation(permutation_df, n_iter, n_subsample)
  permutation_results <- permutation_results[permutation_results$category == "Co-promoted",]
  
  plot_data_permutation <- data.frame(mean = NA, sd = NA, n_pathogens = unique(permutation_results$n_pathogens))
  
  for (n in permutation_results$n_pathogens) {
    plot_data_permutation$mean[plot_data_permutation$n_pathogens == n] <- mean(permutation_results$mean[permutation_results$n_pathogens == n])
    plot_data_permutation$sd[plot_data_permutation$n_pathogens == n] <- sd(permutation_results$mean[permutation_results$n_pathogens == n])
  }
  
  p <- ggplot() +
    geom_line(data=plot_data, 
              aes(x = n_pathogens, y = value), 
              linewidth=1, 
              color = "#000000") +
    geom_ribbon(data = plot_data_permutation, 
                aes(x = n_pathogens,
                    ymin = mean - sd,
                    ymax = mean + sd),
                fill = "#636363", alpha = 0.2) +
    geom_line(data = plot_data_permutation,
              aes(x = n_pathogens, y = mean),
              color = "#636363", linewidth = 1, linetype = "dashed") +
    theme_minimal() +
    scale_x_continuous(breaks = 0:15) +
    scale_y_continuous(labels = percent) +
    ylab("Proportion") +
    xlab(expression("Number of Pathogens (" >= ")")) +
    theme(axis.text.x=element_text(size = rel(1.5)),
          axis.text.y=element_text(size = rel(1.3)),
          axis.title=element_text(size=rel(1.3)),
          legend.title=element_text(size=rel(1.2)), 
          legend.text=element_text(size=rel(1)),
          plot.tag = element_text(size=rel(2.5), face = "bold"),
          plot.title = element_text(size = rel(1.5)))
  
  return(p)
}

fig_s5 <- plot_n_pathogens(metadata, "Established", 1000)

pdf("figure_s5.pdf", height = 6, width = 8)
plot(fig_s5)
dev.off()

#----------------- FIGURE S5a -------------------------------------------------

library(tidyr)
read_genetic_compatibility <- function(filename, status, category) {
  df <- data.frame(fread(filename)) %>% 
    tibble::column_to_rownames('X0') 
  
  status_ids <- metadata[, colnames(df)] %>% 
    select(which(.["category", ] == status)) %>%
    colnames()
  
  df <- df[,status_ids]
  
  value <- c()
  species <- c()
  for (sp in rownames(df)) {
    value <- c(value, as.numeric(df[sp,]))
    species <- c(species, rep(sp, ncol(df)))
  }
  
  output <- data.frame(Value = value, 
                       Species = species, 
                       Category = paste(c(category, " (n=", ncol(df), ")"), collapse=""))
  
  return(output)
}

gcomp_co_prom <- read_genetic_compatibility("gcomp_co_prom.tsv", "Established", "Co-promoted")
gcomp_hg_prom <- read_genetic_compatibility("gcomp_hg_prom.tsv", "Established", "HG-promoted")
gcomp_ww_prom <- read_genetic_compatibility("gcomp_ww_prom.tsv", "Established", "WW-promoted")
gcomp_non_prom <- read_genetic_compatibility("gcomp_non_prom.tsv", "Established", "Non-promoted")

plot_data <- rbind(gcomp_co_prom,
                   gcomp_hg_prom,
                   gcomp_ww_prom,
                   gcomp_non_prom)

fig_s6a <- ggplot(plot_data, aes(x = Species, y = Value, fill = factor(Category, levels = c("Co-promoted (n=41)", "HG-promoted (n=42)", "WW-promoted (n=20)", "Non-promoted (n=229)")))) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = NULL,
       y = "Gene-genome 5mer distance",
       tag = "a",
       title = "Established",
       fill = NULL) +
  scale_fill_manual(values=c(
    "Co-promoted (n=41)"='#000000',
    "HG-promoted (n=42)" ="#D55E00",
    "WW-promoted (n=20)" = "#0072B2",
    "Non-promoted (n=229)"="#F0E442"
  ), name = "") +
  theme(axis.text.x=element_text(size = rel(1.3), angle = 45, vjust = 1, hjust = 1, face = "italic"),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1)),
        plot.tag = element_text(size=rel(2.5), face = "bold"),
        plot.title = element_text(size = rel(2)),
        plot.margin = ggplot2::margin(t = 10,
                                      r = 10,
                                      b = 10,
                                      l = 60))

# Wilcoxon's signed rank test

combinations <- data.frame(crossing(x = unique(plot_data$Category), 
                                    y  = unique(plot_data$Category))) %>%
  filter(x != y)

test_data <- as.data.frame(matrix(nrow=nrow(combinations), ncol=15))
colnames(test_data) <- unique(plot_data$Species)
test_data <- cbind(combinations, test_data)

for (i in 1:nrow(test_data)) {
  for (j in 3:ncol(test_data)) {
    dist1 <- plot_data$Value[plot_data$Category == test_data$x[i] & plot_data$Species == colnames(test_data)[j]]
    dist2 <- plot_data$Value[plot_data$Category == test_data$y[i] & plot_data$Species == colnames(test_data)[j]]
    test_data[i, j] <- median(dist1) - median(dist2)
  }
}                         

pval <- c()
for (i in 1:nrow(test_data)) {
  pval <- c(pval, wilcox.test(x = as.numeric(test_data[i, 3:ncol(test_data)]), alternative = "less")$p.value)
}

test_results <- data.frame(x = test_data$x, 
                           y = test_data$y, 
                           p_value = pval)

#----------------- FIGURE S5b -------------------------------------------------

gcomp_hg_prom <- read_genetic_compatibility("gcomp_hg_prom.tsv", "Latent", "HG-promoted")
gcomp_ww_prom <- read_genetic_compatibility("gcomp_ww_prom.tsv", "Latent", "WW-promoted")
gcomp_non_prom <- read_genetic_compatibility("gcomp_non_prom.tsv", "Latent", "Non-promoted")

plot_data <- rbind(gcomp_hg_prom,
                   gcomp_ww_prom,
                   gcomp_non_prom)

fig_s6b <- ggplot(plot_data, aes(x = Species, y = Value, fill = factor(Category, levels = c("HG-promoted (n=70)", "WW-promoted (n=28)", "Non-promoted (n=1159)")))) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = NULL,
       y = "Gene-genome 5mer distance",
       tag = "b",
       title = "Latent",
       fill = NULL) +
  scale_fill_manual(values=c(
    "HG-promoted (n=70)" ="#D55E00",
    "WW-promoted (n=28)" = "#0072B2",
    "Non-promoted (n=1159)"="#F0E442"
  ), name = "") +
  theme(axis.text.x=element_text(size = rel(1.3), angle = 45, vjust = 1, hjust = 1, face = "italic"),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1)),
        plot.tag = element_text(size=rel(3), face = "bold"),
        plot.title = element_text(size = rel(2)),
        plot.margin = ggplot2::margin(t = 10,
                                      r = 10,
                                      b = 10,
                                      l = 60))

# Wilcoxon's signed rank test

combinations <- data.frame(crossing(x = unique(plot_data$Category), 
                                    y  = unique(plot_data$Category))) %>%
  filter(x != y)

test_data <- as.data.frame(matrix(nrow=nrow(combinations), ncol=15))
colnames(test_data) <- unique(plot_data$Species)
test_data <- cbind(combinations, test_data)

for (i in 1:nrow(test_data)) {
  for (j in 3:ncol(test_data)) {
    dist1 <- plot_data$Value[plot_data$Category == test_data$x[i] & plot_data$Species == colnames(test_data)[j]]
    dist2 <- plot_data$Value[plot_data$Category == test_data$y[i] & plot_data$Species == colnames(test_data)[j]]
    test_data[i, j] <- median(dist1) - median(dist2)
  }
}                         

pval <- c()
for (i in 1:nrow(test_data)) {
  pval <- c(pval, wilcox.test(x = as.numeric(test_data[i, 3:ncol(test_data)]), alternative = "less")$p.value)
}

test_results <- data.frame(x = test_data$x, 
                           y = test_data$y, 
                           p_value = pval)


pdf("figure_s6.pdf", height = 14, width = 12)
fig_s6a /
  fig_s6b
dev.off()

#----------------- FIGURE S6 --------------------------------------------------

read_genetic_compatibility <- function(filename, status, category) {
  df <- data.frame(fread(filename)) %>% 
    tibble::column_to_rownames('X0') 
  
  status_ids <- metadata[, colnames(df)] %>% 
    select(which(.["category", ] == status)) %>%
    colnames()
  
  df <- df[,status_ids]
  
  value <- c()
  species <- c()
  for (sp in rownames(df)) {
    value <- c(value, as.numeric(df[sp,]))
    species <- c(species, rep(sp, ncol(df)))
  }
  
  output <- data.frame(Value = value, 
                       Species = species, 
                       Category = paste(c(category, " (n=", ncol(df), ")"), collapse=""))
  
  output$Species <- sub(" .*", "", output$Species)
  
  return(output)
}

plot_gcomp <- function(categ, title, tag) {
  gcomp_co_prom <- read_genetic_compatibility(paste(c("gcomp_co_promoted_", categ, ".txt"), collapse = ""), "Established", "Co-promoted")
  gcomp_hg_prom <- read_genetic_compatibility(paste(c("gcomp_hg_promoted_", categ, ".txt"), collapse = ""), "Established", "HG-promoted")
  gcomp_ww_prom <- read_genetic_compatibility(paste(c("gcomp_ww_promoted_", categ, ".txt"), collapse = ""), "Established", "WW-promoted")
  gcomp_non_prom <- read_genetic_compatibility(paste(c("gcomp_non_promoted_", categ, ".txt"), collapse = ""), "Established", "Non-promoted")
  
  plot_data <- rbind(gcomp_co_prom, gcomp_hg_prom, gcomp_ww_prom, gcomp_non_prom)
  
  p <- ggplot(plot_data, aes(x = Species, y = Value, fill = factor(Category, levels = c("Co-promoted (n=41)", "HG-promoted (n=42)", "WW-promoted (n=20)", "Non-promoted (n=229)")))) +
    geom_boxplot() +
    theme_minimal() +
    ylim(c(0.02, 0.09)) +
    labs(x = NULL,
         y = "Gene-genome 5mer distance",
         tag = tag,
         title = title,
         fill = NULL) +
    scale_fill_manual(values=c(
      "Co-promoted (n=41)"='#000000',
      "HG-promoted (n=42)" ="#D55E00",
      "WW-promoted (n=20)" = "#0072B2",
      "Non-promoted (n=229)"="#F0E442"
    ), name = "") +
    theme(axis.text.x=element_text(size = rel(1.3), angle = 45, vjust = 1, hjust = 1, face = "italic"),
          axis.text.y=element_text(size = rel(1.3)),
          axis.title=element_text(size=rel(1.3)),
          legend.title=element_text(size=rel(1.2)), 
          legend.text=element_text(size=rel(1)),
          plot.tag = element_text(size=rel(2.5), face = "bold"),
          plot.title = element_text(size = rel(2)),
          plot.margin = ggplot2::margin(t = 10,
                                        r = 10,
                                        b = 10,
                                        l = 60))
  
  return(p)
}


fig_s7a <- plot_gcomp("a-f", "Abundant bacterial genera: A-F", "a")
fig_s7b <- plot_gcomp("g-p", "Abundant bacterial genera: G-P", "b")
fig_s7c <- plot_gcomp("r-z", "Abundant bacterial genera: R-Z", "c")

pdf("figure_s7.pdf", width=20, height = 20)
fig_s7a /
  fig_s7b /
  fig_s7c
dev.off()

#----------------- FIGURE S8 --------------------------------------------------

read_genetic_compatibility <- function(filename, rel_ab, status, category) {
  df <- data.frame(fread(filename)) %>% 
    tibble::column_to_rownames('X0') 
  
  status_ids <- metadata[, colnames(df)] %>% 
    select(which(.["category", ] == status)) %>%
    colnames()
  
  df <- df[,status_ids]
  prevalence <- rel_ab[status_ids]
  
  genus <- unique(sub(" .*", "", rownames(df)))
  
  tmp <- as.data.frame(matrix(nrow = length(genus), ncol = ncol(df)))
  colnames(tmp) <- colnames(df)
  rownames(tmp) <- genus
  
  for (g in genus) {
    tmp[g,] <- apply(df[sub(" .*", "", rownames(df)) == g,], 2, median)
    
  }
  
  output <- data.frame(comp = apply(tmp, 2, mean), 
                       prevalence = prevalence,
                       category = category)
  
  return(output)
}

counts <- rbind(data_hm, data_ww)
rel_ab <- colSums(counts >= 3)/nrow(counts)

gcomp_co_prom <- read_genetic_compatibility("genus_gcomp_co_promoted.tsv", rel_ab, "Established", "Co-promoted")
gcomp_hg_prom <- read_genetic_compatibility("genus_gcomp_hg_promoted.tsv", rel_ab, "Established", "HG-promoted")
gcomp_ww_prom <- read_genetic_compatibility("genus_gcomp_ww_promoted.tsv", rel_ab, "Established", "WW-promoted")
gcomp_non_prom <- read_genetic_compatibility("genus_gcomp_non_promoted.tsv", rel_ab, "Established", "Non-promoted")

cor_table <- rbind(gcomp_co_prom, gcomp_hg_prom, gcomp_ww_prom, gcomp_non_prom)
cor.test(cor_table$comp, cor_table$prevalence, method = "spearman")

fig_s8 <- ggplot(data=cor_table, aes(x=comp, y = prevalence, color=factor(category, levels = c("Co-promoted", "HG-promoted", "WW-promoted", "Non-promoted")))) +
  annotate("segment", x=min(cor_table$comp)*0.95, xend=max(cor_table$comp)*1.05, y=min(cor_table$prevalence), yend=min(cor_table$prevalence),
           arrow=arrow(length=unit(0.3, "cm"), type="closed")) + 
  annotate("segment", x=min(cor_table$comp)*0.95, xend=min(cor_table$comp)*0.95, y=min(cor_table$prevalence), yend=max(cor_table$prevalence)*1.05,
           arrow=arrow(length=unit(0.3, "cm"), type="closed")) +
  geom_point(size=3) +
  ylab("Prevalence in metagenomes") +
  xlab("Genetic incompatibility with resident bacteria\n(distance between 5mer distributions)") +
  theme_minimal() +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values=arg_category_palette, name = "") +
  theme(axis.text.x=element_text(size = rel(1.3)),
        axis.text.y=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1.2)),
        plot.tag = element_text(size=rel(2.5), face = "bold"),
        plot.title = element_text(size = rel(1.5))) 

pdf("figure_s8.pdf", height = 7, width = 9.5)
plot(fig_s8)
dev.off()
